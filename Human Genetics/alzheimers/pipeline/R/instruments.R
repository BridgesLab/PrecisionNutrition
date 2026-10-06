# instruments.R — instrument selection, harmonisation and the SNP attrition log.
#
# Three instrument sets per exposure, one per method family:
#   classical : p < 5e-8, clump r2 < 0.001 / 10 Mb      (the prior univariable pipeline's criteria)
#   mraid     : HapMap3, p < 5e-8, light clump r2 < 0.5  (MRAID wants correlated candidates; the
#               light clump only brings very polygenic exposures under its 10,000-SNP cap)
#   cause     : built inside cause.R (genome-wide parameters + p < 1e-3, r2 < 0.01 pruning)
#
# Every function returns the SNPs it dropped, by step, so the attrition table can say exactly where
# each instrument was lost. Clumping runs locally with plink2 on the same 1000G EUR panel the
# OpenGWAS LD server uses (deviation from the prior pipeline, forced by the /ld/clump outage).

suppressPackageStartupMessages({ library(data.table); library(dplyr) })

# One attrition record. `snps` are the variants REMOVED at this step.
attr_row <- function(assoc, set, step, n_after, snps = character())
  tibble(association = assoc, instrument_set = set, step = step, n_after = as.integer(n_after),
         n_removed = length(snps),
         removed = if (length(snps)) paste(sort(unique(snps)), collapse = ";") else NA_character_)

# plink2 clumping. Returns the index SNPs. p1 = p2 = threshold, as ieugwasr::ld_clump() does.
plink_clump <- function(d, bfile, p_thresh, r2, kb, plink2 = "plink2", work = tempdir()) {
  d <- d[p < p_thresh]
  if (!nrow(d)) return(character())
  f <- tempfile("clump_", tmpdir = work)
  fwrite(d[, .(ID = SNP, P = p)], paste0(f, ".tsv"), sep = "\t")
  # system2() does not quote arguments, and this project's path contains a space.
  args <- c("--bfile", shQuote(bfile), "--clump", shQuote(paste0(f, ".tsv")),
            "--clump-id-field", "ID", "--clump-p-field", "P", "--clump-p1", p_thresh,
            "--clump-p2", p_thresh, "--clump-r2", r2, "--clump-kb", kb, "--out", shQuote(f),
            "--silent", "--threads", 2)
  status <- system2(plink2, args, stdout = FALSE, stderr = FALSE)
  out <- paste0(f, ".clumps")
  if (status != 0 || !file.exists(out)) stop("plink2 --clump failed (status ", status, ")")
  fread(out, select = "ID")$ID
}

# Classical set: genome-wide significant, independent.
select_classical <- function(g, key, cfg, bfile, plink2) {
  sig <- g[p < cfg$p]
  idx <- plink_clump(sig, bfile, cfg$p, cfg$r2, cfg$kb, plink2)
  inst <- sig[SNP %chin% idx]
  list(inst = inst,
       attrition = bind_rows(
         attr_row(key, "classical", paste0("p < ", cfg$p), nrow(sig)),
         attr_row(key, "classical", sprintf("clump r2 < %s, %s kb (local plink2, 1000G EUR)",
                                            cfg$r2, cfg$kb), nrow(inst), setdiff(sig$SNP, idx))))
}

# MRAID candidate set: HapMap3, genome-wide significant, light clump; if the result is still above
# the cap, step down through `r2_ladder` and record which r2 was used.
select_mraid <- function(g, key, cfg, bfile, plink2, hm3) {
  sig <- g[p < cfg$p]
  h <- sig[SNP %chin% hm3]
  att <- bind_rows(attr_row(key, "mraid", paste0("p < ", cfg$p), nrow(sig)),
                   attr_row(key, "mraid", "HapMap3", nrow(h), setdiff(sig$SNP, h$SNP)))
  for (r2 in cfg$r2_ladder) {
    idx <- plink_clump(h, bfile, cfg$p, r2, cfg$kb, plink2)
    if (length(idx) <= cfg$cap) break
  }
  inst <- h[SNP %chin% idx]
  att <- bind_rows(att, attr_row(key, "mraid",
    sprintf("clump r2 < %s, %s kb%s", r2, cfg$kb,
            if (r2 != cfg$r2_ladder[1]) sprintf(" (stepped down from %s to meet cap %d)",
                                                cfg$r2_ladder[1], cfg$cap) else ""),
    nrow(inst), setdiff(h$SNP, idx)))
  if (nrow(inst) > cfg$cap) warning(key, ": MRAID candidate set still above cap (", nrow(inst), ")")
  list(inst = inst, attrition = att, r2_used = r2)
}

# ---- outcome versions ---------------------------------------------------------------------------

# SlopeHunter-adjusted outcome, genome-wide: beta - b_SH * beta_life, aligned to the outcome's
# effect allele. Variants without a lifespan association are dropped (logged by the caller).
# Lifespan is a Martingale residual increasing with MORTALITY; b_SH was fit on the same axis, so
# no sign change is applied here (see exposure_ad_indexevent.qmd).
adjust_outcome_gw <- function(y, life, b_sh) {
  m <- merge(y, life[, .(SNP, ea_l = ea, oa_l = oa, beta_l = beta, se_l = se)], by = "SNP")
  same <- m$ea == m$ea_l & m$oa == m$oa_l
  flip <- m$ea == m$oa_l & m$oa == m$ea_l
  m <- m[same | flip]
  m[, beta_l := fifelse(ea == ea_l, beta_l, -beta_l)]
  m[, `:=`(beta = beta - b_sh * beta_l, se = sqrt(se^2 + b_sh^2 * se_l^2))]
  m[, p := 2 * pnorm(-abs(beta / se))]
  m[, c("ea_l", "oa_l", "beta_l", "se_l") := NULL]
  setattr(m, "n_input", nrow(y))
  m[]
}

# ---- harmonisation --------------------------------------------------------------------------------

# Exposure instruments x outcome, through TwoSampleMR::harmonise_data(action = 2): the prior
# pipeline's rule (infer palindromes from frequency, drop the ambiguous ones). Returns the
# harmonised frame plus attrition for "missing in outcome" and each harmonisation drop reason.
harmonise_pair <- function(inst, out, key, set, action = 2, exp_name = "exposure",
                           out_name = "outcome") {
  o <- out[SNP %chin% inst$SNP]
  att <- attr_row(key, set, "present in outcome GWAS", nrow(o), setdiff(inst$SNP, o$SNP))
  if (!nrow(o)) return(list(dat = NULL, attrition = att))
  ex <- data.frame(SNP = inst$SNP, beta.exposure = inst$beta, se.exposure = inst$se,
                   effect_allele.exposure = inst$ea, other_allele.exposure = inst$oa,
                   eaf.exposure = inst$eaf, pval.exposure = inst$p, samplesize.exposure = inst$n,
                   exposure = exp_name, id.exposure = exp_name, mr_keep.exposure = TRUE)
  ou <- data.frame(SNP = o$SNP, beta.outcome = o$beta, se.outcome = o$se,
                   effect_allele.outcome = o$ea, other_allele.outcome = o$oa,
                   eaf.outcome = o$eaf, pval.outcome = o$p, samplesize.outcome = o$n,
                   outcome = out_name, id.outcome = out_name, mr_keep.outcome = TRUE)
  h <- suppressMessages(TwoSampleMR::harmonise_data(ex, ou, action = action))
  kept <- h$SNP[h$mr_keep]
  pal <- h$SNP[!h$mr_keep & h$palindromic %in% TRUE & h$ambiguous %in% TRUE]
  # Anything else that went in and did not come out kept: incompatible alleles, which
  # harmonise_data() either flags or removes outright.
  other <- setdiff(o$SNP, c(kept, pal))
  h <- h[h$mr_keep, ]
  att <- bind_rows(att,
    attr_row(key, set, "harmonise: ambiguous palindromic dropped", nrow(o) - length(pal), pal),
    attr_row(key, set, "harmonise: allele mismatch dropped", nrow(h), other))
  list(dat = h, attrition = att)
}

# plink2 binary: PLINK2_BIN env var (how Great Lakes jobs pass the module's copy) > config > PATH.
# The Intel builds in /usr/local/bin and genetics.binaRies do not run on this Mac (no Rosetta), so
# the config points at the arm64 build in ~/bin.
resolve_plink2 <- function(cfg_path = NULL) {
  for (cand in c(Sys.getenv("PLINK2_BIN"), cfg_path, unname(Sys.which("plink2")))) {
    if (!is.null(cand) && nzchar(cand)) {
      cand <- path.expand(cand)
      ok <- tryCatch(any(grepl("PLINK v2", suppressWarnings(system2(cand, "--version", stdout = TRUE, stderr = TRUE)))),
                     error = function(e) FALSE)
      if (ok) return(cand)
    }
  }
  stop("no working plink2 found: set PLINK2_BIN or paths$plink2 in pipeline/config.yml")
}

# apss.R — MR-APSS (Hu et al. 2022, PNAS), the pipeline's MAIN method.
#
# Category 1 (models correlated pleiotropy through a background genetic-correlation term, Omega)
# AND models sample overlap / population stratification through C, the matrix of LDSC
# intercepts. So, unlike MRAID, one method covers every association whether or not the two GWAS
# share UK Biobank. Mirrors calcium-cholesterol/robust_mr_apss.qmd, which ran successfully:
#   - QC to MR-APSS's published criteria (HapMap3, MAF > 0.05, unambiguous strand, no MHC,
#     chi2 < max(N/1000, 80)); see calcium-cholesterol/R/robust_mr_helpers.R::qc_gwas()
#   - C and Omega by bivariate LDSC (MRAPSS::est_paras) on the eur_w_ld_chr LD scores
#   - instruments at p < 5e-5 with Cor.SelectionBias = TRUE (winner's-curse correction), clumped
#     at r2 < 0.001 / 1 Mb, here with local plink2 (the calcium arm's plink 1.9 clumper is an
#     Intel binary that no longer runs on this Mac)
# Scale: per SD of exposure on per SD of outcome (MR-APSS works on z / sqrt(N)).

suppressPackageStartupMessages({ library(data.table); library(dplyr) })

# Standard format -> MR-APSS input (SNP A1 A2 Z P N chi2), with a QC log.
apss_input <- function(g, hm3, n_fallback, label, maf_min = 0.05) {
  d <- copy(g)[, n := fifelse(is.na(n), n_fallback, n)]
  steps <- list(); note <- function(s) steps[[length(steps) + 1]] <<- tibble(dataset = label, step = s, n_variants = nrow(d))
  note("pipeline QC output")
  d <- d[!paste0(ea, oa) %chin% c("AT", "TA", "CG", "GC")];      note("unambiguous strand")
  d <- d[SNP %chin% hm3];                                         note("HapMap3")
  d <- d[is.na(eaf) | pmin(eaf, 1 - eaf) > maf_min];              note(paste0("MAF > ", maf_min))
  d <- d[!(chr == "6" & pos > 26e6 & pos < 34e6)];                note("MHC removed")
  d[, z := beta / se]
  chi2_max <- max(stats::median(d$n) / 1000, 80)
  d <- d[z^2 < chi2_max];                                         note(sprintf("chi2 < %.0f", chi2_max))
  out <- as.data.frame(d[, .(SNP, A1 = ea, A2 = oa, Z = z, P = p, N = n, chi2 = z^2)])
  attr(out, "qc_log") <- bind_rows(steps)
  out
}

run_apss <- function(g_exp, g_out, key, hm3, ldsc, n_exp, n_out, bfile, plink2, cfg) {
  e <- apss_input(g_exp, hm3, n_exp, "exposure")
  o <- apss_input(g_out, hm3, n_out, "outcome")
  # est_paras_jk() == MRAPSS::est_paras() plus the LDSC jackknife delete-one values (latent.R),
  # needed for the background-slope interval. C and Omega are identical to the package's.
  paras <- est_paras_jk(dat1 = e, dat2 = o, trait1.name = "exposure", trait2.name = "outcome",
                        ld = ldsc$ld, M = ldsc$M)
  label <- cfg$label %||% "MR-APSS"
  set <- if (label == "MR-APSS") "mrapss" else paste0("mrapss_", format(cfg$iv_threshold, scientific = TRUE))
  d <- as.data.table(paras$dat)
  att <- bind_rows(
    attr(e, "qc_log") |> transmute(association = key, instrument_set = set,
                                   step = paste("exposure QC:", step), n_after = as.integer(n_variants)),
    attr_row(key, set, "merged with outcome (est_paras)", nrow(d)))
  cand <- d[pval.exp < cfg$iv_threshold]
  idx <- plink_clump(data.table(SNP = cand$SNP, p = cand$pval.exp), bfile, cfg$iv_threshold,
                     cfg$clump_r2, cfg$clump_kb, plink2)
  att <- bind_rows(att,
    attr_row(key, set, paste0("p < ", cfg$iv_threshold), nrow(cand)),
    attr_row(key, set, sprintf("clump r2 < %s, %s kb", cfg$clump_r2, cfg$clump_kb),
             length(idx), setdiff(cand$SNP, idx)))
  MRdat <- sanitise_mrdat(as.data.frame(cand[SNP %chin% idx]), cfg$iv_threshold)
  fit <- MRAPSS::MRAPSS(MRdat, exposure = "exposure", outcome = "outcome", C = paras$C,
                        Omega = paras$Omega, Cor.SelectionBias = isTRUE(cfg$cor_selection_bias))
  td <- tidy_apss(fit, exposure = key, outcome = key, arm = key)
  bg <- tidy_background(paras, exposure = key, outcome = key, arm = key)
  res <- tibble(method = label, scale = "SD(outcome) per SD(exposure)",
                b = td$b, se = td$se, ci_lo = td$lci, ci_hi = td$uci, p = td$pval,
                n_snp = td$n_iv, apss_n_valid_iv = td$n_valid_iv, apss_pi0 = td$pi0,
                apss_tau2 = td$tau2, apss_C12 = bg$C12, apss_rg = bg$rg,
                apss_overlap_flag = bg$overlap_flag,
                apss_threshold = cfg$iv_threshold) |>
    bind_cols(apss_background_slope(paras, n_exp, n_out))
  inst <- as_tibble(MRdat) |> transmute(SNP, effect_allele = A1, other_allele = A2,
                                        beta_exposure = b.exp, se_exposure = se.exp,
                                        p_exposure = pval.exp, beta_outcome = b.out,
                                        se_outcome = se.out)
  list(result = res, instruments = inst, attrition = att, background = bg, set = set)
}

analyse_apss <- function(g_exp, g_out, assoc_id, exposure, outcome, arm, hm3, ldsc, n_exp, n_out,
                         bfile, plink2, cfg) {
  r <- run_apss(g_exp, g_out, assoc_id, hm3, ldsc, n_exp, n_out, bfile, plink2, cfg)
  list(rows = bind_cols(meta_cols(assoc_id, exposure, outcome, arm), r$result),
       instruments = r$instruments |> mutate(assoc_id = assoc_id, instrument_set = r$set, .before = 1),
       attrition = r$attrition,
       background = r$background |> mutate(assoc_id = assoc_id, .before = 1))
}

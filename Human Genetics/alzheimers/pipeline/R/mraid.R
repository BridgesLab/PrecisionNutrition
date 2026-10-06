# mraid.R — MRAID (Yuan et al. 2022, Sci Adv; PMID 35235357) on correlated candidate instruments.
#
# Inputs MRAID needs: exposure and outcome z-scores for the candidate SNPs, and the LD matrix among
# them (from the 1000G EUR panel for both, as the paper does when individual data are unavailable).
# MRAID assumes the two GWAS are from NON-OVERLAPPING samples; the pipeline only makes it primary
# where that holds (see primary_method() in pipeline.R).
#
# Scale: MRAID works on standardised effects (beta = Z / sqrt(N - 1)), so its causal effect is in
# SD(outcome) per SD(exposure), observed scale for binary traits. Compare it with the classical
# estimates on sign and p, not magnitude.

suppressPackageStartupMessages({ library(data.table); library(dplyr) })

# Signed LD (r) among `snps` from the plink panel, one block per chromosome. Cross-chromosome LD is
# set to 0 rather than estimated from 503 people, where it would be pure noise. Returns the matrix
# ordered as `snps` (SNPs missing from the panel are dropped and reported via attr "missing").
# r is computed on allele counts of bim column 5; callers must align z-scores to that allele.
ld_matrix <- function(snps, bim, bfile, plink2, work = tempdir()) {
  b <- bim[SNP %chin% snps]
  R <- matrix(0, nrow(b), nrow(b), dimnames = list(b$SNP, b$SNP))
  for (ch in unique(b$chr)) {
    ids <- b[chr == ch, SNP]
    if (length(ids) == 1) { R[ids, ids] <- 1; next }
    f <- tempfile(paste0("ld_chr", ch, "_"), tmpdir = work)
    writeLines(ids, paste0(f, ".snps"))
    status <- system2(plink2, c("--bfile", shQuote(bfile), "--extract", shQuote(paste0(f, ".snps")),
                                "--r-unphased", "square", "--out", shQuote(f), "--silent",
                                "--threads", 2), stdout = FALSE, stderr = FALSE)
    mf <- paste0(f, ".unphased.vcor1")
    if (status != 0 || !file.exists(mf)) stop("plink2 --r-unphased failed on chr ", ch)
    order_ids <- fread(paste0(mf, ".vars"), header = FALSE)[[1]]
    m <- as.matrix(fread(mf, header = FALSE))
    m[!is.finite(m)] <- 0                       # monomorphic in the panel
    R[order_ids, order_ids] <- m
    unlink(list.files(dirname(f), pattern = basename(f), full.names = TRUE))
  }
  diag(R) <- 1
  keep <- intersect(snps, rownames(R))
  out <- R[keep, keep, drop = FALSE]
  attr(out, "missing") <- setdiff(snps, keep)
  out
}

# Thin wrapper over MRAID's internal sampler. Identical to MRAID::MRAID() (package version 1.0),
# including its IVW starting value and its re-run with maxvar = 0 when the first pass returns
# exactly 0, but it also returns the posterior SD that MRAID() computes and then discards.
mraid_fit <- function(zx, zy, R, n1, n2, gibbs = 1000, burnin = 0.2, seed = 20261004) {
  set.seed(seed)
  bx <- zx / sqrt(n1 - 1); by <- zy / sqrt(n2 - 1)
  # IVW start on a near-independent subset (|r| <= 0.25), same greedy rule as MRAID(): repeatedly
  # drop the SNP with the most |r| > 0.25 partners (first such SNP on ties). MRAID() recomputes the
  # counts on a copy of the matrix after every drop, which is O(p^3) and takes hours at p ~ 5,000;
  # here the counts are updated in place, which selects exactly the same SNPs.
  hit <- abs(R) > 0.25; diag(hit) <- FALSE
  alive <- rep(TRUE, ncol(R))
  cnt <- colSums(hit)
  while (any(cnt[alive] > 0)) {
    worst <- which(alive)[which.max(cnt[alive])]
    alive[worst] <- FALSE
    cnt <- cnt - hit[, worst]
  }
  keep <- which(alive)
  alpha0 <- stats::coef(stats::lm(by[keep] ~ -1 + bx[keep]))[[1]]
  run <- function(maxvar) MRAID::MRAID_CPP(
    bx, by, R, R, n1, n2, Gibbsnumberin = gibbs, burninproportion = burnin,
    initial_betain = rep(0, length(bx)), pi_beta_shape_in = 0.5, pi_beta_scale_in = 4.5,
    pi_c_shape_in = 0.5, pi_c_scale_in = 9.5, pi_1_shape_in = 0.5, pi_1_scale_in = 1.5,
    pi_0_shape_in = 0.05, pi_0_scale_in = 9.95, maxvarin = maxvar, alphain = alpha0)
  re <- run(1000)
  if (re$alpha == 0) re <- run(0)
  tibble(method = "MRAID", scale = "SD(outcome) per SD(exposure)",
         b = re$alpha, se = re$sd, p = 2 * stats::pnorm(-abs(re$alpha / re$sd)),
         n_snp = length(bx), mraid_rho = re$rho, mraid_sigmabeta = re$sigmabeta,
         mraid_sigmaeta = re$sigmaeta, mraid_sigma2x = re$sigma2x, mraid_sigma2y = re$sigma2y,
         mraid_n_independent_start = length(keep))
}

# Candidate set -> harmonised to each other (prior rule, action = 2) -> aligned to the panel's
# column-5 allele -> LD matrix -> MRAID. Returns the estimate, the instrument table and attrition.
run_mraid <- function(inst, out, key, bim, bfile, plink2, n_exp, n_out, cfg) {
  h <- harmonise_pair(inst, out, key, "mraid")
  if (is.null(h$dat) || nrow(h$dat) < 3)
    return(list(result = tibble(method = "MRAID", note = "not run: < 3 harmonised SNPs"),
                instruments = NULL, attrition = h$attrition))
  dat <- as.data.table(h$dat)[, .(SNP, ea = effect_allele.exposure, oa = other_allele.exposure,
                                  bx = beta.exposure, sx = se.exposure, by = beta.outcome,
                                  sy = se.outcome, px = pval.exposure)]
  dat <- merge(dat, bim[, .(SNP, chr, pos, a1, a2)], by = "SNP")
  same <- dat$ea == dat$a1 & dat$oa == dat$a2
  flip <- dat$ea == dat$a2 & dat$oa == dat$a1
  mismatch <- dat$SNP[!(same | flip)]
  dat <- dat[same | flip]
  dat[ea == a2, `:=`(bx = -bx, by = -by)]
  dat[, `:=`(zx = bx / sx, zy = by / sy)]
  setorder(dat, chr, pos)
  R <- ld_matrix(dat$SNP, bim, bfile, plink2)
  dat <- dat[SNP %chin% rownames(R)]
  att <- bind_rows(h$attrition,
    attr_row(key, "mraid", "aligned to LD panel alleles", nrow(dat), c(mismatch, attr(R, "missing"))))
  res <- mraid_fit(dat$zx, dat$zy, R[dat$SNP, dat$SNP], n_exp, n_out,
                   gibbs = cfg$gibbs, burnin = cfg$burnin)
  list(result = res, instruments = dat[, .(SNP, chr, pos, effect_allele = a1, other_allele = a2,
                                           beta_exposure = bx, se_exposure = sx, p_exposure = px,
                                           beta_outcome = by, se_outcome = sy)],
       attrition = att)
}

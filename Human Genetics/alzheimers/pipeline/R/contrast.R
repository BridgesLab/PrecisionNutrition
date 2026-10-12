# contrast.R — fasting vs 2 h glucose -> AD, Step 0: instruments, overlap and the conditional-F gate.
#
# OUTCOME-BLIND. Nothing in this file reads an AD GWAS. Instruments are selected on a T2D GWAS and
# their effects on fasting glucose (FG) and 2 h glucose (2hGlu) come from the glycemic GWAS, so the
# selection p-value and the exposure betas are from different studies (three-sample design).
# The decision rule (config `contrast: cond_f_min`) is applied here and written to
# results/pipeline/contrast/ before any outcome is merged:
#   conditional F for 2hGlu >= cond_f_min  -> MVMR is primary
#   otherwise                              -> contrast regression is primary, MVMR shown as underpowered
#
# Scale: betas stay in each GWAS's own unit (mmol/L for the MAGIC files). Conditional F does not
# depend on the scale; the strength of the (2hGlu - FG) contrast does, so it is given both per mmol/L
# and per SD (config `contrast: sd_mmol`).

suppressPackageStartupMessages({ library(data.table); library(dplyr) })

# T2D-selected instruments: genome-wide significant, clumped. `restrict` (rsIDs) limits the
# candidates first, e.g. to the published index variants that carry a cluster assignment.
contrast_select <- function(g_t2d, key, cfg, bfile, plink2, restrict = NULL) {
  att <- NULL
  if (!is.null(restrict)) {
    n0 <- nrow(g_t2d)
    g_t2d <- g_t2d[SNP %chin% restrict]
    att <- attr_row(key, "contrast", sprintf("restricted to %d listed index variants", length(restrict)),
                    nrow(g_t2d), setdiff(restrict, g_t2d$SNP)) |>
      mutate(n_removed = n0 - nrow(g_t2d))          # `removed` lists the index variants NOT found
  }
  sel <- select_classical(g_t2d, key, cfg, bfile, plink2)
  sel$attrition <- bind_rows(att, sel$attrition |> mutate(instrument_set = "contrast"))
  sel
}

# Instruments x FG x 2hGlu on one allele: the T2D risk-increasing allele. Reuses harmonise_pair()
# (TwoSampleMR action 2), so palindromic and mismatched variants are dropped and logged as elsewhere.
contrast_harmonise <- function(inst, g_fg, g_2h, key) {
  one <- function(g, nm) {
    h <- harmonise_pair(inst, g, key, paste0("contrast: ", nm), exp_name = "t2d", out_name = nm)
    d <- if (is.null(h$dat)) NULL else
      as_tibble(h$dat) |> transmute(SNP, ea = effect_allele.exposure, oa = other_allele.exposure,
                                    eaf = eaf.exposure, beta_t2d = beta.exposure, se_t2d = se.exposure,
                                    p_t2d = pval.exposure, beta = beta.outcome, se = se.outcome,
                                    p = pval.outcome)
    list(dat = d, attrition = h$attrition)
  }
  f <- one(g_fg, "fg"); t <- one(g_2h, "glu2h")
  if (is.null(f$dat) || is.null(t$dat)) stop(key, ": no instruments found in an exposure GWAS")
  d <- inner_join(f$dat |> rename(beta_fg = beta, se_fg = se, p_fg = p),
                  t$dat |> select(SNP, beta_2h = beta, se_2h = se, p_2h = p), by = "SNP") |>
    left_join(as_tibble(inst) |> select(SNP, chr, pos), by = "SNP")
  flip <- d$beta_t2d < 0
  d <- d |> mutate(ea0 = ea, ea = if_else(flip, oa, ea), oa = if_else(flip, ea0, oa),
                   eaf = if_else(flip, 1 - eaf, eaf),
                   across(c(beta_t2d, beta_fg, beta_2h), \(b) if_else(flip, -b, b))) |>
    select(SNP, chr, pos, ea, oa, eaf, starts_with("beta_"), starts_with("se_"), starts_with("p_"))
  list(dat = d, attrition = bind_rows(f$attrition, t$attrition,
         attr_row(key, "contrast", "in both exposure GWAS, aligned to the T2D risk allele", nrow(d),
                  setdiff(union(f$dat$SNP, t$dat$SNP), d$SNP))))
}

# Estimation-error correlations between the three GWAS (null-variant z-score correlation,
# classical.R::error_cor(); approximates the cross-trait LDSC intercept). FG-2hGlu is the one the
# conditional F and the contrast's SE need; the two T2D rows say how far selection on T2D is
# independent of the noise in the exposure betas (winner's curse).
contrast_overlap <- function(g_t2d, g_fg, g_2h, set) {
  r <- function(a, b) error_cor(a, b)[1, 2]
  tibble(set = set, pair = c("fg-glu2h", "t2d-fg", "t2d-glu2h"),
         error_cor = c(r(g_fg, g_2h), r(g_t2d, g_fg), r(g_t2d, g_2h)))
}

# Sanderson-Windmeijer conditional F (MVMR::strength_mvmr), with the FG/2hGlu error covariance
# rho * se_fg * se_2h for each SNP. BYG is a placeholder (strength does not use the outcome); the T2D
# beta fills it so that no outcome data are needed.
conditional_f <- function(d, rho) {
  se <- cbind(d$se_fg, d$se_2h)
  inp <- MVMR::format_mvmr(BXGs = cbind(d$beta_fg, d$beta_2h), BYG = d$beta_t2d, seBXGs = se,
                           seBYG = d$se_t2d, RSID = d$SNP)
  gencov <- MVMR::phenocov_mvmr(matrix(c(1, rho, rho, 1), 2), se)
  utils::capture.output(s <- MVMR::strength_mvmr(inp, gencov))      # it prints its own table
  c(fg = s[1, 1], glu2h = s[1, 2])
}

# Mean z^2 of the per-SNP contrast (2hGlu/s_2h - FG/s_fg): the first-stage strength of the regressor
# the contrast regression uses. s = 1 gives the mmol/L contrast.
contrast_f <- function(d, rho, s_fg = 1, s_2h = 1) {
  x <- d$beta_2h / s_2h - d$beta_fg / s_fg
  v <- (d$se_2h / s_2h)^2 + (d$se_fg / s_fg)^2 - 2 * rho * d$se_2h * d$se_fg / (s_2h * s_fg)
  mean(x^2 / v)
}

contrast_strength <- function(d, rho, set, cfg) {
  cf <- conditional_f(d, rho)
  sd <- cfg$sd_mmol
  primary <- if (cf[["glu2h"]] >= cfg$cond_f_min) "MVMR" else "contrast regression"
  tibble(set = set, n_snp = nrow(d),
         n_nominal_fg = sum(d$p_fg < 0.05), n_nominal_2h = sum(d$p_2h < 0.05),
         n_gws_fg = sum(d$p_fg < 5e-8), n_gws_2h = sum(d$p_2h < 5e-8),
         mean_F_fg = mean((d$beta_fg / d$se_fg)^2), mean_F_2h = mean((d$beta_2h / d$se_2h)^2),
         cor_beta_fg_2h = cor(d$beta_fg, d$beta_2h), error_cor_fg_2h = rho,
         cond_F_fg = cf[["fg"]], cond_F_2h = cf[["glu2h"]],
         contrast_F_mmol = contrast_f(d, rho),
         contrast_F_sd = contrast_f(d, rho, sd$fg, sd$glu2h),
         cond_f_min = cfg$cond_f_min, primary_analysis = primary)
}

# What the target calls: one instrument set end to end.
run_contrast_step0 <- function(g_t2d, g_fg, g_2h, set, cfg, bfile, plink2, restrict = NULL) {
  sel <- contrast_select(g_t2d, set, cfg$instruments, bfile, plink2, restrict)
  h <- contrast_harmonise(sel$inst, g_fg, g_2h, set)
  ov <- contrast_overlap(g_t2d, g_fg, g_2h, set)
  rho <- ov$error_cor[ov$pair == "fg-glu2h"]
  list(instruments = h$dat |> mutate(set = set, .before = 1),
       strength = contrast_strength(h$dat, rho, set, cfg),
       overlap = ov,
       attrition = bind_rows(sel$attrition, h$attrition))
}

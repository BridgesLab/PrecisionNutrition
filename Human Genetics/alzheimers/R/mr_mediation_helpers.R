# mr_mediation_helpers.R — two-step + MVMR mediation engine for exposure -> mediator -> AD triads.
#
# Consumed by mr_mediation.qmd. Everything that touches OpenGWAS is in section 1; everything after
# that is pure computation on harmonised tables, so it can be tested offline on simulated data.
#
# Estimands (Carter et al. 2021, Eur J Epidemiol, PMID 33961203):
#   total    = X -> Y                 (univariable, method_XY_total, X's own instruments)
#   xm       = X -> M                 (univariable, method_XM, X's own instruments)
#   direct   = X -> Y | M             (MVMR, method_MVMR, re-clumped X u M instruments)
#   my       = M -> Y | X             (same MVMR fit)
#   indirect = xm * my                (product of coefficients)
#
# Survival-collider adjustment is applied at the SNP level, before any leg is fit:
#   beta_adj = beta_Y - b_SH * beta_lifespan,   se_adj = sqrt(se_Y^2 + b_SH^2 * se_lifespan^2)
# This is the SNP-level form of correct_estimate() in R/fit_selection_slope.R (identical for IVW,
# and the only form that also works for MVMR and the robust methods). The uncertainty in b_SH is
# NOT folded into se_adj: it is shared by every SNP, so it is propagated by redrawing b_SH in the
# Monte Carlo (section 5) instead.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- 0. Defaults inherited from the prior univariable pipeline ----------------------------------
# Source: exposure_survival_screen.qmd / exposure_ad_indexevent.qmd (P_THRESH, R2, KB; harmonise
# action = 2; extract_outcome_data() called with package defaults, i.e. LD proxies ON at r2 >= 0.8).
# Clumping goes through the OpenGWAS LD server (1000G EUR), as extract_instruments() does by default.
# The prior analyses applied NO F-statistic filter (mean F was reported only), so f_min = NA here;
# setting it is a deviation and is logged.
INSTRUMENT_DEFAULTS <- list(
  p_thresh            = 5e-8,
  clump_r2            = 0.001,
  clump_kb            = 10000,
  ld_pop              = "EUR",
  harmonise_action    = 2,       # infer palindromes from EAF, drop the ambiguous ones
  proxies             = TRUE,
  proxy_rsq           = 0.8,
  proxy_palindromes   = 1,
  proxy_maf_threshold = 0.3,
  f_min               = NA_real_
)

# Merge a per-row override (named list, possibly empty) onto the defaults; return the parameters
# and a log of what changed, so every deviation is recorded rather than remembered.
resolve_instrument_params <- function(override = list(), triad_id = NA_character_) {
  override <- override %||% list()
  bad <- setdiff(names(override), names(INSTRUMENT_DEFAULTS))
  if (length(bad)) stop("unknown instrument parameter(s) for ", triad_id, ": ",
                        paste(bad, collapse = ", "))
  par <- utils::modifyList(INSTRUMENT_DEFAULTS, override)
  changed <- names(override)[!map2_lgl(override, INSTRUMENT_DEFAULTS[names(override)], identical)]
  log <- if (length(changed))
    tibble(triad_id = triad_id, outcome_id = NA_character_, stage = "instrument_params",
           type = "deviation",
           detail = paste0(changed, ": ", map_chr(INSTRUMENT_DEFAULTS[changed], format),
                           " -> ", map_chr(par[changed], format)),
           n = NA_integer_, snps = NA_character_)
  else empty_log()
  list(par = par, log = log)
}

# Every step returns its log rows instead of writing to a global, so knitr caching and re-runs
# can never silently drop or duplicate log entries.
empty_log <- function()
  tibble(triad_id = character(), outcome_id = character(), stage = character(),
         type = character(), detail = character(), n = integer(), snps = character())

log_row <- function(triad_id, outcome_id, stage, type, detail, snps = character()) {
  # SNP-list events with no SNPs are not events; deviations and notes are always kept.
  if (!length(snps) && !type %in% c("deviation", "note")) return(empty_log())
  tibble(triad_id = triad_id %||% NA_character_, outcome_id = outcome_id %||% NA_character_,
         stage = stage, type = type, detail = detail,
         n = as.integer(length(snps)),
         snps = if (length(snps)) paste(sort(unique(snps)), collapse = ";") else NA_character_)
}

# ---- 1. OpenGWAS access (cached to disk) ---------------------------------------------------------

# Retry a (possibly rate-limited) OpenGWAS call with backoff. `fn` is a thunk. Same as Step 2,
# except the last error is kept (og_last_error()) so a failure can say WHY: a 502 from the server
# and an expired token otherwise look identical to the caller.
.og_state <- new.env(parent = emptyenv())
og_last_error <- function() .og_state$last %||% "no error captured (call returned an empty result)"
og_retry <- function(fn, tries = 4, base = 3) {
  .og_state$last <- NULL
  for (i in seq_len(tries)) {
    r <- tryCatch(fn(), error = function(e) {
      msg <- conditionMessage(e)
      code <- regmatches(msg, regexpr("Status code from OpenGWAS API: [0-9]+", msg))
      .og_state$last <- if (length(code)) code else substr(gsub("\\s+", " ", msg), 1, 200)
      NULL })
    if (!is.null(r) && (!is.data.frame(r) || nrow(r) > 0)) return(r)
    if (i < tries) Sys.sleep(base * i)
  }
  NULL
}

# Disk cache keyed on a human-readable name plus a hash of the inputs. knitr caches key on chunk
# code only (the calcium arm was burned by this), so the key must include every parameter.
# A NULL / empty result is never written, so a failed pull is retried on the next render.
cache_rds <- function(dir, name, key, expr, refresh = FALSE) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  h <- substr(digest_str(key), 1, 10)
  path <- file.path(dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", name), "_", h, ".rds"))
  if (!refresh && file.exists(path)) return(readRDS(path))
  val <- expr
  if (!is.null(val) && (!is.data.frame(val) || nrow(val) > 0)) saveRDS(val, path)
  val
}

digest_str <- function(x) {
  f <- tempfile(); on.exit(unlink(f))
  writeLines(paste(deparse(x), collapse = ""), f)
  unname(tools::md5sum(f))
}

pull_instruments <- function(id, par, cache_dir) {
  cache_rds(cache_dir, paste0("inst_", id), list(id, par$p_thresh, par$clump_r2, par$clump_kb),
            og_retry(function() TwoSampleMR::extract_instruments(
              id, p1 = par$p_thresh, clump = TRUE, r2 = par$clump_r2, kb = par$clump_kb)))
}

pull_snps <- function(snps, id, par, cache_dir) {
  snps <- sort(unique(snps))
  cache_rds(cache_dir, paste0("lookup_", id), list(snps, id, par[grep("^proxy", names(par))]),
            og_retry(function() TwoSampleMR::extract_outcome_data(
              snps = snps, outcomes = id, proxies = par$proxies, rsq = par$proxy_rsq,
              palindromes = par$proxy_palindromes, maf_threshold = par$proxy_maf_threshold)))
}

# Re-clump X u M with the same parameters and extract both exposures for the survivors, aligned
# to X's alleles. This is TwoSampleMR::mv_extract_exposures(), which does exactly that.
pull_mv_exposures <- function(x_id, m_id, par, cache_dir) {
  cache_rds(cache_dir, paste0("mvexp_", x_id, "_", m_id), list(x_id, m_id, par),
            og_retry(function() TwoSampleMR::mv_extract_exposures(
              c(x_id, m_id), clump_r2 = par$clump_r2, clump_kb = par$clump_kb,
              harmonise_strictness = par$harmonise_action, find_proxies = par$proxies,
              pval_threshold = par$p_thresh, pop = par$ld_pop)))
}

# Which SNPs came back as an LD proxy rather than themselves. extract_outcome_data() marks these
# with proxy.outcome = TRUE when proxies are on.
proxied_snps <- function(d) {
  if (is.null(d) || !"proxy.outcome" %in% names(d)) return(character())
  d$SNP[!is.na(d$proxy.outcome) & d$proxy.outcome]
}

# ---- 2. Outcome tables (one per outcome, over the union of every triad's instruments) ------------

# Y effects plus the lifespan (selection-axis) effects aligned to Y's effect allele. Built ONCE per
# outcome; every triad then reads from it, so all triads share one version of the outcome.
build_outcome_table <- function(outcome_id, snps, life_id, par, cache_dir) {
  lg <- empty_log()
  y <- pull_snps(snps, outcome_id, par, cache_dir)
  if (is.null(y)) stop("outcome pull failed: ", outcome_id)
  y <- y |> filter(mr_keep.outcome)
  lg <- bind_rows(lg,
    log_row(NA, outcome_id, "outcome_table", "missing_in_outcome",
            "instrument SNPs with no outcome association (after proxy search)",
            setdiff(snps, y$SNP)),
    log_row(NA, outcome_id, "outcome_table", "proxy_used", "outcome looked up via LD proxy",
            proxied_snps(y)))

  life <- pull_snps(y$SNP, life_id, par, cache_dir)
  life_h <- if (is.null(life)) NULL else {
    yx <- TwoSampleMR::convert_outcome_to_exposure(y)
    TwoSampleMR::harmonise_data(yx, life, action = par$harmonise_action) |> filter(mr_keep)
  }
  lg <- bind_rows(lg,
    log_row(NA, outcome_id, "outcome_table", "missing_in_lifespan",
            paste0("no ", life_id, " association, or dropped aligning lifespan to the outcome ",
                   "alleles; these SNPs are DROPPED from every SlopeHunter-adjusted arm"),
            setdiff(y$SNP, life_h$SNP %||% character())),
    log_row(NA, outcome_id, "outcome_table", "proxy_used", "lifespan looked up via LD proxy",
            proxied_snps(life)))

  ot <- y |>
    transmute(SNP, effect_allele = effect_allele.outcome, other_allele = other_allele.outcome,
              eaf = eaf.outcome, beta_y = beta.outcome, se_y = se.outcome,
              pval_y = pval.outcome, samplesize_y = samplesize.outcome)
  if (!is.null(life_h))
    ot <- ot |> left_join(life_h |> transmute(SNP, beta_life = beta.outcome, se_life = se.outcome),
                          by = "SNP")
  else ot <- ot |> mutate(beta_life = NA_real_, se_life = NA_real_)
  list(table = ot |> mutate(outcome_id = outcome_id, life_id = life_id), log = lg)
}

# Align the outcome table to an exposure frame's alleles. Rather than re-implement allele logic,
# harmonise a copy whose beta is 1: after harmonise_data() it holds the sign (+1 kept, -1 flipped)
# that TwoSampleMR applied to that SNP, and the same sign is then applied to BOTH beta_y and
# beta_life. Using one call guarantees Y and lifespan get identical keep/flip decisions.
align_outcome <- function(exp_dat, ot, par) {
  unit <- ot |>
    transmute(SNP, beta.outcome = 1, se.outcome = 1, effect_allele.outcome = effect_allele,
              other_allele.outcome = other_allele, eaf.outcome = eaf, pval.outcome = pval_y,
              outcome = outcome_id, id.outcome = outcome_id, mr_keep.outcome = TRUE)
  h <- suppressMessages(TwoSampleMR::harmonise_data(exp_dat, as.data.frame(unit),
                                                    action = par$harmonise_action))
  h <- h |> filter(mr_keep) |> select(SNP, flip = beta.outcome) |> distinct(SNP, .keep_all = TRUE)
  ot |> inner_join(h, by = "SNP") |>
    mutate(by = beta_y * flip, sy = se_y, bl = beta_life * flip, sl = se_life) |>
    select(SNP, by, sy, bl, sl, samplesize_y)
}

# SNP-level SlopeHunter adjustment. b_sh = 0 returns the raw outcome unchanged.
adjust_outcome <- function(by, sy, bl, sl, b_sh) {
  if (b_sh == 0) return(list(b = by, s = sy))
  list(b = by - b_sh * bl, s = sqrt(sy^2 + b_sh^2 * sl^2))
}

# ---- 3. Per-triad analysis frames ----------------------------------------------------------------

# Turn the raw pulls for one triad x outcome into the three aligned SNP tables the legs need:
#   uv : X's own instruments  -> bx, sx, by, sy, bl, sl            (total effect)
#   xm : X's own instruments in the M GWAS (TwoSampleMR harmonised frame)  (X -> M)
#   mv : re-clumped X u M     -> bx, sx, bm, sm, by, sy, bl, sl    (MVMR)
# All three are aligned to X's effect allele, which is what lets the SNP-level bootstrap in
# section 5 reuse one noise draw per SNP x GWAS across tables.
build_triad_frames <- function(pulls, ot, par, adjusted, common_snp_set, triad_id, outcome_id) {
  lg <- empty_log()
  x_inst <- pulls$x_inst

  uv <- align_outcome(x_inst, ot, par) |>
    inner_join(x_inst |> transmute(SNP, bx = beta.exposure, sx = se.exposure,
                                   pval_x = pval.exposure, samplesize_x = samplesize.exposure),
               by = "SNP")

  mvx <- pulls$mv_exp |> filter(id.exposure == pulls$x_id)
  mvm <- pulls$mv_exp |> filter(id.exposure == pulls$m_id)
  mv <- align_outcome(mvx, ot, par) |>
    inner_join(mvx |> transmute(SNP, bx = beta.exposure, sx = se.exposure, pval_x = pval.exposure),
               by = "SNP") |>
    inner_join(mvm |> transmute(SNP, bm = beta.exposure, sm = se.exposure, pval_m = pval.exposure),
               by = "SNP") |>
    mutate(in_x_inst = SNP %in% pulls$x_inst$SNP, in_m_inst = SNP %in% pulls$m_inst$SNP)

  lg <- bind_rows(lg,
    log_row(triad_id, outcome_id, "harmonise_total", "dropped",
            "X instruments lost aligning to the outcome (missing or ambiguous palindrome)",
            setdiff(x_inst$SNP, uv$SNP)),
    log_row(triad_id, outcome_id, "harmonise_mvmr", "dropped",
            "re-clumped X u M SNPs lost aligning to the outcome",
            setdiff(unique(pulls$mv_exp$SNP), mv$SNP)))

  # SlopeHunter arms need lifespan for every SNP. Missing lifespan -> drop (proxies were already
  # tried when the outcome table was built). With common_snp_set the raw arm drops the same SNPs,
  # so the adjusted-vs-raw contrast isolates the adjustment rather than a change of SNP set.
  if (adjusted || common_snp_set) {
    lost_uv <- uv$SNP[is.na(uv$bl)]; lost_mv <- mv$SNP[is.na(mv$bl)]
    why <- if (adjusted) "no lifespan effect; dropped from the SlopeHunter-adjusted arm"
           else "no lifespan effect; dropped from the RAW arm too (common_snp_set = TRUE)"
    lg <- bind_rows(lg,
      log_row(triad_id, outcome_id, "slopehunter_total", "dropped", why, lost_uv),
      log_row(triad_id, outcome_id, "slopehunter_mvmr",  "dropped", why, lost_mv))
    uv <- uv |> filter(!is.na(bl)); mv <- mv |> filter(!is.na(bl))
  }

  # X -> M uses its own lookup (X instruments in the M GWAS), independent of the outcome.
  xm <- suppressMessages(TwoSampleMR::harmonise_data(pulls$x_inst, pulls$x_in_m,
                                                     action = par$harmonise_action)) |>
    filter(mr_keep)
  lg <- bind_rows(lg,
    log_row(triad_id, NA, "harmonise_xm", "dropped", "X instruments lost in the mediator GWAS",
            setdiff(x_inst$SNP, xm$SNP)),
    log_row(triad_id, NA, "harmonise_xm", "proxy_used", "mediator looked up via LD proxy",
            proxied_snps(pulls$x_in_m)))

  # Optional F filter (a deviation: the prior pipeline had none).
  if (!is.na(par$f_min)) {
    weak <- uv$SNP[(uv$bx / uv$sx)^2 < par$f_min]
    lg <- bind_rows(lg, log_row(triad_id, outcome_id, "f_filter", "deviation",
                                paste0("univariable F < ", par$f_min, " removed"), weak))
    uv <- uv |> filter(!SNP %in% weak); xm <- xm |> filter(!SNP %in% weak)
  }
  list(uv = uv, xm = xm, mv = mv, log = lg)
}

# ---- 4. Estimators ---------------------------------------------------------------------------------

# Weighted regression through the origin, the engine shared by IVW and MV-IVW. `se_mode`:
#   "fe"     fixed effect (residual SE fixed at 1)
#   "mre"    multiplicative random effects, residual SE floored at 1 (as MendelianRandomization
#            mr_ivw/mr_mvivw do; NOTE TwoSampleMR 0.7.5's mr_ivw_mre does NOT floor, so under-
#            dispersed instruments get an SE smaller than fixed effects there)
#   "presso" residual SE unfloored, which is what MR-PRESSO's own lm() reports
# Returns the coefficient vector and its covariance matrix (needed jointly for direct and M|X).
ivw_core <- function(bx, by, sy, se_mode = "mre") {
  bx <- as.matrix(bx)
  w <- 1 / sy^2
  fit <- stats::lm(by ~ 0 + bx, weights = w)
  est <- stats::coef(fit); names(est) <- colnames(bx) %||% paste0("b", seq_len(ncol(bx)))
  sigma <- if (nrow(bx) > ncol(bx)) summary(fit)$sigma else 1
  scale <- switch(se_mode, fe = 1 / sigma, mre = 1 / min(1, sigma), presso = 1)
  V <- stats::vcov(fit) * scale^2
  dimnames(V) <- list(names(est), names(est))
  Q <- sum(w * (by - bx %*% est)^2)
  list(est = est, vcov = V, Q = Q, Q_df = nrow(bx) - ncol(bx), n_snp = nrow(bx))
}

# Univariable leg. `dat` is a TwoSampleMR-harmonised frame (beta.exposure, beta.outcome, ...).
# Method codes follow the triad config. MR-PRESSO is run once to find outliers; the result carries
# them so that slope-grid refits reuse the same outlier set (refitting PRESSO's 1000-draw outlier
# search at every grid point is not feasible and would make the curve jumpy).
UV_METHODS <- c("ivw_fe", "ivw_mre", "mr_egger", "weighted_median", "mr_raps", "mr_presso")

run_uv_leg <- function(dat, method, presso_outliers = NULL, presso_nboot = 1000) {
  method <- tolower(method)
  if (!method %in% UV_METHODS) stop("unknown univariable method '", method, "'; allowed: ",
                                    paste(UV_METHODS, collapse = ", "))
  bx <- dat$beta.exposure; by <- dat$beta.outcome; sy <- dat$se.outcome
  core <- ivw_core(bx, by, sy, "mre")
  base <- list(method = method, method_run = method, n_snp = nrow(dat), Q = core$Q,
               Q_p = stats::pchisq(core$Q, core$Q_df, lower.tail = FALSE),
               n_outlier = NA_integer_, presso_global_p = NA_real_, outliers = NULL)

  out <- switch(method,
    ivw_fe  = { f <- ivw_core(bx, by, sy, "fe");  list(b = f$est[[1]], se = sqrt(f$vcov[1, 1])) },
    ivw_mre = list(b = core$est[[1]], se = sqrt(core$vcov[1, 1])),
    mr_presso = {
      if (is.null(presso_outliers)) {
        pr <- if (nrow(dat) >= 4) tryCatch(suppressWarnings(MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
          SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
          data = as.data.frame(dat), NbDistribution = presso_nboot, SignifThreshold = 0.05)),
          error = function(e) NULL) else NULL
        if (is.null(pr)) {
          base$method_run <- "ivw_mre (MR-PRESSO failed or < 4 SNPs)"
          presso_outliers <- character()
        } else {
          idx <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
          presso_outliers <- if (is.numeric(idx)) dat$SNP[idx] else character()
          gp <- pr$`MR-PRESSO results`$`Global Test`$Pvalue
          base$presso_global_p <- suppressWarnings(as.numeric(sub("^<", "", gp)))
        }
      }
      keep <- !dat$SNP %in% presso_outliers
      f <- ivw_core(bx[keep], by[keep], sy[keep], "presso")
      base$n_outlier <- length(presso_outliers); base$outliers <- presso_outliers
      base$n_snp <- sum(keep)
      list(b = f$est[[1]], se = sqrt(f$vcov[1, 1]))
    },
    {
      tsm <- c(mr_egger = "mr_egger_regression", weighted_median = "mr_weighted_median",
               mr_raps = "mr_raps")[[method]]
      r <- tryCatch(suppressMessages(TwoSampleMR::mr(dat, method_list = tsm)),
                    error = function(e) NULL)
      if (is.null(r) || nrow(r) == 0 || is.na(r$b[1])) {
        base$method_run <- paste0("ivw_mre (", method, " failed)")
        list(b = core$est[[1]], se = sqrt(core$vcov[1, 1]))
      } else list(b = r$b[1], se = r$se[1])
    })
  c(base, out)
}

# Multivariable leg: columns of bx are (X, M); returns est = c(direct, my) and their 2x2 vcov.
MV_METHODS <- c("mv_ivw_mre", "mv_ivw_fe", "mv_presso", "grapple")

run_mv_leg <- function(mv, by, sy, method, presso_outliers = NULL, presso_nboot = 1000,
                       grapple_loss = "huber") {
  method <- tolower(method)
  if (!method %in% MV_METHODS) stop("unknown MVMR method '", method, "'; allowed: ",
                                    paste(MV_METHODS, collapse = ", "))
  BX <- cbind(direct = mv$bx, my = mv$bm)
  base <- list(method = method, method_run = method, n_snp = nrow(mv),
               n_outlier = NA_integer_, presso_global_p = NA_real_, outliers = NULL)
  fit <- switch(method,
    mv_ivw_mre = ivw_core(BX, by, sy, "mre"),
    mv_ivw_fe  = ivw_core(BX, by, sy, "fe"),
    mv_presso = {
      if (is.null(presso_outliers)) {
        df <- data.frame(by = by, sy = sy, bx = mv$bx, sx = mv$sx, bm = mv$bm, sm = mv$sm)
        pr <- tryCatch(suppressWarnings(MRPRESSO::mr_presso(
          BetaOutcome = "by", BetaExposure = c("bx", "bm"), SdOutcome = "sy",
          SdExposure = c("sx", "sm"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = df,
          NbDistribution = presso_nboot, SignifThreshold = 0.05)), error = function(e) NULL)
        if (is.null(pr)) { base$method_run <- "mv_ivw_mre (MV MR-PRESSO failed)"; idx <- NULL }
        else {
          idx <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
          gp <- pr$`MR-PRESSO results`$`Global Test`$Pvalue
          base$presso_global_p <- suppressWarnings(as.numeric(sub("^<", "", gp)))
        }
        presso_outliers <- if (is.numeric(idx)) mv$SNP[idx] else character()
      }
      keep <- !mv$SNP %in% presso_outliers
      base$n_outlier <- length(presso_outliers); base$outliers <- presso_outliers
      base$n_snp <- sum(keep)
      ivw_core(BX[keep, , drop = FALSE], by[keep], sy[keep], "presso")
    },
    grapple = {
      # GRAPPLE expects getInput()-style columns. Instruments are already selected upstream, so
      # p.thres = 1 keeps every row. Not exercised locally (GRAPPLE not installed when written):
      # check the first render's output against mv_ivw_mre before trusting it.
      if (!requireNamespace("GRAPPLE", quietly = TRUE)) {
        base$method_run <- "mv_ivw_mre (GRAPPLE not installed)"
        ivw_core(BX, by, sy, "mre")
      } else {
        gd <- data.frame(SNP = mv$SNP, gamma_exp1 = mv$bx, gamma_exp2 = mv$bm, gamma_out1 = by,
                         se_exp1 = mv$sx, se_exp2 = mv$sm, se_out1 = sy,
                         selection_pvals = pmin(mv$pval_x, mv$pval_m))
        g <- tryCatch(GRAPPLE::grappleRobustEst(gd, p.thres = 1, plot.it = FALSE,
                                                loss.function = grapple_loss),
                      error = function(e) NULL)
        if (is.null(g)) { base$method_run <- "mv_ivw_mre (GRAPPLE failed)"; ivw_core(BX, by, sy, "mre") }
        else {
          est <- setNames(as.numeric(g$beta.hat), colnames(BX))
          V <- as.matrix(g$beta.var); dimnames(V) <- list(names(est), names(est))
          list(est = est, vcov = V, Q = NA_real_, Q_df = NA_integer_, n_snp = nrow(mv))
        }
      }
    })
  c(base, list(est = fit$est, vcov = fit$vcov))
}

# ---- 5. Legs at a given selection slope, and Monte Carlo propagation -----------------------------

# Fit total / X->M / MVMR on the frames at one value of b_SH. X->M does not involve Y, so it is
# passed in pre-fit. `frozen` carries MR-PRESSO outlier sets from the point fit.
legs_at_slope <- function(fr, b_sh, methods, xm_fit, frozen = list(), presso_nboot = 1000,
                          grapple_loss = "huber") {
  ya <- adjust_outcome(fr$uv$by, fr$uv$sy, fr$uv$bl, fr$uv$sl, b_sh)
  uvdat <- data.frame(SNP = fr$uv$SNP, beta.exposure = fr$uv$bx, se.exposure = fr$uv$sx,
                      beta.outcome = ya$b, se.outcome = ya$s, mr_keep = TRUE,
                      id.exposure = "X", id.outcome = "Y", exposure = "X", outcome = "Y")
  tot <- run_uv_leg(uvdat, methods$total, presso_outliers = frozen$total,
                    presso_nboot = presso_nboot)
  ma <- adjust_outcome(fr$mv$by, fr$mv$sy, fr$mv$bl, fr$mv$sl, b_sh)
  mvf <- run_mv_leg(fr$mv, ma$b, ma$s, methods$mvmr, presso_outliers = frozen$mvmr,
                    presso_nboot = presso_nboot, grapple_loss = grapple_loss)
  list(total = tot, xm = xm_fit, mv = mvf)
}

# Flatten a legs object to the 4-vector of estimates and the 4x4 covariance used for sampling.
# Order: total, xm, direct, my. Off-diagonal terms between legs come from `R` (a correlation
# matrix estimated by snp_bootstrap_cor); the direct/my block is the MVMR method's own vcov.
LEG_NAMES <- c("total", "xm", "direct", "my")

legs_mu_sigma <- function(legs, R = diag(4)) {
  mu <- c(legs$total$b, legs$xm$b, legs$mv$est[["direct"]], legs$mv$est[["my"]])
  se <- c(legs$total$se, legs$xm$se, sqrt(diag(legs$mv$vcov)))
  S <- diag(se) %*% R %*% diag(se)
  S[3:4, 3:4] <- legs$mv$vcov
  dimnames(S) <- list(LEG_NAMES, LEG_NAMES)
  list(mu = setNames(mu, LEG_NAMES), Sigma = make_psd(S))
}

make_psd <- function(S, eps = 1e-12) {
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  if (min(e$values) >= eps) return(S)
  out <- e$vectors %*% diag(pmax(e$values, eps)) %*% t(e$vectors)
  dimnames(out) <- dimnames(S); out
}

# Correlation between the four leg estimates. Total, X->M and MVMR reuse the same SNPs and the same
# GWAS (X instruments appear in all three; Y in two), so their sampling errors are correlated and
# the ratio indirect/total cannot be sampled as if they were independent. Estimated by a parametric
# SNP-level bootstrap using fast IVW-type fits: one noise draw per SNP x GWAS, shared across tables
# (possible because every table is aligned to X's effect allele). Only the CORRELATION is taken
# from here; each leg keeps its own method's SE.
snp_bootstrap_cor <- function(fr, b_sh, n_boot = 300, seed = 1) {
  set.seed(seed)
  snps <- unique(c(fr$uv$SNP, fr$mv$SNP, fr$xm$SNP))
  idx <- function(s) match(s, snps)
  draws <- replicate(n_boot, {
    z <- matrix(stats::rnorm(length(snps) * 4), ncol = 4,
                dimnames = list(NULL, c("x", "m", "y", "life")))
    u <- fr$uv; i <- idx(u$SNP)
    ya <- adjust_outcome(u$by + u$sy * z[i, "y"], u$sy, u$bl + u$sl * z[i, "life"], u$sl, b_sh)
    tot <- ivw_core(u$bx + u$sx * z[i, "x"], ya$b, ya$s)$est[[1]]
    x <- fr$xm; i <- idx(x$SNP)
    xm <- ivw_core(x$beta.exposure + x$se.exposure * z[i, "x"],
                   x$beta.outcome + x$se.outcome * z[i, "m"], x$se.outcome)$est[[1]]
    v <- fr$mv; i <- idx(v$SNP)
    ma <- adjust_outcome(v$by + v$sy * z[i, "y"], v$sy, v$bl + v$sl * z[i, "life"], v$sl, b_sh)
    mv <- ivw_core(cbind(v$bx + v$sx * z[i, "x"], v$bm + v$sm * z[i, "m"]), ma$b, ma$s)$est
    c(tot, xm, mv)
  })
  R <- stats::cor(t(draws))
  dimnames(R) <- list(LEG_NAMES, LEG_NAMES)
  R
}

# Default leg sampler: frequentist sampling distributions, MVN(mu, Sigma).
# THIS IS THE BAYESIAN HOOK. Any function(mu, Sigma, n) returning an n x 4 matrix with columns
# LEG_NAMES can replace it: e.g. posterior draws of each leg from brms fits of the SNP-level
# regressions (brm(by | se(sy) ~ 0 + bx)), or MVN draws from a normal-normal posterior with
# informative priors on the legs. Everything downstream (indirect, proportion, gates) is computed
# per draw and does not care where the draws came from.
sample_legs_mvn <- function(mu, Sigma, n) {
  out <- MASS::mvrnorm(n, mu, Sigma)
  if (n == 1) out <- matrix(out, nrow = 1)
  colnames(out) <- LEG_NAMES
  out
}

# Draws of the selection slope. Uses saved SlopeHunter bootstrap replicates when they exist
# (the fits CSV only stores b_SH and its bootstrap SE, so by default this falls back to
# N(b_SH, se_SH^2) and says so).
draw_slope <- function(n, b_sh, se_sh, boot = NULL) {
  if (!is.null(boot) && length(boot) >= 50)
    return(list(draws = sample(boot, n, replace = TRUE), source = "SlopeHunter bootstrap replicates"))
  list(draws = stats::rnorm(n, b_sh, se_sh), source = "normal approximation N(b_SH, se_SH^2)")
}

# Full Monte Carlo for one triad x outcome x adjustment arm.
# Unadjusted: legs fit once, draws = MVN(mu, Sigma).
# Adjusted:   legs are refit on a grid of b_SH values (+/- grid_sd SE). For each draw, b_SH_d is
#             sampled, the leg estimates, SEs and MVMR covariance are interpolated at b_SH_d (the
#             "re-adjust the outcome and re-run the legs" step, done once per grid point rather
#             than once per draw), and the legs are then sampled around those values.
run_mc <- function(fr, methods, xm_fit, b_sh = 0, se_sh = 0, slope_boot = NULL,
                   n_draws = 5000, n_grid = 21, grid_sd = 4, n_boot_cor = 300, seed = 1,
                   sampler = sample_legs_mvn, presso_nboot = 1000, grapple_loss = "huber") {
  set.seed(seed)
  point <- legs_at_slope(fr, b_sh, methods, xm_fit, presso_nboot = presso_nboot,
                         grapple_loss = grapple_loss)
  frozen <- list(total = point$total$outliers, mvmr = point$mv$outliers)
  R <- snp_bootstrap_cor(fr, b_sh, n_boot = n_boot_cor, seed = seed)
  ps <- legs_mu_sigma(point, R)

  adjusted <- b_sh != 0 && is.finite(se_sh) && se_sh > 0
  if (!adjusted) {
    d <- sampler(ps$mu, ps$Sigma, n_draws)
    slope <- rep(b_sh, n_draws); slope_src <- "fixed (no adjustment)"; n_clip <- 0L
  } else {
    grid <- b_sh + se_sh * seq(-grid_sd, grid_sd, length.out = n_grid)
    gfit <- map(grid, function(g) {
      lp <- legs_mu_sigma(legs_at_slope(fr, g, methods, xm_fit, frozen = frozen,
                                        presso_nboot = presso_nboot,
                                        grapple_loss = grapple_loss), R)
      c(lp$mu, Sigma = as.vector(lp$Sigma))
    }) |> do.call(what = rbind)
    sl <- draw_slope(n_draws, b_sh, se_sh, slope_boot)
    slope <- sl$draws; slope_src <- sl$source
    n_clip <- sum(slope < min(grid) | slope > max(grid))
    slope_c <- pmin(pmax(slope, min(grid)), max(grid))
    interp <- apply(gfit, 2, function(col) stats::approx(grid, col, xout = slope_c)$y)
    d <- t(vapply(seq_len(n_draws), function(k) {
      S <- make_psd(matrix(interp[k, -(1:4)], 4, 4))
      sampler(interp[k, 1:4], S, 1)[1, ]
    }, numeric(4)))
    colnames(d) <- LEG_NAMES
  }
  draws <- as_tibble(d) |>
    mutate(draw = row_number(), b_SH = slope, indirect = xm * my, proportion = indirect / total,
           difference = total - direct, .before = 1)
  list(point = point, R = R, draws = draws, slope_source = slope_src, n_slope_clipped = n_clip)
}

# ---- 6. Diagnostics -------------------------------------------------------------------------------

# Sanderson-Windmeijer conditional F (two-sample form) and the MVMR Q statistic, via the MVMR
# package. gencov = 0 assumes X and M come from non-overlapping samples; if they overlap, pass a
# covariance from MVMR::phenocov_mvmr() instead.
mvmr_strength <- function(mv, by, sy, gencov = 0) {
  fd <- MVMR::format_mvmr(BXGs = cbind(mv$bx, mv$bm), BYG = by,
                          seBXGs = cbind(mv$sx, mv$sm), seBYG = sy, RSID = mv$SNP)
  # Both functions print to stdout; capture it so the rendered page stays clean.
  # gencov = 0 is deliberate, so its "covariance not specified" warning is silenced too.
  quiet <- function(expr) { utils::capture.output(r <- suppressWarnings(suppressMessages(expr))); r }
  sw <- tryCatch(quiet(MVMR::strength_mvmr(fd, gencov)),   error = function(e) NULL)
  pq <- tryCatch(quiet(MVMR::pleiotropy_mvmr(fd, gencov)), error = function(e) NULL)
  tibble(cond_F_x = if (is.null(sw)) NA_real_ else as.numeric(sw[1, 1]),
         cond_F_m = if (is.null(sw)) NA_real_ else as.numeric(sw[1, 2]),
         mvmr_Q   = if (is.null(pq)) NA_real_ else as.numeric(pq$Qstat),
         mvmr_Q_p = if (is.null(pq)) NA_real_ else as.numeric(pq$Qpval))
}

# Steiger directionality on the X -> M leg: SNPs that explain more variance in M than in X are
# candidates for reverse causation (glucose/insulin -> lean mass) and are reported, plus a
# sensitivity X -> M estimate with them removed.
steiger_xm <- function(xm, method) {
  st <- tryCatch(suppressMessages(TwoSampleMR::steiger_filtering(xm)), error = function(e) NULL)
  if (is.null(st) || !"steiger_dir" %in% names(st))
    return(tibble(steiger_n = NA_integer_, steiger_n_wrong = NA_integer_,
                  xm_b_steiger = NA_real_, xm_se_steiger = NA_real_, steiger_wrong_snps = NA_character_))
  wrong <- st$SNP[!st$steiger_dir & st$steiger_pval < 0.05]
  keep <- st |> filter(!SNP %in% wrong)
  f <- if (nrow(keep) >= 3) run_uv_leg(keep, method) else list(b = NA_real_, se = NA_real_)
  tibble(steiger_n = nrow(st), steiger_n_wrong = length(wrong),
         xm_b_steiger = f$b, xm_se_steiger = f$se,
         steiger_wrong_snps = if (length(wrong)) paste(wrong, collapse = ";") else NA_character_)
}

# How much do the two instrument sets share? X instruments' associations with M, and vice versa.
cross_lookup <- function(xm, mx, p_strict = 5e-8, p_nominal = 1e-5) {
  f <- function(d, lab) tibble(
    direction = lab, n_looked_up = nrow(d),
    n_gw_sig = sum(d$pval.outcome < p_strict, na.rm = TRUE),
    n_suggestive = sum(d$pval.outcome < p_nominal, na.rm = TRUE),
    n_nominal = sum(d$pval.outcome < 0.05, na.rm = TRUE),
    expected_nominal = round(0.05 * nrow(d), 1))
  bind_rows(f(xm, "X instruments in M GWAS"), f(mx, "M instruments in X GWAS"))
}

# ---- 7. Summaries and gates ------------------------------------------------------------------------

mc_summary <- function(x, prefix) {
  x <- x[is.finite(x)]
  p <- if (length(x)) max(2 * min(mean(x > 0), mean(x < 0)), 1 / length(x)) else NA_real_
  tibble(!!paste0(prefix, "_mc_median") := stats::median(x),
         !!paste0(prefix, "_mc_lo") := unname(stats::quantile(x, 0.025)),
         !!paste0(prefix, "_mc_hi") := unname(stats::quantile(x, 0.975)),
         !!paste0(prefix, "_mc_p") := p)
}

# Pre-specified gates. A triad that fails is still reported; it just gets no percent mediated.
apply_gates <- function(row, f_min = 10, equiv_bound = NA_real_, equiv_draws = NULL) {
  excl0 <- function(lo, hi) !is.na(lo) && !is.na(hi) && (lo > 0 || hi < 0)
  g_F     <- isTRUE(row$cond_F_x > f_min) && isTRUE(row$cond_F_m > f_min)
  g_total <- excl0(row$total_mc_lo, row$total_mc_hi)
  g_sign  <- isTRUE(sign(row$indirect_b) == sign(row$total_b))
  flag_direct_opp <- isTRUE(sign(row$direct_b) != sign(row$total_b))
  prop_ok <- isTRUE(row$proportion_b >= 0 && row$proportion_b <= 1)
  # Equivalence (TOST at alpha = 0.05): indirect's 90% interval inside (-bound, +bound).
  equiv <- if (is.na(equiv_bound)) "bound not set"
           else if (is.null(equiv_draws)) "no draws"
           else {
             q <- stats::quantile(equiv_draws, c(0.05, 0.95), na.rm = TRUE)
             if (q[1] > -equiv_bound && q[2] < equiv_bound) "indirect equivalent to 0"
             else if (excl0(row$indirect_mc_lo, row$indirect_mc_hi)) "indirect non-zero"
             else "inconclusive"
           }
  report_prop <- g_F && g_total && g_sign && !flag_direct_opp && prop_ok
  tibble(gate_condF = g_F, gate_total_ci = g_total, gate_sign_concordant = g_sign,
         flag_direct_opposite = flag_direct_opp, flag_prop_outside_01 = !prop_ok,
         equivalence = equiv, report_proportion = report_prop,
         gate_summary = if (report_prop) "PASS" else paste(c(
           if (!g_F) "conditional F <= 10", if (!g_total) "total CI includes 0",
           if (!g_sign) "indirect/total discordant", if (flag_direct_opp) "direct opposite sign",
           if (!prop_ok) "proportion outside [0,1]"), collapse = "; "))
}

# ---- 8. Orchestration ------------------------------------------------------------------------------

# One triad x outcome x adjustment arm -> one result row, its Monte Carlo draws, and its log.
#   tr       one row of the triad config (triad_id, outcome_id, method_*, mediator_type, ...)
#   pulls    list(x_id, m_id, x_inst, m_inst, mv_exp, x_in_m, m_in_x) for this triad
#   ot       the outcome table for tr$outcome_id (build_outcome_table()$table)
#   bsh      one row: b_SH, se_SH, trusted, reason, and optionally a list-column `boot`
#   st       settings list: n_draws, n_grid, grid_sd, n_boot_cor, seed, cond_f_min,
#            k_min_x_in_mediator, common_snp_set, presso_nboot, grapple_loss, sampler
run_triad_arm <- function(tr, pulls, ot, adjustment, bsh, par, st) {
  ids <- tibble(triad_id = tr$triad_id, outcome_id = tr$outcome_id, adjustment = adjustment)
  stub <- function(status, log = empty_log())
    list(result = ids |> mutate(status = status), draws = NULL, log = log)

  adjusted <- adjustment == "slopehunter"
  if (adjusted && !isTRUE(bsh$trusted))
    return(stub(paste0("NOT RUN: b_SH not applied (", bsh$reason %||% "no fit", ")")))

  n_x_in_m <- n_distinct(pulls$x_in_m$SNP %||% character())
  if (identical(tr$mediator_type, "eqtl") && n_x_in_m < st$k_min_x_in_mediator)
    return(stub(paste0("SKIPPED: ", n_x_in_m, " X instruments present in the mediator dataset (< ",
                       st$k_min_x_in_mediator, "); most eQTL resources are cis-only")))

  fr <- build_triad_frames(pulls, ot, par, adjusted, st$common_snp_set, tr$triad_id, tr$outcome_id)
  if (min(nrow(fr$uv), nrow(fr$xm), nrow(fr$mv)) < 3)
    return(stub(sprintf("NOT RUN: too few SNPs (total %d, X->M %d, MVMR %d)",
                        nrow(fr$uv), nrow(fr$xm), nrow(fr$mv)), fr$log))

  methods <- list(total = tr$method_XY_total, mvmr = tr$method_MVMR)
  xm_fit <- run_uv_leg(fr$xm, tr$method_XM, presso_nboot = st$presso_nboot)
  b_sh  <- if (adjusted) bsh$b_SH else 0
  se_sh <- if (adjusted) bsh$se_SH else 0
  boot  <- if (adjusted && "boot" %in% names(bsh)) bsh$boot[[1]] else NULL

  mc <- run_mc(fr, methods, xm_fit, b_sh, se_sh, slope_boot = boot, n_draws = st$n_draws,
               n_grid = st$n_grid, grid_sd = st$grid_sd, n_boot_cor = st$n_boot_cor,
               seed = st$seed, sampler = st$sampler %||% sample_legs_mvn,
               presso_nboot = st$presso_nboot, grapple_loss = st$grapple_loss)
  P <- mc$point
  ma <- adjust_outcome(fr$mv$by, fr$mv$sy, fr$mv$bl, fr$mv$sl, b_sh)

  tot_b <- P$total$b; xm_b <- P$xm$b; dir_b <- P$mv$est[["direct"]]; my_b <- P$mv$est[["my"]]
  se <- c(total = P$total$se, xm = P$xm$se, sqrt(diag(P$mv$vcov)))
  ind_b <- xm_b * my_b

  row <- ids |> mutate(
    status = "run",
    n_snp_total = P$total$n_snp, n_snp_xm = P$xm$n_snp, n_snp_mvmr = P$mv$n_snp,
    n_mvmr_x_inst = sum(fr$mv$in_x_inst), n_mvmr_m_inst = sum(fr$mv$in_m_inst),
    method_total = P$total$method_run, method_xm = P$xm$method_run, method_mvmr = P$mv$method_run,
    total_b = tot_b,  total_se = se[["total"]],
    xm_b = xm_b,      xm_se = se[["xm"]],
    direct_b = dir_b, direct_se = se[["direct"]],
    my_b = my_b,      my_se = se[["my"]],
    cov_direct_my = P$mv$vcov["direct", "my"],
    indirect_b = ind_b,
    # Delta method, quick check only: ignores b_SH uncertainty and the leg correlations.
    indirect_se_delta = sqrt(xm_b^2 * se[["my"]]^2 + my_b^2 * se[["xm"]]^2),
    proportion_b = ind_b / tot_b,
    difference_b = tot_b - dir_b,
    # Consistency: direct + indirect should reproduce total. They are estimated on different SNP
    # sets with different methods, so exact agreement is not expected.
    consistency_gap = (dir_b + ind_b) - tot_b,
    consistency_gap_in_total_se = ((dir_b + ind_b) - tot_b) / se[["total"]],
    total_Q_p = P$total$Q_p, xm_Q_p = P$xm$Q_p,
    total_n_outlier = P$total$n_outlier, mvmr_n_outlier = P$mv$n_outlier,
    b_SH = b_sh, se_SH = se_sh, slope_source = mc$slope_source,
    n_slope_clipped = mc$n_slope_clipped,
    r_total_direct = mc$R["total", "direct"]
  )
  mcs <- map2(c("total", "xm", "direct", "my", "indirect", "proportion", "difference"),
              c("total", "xm", "direct", "my", "indirect", "proportion", "difference"),
              \(col, pre) mc_summary(mc$draws[[col]], pre)) |> bind_cols()
  row <- bind_cols(row, mcs, mvmr_strength(fr$mv, ma$b, ma$s, gencov = st$gencov %||% 0))
  gates <- apply_gates(row, st$cond_f_min, tr$equiv_bound_indirect %||% NA_real_,
                       mc$draws$indirect)
  row <- bind_cols(row, gates) |>
    mutate(across(starts_with("proportion_"), \(x) if_else(report_proportion, x, NA_real_)))

  list(result = row,
       draws = bind_cols(ids, mc$draws),
       log = bind_rows(fr$log,
         if (mc$n_slope_clipped > 0)
           log_row(tr$triad_id, tr$outcome_id, "monte_carlo", "note",
                   paste(mc$n_slope_clipped, "b_SH draws fell outside the refit grid and were",
                         "clamped to its edge")) else NULL))
}

# Per-triad diagnostics that do not depend on the outcome: the X -> M assumption checks that the
# prior univariable work did not cover, Steiger on X -> M, the reverse M -> X estimate, and the
# instrument cross-lookup.
triad_diagnostics <- function(tr, pulls, par, presso_nboot = 1000) {
  xm <- suppressMessages(TwoSampleMR::harmonise_data(pulls$x_inst, pulls$x_in_m,
                                                     action = par$harmonise_action)) |>
    filter(mr_keep)
  mx <- if (is.null(pulls$m_in_x)) NULL else
    suppressMessages(TwoSampleMR::harmonise_data(pulls$m_inst, pulls$m_in_x,
                                                 action = par$harmonise_action)) |>
      filter(mr_keep)
  egg <- tryCatch(TwoSampleMR::mr_pleiotropy_test(xm), error = function(e) NULL)
  pr  <- run_uv_leg(xm, "mr_presso", presso_nboot = presso_nboot)
  ivw <- run_uv_leg(xm, "ivw_mre")
  rev <- if (!is.null(mx) && nrow(mx) >= 3) run_uv_leg(mx, "ivw_mre") else list(b = NA, se = NA, n_snp = nrow(mx %||% data.frame()))
  checks <- tibble(
    triad_id = tr$triad_id,
    xm_n_snp = nrow(xm), xm_mean_F = mean((xm$beta.exposure / xm$se.exposure)^2),
    xm_ivw_mre_b = ivw$b, xm_ivw_mre_se = ivw$se, xm_Q_p = ivw$Q_p,
    xm_egger_intercept = egg$egger_intercept %||% NA_real_, xm_egger_int_p = egg$pval %||% NA_real_,
    xm_presso_global_p = pr$presso_global_p, xm_presso_n_outlier = pr$n_outlier,
    xm_presso_b = pr$b,
    reverse_mx_n_snp = rev$n_snp, reverse_mx_b = rev$b, reverse_mx_se = rev$se,
    reverse_mx_p = 2 * stats::pnorm(-abs(rev$b / rev$se))
  )
  list(checks = bind_cols(checks, steiger_xm(xm, tr$method_XM)),
       cross = cross_lookup(xm, mx %||% xm[0, ]) |> mutate(triad_id = tr$triad_id, .before = 1),
       xm_dat = xm)
}

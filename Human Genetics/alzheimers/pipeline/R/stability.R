# stability.R — pre-specified stability flags for the category-2 comparison (MRBEE, MR-RAPS).
#
# Rule, fixed before results were examined: MRBEE and MR-RAPS are ALWAYS reported; a row is flagged
# "unstable" only on process grounds, never on the size or direction of its estimate:
#   - the fit errored or returned no estimate
#   - the SE is not finite or not positive
#   - MR-RAPS: the over-dispersion optimiser FAILED ("Did not find a solution for tau2" /
#     "Cannot find solution with finite over.dispersion"). NOT flagged: "Estimated overdispersion
#     is negative. Using tau2 = 0", which is RAPS correctly finding no heterogeneity (e.g. 2hGlu ->
#     Kunkle, Q p = 0.60); that warning is kept in stability_note but the fit counts as stable.
# The two methods are refit on the classical target's stored instruments (seconds), so the flags
# need no rerun of the full classical suite; the estimates reported remain the suite's own.

suppressPackageStartupMessages({ library(dplyr) })

capture_fit <- function(expr) {
  warns <- character()
  val <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                             warning = function(w) { warns <<- c(warns, conditionMessage(w))
                                                     invokeRestart("muffleWarning") })
  list(value = val, warnings = warns)
}

stability_flags <- function(cl) {
  ids <- cl$rows |> distinct(assoc_id)
  if (is.null(cl$instruments) || nrow(cl$instruments) < 3)
    return(tibble(assoc_id = ids$assoc_id, method = c("MRBEE", "MR-RAPS"),
                  stable = FALSE, stability_note = "not run: < 3 instruments"))
  i <- cl$instruments
  dat <- data.frame(SNP = i$SNP, beta.exposure = i$beta_exposure, se.exposure = i$se_exposure,
                    beta.outcome = i$beta_outcome, se.outcome = i$se_outcome, mr_keep = TRUE,
                    id.exposure = "x", id.outcome = "y", exposure = "x", outcome = "y")
  r12 <- cl$rows$mrbee_error_cor[cl$rows$method == "MRBEE"][1]

  judge <- function(method, fit, b_se, extra_bad = character()) {
    v <- fit$value
    reasons <- c(
      if (inherits(v, "error")) paste("error:", conditionMessage(v)),
      if (!inherits(v, "error") && (!is.finite(b_se[1]) || !is.finite(b_se[2]) || b_se[2] <= 0))
        "no finite estimate/SE",
      extra_bad)
    tibble(assoc_id = ids$assoc_id, method = method, stable = !length(reasons),
           stability_note = if (length(reasons)) paste(unique(reasons), collapse = "; ") else
             if (length(fit$warnings)) paste("warnings:", paste(unique(fit$warnings), collapse = "; ")) else "")
  }

  raps <- capture_fit(suppressMessages(TwoSampleMR::mr(dat, method_list = "mr_raps")))
  rb <- if (!inherits(raps$value, "error") && nrow(raps$value)) c(raps$value$b[1], raps$value$se[1]) else c(NA, NA)
  tau0 <- grep("Did not find a solution for tau2|Cannot find solution with finite over.dispersion",
               raps$warnings, value = TRUE)

  mb <- if (is.na(r12)) list(value = simpleError("no error-correlation estimate"), warnings = character()) else
    capture_fit(MRBEE::MRBEE.IMRP.UV(by = dat$beta.outcome, bx = dat$beta.exposure,
                                     byse = dat$se.outcome, bxse = dat$se.exposure,
                                     Rxy = matrix(c(1, r12, r12, 1), 2)))
  mbb <- if (!inherits(mb$value, "error")) c(as.numeric(mb$value$theta)[1], sqrt(as.numeric(mb$value$vartheta)[1])) else c(NA, NA)

  bind_rows(judge("MR-RAPS", raps, rb, if (length(tau0)) "over-dispersion optimiser failed (tau2 forced to 0)"),
            judge("MRBEE", mb, mbb))
}

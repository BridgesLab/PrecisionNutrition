# plots.R — per-association diagnostic figures, in the style of the calcium arm
# (calcium-cholesterol/mr-calcium-tc.qmd: scatter, "volcano" funnel, leave-one-out).
#
# Built from the classical and MR-APSS targets only (never from `results`), so these targets do not
# rebuild every time an unrelated association changes. Each returns the paths it wrote, for a
# format = "file" target.

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

FIG_COLOURS <- c("IVW-MRE" = "#00274c", "MR-Egger" = "#d86018", "Weighted median" = "#75988d",
                 "MR-PRESSO" = "#a5a508", "MRBEE" = "#702082", "MR-RAPS" = "#ffcb05",
                 "MR-APSS" = "#ffcb05", "CAUSE" = "#d86018")

# Instruments oriented to the exposure-raising allele, as TwoSampleMR and MR plots conventionally do.
oriented <- function(i) {
  s <- sign(i$beta_exposure); s[s == 0] <- 1
  tibble(SNP = i$SNP, bx = i$beta_exposure * s, sx = i$se_exposure,
         by = i$beta_outcome * s, sy = i$se_outcome)
}

pretty_id <- function(assoc_id) gsub("__", "  →  ", sub("__raw$", "", sub("__slopehunter$", "  [SlopeHunter]", assoc_id)))

# Scatter of SNP effects with 95% CIs on both axes and one line per method. Egger is drawn with its
# intercept; the others pass through the origin.
plot_scatter <- function(d, fits, xlab, ylab, title) {
  ggplot(d, aes(bx, by)) +
    geom_hline(yintercept = 0, colour = "grey75") + geom_vline(xintercept = 0, colour = "grey75") +
    geom_errorbar(aes(ymin = by - 1.96 * sy, ymax = by + 1.96 * sy), alpha = 0.25, width = 0) +
    geom_errorbar(aes(xmin = bx - 1.96 * sx, xmax = bx + 1.96 * sx), alpha = 0.25, width = 0,
                  orientation = "y") +
    geom_point(size = 1) +
    geom_abline(data = fits, aes(intercept = intercept, slope = b, colour = method), linewidth = 0.9) +
    scale_colour_manual(values = FIG_COLOURS) +
    guides(colour = guide_legend(nrow = 2)) +
    labs(x = xlab, y = ylab, colour = "MR method", title = title) +
    theme_classic(base_size = 14) + theme(legend.position = "bottom")
}

# "Volcano" funnel: single-SNP Wald ratios against their precision, the IVW-MRE estimate, and the
# pseudo-95% cone IVW +/- 1.96/precision. Asymmetry suggests directional pleiotropy. The x-axis is
# limited to the central 99% of ratios so one wild SNP cannot flatten the plot; the caption counts
# how many fall outside.
plot_funnel <- function(d, b_ivw, title) {
  w <- d |> mutate(b = by / bx, se = sy / abs(bx), precision = 1 / se)
  lims <- stats::quantile(w$b, c(0.005, 0.995), na.rm = TRUE)
  lims <- range(c(lims, b_ivw))
  pmax <- max(w$precision) * 1.1
  cone <- tibble(precision = seq(pmax / 1000, pmax, length.out = 1000)) |>
    mutate(lower = b_ivw - 1.96 / precision, upper = b_ivw + 1.96 / precision)
  n_out <- sum(w$b < lims[1] | w$b > lims[2])
  ggplot(w, aes(b, precision)) +
    geom_point(size = 1) +
    geom_vline(xintercept = b_ivw, colour = "#ff7f0e", linewidth = 1) +
    geom_line(data = cone, aes(lower, precision), linetype = "dashed") +
    geom_line(data = cone, aes(upper, precision), linetype = "dashed") +
    coord_cartesian(xlim = lims, ylim = c(0, pmax)) +
    labs(x = "Single-SNP estimate (Wald ratio)", y = "Precision (1/SE)", title = title,
         caption = if (n_out) sprintf("%d SNP(s) beyond the x-axis limits", n_out) else NULL) +
    theme_classic(base_size = 14)
}

# Leave-one-out IVW-MRE (residual SE floored at 1, as in the pipeline's IVW-MRE). Points sorted by
# estimate; the band is the all-SNP estimate's 95% CI; the 5 SNPs that move it most are labelled.
loo_ivw <- function(d) {
  fit1 <- function(k) {
    s <- summary(stats::lm(d$by[k] ~ -1 + d$bx[k], weights = 1 / d$sy[k]^2))
    c(b = s$coefficients[1, 1], se = s$coefficients[1, 2] / s$sigma * max(1, s$sigma))
  }
  all <- fit1(seq_len(nrow(d)))
  out <- t(vapply(seq_len(nrow(d)), \(j) fit1(-j), numeric(2)))
  tibble(SNP = d$SNP, b_loo = out[, "b"], se_loo = out[, "se"], b_all = all[["b"]],
         se_all = all[["se"]], shift = out[, "b"] - all[["b"]]) |> arrange(desc(abs(shift)))
}

plot_loo <- function(loo, title) {
  d <- loo |> arrange(b_loo) |> mutate(rank = row_number(), top = rank(-abs(shift)) <= 5)
  ggplot(d, aes(rank, b_loo)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = d$b_all[1] - 1.96 * d$se_all[1],
             ymax = d$b_all[1] + 1.96 * d$se_all[1], alpha = 0.12, fill = "#00274c") +
    geom_hline(yintercept = d$b_all[1], colour = "#00274c") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 1) +
    geom_text(data = d |> filter(top), aes(label = SNP), size = 3, vjust = -0.8, check_overlap = TRUE) +
    labs(x = "SNP left out (sorted by resulting estimate)", y = "IVW-MRE without that SNP",
         title = title, caption = "Band: all-SNP IVW-MRE 95% CI. Labels: the 5 most influential SNPs.") +
    theme_classic(base_size = 14) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

save_fig <- function(p, path, w = 7, h = 5.5) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggsave(path, p, width = w, height = h, dpi = 200)
  path
}

# Classical diagnostics for one association: scatter, funnel, LOO plot and LOO table.
diagnostics_classical <- function(cl, out_dir) {
  # Attached here, not only at the top of this file: these targets run on crew workers, which load
  # just the packages in tar_option_set(packages = ...), not this file's library() calls.
  suppressPackageStartupMessages(library(ggplot2))
  if (is.null(cl$instruments) || nrow(cl$instruments) < 3) return(character())
  id <- cl$rows$assoc_id[1]; dir <- file.path(out_dir, id); title <- pretty_id(id)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  d <- oriented(cl$instruments)
  r <- cl$rows
  fits <- r |>
    filter(method %in% c("IVW-MRE", "MR-Egger", "Weighted median", "MR-PRESSO", "MRBEE"), is.finite(b)) |>
    transmute(method, b, intercept = if_else(method == "MR-Egger", coalesce(egger_intercept, 0), 0))
  b_ivw <- r$b[r$method == "IVW-MRE"]
  loo <- loo_ivw(d)
  readr::write_csv(loo, file.path(dir, "leave_one_out.csv"))
  c(save_fig(plot_scatter(d, fits, "SNP effect on exposure", "SNP effect on outcome", title),
             file.path(dir, "scatter_classical.png")),
    save_fig(plot_funnel(d, b_ivw, title), file.path(dir, "funnel.png")),
    save_fig(plot_loo(loo, title), file.path(dir, "leave_one_out.png")),
    file.path(dir, "leave_one_out.csv"))
}

# MR-APSS scatter on its own (per-SD) scale and instruments, with the MR-APSS slope.
diagnostics_apss <- function(ap, out_dir) {
  suppressPackageStartupMessages(library(ggplot2))
  if (is.null(ap$instruments) || nrow(ap$instruments) < 3) return(character())
  id <- ap$rows$assoc_id[1]
  d <- oriented(ap$instruments)
  fits <- ap$rows |> transmute(method = "MR-APSS", b, intercept = 0)
  save_fig(plot_scatter(d, fits, "SNP effect on exposure (per SD)", "SNP effect on outcome (per SD)",
                        paste0(pretty_id(id), "  (MR-APSS instruments, p < 5e-5)")),
           file.path(out_dir, id, "scatter_mrapss.png"))
}

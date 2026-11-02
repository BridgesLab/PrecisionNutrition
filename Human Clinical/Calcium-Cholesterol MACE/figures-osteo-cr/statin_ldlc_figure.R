library(ggplot2)
library(dplyr)
library(patchwork)

# ── Panel A: Primary model — all covariates ─────────────────────────────────
# Osteoporosis Fine-Gray 3-state, time-averaged LDL-C, fully adjusted
# From ANALYSIS_NOTES.md (run 2026-04-06)

primary <- tibble::tribble(
  ~term,                        ~subhr,  ~lo,    ~hi,    ~highlight,
  "Time-avg LDL-C (per 10 mg/dL)", 0.985,  0.962,  1.009,  "LDL-C",
  "Ever statin use",            0.646,  0.530,  0.787,  "Statin",
  "Diabetes",                   0.782,  0.650,  0.941,  "Other",
  "Hypertension",               0.992,  0.762,  1.292,  "Other",
  "Black (vs White)",           0.662,  0.461,  0.951,  "Other",
  "Hispanic (vs White)",        0.584,  0.317,  1.076,  "Other",
  "Asian (vs White)",           0.873,  0.619,  1.232,  "Other",
  "Female (vs Male)",           6.035,  4.941,  7.372,  "Other",
  "Age at baseline (per yr)",   1.055,  1.048,  1.061,  "Other"
) |>
  mutate(
    term = factor(term, levels = rev(term)),
    p_label = case_when(
      highlight == "LDL-C" ~ "p = 0.21",
      highlight == "Statin" ~ "p < 0.001",
      term == "Age at baseline (per yr)" ~ "p < 0.001",
      term == "Female (vs Male)" ~ "p < 0.001",
      term == "Black (vs White)" ~ "p = 0.026",
      term == "Diabetes" ~ "p = 0.009",
      TRUE ~ ""
    )
  )

col_ldlc  <- "#00274C"   # Michigan Blue
col_statin <- "#CC4400"  # Burnt orange
col_other  <- "#888888"

pA <- ggplot(primary, aes(x = subhr, xmin = lo, xmax = hi, y = term,
                           color = highlight)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_errorbarh(height = 0.25, linewidth = 0.8) +
  geom_point(size = 3.2, shape = 15) +
  geom_text(aes(label = sprintf("%.3f", subhr)), hjust = -0.35, size = 3,
            color = "black", fontface = "plain") +
  scale_color_manual(values = c("LDL-C" = col_ldlc, "Statin" = col_statin,
                                 "Other" = col_other),
                     guide = "none") +
  scale_x_log10(breaks = c(0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 5.0, 10.0),
                labels = c("0.3", "0.5", "0.7", "1.0", "1.5", "2.0", "5.0", "10.0"),
                limits = c(0.22, 16)) +
  labs(
    title = "A   Primary model: osteoporosis",
    subtitle = "Fine-Gray 3-state (MACE + death competing)  |  n = 5,784  |  595 events",
    x = "Subdistribution hazard ratio (95% CI)",
    y = NULL
  ) +
  # Highlight boxes
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = 8.5, ymax = 9.5, fill = "#00274C", alpha = 0.08) +  # LDL-C row
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = 7.5, ymax = 8.5, fill = "#CC4400", alpha = 0.08) +  # Statin row
  annotate("label", x = 12, y = 9, label = "NULL\n(p = 0.21)",
           color = col_ldlc, size = 3, fontface = "bold", label.size = 0.4,
           fill = "#EEF4FF") +
  annotate("label", x = 12, y = 8, label = "PROTECTIVE\n(p < 0.001)",
           color = col_statin, size = 3, fontface = "bold", label.size = 0.4,
           fill = "#FFF0EE") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.y = element_text(color = c(rep("gray40", 7),
                                          col_statin, col_ldlc),
                               face = c(rep("plain", 7), "bold", "bold")),
    panel.grid.major.x = element_line(color = "gray92")
  )


# ── Panel B: LDL-C Sub-HR across all model specifications ──────────────────
models <- tibble::tribble(
  ~model,                             ~subhr,  ~lo,    ~hi,    ~type,
  "Fine-Gray: time-avg (PRIMARY)",    0.985,  0.962,  1.009,  "Fine-Gray",
  "Fine-Gray: baseline LDL-C",        0.989,  0.969,  1.010,  "Fine-Gray",
  "Fine-Gray: step/LOCF",             0.987,  0.964,  1.011,  "Fine-Gray",
  "Fine-Gray: statin-naive",          0.999,  0.954,  1.046,  "Fine-Gray",
  "Fine-Gray: all patients",          0.992,  0.980,  1.004,  "Fine-Gray",
  "Fine-Gray 2-state: MACE competing",0.983,  0.960,  1.007,  "Fine-Gray",
  "Cox: time-avg LDL-C",             0.981,  0.959,  1.004,  "Cox",
  "Cox: baseline LDL-C",             0.987,  0.967,  1.007,  "Cox",
  "Cox: time-varying LOCF",          0.984,  0.957,  1.011,  "Cox"
) |>
  mutate(model = factor(model, levels = rev(model)),
         is_primary = model == "Fine-Gray: time-avg (PRIMARY)")

pB <- ggplot(models, aes(x = subhr, xmin = lo, xmax = hi, y = model,
                          color = type, shape = is_primary, size = is_primary)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_errorbarh(height = 0.25, linewidth = 0.8, aes(size = NULL)) +
  geom_point() +
  scale_color_manual(values = c("Fine-Gray" = "#1A5C96", "Cox" = "#5A5A8A"),
                     name = "Framework") +
  scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 15), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2.8), guide = "none") +
  scale_x_continuous(limits = c(0.88, 1.12),
                     breaks = c(0.90, 0.95, 1.00, 1.05, 1.10)) +
  annotate("rect", xmin = 0.975, xmax = 1.025, ymin = -Inf, ymax = Inf,
           fill = "gray90", alpha = 0.4) +
  annotate("text", x = 0.893, y = 1.3, label = "← lower\nosteoporosis\nrisk",
           size = 2.8, color = "gray50", hjust = 0, lineheight = 0.9) +
  annotate("text", x = 1.107, y = 1.3, label = "higher\nosteoporosis\nrisk →",
           size = 2.8, color = "gray50", hjust = 1, lineheight = 0.9) +
  labs(
    title = "B   LDL-C Sub-HR across 9 model specifications",
    subtitle = "All consistent with null (shaded band = ±0.025 around 1.0)",
    x = "LDL-C Sub-HR per 10 mg/dL (95% CI)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    panel.grid.major.x = element_line(color = "gray92")
  )


# ── Combine and save ────────────────────────────────────────────────────────
combined <- pA + pB +
  plot_annotation(
    title = "Statin Bone Protection Is Independent of LDL-C",
    subtitle = "Michigan Genomics Initiative EHR cohort (Fine-Gray competing risks, n ≈ 5,784; 595 osteoporosis events)",
    caption = "Primary model adjusted for: age, sex, race/ethnicity, ever statin, diabetes, hypertension.  LDL-C null across all 9 specifications; statin protective in primary model after LDL-C adjustment.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = "#00274C"),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      plot.caption = element_text(size = 8, color = "gray40")
    )
  )

out <- "/Users/davebrid/Documents/GitHub/PrecisionNutrition/Human Clinical/Calcium-Cholesterol MACE/figures-osteo-cr/statin_vs_ldlc_independence.png"
ggsave(out, combined, width = 13, height = 6, dpi = 150, bg = "white")
cat("Saved:", out, "\n")

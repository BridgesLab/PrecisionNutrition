library(ggplot2)
library(dplyr)
library(patchwork)

# Drug-target MR estimates from ANALYSIS.md
# Serum calcium outcome
ca <- tibble::tribble(
  ~target,     ~label,               ~b,      ~lo,     ~hi,     ~n_snps, ~f_stat,
  "HMGCR",     "HMGCR\n(Statins)",   0.185,   0.155,   0.215,   37,      210,
  "PCSK9",     "PCSK9\n(PCSK9i)",    0.004,  -0.040,   0.047,    3,     2100,
  "NPC1L1",    "NPC1L1\n(Ezetimibe)",-0.040, -0.192,   0.112,    1,       NA,
) |> mutate(
  target = factor(target, levels = c("NPC1L1","PCSK9","HMGCR")),
  sig = abs(b) / ((hi - lo) / (2 * 1.96)) > 1.96
)

# Heel eBMD outcome
bmd <- tibble::tribble(
  ~target,     ~label,                ~b,      ~lo,     ~hi,     ~n_snps, ~f_stat,
  "HMGCR",     "HMGCR\n(Statins)",  -0.115,  -0.149,  -0.081,   37,      210,
  "PCSK9",     "PCSK9\n(PCSK9i)",   -0.014,  -0.069,   0.041,    3,     2100,
  "NPC1L1",    "NPC1L1\n(Ezetimibe)", 0.012, -0.078,   0.102,    1,       NA,
) |> mutate(
  target = factor(target, levels = c("NPC1L1","PCSK9","HMGCR")),
  sig = abs(b) / ((hi - lo) / (2 * 1.96)) > 1.96
)

col_hmgcr  <- "#006600"
col_pcsk9  <- "#AA0000"
col_npc1l1 <- "#774400"
cols <- c("HMGCR" = col_hmgcr, "PCSK9" = col_pcsk9, "NPC1L1" = col_npc1l1)

# ── Panel A: Drug-target → Serum Calcium ────────────────────────────────────
pA <- ggplot(ca, aes(x = b, xmin = lo, xmax = hi, y = target, color = target)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2.5, ymax = 3.5,
           fill = col_hmgcr, alpha = 0.06) +
  geom_errorbarh(height = 0.22, linewidth = 1.0) +
  geom_point(aes(size = sig, shape = sig)) +
  geom_text(aes(label = sprintf("β = %+.3f", b)), hjust = -0.25,
            size = 3.5, color = "black", fontface = "plain") +
  scale_color_manual(values = cols, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 15), guide = "none") +
  scale_size_manual(values  = c(`TRUE` = 4.5, `FALSE` = 3.0), guide = "none") +
  scale_x_continuous(limits = c(-0.32, 0.35),
                     breaks = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3)) +
  scale_y_discrete(labels = c(
    "NPC1L1" = "NPC1L1\n(Ezetimibe)",
    "PCSK9"  = "PCSK9\n(PCSK9i)",
    "HMGCR"  = "HMGCR\n(Statins)"
  )) +
  annotate("label", x = 0.30, y = 3, label = "p < 10⁻²⁰",
           color = col_hmgcr, size = 3, fontface = "bold",
           fill = "#EEFAEE", label.size = 0.3) +
  annotate("label", x = 0.30, y = 2, label = "p = 0.87\n(null)",
           color = col_pcsk9, size = 3, fontface = "bold",
           fill = "#FAEAEA", label.size = 0.3) +
  labs(
    title = "A   Drug-target MR → Serum Calcium",
    subtitle = "Outcomes: UKB + MGI serum calcium (n ≈ 102k)",
    x = "Causal effect on serum calcium (SD units per SD LDL-C equivalent)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.y   = element_text(
      color = c(col_npc1l1, col_pcsk9, col_hmgcr),
      face  = c("plain", "plain", "bold"), size = 12
    ),
    panel.grid.major.x = element_line(color = "gray92")
  )

# ── Panel B: Drug-target → Heel eBMD ────────────────────────────────────────
pB <- ggplot(bmd, aes(x = b, xmin = lo, xmax = hi, y = target, color = target)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2.5, ymax = 3.5,
           fill = col_hmgcr, alpha = 0.06) +
  geom_errorbarh(height = 0.22, linewidth = 1.0) +
  geom_point(aes(size = sig, shape = sig)) +
  geom_text(aes(label = sprintf("β = %+.3f", b)), hjust = -0.25,
            size = 3.5, color = "black", fontface = "plain") +
  scale_color_manual(values = cols, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 15), guide = "none") +
  scale_size_manual(values  = c(`TRUE` = 4.5, `FALSE` = 3.0), guide = "none") +
  scale_x_continuous(limits = c(-0.22, 0.17),
                     breaks = c(-0.2, -0.1, 0, 0.1)) +
  scale_y_discrete(labels = c(
    "NPC1L1" = "NPC1L1\n(Ezetimibe)",
    "PCSK9"  = "PCSK9\n(PCSK9i)",
    "HMGCR"  = "HMGCR\n(Statins)"
  )) +
  annotate("label", x = 0.12, y = 3, label = "p < 10⁻¹²",
           color = col_hmgcr, size = 3, fontface = "bold",
           fill = "#EEFAEE", label.size = 0.3) +
  annotate("label", x = 0.12, y = 2, label = "p = 0.61\n(null)",
           color = col_pcsk9, size = 3, fontface = "bold",
           fill = "#FAEAEA", label.size = 0.3) +
  labs(
    title = "B   Drug-target MR → Heel eBMD",
    subtitle = "Outcome: Morris 2019 heel eBMD (n = 426,824)",
    x = "Causal effect on heel eBMD (SD units per SD LDL-C equivalent)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.y   = element_text(
      color = c(col_npc1l1, col_pcsk9, col_hmgcr),
      face  = c("plain", "plain", "bold"), size = 12
    ),
    panel.grid.major.x = element_line(color = "gray92")
  )

combined <- pA + pB +
  plot_annotation(
    title    = "HMGCR (Mevalonate Pathway) Drives Both Outcomes; PCSK9 Does Not",
    subtitle = "Drug-target MR using cis-acting SNPs as instruments.  NPC1L1 instrument n=1 (imprecise).",
    caption  = "Instrument F-statistics: HMGCR >200, PCSK9 >2,000, NPC1L1 ~25.  HMGCR effect 3–5× larger than genome-wide MR for both outcomes.",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, color = "#00274C"),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      plot.caption  = element_text(size = 8,  color = "gray40")
    )
  )

FIG_BASE <- "/Users/davebrid/Documents/GitHub/PrecisionNutrition"
out <- file.path(FIG_BASE, "Human Genetics/calcium-cholesterol/figures/drug-target-mr-bmd-ca-1.png")
ggsave(out, combined, width = 13, height = 5, dpi = 150, bg = "white")
cat("Saved:", out, "\n")

---
title: "Layer 3 — Formal direct (druggable) effect: cis-MR, MVMR, and colocalization"
author: "Dave Bridges and Katie Kittell"
date: today
format:
  html:
    toc: true
    toc-location: right
    keep-md: true
    code-fold: true
    code-summary: "Show the code"
knitr:
  opts_chunk:
    fig.path: "figures/"
    dev: ["png", "pdf"]
execute:
  echo: true
  warning: false
  message: false
---


::: {.cell}

```{.r .cell-code}
library(TwoSampleMR)
library(ieugwasr)
library(MVMR)
library(coloc)
library(knitr)
source(here::here("R", "helpers.R"))
cfg <- load_config()
set.seed(cfg$seed)
```
:::


## Purpose

This is the layer that produces the **druggable** estimate. For each receptor gene we
**switch the exposure from smoking behavior to the gene's cis-eQTL**, estimate the
behavior-independent (direct) effect on AD via cis-MR and MVMR (receptor expression |
smoking heaviness), and confirm each locus with colocalization. Isolating the receptor
axis also removes the survival/selection collider, which lives on the behavioral axis.

Per-locus **direction logic** (APPROACH.md §7) is recorded explicitly: sign A = allele→
expression, sign B = allele→AD; higher function→lower AD ⇒ agonist/PAM indicated.

## nAChR gene coordinates (GRCh37)


::: {.cell}

```{.r .cell-code}
# hg19/GRCh37 cis windows (gene body +/- cis_window_kb) for the receptor panel.
gene_coords <- tribble(
  ~gene,     ~chr, ~start,     ~end,
  "CHRNA5",   15,  78857862,   78885393,
  "CHRNA3",   15,  78885394,   78919647,
  "CHRNB4",   15,  78919290,   78937340,
  "CHRNA4",   20,  61975397,   62019864,
  "CHRNB2",    1, 154543992,  154579816,
  "CHRNA7",   15,  32322677,   32464722,
  "CHRNA6",    8,  42634258,   42648862,
  "CHRNB3",    8,  42551452,   42591443
) |>
  mutate(win_start = start - cfg$cis_window_kb * 1000,
         win_end   = end   + cfg$cis_window_kb * 1000)
kable(gene_coords, caption = "nAChR cis windows (GRCh37)")
```

::: {.cell-output-display}


Table: nAChR cis windows (GRCh37)

|gene   | chr|     start|       end| win_start|   win_end|
|:------|---:|---------:|---------:|---------:|---------:|
|CHRNA5 |  15|  78857862|  78885393|  78757862|  78985393|
|CHRNA3 |  15|  78885394|  78919647|  78785394|  79019647|
|CHRNB4 |  15|  78919290|  78937340|  78819290|  79037340|
|CHRNA4 |  20|  61975397|  62019864|  61875397|  62119864|
|CHRNB2 |   1| 154543992| 154579816| 154443992| 154679816|
|CHRNA7 |  15|  32322677|  32464722|  32222677|  32564722|
|CHRNA6 |   8|  42634258|  42648862|  42534258|  42748862|
|CHRNB3 |   8|  42551452|  42591443|  42451452|  42691443|


:::
:::


## cis-eQTL instruments (eQTLGen blood) + cis-MR per gene


::: {.cell}

```{.r .cell-code}
# eQTLGen cis-eQTL datasets in OpenGWAS follow eqtl-a-ENSG... ; we resolve by gene symbol.
# Strategy: pull AD region assoc (associations) within each cis window from the AD GWAS,
# and the gene's eQTL via tophits/associations. Where a dedicated eQTL dataset is absent,
# the gene is reported as "no cis instrument available".

# eQTLGen OpenGWAS ids for the receptor panel (ENSG-based). NA = not in OpenGWAS eQTLGen.
eqtl_ids <- c(
  CHRNA5 = "eqtl-a-ENSG00000169684",
  CHRNA3 = "eqtl-a-ENSG00000080644",
  CHRNB4 = "eqtl-a-ENSG00000117971",
  CHRNA4 = "eqtl-a-ENSG00000101204",
  CHRNB2 = "eqtl-a-ENSG00000160716",
  CHRNA7 = "eqtl-a-ENSG00000175344",
  CHRNA6 = "eqtl-a-ENSG00000147434",
  CHRNB3 = "eqtl-a-ENSG00000147432")

run_cis_mr <- function(gene) {
  id <- eqtl_ids[[gene]]
  win <- gene_coords |> filter(gene == !!gene)
  out <- list(gene = gene, status = "ok")
  exp <- tryCatch(
    extract_instruments(id, p1 = 5e-8, clump = TRUE),
    error = function(e) NULL)
  if (is.null(exp) || nrow(exp) == 0) {
    # relax to cis-window tophits at 1e-5 if no genome-wide cis hit
    exp <- tryCatch(extract_instruments(id, p1 = 1e-5, clump = TRUE),
                    error = function(e) NULL)
  }
  if (is.null(exp) || nrow(exp) == 0)
    return(tibble(gene = gene, status = "no cis instrument", nsnp = 0,
                  b = NA, se = NA, pval = NA, signA = NA))
  exp$exposure <- gene
  ado <- extract_outcome_data(exp$SNP, cfg$opengwas$ad_primary, proxies = TRUE, rsq = 0.8)
  hd <- harmonise_data(exp, ado, action = 2) |> filter(mr_keep)
  if (nrow(hd) == 0)
    return(tibble(gene = gene, status = "no harmonised SNP", nsnp = 0,
                  b = NA, se = NA, pval = NA, signA = NA))
  m <- mr(hd, method_list = if (nrow(hd) >= 2) "mr_ivw" else "mr_wald_ratio")
  tibble(gene = gene, status = "ok", nsnp = nrow(hd),
         b = m$b[1], se = m$se[1], pval = m$pval[1],
         # signA: direction of lead allele on expression (beta.exposure of strongest SNP)
         signA = sign(hd$beta.exposure[which.max(abs(hd$beta.exposure))]))
}

cis_mr <- purrr::map_dfr(names(eqtl_ids), run_cis_mr)
write_csv(cis_mr, here::here("results", "cis_mr_per_gene.csv"))
cis_mr |> kable(caption = "Cis-eQTL → AD MR per receptor gene", digits = 3)
```

::: {.cell-output-display}


Table: Cis-eQTL → AD MR per receptor gene

|gene   |status            | nsnp|     b|    se|  pval| signA|
|:------|:-----------------|----:|-----:|-----:|-----:|-----:|
|CHRNA5 |no cis instrument |    0|    NA|    NA|    NA|    NA|
|CHRNA3 |no cis instrument |    0|    NA|    NA|    NA|    NA|
|CHRNB4 |no cis instrument |    0|    NA|    NA|    NA|    NA|
|CHRNA4 |no cis instrument |    0|    NA|    NA|    NA|    NA|
|CHRNB2 |ok                |    1| 0.002| 0.039| 0.955|     1|
|CHRNA7 |no cis instrument |    0|    NA|    NA|    NA|    NA|
|CHRNA6 |no cis instrument |    0|    NA|    NA|    NA|    NA|
|CHRNB3 |no cis instrument |    0|    NA|    NA|    NA|    NA|


:::
:::


## MVMR: receptor expression | smoking heaviness


::: {.cell}

```{.r .cell-code}
# For each gene with a usable cis instrument, build a 2-exposure MVMR:
#   X1 = receptor cis-expression, X2 = smoking heaviness (CPD).
# Instruments = receptor cis-SNPs + genome-wide smoking SNPs (non-receptor loci).
run_mvmr <- function(gene) {
  id <- eqtl_ids[[gene]]
  exp_r <- tryCatch(extract_instruments(id, p1 = 1e-5, clump = TRUE), error = function(e) NULL)
  if (is.null(exp_r) || nrow(exp_r) == 0) return(NULL)
  exp_r$exposure <- "receptor"; exp_r$id.exposure <- "receptor"

  exp_s <- read_csv(here::here("data", "Instruments - cig.liu - Sleep.csv"),
                    show_col_types = FALSE) |>
    mutate(id.exposure = "smoking", exposure = "smoking")

  snps <- union(exp_r$SNP, exp_s$SNP)
  # Outcome and cross-exposure effects
  smk_full <- tryCatch(extract_outcome_data(snps, cfg$opengwas$smk_cpd), error = function(e) NULL)
  rec_full <- tryCatch(extract_outcome_data(snps, id), error = function(e) NULL)
  ado      <- tryCatch(extract_outcome_data(snps, cfg$opengwas$ad_primary, proxies = TRUE, rsq = 0.8),
                       error = function(e) NULL)
  if (any(vapply(list(smk_full, rec_full, ado), is.null, logical(1)))) return(NULL)

  mvdat <- tryCatch(
    mv_harmonise_data(
      exposure_dat = bind_rows(
        exp_r |> transmute(SNP, exposure, id.exposure, beta.exposure, se.exposure,
                           effect_allele.exposure, other_allele.exposure, eaf.exposure,
                           pval.exposure),
        # use smoking GWAS betas for the smoking exposure across all SNPs
        rec_full |> transmute(SNP, exposure = "smoking", id.exposure = "smoking")  # placeholder
      ),
      outcome_dat = ado),
    error = function(e) NULL)
  # NOTE: full MVMR harmonisation requires aligned beta matrices; this helper documents
  # the intended construction. Results below use MVMR::format_mvmr on aligned vectors.
  NULL
}

# Aligned-vector MVMR using TwoSampleMR's mv_extract_exposures for robustness.
run_mvmr2 <- function(gene) {
  id <- eqtl_ids[[gene]]
  mvexp <- tryCatch(
    mv_extract_exposures(c(id, cfg$opengwas$smk_cpd), pval_threshold = 1e-5),
    error = function(e) NULL)
  if (is.null(mvexp) || nrow(mvexp) == 0) return(tibble(gene = gene, status = "no MVMR exposures"))
  mvout <- extract_outcome_data(unique(mvexp$SNP), cfg$opengwas$ad_primary, proxies = TRUE, rsq = 0.8)
  mvdat <- mv_harmonise_data(mvexp, mvout)
  res <- mv_multiple(mvdat)$result

  # Conditional F via MVMR package
  fmt <- tryCatch({
    bx <- mvdat$exposure_beta; sx <- mvdat$exposure_se
    F.data <- format_mvmr(BXGs = bx, BYG = mvdat$outcome_beta,
                          seBXGs = sx, seBYG = mvdat$outcome_se,
                          RSID = rownames(bx))
    sw <- strength_mvmr(F.data, gen_cov = 0)
    as.numeric(sw$exposure1)
  }, error = function(e) NA)

  res |>
    mutate(gene = gene,
           is_receptor = grepl(id, id.exposure, fixed = TRUE),
           cond_F = ifelse(row_number() == 1, fmt, NA)) |>
    dplyr::select(gene, exposure, b, se, pval, cond_F)
}

mvmr_res <- purrr::map_dfr(names(eqtl_ids), function(g)
  tryCatch(run_mvmr2(g), error = function(e) tibble(gene = g, exposure = NA,
                                                    b = NA, se = NA, pval = NA, cond_F = NA)))
write_csv(mvmr_res, here::here("results", "mvmr_direct_effect.csv"))
mvmr_res |> kable(caption = "MVMR direct effect (receptor expression | smoking)", digits = 3)
```

::: {.cell-output-display}


Table: MVMR direct effect (receptor expression | smoking)

|gene   |status            |exposure                                               |      b|    se|  pval|cond_F |
|:------|:-----------------|:------------------------------------------------------|------:|-----:|-----:|:------|
|CHRNA5 |no MVMR exposures |NA                                                     |     NA|    NA|    NA|NA     |
|CHRNA3 |no MVMR exposures |NA                                                     |     NA|    NA|    NA|NA     |
|CHRNB4 |no MVMR exposures |NA                                                     |     NA|    NA|    NA|NA     |
|CHRNA4 |no MVMR exposures |NA                                                     |     NA|    NA|    NA|NA     |
|CHRNB2 |NA                |ENSG00000160716 &#124;&#124; id:eqtl-a-ENSG00000160716 |  0.000| 0.077| 0.998|NA     |
|CHRNB2 |NA                |Cigarettes smoked per day &#124;&#124; id:ieu-b-142    | -0.119| 0.075| 0.112|NA     |
|CHRNA7 |no MVMR exposures |NA                                                     |     NA|    NA|    NA|NA     |
|CHRNA6 |no MVMR exposures |NA                                                     |     NA|    NA|    NA|NA     |
|CHRNB3 |no MVMR exposures |NA                                                     |     NA|    NA|    NA|NA     |


:::
:::


## Colocalization (cis-eQTL vs AD)


::: {.cell}

```{.r .cell-code}
# coloc.abf per gene over the cis window. Pull regional assoc for AD and eQTL.
run_coloc <- function(gene) {
  id <- eqtl_ids[[gene]]
  win <- gene_coords |> filter(gene == !!gene)
  region <- sprintf("%d:%d-%d", win$chr, win$win_start, win$win_end)
  ad_reg  <- tryCatch(associations(region, cfg$opengwas$ad_primary), error = function(e) NULL)
  eq_reg  <- tryCatch(associations(region, id), error = function(e) NULL)
  if (is.null(ad_reg) || is.null(eq_reg) || nrow(ad_reg) == 0 || nrow(eq_reg) == 0)
    return(tibble(gene = gene, status = "no regional data",
                  PP.H3 = NA, PP.H4 = NA, nsnp = 0))
  shared <- intersect(ad_reg$rsid, eq_reg$rsid)
  ad_reg <- ad_reg |> filter(rsid %in% shared) |> distinct(rsid, .keep_all = TRUE)
  eq_reg <- eq_reg |> filter(rsid %in% shared) |> distinct(rsid, .keep_all = TRUE) |>
    arrange(match(rsid, ad_reg$rsid))
  d1 <- list(beta = ad_reg$beta, varbeta = ad_reg$se^2, type = "cc", snp = ad_reg$rsid)
  d2 <- list(beta = eq_reg$beta, varbeta = eq_reg$se^2, type = "quant",
             snp = eq_reg$rsid, sdY = 1)
  cr <- tryCatch(coloc.abf(d1, d2), error = function(e) NULL)
  if (is.null(cr)) return(tibble(gene = gene, status = "coloc failed",
                                 PP.H3 = NA, PP.H4 = NA, nsnp = length(shared)))
  tibble(gene = gene, status = "ok",
         PP.H3 = cr$summary["PP.H3.abf"], PP.H4 = cr$summary["PP.H4.abf"],
         nsnp = length(shared))
}

coloc_res <- purrr::map_dfr(names(eqtl_ids), function(g)
  tryCatch(run_coloc(g), error = function(e) tibble(gene = g, status = "error",
                                                    PP.H3 = NA, PP.H4 = NA, nsnp = 0)))
```

::: {.cell-output .cell-output-stdout}

```
PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
 2.24e-62  3.91e-64  9.79e-01  1.70e-02  4.31e-03 
[1] "PP abf for shared variant: 0.431%"
```


:::

```{.r .cell-code}
write_csv(coloc_res, here::here("results", "coloc_per_gene.csv"))
coloc_res |> kable(caption = "Colocalization (cis-eQTL vs AD) per gene", digits = 3)
```

::: {.cell-output-display}


Table: Colocalization (cis-eQTL vs AD) per gene

|gene   |status           | PP.H3| PP.H4| nsnp|
|:------|:----------------|-----:|-----:|----:|
|CHRNA5 |no regional data |    NA|    NA|    0|
|CHRNA3 |no regional data |    NA|    NA|    0|
|CHRNB4 |no regional data |    NA|    NA|    0|
|CHRNA4 |no regional data |    NA|    NA|    0|
|CHRNB2 |ok               | 0.017| 0.004|  614|
|CHRNA7 |no regional data |    NA|    NA|    0|
|CHRNA6 |no regional data |    NA|    NA|    0|
|CHRNB3 |no regional data |    NA|    NA|    0|


:::
:::


## Target nomination table


::: {.cell}

```{.r .cell-code}
# signB = direction of cis-MR effect on AD (beta sign). Modality logic per APPROACH.md §7:
#   higher function -> lower AD  => agonist/PAM ; higher function -> higher AD => antagonist.
targets <- cis_mr |>
  dplyr::select(gene, cis_b = b, cis_se = se, cis_p = pval, signA, nsnp) |>
  left_join(coloc_res |> dplyr::select(gene, PP.H4), by = "gene") |>
  left_join(mvmr_res |> filter(!is.na(cond_F)) |>
              group_by(gene) |> summarise(cond_F = first(cond_F), .groups = "drop"),
            by = "gene") |>
  mutate(
    signB = sign(cis_b),
    # express expression-increasing-allele effect on AD: align so signA = +1 (higher expr)
    expr_to_AD = signB * signA,                     # +1: higher expr -> higher AD
    implied_modality = case_when(
      is.na(expr_to_AD) ~ NA_character_,
      expr_to_AD < 0 ~ "agonist / PAM (higher function protective)",
      expr_to_AD > 0 ~ "antagonist (higher function harmful)",
      TRUE ~ "ambiguous"),
    chrfam7a_flag = gene == "CHRNA7",
    nominated = !is.na(PP.H4) & PP.H4 >= cfg$thresholds$coloc_h4_min &
                !is.na(cis_p) & cis_p < 0.05)

write_tsv(targets, here::here("results", "targets.tsv"))
targets |> kable(caption = "Target nomination (Layer 3)", digits = 3)
```

::: {.cell-output-display}


Table: Target nomination (Layer 3)

|gene   | cis_b| cis_se| cis_p| signA| nsnp| PP.H4|cond_F | signB| expr_to_AD|implied_modality                     |chrfam7a_flag |nominated |
|:------|-----:|------:|-----:|-----:|----:|-----:|:------|-----:|----------:|:------------------------------------|:-------------|:---------|
|CHRNA5 |    NA|     NA|    NA|    NA|    0|    NA|NA     |    NA|         NA|NA                                   |FALSE         |FALSE     |
|CHRNA3 |    NA|     NA|    NA|    NA|    0|    NA|NA     |    NA|         NA|NA                                   |FALSE         |FALSE     |
|CHRNB4 |    NA|     NA|    NA|    NA|    0|    NA|NA     |    NA|         NA|NA                                   |FALSE         |FALSE     |
|CHRNA4 |    NA|     NA|    NA|    NA|    0|    NA|NA     |    NA|         NA|NA                                   |FALSE         |FALSE     |
|CHRNB2 | 0.002|  0.039| 0.955|     1|    1| 0.004|NA     |     1|          1|antagonist (higher function harmful) |FALSE         |FALSE     |
|CHRNA7 |    NA|     NA|    NA|    NA|    0|    NA|NA     |    NA|         NA|NA                                   |TRUE          |FALSE     |
|CHRNA6 |    NA|     NA|    NA|    NA|    0|    NA|NA     |    NA|         NA|NA                                   |FALSE         |FALSE     |
|CHRNB3 |    NA|     NA|    NA|    NA|    0|    NA|NA     |    NA|         NA|NA                                   |FALSE         |FALSE     |


:::
:::


## Decision


::: {.cell}

```{.r .cell-code}
nom <- targets |> filter(nominated)
log_decision("L3", "Genes nominated (coloc H4>=0.8 & cis-MR p<0.05)",
             sprintf("%d gene(s)", nrow(nom)),
             ifelse(nrow(nom) > 0, paste(nom$gene, collapse = ","), "none survive"))
if (nrow(nom) > 0) {
  cat("Nominated targets:\n"); print(nom |> dplyr::select(gene, expr_to_AD, implied_modality, PP.H4, cond_F))
} else cat("No genes survive Layer 3 nomination thresholds.\n")
```

::: {.cell-output .cell-output-stdout}

```
No genes survive Layer 3 nomination thresholds.
```


:::
:::


::: callout-note
**Honesty flags.** eQTL captures *expression*, not channel *function* — a real functional
effect can be invisible here. MVMR with shared cis instruments is prone to conditional
weak-instrument bias; report conditional F honestly. CHRNA7 results are sensitivity-
qualified due to the CHRFAM7A duplication (see Layer 4 / APPROACH.md §9).
:::

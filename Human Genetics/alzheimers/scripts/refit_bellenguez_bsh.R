# refit_bellenguez_bsh.R — fit the Bellenguez 2022 SlopeHunter selection slope (b_SH) in this folder.
#
# The primary b_SH (-0.944) was fit by the smoking pipeline (smoking/04e_overlap_indexevent.qmd,
# smoking/R/prep_indexevent_genomewide.R) and the AD pipeline read it from
# smoking/results/indexevent_slopehunter_fits.csv. This script refits it with this folder's engine
# (R/fit_selection_slope.R), using the recipe already used here for Kunkle and Uffelmann
# (README, "Fitting b_SH"). That way the primary correction has its provenance in this folder.
#
# Inputs (gitignored, data/cache/sumstats/):
#   AD_bellenguez_harmonised.tsv.gz      GWAS Catalog harmonised Bellenguez 2022 (GCST90027158, Build38):
#     https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/
#       GCST90027158/harmonised/35379992-GCST90027158-MONDO_0004975-Build38.f.tsv.gz
#     (smoking/scripts/fetch_sumstats.sh downloads it to the same name)
#   LIFESPAN_pilling_harmonised.tsv.gz   parental lifespan (GCST006697), as used by every b_SH fit
#   data/cache/EUR.{bed,bim,fam}         1000G EUR LD panel
# Output:
#   R/fits/bellenguez_lifespan_slopehunter_fits.csv   (same layout as the other R/fits/ files)
#   prints a comparison with the smoking fit
#
# Run from Human Genetics/alzheimers:   Rscript scripts/refit_bellenguez_bsh.R

source("R/fit_selection_slope.R")
cfg <- yaml::read_yaml("pipeline/config.yml")

ss       <- "data/cache/sumstats"
ad_file  <- file.path(ss, "AD_bellenguez_harmonised.tsv.gz")
life_file <- file.path(ss, "LIFESPAN_pilling_harmonised.tsv.gz")
clump_in <- file.path(ss, "bellenguez_clump_input.tsv")
clump_out <- file.path(ss, "bellenguez_clumped")
out_csv  <- "R/fits/bellenguez_lifespan_slopehunter_fits.csv"
plink2   <- Sys.getenv("PLINK2_BIN", path.expand(cfg$paths$plink2))

stopifnot(file.exists(ad_file), file.exists(life_file), file.exists(plink2))
message("md5 ", basename(ad_file), ": ", unname(tools::md5sum(ad_file)))
message("md5 ", basename(life_file), ": ", unname(tools::md5sum(life_file)))

# 1. merge + allele-align; clump input is every merged SNP with its lifespan p-value
m <- prep_selection_merge(ad_file, life_file)
data.table::fwrite(m[, .(ID = SNP, P = p_life)], clump_in, sep = "\t")

# 2. LD-clump on the lifespan axis (same settings as Kunkle, Uffelmann and the smoking fit)
status <- system2(plink2, c("--bfile", cfg$paths$ld_bfile, "--clump", clump_in,
                            "--clump-p1", "1e-3", "--clump-p2", "0.01", "--clump-r2", "0.01",
                            "--clump-kb", "1000", "--out", clump_out))
stopifnot(status == 0)

# 3. fit SlopeHunter (full and APOE-excluded PRIMARY sets; seed 20260604)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
fits <- build_selection_slope(ad_file = ad_file, life_file = life_file,
                              clumped_snps = paste0(clump_out, ".clumps"), out_csv = out_csv)

# 4. compare with the smoking fit the pipeline has used so far
old <- readr::read_csv("smoking/results/indexevent_slopehunter_fits.csv", show_col_types = FALSE)
cmp <- dplyr::bind_rows(dplyr::mutate(old, source = "smoking pipeline"),
                        dplyr::mutate(fits, source = "this folder (refit)")) |>
  dplyr::filter(grepl("PRIMARY", set)) |>
  dplyr::select(source, set, n_snp, n_fit, b_SH, se_SH, ci_lo, ci_hi)
print(as.data.frame(cmp), digits = 10)
message("Identical to the smoking fit: ",
        isTRUE(all.equal(cmp$b_SH[1], cmp$b_SH[2], tolerance = 1e-10)) &&
          isTRUE(all.equal(cmp$se_SH[1], cmp$se_SH[2], tolerance = 1e-10)) &&
          cmp$n_fit[1] == cmp$n_fit[2])

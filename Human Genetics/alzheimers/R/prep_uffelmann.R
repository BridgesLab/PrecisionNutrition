# prep_uffelmann.R — bring the Uffelmann 2026 (PGC-ALZ3) AD summary statistics into the pipeline.
#
# Inputs (CC BY-NC-SA 4.0; downloaded by hand from the CNCR SURF share into data/cache/sumstats/,
# which is gitignored; 23andMe excluded by the authors):
#   newc_combined_eur_neff0.6_nsumstats1.txt.gz        cohorts new to PGC-ALZ (independent of Bellenguez)
#   twent_case_control_eur_neff0.6_nsumstats1.txt.gz   main analysis, clinical cases only
#   twent_combined_eur_neff0.6_nsumstats1.txt.gz       main analysis, clinical + proxy cases (combined)
#
# 1. Rewrites each file in the GWAS-Catalog-harmonised layout read_gwas_tsv() and
#    prep_selection_merge() expect (variant_id = rsID, GRCh37).
# 2. Fits the SlopeHunter survival slope against parental lifespan exactly as for Kunkle
#    (clump on lifespan p: plink2 --clump-p1 1e-3 --clump-p2 0.01 --clump-r2 0.01 --clump-kb 1000),
#    and prints the orientation and identifiability checks used to decide `trusted` in config.yml.
#
# Run from alzheimers/:  Rscript R/prep_uffelmann.R

suppressPackageStartupMessages({ library(data.table); library(dplyr) })
source("R/fit_selection_slope.R")

sumstats <- "data/cache/sumstats"
life_file <- file.path(sumstats, "LIFESPAN_pilling_harmonised.tsv.gz")
plink2 <- path.expand(Sys.getenv("PLINK2_BIN", "~/bin/plink2"))
files <- c(uffelmann_newc = "newc_combined_eur_neff0.6_nsumstats1.txt.gz",
           uffelmann_cc   = "twent_case_control_eur_neff0.6_nsumstats1.txt.gz",
           uffelmann_all  = "twent_combined_eur_neff0.6_nsumstats1.txt.gz")
# Optional: Rscript R/prep_uffelmann.R uffelmann_all   (process only the named files)
if (length(args <- commandArgs(trailingOnly = TRUE))) files <- files[args]

for (k in names(files)) {
  message("== ", k)
  out <- file.path(sumstats, sprintf("AD_%s_harmonised.tsv.gz", k))
  if (!file.exists(out)) {
    d <- fread(file.path(sumstats, files[[k]]),
               select = c("rsid", "chromosome", "base_pair_location", "effect_allele", "other_allele",
                          "effect_allele_frequency", "beta", "standard_error", "p_value", "neff"),
               colClasses = list(character = "p_value"), showProgress = FALSE)
    setnames(d, c("rsid", "neff"), c("variant_id", "n"))
    d <- d[grepl("^rs[0-9]+$", variant_id)]
    fwrite(d, out, sep = "\t", na = "NA")
    message("  wrote ", out, " (", nrow(d), " variants; median neff ", round(median(d$n)), ")")
  }

  print(check_ad_orientation_anchors(out)); message("  orientation: ", attr(check_ad_orientation_anchors(out), "verdict"))

  m <- prep_selection_merge(out, life_file)
  clump_in <- file.path(sumstats, sprintf("%s_clump_input.tsv", k))
  fwrite(m[, .(ID = SNP, P = p_life)], clump_in, sep = "\t")
  system2(plink2, c("--bfile", "data/cache/EUR", "--clump", clump_in, "--clump-kb", "1000",
                    "--clump-p1", "1e-3", "--clump-p2", "0.01", "--clump-r2", "0.01",
                    "--out", file.path(sumstats, paste0(k, "_clumped"))), stdout = FALSE)
  clumped <- fread(file.path(sumstats, paste0(k, "_clumped.clumps")))$ID
  mc <- m[SNP %in% clumped]
  print(bsh_identifiability(mc))
  fits <- fit_selection_slope(mc)
  readr::write_csv(fits, sprintf("R/fits/%s_lifespan_slopehunter_fits.csv", k))
  print(fits)
  print(bsh_null_permutation(mc[apoe == FALSE]))
}

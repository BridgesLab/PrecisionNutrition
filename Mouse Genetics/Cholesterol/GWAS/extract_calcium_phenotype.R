## Extract serum calcium phenotypes from Svenson DO HFD v12 for GEMMA
##
## Writes one file per phenotype, one value per line, in the same sample
## order as Genotypes_all.fam (= Genotypes_all.bimbam).  Missing values
## are written as "NA" (GEMMA accepts this for LMM).
##
## Run on Great Lakes:
##   module load R/4.3.1
##   Rscript extract_calcium_phenotype.R

GEMMA_DIR <- "/home/davebrid/bridgeslab/GEMMA/DO_HFD"

# ── Load phenotype dataset ───────────────────────────────────────────────────
phen <- readRDS(file.path(GEMMA_DIR, "dataset.phenotype.Svenson_DO_HFD.v12.Rds"))

# Rank-normalised matrix (rows = samples, cols = phenotypes)
pheno_rz  <- phen$data$rz
pheno_raw <- phen$data$raw
cat("Phenotype matrix dimensions:", dim(pheno_rz), "\n")

# ── Sample order from .fam file (must match bimbam genotype order) ───────────
fam <- read.table(
  file.path(GEMMA_DIR, "Genotypes_all.fam"),
  header = FALSE,
  col.names = c("fid", "iid", "pid", "mid", "sex", "pheno")
)
cat("FAM samples:", nrow(fam), "\n")

# Align: phenotype rows -> fam order (some mice may lack phenotype data -> NA)
idx <- match(fam$iid, rownames(pheno_rz))
cat("Matched:", sum(!is.na(idx)), "of", nrow(fam), "samples\n")
if (sum(!is.na(idx)) < nrow(fam) * 0.9)
  warning("More than 10% of FAM samples are missing from phenotype dataset")

# ── Helper: extract one phenotype in fam order and write ────────────────────
write_pheno <- function(mat, col, outfile) {
  if (!col %in% colnames(mat)) stop("Column not found: ", col)
  vals <- mat[idx, col]          # NA for unmatched samples
  writeLines(ifelse(is.na(vals), "NA", as.character(vals)), outfile)
  cat("Wrote", outfile, "— n non-NA:", sum(!is.na(vals)), "\n")
}

# ── Write calcium phenotype files ───────────────────────────────────────────
# Primary: 19-week calcium (taken alongside HFD) — rank-normalised
write_pheno(pheno_rz,  "calcium2",        file.path(GEMMA_DIR, "Calcium2_all.txt"))

# Secondary / sensitivity
write_pheno(pheno_rz,  "calcium1",        file.path(GEMMA_DIR, "Calcium1_all.txt"))
write_pheno(pheno_rz,  "non.fast.calcium", file.path(GEMMA_DIR, "Calcium_nonfasting_all.txt"))

# Raw values alongside (useful for effect-size interpretation)
write_pheno(pheno_raw, "calcium2",        file.path(GEMMA_DIR, "Calcium2_raw_all.txt"))
write_pheno(pheno_raw, "calcium1",        file.path(GEMMA_DIR, "Calcium1_raw_all.txt"))

cat("\nDone. Primary GWAS phenotype: Calcium2_all.txt (19-week, rank-normalised)\n")

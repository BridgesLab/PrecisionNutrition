# bsh_propagate_existing.R — apply the b_SH uncertainty propagation (pipeline/R/selection.R) to the
# result CSVs already on disk, and write results/pipeline/bsh_grid.csv.
#
# The pipeline does this itself when it writes these files (csv_results, csv_contrast_results,
# csv_contrast_groups, csv_cis_results, csv_pf and csv_bsh_grid in _targets.R). This script exists
# because the CSVs in git come from Great Lakes runs and the local targets store is older, so
# rebuilding csv_results locally would write stale results. It is idempotent (already-propagated
# rows are restored first) and rewrites only the cells whose values change: the SlopeHunter rows'
# se / p / CI, plus the added *_bsh columns.
#
# Run from Human Genetics/alzheimers:   Rscript pipeline/bsh_propagate_existing.R

suppressPackageStartupMessages({ library(dplyr); library(readr) })
source("pipeline/R/selection.R")

cfg     <- yaml::read_yaml("pipeline/config.yml")
res_dir <- cfg$paths$results_dir
con_dir <- file.path(res_dir, "contrast")
par <- bsh_params(load_bsh_table(cfg$slopehunter$fits),
                  read_csv(file.path(res_dir, "bsh_null_checks.csv"), show_col_types = FALSE))
print(par)

tables <- tibble::tribble(
  ~table,                    ~path,                                                   ~keys,
  "results",                 file.path(res_dir, "results.csv"),                       "results",
  "contrast",                file.path(con_dir, "results.csv"),                       "contrast",
  "cluster",                 file.path(con_dir, "cluster_results.csv"),               "cluster",
  "drug_target",             file.path(con_dir, "drug_target_results.csv"),           "drug_target",
  "post_freeze_drug_target", file.path(con_dir, "post_freeze_drug_target_results.csv"), "drug_target")

num_cols <- c("b", "se", "p", "ci_lo", "ci_hi", "p_gamma_below_0", BSH_ADDED)

# Every column as text, so cells that do not change are written back byte for byte.
read_text <- function(path) read_csv(path, col_types = cols(.default = col_character()), na = character(),
                                     name_repair = "minimal", progress = FALSE)
to_value <- function(x, logical = FALSE) {
  x[x %in% c("NA", "")] <- NA
  if (logical) as.logical(x) else suppressWarnings(as.numeric(x))
}

grids <- purrr::pmap(tables, function(table, path, keys) {
  txt  <- read_text(path)
  work <- txt
  for (col in intersect(num_cols, names(work))) work[[col]] <- to_value(work[[col]])
  if ("ci_excludes_0" %in% names(work)) work$ci_excludes_0 <- to_value(work$ci_excludes_0, logical = TRUE)
  new <- propagate_bsh(work, BSH_KEYS[[keys]], par)
  out <- new
  for (col in names(new)) {
    v <- new[[col]]
    if (is.character(v)) next
    s <- ifelse(is.na(v), NA_character_, as.character(v))
    if (col %in% names(txt)) {
      # Unchanged within 1e-12 (relative) keeps the original text: the *_cond_bsh columns are
      # written with 15 digits, so a re-run restores the conditional values with ~1e-16 noise.
      o <- to_value(txt[[col]], logical = is.logical(v))
      same <- (is.na(o) & is.na(v)) |
        (!is.na(o) & !is.na(v) & (o == v | abs(o - v) <= 1e-12 * pmax(abs(o), abs(v))))
      s[same] <- txt[[col]][same]
    }
    out[[col]] <- s
  }
  out[] <- lapply(out, \(x) replace(x, x %in% "NA", NA))
  write_csv(out, path, na = "NA")
  n <- sum(!is.na(new$dtheta_dbsh))
  message(sprintf("%-24s %4d SlopeHunter rows propagated  (%s)", table, n, path))
  bsh_sensitivity_grid(work, BSH_KEYS[[keys]], par, table)
})

grid <- bind_rows(grids)
write_csv(grid, file.path(res_dir, "bsh_grid.csv"))
message("bsh_grid.csv: ", nrow(grid), " rows")

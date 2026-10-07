# _targets.R — univariable robust-MR pipeline + the existing AD notebooks.
#
# Run from this folder:
#   targets::tar_make()                         # everything (CAUSE targets need Great Lakes)
#   targets::tar_make(names = !starts_with("cause_") & !starts_with("nb_"))   # local, no CAUSE
#   targets::tar_visnetwork()                   # dependency graph
#   targets::tar_read(results)                  # final results table
# On Great Lakes:  PIPELINE_PROFILE=greatlakes  (see pipeline/README.md)
#
# What is where:
#   pipeline/config.yml   datasets, associations, thresholds, Great Lakes resources
#   pipeline/R/*.R        read/QC, instruments, classical suite, MRAID, CAUSE, SlopeHunter, tables
#   R/fit_selection_slope.R                       SlopeHunter engine (Step 2)
#   ../calcium-cholesterol/R/robust_mr_helpers.R  CAUSE posterior summaries (calcium arm)

library(targets)
library(tarchetypes)
library(crew)

tar_source(c("pipeline/R", "R/fit_selection_slope.R", "../calcium-cholesterol/R/robust_mr_helpers.R"))

cfg     <- yaml::read_yaml("pipeline/config.yml")
ds      <- dataset_table(cfg)
assoc   <- association_table(cfg, ds)
profile <- Sys.getenv("PIPELINE_PROFILE", "local")

# Config sections as separate globals, so editing one section only invalidates its own targets.
CFG_QC    <- cfg$qc
CFG_INST  <- cfg$instruments
CFG_METH  <- cfg$methods
CFG_SH    <- cfg$slopehunter
LD_BFILE  <- cfg$paths$ld_bfile
GWAS_DIR  <- cfg$paths$gwas_dir
RES_DIR   <- cfg$paths$results_dir
FIG_DIR   <- file.path(cfg$paths$results_dir, "figures")
HM3_PATH  <- cfg$paths$hapmap3
PLINK2    <- cfg$paths$plink2     # resolved inside commands, so the machine-specific path is not a dependency
SEED      <- cfg$seed
CFG_APSS  <- cfg$mrapss
CFG_APSS_GW <- utils::modifyList(cfg$mrapss, cfg$mrapss_sensitivity)   # p < 5e-8 sensitivity
LDSC_DIR  <- cfg$paths$ldscores
RUN_MRAID <- isTRUE(cfg$analysis$run_mraid)
MRAID_MAX <- cfg$analysis$mraid_rule$max_candidates %||% Inf
mraid_eligible <- function(a) !a$classical_only &
  (a$overlap == "none" | !isTRUE(cfg$analysis$mraid_rule$require_no_overlap))

# ---- controllers: local workers always; Slurm on Great Lakes ---------------------------------
local_ctl <- crew_controller_local(name = "local", workers = 2, seconds_idle = 60)
if (profile == "greatlakes") {
  gl <- cfg$greatlakes
  slurm_ctl <- crew.cluster::crew_controller_slurm(
    name = "slurm", workers = gl$workers, seconds_idle = 300,
    options_cluster = crew.cluster::crew_options_slurm(
      script_lines = c(sprintf("#SBATCH --account=%s", gl$account), gl$script_lines),
      memory_gigabytes_required = gl$memory_gb, cpus_per_task = gl$cpus,
      time_minutes = gl$time_minutes, partition = gl$partition,
      log_output = "logs/crew_%A.out", log_error = "logs/crew_%A.err"))
  controller <- crew_controller_group(local_ctl, slurm_ctl)
  heavy <- tar_resources(crew = tar_resources_crew(controller = "slurm"))
} else {
  controller <- local_ctl
  heavy <- tar_resources(crew = tar_resources_crew(controller = "local"))
}

tar_option_set(
  packages = c("data.table", "dplyr", "tidyr", "purrr", "TwoSampleMR"),
  format = "qs", memory = "transient", garbage_collection = TRUE,
  controller = controller, seed = SEED)

# ---- static branching tables ------------------------------------------------------------------
sym_vec <- function(x) rlang::syms(x)

ds_values <- ds |> dplyr::select(key, id, source, file, n, p_from_beta)

trusted_out <- names(purrr::keep(CFG_SH$fits, \(f) isTRUE(f$trusted)))
sh_values <- tibble::tibble(out = trusted_out, g = sym_vec(paste0("gwas_", trusted_out)))

exp_all   <- unique(assoc$exposure)
exp_mraid <- if (RUN_MRAID) unique(assoc$exposure[mraid_eligible(assoc)]) else character()

assoc <- assoc |>
  dplyr::left_join(ds |> dplyr::select(exposure = key, n_exp = n), by = "exposure") |>
  dplyr::left_join(ds |> dplyr::select(outcome = key, n_out = n), by = "outcome") |>
  dplyr::mutate(g_exp = sym_vec(paste0("gwas_", exposure)),
                g_out = sym_vec(ifelse(arm == "slopehunter", paste0("gwas_sh_", outcome),
                                       paste0("gwas_", outcome))),
                inst_c = sym_vec(paste0("inst_classical_", exposure)),
                inst_m = sym_vec(paste0("inst_mraid_", exposure)))
cols <- c("assoc_id", "exposure", "outcome", "arm", "n_exp", "n_out", "g_exp", "g_out")
v_classical <- assoc[, c(cols, "inst_c")]
v_apss      <- assoc[assoc$primary_method == "MR-APSS", cols]
v_cause     <- assoc[assoc$sensitivity_method == "CAUSE", cols]
v_mraid     <- if (RUN_MRAID) assoc[mraid_eligible(assoc), c(cols, "inst_m")] else assoc[0, c(cols, "inst_m")]

# ---- targets ----------------------------------------------------------------------------------
t_reference <- list(
  tar_target(config_file, "pipeline/config.yml", format = "file"),
  tar_target(bim_file, paste0(LD_BFILE, ".bim"), format = "file"),
  tar_target(ld_bim, read_bim(sub("\\.bim$", "", bim_file))),
  tar_target(hm3_file, HM3_PATH, format = "file"),
  tar_target(hm3, data.table::fread(hm3_file)$SNP),
  tar_target(ldsc_dir, LDSC_DIR, format = "file"),
  tar_target(ldsc, read_ldscores(ldsc_dir))
)

# One download + read/QC per dataset. Local files are tracked as files; OpenGWAS ones are fetched.
t_datasets <- tar_map(
  values = ds_values, names = "key", unlist = FALSE,
  tar_target(gwas_file, if (source == "opengwas_vcf") fetch_opengwas_vcf(id, GWAS_DIR) else file,
             format = "file"),
  tar_target(gwas, load_gwas(gwas_file, list(key = key, source = source, n = n, p_from_beta = p_from_beta),
                             ld_bim$SNP, CFG_QC)),
  tar_target(gwas_qc, attr(gwas, "qc_log"))
)

# SlopeHunter: Step 2's fits, the genome-wide adjusted outcome, and Bellenguez's open checks.
t_selection <- list(
  tar_target(bsh_files, purrr::keep(purrr::map_chr(CFG_SH$fits, "fits_csv"), file.exists), format = "file"),
  tar_target(bsh_table, { bsh_files; load_bsh_table(CFG_SH$fits) }),
  tar_map(values = sh_values, names = "out", unlist = FALSE,
    tar_target(gwas_sh, adjust_outcome_gw(g, gwas_lifespan, bsh_table$b_SH[bsh_table$dataset == out])),
    tar_target(sh_merge, selection_merge(g, gwas_lifespan, LD_BFILE, resolve_plink2(PLINK2)),
               resources = heavy),
    tar_target(bsh_check, bsh_checks(sh_merge, bsh_table$b_SH[bsh_table$dataset == out],
                                     n_perm = CFG_SH$null_permutations) |>
                 dplyr::mutate(dataset = out, .before = 1), resources = heavy)
  )
)

# Instrument sets per exposure.
t_instruments <- list(
  tar_map(values = tibble::tibble(key = exp_all, g = sym_vec(paste0("gwas_", exp_all))),
          names = "key", unlist = FALSE,
          tar_target(inst_classical, select_classical(g, key, CFG_INST$classical, LD_BFILE,
                                                      resolve_plink2(PLINK2)))),
  if (length(exp_mraid))
    tar_map(values = tibble::tibble(key = exp_mraid, g = sym_vec(paste0("gwas_", exp_mraid))),
            names = "key", unlist = FALSE,
            tar_target(inst_mraid, select_mraid(g, key, CFG_INST$mraid, LD_BFILE,
                                                resolve_plink2(PLINK2), hm3)))
)

# Per association: classical suite (+ MRBEE/RAPS stability flags) on every row; MR-APSS (main) and
# CAUSE (sensitivity) on every non-lifespan row. MRAID only if analysis$run_mraid is true.
m_classical <- tar_map(values = v_classical, names = "assoc_id", unlist = FALSE,
  tar_target(classical, analyse_classical(inst_c, g_exp, g_out, assoc_id, exposure, outcome, arm,
                                          CFG_METH$presso_nboot, SEED)),
  tar_target(stability, stability_flags(classical)),
  # Scatter, funnel and leave-one-out plots + LOO table (results/pipeline/figures/<assoc_id>/).
  tar_target(diag_classical, diagnostics_classical(classical, FIG_DIR), format = "file"))
m_apss <- tar_map(values = v_apss, names = "assoc_id", unlist = FALSE,
  tar_target(apss, analyse_apss(g_exp, g_out, assoc_id, exposure, outcome, arm, hm3, ldsc, n_exp,
                                n_out, LD_BFILE, resolve_plink2(PLINK2), CFG_APSS),
             resources = heavy),
  tar_target(diag_apss, diagnostics_apss(apss, FIG_DIR), format = "file"),
  # Threshold sensitivity: genome-wide-significant instruments only (config mrapss_sensitivity).
  tar_target(apss_gw, analyse_apss(g_exp, g_out, assoc_id, exposure, outcome, arm, hm3, ldsc, n_exp,
                                   n_out, LD_BFILE, resolve_plink2(PLINK2), CFG_APSS_GW),
             resources = heavy))
m_mraid <- if (RUN_MRAID) tar_map(values = v_mraid, names = "assoc_id", unlist = FALSE,
  tar_target(mraid, analyse_mraid(inst_m, g_out, assoc_id, exposure, outcome, arm, ld_bim, LD_BFILE,
                                  resolve_plink2(PLINK2), n_exp, n_out, CFG_METH[["mraid"]],
                                  max_candidates = MRAID_MAX),
             resources = heavy)) else list(mraid = list())
m_cause <- tar_map(values = v_cause, names = "assoc_id", unlist = FALSE,
  tar_target(cause, analyse_cause(g_exp, g_out, assoc_id, exposure, outcome, arm, hm3, n_exp, n_out,
                                  LD_BFILE, resolve_plink2(PLINK2), CFG_METH[["cause"]], SEED),
             resources = heavy,
             # Locally CAUSE cannot run (loo >= 2.10); record the error and keep going.
             error = if (profile == "greatlakes") "stop" else "null"))

pluck_all <- function(x, field) dplyr::bind_rows(purrr::map(purrr::compact(x), field))

t_tables <- list(
  tar_combine(classical_all, m_classical[["classical"]], command = list(!!!.x)),
  tar_combine(stability_all, m_classical[["stability"]], command = dplyr::bind_rows(!!!.x)),
  tar_combine(apss_all, m_apss[["apss"]], command = list(!!!.x)),
  tar_combine(apss_gw_all, m_apss[["apss_gw"]], command = list(!!!.x)),
  if (RUN_MRAID) tar_combine(mraid_all, m_mraid[["mraid"]], command = list(!!!.x))
  else tar_target(mraid_all, list()),
  tar_combine(cause_all, m_cause[["cause"]], command = list(!!!.x)),
  tar_target(results, assemble_results(assoc, pluck_all(classical_all, "rows"),
                                       dplyr::bind_rows(pluck_all(apss_all, "rows"),
                                                        pluck_all(apss_gw_all, "rows"),
                                                        pluck_all(cause_all, "rows"),
                                                        pluck_all(mraid_all, "rows")),
                                       stability_all)),
  # Instrument table for the MAIN method of each association (MR-APSS; the classical set for
  # classical-only rows). Every association's classical set is in instruments_classical.
  tar_target(instruments_primary, dplyr::bind_rows(pluck_all(apss_all, "instruments"),
                                                   pluck_all(classical_all, "instruments") |>
                                                     dplyr::filter(assoc_id %in% assoc$assoc_id[assoc$primary_method == "none"]))),
  tar_target(instruments_classical, pluck_all(classical_all, "instruments")),
  tar_combine(dataset_qc, t_datasets[["gwas_qc"]], command = dplyr::bind_rows(!!!.x)),
  tar_target(attrition, dplyr::bind_rows(
    dataset_qc |> dplyr::transmute(association = dataset, instrument_set = "dataset QC", step,
                                    n_after = as.integer(n_variants)),
    pluck_all(classical_all, "attrition"), pluck_all(apss_all, "attrition"),
    pluck_all(apss_gw_all, "attrition"),
    pluck_all(cause_all, "attrition"), pluck_all(mraid_all, "attrition")) |> dplyr::distinct()),
  tar_target(apss_background, pluck_all(apss_all, "background")),
  tar_combine(bsh_checks_all, t_selection[[3]][["bsh_check"]], command = dplyr::bind_rows(!!!.x)),
  tar_target(step2_gap, step2_gap_fill(results, bsh_table)),
  tar_target(datasets_used, ds |> dplyr::mutate(cohorts = purrr::map_chr(cohorts, paste, collapse = "+"))),
  tar_target(cause_weights, pluck_all(cause_all, "weights")),
  # Latent residual heritable confounding (latent.R): CAUSE q*eta and MR-APSS background slope,
  # one row per association x method. The MR-APSS background (Omega, C) does not depend on the
  # instrument threshold, so it is taken from the main (p < 5e-5) fit.
  tar_target(latent_confounding, latent_table(pluck_all(cause_all, "rows"), pluck_all(apss_all, "rows"))),
  # CAUSE grid posteriors, kept so the latent terms can be recomputed without another CAUSE run.
  tar_target(cause_posteriors, purrr::set_names(purrr::map(purrr::compact(cause_all), "posterior"),
                                                purrr::map_chr(purrr::compact(cause_all), \(x) x$rows$assoc_id[1]))),

  # CSV outputs (tracked in git; SNP lists are semicolon-joined in `removed`).
  tar_target(csv_results, write_csv_target(results, file.path(RES_DIR, "results.csv")), format = "file"),
  tar_target(csv_instruments, write_csv_target(instruments_primary, file.path(RES_DIR, "instruments_primary.csv")), format = "file"),
  tar_target(csv_instruments_classical, write_csv_target(instruments_classical, file.path(RES_DIR, "instruments_classical.csv")), format = "file"),
  tar_combine(figures_all, m_classical[["diag_classical"]], m_apss[["diag_apss"]], command = c(!!!.x), format = "file"),
  tar_target(csv_attrition, write_csv_target(attrition, file.path(RES_DIR, "instrument_attrition.csv.gz")), format = "file"),
  tar_target(csv_apss_background, write_csv_target(apss_background, file.path(RES_DIR, "mrapss_background.csv")), format = "file"),
  tar_target(csv_associations, write_csv_target(assoc |> dplyr::select(assoc_id, block, exposure, outcome, arm, overlap, primary_method, sensitivity_method),
                                                file.path(RES_DIR, "associations.csv")), format = "file"),
  tar_target(csv_datasets, write_csv_target(datasets_used, file.path(RES_DIR, "datasets.csv")), format = "file"),
  tar_target(csv_step2_gap, write_csv_target(step2_gap, file.path(RES_DIR, "step2_lifespan_gap_fill.csv")), format = "file"),
  tar_target(csv_bsh_checks, write_csv_target(bsh_checks_all, file.path(RES_DIR, "bsh_null_checks.csv")), format = "file"),
  tar_target(csv_latent, write_csv_target(latent_confounding, file.path(RES_DIR, "latent_confounding.csv")), format = "file"),
  tar_target(rds_cause_posteriors, { dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)
                                     f <- file.path(RES_DIR, "cause_posteriors.rds"); saveRDS(cause_posteriors, f); f },
             format = "file")
)

# ---- existing notebooks -----------------------------------------------------------------------
# Rendered as pipeline steps so their order is enforced (Step 1 writes the CSV Step 2 reads). They
# still pull from the OpenGWAS API themselves; while /tophits is down they error, and error =
# "continue" lets the rest of the pipeline finish. The MVMR notebook joins here when it runs.
if (!nzchar(Sys.which("quarto")) && file.exists("/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto"))
  Sys.setenv(QUARTO_PATH = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto")
render_nb <- function(path, outputs) { quarto::quarto_render(path, quiet = TRUE); c(path, outputs) }
t_notebooks <- list(
  tar_target(nb_step1_survival_screen,
             render_nb("exposure_survival_screen.qmd",
                       c("exposure_survival_screen.csv", "exposure_survival_verdict.csv", "ad_longevity_arm.csv")),
             format = "file", error = "continue", deployment = "main"),
  tar_target(nb_step2_indexevent,
             { nb_step1_survival_screen; bsh_files
               render_nb("exposure_ad_indexevent.qmd", "exposure_ad_indexevent_corrected.csv") },
             format = "file", error = "continue", deployment = "main"),
  # Results + diagnostics report. Reads only results/pipeline/, so it re-renders whenever the
  # tables or figures change. Local only (nb_ targets are excluded on Great Lakes).
  tar_target(nb_report,
             { csv_results; csv_associations; csv_datasets; csv_step2_gap; csv_bsh_checks
               csv_apss_background; csv_instruments_classical; figures_all; csv_latent
               render_nb("robust_mr_report.qmd", "robust_mr_report.html") },
             format = "file", error = "continue", deployment = "main"),
  tar_target(nb_mediation,
             { bsh_files; render_nb("mr_mediation.qmd", "mr_mediation_results.csv") },
             format = "file", error = "continue", deployment = "main")
)

list(t_reference, t_datasets, t_selection, t_instruments, m_classical, m_apss,
     if (RUN_MRAID) m_mraid, m_cause, t_tables,
     t_notebooks)

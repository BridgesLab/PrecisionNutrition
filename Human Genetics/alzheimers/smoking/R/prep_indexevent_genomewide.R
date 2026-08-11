# prep_indexevent_genomewide.R — Step 1 of the properly-powered index-event correction.
#
# Builds a genome-wide, allele-aligned, LD-prunable merge of:
#   incidence/selection axis = parental lifespan (Pilling 2017, GCST006697, GRCh37)
#   prognosis                = AD (Bellenguez 2022, GCST90027158, GRCh38)
# Merged on rsID (build-agnostic). AD effects are aligned to the LIFESPAN effect allele, since
# SlopeHunter treats incidence as x. Palindromic and multi-allelic/indel SNPs are dropped
# (strand ambiguity). Writes a clump-ready table; APOE exclusion happens after clumping so the
# region can be reported both ways.
#
# Fixes the two defects of the original underpowered run (n=10 tophits):
#   (1) APOE (rs429358) dominated and is genuine AD-longevity pleiotropy, not a collider;
#   (2) a 15q25 SNP in the fitting set was circular (it is the locus under test).

suppressPackageStartupMessages({library(data.table); library(here)})

SUMSTAT_DIR <- here::here("data", "cache", "sumstats")
AD_FILE   <- file.path(SUMSTAT_DIR, "AD_bellenguez_harmonised.tsv.gz")
LIFE_FILE <- file.path(SUMSTAT_DIR, "LIFESPAN_pilling_harmonised.tsv.gz")

# Pre-filter the incidence trait: SlopeHunter fits on SNPs with incidence p < xp_thresh
# (default 1e-3), so keeping p < 0.01 is generous and makes the merge tractable.
LIFE_P_KEEP <- 0.01

is_palindromic <- function(a1, a2) {
  p <- paste0(a1, a2)
  p %in% c("AT", "TA", "CG", "GC")
}
is_snp <- function(a1, a2) nchar(a1) == 1L & nchar(a2) == 1L &
  a1 %in% c("A", "C", "G", "T") & a2 %in% c("A", "C", "G", "T")

prep_indexevent_data <- function() {
  cols <- c("variant_id", "chromosome", "base_pair_location", "effect_allele",
            "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")

  message("Reading lifespan (incidence) sumstats ...")
  life <- fread(cmd = paste("gunzip -c", shQuote(LIFE_FILE)), select = cols,
                showProgress = FALSE)
  setnames(life, c("SNP","chr_life","pos_life","ea_life","oa_life","eaf_life",
                   "beta_life","se_life","p_life"))
  life <- life[!is.na(p_life) & p_life < LIFE_P_KEEP &
                 !is.na(beta_life) & !is.na(se_life) & se_life > 0]
  life <- life[is_snp(ea_life, oa_life)]
  message("  lifespan SNPs with p < ", LIFE_P_KEEP, ": ", nrow(life))

  message("Reading AD (prognosis) sumstats ...")
  ad <- fread(cmd = paste("gunzip -c", shQuote(AD_FILE)), select = cols, showProgress = FALSE)
  setnames(ad, c("SNP","chr_ad","pos_ad","ea_ad","oa_ad","eaf_ad","beta_ad","se_ad","p_ad"))
  ad <- ad[!is.na(beta_ad) & !is.na(se_ad) & se_ad > 0]
  ad <- ad[is_snp(ea_ad, oa_ad)]
  message("  AD SNPs usable: ", nrow(ad))

  message("Merging on rsID ...")
  m <- merge(life, ad, by = "SNP")
  message("  merged: ", nrow(m))

  # Drop palindromic (strand-ambiguous) SNPs.
  m <- m[!is_palindromic(ea_life, oa_life)]

  # Align AD effect to the LIFESPAN effect allele; drop allele mismatches.
  same <- m$ea_ad == m$ea_life & m$oa_ad == m$oa_life
  flip <- m$ea_ad == m$oa_life & m$oa_ad == m$ea_life
  m <- m[same | flip]
  m[, beta_ad_aligned := fifelse(ea_ad == ea_life, beta_ad, -beta_ad)]
  message("  after allele alignment (palindromic/mismatch dropped): ", nrow(m))

  # Common variants only.
  m <- m[is.na(eaf_life) | (eaf_life > 0.01 & eaf_life < 0.99)]

  # APOE region flags in BOTH builds (AD file is GRCh38, lifespan GRCh37) — generous windows.
  m[, apoe := (chr_ad == 19 & pos_ad > 43.9e6 & pos_ad < 45.9e6) |
              (chr_life == 19 & pos_life > 44.4e6 & pos_life < 46.5e6)]
  # 15q25 CHRNA5/A3/B4 — the locus under test; including it in the fit is circular.
  m[, locus15q25 := (chr_life == 15 & pos_life > 78.5e6 & pos_life < 79.2e6)]

  message("  final: ", nrow(m), " | APOE-region: ", sum(m$apoe),
          " | 15q25: ", sum(m$locus15q25))
  m[]
}

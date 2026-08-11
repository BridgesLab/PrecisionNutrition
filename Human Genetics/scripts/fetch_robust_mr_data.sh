#!/usr/bin/env bash
# fetch_robust_mr_data.sh -------------------------------------------------------
# Download the genome-wide summary statistics and reference panels needed by the
# robust MR arm (CAUSE + MR-APSS). Unlike the existing OpenGWAS-driven scripts in
# this repo, both methods need FULL genome-wide statistics, not clumped instruments.
#
# Usage:  bash scripts/fetch_robust_mr_data.sh [target_dir]
# Default target_dir is ./raw_data
#
# Disk: ~12 GB downloaded, ~25 GB after decompression of the reference panel.
#
# NOTE ON URLS: upstream hosts (GEFOS, GLGC, Broad) have all reorganised in recent
# years. Every download below is attempted against a primary and one or more
# fallbacks, and the script reports what it could not get rather than failing
# silently. Anything it cannot fetch is listed at the end with a manual link.

set -uo pipefail

TARGET="${1:-raw_data}"
RAW="${TARGET}/robust_mr"
REF="${TARGET}/reference"
mkdir -p "${RAW}" "${REF}"

MISSING=()

# --retry-all-errors needs curl >= 7.71. macOS has shipped older ones; probe once
# rather than having every download die on an unknown-option error.
if curl --help all 2>/dev/null | grep -q -- '--retry-all-errors'; then
  RETRY_ALL="--retry-all-errors"
else
  RETRY_ALL=""
fi

# try_get <destination> <url> [url ...]
#
# Resumable: partial downloads are kept as <dest>.part and continued on the next
# run. Safe to Ctrl-C, close the laptop, or re-run after the wifi drops — nothing
# already fetched is re-fetched. A completed file is never re-downloaded.
try_get () {
  local dest="$1"; shift
  if [[ -s "${dest}" ]]; then
    echo "  [skip] $(basename "${dest}") already complete"
    return 0
  fi
  for url in "$@"; do
    # Plain string, not an array: macOS ships bash 3.2, where expanding an empty
    # array under `set -u` is an "unbound variable" error. Left deliberately
    # unquoted below so it word-splits into two arguments (or none).
    local resume=""
    if [[ -s "${dest}.part" ]]; then
      echo "  [cont] $(basename "${dest}") from $(du -h "${dest}.part" | cut -f1)"
      resume="--continue-at -"
    else
      echo "  [get ] ${url}"
    fi
    # Retries + a stall detector ride out home-wifi hiccups rather than failing
    # over to a fallback URL that is probably also fine.
    # shellcheck disable=SC2086
    if curl -fL ${resume} \
            --retry 10 --retry-delay 5 ${RETRY_ALL} \
            --connect-timeout 30 --speed-limit 1024 --speed-time 120 \
            --progress-bar -o "${dest}.part" "${url}"; then
      mv "${dest}.part" "${dest}"
      echo "  [ok  ] $(basename "${dest}")"
      return 0
    fi
    # Keep .part on transient failure so the next run resumes; only discard it
    # when moving on to a different URL, where a resume would corrupt the file.
    if [[ $# -gt 1 ]]; then rm -f "${dest}.part"; fi
  done
  echo "  [FAIL] $(basename "${dest}") — re-run this script to resume"
  MISSING+=("$(basename "${dest}")  <-  $1")
  return 1
}

echo "== 1. Exposure: LDL-C, GLGC 2021 (Graham et al., Nature) — EUR meta-analysis"
echo "   Hosted by the Willer lab at UMich."
GLGC="https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific"
try_get "${RAW}/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" \
  "${GLGC}/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"
try_get "${RAW}/TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" \
  "${GLGC}/TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"
echo "   If the filenames have changed, list the directory:"
echo "     curl -s ${GLGC}/ | grep -o 'LDL[^\"]*gz'"

echo
echo "== 2. Exposure (overlap-free arm): LDL-C, GLGC 2013 (Willer et al.), pre-UK Biobank"
if [[ -s "${TARGET}/jointGwasMc_LDL.txt.gz" ]]; then
  echo "  [skip] jointGwasMc_LDL.txt.gz already in raw_data/ — this is the pre-UKB arm"
else
  try_get "${TARGET}/jointGwasMc_LDL.txt.gz" \
    "https://csg.sph.umich.edu/willer/public/lipids2013/jointGwasMc_LDL.txt.gz"
fi

echo
echo "== 3. Outcome: heel eBMD, Morris 2019 (UK Biobank, N=426,824)"
echo "   GWAS Catalog accession GCST006979. Take the build37 '.f' file, NOT the"
echo "   '.h' harmonised one — .h is remapped to GRCh38 and everything else in"
echo "   this pipeline (GLGC POS_b37, the clumping panel, the MHC window) is b37."
EBI="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006979"
try_get "${RAW}/30598549-GCST006979-EFO_0009270-build37.f.tsv.gz" \
  "${EBI}/30598549-GCST006979-EFO_0009270-build37.f.tsv.gz"
echo "   Directory listing if the filename has changed:"
echo "     ${EBI}/"

echo
echo "== 4. LDSC reference: HapMap3 snplist + European LD scores"
echo "   Required by MR-APSS est_paras() for the C matrix (the overlap correction)."
echo "   The Broad's alkesgroup host retired these; Zenodo mirrors are primary now."
# HapMap3 snplist. This Zenodo record is the one the maintained LDSC tutorials
# point at; note it is .gz here, not the .bz2 the old Broad host served.
try_get "${REF}/w_hm3.snplist.gz" \
  "https://zenodo.org/records/7773502/files/w_hm3.snplist.gz?download=1"
[[ -s "${REF}/w_hm3.snplist.gz" ]] && gunzip -kf "${REF}/w_hm3.snplist.gz"

# LD scores. The original eur_w_ld_chr host is retired and the Zenodo mirrors use
# inconsistent filenames, so fall back to the S-LDSC 1000G Phase 3 EUR bundle,
# which is actively maintained. read_ldscores() in R accepts either layout, so it
# does not matter which of these you end up with.
try_get "${REF}/1000G_Phase3_ldscores.tgz" \
  "https://zenodo.org/records/10515792/files/1000G_Phase3_ldscores.tgz?download=1"
if [[ -s "${REF}/1000G_Phase3_ldscores.tgz" ]]; then
  mkdir -p "${REF}/eur_w_ld_chr"
  tar -xzf "${REF}/1000G_Phase3_ldscores.tgz" -C "${REF}/eur_w_ld_chr" --strip-components=1 || \
    tar -xzf "${REF}/1000G_Phase3_ldscores.tgz" -C "${REF}/eur_w_ld_chr" || true
fi
echo "   If that 404s, browse these and take whatever LD-score archive is listed:"
echo "     https://zenodo.org/records/10515792   (S-LDSC reference files, maintained)"
echo "     https://zenodo.org/records/18749273   (Feb 2026, 1000G Phase 3 EUR LD scores)"
echo "     https://zenodo.org/records/8182036    (Jul 2023, copy of the alkesgroup files)"
echo "   Unpack so that raw_data/reference/eur_w_ld_chr/ contains per-chromosome"
echo "   *.l2.ldscore.gz and *.l2.M_5_50 files. Either naming convention is fine:"
echo "     1.l2.ldscore.gz          (classic eur_w_ld_chr)"
echo "     LDscore.1.l2.ldscore.gz  (S-LDSC bundles)"

echo
echo "== 5. 1000 Genomes EUR plink panel (LD clumping for both methods)"
LOCAL_PANEL="${TARGET}/../alzheimers/data/cache/EUR.bed"
if [[ -s "${LOCAL_PANEL}" ]]; then
  echo "  [skip] already present at alzheimers/data/cache/EUR.{bed,bim,fam}"
  echo "         config_robust_mr.yml points at it — nothing to download."
else
  mkdir -p "${REF}/1kg_eur"
  try_get "${REF}/1kg_eur/EUR.tgz" "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"
  if [[ -s "${REF}/1kg_eur/EUR.tgz" ]]; then
    tar -xzf "${REF}/1kg_eur/EUR.tgz" -C "${REF}/1kg_eur" || true
  fi
fi

echo
echo "== 6. OPTIONAL — CAUSE-native LD blocks (1000G EUR, r2 and snpdata RDS)"
echo "   Only needed if you prefer cause::ld_prune() over plink clumping."
echo "   https://zenodo.org/record/1464357   (chr1-22 *_AF0.05_0.1.RDS + *_snpdata.RDS)"

echo
echo "=============================================================="
if [[ ${#MISSING[@]} -eq 0 ]]; then
  echo "All downloads succeeded."
else
  echo "Could not fetch ${#MISSING[@]} file(s) — get these manually:"
  printf '  %s\n' "${MISSING[@]}"
fi
echo
echo "Next: open robust_mr_prep.qmd and run the header-inspection chunk BEFORE"
echo "anything else. Correct the cols: blocks in config_robust_mr.yml to match."
echo "=============================================================="

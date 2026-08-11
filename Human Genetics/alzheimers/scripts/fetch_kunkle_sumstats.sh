#!/usr/bin/env bash
# fetch_kunkle_sumstats.sh — turnkey inputs for the Kunkle x lifespan SlopeHunter slope (b_SH).
#
# The Bellenguez b_SH is reused from the smoking pipeline, but Kunkle (clinical, no proxy) can
# carry a different selection slope and has no cached fit. This downloads + normalizes the three
# inputs R/fit_selection_slope.R needs, then prints the fit commands.
#
# Downloads (idempotent — existing files are skipped):
#   1. Kunkle 2019 IGAP Stage 1 AD   GCST007511  (prognosis trait)      ~543 MB raw .txt
#      -> normalized to the GWAS-Catalog-harmonised schema prep_selection_merge() expects.
#   2. Pilling 2017 parental lifespan GCST006697 (incidence/selection)  ~265 MB
#   3. 1000 Genomes EUR LD panel     (MRC-IEU mirror)                   ~1.5 GB -> ~2.6 GB extracted
#
# 2 and 3 are shared with smoking/scripts/fetch_sumstats.sh; if that cache exists they are
# symlinked rather than re-downloaded.
#
# DATA TERMS: the Kunkle/IGAP summary statistics are released for research use — by running this
# you accept the IGAP data terms (see the README fetched alongside the file). This pulls the
# publicly hosted EBI GWAS-Catalog copy; no login is required for that mirror.
#
# Usage:  bash scripts/fetch_kunkle_sumstats.sh
# Requires: curl, gzip, awk, tar (and plink2 for the clumping step that follows).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CACHE="$ROOT/data/cache"
SUMSTATS="$CACHE/sumstats"
SMOKING_CACHE="$ROOT/smoking/data/cache"          # to reuse shared lifespan + LD panel
mkdir -p "$SUMSTATS"

GWAS_FTP="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
KUNKLE_URL="$GWAS_FTP/GCST007001-GCST008000/GCST007511/Kunkle_etal_Stage1_results.txt"
KUNKLE_README="$GWAS_FTP/GCST007001-GCST008000/GCST007511/Kunkle_etal_2019_IGAP_summary_statistics_README_0.docx"
LIFE_URL="$GWAS_FTP/GCST006001-GCST007000/GCST006697/harmonised/29227965-GCST006697-EFO_0007796-build37.f.tsv.gz"
LD_URL="http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"

KUNKLE_RAW="$SUMSTATS/Kunkle_etal_Stage1_results.txt"
KUNKLE_OUT="$SUMSTATS/AD_kunkle_harmonised.tsv.gz"
LIFE_OUT="$SUMSTATS/LIFESPAN_pilling_harmonised.tsv.gz"
LD_TGZ="$CACHE/1kg.v3.tgz"

fetch () {  # fetch <url> <dest> <label> <min_bytes>
  local url="$1" dest="$2" label="$3" min="$4"
  if [[ -f "$dest" && $(stat -f%z "$dest" 2>/dev/null || stat -c%s "$dest") -ge "$min" ]]; then
    echo "  [skip] $label already present ($(du -h "$dest" | cut -f1))"
    return
  fi
  echo "  [get ] $label ..."
  curl -fSL --retry 3 --retry-delay 5 -o "$dest" "$url"
  echo "  [ok  ] $label -> $(du -h "$dest" | cut -f1)"
}

reuse_or_fetch () {  # reuse_or_fetch <shared_path> <dest> <url> <label> <min_bytes>
  local shared="$1" dest="$2" url="$3" label="$4" min="$5"
  if [[ -f "$dest" ]]; then echo "  [skip] $label already present"; return; fi
  if [[ -f "$shared" ]]; then
    echo "  [link] $label reused from smoking cache"; ln -s "$shared" "$dest"; return
  fi
  fetch "$url" "$dest" "$label" "$min"
}

echo "Fetching Kunkle x lifespan b_SH inputs into $CACHE"
echo

echo "1/3 Kunkle 2019 IGAP Stage 1 AD summary statistics (GRCh37)"
if [[ -f "$KUNKLE_OUT" && $(stat -f%z "$KUNKLE_OUT" 2>/dev/null || stat -c%s "$KUNKLE_OUT") -ge 50000000 ]]; then
  echo "  [skip] AD_kunkle_harmonised.tsv.gz already present ($(du -h "$KUNKLE_OUT" | cut -f1))"
else
  fetch "$KUNKLE_URL"    "$KUNKLE_RAW"                                   "Kunkle_etal_Stage1_results.txt" 400000000
  fetch "$KUNKLE_README" "$SUMSTATS/Kunkle_README.docx"                 "Kunkle README (data terms)"     10000
  echo "  [norm] mapping IGAP columns -> harmonised schema (adds NA effect_allele_frequency) ..."
  # Header-name-driven (order/delimiter agnostic): Chromosome, Position, MarkerName(rsID),
  # Effect_allele, Non_Effect_allele, Beta, SE, Pvalue. AD-side EAF is not used downstream.
  awk '
    NR==1 {
      for (i=1;i<=NF;i++) { k=tolower($i); gsub(/[^a-z0-9]/,"",k); col[k]=i }
      c_chr=col["chromosome"]; c_pos=col["position"]; c_snp=col["markername"];
      c_ea=col["effectallele"]; c_oa=col["noneffectallele"];
      c_b=col["beta"]; c_se=col["se"]; c_p=col["pvalue"];
      if (!(c_chr&&c_pos&&c_snp&&c_ea&&c_oa&&c_b&&c_se&&c_p)) {
        print "ERROR: a required Kunkle column was not found in the header" > "/dev/stderr"; exit 1 }
      print "variant_id\tchromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tbeta\tstandard_error\tp_value";
      next
    }
    { print $c_snp"\t"$c_chr"\t"$c_pos"\t"$c_ea"\t"$c_oa"\tNA\t"$c_b"\t"$c_se"\t"$c_p }
  ' "$KUNKLE_RAW" | gzip -c > "$KUNKLE_OUT"
  echo "  [ok  ] AD_kunkle_harmonised.tsv.gz -> $(du -h "$KUNKLE_OUT" | cut -f1)"
  rm -f "$KUNKLE_RAW"   # drop the 543 MB raw; the normalized gz is the input
fi

echo "2/3 Pilling 2017 parental-longevity summary statistics (Martingale scale; +beta = shorter life)"
reuse_or_fetch "$SMOKING_CACHE/sumstats/LIFESPAN_pilling_harmonised.tsv.gz" "$LIFE_OUT" \
               "$LIFE_URL" "LIFESPAN_pilling_harmonised.tsv.gz" 200000000

echo "3/3 1000 Genomes EUR LD reference panel (for plink2 --clump)"
if [[ -f "$CACHE/EUR.bed" ]]; then
  echo "  [skip] EUR panel already present"
elif [[ -f "$SMOKING_CACHE/EUR.bed" ]]; then
  echo "  [link] EUR panel reused from smoking cache"
  for ext in bed bim fam; do ln -sf "$SMOKING_CACHE/EUR.$ext" "$CACHE/EUR.$ext"; done
else
  fetch "$LD_URL" "$LD_TGZ" "1kg.v3.tgz" 1400000000
  echo "  [tar ] extracting EUR.{bed,bim,fam} ..."
  tar -xzf "$LD_TGZ" -C "$CACHE"
fi

cat <<'NEXT'

Done. Next steps (produces R/fits/kunkle_lifespan_slopehunter_fits.csv, consumed by
exposure_ad_indexevent.qmd):

  # 1. merge + allele-align, write the plink clump input (ID = rsID, P = lifespan p)
  Rscript -e 'source("R/fit_selection_slope.R");
              m <- prep_selection_merge("data/cache/sumstats/AD_kunkle_harmonised.tsv.gz",
                                        "data/cache/sumstats/LIFESPAN_pilling_harmonised.tsv.gz");
              saveRDS(m, "data/cache/sumstats/merged_kunkle_lifespan.rds");
              data.table::fwrite(m[, .(ID = SNP, P = p_life)],
                                 "data/cache/sumstats/kunkle_clump_input.tsv", sep = "\t")'

  # 2. LD-clump to an independent SNP set
  cd data/cache && plink2 --bfile EUR --clump sumstats/kunkle_clump_input.tsv \
       --clump-p1 1e-3 --clump-p2 0.01 --clump-r2 0.01 --clump-kb 1000 \
       --out sumstats/kunkle_clumped && cd ../..

  # 3. fit SlopeHunter on the clumped set -> fits CSV (one call: preps, filters, fits, writes)
  Rscript -e 'source("R/fit_selection_slope.R");
              dir.create("R/fits", showWarnings = FALSE, recursive = TRUE);
              fits <- build_selection_slope(
                ad_file      = "data/cache/sumstats/AD_kunkle_harmonised.tsv.gz",
                life_file    = "data/cache/sumstats/LIFESPAN_pilling_harmonised.tsv.gz",
                clumped_snps = "data/cache/sumstats/kunkle_clumped.clumps",
                out_csv      = "R/fits/kunkle_lifespan_slopehunter_fits.csv");
              print(fits)'

Then re-render exposure_ad_indexevent.qmd — the Kunkle rows in the correction table will populate.
NEXT

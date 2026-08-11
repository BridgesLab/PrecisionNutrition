#!/usr/bin/env bash
# fetch_sumstats.sh — download the large external inputs required by Layer 4e
# (the genome-wide index-event / SlopeHunter correction).
#
# These files are NOT pullable from the OpenGWAS API (which only serves clumped tophits and
# regional queries). The underpowered n=10 version of Task D was a direct consequence of that
# limitation; the genome-wide redo needs full summary statistics plus an LD reference panel.
#
# Downloads (~3.5 GB total, idempotent — existing files are skipped):
#   1. Bellenguez 2022 AD          GCST90027158  (prognosis trait)      ~575 MB
#   2. Pilling 2017 parental lifespan GCST006697 (incidence/selection)  ~265 MB
#   3. 1000 Genomes EUR LD panel   (MRC-IEU mirror)                     ~1.5 GB -> ~2.6 GB extracted
#
# Consumed by:
#   R/prep_indexevent_genomewide.R   (merge + allele-align -> merged_AD_lifespan.rds)
#   plink2 --clump                   (LD-independent SNP set)
#   R/run_indexevent_genomewide.R    (SlopeHunter fit + correction)
#
# Usage:  bash scripts/fetch_sumstats.sh
# Requires: curl, tar, and (for the clumping step that follows) plink2.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CACHE="$ROOT/data/cache"
SUMSTATS="$CACHE/sumstats"
mkdir -p "$SUMSTATS"

GWAS_FTP="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
AD_URL="$GWAS_FTP/GCST90027001-GCST90028000/GCST90027158/harmonised/35379992-GCST90027158-MONDO_0004975-Build38.f.tsv.gz"
LIFE_URL="$GWAS_FTP/GCST006001-GCST007000/GCST006697/harmonised/29227965-GCST006697-EFO_0007796-build37.f.tsv.gz"
LD_URL="http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"

AD_OUT="$SUMSTATS/AD_bellenguez_harmonised.tsv.gz"
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

echo "Fetching Layer 4e inputs into $CACHE"
echo
echo "1/3 Bellenguez 2022 AD summary statistics (GRCh38, harmonised)"
fetch "$AD_URL" "$AD_OUT" "AD_bellenguez_harmonised.tsv.gz" 400000000

echo "2/3 Pilling 2017 parental-longevity summary statistics (harmonised)"
echo "     NB: Martingale-residual scale — POSITIVE beta = SHORTER lifespan."
echo "     (verified against ApoE4 rs429358 beta=+0.057 and CHRNA5 rs16969968 beta=+0.025)"
fetch "$LIFE_URL" "$LIFE_OUT" "LIFESPAN_pilling_harmonised.tsv.gz" 200000000

echo "3/3 1000 Genomes EUR LD reference panel (for plink2 --clump)"
fetch "$LD_URL" "$LD_TGZ" "1kg.v3.tgz" 1400000000
if [[ ! -f "$CACHE/EUR.bed" ]]; then
  echo "  [tar ] extracting EUR.{bed,bim,fam} ..."
  tar -xzf "$LD_TGZ" -C "$CACHE"
fi
echo "  [ok  ] EUR panel: $(ls "$CACHE"/EUR.* 2>/dev/null | wc -l | tr -d ' ') files"

cat <<'NEXT'

Done. Next steps (Layer 4e):

  # 1. merge + allele-align the two sumstats (writes merged_AD_lifespan.rds + clump_input.tsv)
  Rscript -e 'source("R/prep_indexevent_genomewide.R"); m <- prep_indexevent_data();
              saveRDS(m, "data/cache/sumstats/merged_AD_lifespan.rds");
              data.table::fwrite(m[, .(ID = SNP, P = p_life)],
                                 "data/cache/sumstats/clump_input.tsv", sep = "\t")'

  # 2. LD-clump to an independent SNP set (~2,346 SNPs)
  cd data/cache && plink2 --bfile EUR --clump sumstats/clump_input.tsv \
       --clump-p1 1e-3 --clump-p2 0.01 --clump-r2 0.01 --clump-kb 1000 \
       --out sumstats/clumped

  # 3. fit SlopeHunter + apply the correction (writes results/indexevent_*.csv)
  Rscript -e 'source("R/run_indexevent_genomewide.R");
              idx <- data.table::fread("data/cache/sumstats/clumped.clumps")$ID;
              r <- run_indexevent_genomewide(idx);
              readr::write_csv(r$fits, "results/indexevent_slopehunter_fits.csv");
              readr::write_csv(r$corrected, "results/indexevent_correction_genomewide.csv")'

Then render 04e_overlap_indexevent.qmd (its live chunks read the results/ CSVs and need no network).
NEXT

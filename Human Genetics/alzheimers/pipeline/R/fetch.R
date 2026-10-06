# fetch.R — download full GWAS summary statistics from OpenGWAS.
#
# Full sumstats come from the file-download endpoint (ieugwasr::gwasinfo_files()), which returns
# short-lived signed URLs. This is separate from the /tophits and /ld/clump endpoints, so it keeps
# working when those are down (as they were from 2026-10-01). OpenGWAS limits downloads to
# 20 datasets per 24 h per account.

# Download <id>.vcf.gz (+ .tbi and the OpenGWAS QC report) into `dir`. Skips a file that is already
# present and passes `gzip -t`, so re-running is cheap and an interrupted download is redone.
# Returns the path of the .vcf.gz (for a targets format = "file" target).
fetch_opengwas_vcf <- function(id, dir = "data/cache/gwas") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  vcf <- file.path(dir, paste0(id, ".vcf.gz"))
  if (file.exists(vcf) && gz_ok(vcf)) return(vcf)

  urls <- ieugwasr::gwasinfo_files(id)[[1]]
  if (!length(urls)) stop("gwasinfo_files() returned no files for ", id)
  name_of <- function(u) basename(sub("[?].*$", "", u))
  old <- options(timeout = 7200); on.exit(options(old))
  for (u in urls) {
    dest <- file.path(dir, name_of(u))
    tmp <- paste0(dest, ".part")
    utils::download.file(u, tmp, mode = "wb", quiet = TRUE, method = "libcurl")
    file.rename(tmp, dest)
  }
  if (!gz_ok(vcf)) stop("downloaded ", vcf, " is not a valid gzip file")
  vcf
}

gz_ok <- function(path) system2("gzip", c("-t", shQuote(path)), stdout = FALSE, stderr = FALSE) == 0

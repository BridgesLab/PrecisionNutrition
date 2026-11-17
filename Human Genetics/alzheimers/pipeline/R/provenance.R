# provenance.R — a stamp beside every output file, recording what produced it.
#
# stamp_output("results/pipeline/results.csv", step) writes results/pipeline/stamps/results.csv.txt:
# the file's md5, the step that wrote it, when and where, the git commit (if the folder is a git
# checkout; on Great Lakes it is usually a copy), a digest of the pipeline code, and the versions of
# the packages that produce the estimates. One stamp per file means no shared file for parallel
# workers to race on, and building a single target by name still stamps what it wrote.
# History across runs is kept by hand in results/pipeline/RUN_INFO.md.
#
# To check whether a Great Lakes copy matched a commit, run pipeline_code_digest() in a checkout of
# that commit and compare it with the stamp's code_digest.

STAMP_PACKAGES <- c("TwoSampleMR", "MRAPSS", "cause", "MRBEE", "mr.raps", "MRPRESSO", "MRAID",
                    "SlopeHunter", "targets", "data.table", "dplyr")

# md5 of the code that defines the analysis, combined into one digest. Paths relative to the
# project folder (the one holding _targets.R).
pipeline_code_digest <- function() {
  files <- sort(c("_targets.R", "pipeline/config.yml", list.files("pipeline/R", "\\.R$", full.names = TRUE),
                  list.files("R", "\\.R$", full.names = TRUE)))
  files <- files[file.exists(files)]
  sums <- unname(tools::md5sum(files))
  tmp <- tempfile(); on.exit(unlink(tmp))
  writeLines(paste(sums, files), tmp)
  unname(tools::md5sum(tmp))
}

stamp_git <- function(...) {
  if (!nzchar(Sys.which("git"))) return(character())
  suppressWarnings(tryCatch(system2("git", c(...), stdout = TRUE, stderr = FALSE), error = \(e) character()))
}

stamp_output <- function(path, step) {
  dir <- file.path(dirname(path), "stamps")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  commit <- stamp_git("rev-parse", "--short", "HEAD")
  dirty  <- stamp_git("status", "--porcelain", "--", "_targets.R", "pipeline", "R")
  pkgs   <- vapply(STAMP_PACKAGES, \(p) tryCatch(as.character(utils::packageVersion(p)), error = \(e) "not installed"), "")
  writeLines(c(
    paste("file:", path),
    paste("md5:", unname(tools::md5sum(path))),
    paste("step:", step),
    paste("written:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("host:", Sys.info()[["nodename"]]),
    paste("git_commit:", if (length(commit)) commit else "not a git checkout"),
    paste("uncommitted_code:", if (!length(commit)) "unknown" else if (length(dirty)) paste(trimws(dirty), collapse = "; ") else "none"),
    paste("code_digest:", pipeline_code_digest()),
    paste("R:", R.version$version.string),
    paste("packages:", paste(names(pkgs), pkgs, sep = " ", collapse = "; "))),
    file.path(dir, paste0(basename(path), ".txt")))
  invisible(path)
}

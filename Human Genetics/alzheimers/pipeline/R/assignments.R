# assignments.R — outcome-blind groupings of the T2D instruments, frozen before AD is merged.
#
# Two groupings, both built from exposure-side information only:
#   clusters      Suzuki 2024 (Nature 627:347) hard clusters of the 1,289 T2D index variants
#                 (8 clusters from 37 cardiometabolic traits; no neurological trait)
#   drug targets  genes that are the mechanism-of-action target of a drug in phase 2+ or approved for T2D
#                 (Open Targets Platform, ChEMBL-derived), grouped into drug classes by the regexes in
#                 config `contrast: drug_targets: classes`
# write_frozen() writes each table to assignments/ and records its checksum in
# assignments/MANIFEST.csv; assert_frozen() is the gate later steps call before touching an outcome.

suppressPackageStartupMessages({ library(data.table); library(dplyr) })

# ---- chr:pos summary statistics (T2DGGI) ----------------------------------------------------------

# T2DGGI files have chromosome/position (hg19) and no rsID, so variants are named from the LD panel's
# bim by chr:pos + alleles. Column names are matched by pattern because the README describes the
# columns without naming them; an unmatched column stops with the header shown.
read_gwas_chrpos <- function(path, bim, columns = NULL) {
  cmd <- if (grepl("\\.zip$", path)) paste("unzip -p", shQuote(path)) else paste("gunzip -c", shQuote(path))
  d <- fread(cmd = cmd, showProgress = FALSE)
  pick <- function(what, pattern) {
    if (!is.null(columns[[what]])) return(columns[[what]])       # config `columns:` overrides the pattern
    hit <- grep(pattern, names(d), ignore.case = TRUE, value = TRUE)
    if (length(hit) != 1) stop(sprintf("%s: %d columns match %s ('%s'); set `columns: {%s: ...}` in config. Header: %s",
                                       basename(path), length(hit), what, pattern, what, paste(names(d), collapse = ", ")))
    hit
  }
  cols <- c(chr = pick("chr", "^chr"), pos = pick("pos", "^(pos|bp|base)"),
            ea = pick("ea", "^(effect.?allele$|ea$|a1$)"), oa = pick("oa", "^(other.?allele$|non.?effect.?allele$|oa$|a2$)"),
            eaf = pick("eaf", "(eaf|freq)"), beta = pick("beta", "beta"),
            se = pick("se", "(^se$|stderr|standard)"), p = pick("p", "(^p$|pval|p.value)"),
            n = pick("neff", "(^neff|effective)"))
  d <- d[, unname(cols), with = FALSE]; setnames(d, names(cols))
  n_total <- nrow(d)
  d[, `:=`(chr = sub("^chr", "", as.character(chr)), pos = as.integer(pos), ea = toupper(ea), oa = toupper(oa))]
  b <- copy(bim)[, chr := as.character(chr)]
  d <- merge(d, b[, .(chr, pos, SNP, a1, a2)], by = c("chr", "pos"))
  d <- d[(ea == a1 & oa == a2) | (ea == a2 & oa == a1)]
  d <- d[, .(SNP, chr, pos, ea, oa, eaf = as.numeric(eaf), beta = as.numeric(beta), se = as.numeric(se),
             p = suppressWarnings(as.numeric(p)), n = as.numeric(n))]
  setattr(d, "n_total", n_total)
  d[]
}

load_gwas_chrpos <- function(path, ds, bim, qc_cfg)
  qc_sumstats(read_gwas_chrpos(path, bim, ds$columns), ds$key, bim$SNP, maf_min = qc_cfg$maf_min,
              n_fixed = ds$n %||% NA_real_)

# ---- Suzuki 2024 clusters ---------------------------------------------------------------------------

# Index variants with their cluster, from the paper's Supplementary Tables workbook. The sheet is set
# in config; columns are found by pattern (rsID, or chr + position mapped through the bim). Returns
# SNP, cluster, and the paper's locus label where there is one.
read_suzuki_clusters <- function(path, sheet, bim, skip = 0) {
  x <- as.data.table(readxl::read_excel(path, sheet = sheet, skip = skip, .name_repair = "minimal"))
  nm <- names(x)
  one <- function(pattern) { h <- grep(pattern, nm, ignore.case = TRUE); if (length(h)) h[1] else NA }
  i_cl <- one("cluster")
  if (is.na(i_cl)) stop("no 'cluster' column in sheet ", sheet, "; columns: ", paste(nm, collapse = ", "))
  i_rs <- which(vapply(x, \(v) mean(grepl("^rs[0-9]+$", as.character(v))) > 0.8, logical(1)))[1]
  out <- data.table(cluster = as.character(x[[i_cl]]),
                    locus = if (!is.na(one("locus"))) as.character(x[[one("locus")]]) else NA_character_)
  if (!is.na(i_rs)) {
    out[, SNP := as.character(x[[i_rs]])]
  } else {
    i_chr <- one("^chr"); i_pos <- one("pos")
    if (is.na(i_chr) || is.na(i_pos)) stop("no rsID or chr/position columns in sheet ", sheet,
                                           "; columns: ", paste(nm, collapse = ", "))
    out[, `:=`(chr = sub("^chr", "", as.character(x[[i_chr]])), pos = as.integer(x[[i_pos]]))]
    out <- merge(out, copy(bim)[, .(chr = as.character(chr), pos, SNP)], by = c("chr", "pos"), all.x = TRUE)
  }
  out <- out[!is.na(cluster) & nzchar(cluster)]
  message(sprintf("Suzuki clusters: %d rows, %d with an rsID in the LD panel, %d clusters",
                  nrow(out), sum(out$SNP %chin% bim$SNP), uniqueN(out$cluster)))
  as_tibble(out[, .(SNP, cluster, locus)])
}

# ---- drug targets -----------------------------------------------------------------------------------

ot_query <- function(query) {
  r <- httr2::request("https://api.platform.opentargets.org/api/v4/graphql") |>
    httr2::req_body_json(list(query = query)) |> httr2::req_timeout(180) |> httr2::req_perform() |>
    httr2::resp_body_json()
  if (!is.null(r$errors)) stop("Open Targets: ", r$errors[[1]]$message)
  r$data
}

# Every drug Open Targets lists for the disease, one row per drug x mechanism x target gene, with the
# platform's data release. This is a snapshot: the target is built once (cue "never") and the table
# is written to assignments/, so a later API or release change cannot move the grouping.
fetch_drug_targets <- function(disease_id) {
  d <- ot_query(sprintf('{ meta { dataVersion { year month } }
    disease(efoId: "%s") { drugAndClinicalCandidates { rows { maxClinicalStage
      drug { id name drugType mechanismsOfAction { rows { mechanismOfAction actionType
        targets { id approvedSymbol } } } } } } } }', disease_id))
  ver <- paste(d$meta$dataVersion$year, d$meta$dataVersion$month, sep = ".")
  purrr::map_dfr(d$disease$drugAndClinicalCandidates$rows, \(r)
    purrr::map_dfr(r$drug$mechanismsOfAction$rows, \(m)
      purrr::map_dfr(m$targets, \(t) tibble(
        drug_id = r$drug$id, drug = r$drug$name, drug_type = r$drug$drugType,
        stage_for_disease = r$maxClinicalStage, mechanism = m$mechanismOfAction,
        action = m$actionType %||% NA_character_, gene = t$approvedSymbol, ensembl_id = t$id)))) |>
    mutate(disease_id = disease_id, opentargets_release = ver)
}

# Mechanisms at the configured clinical stages -> drug class (first matching regex in `classes`; NA = not a
# glucose-lowering class or too unspecific to instrument, kept in the table with class NA).
classify_drug_targets <- function(drugs, cfg) {
  cls <- cfg$classes
  drugs |> filter(stage_for_disease %in% cfg$stages) |>
    mutate(drug_class = purrr::map_chr(mechanism, \(m) {
      hit <- names(cls)[purrr::map_lgl(cls, \(rx) grepl(rx, m, ignore.case = TRUE))]
      if (length(hit)) hit[1] else NA_character_ })) |>
    group_by(gene, ensembl_id, drug_class, mechanism, action) |>
    summarise(n_drugs = n_distinct(drug_id), stages = paste(sort(unique(stage_for_disease)), collapse = "; "),
              drugs = paste(sort(unique(drug)), collapse = "; "),
              opentargets_release = first(opentargets_release), .groups = "drop")
}

# GRCh37 gene coordinates (the build of every GWAS here) from Ensembl's GRCh37 REST server.
gene_coords_b37 <- function(symbols) {
  r <- httr2::request("https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens") |>
    httr2::req_headers(Accept = "application/json") |>
    httr2::req_body_json(list(symbols = as.list(unique(symbols)))) |> httr2::req_timeout(120) |>
    httr2::req_perform() |> httr2::resp_body_json()
  purrr::imap_dfr(r, \(g, s) tibble(gene = s, chr = g$seq_region_name, start = g$start, end = g$end))
}

# Instrument -> drug-target gene when the variant lies within `window_kb` of the gene body. A variant
# near two targets takes the nearer one (ABCC8/KCNJ11 are adjacent and share a class).
assign_drug_targets <- function(inst, targets, coords, window_kb) {
  tg <- targets |> filter(!is.na(drug_class)) |> distinct(gene, drug_class) |>
    inner_join(coords |> filter(chr %in% c(1:22, "X")), by = "gene")
  w <- window_kb * 1000
  purrr::pmap_dfr(inst |> distinct(SNP, chr, pos), \(SNP, chr, pos) {
    near <- tg[tg$chr == as.character(chr) & pos >= tg$start - w & pos <= tg$end + w, ]
    if (!nrow(near)) return(NULL)
    near |> mutate(SNP = SNP, distance_kb = pmax(start - pos, pos - end, 0) / 1000) |>
      slice_min(distance_kb, n = 1, with_ties = FALSE) |>
      select(SNP, drug_target_gene = gene, drug_class, distance_kb)
  })
}

# ---- drug-target (cis) instruments ---------------------------------------------------------------------

# Drug-target MR instruments, as in the lipid-target literature: variants in or near each target gene
# (gene body +- window_kb) associated with the selection trait (T2D) at a RELAXED threshold, lightly
# clumped. Clumping is done once per drug class, so adjacent targets of one class (ABCC8/KCNJ11) share
# one set and no variant is counted twice; each kept variant is labelled with its nearest target gene.
# Selection is on T2D and the FG / 2hGlu effects come from MAGIC, as for the genome-wide sets.
select_cis <- function(g_t2d, targets, coords, cfg, bfile, plink2) {
  tg <- targets |> filter(!is.na(drug_class)) |> distinct(gene, drug_class) |>
    inner_join(coords |> filter(chr %in% as.character(1:22)), by = "gene")
  w <- cfg$window_kb * 1000
  purrr::map_dfr(split(tg, tg$drug_class), \(cl) {
    cand <- rbindlist(purrr::pmap(cl[c("chr", "start", "end")], \(chr, start, end) {
      g_chr <- chr; lo <- start - w; hi <- end + w
      g_t2d[chr == g_chr & pos >= lo & pos <= hi & p < cfg$p] }))
    cand <- unique(cand, by = "SNP")
    if (!nrow(cand)) return(NULL)
    idx <- plink_clump(cand, bfile, cfg$p, cfg$r2, cfg$clump_kb, plink2)
    k <- cand[SNP %chin% idx]
    near <- purrr::map2_chr(k$chr, k$pos, \(ch, bp) {
      c2 <- cl[cl$chr == ch, ]; c2$gene[which.min(pmax(c2$start - bp, bp - c2$end, 0))] })
    as_tibble(k) |> mutate(drug_class = cl$drug_class[1], drug_target_gene = near, .before = 1)
  })
}

# What the target calls: cis instruments -> FG and 2hGlu effects on the T2D risk allele, plus a per-gene
# strength summary (mean z^2 on each trait). Outcome-blind.
run_cis_step0 <- function(g_t2d, g_fg, g_2h, set, targets, coords, cfg, bfile, plink2) {
  sel <- select_cis(g_t2d, targets, coords, cfg, bfile, plink2)
  if (!nrow(sel)) stop(set, ": no cis instruments at p < ", cfg$p)
  h <- contrast_harmonise(as.data.table(sel), g_fg, g_2h, paste0("cis_", set))
  d <- h$dat |> inner_join(sel |> select(SNP, drug_class, drug_target_gene), by = "SNP") |>
    mutate(set = set, .before = 1)
  summ <- d |> group_by(set, drug_class, drug_target_gene) |>
    summarise(n_snp = n(), min_p_t2d = min(p_t2d), mean_F_t2d = mean((beta_t2d / se_t2d)^2),
              mean_F_fg = mean((beta_fg / se_fg)^2), mean_F_2h = mean((beta_2h / se_2h)^2),
              n_nominal_fg = sum(p_fg < 0.05), n_nominal_2h = sum(p_2h < 0.05), .groups = "drop") |>
    full_join(targets |> filter(!is.na(drug_class)) |> distinct(drug_class, drug_target_gene = gene) |>
                mutate(set = set), by = c("set", "drug_class", "drug_target_gene")) |>
    mutate(n_selected = purrr::map_int(drug_target_gene, \(g) sum(sel$drug_target_gene == g)),
           n_snp = coalesce(n_snp, 0L)) |>
    arrange(drug_class, drug_target_gene)
  list(instruments = d, summary = summ, attrition = h$attrition)
}

# One row per instrument: cluster and drug-target group (NA = unassigned; the unassigned fraction is
# reported, not hidden).
build_assignments <- function(inst, clusters, drug_map) {
  out <- inst |> select(SNP, chr, pos, ea, oa)
  out <- if (is.null(clusters)) mutate(out, cluster = NA_character_, locus = NA_character_)
         else left_join(out, distinct(clusters, SNP, .keep_all = TRUE), by = "SNP")
  left_join(out, drug_map, by = "SNP")
}

# ---- freeze -------------------------------------------------------------------------------------------

git_head <- function() tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)[1],
                                error = function(e) NA_character_)

# Write a table to assignments/<name>.csv and record it in assignments/MANIFEST.csv. `frozen_at` only
# changes when the content does, so re-running the pipeline does not re-stamp an unchanged file.
write_frozen <- function(x, name, dir = "assignments") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(dir, paste0(name, ".csv")); mf <- file.path(dir, "MANIFEST.csv")
  readr::write_csv(x, f)
  row <- tibble(file = basename(f), md5 = unname(tools::md5sum(f)), n_rows = nrow(x),
                frozen_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), git_head_when_written = git_head())
  man <- if (file.exists(mf)) readr::read_csv(mf, col_types = "ccicc") else row[0, ]
  old <- man[man$file == row$file, ]
  if (nrow(old) && identical(old$md5, row$md5)) row <- old
  readr::write_csv(bind_rows(man[man$file != row$file, ], row) |> arrange(file), mf)
  c(f, mf)
}

# Gate for every step that reads an outcome: the files must match the manifest AND be committed, so
# the git history timestamps the freeze. Returns the manifest.
assert_frozen <- function(dir = "assignments") {
  man <- readr::read_csv(file.path(dir, "MANIFEST.csv"), col_types = "ccicc")
  now <- unname(tools::md5sum(file.path(dir, man$file)))
  if (any(is.na(now)) || any(now != man$md5))
    stop("assignments/ does not match MANIFEST.csv: ", paste(man$file[is.na(now) | now != man$md5], collapse = ", "))
  dirty <- system2("git", c("status", "--porcelain", "--", shQuote(dir)), stdout = TRUE)
  if (length(dirty)) stop("assignments/ has uncommitted changes; commit the freeze before unblinding:\n",
                          paste(dirty, collapse = "\n"))
  man
}

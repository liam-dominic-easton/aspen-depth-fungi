## helpers.R — Shared functions for the depth-paper analysis scripts.
##
## Loads phyloseq objects, environmental data, and the MaarjAM virtual-taxon
## (VT) mapping; provides distance / PCoA / PERMANOVA helpers used across the
## analysis scripts in this package.
##
## All paths are resolved relative to PACKAGE_ROOT, the repository root (the
## directory that contains the `scripts/`, `data/`, and `outputs/` folders).
## Run any script from PACKAGE_ROOT as the working directory.

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(readr)
  library(vegan)
})

## ---- locate package root ----------------------------------------------------
find_package_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "scripts", "helpers.R")) &&
        file.exists(file.path(current, "data", "phyloseq"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate the repository root. Run scripts from the folder that contains scripts/ and data/.")
}

PACKAGE_ROOT <- find_package_root()
DATA_DIR     <- file.path(PACKAGE_ROOT, "data")
PHY_DIR      <- file.path(DATA_DIR, "phyloseq")
ENV_DIR      <- file.path(DATA_DIR, "environmental")
TAX_DIR      <- file.path(DATA_DIR, "taxonomy")
PRECOMP_DIR  <- file.path(DATA_DIR, "precomputed")
## OUTPUT_TAG routes a sensitivity run to its own folders, e.g.
## OUTPUT_TAG=no_MOO3B -> outputs/tables_no_MOO3B, outputs/figures_no_MOO3B
OUTPUT_TAG   <- Sys.getenv("OUTPUT_TAG", "")
.tag         <- if (nzchar(OUTPUT_TAG)) paste0("_", OUTPUT_TAG) else ""
TABLES_DIR   <- file.path(PACKAGE_ROOT, "outputs", paste0("tables", .tag))
FIGS_DIR     <- file.path(PACKAGE_ROOT, "outputs", paste0("figures", .tag))
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGS_DIR,   showWarnings = FALSE, recursive = TRUE)

## ---- sample exclusion switch -------------------------------------------------
## Comma-separated sample IDs to drop everywhere, e.g. EXCLUDE_SAMPLES=MOO3B
## for a sensitivity run without that sample.
EXCLUDE_SAMPLES <- setdiff(trimws(strsplit(Sys.getenv("EXCLUDE_SAMPLES", ""), ",")[[1]]), "")

## Number of permutations for all permutation tests (9,999 by default, so that
## BH-adjusted p-values are not bounded by permutation resolution).
NPERM <- as.integer(Sys.getenv("NPERM", "9999"))

## ---- biome assignments per node ---------------------------------------------
BOREAL_NODES   <- c("HIG", "NOT", "RED", "CAL", "GLA")
PARKLAND_NODES <- c("EDG", "DUN", "YRK", "MOO")

## ---- phyloseq sample_data → data.frame --------------------------------------
ps_sample_df <- function(ps) {
  df <- data.frame(sample_data(ps), stringsAsFactors = FALSE, check.names = FALSE)
  for (col in names(df)) {
    if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
  }
  df
}

## ---- load a phyloseq library, restricted to the clean sample set ------------
##
## Each library has two phyloseq files:
##   ps_<lib>_all.rds     — the full set of samples after DADA2
##   ps_<lib>_clean.rds   — the QC-passing subset used in the manuscript
## load_ps() returns the intersection: full taxonomy & abundance tables, but
## limited to the 52-sample analysis set.
load_ps <- function(library_name = c("ITS", "18S")) {
  library_name <- match.arg(library_name)
  files <- switch(library_name,
    "ITS" = list(all = "ps_fungi_all.rds", clean = "ps_fungi_clean.rds"),
    "18S" = list(all = "ps_amf_all.rds", clean = "ps_amf_euk_clean.rds")
  )
  ps_all   <- readRDS(file.path(PHY_DIR, files$all))
  ps_clean <- readRDS(file.path(PHY_DIR, files$clean))
  ps <- prune_samples(sample_names(ps_all) %in% sample_names(ps_clean), ps_all)
  if (length(EXCLUDE_SAMPLES)) {
    key <- get_sample_key(ps)
    drop <- key$sample_name[key$sample_id_clean %in% EXCLUDE_SAMPLES]
    ps <- prune_samples(!sample_names(ps) %in% drop, ps)
    ps <- prune_taxa(taxa_sums(ps) > 0, ps)
    cat("  excluded:", paste(EXCLUDE_SAMPLES, collapse = ", "), "\n")
  }
  cat(library_name, ":", phyloseq::nsamples(ps), "samples loaded\n")
  ps
}

## ---- get sample-name → site-id key from phyloseq metadata -------------------
get_sample_key <- function(ps) {
  sd <- ps_sample_df(ps) |>
    tibble::rownames_to_column("sample_name") |>
    mutate(
      node            = toupper(as.character(node)),
      stand           = as.character(stand),
      depth           = toupper(as.character(depth)),
      sample_id_clean = paste0(node, stand, depth)
    )
  sd[, c("sample_name", "sample_id_clean", "node", "stand", "depth")]
}

## ---- match phyloseq samples to environmental data ---------------------------
match_ps_to_env <- function(ps, env) {
  key <- get_sample_key(ps)
  inner_join(key, env, by = "sample_id_clean", suffix = c(".ps", ""))
}

## ---- load environmental data with biome and depth labels --------------------
load_env <- function() {
  env <- read_csv(file.path(ENV_DIR, "depth_analysis_data.csv"), show_col_types = FALSE)
  env <- env |> filter(!sample_id_clean %in% EXCLUDE_SAMPLES)
  env |>
    mutate(
      ## Mean annual daily temperature derived from the supplied CFS/CIPHA
      ## 2008–2018 mean daily maximum and minimum temperature fields.
      mat_10_year = (x10_year_t_max_mean + x10_year_t_min_mean) / 2,
      biome = case_when(
        node %in% BOREAL_NODES   ~ "boreal",
        node %in% PARKLAND_NODES ~ "parkland",
        TRUE                     ~ NA_character_
      ),
      depth_factor = factor(depth, levels = c("A", "B"),
                            labels = c("surface", "subsoil"))
    )
}

## ---- build VT-level abundance matrix for AMF (18S) --------------------------
##
## The 18S library was sequenced forward-only (~254 bp) because the NS31/AML2
## amplicon (~550 bp) exceeded the merged-read length. The first 22 bp of each
## ASV fall in a conserved region not present in MaarjAM reference sequences
## (5'-AGCAGCCGCGGTAATTCCAGCT); we trim that leader before VT assignment and
## fall back to it if the otu_table sequences contain it untrimmed.
build_vt_matrix <- function(ps_18s) {
  amf_otu <- as.data.frame(as(otu_table(ps_18s), "matrix"))
  if (taxa_are_rows(ps_18s)) amf_otu <- as.data.frame(t(amf_otu))

  asv_map <- read_csv(file.path(TAX_DIR, "easton_18s_asv_map.csv"),
                      show_col_types = FALSE)
  vtx_map <- read_csv(file.path(TAX_DIR, "maarjam_vsearch_mapping_final.csv"),
                      show_col_types = FALSE)
  seq_to_amf <- left_join(asv_map,
                          vtx_map[, c("asv_id", "is_amf", "vtx_vsearch")],
                          by = "asv_id")

  otu_seqs     <- colnames(amf_otu)
  matched_seqs <- intersect(otu_seqs, seq_to_amf$sequence)
  if (length(matched_seqs) == 0) {
    leader           <- "AGCAGCCGCGGTAATTCCAGCT"
    otu_seqs_trimmed <- ifelse(startsWith(otu_seqs, leader),
                               substr(otu_seqs, nchar(leader) + 1, nchar(otu_seqs)),
                               otu_seqs)
    vtx_lookup       <- setNames(seq_to_amf$vtx_vsearch, seq_to_amf$sequence)
    vtx_flags        <- vtx_lookup[otu_seqs_trimmed]
  } else {
    vtx_lookup <- setNames(seq_to_amf$vtx_vsearch, seq_to_amf$sequence)
    vtx_flags  <- vtx_lookup[otu_seqs]
  }

  has_vtx   <- !is.na(vtx_flags) & vtx_flags != ""
  vtx_cols  <- which(has_vtx)
  unique_vt <- sort(unique(vtx_flags[has_vtx]))

  vt_mat <- matrix(0, nrow = nrow(amf_otu), ncol = length(unique_vt),
                   dimnames = list(rownames(amf_otu), unique_vt))
  for (i in seq_along(vtx_cols)) {
    col_idx <- vtx_cols[i]
    vt_code <- vtx_flags[col_idx]
    vt_mat[, vt_code] <- vt_mat[, vt_code] + amf_otu[, col_idx]
  }
  cat("VT matrix:", nrow(vt_mat), "samples x", ncol(vt_mat), "VT\n")
  vt_mat
}

## ---- Aitchison distance from a samples-by-taxa matrix -----------------------
aitchison_distance_matrix <- function(mat, pseudocount = 0.5) {
  mat_pseudo <- mat + pseudocount
  gm  <- exp(rowMeans(log(mat_pseudo)))
  clr <- log(mat_pseudo / gm)
  stats::dist(clr)
}

## ---- Aitchison distance from a phyloseq object ------------------------------
aitchison_distance <- function(ps, pseudocount = 0.5) {
  otu <- as.matrix(as(otu_table(ps), "matrix"))
  if (taxa_are_rows(ps)) otu <- t(otu)
  otu_pseudo <- otu + pseudocount
  gm  <- exp(rowMeans(log(otu_pseudo)))
  clr <- log(otu_pseudo / gm)
  stats::dist(clr)
}

## ---- PCoA on Aitchison distances --------------------------------------------
aitchison_pcoa <- function(ps, pseudocount = 0.5) {
  d <- aitchison_distance(ps, pseudocount = pseudocount)
  list(ordination = ape::pcoa(d), distance = d)
}

## ---- thin PERMANOVA wrapper -------------------------------------------------
##
## by = "margin" requests Type-III sums of squares; nperm controls the number
## of permutations; strata is passed to permute::how() for restricted
## permutations within blocks (e.g. depth labels swapped within stands).
run_permanova <- function(dist_obj, metadata, formula_rhs, nperm = 999,
                          by = "margin", strata = NULL) {
  fml <- as.formula(paste("dist_obj ~", formula_rhs))
  set.seed(42)
  if (!is.null(strata)) {
    perm <- permute::how(blocks = strata, nperm = nperm)
    res  <- adonis2(fml, data = metadata, permutations = perm, by = by)
  } else {
    res  <- adonis2(fml, data = metadata, permutations = nperm, by = by)
  }
  res |> as.data.frame() |> tibble::rownames_to_column("term")
}

## ---- subset a distance object by sample names -------------------------------
subset_dist <- function(dist_obj, keep_names) {
  all_names <- labels(dist_obj)
  keep_idx  <- all_names %in% keep_names
  as.dist(as.matrix(dist_obj)[keep_idx, keep_idx])
}


## ============================================================================
## Design-based permutation inference
## ----------------------------------------------------------------------------
## Samples are nested: 2 depths within 27 stands within 9 nodes. Each predictor
## is tested by permuting the units at which it varies:
##   * node-level terms (latitude, biome, climate): values permuted among nodes
##     (optionally within strata, e.g. MAT within biome);
##   * stand-level covariates that vary mainly among nodes (pH, sand, basal
##     area): node means permuted among nodes, as for climate;
##   * within-stand terms (depth): stand conditioned out (dbRDA + Condition) and
##     residuals permuted within stands (vegan, see test_within_stand()).
## Depth x moderator interactions use paired surface->subsoil CLR shifts
## (paired_shift_test()), with the moderator permuted at its own level.
## ============================================================================

## align metadata rows to a dist object, failing loudly on any mismatch
align_md <- function(md, d, id_col = "sample_name") {
  md <- md[match(labels(d), md[[id_col]]), , drop = FALSE]
  stopifnot(identical(md[[id_col]], labels(d)))
  md
}

gower_centre <- function(d) {
  A <- -0.5 * as.matrix(d)^2
  n <- nrow(A)
  A - outer(rowMeans(A), rep(1, n)) - outer(rep(1, n), colMeans(A)) + mean(A)
}

.hat <- function(X) {
  q <- qr(X)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  Q %*% t(Q)
}

## marginal pseudo-F and R2 for `focal` given the other `terms`
## (identical to adonis2(by = "margin") for main-effect models)
marginal_F <- function(G, md, terms, focal) {
  Xf <- model.matrix(reformulate(terms), md)
  rest <- setdiff(terms, focal)   # focal may be several columns (tested jointly)
  Xr <- if (length(rest)) model.matrix(reformulate(rest), md) else matrix(1, nrow(md), 1)
  Hf <- .hat(Xf); Hr <- .hat(Xr)
  df_t <- qr(Xf)$rank - qr(Xr)$rank
  df_r <- nrow(md) - qr(Xf)$rank
  ss_t <- sum((Hf - Hr) * G); ss_r <- sum((diag(nrow(md)) - Hf) * G)
  c(F = (ss_t / df_t) / (ss_r / df_r), R2 = ss_t / sum(diag(G)), df = df_t)
}

## permutation p for `focal` (one or more columns, permuted jointly), at the
## level of `unit`:
## unit = "node":  focal replaced by node means (first value if non-numeric);
##                 nodes permuted among nodes, within `strata` if given
## unit = "stand": focal values kept per (stand, depth); stands permuted among
##                 stands with the same set of depths (paired vs unpaired)
unit_perm_test <- function(G, md, terms, focal, unit = c("node", "stand"),
                           strata = NULL, nperm = NPERM, seed = 42) {
  unit <- match.arg(unit)
  agg <- function(x) if (is.numeric(x)) mean(x) else x[1]
  if (unit == "node") {
    u <- md$node
    units <- sort(unique(u))
    tabs <- lapply(focal, function(f) sapply(split(md[[f]], u), agg)[units])
    names(tabs) <- focal
    for (f in focal) md[[f]] <- unname(tabs[[f]][u])
    grp <- if (is.null(strata)) setNames(rep("all", length(units)), units)
           else sapply(split(md[[strata]], u), function(x) x[1])[units]
    apply_map <- function(map) { mp <- md; for (f in focal) mp[[f]] <- unname(tabs[[f]][map[u]]); mp }
  } else {
    u <- md$stand_id
    units <- sort(unique(u))
    key <- paste(u, md$depth)
    tabs <- lapply(focal, function(f) setNames(md[[f]], key)); names(tabs) <- focal
    grp <- sapply(split(md$depth, u), function(x) paste(sort(unique(x)), collapse = ""))[units]
    apply_map <- function(map) { mp <- md; for (f in focal) mp[[f]] <- unname(tabs[[f]][paste(map[u], md$depth)]); mp }
  }
  draw_map <- function() {
    map <- setNames(units, units)
    for (g in unique(grp)) { i <- which(grp == g); map[i] <- units[i][sample.int(length(i))] }
    map
  }
  obs <- marginal_F(G, md, terms, focal)
  set.seed(seed)
  perm_F <- vapply(seq_len(nperm), function(i) marginal_F(G, apply_map(draw_map()), terms, focal)[["F"]],
                   numeric(1))
  c(obs, p = (sum(perm_F >= obs[["F"]] - 1e-10) + 1) / (nperm + 1), nperm = nperm)
}

## within-stand test: stand conditioned out, residuals permuted within stands
test_within_stand <- function(d, md, term, covars = character(0), nperm = NPERM, seed = 42) {
  md <- align_md(md, d)
  cond <- paste(c("stand_id", covars), collapse = " + ")
  f <- as.formula(paste("d ~", term, "+ Condition(", cond, ")"))
  set.seed(seed)
  a <- anova(vegan::dbrda(f, data = md),
             permutations = permute::how(blocks = factor(md$stand_id), nperm = nperm))
  c(F = a$F[1], p = a$`Pr(>F)`[1], nperm = nperm)
}

## paired surface->subsoil shifts: rows = paired stands, columns = CLR taxa
paired_shifts <- function(clr, md) {
  md <- md[md$sample_name %in% rownames(clr), ]
  pairs <- md |> group_by(stand_id) |> filter(n_distinct(depth) == 2) |> ungroup()
  st <- pairs |> group_by(stand_id) |>
    summarise(node = first(node), biome = first(biome),
              across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
    arrange(stand_id)
  D <- t(vapply(st$stand_id, function(s) {
    clr[pairs$sample_name[pairs$stand_id == s & pairs$depth == "B"], ] -
      clr[pairs$sample_name[pairs$stand_id == s & pairs$depth == "A"], ]
  }, numeric(ncol(clr))))
  rownames(D) <- st$stand_id
  list(D = D, st = st)
}

## does the paired depth shift vary with a moderator? (depth x moderator)
## moderator values are stand means across the two depths
paired_shift_test <- function(ps_obj, moderator, unit = c("node", "stand"),
                              covars = character(0), strata = NULL, nperm = NPERM) {
  unit <- match.arg(unit)
  D <- ps_obj$D; st <- ps_obj$st
  st$depth <- "pair"
  G <- tcrossprod(scale(D, scale = FALSE))
  unit_perm_test(G, st, c(moderator, covars), moderator, unit = unit,
                 strata = strata, nperm = nperm)
}

clr_matrix <- function(mat, pseudocount = 0.5) {
  L <- log(mat + pseudocount)
  L - rowMeans(L)
}


## ---- indicator-species permutation settings ---------------------------------
## Depth indicators: depth labels permuted within stands (paired design), with
## enough permutations that BH-adjusted p-values are not bounded by resolution.
NPERM_IND <- as.integer(Sys.getenv("NPERM_IND", "99999"))

## Biome indicators: biome is a node-level label, so p-values come from the
## exact null over all C(9, 4) = 126 assignments of "parkland" to 4 of 9 nodes.
indval_biome_exact <- function(rel, node, biome) {
  nodes <- sort(unique(node))
  park_obs <- sort(unique(node[biome == "parkland"]))
  stat_for <- function(park) {
    cl <- factor(ifelse(node %in% park, "parkland", "boreal"))
    apply(indicspecies::strassoc(rel, cl, func = "IndVal.g"), 1, max)
  }
  obs <- stat_for(park_obs)
  combs <- combn(nodes, length(park_obs), simplify = FALSE)
  null <- vapply(combs, stat_for, numeric(ncol(rel)))
  data.frame(taxon = colnames(rel), stat_node = obs,
             p_node_exact = rowMeans(null >= obs - 1e-12),
             n_assignments = length(combs))
}

cat("helpers.R loaded; package root =", PACKAGE_ROOT, "\n")

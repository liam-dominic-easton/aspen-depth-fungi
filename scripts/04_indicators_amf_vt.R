## 04_indicators_amf_vt.R — AMF (18S, VT-level) indicators.
##
## Outputs (in outputs/tables/):
##   indicator_vtx_by_depth.csv         — AMF VT depth indicators (none after BH; quoted in the text)
##   indicator_vtx_by_biome.csv         — AMF VT biome indicators (none after BH; not reported)
##   (plus *_full.csv versions, vtx_taxonomy_lookup.csv and vtx_abundance_by_depth.csv)
##
## AMF VT indicators use min_prev = 3 samples (low because there are relatively
## few unique VT in this dataset).
## IndVal.g. Depth: within-stand permutations (NPERM_IND). Biome: exact node-level
## null (126 assignments).

source("scripts/helpers.R")
library(indicspecies)
library(glue)

# ---------- Helper: generic indicator analysis ----------
run_indicator <- function(otu_mat, grouping, label, min_prev = 5,
                          control = how(nperm = NPERM_IND)) {
  valid <- !is.na(grouping) & grouping != ""
  otu_sub <- otu_mat[valid, , drop = FALSE]
  grp <- factor(grouping[valid])

  prev <- colSums(otu_sub > 0)
  otu_sub <- otu_sub[, prev >= min_prev, drop = FALSE]

  if (ncol(otu_sub) < 2) { cat(glue("[{label}] Too few taxa after prevalence filter"), "\n"); return(NULL) }

  rel <- sweep(otu_sub, 1, rowSums(otu_sub), "/")
  rel[is.nan(as.matrix(rel))] <- 0

  cat(glue("[{label}] {nrow(rel)} samples, {ncol(rel)} taxa, {nlevels(grp)} groups"), "\n")

  set.seed(42)
  ind <- tryCatch(
    multipatt(rel, grp, func = "IndVal.g", duleg = TRUE, control = control),
    error = function(e) { cat(glue("[{label}] Error: {e$message}"), "\n"); NULL }
  )
  if (is.null(ind)) return(NULL)

  sign_tbl <- as.data.frame(ind$sign) %>% tibble::rownames_to_column("taxon") %>%
    mutate(p_bh = p.adjust(p.value, method = "BH"))

  group_names <- levels(grp)
  mp_cols <- paste0("s.", group_names)
  lookup <- setNames(group_names, mp_cols)
  gcols <- intersect(mp_cols, colnames(sign_tbl))

  if (length(gcols) > 0) {
    sign_tbl$indicator_for <- apply(sign_tbl[, gcols, drop = FALSE], 1, function(x) {
      paste(unname(lookup[gcols[x == 1]]), collapse = "+")
    })
  } else {
    sign_tbl$indicator_for <- NA_character_
  }

  sig <- sign_tbl %>% filter(p.value < 0.05) %>% arrange(p.value)
  cat(glue("[{label}] Significant: {nrow(sig)}/{nrow(sign_tbl)}"), "\n")

  list(full = sign_tbl, significant = sig, model = ind, rel = rel, valid = valid)
}

# ---------- Load data ----------
cat("Loading data...\n")
ps_18s <- load_ps("18S")
env    <- load_env()

# =============================================================================
# 1. AMF indicator analysis at VTX level
# =============================================================================
cat("\n=== AMF (18S) VTX-level indicator analysis ===\n")

# Bridge: phyloseq sequences → ASV_ID → VTX
asv_map <- read.csv(
  file.path(TAX_DIR, "easton_18s_asv_map.csv"),
  stringsAsFactors = FALSE
)
vtx_map <- read.csv(
  file.path(TAX_DIR, "maarjam_vsearch_mapping_final.csv"),
  stringsAsFactors = FALSE
)

cat(glue("ASV map: {nrow(asv_map)} entries"), "\n")
cat(glue("VTX map: {nrow(vtx_map)} entries, {sum(vtx_map$is_amf == 'True' | vtx_map$is_amf == TRUE)} AMF"), "\n")

# Build sequence → VTX lookup
seq_to_asv <- setNames(asv_map$asv_id, asv_map$sequence)
asv_to_vtx <- setNames(vtx_map$vtx_vsearch, vtx_map$asv_id)

# Parse MaarjAM FASTA headers for VTX → Family/Genus lookup
cat("\n--- Parsing MaarjAM FASTA for VTX taxonomy ---\n")
fasta_path <- file.path(DATA_DIR, "reference", "maarjam_SSU_type_2021", "maarjam_database_SSU.fasta")
fasta_lines <- readLines(fasta_path, warn = FALSE)
header_lines <- fasta_lines[grepl("^>", fasta_lines)]
cat(glue("MaarjAM headers: {length(header_lines)}"), "\n")

# Parse: >gi|ID|gb|ACC| Family Genus ... VTXnnnnn
vtx_taxonomy <- data.frame(vtx = character(), family = character(), genus = character(),
                            stringsAsFactors = FALSE)

for (h in header_lines) {
  # Extract VTX code
  vtx_match <- regmatches(h, regexpr("VTX\\d+", h))
  if (length(vtx_match) == 0) next

  # Extract text after last |
  desc <- sub("^.*\\|\\s*", "", h)
  words <- strsplit(trimws(desc), "\\s+")[[1]]

  family <- NA_character_
  genus <- NA_character_

  if (length(words) >= 1 && grepl("aceae$|ineae$", words[1])) {
    family <- words[1]
  }
  if (length(words) >= 2 && !grepl("^(sp\\.|18S|rDNA|TYPE|VTX)", words[2])) {
    genus <- words[2]
  }

  vtx_taxonomy <- rbind(vtx_taxonomy, data.frame(
    vtx = vtx_match, family = family, genus = genus, stringsAsFactors = FALSE
  ))
}

# Deduplicate: take most common family/genus per VTX
vtx_tax_lookup <- vtx_taxonomy %>%
  filter(!is.na(family) | !is.na(genus)) %>%
  group_by(vtx) %>%
  summarise(
    family = names(sort(table(family[!is.na(family)]), decreasing = TRUE))[1],
    genus = names(sort(table(genus[!is.na(genus)]), decreasing = TRUE))[1],
    n_seqs = n(),
    .groups = "drop"
  )

cat(glue("VTX with taxonomy: {nrow(vtx_tax_lookup)} VTX"), "\n")
cat(glue("Families: {length(unique(vtx_tax_lookup$family[!is.na(vtx_tax_lookup$family)]))}"), "\n")
cat(glue("Genera: {length(unique(vtx_tax_lookup$genus[!is.na(vtx_tax_lookup$genus)]))}"), "\n")

# Save the lookup
write_csv(vtx_tax_lookup, file.path(TABLES_DIR, "vtx_taxonomy_lookup.csv"))
cat("Saved: vtx_taxonomy_lookup.csv\n")

# Build VTX-level abundance matrix
cat("\n--- Building VTX abundance matrix ---\n")
amf_18s_otu <- as.data.frame(as(otu_table(ps_18s), "matrix"))
if (taxa_are_rows(ps_18s)) amf_18s_otu <- as.data.frame(t(amf_18s_otu))

# Map each column (sequence) to VTX
seq_names <- colnames(amf_18s_otu)
vtx_assignments <- sapply(seq_names, function(s) {
  aid <- seq_to_asv[s]
  if (is.na(aid)) return(NA_character_)
  vt <- asv_to_vtx[aid]
  if (is.na(vt) || vt == "") return(NA_character_)
  vt
})

amf_cols <- !is.na(vtx_assignments)
cat(glue("18S ASVs with VTX: {sum(amf_cols)} of {length(seq_names)}"), "\n")

if (sum(amf_cols) > 0) {
  amf_otu_sub <- amf_18s_otu[, amf_cols, drop = FALSE]
  amf_vtx <- vtx_assignments[amf_cols]

  # Aggregate by VTX
  unique_vtx <- unique(amf_vtx)
  vtx_mat <- sapply(unique_vtx, function(v) {
    cols <- which(amf_vtx == v)
    if (length(cols) == 1) return(amf_otu_sub[, cols])
    rowSums(amf_otu_sub[, cols, drop = FALSE])
  })
  colnames(vtx_mat) <- unique_vtx
  rownames(vtx_mat) <- rownames(amf_18s_otu)

  cat(glue("VTX matrix: {nrow(vtx_mat)} samples x {ncol(vtx_mat)} VTX"), "\n")

  # Match to env
  matched_18s <- match_ps_to_env(ps_18s, env)
  vtx_matched <- vtx_mat[rownames(vtx_mat) %in% matched_18s$sample_name, , drop = FALSE]
  vtx_matched <- vtx_matched[match(
    matched_18s$sample_name[matched_18s$sample_name %in% rownames(vtx_matched)],
    rownames(vtx_matched)
  ), , drop = FALSE]
  matched_18s_sub <- matched_18s %>% filter(sample_name %in% rownames(vtx_matched))

  cat(glue("Matched VTX matrix: {nrow(vtx_matched)} samples x {ncol(vtx_matched)} VTX"), "\n")

  # Run indicator analysis by depth
  cat("\n--- AMF VTX indicators by DEPTH ---\n")
  ind_vtx_depth <- run_indicator(vtx_matched, as.character(matched_18s_sub$depth_factor), "AMF-VTX-depth", min_prev = 3,
                                 control = how(blocks = factor(paste0(matched_18s_sub$node, matched_18s_sub$stand)),
                                               nperm = NPERM_IND))

  if (!is.null(ind_vtx_depth)) {
    # Annotate with taxonomy
    ind_vtx_depth$significant <- ind_vtx_depth$significant %>%
      left_join(vtx_tax_lookup, by = c("taxon" = "vtx"))

    ind_vtx_depth$full <- ind_vtx_depth$full %>%
      left_join(vtx_tax_lookup, by = c("taxon" = "vtx"))

    cat("\nSignificant VTX indicators:\n")
    if (nrow(ind_vtx_depth$significant) > 0) {
      print(ind_vtx_depth$significant %>% select(taxon, indicator_for, stat, p.value, family, genus))
    } else {
      cat("  (none)\n")
    }

    write_csv(ind_vtx_depth$significant, file.path(TABLES_DIR, "indicator_vtx_by_depth.csv"))
    write_csv(ind_vtx_depth$full, file.path(TABLES_DIR, "indicator_vtx_by_depth_full.csv"))
    cat("Saved: indicator_vtx_by_depth.csv, indicator_vtx_by_depth_full.csv\n")
  }

  # VTX indicators by biome
  cat("\n--- AMF VTX indicators by BIOME ---\n")
  ## biome is node-level: p-values come from the exact node-level null over the
  ## 126 assignments of "parkland" to 4 of 9 nodes (permuting the 52 samples
  ## freely would pseudoreplicate)
  ind_vtx_biome <- run_indicator(vtx_matched, matched_18s_sub$biome, "AMF-VTX-biome", min_prev = 3,
                                 control = how(nperm = 99))

  if (!is.null(ind_vtx_biome)) {
    ex <- indval_biome_exact(ind_vtx_biome$rel, matched_18s_sub$node[ind_vtx_biome$valid],
                             matched_18s_sub$biome[ind_vtx_biome$valid])
    ind_vtx_biome$full <- ind_vtx_biome$full |> left_join(ex, by = "taxon") |>
      mutate(p_node_bh = p.adjust(p_node_exact, "BH"))
    ind_vtx_biome$significant <- ind_vtx_biome$full |> filter(p_node_exact < 0.05) |> arrange(p_node_exact)
    ind_vtx_biome$significant <- ind_vtx_biome$significant %>%
      left_join(vtx_tax_lookup, by = c("taxon" = "vtx"))

    cat("\nSignificant VTX biome indicators:\n")
    if (nrow(ind_vtx_biome$significant) > 0) {
      print(ind_vtx_biome$significant %>% select(taxon, indicator_for, stat, p_node_exact, family, genus))
    } else {
      cat("  (none)\n")
    }

    write_csv(ind_vtx_biome$significant, file.path(TABLES_DIR, "indicator_vtx_by_biome.csv"))
    write_csv(ind_vtx_biome$full, file.path(TABLES_DIR, "indicator_vtx_by_biome_full.csv"))
    cat("Saved: indicator_vtx_by_biome.csv\n")
  }

  # VTX abundance by depth summary
  cat("\n--- VTX abundance by depth ---\n")
  surf_idx <- matched_18s_sub$depth_factor == "surface"
  sub_idx  <- matched_18s_sub$depth_factor == "subsoil"

  vtx_depth <- data.frame(
    vtx = colnames(vtx_matched),
    mean_surface = colMeans(vtx_matched[surf_idx, , drop = FALSE]),
    mean_subsoil = colMeans(vtx_matched[sub_idx, , drop = FALSE]),
    prev_surface = colSums(vtx_matched[surf_idx, , drop = FALSE] > 0),
    prev_subsoil = colSums(vtx_matched[sub_idx, , drop = FALSE] > 0),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      total = mean_surface + mean_subsoil,
      depth_ratio = ifelse(mean_surface > 0, mean_subsoil / mean_surface, NA)
    ) %>%
    left_join(vtx_tax_lookup, by = c("vtx" = "vtx")) %>%
    arrange(desc(total))

  write_csv(vtx_depth, file.path(TABLES_DIR, "vtx_abundance_by_depth.csv"))
  cat("Saved: vtx_abundance_by_depth.csv\n")
  cat("\nTop 15 VTX by abundance:\n")
  print(vtx_depth %>% select(vtx, family, genus, mean_surface, mean_subsoil, depth_ratio) %>% head(15))
}

cat("\n=== AMF indicator script complete ===\n")

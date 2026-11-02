## 05_nestedness_betadisper.R — beta-dispersion + AMF nestedness.
##
## Outputs (in outputs/tables/ and outputs/figures/):
##   betadisper_by_depth.csv           — beta-dispersion test per library (P values quoted in the text)
##   betadisper_by_depth_points.csv    — distance of each sample to its depth-group median
##   nestedness_results.csv            — NODF (overall, rows, columns) + paired overlap and its null (Table S6)
##   vt_depth_overlap.csv              — per-stand paired VT overlap and nulls (summarized in Table S6)
##   vt_richness_by_depth.csv          — VT richness per sample
##   betadisper_by_depth.png           — dispersion plot (not in the manuscript)
##
## Nestedness is tested only for AMF VT because VT richness shows the
## strongest depth decline; ITS guild richness differences are weaker.
## NODF is computed with vegan::nestednodf using a binary VT × sample matrix
## ordered by depth, with null model "r1" and NPERM simulations.

source("scripts/helpers.R")
library(readr)
library(glue)
library(ggplot2)

cat("=== Nestedness & VT richness by depth ===\n")

# =============================================================================
# 1. Load 18S phyloseq and build VT-level abundance matrix
# =============================================================================
cat("Loading 18S phyloseq...\n")
ps_18s <- load_ps("18S")
env    <- load_env()

# Get OTU table (samples × sequences)
amf_otu <- as.data.frame(as(otu_table(ps_18s), "matrix"))
if (taxa_are_rows(ps_18s)) amf_otu <- as.data.frame(t(amf_otu))
cat("OTU table: ", nrow(amf_otu), " samples × ", ncol(amf_otu), " ASVs\n")

# Bridge sequences → ASV_IDs → VTX codes
PREV_DATA <- TAX_DIR

asv_map <- read_csv(file.path(PREV_DATA, "easton_18s_asv_map.csv"),
                    show_col_types = FALSE)
vtx_map <- read_csv(file.path(PREV_DATA, "maarjam_vsearch_mapping_final.csv"),
                    show_col_types = FALSE)

seq_to_amf <- asv_map %>%
  left_join(vtx_map %>% select(asv_id, is_amf, vtx_vsearch), by = "asv_id")

otu_seqs <- colnames(amf_otu)

# Try direct match first, then trimmed match
matched_seqs <- intersect(otu_seqs, seq_to_amf$sequence)

if (length(matched_seqs) == 0) {
  leader <- "AGCAGCCGCGGTAATTCCAGCT"
  otu_seqs_trimmed <- ifelse(
    startsWith(otu_seqs, leader),
    substr(otu_seqs, nchar(leader) + 1, nchar(otu_seqs)),
    otu_seqs
  )
  vtx_lookup <- setNames(seq_to_amf$vtx_vsearch, seq_to_amf$sequence)
  vtx_flags <- vtx_lookup[otu_seqs_trimmed]
} else {
  vtx_lookup <- setNames(seq_to_amf$vtx_vsearch, seq_to_amf$sequence)
  vtx_flags <- vtx_lookup[otu_seqs]
}

has_vtx <- !is.na(vtx_flags) & vtx_flags != ""
cat("ASVs with VT assignment: ", sum(has_vtx), "\n")
cat("Unique VT: ", length(unique(vtx_flags[has_vtx])), "\n")

# Aggregate reads to VT level
vtx_cols <- which(has_vtx)
unique_vt <- sort(unique(vtx_flags[has_vtx]))

vt_matrix <- matrix(0, nrow = nrow(amf_otu), ncol = length(unique_vt),
                    dimnames = list(rownames(amf_otu), unique_vt))

for (i in seq_along(vtx_cols)) {
  col_idx <- vtx_cols[i]
  vt_code <- vtx_flags[col_idx]
  vt_matrix[, vt_code] <- vt_matrix[, vt_code] + amf_otu[, col_idx]
}

cat("VT abundance matrix: ", nrow(vt_matrix), " samples × ", ncol(vt_matrix), " VT\n")

# =============================================================================
# 2. Match samples to env, split by depth
# =============================================================================
md <- match_ps_to_env(ps_18s, env)

# Keep only samples in VT matrix
md <- md %>% filter(sample_name %in% rownames(vt_matrix))

# Create paired stand identifier (node + stand number)
md <- md %>% mutate(stand_id = paste0(node, stand))

# Split by depth
surface_ids <- md %>% filter(depth_factor == "surface") %>% pull(sample_name)
subsoil_ids <- md %>% filter(depth_factor == "subsoil") %>% pull(sample_name)

cat("Surface samples: ", length(surface_ids), "\n")
cat("Subsoil samples: ", length(subsoil_ids), "\n")

# =============================================================================
# 3. VT richness by depth — descriptive paired summary
# =============================================================================
cat("\n=== VT richness by depth ===\n")

vt_richness_df <- md %>%
  mutate(vt_richness = rowSums(vt_matrix[sample_name, ] > 0)) %>%
  select(sample_name, sample_id_clean, node, stand, depth_factor, stand_id, vt_richness)

# Summary by depth
richness_summary <- vt_richness_df %>%
  group_by(depth_factor) %>%
  summarise(
    n = n(),
    mean_vt = round(mean(vt_richness), 1),
    sd_vt = round(sd(vt_richness), 1),
    median_vt = median(vt_richness),
    min_vt = min(vt_richness),
    max_vt = max(vt_richness),
    .groups = "drop"
  )
cat("\nVT richness summary:\n")
print(as.data.frame(richness_summary))

# Paired test: match surface-subsoil within each stand
paired <- vt_richness_df %>%
  tidyr::pivot_wider(
    id_cols = stand_id,
    names_from = depth_factor,
    values_from = vt_richness
  ) %>%
  filter(!is.na(surface) & !is.na(subsoil))

cat("\nPaired stands: ", nrow(paired), "\n")
cat("Surface > subsoil: ", sum(paired$surface > paired$subsoil), "/", nrow(paired), "\n")
cat("Surface = subsoil: ", sum(paired$surface == paired$subsoil), "/", nrow(paired), "\n")
cat("Surface < subsoil: ", sum(paired$surface < paired$subsoil), "/", nrow(paired), "\n")

## Inference on AMF richness is the Bayesian paired model (06); this section
## is descriptive.
cat("Mean VT richness, paired stands: surface", round(mean(paired$surface), 2),
    "| subsoil", round(mean(paired$subsoil), 2), "\n")
cat("Mean VT richness, all subsoil samples:",
    round(mean(vt_richness_df$vt_richness[vt_richness_df$depth_factor == "subsoil"]), 2), "\n")

# Save
write_csv(vt_richness_df, file.path(TABLES_DIR, "vt_richness_by_depth.csv"))
cat("Saved: vt_richness_by_depth.csv\n")

# =============================================================================
# 4. Nestedness test (NODF) — are subsoil VT nested within surface VT?
# =============================================================================
cat("\n=== Nestedness analysis ===\n")
nodf_empirical_p <- NA_real_

# Create presence/absence VT matrix
vt_pa <- (vt_matrix > 0) * 1

# Site-level nestedness: all samples
cat("\n--- All-sample nestedness (NODF) ---\n")
nodf_all <- tryCatch({
  vegan::nestednodf(vt_pa)
}, error = function(e) {
  cat("nestednodf error: ", e$message, "\n")
  NULL
})

if (!is.null(nodf_all)) {
  cat("NODF (columns/species): ", round(nodf_all$statistic["N.columns"], 2), "\n")
  cat("NODF (rows/sites): ", round(nodf_all$statistic["N.rows"], 2), "\n")
  cat("NODF (overall): ", round(nodf_all$statistic["NODF"], 2), "\n")

  # Null model test using oecosimu
  cat("\nRunning null model test (", NPERM, " simulations)...\n")
  set.seed(42)
  nodf_test <- tryCatch({
    vegan::oecosimu(vt_pa, vegan::nestednodf, method = "r1", nsimul = NPERM,
                    statistic = "statistic")
  }, error = function(e) {
    cat("oecosimu error: ", e$message, "\n")
    NULL
  })

  if (!is.null(nodf_test)) {
    observed_nodf <- unname(nodf_test$oecosimu$statistic["NODF"])
    simulated_nodf <- nodf_test$oecosimu$simulated[3, ]
    nodf_empirical_p <- (1 + sum(simulated_nodf >= observed_nodf)) /
      (1 + length(simulated_nodf))
    cat("Null model test (method = r1):\n")
    cat("  Observed NODF: ", round(nodf_test$oecosimu$statistic["NODF"], 2), "\n")
    cat("  Mean null NODF: ", round(mean(nodf_test$oecosimu$simulated[3, ]), 2), "\n")
    cat("  empirical upper-tail p: ", round(nodf_empirical_p, 4), "\n")
    cat("  z: ", round(nodf_test$oecosimu$z["NODF"], 3), "\n")
  }
}

# Depth-stratified nestedness: test whether subsoil is nested within surface
# For each stand, check: how many subsoil VT are also found at surface?
cat("\n--- Paired depth nestedness check ---\n")

paired_nest <- paired %>%
  rowwise() %>%
  mutate(stand_id = stand_id) %>%
  ungroup()

# For each paired stand, compute VT overlap
nest_stats <- list()
for (i in seq_len(nrow(paired))) {
  sid <- paired$stand_id[i]
  surf_name <- md %>% filter(stand_id == sid, depth_factor == "surface") %>% pull(sample_name)
  sub_name  <- md %>% filter(stand_id == sid, depth_factor == "subsoil") %>% pull(sample_name)

  if (length(surf_name) == 1 && length(sub_name) == 1) {
    surf_vt <- unique_vt[vt_pa[surf_name, ] > 0]
    sub_vt  <- unique_vt[vt_pa[sub_name, ] > 0]

    shared <- length(intersect(surf_vt, sub_vt))
    only_surf <- length(setdiff(surf_vt, sub_vt))
    only_sub  <- length(setdiff(sub_vt, surf_vt))

    # Nestedness: proportion of subsoil VT that are also in surface
    prop_nested <- if (length(sub_vt) > 0) shared / length(sub_vt) else NA

    nest_stats[[i]] <- data.frame(
      stand_id = sid,
      n_surface_vt = length(surf_vt),
      n_subsoil_vt = length(sub_vt),
      n_shared = shared,
      n_only_surface = only_surf,
      n_only_subsoil = only_sub,
      prop_subsoil_nested = round(prop_nested, 3),
      stringsAsFactors = FALSE
    )
  }
}

nest_df <- bind_rows(nest_stats)
cat("\nPaired depth VT overlap:\n")
cat("  Mean surface VT: ", round(mean(nest_df$n_surface_vt), 1), "\n")
cat("  Mean subsoil VT: ", round(mean(nest_df$n_subsoil_vt), 1), "\n")
cat("  Mean shared VT:  ", round(mean(nest_df$n_shared), 1), "\n")
cat("  Mean subsoil-only VT: ", round(mean(nest_df$n_only_subsoil), 1), "\n")
cat("  Mean proportion subsoil nested in surface: ",
    round(mean(nest_df$prop_subsoil_nested, na.rm = TRUE), 3), "\n")

# =============================================================================
# 5. Compile and save results
# =============================================================================
## Null for the paired overlap: the share of a subsoil
## sample's VT also found in the surface sample of OTHER stands in the same node,
## and of any other stand. The observed same-stand overlap is compared with these.
surf_md <- md %>% filter(depth_factor == "surface")
ov <- function(sub_name, surf_name) {
  sb <- vt_pa[sub_name, ] > 0
  if (sum(sb) == 0) return(NA_real_)
  sum(sb & vt_pa[surf_name, ] > 0) / sum(sb)
}
null_rows <- lapply(nest_df$stand_id, function(sid) {
  sub_name <- md %>% filter(stand_id == sid, depth_factor == "subsoil") %>% pull(sample_name)
  nd <- md$node[md$stand_id == sid][1]
  same_node <- surf_md %>% filter(node == nd, stand_id != sid) %>% pull(sample_name)
  others <- surf_md %>% filter(stand_id != sid) %>% pull(sample_name)
  data.frame(stand_id = sid,
             overlap_same_node_other_stands = if (length(same_node)) mean(sapply(same_node, ov, sub_name = sub_name)) else NA,
             overlap_any_other_stand = mean(sapply(others, ov, sub_name = sub_name)))
})
nest_df <- nest_df %>% left_join(bind_rows(null_rows), by = "stand_id")
cat("  Overlap null, other stands in same node: ",
    round(mean(nest_df$overlap_same_node_other_stands, na.rm = TRUE), 3), "\n")
cat("  Overlap null, any other stand: ", round(mean(nest_df$overlap_any_other_stand), 3), "\n")

nodf_stat <- if (!is.null(nodf_test)) nodf_test$oecosimu else NULL
nodf_null_mean <- if (!is.null(nodf_stat)) rowMeans(nodf_stat$simulated) else rep(NA, 3)
nodf_p <- if (!is.null(nodf_stat)) {
  obs <- nodf_stat$statistic
  (1 + rowSums(nodf_stat$simulated >= obs)) / (1 + ncol(nodf_stat$simulated))
} else rep(NA, 3)
nestedness_results <- data.frame(
  metric = c("nodf_overall", "nodf_rows_samples", "nodf_columns_vt",
             "nodf_overall_null_mean", "nodf_rows_null_mean", "nodf_columns_null_mean",
             "nodf_overall_z", "nodf_rows_z", "nodf_columns_z",
             "nodf_overall_p", "nodf_rows_p", "nodf_columns_p", "nodf_nsimul",
             "mean_surface_vt_paired", "mean_subsoil_vt_paired", "mean_subsoil_vt_all27",
             "mean_shared_vt", "mean_prop_subsoil_nested",
             "mean_overlap_null_same_node", "mean_overlap_null_any_stand",
             "n_paired_stands",
             ## matched comparison: MOO3 has no other surface
             ## sample in its node, so the same-node null uses 24 stands
             "mean_prop_subsoil_nested_matched", "n_stands_same_node_null"),
  value = c(
    nodf_all$statistic[c("NODF", "N.rows", "N.columns")],
    nodf_null_mean[c(3, 2, 1)],
    if (!is.null(nodf_stat)) nodf_stat$z[c("NODF", "N.rows", "N.columns")] else rep(NA, 3),
    nodf_p[c(3, 2, 1)], NPERM,
    mean(nest_df$n_surface_vt), mean(nest_df$n_subsoil_vt),
    mean(vt_richness_df$vt_richness[vt_richness_df$depth_factor == "subsoil"]),
    mean(nest_df$n_shared), mean(nest_df$prop_subsoil_nested, na.rm = TRUE),
    mean(nest_df$overlap_same_node_other_stands, na.rm = TRUE),
    mean(nest_df$overlap_any_other_stand), nrow(paired),
    mean(nest_df$prop_subsoil_nested[!is.na(nest_df$overlap_same_node_other_stands)], na.rm = TRUE),
    sum(!is.na(nest_df$overlap_same_node_other_stands))),
  stringsAsFactors = FALSE
)

write_csv(nestedness_results, file.path(TABLES_DIR, "nestedness_results.csv"))
cat("\nSaved: nestedness_results.csv\n")

# Also save the per-stand overlap table
write_csv(nest_df, file.path(TABLES_DIR, "vt_depth_overlap.csv"))
cat("Saved: vt_depth_overlap.csv\n")

# =============================================================================
# 6. Beta dispersion by depth for the two active community matrices
# =============================================================================
dispersion_test <- function(library_label, distance_object, metadata) {
  ids <- labels(distance_object)
  metadata <- metadata[match(ids, metadata$sample_name), , drop = FALSE]
  stopifnot(identical(ids, metadata$sample_name))
  group <- droplevels(metadata$depth_factor)
  fit <- vegan::betadisper(distance_object, group, bias.adjust = TRUE)
  ## depth labels permuted within stands, NPERM permutations
  set.seed(42)
  ctrl <- permute::how(blocks = factor(paste0(metadata$node, metadata$stand)), nperm = NPERM)
  test <- vegan::permutest(fit, permutations = ctrl)
  p_value <- unname(test$tab[1, "Pr(>F)"])
  list(
    points = tibble::tibble(
      library = library_label,
      sample_name = ids,
      depth = as.character(group),
      distance_to_centroid = fit$distances
    ),
    summary = tibble::tibble(
      library = library_label,
      F = unname(test$tab[1, "F"]),
      p_value = p_value,
      n = length(group)
    )
  )
}

ps_its <- load_ps("ITS")
md_its <- match_ps_to_env(ps_its, env)
its_disp <- dispersion_test("ITS", aitchison_distance(ps_its), md_its)
amf_disp <- dispersion_test("AMF VT", aitchison_distance_matrix(vt_matrix), md)

dispersion_points <- bind_rows(its_disp$points, amf_disp$points)
dispersion_summary <- bind_rows(its_disp$summary, amf_disp$summary)
write_csv(dispersion_points, file.path(TABLES_DIR, "betadisper_by_depth_points.csv"))
write_csv(dispersion_summary, file.path(TABLES_DIR, "betadisper_by_depth.csv"))

labels <- dispersion_summary |>
  mutate(label = paste0("permutation p = ", format.pval(p_value, digits = 2)),
         depth = "surface")

## depth colours as in Figures 2 and S1 (surface brown, subsoil teal), surface first
dispersion_plot <- ggplot(dispersion_points |>
                            mutate(depth = factor(depth, levels = c("surface", "subsoil"))),
                          aes(x = depth, y = distance_to_centroid, fill = depth)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.70) +
  geom_text(data = labels, aes(x = depth, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.4, size = 3.2) +
  facet_wrap(~library, scales = "free_y") +
  scale_fill_manual(values = c(surface = "#8C510A", subsoil = "#01665E")) +
  labs(x = NULL, y = "Distance to depth-group spatial median") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave(file.path(FIGS_DIR, "betadisper_by_depth.png"), dispersion_plot,
       width = 7.2, height = 3.5, dpi = 300)
cat("Saved: betadisper_by_depth.csv and betadisper_by_depth.png\n")

cat("\n=== Script 05 complete ===\n")

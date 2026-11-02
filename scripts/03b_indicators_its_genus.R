## 03b_indicators_its_genus.R — genus-level ITS indicator analysis (Table 2).
##
## Outputs (in outputs/tables/):
##   indicator_genera_by_depth.csv       — genus-aggregated depth indicators with
##                                         unadjusted p < 0.05; rows with
##                                         BH-adjusted p < 0.05 form Table 2
##   indicator_genera_by_depth_full.csv  — all tested genus-level groups
##   ecm_genus_env_by_depth.csv          — depth-stratified ECM genus-environment correlations
##   genus_depth_preferences.csv         — ECM and saprotroph genus depth ratios
##
## Aggregates ITS ASVs to genus level (using the phyloseq taxonomy), then runs
## indicspecies::multipatt (IndVal.g) on relative abundance with a prevalence
## filter (>= 5 samples); depth labels are permuted within stands (NPERM_IND
## permutations) and p-values are BH-adjusted across all tested groups.
## All ITS ASVs enter the general indicator analysis; FUNGuild subsets are
## used only for the separate ECM/saprotroph summaries.

source("scripts/helpers.R")
library(indicspecies)
library(ggplot2)
library(glue)

cat("Loading data...\n")
ps_its <- load_ps("ITS")
env    <- load_env()

matched <- match_ps_to_env(ps_its, env)
cat("Matched ITS samples:", nrow(matched), "\n")

its_otu <- as.data.frame(as(otu_table(ps_its), "matrix"))
if (taxa_are_rows(ps_its)) its_otu <- as.data.frame(t(its_otu))

its_tax <- as.data.frame(as(tax_table(ps_its), "matrix"))

# Load FUNGuild
guild_path <- file.path(TAX_DIR, "otu_table_for_FUNGuild.guilds.txt")
if (!file.exists(guild_path)) stop("FUNGuild file not found")
funguild <- read.delim(guild_path, stringsAsFactors = FALSE)
cat("FUNGuild loaded:", nrow(funguild), "rows\n")

# Subset OTU to matched samples
matched_names <- matched$sample_name
otu_matched <- its_otu[rownames(its_otu) %in% matched_names, ]
otu_matched <- otu_matched[match(
  matched$sample_name[matched$sample_name %in% rownames(otu_matched)],
  rownames(otu_matched)
), ]
matched_sub <- matched %>% filter(sample_name %in% rownames(otu_matched))

# ---------- Identify guild column and ECM/saprotroph ASVs ----------
guild_col <- intersect(c("Guild", "guild"), names(funguild))
gcol <- guild_col[1]
genus_col <- intersect(c("Genus", "genus"), names(its_tax))

ecm_rows <- funguild %>% filter(grepl("Ectomycorrhiz", !!sym(gcol), ignore.case = TRUE))
sapro_rows <- funguild %>% filter(grepl("Saprotroph", !!sym(gcol), ignore.case = TRUE))

ecm_seqs  <- intersect(ecm_rows[["OTU_ID"]], colnames(otu_matched))
sapro_seqs <- intersect(sapro_rows[["OTU_ID"]], colnames(otu_matched))
cat("ECM ASVs:", length(ecm_seqs), "| Saprotroph ASVs:", length(sapro_seqs), "\n")

# ---------- Helper: aggregate OTU to genus ----------
aggregate_to_genus <- function(otu_mat, asv_ids, tax_df, genus_column) {
  asv_genus <- tax_df[asv_ids, genus_column]
  names(asv_genus) <- asv_ids
  asv_genus[is.na(asv_genus) | asv_genus == "" | asv_genus == "unidentified"] <- "Unknown"
  asv_genus <- gsub("^[gG]__", "", asv_genus)
  asv_genus[asv_genus == ""] <- "Unknown"

  genera <- unique(asv_genus)
  genus_mat <- sapply(genera, function(g) {
    asvs <- names(asv_genus)[asv_genus == g]
    if (length(asvs) == 1) return(otu_mat[, asvs])
    rowSums(otu_mat[, asvs, drop = FALSE])
  })
  colnames(genus_mat) <- genera
  genus_mat
}

ecm_genus_mat <- aggregate_to_genus(otu_matched, ecm_seqs, its_tax, genus_col[1])
sapro_genus_mat <- aggregate_to_genus(otu_matched, sapro_seqs, its_tax, genus_col[1])
cat("ECM genera:", ncol(ecm_genus_mat), "| Saprotroph genera:", ncol(sapro_genus_mat), "\n")

# =============================================================================
# 3a. Depth-stratified ECM genus-environment correlations
# =============================================================================
cat("\n=== 3a. Depth-stratified ECM genus-env correlations ===\n")

# Top 20 named ECM genera
ranked_ecm <- sort(colSums(ecm_genus_mat), decreasing = TRUE)
named_ecm <- names(ranked_ecm)[names(ranked_ecm) != "Unknown"]
top20_ecm <- named_ecm[1:min(20, length(named_ecm))]

env_vars <- intersect(c("sand", "cmi_10_year", "x10_year_t_max_mean", "p_h", "em_colonization"), names(matched_sub))
cor_by_depth <- list()

for (depth_level in c("surface", "subsoil")) {
  idx <- matched_sub$depth_factor == depth_level
  for (g in top20_ecm) {
    for (ev in env_vars) {
      x <- matched_sub[[ev]][idx]
      y <- ecm_genus_mat[idx, g]
      if (sum(!is.na(x) & y > 0) >= 5) {
        ct <- cor.test(x, y, method = "spearman", exact = FALSE)
        cor_by_depth[[length(cor_by_depth) + 1]] <- data.frame(
          genus = g, env = ev, depth = depth_level,
          rho = round(ct$estimate, 3),
          p = signif(ct$p.value, 3),
          n = sum(!is.na(x)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

if (length(cor_by_depth) > 0) {
  cor_depth_df <- bind_rows(cor_by_depth)
  write_csv(cor_depth_df, file.path(TABLES_DIR, "ecm_genus_env_by_depth.csv"))
  cat("Saved: ecm_genus_env_by_depth.csv\n")

  sig <- cor_depth_df %>% filter(p < 0.05)
  cat("Significant (p<0.05):", nrow(sig), "of", nrow(cor_depth_df), "\n")

  # Summary: which genera have stronger subsoil correlations?
  wide <- cor_depth_df %>%
    tidyr::pivot_wider(names_from = depth, values_from = c(rho, p),
                       names_sep = "_")
  wide$abs_diff <- abs(wide$rho_subsoil) - abs(wide$rho_surface)
  stronger_sub <- wide %>% filter(abs_diff > 0)
  cat("Genus-env pairs with |rho| stronger in subsoil:", nrow(stronger_sub),
      "of", nrow(wide), "\n")
}

# =============================================================================
# 3b. Indicator genera by depth
# =============================================================================
cat("\n=== 3b. Indicator genera by depth ===\n")

# Aggregate every ITS ASV to genus for the general indicator analysis, so that
# pathogen-only and other fungal genera are included alongside ECM and
# saprotroph genera.
all_genus_mat <- aggregate_to_genus(
  otu_matched, colnames(otu_matched), its_tax, genus_col[1]
)

# Filter: present in >= 5 samples, exclude "Unknown"
prev <- colSums(all_genus_mat > 0)
keep <- prev >= 5 & colnames(all_genus_mat) != "Unknown"
genus_filt <- all_genus_mat[, keep, drop = FALSE]

# Relative abundance among all ITS reads. The denominator is calculated before
# removing unidentified and low-prevalence genera.
rel <- sweep(genus_filt, 1, rowSums(otu_matched), "/")
rel[is.nan(rel)] <- 0

grp <- factor(matched_sub$depth_factor)
cat("Running multipatt on", ncol(genus_filt), "genera,", nrow(genus_filt), "samples\n")

set.seed(42)
ind_genus <- tryCatch(
  multipatt(rel, grp, func = "IndVal.g", duleg = TRUE,
            control = how(blocks = factor(paste0(matched_sub$node, matched_sub$stand)), nperm = NPERM_IND)),
  error = function(e) { cat("Error:", e$message, "\n"); NULL }
)

if (!is.null(ind_genus)) {
  sign_tbl <- as.data.frame(ind_genus$sign) %>%
    tibble::rownames_to_column("genus") %>%
    mutate(p_bh = p.adjust(p.value, method = "BH"))

  group_names <- levels(grp)
  mp_cols <- paste0("s.", group_names)
  lookup <- setNames(group_names, mp_cols)
  gcols <- intersect(mp_cols, colnames(sign_tbl))

  if (length(gcols) > 0) {
    sign_tbl$indicator_for <- apply(sign_tbl[, gcols, drop = FALSE], 1, function(x) {
      paste(unname(lookup[gcols[x == 1]]), collapse = "+")
    })
  }

  sig_genera <- sign_tbl %>% filter(p.value < 0.05) %>% arrange(p.value)
  cat("Significant indicator genera:", nrow(sig_genera), "\n")
  if (nrow(sig_genera) > 0) {
    cat("Surface indicators:\n")
    print(sig_genera %>% filter(indicator_for == "surface") %>%
            select(genus, stat, p.value))
    cat("Subsoil indicators:\n")
    print(sig_genera %>% filter(indicator_for == "subsoil") %>%
            select(genus, stat, p.value))
  }

  write_csv(sign_tbl, file.path(TABLES_DIR, "indicator_genera_by_depth_full.csv"))
  write_csv(sig_genera, file.path(TABLES_DIR, "indicator_genera_by_depth.csv"))
  cat("Saved: indicator_genera_by_depth.csv\n")
}

# =============================================================================
# 3c. Genus depth-preference ratios
# =============================================================================
cat("\n=== 3c. Genus depth-preference ratios ===\n")

surf_idx <- matched_sub$depth_factor == "surface"
sub_idx  <- matched_sub$depth_factor == "subsoil"

# ECM genera
ecm_depth_pref <- data.frame(
  genus = colnames(ecm_genus_mat),
  guild = "ECM",
  mean_surface = colMeans(ecm_genus_mat[surf_idx, , drop = FALSE]),
  mean_subsoil = colMeans(ecm_genus_mat[sub_idx, , drop = FALSE]),
  stringsAsFactors = FALSE
) %>%
  filter(genus != "Unknown") %>%
  mutate(
    depth_ratio = mean_subsoil / (mean_surface + 1),
    preference = case_when(
      depth_ratio > 2   ~ "subsoil_specialist",
      depth_ratio < 0.5 ~ "surface_specialist",
      TRUE              ~ "generalist"
    )
  ) %>%
  arrange(desc(depth_ratio))

# Saprotroph genera
sapro_depth_pref <- data.frame(
  genus = colnames(sapro_genus_mat),
  guild = "Saprotroph",
  mean_surface = colMeans(sapro_genus_mat[surf_idx, , drop = FALSE]),
  mean_subsoil = colMeans(sapro_genus_mat[sub_idx, , drop = FALSE]),
  stringsAsFactors = FALSE
) %>%
  filter(genus != "Unknown") %>%
  mutate(
    depth_ratio = mean_subsoil / (mean_surface + 1),
    preference = case_when(
      depth_ratio > 2   ~ "subsoil_specialist",
      depth_ratio < 0.5 ~ "surface_specialist",
      TRUE              ~ "generalist"
    )
  ) %>%
  arrange(desc(depth_ratio))

depth_pref <- bind_rows(ecm_depth_pref, sapro_depth_pref)
write_csv(depth_pref, file.path(TABLES_DIR, "genus_depth_preferences.csv"))
cat("Saved: genus_depth_preferences.csv\n")

cat("\nECM depth specialists:\n")
print(ecm_depth_pref %>% filter(preference != "generalist") %>%
        select(genus, mean_surface, mean_subsoil, depth_ratio, preference) %>%
        head(20))

cat("\nSaprotroph depth specialists:\n")
print(sapro_depth_pref %>% filter(preference != "generalist") %>%
        select(genus, mean_surface, mean_subsoil, depth_ratio, preference) %>%
        head(20))

cat("\n=== 03b_indicators_its_genus.R complete ===\n")

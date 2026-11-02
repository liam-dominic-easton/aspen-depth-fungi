## 03_indicators_its.R — ITS (fungi) ASV-level indicator species analysis.
##
## Outputs (in outputs/tables/):
##   indicator_by_depth.csv, indicator_by_depth_full.csv  — ASV depth indicators
##                                                          (the full table is read by 15)
##   indicator_by_biome.csv, indicator_by_biome_full.csv  — biome (boreal vs parkland) indicators
##   indicator_by_depth_biome.csv                         — depth x biome indicators
##   ecm_genus_env_correlations.csv                       — top ECM genera vs environment
##   saprotroph_by_depth.csv                              — saprotroph genera by depth
##
## Genus-level depth indicators (Table 2) are computed in 03b_indicators_its_genus.R.
## Indicators are computed with indicspecies::multipatt at IndVal.g on
## relative-abundance data; depth labels are permuted within stands
## (NPERM_IND permutations). Taxa with prevalence < 5 samples are dropped
## before testing.

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
cat("Taxonomy columns:", paste(names(its_tax), collapse = ", "), "\n")

# Load FUNGuild for guild assignments
guild_path <- file.path(TAX_DIR, "otu_table_for_FUNGuild.guilds.txt")
if (file.exists(guild_path)) {
  funguild <- read.delim(guild_path, stringsAsFactors = FALSE)
  cat("FUNGuild loaded:", nrow(funguild), "rows\n")
  cat("FUNGuild columns:", paste(names(funguild), collapse = ", "), "\n")
} else {
  cat("WARNING: FUNGuild file not found\n")
  funguild <- NULL
}

# ---------- Helper: Indicator species ----------
run_indicator <- function(otu_mat, grouping, label, min_prev = 5,
                          control = how(nperm = NPERM_IND)) {
  valid <- !is.na(grouping) & grouping != ""
  otu_sub <- otu_mat[valid, , drop = FALSE]
  grp <- factor(grouping[valid])

  prev <- colSums(otu_sub > 0)
  otu_sub <- otu_sub[, prev >= min_prev, drop = FALSE]

  if (ncol(otu_sub) < 2) { cat(glue("[{label}] Too few taxa"), "\n"); return(NULL) }

  rel <- sweep(otu_sub, 1, rowSums(otu_sub), "/")
  rel[is.nan(as.matrix(rel))] <- 0

  cat(glue("[{label}] {nrow(rel)} samples, {ncol(rel)} taxa, {nlevels(grp)} groups"), "\n")

  set.seed(42)
  ind <- tryCatch(
    multipatt(rel, grp, func = "IndVal.g", duleg = TRUE, control = control),
    error = function(e) { cat(glue("[{label}] Error: {e$message}"), "\n"); NULL }
  )
  if (is.null(ind)) return(NULL)

  sign_tbl <- as.data.frame(ind$sign) %>%
    tibble::rownames_to_column("taxon") %>%
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

  # Add taxonomy
  tax_info <- its_tax %>% tibble::rownames_to_column("taxon")
  sign_tbl <- sign_tbl %>% left_join(tax_info, by = "taxon")

  sig <- sign_tbl %>% filter(p.value < 0.05) %>% arrange(p.value)
  cat(glue("[{label}] Significant: {nrow(sig)}/{nrow(sign_tbl)}"), "\n")

  list(full = sign_tbl, significant = sig, model = ind, rel = rel, valid = valid)
}

# ---------- Subset OTU to matched samples ----------
matched_names <- matched$sample_name
otu_matched <- its_otu[rownames(its_otu) %in% matched_names, ]
otu_matched <- otu_matched[match(matched$sample_name[matched$sample_name %in% rownames(otu_matched)], rownames(otu_matched)), ]
matched_sub <- matched %>% filter(sample_name %in% rownames(otu_matched))

# =============================================================================
# 1. Indicator species by DEPTH
# =============================================================================
cat("\n=== Indicators by DEPTH ===\n")
## depth labels permuted within stands (paired design), NPERM_IND permutations
ind_depth <- run_indicator(otu_matched, as.character(matched_sub$depth_factor), "depth",
                           control = how(blocks = factor(paste0(matched_sub$node, matched_sub$stand)),
                                         nperm = NPERM_IND))
if (!is.null(ind_depth)) {
  write_csv(ind_depth$full, file.path(TABLES_DIR, "indicator_by_depth_full.csv"))
  write_csv(ind_depth$significant, file.path(TABLES_DIR, "indicator_by_depth.csv"))
  cat("Saved: indicator_by_depth.csv\n")

  for (d in c("surface", "subsoil")) {
    n <- sum(ind_depth$significant$indicator_for == d)
    cat(glue("  {d} indicators: {n}"), "\n")
  }
}

# =============================================================================
# 2. Indicator species by BIOME
# =============================================================================
cat("\n=== Indicators by BIOME ===\n")
## biome is node-level: exact node-level null over 126 assignments (not reported
## in the manuscript; kept for completeness)
ind_biome <- run_indicator(otu_matched, matched_sub$biome, "biome", control = how(nperm = 99))
if (!is.null(ind_biome)) {
  ex <- indval_biome_exact(ind_biome$rel, matched_sub$node[ind_biome$valid], matched_sub$biome[ind_biome$valid])
  ind_biome$full <- ind_biome$full |> left_join(ex, by = "taxon") |>
    mutate(p_node_bh = p.adjust(p_node_exact, "BH"))
  ind_biome$significant <- ind_biome$full |> filter(p_node_exact < 0.05) |> arrange(p_node_exact)
  write_csv(ind_biome$full, file.path(TABLES_DIR, "indicator_by_biome_full.csv"))
  write_csv(ind_biome$significant, file.path(TABLES_DIR, "indicator_by_biome.csv"))
  cat("Saved: indicator_by_biome.csv\n")
}

# =============================================================================
# 3. Indicator species by DEPTH × BIOME
# =============================================================================
cat("\n=== Indicators by DEPTH × BIOME ===\n")
depth_biome <- paste(matched_sub$depth_factor, matched_sub$biome, sep = "_")
## not reported in the manuscript (free permutations; biome part pseudoreplicated)
ind_db <- run_indicator(otu_matched, depth_biome, "depth×biome", control = how(nperm = 999))
if (!is.null(ind_db)) {
  write_csv(ind_db$significant, file.path(TABLES_DIR, "indicator_by_depth_biome.csv"))
  cat("Saved: indicator_by_depth_biome.csv\n")
}

# =============================================================================
# 4. ECM genus analysis (using FUNGuild)
# =============================================================================
cat("\n=== ECM genus analysis ===\n")

if (!is.null(funguild)) {
  # Find ECM taxa in FUNGuild output
  guild_col <- intersect(c("Guild", "guild"), names(funguild))
  if (length(guild_col) > 0) {
    gcol <- guild_col[1]
    ecm_rows <- funguild %>%
      filter(grepl("Ectomycorrhiz", !!sym(gcol), ignore.case = TRUE))

    # The FUNGuild Taxon column should match OTU sequences
    taxon_col <- intersect(c("Taxon", "taxonomy", "OTU_ID"), names(funguild))
    cat("ECM rows in FUNGuild:", nrow(ecm_rows), "\n")

    # Match ECM ASVs to OTU table via OTU_ID (full sequences)
    ecm_seqs <- intersect(ecm_rows[["OTU_ID"]], colnames(otu_matched))

    cat("ECM ASVs matched to OTU table:", length(ecm_seqs), "\n")

    if (length(ecm_seqs) > 5) {
      ecm_otu <- otu_matched[, ecm_seqs, drop = FALSE]

      # Map to genus from taxonomy table
      genus_col <- intersect(c("Genus", "genus"), names(its_tax))
      if (length(genus_col) > 0) {
        asv_genus <- its_tax[ecm_seqs, genus_col[1]]
        names(asv_genus) <- ecm_seqs

        # Clean: replace NA/empty with "Unknown" and strip FUNGuild prefixes
        asv_genus[is.na(asv_genus) | asv_genus == "" | asv_genus == "unidentified"] <- "Unknown"
        asv_genus <- gsub("^[gG]__", "", asv_genus)
        asv_genus[asv_genus == ""] <- "Unknown"

        cat("ECM genus assignments:", length(unique(asv_genus)), "genera (",
            sum(asv_genus == "Unknown"), "unknown ASVs)\n")

        # Aggregate by genus
        genera <- unique(asv_genus)
        genus_mat <- sapply(genera, function(g) {
          asvs <- names(asv_genus)[asv_genus == g]
          if (length(asvs) == 1) return(ecm_otu[, asvs])
          rowSums(ecm_otu[, asvs, drop = FALSE])
        })

        # Top 20 by abundance (exclude Unknown for ranking, keep for completeness)
        ranked <- sort(colSums(genus_mat), decreasing = TRUE)
        named_genera <- names(ranked)[names(ranked) != "Unknown"]
        top20 <- named_genera[1:min(20, length(named_genera))]
        cat("Top 20 ECM genera:", paste(top20, collapse = ", "), "\n")

        # Environmental correlations
        env_vars <- intersect(c("sand", "latitude", "p_h", "em_colonization"), names(matched_sub))
        cor_results <- list()

        for (g in top20) {
          for (ev in env_vars) {
            x <- matched_sub[[ev]]
            y <- genus_mat[, g]
            if (sum(!is.na(x) & y > 0) >= 5) {
              ct <- cor.test(x, y, method = "spearman", exact = FALSE)
              cor_results[[length(cor_results) + 1]] <- data.frame(
                genus = g, env = ev, rho = round(ct$estimate, 3),
                p = signif(ct$p.value, 3), stringsAsFactors = FALSE
              )
            }
          }
        }

        if (length(cor_results) > 0) {
          cor_df <- bind_rows(cor_results)
          write_csv(cor_df, file.path(TABLES_DIR, "ecm_genus_env_correlations.csv"))
          cat("Saved: ecm_genus_env_correlations.csv\n")
          sig <- cor_df %>% filter(p < 0.05)
          if (nrow(sig) > 0) { cat("Significant:\n"); print(sig) }
        }
      }
    }
  }
}

# =============================================================================
# 5. Saprotroph depth analysis
# =============================================================================
cat("\n=== Saprotroph by depth ===\n")

if (!is.null(funguild) && length(guild_col) > 0) {
  sapro_rows <- funguild %>%
    filter(grepl("Saprotroph", !!sym(gcol), ignore.case = TRUE))

  sapro_seqs <- intersect(sapro_rows[["OTU_ID"]], colnames(otu_matched))
  cat("Saprotroph ASVs matched:", length(sapro_seqs), "\n")

  if (length(sapro_seqs) > 5) {
    genus_col <- intersect(c("Genus", "genus"), names(its_tax))
    if (length(genus_col) > 0) {
      sapro_otu <- otu_matched[, sapro_seqs, drop = FALSE]
      asv_genus <- its_tax[sapro_seqs, genus_col[1]]
      names(asv_genus) <- sapro_seqs

      # Clean NA/empty genus
      asv_genus[is.na(asv_genus) | asv_genus == "" | asv_genus == "unidentified"] <- "Unknown"
      asv_genus <- gsub("^[gG]__", "", asv_genus)
      asv_genus[asv_genus == ""] <- "Unknown"

      cat("Saprotroph genus assignments:", length(unique(asv_genus)), "genera (",
          sum(asv_genus == "Unknown"), "unknown ASVs)\n")

      genera <- unique(asv_genus)
      genus_mat <- sapply(genera, function(g) {
        asvs <- names(asv_genus)[asv_genus == g]
        if (length(asvs) == 1) return(sapro_otu[, asvs])
        rowSums(sapro_otu[, asvs, drop = FALSE])
      })

      surf_idx <- matched_sub$depth_factor == "surface"
      sub_idx  <- matched_sub$depth_factor == "subsoil"

      sapro_depth <- data.frame(
        genus = colnames(genus_mat),
        mean_surface = colMeans(genus_mat[surf_idx, , drop = FALSE]),
        mean_subsoil = colMeans(genus_mat[sub_idx, , drop = FALSE]),
        stringsAsFactors = FALSE
      ) %>%
        mutate(total = mean_surface + mean_subsoil) %>%
        arrange(desc(total))

      write_csv(sapro_depth, file.path(TABLES_DIR, "saprotroph_by_depth.csv"))
      cat("Saved: saprotroph_by_depth.csv\n")
      cat("Top 10:\n")
      print(head(sapro_depth, 10))
    }
  }
}

cat("\n=== 03_indicators_its.R complete ===\n")

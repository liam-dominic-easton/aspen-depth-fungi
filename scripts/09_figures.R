## 09_figures.R - Main-text Figures 1, 3, 5 and 6 for the manuscript.
##
## Inputs:
##   outputs/tables/bayesian_depth_comparisons.csv       (Figure 5; from 06)
##   outputs/tables/root_biomass_bayesian_models.csv     (Figure 6; from 07)
##   data/phyloseq/ps_*.rds                              (Figure 3)
##   data/taxonomy/maarjam_vsearch_mapping_final.csv     (Figure 3, 18S panel)
##   data/environmental/depth_analysis_data.csv          (Figures 1 and 3)
##   data/reference/natural_earth/ne_50m_admin_1_states_provinces_lines.{shp,dbf}
##                                                        (Figure 1A provincial boundaries)
##
## Outputs (in outputs/figures/):
##   study_area_climate.png              - Figure 1
##   pcoa_with_envfit.png                - Figure 3
##   depth_comparison_forest.png         - Figure 5
##   root_biomass_model_comparison.png   - Figure 6
## (Figure 2: 18_latitude_figure.R; Figure 4: 20_paired_stand_figure.R; the
##  variance-partitioning diagram, Figure S1, is made by 12_guild_balance_roots_varpart.R.)
##
## Figures 5 and 6 read the model tables written by
## scripts/06_bayesian_depth_comparisons.R and scripts/07_root_biomass_bayesian.R,
## so run those scripts first.

source("scripts/helpers.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(maps)
})

# ============================================================================
# 0. Figure 1: study area and supplied CFS/CIPHA climate gradient
# ============================================================================
env_fig <- load_env() |>
  distinct(node_stand, .keep_all = TRUE) |>
  mutate(biome = factor(biome, levels = c("boreal", "parkland")))

na_map <- map_data("world", region = c("Canada", "USA"))  # country outlines for context
biome_cols <- c(boreal = "#2166AC", parkland = "#B2182B")

## Canadian provincial/territorial boundaries: Natural Earth 1:50m
## admin-1 boundary lines, v5.1.0 (public domain; naturalearthdata.com, file
## ne_50m_admin_1_states_provinces_lines.zip from the NACIS CDN, naciscdn.org).
## Unzip into data/reference/natural_earth/. Read with base R, so sf is not
## required: .dbf for the country code, .shp for the line coordinates.
NE_ADMIN1 <- file.path(DATA_DIR, "reference", "natural_earth",
                       "ne_50m_admin_1_states_provinces_lines")
read_dbf_field <- function(path, field) {
  con <- file(path, "rb"); on.exit(close(con))
  le <- function(b) readBin(b, "integer", size = length(b), endian = "little",
                            signed = length(b) > 2)
  hdr <- readBin(con, "raw", 32)
  n_rec <- le(hdr[5:8]); hdr_len <- le(hdr[9:10]); rec_len <- le(hdr[11:12])
  n_fld <- (hdr_len - 33) %/% 32
  fd <- matrix(readBin(con, "raw", hdr_len - 32)[seq_len(32 * n_fld)], nrow = 32)
  nm <- apply(fd[1:11, , drop = FALSE], 2, function(b) rawToChar(b[b != as.raw(0)]))
  len <- as.integer(fd[17, ]); start <- 2 + c(0, cumsum(len))  # byte 1 = deletion flag
  j <- match(field, nm)
  rec <- matrix(readBin(con, "raw", n_rec * rec_len), nrow = rec_len)
  trimws(apply(rec[start[j]:(start[j] + len[j] - 1), , drop = FALSE], 2, rawToChar))
}
read_shp_lines <- function(path) {   # polyline shapefile -> long/lat per record and part
  con <- file(path, "rb"); on.exit(close(con))
  readBin(con, "raw", 100)
  out <- list()
  repeat {
    h <- readBin(con, "integer", 2, size = 4, endian = "big")
    if (length(h) < 2) break
    b <- readBin(con, "raw", 2 * h[2])
    i4 <- function(from, n = 1) readBin(b[from:(from + 4 * n - 1)], "integer", n,
                                        size = 4, endian = "little")
    if (i4(1) == 0) next                                    # null shape
    np <- i4(37); npt <- i4(41); parts <- i4(45, np)
    xy <- matrix(readBin(b[(45 + 4 * np):length(b)], "double", 2 * npt, size = 8,
                         endian = "little"), ncol = 2, byrow = TRUE)
    out[[length(out) + 1]] <- data.frame(record = h[1], part = findInterval(seq_len(npt) - 1, parts),
                                         long = xy[, 1], lat = xy[, 2])
  }
  do.call(rbind, out)
}
have_ne <- file.exists(paste0(NE_ADMIN1, c(".shp", ".dbf")))
if (all(have_ne)) {
  ne_country <- read_dbf_field(paste0(NE_ADMIN1, ".dbf"), "ADM0_A3")
  prov_lines <- read_shp_lines(paste0(NE_ADMIN1, ".shp")) |>
    filter(record %in% which(ne_country == "CAN")) |>
    mutate(group = paste(record, part))
}

## label positions chosen by hand to sit inside each unit and clear of the stands
prov_labs <- tibble::tribble(
  ~long,  ~lat,  ~label,
  -125.2, 54.6,  "British\nColumbia",
  -113.6, 59.2,  "Alberta",
  -105.9, 56.0,  "Saskatchewan",
   -98.0, 55.0,  "Manitoba",
   -89.6, 50.4,  "Ontario")
country_labs <- tibble::tribble(
  ~long,  ~lat,  ~label,
   -89.9, 53.7,  "CANADA",
  -108.0, 45.0,  "USA")

p_map <- ggplot(na_map, aes(long, lat, group = group)) +
  geom_polygon(fill = "grey95", color = "grey50", linewidth = 0.3)
if (all(have_ne)) {
  p_map <- p_map + geom_path(data = prov_lines, color = "grey60", linewidth = 0.25)
}
p_map <- p_map +
  geom_text(data = country_labs, aes(long, lat, label = label), inherit.aes = FALSE,
            size = 3.4, fontface = "bold", color = "grey35") +
  geom_text(data = prov_labs, aes(long, lat, label = label), inherit.aes = FALSE,
            size = 2.8, fontface = "italic", color = "grey35", lineheight = 0.85) +
  geom_point(data = env_fig, aes(longitude, latitude, color = biome),
             inherit.aes = FALSE, size = 2.2) +
  coord_quickmap(xlim = c(-135, -85), ylim = c(41, 62), expand = FALSE) +
  scale_color_manual(values = biome_cols) +
  labs(title = "A", x = "Longitude", y = "Latitude", color = "Biome") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), panel.grid = element_blank(),
        ## biome key inside the map (Pacific, lower left) so the map can use
        ## the full figure height; it also serves panels B and C
        legend.position = "inside", legend.position.inside = c(0.015, 0.025),
        legend.justification.inside = c(0, 0),
        legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.3),
        legend.key.size = unit(0.4, "cm"), legend.margin = margin(3, 5, 3, 4))

## B and C: confidence band and line drawn first, opaque points on top, so a
## point's colour does not change where it overlaps the band
p_cmi <- ggplot(env_fig, aes(latitude, cmi_10_year)) +
  geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey25",
              fill = "grey80", linewidth = 0.7) +
  geom_point(aes(color = biome), size = 2.2) +
  scale_color_manual(values = biome_cols, guide = "none") +
  labs(title = "B", x = "Latitude", y = "Climate moisture index (cm per year)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

p_mat <- ggplot(env_fig, aes(latitude, mat_10_year)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey25",
              fill = "grey80", linewidth = 0.7) +
  geom_point(aes(color = biome), size = 2.2) +
  scale_color_manual(values = biome_cols, guide = "none") +
  labs(title = "C", x = "Latitude", y = "Mean annual temperature (\u00B0C)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

p_fig1 <- p_map + p_cmi + p_mat + plot_layout(widths = c(2.4, 1, 1))
if (all(have_ne)) {
  ggsave(file.path(FIGS_DIR, "study_area_climate.png"), p_fig1,
         width = 10.5, height = 4.4, dpi = 300)
  cat("Saved: outputs/figures/study_area_climate.png\n")
} else {
  warning("Natural Earth admin-1 lines not found at ", NE_ADMIN1,
          ".{shp,dbf}; Figure 1 not saved (existing PNG left unchanged).")
}

# ============================================================================
# 1. Figure 5: subsoil-minus-surface differences (Bayesian paired models, 06)
# ============================================================================
# contrast_* in the table = surface - subsoil on the response scale, computed
# per posterior draw for a typical stand (population level). It is plotted with
# the sign reversed, as subsoil - surface, to match Figure 6 (depth coefficient,
# subsoil relative to surface): the difference is negated and the 90% HDI
# becomes [-hi, -lo]; Pr(direction) is unchanged. Proportions in percentage
# points; log(saprotroph:ECM) on the log scale; richness in taxa; roots in g
# and ECM tips per m. Intervals are 90% HDIs; black points Pr(direction) >= 0.95, grey points lower.
# Pathogen rows are omitted (see guild-ambiguity tables).
cat("\n=== Figure 5: guild depth comparisons ===\n")
depth_df <- read_csv(file.path(TABLES_DIR, "bayesian_depth_comparisons.csv"), show_col_types = FALSE)
panel_of <- c("ECM proportion" = "Guild shares (percentage points)",
              "Saprotroph proportion" = "Guild shares (percentage points)",
              "Saprotroph:ECM ratio (log)" = "Guild balance (log ratio)",
              "ECM richness" = "Richness (taxa)", "Saprotroph richness" = "Richness (taxa)",
              "AMF richness" = "Richness (taxa)",
              "Fine root biomass" = "Fine-root biomass (g)", "ECM tips per m root" = "ECM tips per m root")
## display labels matching the manuscript text
display_of <- c("ECM proportion" = "ECM share", "Saprotroph proportion" = "Saprotroph share",
                "Saprotroph:ECM ratio (log)" = "log(saprotroph:ECM)", "AMF richness" = "AM fungal richness")
relabel <- function(x) ifelse(x %in% names(display_of), unname(display_of[x]), x)
bt <- depth_df |>
  filter(outcome %in% names(panel_of)) |>
  mutate(panel = unname(panel_of[outcome]), outcome = relabel(outcome),
         bold = effect_pr_direction >= 0.95,
         pr_label = sprintf("Pr = %.2f", effect_pr_direction),
         diff_mean = -contrast_mean,          # subsoil - surface
         diff_lo   = -contrast_hdi_hi,
         diff_hi   = -contrast_hdi_lo)
print(bt |> select(outcome, diff_mean, diff_lo, diff_hi, effect_pr_direction))

make_panel <- function(data, title) {
  data <- data |> mutate(outcome = factor(outcome, levels = rev(unique(outcome))))
  ggplot(data, aes(x = diff_mean, y = outcome)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = diff_lo, xmax = diff_hi),
                  width = 0.3, linewidth = 0.6, orientation = "y") +
    geom_point(aes(fill = bold), shape = 21, size = 3) +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "grey60"), guide = "none") +
    labs(title = title, x = "Difference, subsoil \u2212 surface", y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          panel.grid.minor = element_blank(), plot.margin = margin(5, 60, 5, 5)) +
    coord_cartesian(clip = "off") +
    geom_text(aes(label = pr_label), x = Inf, hjust = -0.05, size = 2.8, color = "grey30")
}
panels <- unique(unname(panel_of))
heights <- sapply(panels, function(pn) max(1, sum(bt$panel == pn)))
p_depth <- wrap_plots(lapply(panels, function(pn) make_panel(filter(bt, panel == pn), pn)), ncol = 1) +
  plot_layout(heights = heights, axis_titles = "collect") +   # one shared x-axis title
  plot_annotation(caption = paste("Subsoil minus surface for a typical stand; positive = higher in the subsoil.",
                                  "90% HDI; black: Pr(direction) ≥ 0.95; grey: lower."))
ggsave(file.path(FIGS_DIR, "depth_comparison_forest.png"), p_depth, width = 7.5, height = 9.5, dpi = 300)
cat("Saved: outputs/figures/depth_comparison_forest.png\n")

# ============================================================================
# 2. Figure 6: depth effect and within-depth root slope (07)
# ============================================================================
# Depth and fine-root biomass do not overlap between depths, so the depth
# coefficient conditional on roots is not identifiable and is not plotted.
# Rows: depth effect from the depth-only model (subsoil - surface), and the
# root slope from the depth + roots model (per SD of log fine-root biomass,
# i.e. among samples from the same depth). Outcomes standardized (SD units).
cat("\n=== Figure 6: root-biomass models ===\n")
root_df <- read_csv(file.path(TABLES_DIR, "root_biomass_bayesian_models.csv"), show_col_types = FALSE)
make_root_panel <- function(r, title) {
  pd <- tibble(predictor = c("Depth: subsoil vs surface (depth-only model)",
                             "Root biomass within depth (depth + roots model)"),
               est = c(r$m1_depth_mean, r$m3_roots_within_depth_mean),
               lo = c(r$m1_depth_hdi_lo, r$m3_roots_within_depth_hdi_lo),
               hi = c(r$m1_depth_hdi_hi, r$m3_roots_within_depth_hdi_hi),
               pr = c(r$m1_depth_pr_direction, r$m3_roots_within_depth_pr_direction)) |>
    mutate(bold = pr >= 0.95, pr_label = sprintf("Pr = %.2f", pr),
           predictor = factor(predictor, levels = rev(predictor)))
  ggplot(pd, aes(x = est, y = predictor)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.6, orientation = "y") +
    geom_point(aes(fill = bold), shape = 21, size = 3) +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "grey60"), guide = "none") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          panel.grid.minor = element_blank(), plot.margin = margin(5, 80, 5, 5)) +
    coord_cartesian(clip = "off") +
    geom_text(aes(label = pr_label), x = Inf, hjust = -0.05, size = 2.8, color = "grey30")
}
ord <- c("ECM proportion", "Saprotroph proportion", "AMF richness")   # same order as Figure 5
p_root <- wrap_plots(lapply(ord, function(o) make_root_panel(filter(root_df, outcome == o), relabel(o))), ncol = 1) +
  plot_annotation(caption = paste0("Coefficients on the standardized outcome scale, 90% HDI; black: Pr(direction) \u2265 0.95; grey: lower.\n",
                                   "The depth effect adjusted for roots is omitted because it relies on extrapolation (roots do not overlap between depths)."),
                  theme = theme(plot.caption = element_text(hjust = 0)))
ggsave(file.path(FIGS_DIR, "root_biomass_model_comparison.png"), p_root, width = 7, height = 6, dpi = 300)
cat("Saved: outputs/figures/root_biomass_model_comparison.png\n")

# ============================================================================
# 3. Combined PCoA + envfit (Figure 3)
# ============================================================================

cat("\n=== Combined PCoA + envfit figure ===\n")

ps_its <- load_ps("ITS")
ps_18s <- load_ps("18S")
env    <- load_env()

libraries <- list(
  its  = list(ps = ps_its, label = "ITS (fungi)"),
  amf  = list(ps = ps_18s, label = "AM fungi (18S)")
)

# Build VT-level matrix for 18S PCoA
cat("Building VT abundance matrix for 18S PCoA...\n")
vt_mat <- build_vt_matrix(ps_18s)

envfit_vars <- c("sand", "cmi_10_year", "mat_10_year", "p_h")
envfit_labels <- c(sand = "Sand", cmi_10_year = "CMI", mat_10_year = "MAT", p_h = "pH")

pcoa_plots <- list()
for (lib_name in names(libraries)) {
  lib <- libraries[[lib_name]]

  if (lib_name == "amf") {
    # VT-level PCoA for 18S
    vt_dist <- aitchison_distance_matrix(vt_mat)
    pcoa_obj <- ape::pcoa(vt_dist)
  } else {
    res <- aitchison_pcoa(lib$ps)
    pcoa_obj <- res$ordination
  }
  md <- match_ps_to_env(lib$ps, env)

  scores <- pcoa_obj$vectors[, 1:2]
  colnames(scores) <- c("PCo1", "PCo2")
  var_exp <- round(pcoa_obj$values$Relative_eig[1:2] * 100, 1)

  plot_data <- cbind(scores, md[match(rownames(scores), md$sample_name), ])

  # Envfit
  env_matrix <- plot_data[, envfit_vars]
  cc <- complete.cases(env_matrix)
  ## Vector directions and r2 from envfit; p-values from design-based
  ## permutation: CMI and MAT permuted among nodes, sand and pH among whole
  ## stands.
  ef <- vegan::envfit(scores[cc, ], env_matrix[cc, ], permutations = 0)
  md_cc <- plot_data[cc, ] |> mutate(stand_id = paste0(node, stand))
  sc_cc <- scores[cc, ]
  r2_of <- function(v) summary(lm(v ~ sc_cc[, 1] + sc_cc[, 2]))$r.squared
  design_p <- function(var, unit) {
    u <- if (unit == "node") md_cc$node else md_cc$stand_id
    x <- md_cc[[var]]
    if (unit == "node") x <- ave(x, u)            # node means
    obs <- r2_of(x)
    units <- sort(unique(u))
    grp <- if (unit == "node") rep("all", length(units)) else
      sapply(split(md_cc$depth, u), function(z) paste(sort(unique(z)), collapse = ""))[units]
    key <- paste(u, md_cc$depth); val <- setNames(x, key)
    set.seed(42)
    null <- replicate(NPERM, {
      map <- setNames(units, units)
      for (g in unique(grp)) { i <- which(grp == g); map[i] <- units[i][sample.int(length(i))] }
      xp <- if (unit == "node") unname(tapply(x, u, mean)[map[u]]) else unname(val[paste(map[u], md_cc$depth)])
      r2_of(xp)
    })
    (sum(null >= obs - 1e-12) + 1) / (NPERM + 1)
  }
  unit_of <- c(sand = "node", p_h = "node", cmi_10_year = "node", mat_10_year = "node")

  ## unit arrows x sqrt(r2), so arrow length is proportional to sqrt(r2)
  ## (scores(ef) would already include sqrt(r2), so the unit arrows are used)
  arrow_coords <- as.data.frame(ef$vectors$arrows)
  colnames(arrow_coords) <- c("PCo1", "PCo2")
  arrow_coords$var <- rownames(arrow_coords)
  arrow_coords$r2 <- ef$vectors$r
  arrow_coords$pval <- sapply(arrow_coords$var, function(v) design_p(v, unit_of[[v]]))
  cat(lib$label, "envfit (design-based p):\n"); print(arrow_coords[, c("var", "r2", "pval")])

  # Scale arrows for readability, then place each label just beyond its tip so
  # the label box cannot obscure the arrowhead.
  axis_range <- max(abs(range(scores[, 1])), abs(range(scores[, 2])))
  arrow_scale <- axis_range * 0.45

  sig_arrows <- arrow_coords %>%
    filter(pval < 0.05) %>%
    mutate(
      xend = PCo1 * sqrt(r2) * arrow_scale,
      yend = PCo2 * sqrt(r2) * arrow_scale,
      arrow_length = pmax(sqrt(xend^2 + yend^2), .Machine$double.eps),
      label_offset = axis_range * 0.055,
      label_x = xend + (xend / arrow_length) * label_offset,
      label_y = yend + (yend / arrow_length) * label_offset,
      label = envfit_labels[var]
    )

  p <- ggplot(plot_data, aes(x = PCo1, y = PCo2)) +
    stat_ellipse(aes(group = depth_factor, linetype = depth_factor),
                 type = "norm", level = 0.95, color = "grey20",
                 linewidth = 0.75,
                 show.legend = c(linetype = TRUE, color = FALSE, shape = FALSE)) +
    geom_point(aes(color = biome, shape = depth_factor), size = 3, alpha = 0.48,
               show.legend = c(color = TRUE, shape = TRUE, linetype = FALSE)) +
    scale_color_manual(values = c("boreal" = "#2166ac", "parkland" = "#b2182b"),
                       name = "Biome") +
    scale_shape_manual(values = c("surface" = 16, "subsoil" = 17),
                       name = "Depth") +
    scale_linetype_manual(values = c("surface" = "solid", "subsoil" = "dashed"),
                          name = "Depth ellipse") +
    guides(
      color = guide_legend(override.aes = list(linetype = 0, linewidth = 0)),
      shape = guide_legend(override.aes = list(linetype = 0, linewidth = 0)),
      linetype = guide_legend(override.aes = list(shape = NA))
    ) +
    labs(
      title = lib$label,
      x = paste0("PCo1 (", var_exp[1], "%)"),
      y = paste0("PCo2 (", var_exp[2], "%)")
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"))

  if (nrow(sig_arrows) > 0) {
    p <- p +
      geom_segment(data = sig_arrows,
                   aes(x = 0, y = 0, xend = xend, yend = yend),
                   arrow = arrow(length = unit(0.12, "cm"),
                                 angle = 25, type = "open"),
                   linewidth = 1.05, color = "black", inherit.aes = FALSE) +
      geom_label_repel(data = sig_arrows,
                       aes(x = label_x, y = label_y, label = label),
                       size = 3.5, color = "black", fontface = "bold",
                       fill = scales::alpha("white", 0.85), label.size = 0.15,
                       point.padding = 0.15, box.padding = 0.45,
                       min.segment.length = 0, segment.color = "grey45",
                       inherit.aes = FALSE)
  }

  pcoa_plots[[lib_name]] <- p
}

combined_pcoa <- pcoa_plots[["its"]] + pcoa_plots[["amf"]] +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(FIGS_DIR, "pcoa_with_envfit.png"),
       combined_pcoa, width = 8.5, height = 5.2, dpi = 300)
cat("Saved: outputs/figures/pcoa_with_envfit.png\n")

cat("\nDone.\n")

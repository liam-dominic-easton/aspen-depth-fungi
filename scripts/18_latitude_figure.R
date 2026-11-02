## 18_latitude_figure.R - Figure 2: guild balance, ECM share and AMF richness
## against latitude at each depth.
##
## Descriptive display of the data behind predictions 1 and 3. Each point is a
## node mean at one depth and each bar spans the stands in that node, so the
## display shows the nine nodes that latitude varies among; surface and subsoil
## are offset slightly. No fitted lines, so the display does not imply more
## certainty than the Bayesian latitude models of 16_latitude_effects.R support.
##
## Inputs:  load_env(); outputs/tables/guild_ambiguity_sample_proportions.csv
## Output:  outputs/figures/latitude_by_depth.png

source("scripts/helpers.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(tidyr) })

env <- load_env() |>
  mutate(depth_label = factor(ifelse(depth == "A", "Surface (0–20 cm)", "Subsoil (20–40 cm)"),
                              levels = c("Surface (0–20 cm)", "Subsoil (20–40 cm)")))

lsr <- read_csv(file.path(TABLES_DIR, "guild_ambiguity_sample_proportions.csv"), show_col_types = FALSE) |>
  filter(method == "priority") |>
  select(sample_id_clean, guild, guild_reads) |>
  pivot_wider(names_from = guild, values_from = guild_reads) |>
  transmute(sample_id_clean, lsr = log((Saprotroph + 0.5) / (ECM + 0.5)))
env <- env |> left_join(lsr, by = "sample_id_clean")
stopifnot(!anyNA(env$lsr), !anyNA(env$prop_ecm), !anyNA(env$vt_richness))

depth_cols <- c("Surface (0–20 cm)" = "#8C510A", "Subsoil (20–40 cm)" = "#01665E")
nodes <- env |>
  pivot_longer(c(lsr, prop_ecm, vt_richness), names_to = "var", values_to = "y") |>
  group_by(node, depth_label, var) |>
  summarise(lat = mean(latitude), mean = mean(y), lo = min(y), hi = max(y), .groups = "drop") |>
  mutate(lat = lat + ifelse(depth_label == "Surface (0–20 cm)", -0.12, 0.12))   # offset depths
panel <- function(v, ylab, tag) {
  ggplot(filter(nodes, var == v), aes(lat, mean, colour = depth_label)) +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.35, alpha = 0.3) +   # faint, so the means stand out
    geom_point(size = 2.3) +
    scale_colour_manual(values = depth_cols, name = NULL) +
    labs(x = "Latitude (°N)", y = ylab, tag = tag) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
}
p <- panel("lsr", "log(saprotroph:ECM reads)", "A") +
  panel("prop_ecm", "ECM share of ITS reads", "B") + scale_y_continuous(labels = scales::percent) +
  panel("vt_richness", "AM fungal richness (VT per sample)", "C") +
  plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(FIGS_DIR, "latitude_by_depth.png"), p, width = 9, height = 3.6, dpi = 300)
cat("Saved: latitude_by_depth.png\n18_latitude_figure.R complete.\n")

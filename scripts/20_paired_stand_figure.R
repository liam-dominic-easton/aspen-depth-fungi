## 20_paired_stand_figure.R - Figure 4: surface and subsoil values for each paired
## stand (ECM share of ITS reads; observed AM fungal richness), lines joining the two
## depths of a stand, grey where the value was lower in the subsoil and orange otherwise,
## with the mean of each depth in black.
## Descriptive display only; depth contrasts are estimated in 06 (Figure 5).
##
## Output: outputs/figures/paired_depth_stands.png

source("scripts/helpers.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(tidyr) })

env <- load_env() |> mutate(stand_id = paste0(node, stand))
paired_ids <- env |> count(stand_id) |> filter(n == 2) |> pull(stand_id)
d <- env |> filter(stand_id %in% paired_ids) |>
  mutate(depth_lab = factor(ifelse(depth == "A", "0–20 cm", "20–40 cm"), levels = c("0–20 cm", "20–40 cm")),
         ecm_pct = 100 * prop_ecm)
stopifnot(n_distinct(d$stand_id) == 25)
## each stand's line is grey if the value was lower in the subsoil and orange otherwise;
## the black line joins the means of the two depths
panel <- function(y, ylab, tag) {
  dir <- d |> select(stand_id, depth, v = all_of(y)) |> pivot_wider(names_from = depth, values_from = v) |>
    transmute(stand_id, change = ifelse(B < A, "lower", "not lower"))
  n_low <- sum(dir$change == "lower")
  lev <- c(lower = sprintf("Lower in subsoil (%d stands)", n_low),
           `not lower` = sprintf("Not lower (%d stands)", nrow(dir) - n_low))
  dd <- d |> left_join(dir, by = "stand_id") |> mutate(change = factor(lev[change], levels = lev))
  mn <- dd |> group_by(depth_lab) |> summarise(m = mean(.data[[y]]), .groups = "drop")
  ggplot(dd, aes(depth_lab, .data[[y]], group = stand_id, colour = change)) +
    geom_line(linewidth = 0.45, alpha = 0.8) + geom_point(size = 1.2, alpha = 0.8) +
    geom_line(data = mn, aes(depth_lab, m, group = 1), inherit.aes = FALSE, linewidth = 1.3, colour = "black") +
    geom_point(data = mn, aes(depth_lab, m), inherit.aes = FALSE, size = 2.8, colour = "black") +
    scale_colour_manual(values = setNames(c("grey65", "#D55E00"), lev), name = NULL) +
    labs(x = NULL, y = ylab, tag = tag) + theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.direction = "vertical")
}
p <- panel("ecm_pct", "ECM share of ITS reads (%)", "A") | panel("vt_richness", "AM fungal richness (VT per sample)", "B")
ggsave(file.path(FIGS_DIR, "paired_depth_stands.png"), p, width = 6.5, height = 4.2, dpi = 300)
cat("Saved: paired_depth_stands.png\n20_paired_stand_figure.R complete.\n")

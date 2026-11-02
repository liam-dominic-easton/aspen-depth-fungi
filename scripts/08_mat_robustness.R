## 08_mat_robustness.R — within-stand surface-vs-subsoil Aitchison
## distance against mean annual temperature (MAT).
##
## Descriptive figure. Inference on depth x MAT is the paired-shift test with
## MAT permuted among nodes (02_permanova_varpart.R,
## depth_x_env_paired_shift.csv); P values are reported in Table S3,
## not on the panels, because the panels show only the size of the shift.
## The line is an OLS fit shown without a confidence band, because the 25
## stands in 9 nodes are not independent.

source("scripts/helpers.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

ps_its <- load_ps("ITS"); ps_18s <- load_ps("18S")
env <- load_env() |> mutate(stand_id = paste0(node, stand))
vt_mat <- build_vt_matrix(ps_18s)

libs <- list(list(label = "ITS", d = aitchison_distance(ps_its), ps = ps_its),
             list(label = "18S (AMF VT)", d = aitchison_distance_matrix(vt_mat), ps = ps_18s))
plots <- lapply(libs, function(lib) {
  dm <- as.matrix(lib$d)
  md <- align_md(match_ps_to_env(lib$ps, env) |> filter(sample_name %in% labels(lib$d)), lib$d)
  paired <- md |> group_by(stand_id) |> filter(n() == 2) |>
    summarise(within_stand_aitchison = dm[sample_name[depth == "A"], sample_name[depth == "B"]],
              MAT = first(mat_10_year), .groups = "drop")
  ggplot(paired, aes(MAT, within_stand_aitchison)) +
    geom_point(size = 2.6, alpha = 0.75) +
    geom_smooth(method = "lm", se = FALSE, color = "#0072B2", formula = y ~ x) +
    labs(title = lib$label,
         x = "Mean annual temperature (°C)",
         y = "Within-stand Aitchison distance\n(surface vs subsoil)") +
    theme_minimal(base_size = 11)
})
ggsave(file.path(FIGS_DIR, "within_stand_mat.png"), wrap_plots(plots, nrow = 1),
       width = 10, height = 4.5, dpi = 300)
cat("Saved:", file.path(FIGS_DIR, "within_stand_mat.png"), "\n")

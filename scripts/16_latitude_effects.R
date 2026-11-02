## 16_latitude_effects.R — Latitude main effects and latitude x depth
## (predictions: (1) latitude, (2) depth, (3) latitude x depth).
##
## Inputs: load_env(); outputs/tables/guild_ambiguity_sample_proportions.csv (10);
##         data/raw_environmental/Master_with pH added.xlsx (sheet "Sheet1")
## Outputs (outputs/tables/):
##   latitude_premise_checks.csv  climate, fine roots, and soil C versus latitude
##                                (Table S5)
##   latitude_bayes_models.csv    y ~ depth + latitude (+ depth:latitude)
##                                + (1 | node) + (1 | stand) (Table S4)
## Composition versus latitude (node-level permutation) is in Table 1 and the
## depth x latitude paired-shift test in depth_x_env_paired_shift.csv, both
## from 02_permanova_varpart.R.

source("scripts/helpers.R")
source("scripts/bayes_helpers.R")
suppressPackageStartupMessages(library(tidyr))

env <- load_env() |> mutate(depth_num = ifelse(depth == "A", 0, 1), stand_id = paste0(node, stand))
env$lat_z <- standardize(env$latitude)
lat_sd <- sd(env$latitude)
half_read <- 0.5 / pmax(env$its_total_reads, 1)

guild <- read_csv(file.path(TABLES_DIR, "guild_ambiguity_sample_proportions.csv"), show_col_types = FALSE) |>
  select(sample_id_clean, method, guild, guild_reads) |>
  pivot_wider(names_from = guild, values_from = guild_reads) |>
  mutate(lsr = log((Saprotroph + 0.5) / (ECM + 0.5))) |>
  select(sample_id_clean, method, lsr) |>
  pivot_wider(names_from = method, values_from = lsr, names_prefix = "lsr_")
env <- env |> left_join(guild, by = "sample_id_clean")
stopifnot(!anyNA(env$lsr_priority))

## ---- 1. premise checks --------------------------------------------------------
cat("\n=== 1. Premise checks ===\n")
raw <- readxl::read_excel(file.path(DATA_DIR, "raw_environmental", "Master_with pH added.xlsx"),
                          sheet = "Sheet1")
names(raw) <- trimws(names(raw))
raw <- raw |>
  mutate(stand_id = paste0(Node, Stand), depth = Depth, lat = Latitude,
         mat = (`10 year TMaxMean` + `10 year TMinMean`) / 2, cmi = `CMI_10 year`,
         fine_roots = Fine_roots) |>
  filter(!paste0(stand_id, depth) %in% EXCLUDE_SAMPLES)

cor_row <- function(df, x, y, label, depth_lab) {
  d <- df |> filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
  ct <- suppressWarnings(cor.test(d[[x]], d[[y]], method = "spearman"))
  nm <- d |> group_by(Node) |> summarise(x = mean(.data[[x]]), y = mean(.data[[y]]), .groups = "drop")
  cn <- cor.test(nm$x, nm$y)
  tibble(variable = label, depth = depth_lab, n = nrow(d), stand_rho = unname(ct$estimate),
         stand_p = ct$p.value, n_nodes = nrow(nm), node_r = unname(cn$estimate), node_p = cn$p.value)
}
stands <- raw |> filter(depth == "A")
premise <- bind_rows(
  cor_row(stands, "lat", "mat", "MAT", "stand-level"),
  cor_row(stands, "lat", "cmi", "CMI", "stand-level"),
  bind_rows(lapply(c("A", "B"), function(dp) {
    d <- raw |> filter(depth == dp)
    bind_rows(lapply(c("fine_roots", "MAOM_Percent_C", "POM_Percent_C"),
                     function(v) cor_row(d, "lat", v, v, dp)))
  })))
prof <- raw |> select(stand_id, Node, lat, depth, fine_roots, MAOM_Percent_C, POM_Percent_C) |>
  pivot_wider(names_from = depth, values_from = c(fine_roots, MAOM_Percent_C, POM_Percent_C)) |>
  mutate(log_root_ratio = log(fine_roots_A / fine_roots_B),
         log_maom_c_ratio = log(MAOM_Percent_C_A / MAOM_Percent_C_B),
         log_pom_c_ratio = log(POM_Percent_C_A / POM_Percent_C_B))
premise <- bind_rows(premise,
                     cor_row(prof, "lat", "log_root_ratio", "log surface:subsoil fine roots", "profile"),
                     cor_row(prof, "lat", "log_maom_c_ratio", "log surface:subsoil MAOM-C", "profile"),
                     cor_row(prof, "lat", "log_pom_c_ratio", "log surface:subsoil POM-C", "profile"))
print(premise, n = Inf, width = 200)
write_csv(premise, file.path(TABLES_DIR, "latitude_premise_checks.csv"))

## ---- 2. Bayesian latitude models ----------------------------------------------
cat("\n=== 2. Bayesian latitude models ===\n")
outcomes <- tribble(
  ~col,                  ~how,    ~display,
  "lsr_priority",        "none",  "log(sap:ECM), priority",
  "lsr_exclusive_only",  "none",  "log(sap:ECM), exclusive",
  "lsr_fractional",      "none",  "log(sap:ECM), fractional",
  "prop_ecm",            "logit", "ECM proportion",
  "prop_saprotroph",     "logit", "Saprotroph proportion",
  "vt_richness",         "log05", "AMF richness"
)
rows <- list()
for (i in seq_len(nrow(outcomes))) {
  oc <- outcomes[i, ]
  eps <- if (oc$how == "logit") half_read else 0.001
  y <- transform_y(env[[oc$col]], oc$how, eps = eps)
  s <- sd(y)
  d <- data.frame(yz = standardize(y), depth_num = env$depth_num, lat_z = env$lat_z,
                  node = env$node, stand = env$stand_id)
  for (mod in c("main", "interaction")) {
    cat("Fitting", oc$display, "|", mod, "\n")
    fit <- if (mod == "main") bayes_fit(yz ~ depth_num + lat_z + (1 | node) + (1 | stand), d)
           else bayes_fit(yz ~ depth_num * lat_z + (1 | node) + (1 | stand), d)
    dr <- as_draws_df(fit); dg <- fit_diagnostics(fit)
    for (term in intersect(c("b_depth_num", "b_lat_z", "b_depth_num:lat_z"), names(dr))) {
      per <- if (grepl("lat_z", term)) lat_sd else 1
      rows[[length(rows) + 1]] <- tibble(
        outcome = oc$display, transform = oc$how, model = mod, term = sub("^b_", "", term), n = nrow(d),
        !!!as.list(summ(dr[[term]], "std")),
        !!!as.list(summ(dr[[term]] * s / per, "per_unit")),
        !!!as.list(dg))
    }
  }
}
bayes_df <- bind_rows(rows) |> mutate(across(where(is.numeric), ~ signif(.x, 4)))
print(bayes_df |> select(outcome, model, term, std_mean, std_hdi_lo, std_hdi_hi, std_pr_direction, max_rhat, divergences),
      n = Inf, width = 200)
write_csv(bayes_df, file.path(TABLES_DIR, "latitude_bayes_models.csv"))
cat("\n16_latitude_effects.R complete. per_unit_* is per degree of latitude on the transformed scale.\n")

## 17_latitude_som_test.R — does mineral-soil organic matter account for the
## latitude effect on guild balance? (exploratory; Table S1)
##
## Two explanations predict more ECM relative to saprotrophs toward the north of
## this transect, where mineral-soil organic matter is lower: (i) saprotrophs track
## the supply of organic matter; (ii) cold, N-poor conditions favour ECM fungi
## independently of organic matter. Under (i) the latitude coefficient should
## shrink once MAOM-C (and POM-C) enter the model; under (ii) it should not.
## MAOM-C concentration is correlated with sand, so a sand-adjusted model is
## included as a texture check.
##
## Output: outputs/tables/latitude_som_models.csv

source("scripts/helpers.R")
source("scripts/bayes_helpers.R")
suppressPackageStartupMessages(library(tidyr))

env <- load_env() |>
  mutate(depth_num = ifelse(depth == "A", 0, 1), stand_id = paste0(node, stand))
env$lat_z <- standardize(env$latitude)
env$maom_z <- standardize(log(env$maom_percent_c))
env$pom_z <- standardize(log(env$pom_percent_c))
env$sand_z <- standardize(env$sand)
half_read <- 0.5 / pmax(env$its_total_reads, 1)

guild <- read_csv(file.path(TABLES_DIR, "guild_ambiguity_sample_proportions.csv"), show_col_types = FALSE) |>
  select(sample_id_clean, method, guild, guild_reads) |>
  pivot_wider(names_from = guild, values_from = guild_reads) |>
  mutate(lsr = log((Saprotroph + 0.5) / (ECM + 0.5))) |>
  select(sample_id_clean, method, lsr) |>
  pivot_wider(names_from = method, values_from = lsr, names_prefix = "lsr_")
env <- env |> left_join(guild, by = "sample_id_clean")
stopifnot(!anyNA(env$lsr_priority), !anyNA(env$maom_z), !anyNA(env$pom_z), !anyNA(env$sand_z))

outcomes <- tribble(
  ~col,                  ~how,    ~display,
  "lsr_priority",        "none",  "log(sap:ECM), priority",
  "lsr_exclusive_only",  "none",  "log(sap:ECM), exclusive",
  "lsr_fractional",      "none",  "log(sap:ECM), fractional",
  "prop_ecm",            "logit", "ECM proportion",
  "prop_saprotroph",     "logit", "Saprotroph proportion"
)
models <- list(
  latitude          = yz ~ depth_num + lat_z + (1 | node) + (1 | stand),
  latitude_maom     = yz ~ depth_num + lat_z + maom_z + (1 | node) + (1 | stand),
  latitude_maom_pom = yz ~ depth_num + lat_z + maom_z + pom_z + (1 | node) + (1 | stand),
  latitude_maom_sand = yz ~ depth_num + lat_z + maom_z + sand_z + (1 | node) + (1 | stand)
)

rows <- list()
for (i in seq_len(nrow(outcomes))) {
  oc <- outcomes[i, ]
  eps <- if (oc$how == "logit") half_read else 0.001
  y <- transform_y(env[[oc$col]], oc$how, eps = eps)
  d <- data.frame(yz = standardize(y), depth_num = env$depth_num, lat_z = env$lat_z,
                  maom_z = env$maom_z, pom_z = env$pom_z, sand_z = env$sand_z,
                  node = env$node, stand = env$stand_id)
  for (m in names(models)) {
    cat("Fitting", oc$display, "|", m, "\n")
    fit <- bayes_fit(models[[m]], d)
    dr <- as_draws_df(fit); dg <- fit_diagnostics(fit)
    for (term in intersect(c("b_lat_z", "b_maom_z", "b_pom_z", "b_sand_z"), names(dr))) {
      rows[[length(rows) + 1]] <- tibble(
        outcome = oc$display, model = m, term = sub("^b_", "", term), n = nrow(d),
        !!!as.list(summ(dr[[term]], "std")), !!!as.list(dg))
    }
  }
}
out <- bind_rows(rows) |> mutate(across(where(is.numeric), ~ signif(.x, 4)))
print(out |> select(outcome, model, term, std_mean, std_hdi_lo, std_hdi_hi, std_pr_direction, max_rhat),
      n = Inf, width = 200)
write_csv(out, file.path(TABLES_DIR, "latitude_som_models.csv"))

## correlations among the predictors, for interpretation
cat("\nPredictor correlations (samples):\n")
print(round(cor(env[, c("lat_z", "maom_z", "pom_z", "sand_z", "depth_num")]), 2))
cat("\n17_latitude_som_test.R complete.\n")

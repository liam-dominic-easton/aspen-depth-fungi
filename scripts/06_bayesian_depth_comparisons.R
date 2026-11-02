## 06_bayesian_depth_comparisons.R — Bayesian surface-vs-subsoil comparisons
## (Figure 5).
##
## Output: outputs/tables/bayesian_depth_comparisons.csv (read by 09_figures.R)
##
## Model:  y_z ~ depth + (1 | stand), y_z = standardized transformed outcome;
##         priors and sampler settings in bayes_helpers.R.
## Transforms:
##   guild proportions      logit, clamped at half a read (0.5 / ITS reads)
##   log(saprotroph:ECM)    log((sap reads + 0.5) / (ECM reads + 0.5)), priority rule
##   richness               log(x + 0.5)
##   fine-root biomass      log(x + 0.01), as in the Methods and Figure 6
##   ECM tips per m root    untransformed (em_m)
## Reported per outcome:
##   effect_*    subsoil - surface on the transformed scale
##   contrast_*  surface - subsoil on the response scale, computed per posterior
##               draw at the population level (proportions in percentage points)

source("scripts/helpers.R")
source("scripts/bayes_helpers.R")
suppressPackageStartupMessages(library(tidyr))

dat <- load_env() |> mutate(depth_num = ifelse(depth == "A", 0, 1))
guild <- read_csv(file.path(TABLES_DIR, "guild_ambiguity_sample_proportions.csv"),
                  show_col_types = FALSE) |>
  filter(method == "priority") |>
  select(sample_id_clean, guild, guild_reads, total_its_reads) |>
  pivot_wider(names_from = guild, values_from = guild_reads) |>
  mutate(log_sap_ecm = log((Saprotroph + 0.5) / (ECM + 0.5)))
dat <- dat |> left_join(guild |> select(sample_id_clean, log_sap_ecm, total_its_reads),
                        by = "sample_id_clean")
stopifnot(!anyNA(dat$log_sap_ecm))
half_read <- 0.5 / pmax(dat$total_its_reads, 1)

outcomes <- tribble(
  ~col,                      ~how,     ~display,                    ~scale,
  "prop_ecm",                "logit",  "ECM proportion",            100,
  "prop_saprotroph",         "logit",  "Saprotroph proportion",     100,
  "prop_plant_pathogen",     "logit",  "Plant pathogen proportion", 100,
  "log_sap_ecm",             "none",   "Saprotroph:ECM ratio (log)", 1,
  "richness_ECM",            "log05",  "ECM richness",              1,
  "richness_Saprotroph",     "log05",  "Saprotroph richness",       1,
  "richness_Plant pathogen", "log05",  "Pathogen richness",         1,
  "vt_richness",             "log05",  "AMF richness",              1,
  "fine_roots",              "log001", "Fine root biomass",         1,
  "em_m",                    "none",   "ECM tips per m root",       1
)

res <- bind_rows(lapply(seq_len(nrow(outcomes)), function(i) {
  oc <- outcomes[i, ]
  cat("Fitting:", oc$display, "\n")
  eps <- if (oc$how == "logit") half_read else 0.001
  r <- paired_depth_model(dat[[oc$col]], dat$depth_num, dat$node_stand, oc$how,
                          scale_factor = oc$scale, eps = eps)
  y <- dat[[oc$col]]
  tibble(outcome = oc$display, transform = oc$how,
         surface_mean = mean(y[dat$depth == "A"], na.rm = TRUE) * oc$scale,
         subsoil_mean = mean(y[dat$depth == "B"], na.rm = TRUE) * oc$scale,
         !!!as.list(r))
}))
res <- res |> mutate(across(where(is.numeric), ~ signif(.x, 5)))
print(res, n = Inf, width = 250)
write_csv(res, file.path(TABLES_DIR, "bayesian_depth_comparisons.csv"))
cat("\nSaved: outputs/tables/bayesian_depth_comparisons.csv\n")

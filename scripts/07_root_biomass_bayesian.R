## 07_root_biomass_bayesian.R — root biomass and the depth contrast (Figure 6;
## Table S8).
##
## Output: outputs/tables/root_biomass_bayesian_models.csv (read by 09_figures.R)
##         outputs/tables/root_biomass_loo.csv
##
## Three models per outcome (standardized outcome; bayes_helpers.R settings):
##   M1 depth only       y ~ depth + (1 | stand)
##   M2 roots only       y ~ z_log_roots + (1 | stand)
##   M3 depth + roots    y ~ depth + z_log_roots + (1 | stand)
## z_log_roots = standardized log(fine-root biomass + 0.01).
##
## Interpretation: fine-root biomass does not overlap between
## depths (every surface sample exceeds every subsoil sample; r = -0.90 with
## depth), so the depth coefficient in M3 is an extrapolation and is NOT
## interpreted. The estimable quantities are the M1 depth effect and the M3
## root slope, i.e. the association of roots with the outcome among samples
## from the same depth. M3 depth is kept in the table, flagged, for transparency.

source("scripts/helpers.R")
source("scripts/bayes_helpers.R")

dat <- load_env() |>
  mutate(depth_num = ifelse(depth == "A", 0, 1),
         z_log_roots = standardize(log(fine_roots + 0.01)))
half_read <- 0.5 / pmax(dat$its_total_reads, 1)

cat("Root overlap between depths: min surface =", min(dat$fine_roots[dat$depth == "A"]),
    "; max subsoil =", max(dat$fine_roots[dat$depth == "B"]), "\n")
cat("cor(depth, z_log_roots) =", round(cor(dat$depth_num, dat$z_log_roots), 3), "\n")

outcomes <- tribble(
  ~col,              ~how,    ~display,
  "prop_ecm",        "logit", "ECM proportion",
  "prop_saprotroph", "logit", "Saprotroph proportion",
  "vt_richness",     "log05", "AMF richness"
)

rows <- list(); loos <- list()
for (i in seq_len(nrow(outcomes))) {
  oc <- outcomes[i, ]
  cat("\n=== ", oc$display, " ===\n", sep = "")
  eps <- if (oc$how == "logit") half_read else 0.001
  y <- transform_y(dat[[oc$col]], oc$how, eps = eps)
  d <- data.frame(yz = standardize(y), depth_num = dat$depth_num,
                  z_log_roots = dat$z_log_roots, stand = dat$node_stand)
  m1 <- bayes_fit(yz ~ depth_num + (1 | stand), d)
  m2 <- bayes_fit(yz ~ z_log_roots + (1 | stand), d)
  m3 <- bayes_fit(yz ~ depth_num + z_log_roots + (1 | stand), d)
  p1 <- as_draws_df(m1); p2 <- as_draws_df(m2); p3 <- as_draws_df(m3)
  rows[[i]] <- tibble(
    outcome = oc$display,
    !!!as.list(summ(p1$b_depth_num, "m1_depth")),
    !!!as.list(summ(p2$b_z_log_roots, "m2_roots")),
    !!!as.list(summ(p3$b_z_log_roots, "m3_roots_within_depth")),
    !!!as.list(summ(p3$b_depth_num, "m3_depth_not_identifiable")),
    max_rhat = max(fit_diagnostics(m1)[["max_rhat"]], fit_diagnostics(m2)[["max_rhat"]],
                   fit_diagnostics(m3)[["max_rhat"]]),
    divergences = fit_diagnostics(m1)[["divergences"]] + fit_diagnostics(m2)[["divergences"]] +
      fit_diagnostics(m3)[["divergences"]])
  ## PSIS-LOO with exact refits (reloo) for observations with Pareto k > 0.7.
  ## Model names are read from the "model" column when present, because
  ## loo >= 2.8 no longer sets rownames().
  loo_exact <- function(m) {
    l <- loo(m)
    k <- loo::pareto_k_values(l)
    n_bad <- sum(k > 0.7)
    if (n_bad > 0) l <- loo(m, reloo = TRUE)
    list(loo = l, n_k_gt_07 = n_bad, max_k_psis = max(k))
  }
  L <- list(depth_only = loo_exact(m1), roots_only = loo_exact(m2), depth_plus_roots = loo_exact(m3))
  cmp <- loo_compare(lapply(L, `[[`, "loo"))
  model_names <- if ("model" %in% colnames(cmp)) as.character(cmp[, "model"]) else rownames(cmp)
  loos[[i]] <- tibble(outcome = oc$display, model = model_names,
                      elpd_diff = as.numeric(cmp[, "elpd_diff"]), se_diff = as.numeric(cmp[, "se_diff"]),
                      n_obs_k_gt_07_refit = sapply(L, `[[`, "n_k_gt_07")[model_names],
                      max_pareto_k_psis = sapply(L, `[[`, "max_k_psis")[model_names],
                      loo_version = as.character(packageVersion("loo")))
}
res <- bind_rows(rows) |> mutate(across(where(is.numeric), ~ signif(.x, 4)))
loo_df <- bind_rows(loos) |> mutate(across(where(is.numeric), ~ signif(.x, 4)))
print(res, width = 250); print(loo_df)
write_csv(res, file.path(TABLES_DIR, "root_biomass_bayesian_models.csv"))
write_csv(loo_df, file.path(TABLES_DIR, "root_biomass_loo.csv"))
cat("\nSaved: outputs/tables/root_biomass_bayesian_models.csv, root_biomass_loo.csv\n")

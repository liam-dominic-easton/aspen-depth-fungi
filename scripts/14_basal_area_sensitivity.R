## 14_basal_area_sensitivity.R -- targeted stand-structure sensitivity
## (Table S3).
##
## Basal area is constant within a stand, so it cannot explain the paired depth
## contrast as a main effect; the question is whether the surface->subsoil
## contrast changes with basal area.
##   * composition: paired CLR shifts regressed on basal area, with basal area
##     node means permuted among nodes;
##   * univariate: Bayesian regressions of the paired change in
##     log(saprotroph:ECM) (priority rule) and in AMF richness on standardized
##     basal area.
## Basal area is missing for NOT2 and NOT3.

source("scripts/helpers.R")
source("scripts/bayes_helpers.R")
suppressPackageStartupMessages(library(tidyr))

ps_its <- load_ps("ITS"); ps_18s <- load_ps("18S")
env <- load_env() |> mutate(stand_id = paste0(node, stand))
vt_mat <- build_vt_matrix(ps_18s)
its_otu <- as(otu_table(ps_its), "matrix"); if (taxa_are_rows(ps_its)) its_otu <- t(its_otu)

community <- bind_rows(lapply(list(list("ITS", its_otu, ps_its), list("AMF", vt_mat, ps_18s)), function(x) {
  md <- match_ps_to_env(x[[3]], env) |> filter(sample_name %in% rownames(x[[2]]))
  sh <- paired_shifts(clr_matrix(x[[2]]), md)
  ok <- !is.na(sh$st$basal_area)
  sh$D <- sh$D[ok, , drop = FALSE]; sh$st <- sh$st[ok, ]
  r <- paired_shift_test(sh, "basal_area", unit = "node")
  tibble(library = x[[1]], response = "community composition (paired shift)",
         test = "depth x basal area", n_stands = sum(ok),
         R2 = r[["R2"]], pseudo_F = r[["F"]], p_value = r[["p"]], permutation = "among nodes (node means)")
}))

stand_ba <- env |> distinct(stand_id, basal_area)
guild_change <- read_csv(file.path(TABLES_DIR, "guild_ambiguity_sample_proportions.csv"),
                         show_col_types = FALSE) |>
  filter(method == "priority", guild %in% c("Saprotroph", "ECM")) |>
  select(node_stand, depth, guild, guild_reads) |>
  pivot_wider(names_from = guild, values_from = guild_reads) |>
  mutate(log_sap_ecm = log((Saprotroph + 0.5) / (ECM + 0.5))) |>
  select(node_stand, depth, log_sap_ecm) |>
  pivot_wider(names_from = depth, values_from = log_sap_ecm) |>
  transmute(stand_id = node_stand, depth_change = B - A) |>
  left_join(stand_ba, by = "stand_id")
amf_change <- tibble(sample_name = rownames(vt_mat), richness = rowSums(vt_mat > 0)) |>
  left_join(match_ps_to_env(ps_18s, env), by = "sample_name") |>
  select(stand_id, depth, richness) |>
  pivot_wider(names_from = depth, values_from = richness) |>
  transmute(stand_id, depth_change = B - A) |>
  left_join(stand_ba, by = "stand_id")

uni <- bind_rows(lapply(list(list("ITS", "paired change in log(saprotroph:ECM)", guild_change),
                             list("AMF", "paired change in AMF richness", amf_change)), function(x) {
  d <- x[[3]] |> filter(!is.na(depth_change), !is.na(basal_area))
  dd <- data.frame(yz = standardize(d$depth_change), xz = standardize(d$basal_area))
  fit <- bayes_fit(yz ~ xz, dd)
  tibble(library = x[[1]], response = x[[2]], test = "depth change x basal area", n_stands = nrow(d),
         !!!as.list(summ(as_draws_df(fit)$b_xz, "std_slope")), !!!as.list(fit_diagnostics(fit)))
}))

write_csv(community, file.path(TABLES_DIR, "basal_area_community_sensitivity.csv"))
write_csv(uni, file.path(TABLES_DIR, "basal_area_depth_change_sensitivity.csv"))
print(community, width = Inf); print(uni, width = Inf)

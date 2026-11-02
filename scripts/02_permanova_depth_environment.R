## 02_permanova_depth_environment.R — Community composition: PERMANOVA and
## depth x environment tests.
##
## Outputs (outputs/tables/):
##   table1_permanova.csv             Table 1: depth, latitude or biome, pH
##   depth_x_env_paired_shift.csv     depth x latitude / MAT / CMI / sand (Table S3)
##   depth_x_mat_given_biome.csv      depth x MAT with depth x biome held constant (Table S3)
##   depth_x_temperature_metrics.csv  depth x Tmax / Tmin / MAT (Table S3)
## (Four-set variance partitioning, Table 3, is in 12_guild_balance_roots_varpart.R.)
##
## Communities: ITS (all ASVs) and AMF virtual taxa (VT), CLR with 0.5
## pseudocount, Aitchison distance.
##
## Inference follows the nesting of the design (see helpers.R):
##   depth      stand conditioned out (dbRDA + Condition(stand)); residuals
##              permuted within stands
##   latitude,  permuted among the 9 nodes (they vary only among nodes)
##   biome
##   pH         node means permuted among nodes (pH varies mainly among nodes)
##   depth x    paired surface->subsoil CLR shifts regressed on the moderator,
##   moderator  moderator permuted among nodes (node-level) or stands (sand)
## R2 values are marginal (adonis2 by = "margin" equivalent) and descriptive.
## Unrestricted permutations of samples would pseudoreplicate node- and
## stand-level terms, so each term is permuted at the level at which it varies.

source("scripts/helpers.R")
suppressPackageStartupMessages(library(tidyr))

cat("\n=== Loading data ===\n")
ps_its <- load_ps("ITS")
ps_18s <- load_ps("18S")
env    <- load_env() |> mutate(stand_id = paste0(node, stand))
vt_mat <- build_vt_matrix(ps_18s)

its_otu <- as(otu_table(ps_its), "matrix"); if (taxa_are_rows(ps_its)) its_otu <- t(its_otu)

make_lib <- function(label, counts, ps) {
  d  <- aitchison_distance_matrix(counts)
  md <- align_md(match_ps_to_env(ps, env) |> filter(sample_name %in% rownames(counts)), d)
  list(label = label, d = d, md = md, G = gower_centre(d),
       shifts = paired_shifts(clr_matrix(counts), md))
}
libs <- list(make_lib("ITS", its_otu, ps_its),
             make_lib("18S (AMF VT)", vt_mat, ps_18s))

## ---- 1. Table 1 ---------------------------------------------------------------
cat("\n=== Table 1 ===\n")
table1_rows <- function(lib, drop_nodes = character(0)) {
  keep <- !lib$md$node %in% drop_nodes
  md <- lib$md[keep, ]
  d  <- as.dist(as.matrix(lib$d)[keep, keep])
  G  <- gower_centre(d)
  out <- list()
  for (geo in c("latitude", "biome")) {
    terms <- c("depth_factor", geo, "p_h")
    dep <- test_within_stand(d, md, "depth_factor", covars = "p_h")
    r_d <- marginal_F(G, md, terms, "depth_factor")
    r_g <- unit_perm_test(G, md, terms, geo, unit = "node")
    r_p <- unit_perm_test(G, md, terms, "p_h", unit = "node")
    out[[geo]] <- tibble(
      library = lib$label, model = paste(terms, collapse = " + "),
      term = c("depth", geo, "pH"),
      R2 = c(r_d[["R2"]], r_g[["R2"]], r_p[["R2"]]),
      F = c(dep[["F"]], r_g[["F"]], r_p[["F"]]),
      p = c(dep[["p"]], r_g[["p"]], r_p[["p"]]),
      permutation = c("within stands, stand conditioned", "among nodes", "among nodes (node means)"),
      n = nrow(md))
  }
  bind_rows(out)
}
table1 <- bind_rows(lapply(libs, table1_rows))
print(table1, n = Inf, width = 200)
write_csv(table1, file.path(TABLES_DIR, "table1_permanova.csv"))

## ---- 2. depth x environment (paired shifts) -------------------------------------
cat("\n=== Depth x environment: paired surface->subsoil shifts ===\n")
mods <- tribble(
  ~moderator,     ~unit,
  "latitude",     "node",
  "mat_10_year",  "node",
  "cmi_10_year",  "node",
  "sand",         "node"
)
ints <- bind_rows(lapply(libs, function(lib) bind_rows(lapply(seq_len(nrow(mods)), function(i) {
  r <- paired_shift_test(lib$shifts, mods$moderator[i], unit = mods$unit[i])
  tibble(library = lib$label, moderator = mods$moderator[i], permutation = paste0("among ", mods$unit[i], "s"),
         n_stands = nrow(lib$shifts$D), R2 = r[["R2"]], F = r[["F"]], p = r[["p"]])
})))) |> mutate(p_bh = p.adjust(p, "BH"))
print(ints, n = Inf, width = 200)
write_csv(ints, file.path(TABLES_DIR, "depth_x_env_paired_shift.csv"))

cat("\n=== Depth x MAT holding depth x biome constant ===\n")
mat_b <- bind_rows(lapply(libs, function(lib) {
  r <- paired_shift_test(lib$shifts, "mat_10_year", unit = "node", covars = "biome", strata = "biome")
  tibble(library = lib$label, moderator = "mat_10_year | biome",
         permutation = "among nodes within biome", R2 = r[["R2"]], F = r[["F"]], p = r[["p"]])
}))
print(mat_b, width = 200)
write_csv(mat_b, file.path(TABLES_DIR, "depth_x_mat_given_biome.csv"))

cat("\n=== Depth x temperature metric ===\n")
temps <- c(Tmax = "x10_year_t_max_mean", Tmin = "x10_year_t_min_mean", MAT = "mat_10_year")
tm <- bind_rows(lapply(libs, function(lib) bind_rows(lapply(names(temps), function(nm) {
  r <- paired_shift_test(lib$shifts, temps[[nm]], unit = "node")
  tibble(library = lib$label, metric = nm, R2 = r[["R2"]], F = r[["F"]], p = r[["p"]])
})))) |> mutate(p_bh = p.adjust(p, "BH"))
print(tm, width = 200)
write_csv(tm, file.path(TABLES_DIR, "depth_x_temperature_metrics.csv"))

cat("\n=== Done ===\n")

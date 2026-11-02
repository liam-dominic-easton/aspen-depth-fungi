## 12_guild_balance_roots_varpart.R
## Guild balance, roots, soil organic matter and variance partitioning:
##   1) saprotroph:ECM depth contrasts under all guild-ambiguity treatments;
##   2) root-biomass and soil-organic-matter models (Table S8);
##   3) paired root-depth-ratio tests (Table S7);
##   4) four-set variance partitioning with unique-fraction tests
##      (Table 3, Table S9, Figure S1).

source("scripts/helpers.R")
source("scripts/bayes_helpers.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

guild_file <- file.path(TABLES_DIR, "guild_ambiguity_sample_proportions.csv")
if (!file.exists(guild_file)) {
  stop("Run scripts/10_guild_ambiguity_sensitivity.R first")
}

env <- load_env() |>
  mutate(
    depth_num = if_else(depth == "B", 1, 0),
    node_stand = paste0(node, stand),
    log_fine_roots = log(fine_roots + 0.01),
    ## Organic matter enters as the measured C concentrations of the two
    ## fractions (as in Karst et al. 2021). Fraction masses were not recorded,
    ## so the fractions cannot be summed into a bulk-soil C total.
    log_maom_c = log(maom_percent_c),
    log_pom_c = log(pom_percent_c)
  ) |>
  mutate(
    z_log_fine_roots = as.numeric(scale(log_fine_roots)),
    z_log_maom_c = as.numeric(scale(log_maom_c)),
    z_log_pom_c = as.numeric(scale(log_pom_c))
  )

guild <- read_csv(guild_file, show_col_types = FALSE) |>
  select(sample_id_clean, guild, method, proportion, guild_reads, total_its_reads) |>
  pivot_wider(
    names_from = guild,
    values_from = c(proportion, guild_reads),
    names_sep = "_"
  ) |>
  left_join(env, by = "sample_id_clean") |>
  mutate(
    ## Half-read continuity correction, expressed on the original read scale.
    log_sap_ecm = log((guild_reads_Saprotroph + 0.5) /
                        (guild_reads_ECM + 0.5))
  )

## Models use the package's Bayesian settings (bayes_helpers.R). Outcomes are standardized;
## coefficients are back-transformed to the outcome's transformed scale per SD of
## each (already standardized) predictor.
## Caveats written into the output:
##   * depth and fine roots do not overlap between depths, so the depth
##     coefficient in any model that also contains roots is not identifiable;
##   * organic matter = MAOM-C and POM-C concentrations (a summed C stock or
##     combined C:N would require fraction masses, which were not recorded).
bayes_terms <- function(f, d, outcome, method, model) {
  vars <- all.vars(f)
  d <- d[complete.cases(d[, vars]), ]
  y <- d[[vars[1]]]; s <- sd(y)
  d$yz <- (y - mean(y)) / s
  f2 <- update(f, yz ~ .)
  fit <- bayes_fit(f2, d)
  dr <- as_draws_df(fit)
  terms <- grep("^b_", names(dr), value = TRUE)
  terms <- setdiff(terms, "b_Intercept")
  bind_rows(lapply(terms, function(t) {
    tibble(outcome = outcome, method = method, model = model, term = sub("^b_", "", t),
           !!!as.list(summ(dr[[t]] * s, "estimate")), n = nrow(d),
           note = case_when(
             sub("^b_", "", t) == "depth_num" & "z_log_fine_roots" %in% vars ~ "not identifiable: roots and depth do not overlap",
             TRUE ~ ""))
  })) |> bind_cols(as_tibble(as.list(fit_diagnostics(fit))))
}

## ---- saprotroph:ECM models ---------------------------------------------------------
ratio_models <- bind_rows(lapply(unique(guild$method), function(m) {
  d <- filter(guild, method == m) |> mutate(stand = node_stand)
  forms <- list(
    depth_only  = log_sap_ecm ~ depth_num + (1 | stand),
    depth_roots = log_sap_ecm ~ depth_num + z_log_fine_roots + (1 | stand),
    depth_som   = log_sap_ecm ~ depth_num + z_log_maom_c + z_log_pom_c + (1 | stand),
    full        = log_sap_ecm ~ depth_num + z_log_fine_roots + z_log_maom_c + z_log_pom_c + (1 | stand)
  )
  bind_rows(lapply(names(forms), function(nm) bayes_terms(forms[[nm]], d, "log(saprotroph:ECM)", m, nm)))
}))
write_csv(ratio_models, file.path(TABLES_DIR, "sap_ecm_root_som_models.csv"))

## ---- ECM and saprotroph component models ------------------------------------------
component_models <- bind_rows(lapply(unique(guild$method), function(m) {
  d <- filter(guild, method == m) |> mutate(stand = node_stand)
  eps <- 0.5 / pmax(d$total_its_reads, 1)
  d$y_ECM <- qlogis(pmin(pmax(d$proportion_ECM, eps), 1 - eps))
  d$y_Saprotroph <- qlogis(pmin(pmax(d$proportion_Saprotroph, eps), 1 - eps))
  bind_rows(
    bayes_terms(y_ECM ~ depth_num + z_log_fine_roots + z_log_maom_c + z_log_pom_c + (1 | stand),
                d, "ECM", m, "full"),
    bayes_terms(y_Saprotroph ~ depth_num + z_log_fine_roots + z_log_maom_c + z_log_pom_c + (1 | stand),
                d, "Saprotroph", m, "full"))
}))
write_csv(component_models, file.path(TABLES_DIR, "guild_component_root_som_models.csv"))

## ---- Paired root-ratio test -------------------------------------------------
paired <- guild |>
  select(method, node_stand, depth, log_sap_ecm, fine_roots,
         maom_percent_c, pom_percent_c, sand) |>
  pivot_wider(
    names_from = depth,
    values_from = c(log_sap_ecm, fine_roots, maom_percent_c, pom_percent_c, sand),
    names_sep = "_"
  ) |>
  mutate(
    delta_log_sap_ecm = log_sap_ecm_B - log_sap_ecm_A,
    log_root_ratio = log((fine_roots_A + 0.01) / (fine_roots_B + 0.01)),
    log_maom_c_ratio = log(maom_percent_c_A / maom_percent_c_B),
    log_pom_c_ratio = log(pom_percent_c_A / pom_percent_c_B),
    ## stand texture: does the depth change vary with sand?
    sand_stand_mean = (sand_A + sand_B) / 2,
    ## prediction 2: how much more steeply roots decline than
    ## each organic-matter fraction (log of root ratio / C ratio)
    log_root_vs_maom = log_root_ratio - log_maom_c_ratio,
    log_root_vs_pom = log_root_ratio - log_pom_c_ratio
  )

## stand-level Bayesian regressions of the paired change on each predictor
## (both standardized, so the slope is a correlation-scale effect); Spearman rho
## kept as a descriptive statistic only
paired_tests <- bind_rows(lapply(unique(paired$method), function(m) {
  d <- filter(paired, method == m)
  bind_rows(lapply(c("log_root_ratio", "log_maom_c_ratio", "log_pom_c_ratio",
                                              "log_root_vs_maom", "log_root_vs_pom", "sand_stand_mean"), function(x) {
    dd <- d[complete.cases(d[, c("delta_log_sap_ecm", x)]), ]
    dd <- data.frame(yz = standardize(dd$delta_log_sap_ecm), xz = standardize(dd[[x]]))
    fit <- bayes_fit(yz ~ xz, dd)
    tibble(method = m, predictor = x, n_stands = nrow(dd),
           !!!as.list(summ(as_draws_df(fit)$b_xz, "std_slope")),
           spearman_rho = cor(dd$yz, dd$xz, method = "spearman"),
           !!!as.list(fit_diagnostics(fit)))
  }))
}))
write_csv(paired, file.path(TABLES_DIR, "sap_ecm_paired_depth_changes.csv"))
write_csv(paired_tests, file.path(TABLES_DIR, "sap_ecm_paired_predictor_tests.csv"))

## ---- Four-set variance partitioning ----------------------------------------
ps_its <- load_ps("ITS")
ps_18s <- load_ps("18S")

build_lib <- function(ps, label, vt_matrix = NULL) {
  if (is.null(vt_matrix)) {
    d <- aitchison_distance(ps)
    md <- match_ps_to_env(ps, env)
  } else {
    d <- aitchison_distance_matrix(vt_matrix)
    md <- match_ps_to_env(ps, env) |>
      filter(sample_name %in% rownames(vt_matrix))
  }
  md <- md |>
    filter(sample_name %in% labels(d)) |>
    arrange(match(sample_name, labels(d))) |>
    mutate(stand_id = paste0(node, stand))
  list(label = label, dist = d, md = md)
}

libs <- list(
  build_lib(ps_its, "ITS"),
  build_lib(ps_18s, "18S (AMF VT)", build_vt_matrix(ps_18s))
)

partial_test <- function(d, md, focal, conditioning, restricted) {
  rhs <- paste(focal, "+ Condition(", paste(conditioning, collapse = " + "), ")")
  mod <- capscale(as.formula(paste("d ~", rhs)), data = md, add = TRUE)
  set.seed(42)
  perms <- if (restricted) {
    permute::how(blocks = factor(md$stand_id), nperm = 999)
  } else {
    999
  }
  a <- anova.cca(mod, permutations = perms)
  tibble(
    df = a$Df[1],
    variance = a$Variance[1],
    pseudo_F = a$F[1],
    p_value = a$`Pr(>F)`[1],
    permutation_scheme = if_else(restricted, "within stand", "unrestricted")
  )
}

vp_fractions <- list()
vp_tests <- list()
vp_objects <- list()
for (lib in libs) {
  preds <- c("depth_num", "sand", "p_h", "cmi_10_year", "mat_10_year",
             "log_fine_roots")
  cc <- complete.cases(lib$md[, preds])
  md <- lib$md[cc, ]
  d <- as.dist(as.matrix(lib$dist)[cc, cc])

  X_depth <- data.frame(depth = md$depth_num)
  X_soil <- data.frame(sand = scale(md$sand)[, 1], pH = scale(md$p_h)[, 1])
  X_climate <- data.frame(CMI = scale(md$cmi_10_year)[, 1],
                          MAT = scale(md$mat_10_year)[, 1])
  X_roots <- data.frame(log_roots = scale(md$log_fine_roots)[, 1])

  vp <- varpart(d, X_depth, X_soil, X_climate, X_roots)
  vp_objects[[lib$label]] <- vp
  ind <- vp$part$indfract
  ## vegan varpart4 layout (X1 = Depth, X2 = Soil, X3 = Climate, X4 = Roots):
  ## [e] X1X2, [f] X2X3, [g] X1X3, [h] X1X4, [i] X2X4, [j] X3X4,
  ## [k] X1X2X4, [l] X1X2X3, [m] X2X3X4, [n] X1X3X4, [o] all.
  x <- "\u2229"
  fraction_labels <- c(
    "Depth unique", "Soil unique", "Climate unique", "Roots unique",
    paste("Depth", x, "Soil"), paste("Soil", x, "Climate"), paste("Depth", x, "Climate"),
    paste("Depth", x, "Roots"), paste("Soil", x, "Roots"), paste("Climate", x, "Roots"),
    paste("Depth", x, "Soil", x, "Roots"), paste("Depth", x, "Soil", x, "Climate"),
    paste("Soil", x, "Climate", x, "Roots"), paste("Depth", x, "Climate", x, "Roots"),
    "All four sets", "Residual"
  )
  vp_fractions[[lib$label]] <- tibble(
    library = lib$label,
    fraction_id = rownames(ind),
    fraction = fraction_labels,
    adj_r2 = ind$Adj.R.square
  )

  ## Unique-fraction tests at the level each set varies:
  ##   Depth, Roots  stand conditioned out, residuals permuted within stands
  ##                 (climate is constant within stands and so is absorbed)
  ##   Soil          node means of sand + pH permuted jointly among nodes
  ##   Climate       CMI + MAT permuted jointly among nodes
  G <- gower_centre(d)
  all_terms <- c("depth_num", "sand", "p_h", "cmi_10_year", "mat_10_year", "log_fine_roots")
  wd <- test_within_stand(d, md, "depth_num", covars = c("sand", "p_h", "log_fine_roots"))
  wr <- test_within_stand(d, md, "log_fine_roots", covars = c("depth_num", "sand", "p_h"))
  so <- unit_perm_test(G, md, all_terms, c("sand", "p_h"), unit = "node")
  cl <- unit_perm_test(G, md, all_terms, c("cmi_10_year", "mat_10_year"), unit = "node")
  vp_tests[[lib$label]] <- tibble(
    library = lib$label,
    predictor_set = c("Depth", "Roots", "Soil", "Climate"),
    pseudo_F = c(wd[["F"]], wr[["F"]], so[["F"]], cl[["F"]]),
    p_value = c(wd[["p"]], wr[["p"]], so[["p"]], cl[["p"]]),
    permutation_scheme = c("within stands, stand conditioned", "within stands, stand conditioned",
                           "among nodes (node means)", "among nodes"),
    nperm = NPERM)
}

## Native four-set vegan diagram, combined across libraries for the manuscript.
## The geometry is schematic rather than area-proportional; the labels are the
## adjusted-R² fractions returned by varpart, rounded to three decimals.
format_vp_fraction <- function(x) {
  rounded <- round(x, 3)
  rounded[rounded == 0] <- 0 # avoid displaying negative zero
  sprintf("%.3f", rounded)
}

venn_path <- file.path(FIGS_DIR, "variance_partition_venn.png")
png(venn_path, width = 10.5, height = 5.4, units = "in", res = 400)
par(mfrow = c(1, 2), mar = c(1.2, 1.2, 2.8, 1.2),
    oma = c(0, 0, 1.7, 0), xpd = NA)

vp_cols <- c(Depth = "#E69F00", Soil = "#56B4E9",
             Climate = "#009E73", Roots = "#CC79A7")
panel_titles <- c("ITS" = "A  ITS fungi", "18S (AMF VT)" = "B  AMF")

for (lib_name in names(vp_objects)) {
  vp <- vp_objects[[lib_name]]
  vegan::showvarparts(
    parts = 4,
    labels = format_vp_fraction(vp$part$indfract$Adj.R.square),
    bg = unname(vp_cols), alpha = 78,
    Xnames = names(vp_cols), id.size = 0.92,
    cex = 0.73
  )
  title(main = panel_titles[[lib_name]], line = 0.35,
        cex.main = 1.15, font.main = 2)
}
mtext("Four-set variance partitioning (adjusted R²)", side = 3,
      outer = TRUE, line = 0.35, cex = 1.05, font = 2)
dev.off()
cat("Saved: outputs/figures/variance_partition_venn.png\n")

vp_fraction_table <- bind_rows(vp_fractions)
write_csv(vp_fraction_table, file.path(TABLES_DIR, "varpart_four_set_all_fractions.csv"))
vp_fraction_table |>
  filter(fraction %in% c("Depth unique", "Soil unique", "Climate unique",
                         "Roots unique", "Residual")) |>
  select(library, fraction, adj_r2) |>
  write_csv(file.path(TABLES_DIR, "table2_varpart.csv"))
write_csv(bind_rows(vp_tests), file.path(TABLES_DIR, "varpart_four_set_unique_tests.csv"))

cat("\nsaprotroph:ECM depth coefficients\n")
print(filter(ratio_models, term == "depth_num"), n = Inf, width = Inf)
cat("\nPaired predictor tests\n")
print(paired_tests, n = Inf, width = Inf)
cat("\nUnique variance-fraction tests\n")
print(bind_rows(vp_tests), n = Inf, width = Inf)

## 13_additional_sensitivities.R
## Additional sensitivity analyses:
##   A) guild:ECM log-ratios, so the pathogen depth signal can be separated from
##      the compositional consequence of the ECM decline (Table S12);
##   B) FUNGuild confidence-rank sensitivity (probable + highly probable only;
##      Table S11);
##   C) correlations between root ECM colonization and bulk-soil DNA measures
##      (Table S2);
##   D) understory mycorrhizal composition versus latitude (understory rows of
##      Table S5).

source("scripts/helpers.R")
source("scripts/bayes_helpers.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

## ---- shared FUNGuild flag construction (mirrors script 10) ------------------

fg <- read_tsv(
  file.path(TAX_DIR, "otu_table_for_FUNGuild.guilds.txt"),
  show_col_types = FALSE,
  name_repair = "minimal",
  guess_max = 10000
)

meta_cols <- c(
  "OTU_ID", "taxonomy", "Taxon", "Taxon Level", "Trophic Mode", "Guild",
  "Growth Morphology", "Trait", "Confidence Ranking", "Notes", "Citation/Source"
)

ps_its <- load_ps("ITS")
key <- get_sample_key(ps_its) |>
  transmute(
    sample_name, sample_id_clean, node, stand, depth,
    node_stand = paste0(node, stand),
    depth_num = if_else(depth == "B", 1, 0)
  )
sample_cols <- intersect(setdiff(names(fg), meta_cols), key$sample_name)
stopifnot(length(sample_cols) == nrow(key))

fg_flags <- fg |>
  mutate(
    Guild = str_squish(replace_na(Guild, "")),
    guild_lower = str_to_lower(Guild),
    confidence = str_squish(replace_na(`Confidence Ranking`, "")),
    is_unassigned = Guild %in% c("", "-") |
      str_detect(guild_lower, "^\\s*(unassigned|unknown|undefined|none)\\s*$"),
    flag_ecm = str_detect(guild_lower, "ectomycorrhizal"),
    flag_amf = str_detect(guild_lower, "arbuscular"),
    flag_saprotroph = str_detect(guild_lower, "saprotroph"),
    flag_pathogen = str_detect(guild_lower, "pathogen"),
    flag_endophyte = str_detect(guild_lower, "endophyte"),
    n_focal_flags = flag_ecm + flag_amf + flag_saprotroph + flag_pathogen + flag_endophyte,
    priority_category = case_when(
      is_unassigned ~ "Unassigned",
      flag_ecm ~ "ECM",
      flag_amf ~ "AMF",
      flag_saprotroph ~ "Saprotroph",
      flag_pathogen ~ "Pathogen",
      flag_endophyte ~ "Endophyte",
      nzchar(Guild) ~ "Other",
      TRUE ~ "Unassigned"
    )
  )

focal <- c("ECM", "Saprotroph", "Pathogen")
flag_name <- c(ECM = "flag_ecm", Saprotroph = "flag_saprotroph",
               Pathogen = "flag_pathogen")

## Read counts per guild under each ambiguity treatment, for an arbitrary subset
## of the FUNGuild rows. The ITS denominator always uses all retained reads, so
## proportions stay comparable when rows are dropped.
guild_reads_for <- function(flags) {
  counts_long <- flags |>
    select(OTU_ID, priority_category, n_focal_flags,
           flag_ecm, flag_amf, flag_saprotroph, flag_pathogen, flag_endophyte,
           all_of(sample_cols)) |>
    pivot_longer(all_of(sample_cols), names_to = "sample_name", values_to = "reads") |>
    filter(reads > 0) |>
    left_join(key, by = "sample_name")

  priority <- counts_long |>
    filter(priority_category %in% focal) |>
    group_by(sample_name, guild = priority_category) |>
    summarise(guild_reads = sum(reads), .groups = "drop") |>
    mutate(method = "priority")

  exclusive <- bind_rows(lapply(focal, function(g) {
    counts_long |>
      filter(n_focal_flags == 1, .data[[flag_name[[g]]]]) |>
      group_by(sample_name) |>
      summarise(guild_reads = sum(reads), .groups = "drop") |>
      mutate(guild = g, method = "exclusive_only")
  }))

  fractional <- bind_rows(lapply(focal, function(g) {
    counts_long |>
      filter(n_focal_flags >= 1, .data[[flag_name[[g]]]]) |>
      group_by(sample_name) |>
      summarise(guild_reads = sum(reads / n_focal_flags), .groups = "drop") |>
      mutate(guild = g, method = "fractional")
  }))

  bind_rows(priority, exclusive, fractional) |>
    complete(sample_name = key$sample_name, guild = focal,
             method = c("priority", "exclusive_only", "fractional"),
             fill = list(guild_reads = 0)) |>
    left_join(key, by = "sample_name")
}

## The denominator is all retained ITS reads, taken once from the full table.
all_counts <- as.matrix(fg_flags[, sample_cols])
storage.mode(all_counts) <- "numeric"
sample_totals <- tibble(sample_name = sample_cols,
                        total_its_reads = colSums(all_counts))

## ---- A. guild:ECM log-ratios ------------------------------------------------

ratio_input <- guild_reads_for(fg_flags) |>
  left_join(sample_totals, by = "sample_name") |>
  select(sample_name, sample_id_clean, node_stand, depth, depth_num,
         method, guild, guild_reads, total_its_reads) |>
  pivot_wider(names_from = guild, values_from = guild_reads) |>
  mutate(
    log_sap_ecm = log((Saprotroph + 0.5) / (ECM + 0.5)),
    log_path_ecm = log((Pathogen + 0.5) / (ECM + 0.5)),
    ## Pathogens against the other free-living guild. Ratios taken against ECM
    ## inherit the ECM depth decline; this one does not, so it tests whether the
    ## pathogen signal is separable from that decline.
    log_path_sap = log((Pathogen + 0.5) / (Saprotroph + 0.5))
  )

ratio_labels <- c(
  log_sap_ecm = "log(saprotroph:ECM)",
  log_path_ecm = "log(pathogen:ECM)",
  log_path_sap = "log(pathogen:saprotroph)"
)

ratio_depth_models <- bind_rows(lapply(
  names(ratio_labels), function(ycol) {
    bind_rows(lapply(unique(ratio_input$method), function(m) {
      d <- filter(ratio_input, method == m)
      d$y <- d[[ycol]]
      paired_summary_bayes(d, "y") |>
        mutate(ratio = unname(ratio_labels[[ycol]]), method = m, .before = 1)
    }))
  }))

ratio_depth_models <- ratio_depth_models |>
  arrange(ratio, factor(method, levels = c("priority", "exclusive_only", "fractional")))
write_csv(ratio_depth_models, file.path(TABLES_DIR, "guild_ratio_depth_models.csv"))

## ---- B. FUNGuild confidence-rank sensitivity --------------------------------

keep_conf <- c("Probable", "Highly Probable")
fg_conf <- fg_flags |>
  mutate(across(all_of(sample_cols), ~ if_else(confidence %in% keep_conf, .x, 0)))

cat("ASVs retained at probable/highly probable confidence:",
    sum(fg_flags$confidence %in% keep_conf), "of", nrow(fg_flags), "\n")

conf_props <- guild_reads_for(fg_conf) |>
  left_join(sample_totals, by = "sample_name") |>
  mutate(proportion = guild_reads / total_its_reads)

analyse_one <- function(df) {
  df <- arrange(df, node_stand, depth_num)
  eps <- 0.5 / pmax(df$total_its_reads, 1)
  df$logit_prop <- qlogis(pmin(pmax(df$proportion, eps), 1 - eps))
  paired_summary_bayes(df, "logit_prop", "proportion", scale = 100)
}

conf_inference <- conf_props |>
  group_by(guild, method) |>
  group_modify(~ analyse_one(.x)) |>
  ungroup() |>
  arrange(guild, factor(method, levels = c("priority", "exclusive_only", "fractional")))
write_csv(conf_inference, file.path(TABLES_DIR, "guild_confidence_sensitivity.csv"))

## log(saprotroph:ECM) under the confidence restriction
conf_ratio <- conf_props |>
  select(sample_name, node_stand, depth, depth_num, method, guild, guild_reads) |>
  pivot_wider(names_from = guild, values_from = guild_reads) |>
  mutate(y = log((Saprotroph + 0.5) / (ECM + 0.5)))
conf_ratio_models <- bind_rows(lapply(unique(conf_ratio$method), function(m) {
  paired_summary_bayes(filter(conf_ratio, method == m), "y") |>
    mutate(ratio = "log(saprotroph:ECM)", method = m, confidence = "probable/highly probable", .before = 1)
}))
write_csv(conf_ratio_models, file.path(TABLES_DIR, "guild_confidence_ratio_models.csv"))
print(conf_ratio_models, width = Inf)

## ---- C. root ECM colonization versus bulk-soil DNA --------------------------

env <- load_env()
dna_measures <- c(
  prop_ecm = "ECM proportion",
  prop_saprotroph = "Saprotroph proportion",
  richness_ECM = "ECM richness",
  vt_richness = "AMF VT richness"
)

## Two colonization metrics exist: em_colonization = ECM
## tips per 100 live tips scored (%); em_m = ECM tips per metre of fine root.
## Both are reported, with BH adjustment across the 8 tests within each metric.
colonization_tests <- bind_rows(lapply(c("em_colonization", "em_m"), function(metric) {
  bind_rows(lapply(names(dna_measures), function(v) {
    bind_rows(lapply(c("A", "B"), function(dp) {
      d <- env |> filter(depth == dp) |> select(all_of(c(metric, v))) |>
        filter(complete.cases(across(everything())))
      ct <- cor.test(d[[metric]], d[[v]], method = "pearson")
      tibble(metric = metric, measure = unname(dna_measures[[v]]),
             depth = if (dp == "A") "surface" else "subsoil", n = nrow(d),
             pearson_r = unname(ct$estimate), p_value = ct$p.value)
    }))
  })) |> mutate(p_bh = p.adjust(p_value, "BH"))
}))
write_csv(colonization_tests, file.path(TABLES_DIR, "colonization_vs_soil_dna.csv"))

## ---- D. understory mycorrhizal composition across the latitudinal gradient --
## Does AM-host cover fall away toward the north? This bears on whether
## stand-level fine-root biomass indexes AM host density consistently across
## the gradient.

## Latitude is a node-level predictor: the three stands at a node share the
## same latitude and must not be counted as independent replicates. Average
## the stand-level vegetation measurements within each of the nine nodes.
understory <- env |>
  distinct(node_stand, .keep_all = TRUE) |>
  group_by(node) |>
  summarise(
    latitude = mean(latitude, na.rm = TRUE),
    across(c(proportion_am, proportion_em, total_cover),
           ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

understory_tests <- bind_rows(lapply(
  names(c(proportion_am = 1, proportion_em = 1, total_cover = 1)), function(v) {
    d <- understory |>
      select(latitude, all_of(v)) |>
      filter(complete.cases(across(everything())))
    ct <- cor.test(d$latitude, d[[v]], method = "pearson")
    tibble(
      measure = c(proportion_am = "AM-associated cover (proportion)",
                  proportion_em = "ECM-associated cover (proportion)",
                  total_cover = "Total understory cover")[[v]],
      n_nodes = nrow(d),
      pearson_r = unname(ct$estimate),
      p_value = ct$p.value
    )
  }))
write_csv(understory_tests, file.path(TABLES_DIR, "understory_vs_latitude.csv"))

cat("\nGuild:ECM log-ratio depth models\n")
print(ratio_depth_models, n = Inf, width = Inf)
cat("\nConfidence-rank sensitivity (probable + highly probable)\n")
print(conf_inference, n = Inf, width = Inf)
cat("\nRoot ECM colonization versus bulk-soil DNA measures\n")
print(colonization_tests, n = Inf, width = Inf)
cat("\nUnderstory composition versus latitude\n")
print(understory_tests, n = Inf, width = Inf)

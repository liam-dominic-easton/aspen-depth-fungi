## 10_guild_ambiguity_sensitivity.R
## Test whether focal ITS guild results depend on the priority rule used for
## FUNGuild annotations containing more than one focal guild.
##
## Outputs (outputs/tables/):
##   guild_ambiguity_sample_proportions.csv  per-sample guild reads under three
##                                           schemes; read by 06, 12, 14, 16-18
##   guild_ambiguity_depth_inference.csv     depth contrasts per scheme (Table S10)
##   guild_ambiguity_method_comparison.csv   agreement between schemes
##   guild_ambiguity_asv_summary.csv         guild flags per ASV

source("scripts/helpers.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})
source("scripts/bayes_helpers.R")

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
sample_cols <- setdiff(names(fg), meta_cols)

ps_its <- load_ps("ITS")
key <- get_sample_key(ps_its) |>
  transmute(
    sample_name,
    sample_id_clean,
    node,
    stand,
    depth,
    node_stand = paste0(node, stand),
    depth_num = if_else(depth == "B", 1, 0)
  )
sample_cols <- intersect(sample_cols, key$sample_name)
if (length(sample_cols) != nrow(key)) {
  stop("FUNGuild sample columns do not match all clean ITS samples")
}

fg_flags <- fg |>
  mutate(
    Guild = str_squish(replace_na(Guild, "")),
    guild_lower = str_to_lower(Guild),
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

count_mat <- as.matrix(fg_flags[, sample_cols])
storage.mode(count_mat) <- "numeric"
fg_flags$total_reads_clean <- rowSums(count_mat)

asv_summary <- fg_flags |>
  transmute(
    OTU_ID, Guild, `Confidence Ranking`, total_reads_clean,
    n_focal_flags, priority_category,
    ECM = flag_ecm, AMF = flag_amf, Saprotroph = flag_saprotroph,
    Pathogen = flag_pathogen, Endophyte = flag_endophyte,
    ambiguous = n_focal_flags > 1
  )
write_csv(asv_summary, file.path(TABLES_DIR, "guild_ambiguity_asv_summary.csv"))

cat("FUNGuild ASVs:", nrow(fg_flags), "\n")
cat("Ambiguous focal assignments:", sum(fg_flags$n_focal_flags > 1), "ASVs;",
    sum(fg_flags$total_reads_clean[fg_flags$n_focal_flags > 1]), "reads\n")

counts_long <- fg_flags |>
  select(OTU_ID, priority_category, n_focal_flags,
         flag_ecm, flag_amf, flag_saprotroph, flag_pathogen, flag_endophyte,
         all_of(sample_cols)) |>
  pivot_longer(all_of(sample_cols), names_to = "sample_name", values_to = "reads") |>
  filter(reads > 0) |>
  left_join(key, by = "sample_name")

focal <- c("ECM", "Saprotroph", "Pathogen")
flag_name <- c(ECM = "flag_ecm", Saprotroph = "flag_saprotroph", Pathogen = "flag_pathogen")
sample_totals <- counts_long |>
  group_by(sample_name) |>
  summarise(total_its_reads = sum(reads), .groups = "drop")

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

sample_props <- bind_rows(priority, exclusive, fractional) |>
  complete(sample_name = key$sample_name, guild = focal,
           method = c("priority", "exclusive_only", "fractional"),
           fill = list(guild_reads = 0)) |>
  left_join(sample_totals, by = "sample_name") |>
  left_join(key, by = "sample_name") |>
  mutate(proportion = guild_reads / total_its_reads) |>
  arrange(guild, method, node_stand, depth)

write_csv(sample_props, file.path(TABLES_DIR, "guild_ambiguity_sample_proportions.csv"))

## Verify that the reconstructed priority rule exactly reproduces the guild
## proportions in depth_analysis_data.csv.
env_current <- read_csv(file.path(ENV_DIR, "depth_analysis_data.csv"), show_col_types = FALSE) |>
  select(sample_id_clean, prop_ecm, prop_saprotroph, prop_plant_pathogen)
priority_check <- sample_props |>
  filter(method == "priority") |>
  select(sample_id_clean, guild, proportion) |>
  pivot_wider(names_from = guild, values_from = proportion) |>
  left_join(env_current, by = "sample_id_clean") |>
  mutate(
    diff_ecm = ECM - prop_ecm,
    diff_saprotroph = Saprotroph - prop_saprotroph,
    diff_pathogen = Pathogen - prop_plant_pathogen
  )
max_diff <- max(abs(unlist(priority_check[, c("diff_ecm", "diff_saprotroph", "diff_pathogen")])),
                na.rm = TRUE)
cat("Maximum difference from the priority-rule proportions in depth_analysis_data.csv:", signif(max_diff, 4), "\n")
if (!is.finite(max_diff) || max_diff > 1e-10) {
  stop("Priority-rule reconstruction does not reproduce the guild proportions in depth_analysis_data.csv")
}

## Depth inference uses the Bayesian paired model used throughout
## (bayes_helpers.R), alongside descriptive paired means/median differences.
## Differences are subsoil - surface in percentage points; bayes_effect_* is on
## the logit scale.
analyse_one <- function(df) {
  df <- arrange(df, node_stand, depth_num)
  eps <- 0.5 / pmax(df$total_its_reads, 1)
  df$logit_prop <- qlogis(pmin(pmax(df$proportion, eps), 1 - eps))
  paired_summary_bayes(df, "logit_prop", "proportion", scale = 100)
}

inference <- sample_props |>
  group_by(guild, method) |>
  group_modify(~ analyse_one(.x)) |>
  ungroup() |>
  arrange(guild, factor(method, levels = c("priority", "exclusive_only", "fractional")))
write_csv(inference, file.path(TABLES_DIR, "guild_ambiguity_depth_inference.csv"))

comparison <- sample_props |>
  select(sample_name, guild, method, proportion) |>
  pivot_wider(names_from = method, values_from = proportion) |>
  group_by(guild) |>
  summarise(
    priority_vs_exclusive_spearman = cor(priority, exclusive_only, method = "spearman"),
    priority_vs_fractional_spearman = cor(priority, fractional, method = "spearman"),
    mean_absolute_change_exclusive_pp = 100 * mean(abs(exclusive_only - priority)),
    mean_absolute_change_fractional_pp = 100 * mean(abs(fractional - priority)),
    .groups = "drop"
  )
write_csv(comparison, file.path(TABLES_DIR, "guild_ambiguity_method_comparison.csv"))

print(inference, n = Inf, width = Inf)
print(comparison, n = Inf, width = Inf)
cat("Saved guild-ambiguity sensitivity outputs in", TABLES_DIR, "\n")

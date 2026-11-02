## 15_pathogen_taxon_audit.R -- taxa underlying the pathogen depth signal
## (Table S15). Reads outputs/tables/indicator_by_depth_full.csv from 03.
##
## This exploratory audit applies the manuscript's priority rule, aggregates
## pathogen-assigned ASVs by genus, and tests paired surface/subsoil changes.
## It also joins the complete ASV-level indicator analysis so that raw and
## BH-adjusted indicator evidence are reported rather than inferred.

source("scripts/helpers.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

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
    node_stand = paste0(node, stand)
  )
sample_cols <- intersect(setdiff(names(fg), meta_cols), key$sample_name)
if (length(sample_cols) != nrow(key)) {
  stop("FUNGuild sample columns do not match the clean ITS sample set")
}

its_tax <- as.data.frame(as(tax_table(ps_its), "matrix")) |>
  tibble::rownames_to_column("OTU_ID") |>
  transmute(
    OTU_ID,
    genus = str_remove(as.character(Genus), "^[gG]__"),
    species = str_remove(as.character(Species), "^[sS]__")
  )

fg_flags <- fg |>
  mutate(
    guild_text = str_to_lower(str_squish(replace_na(Guild, ""))),
    is_unassigned = Guild %in% c("", "-") |
      str_detect(guild_text, "^\\s*(unassigned|unknown|undefined|none)\\s*$"),
    flag_ecm = str_detect(guild_text, "ectomycorrhizal"),
    flag_amf = str_detect(guild_text, "arbuscular"),
    flag_saprotroph = str_detect(guild_text, "saprotroph"),
    flag_pathogen = str_detect(guild_text, "pathogen"),
    flag_plant_pathogen = str_detect(guild_text, "plant pathogen"),
    flag_endophyte = str_detect(guild_text, "endophyte"),
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
  ) |>
  left_join(its_tax, by = "OTU_ID") |>
  mutate(
    genus = if_else(is.na(genus) | genus == "" | genus == "unidentified",
                    "Unidentified", genus)
  )

counts_long <- fg_flags |>
  select(OTU_ID, genus, species, priority_category, flag_plant_pathogen,
         all_of(sample_cols)) |>
  pivot_longer(all_of(sample_cols), names_to = "sample_name", values_to = "reads") |>
  mutate(reads = as.numeric(reads)) |>
  left_join(key, by = "sample_name")

sample_totals <- counts_long |>
  group_by(sample_name) |>
  summarise(total_its_reads = sum(reads), .groups = "drop")

pathogen_asvs <- fg_flags |>
  filter(priority_category == "Pathogen") |>
  select(OTU_ID, genus, species, flag_plant_pathogen)

genus_props <- counts_long |>
  filter(priority_category == "Pathogen") |>
  group_by(sample_name, sample_id_clean, node_stand, depth, genus) |>
  summarise(reads = sum(reads), .groups = "drop") |>
  complete(
    nesting(sample_name, sample_id_clean, node_stand, depth),
    genus = unique(pathogen_asvs$genus),
    fill = list(reads = 0)
  ) |>
  left_join(sample_totals, by = "sample_name") |>
  mutate(proportion = reads / total_its_reads)

summarize_genus <- function(df) {
  paired <- df |>
    select(node_stand, depth, proportion) |>
    pivot_wider(names_from = depth, values_from = proportion) |>
    filter(!is.na(A), !is.na(B))
  wt <- suppressWarnings(wilcox.test(
    paired$B, paired$A, paired = TRUE, exact = FALSE
  ))
  tibble(
    paired_stands = nrow(paired),
    surface_mean_percent = 100 * mean(paired$A),
    subsoil_mean_percent = 100 * mean(paired$B),
    change_percentage_points = 100 * mean(paired$B - paired$A),
    stands_increasing = sum(paired$B > paired$A),
    stands_decreasing = sum(paired$B < paired$A),
    ties = sum(paired$B == paired$A),
    wilcoxon_p = wt$p.value
  )
}

genus_summary <- genus_props |>
  group_by(genus) |>
  group_modify(~ summarize_genus(.x)) |>
  ungroup() |>
  mutate(wilcoxon_bh = p.adjust(wilcoxon_p, method = "BH"),
         ## 34 tested groups = named genera + an unidentified pool + higher-rank
         ## placeholders
         group_type = case_when(genus == "Unidentified" ~ "unidentified pool",
                                grepl("incertae", genus, ignore.case = TRUE) ~ "higher-rank placeholder",
                                TRUE ~ "named genus"))

read_and_prevalence <- genus_props |>
  group_by(genus, depth) |>
  summarise(
    prevalence = sum(reads > 0),
    reads = sum(reads),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = depth,
    values_from = c(reads, prevalence),
    names_glue = "{.value}_{depth}"
  )

asv_counts <- pathogen_asvs |>
  group_by(genus) |>
  summarise(
    n_asvs = n(),
    n_plant_pathogen_asvs = sum(flag_plant_pathogen),
    species_labels = paste(sort(unique(na.omit(species))), collapse = "; "),
    .groups = "drop"
  )

indicator_full <- read_csv(
  file.path(TABLES_DIR, "indicator_by_depth_full.csv"),
  show_col_types = FALSE
) |>
  select(taxon, indicator_for, indicator_stat = stat,
         indicator_p = p.value, indicator_bh = p_bh)

indicator_summary <- pathogen_asvs |>
  inner_join(indicator_full, by = c("OTU_ID" = "taxon")) |>
  filter(indicator_for == "subsoil") |>
  group_by(genus) |>
  summarise(
    n_raw_subsoil_indicator_asvs = sum(indicator_p < 0.05),
    minimum_indicator_p = min(indicator_p),
    minimum_indicator_bh = min(indicator_bh),
    .groups = "drop"
  )

output <- genus_summary |>
  left_join(read_and_prevalence, by = "genus") |>
  left_join(asv_counts, by = "genus") |>
  left_join(indicator_summary, by = "genus") |>
  mutate(
    n_raw_subsoil_indicator_asvs = replace_na(n_raw_subsoil_indicator_asvs, 0L)
  ) |>
  arrange(wilcoxon_p)

write_csv(output, file.path(TABLES_DIR, "pathogen_genera_depth_audit.csv"))

cat("Priority-rule pathogen ASVs:", nrow(pathogen_asvs), "in",
    n_distinct(pathogen_asvs$genus), "genera\n")
cat("Plant-pathogen-labelled among them:",
    sum(pathogen_asvs$flag_plant_pathogen), "ASVs\n\n")
print(as_tibble(output), n = Inf, width = Inf)

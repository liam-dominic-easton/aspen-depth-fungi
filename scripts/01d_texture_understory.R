## 01d_texture_understory.R - Soil texture and understory mycorrhizal composition
## from the raw source files.
##
## Inputs:
##   data/raw_environmental/Soil texture.xls     particle-size analysis (laser
##                                               diffraction), one column per
##                                               sample, rows "Sand", "Silt",
##                                               "Clay" (% of the < 2 mm fraction)
##   data/raw_environmental/Plant understorey.xlsx  % cover of each species in
##                                               two 1-m2 quadrats per stand
##   data/reference/meng2023_5_11770.csv         plant mycorrhizal types of Meng
##                                               et al. (2023), Ecol. Lett. 26:1862;
##                                               Dryad doi:10.5061/dryad.x0k6djhq2
##   data/reference/understory_taxon_names.csv   name used for each understory taxon
##                                               in the Meng table (species, or genus /
##                                               family for taxa not identified to
##                                               species or not in the table)
##
## Outputs: overwrites sand, silt, clay, proportion_am, proportion_em and
## proportion_nm in data/environmental/depth_analysis_data.csv;
## outputs/tables/understory_mycorrhizal_types.csv (type assigned to each taxon).
##
## Texture values are matched to samples by name (the raw file writes some
## names with a space, e.g. "RED 2A", and some without, e.g. "RED1A").

source("scripts/helpers.R")
suppressPackageStartupMessages({
  library(readxl)
  library(tidyr)
})

env_file <- file.path(ENV_DIR, "depth_analysis_data.csv")
env <- read_csv(env_file, show_col_types = FALSE)
env_key <- paste0(env$node, env$stand, env$depth)

## ---- texture -----------------------------------------------------------------
tx <- suppressMessages(read_excel(file.path(DATA_DIR, "raw_environmental", "Soil texture.xls"),
                                  sheet = "Data Table", col_names = FALSE))
lab <- trimws(as.character(tx[[1]]))
row_of <- function(k) { i <- which(lab == k); stopifnot(length(i) == 1); i }
texture <- tibble(
  sample = gsub(" ", "", unlist(tx[2, -1])),
  sand   = as.numeric(unlist(tx[row_of("Sand"), -1])),
  silt   = as.numeric(unlist(tx[row_of("Silt"), -1])),
  clay   = as.numeric(unlist(tx[row_of("Clay"), -1])),
  total  = as.numeric(unlist(tx[row_of("Total (%)"), -1]))
)
stopifnot(!anyDuplicated(texture$sample), all(abs(texture$total - 100) < 0.01),
          all(env_key %in% texture$sample))
m <- match(env_key, texture$sample)
changed <- abs(env$sand - texture$sand[m]) > 0.01
cat("Sand values replaced:", sum(changed), "of", nrow(env),
    if (any(changed)) paste0("(", paste(env_key[changed], collapse = ", "), ")") else "", "\n")
env$sand <- texture$sand[m]
env$silt <- texture$silt[m]
env$clay <- texture$clay[m]

## ---- understory mycorrhizal composition ----------------------------------------
## Each taxon takes its mycorrhizal type from Meng et al. (2023), who group
## dual AM + ECM species with ECM. Taxa identified only to genus, or absent
## from the table, take the most common type among their congeners (Poaceae for
## unidentified grasses). Types are grouped as AM, ECM, or other (non-mycorrhizal,
## ericoid, orchid); proportion_nm holds the "other" share. Shares are of the
## summed cover of both quadrats in a stand.
meng <- suppressMessages(read_csv(file.path(DATA_DIR, "reference", "meng2023_5_11770.csv"), show_col_types = FALSE)) |>
  transmute(species = gsub("_", " ", Species), Genus, Family, type = Mycorrhizal_type)
names_tab <- read_csv(file.path(DATA_DIR, "reference", "understory_taxon_names.csv"), show_col_types = FALSE)
modal <- function(x) { t <- sort(table(x), decreasing = TRUE); names(t)[1] }
types <- names_tab |>
  rowwise() |>
  mutate(
    pool = list(switch(level,
      species = meng$type[meng$species == lookup],
      genus   = meng$type[meng$Genus %in% strsplit(lookup, "|", fixed = TRUE)[[1]]],
      family  = meng$type[meng$Family == lookup])),
    n_records = length(pool),
    meng_type = if (n_records) modal(pool) else NA_character_,
    type_counts = paste(names(table(pool)), table(pool), sep = " ", collapse = ", ")
  ) |>
  ungroup() |>
  select(-pool) |>
  mutate(category = case_when(meng_type == "AM" ~ "AM", meng_type == "ECM" ~ "ECM", TRUE ~ "other"))
if (anyNA(types$meng_type)) stop("No Meng et al. record for: ",
                                 paste(types$taxon[is.na(types$meng_type)], collapse = ", "))
write_csv(types, file.path(TABLES_DIR, "understory_mycorrhizal_types.csv"))

us <- suppressMessages(read_excel(file.path(DATA_DIR, "raw_environmental", "Plant understorey.xlsx"),
                                  sheet = "raw data"))
taxa <- setdiff(names(us), c("Node", "Stand", "Square metre"))
missing <- setdiff(taxa, types$taxon)
if (length(missing)) stop("No name mapping for: ", paste(missing, collapse = ", "))
comp <- us |>
  pivot_longer(all_of(taxa), names_to = "taxon", values_to = "cover") |>
  filter(!is.na(cover), cover > 0) |>
  left_join(types |> select(taxon, category), by = "taxon") |>
  group_by(node = Node, stand = as.numeric(Stand)) |>
  summarise(total = sum(cover),
            proportion_am = sum(cover[category == "AM"]) / total,
            proportion_em = sum(cover[category == "ECM"]) / total,
            proportion_nm = sum(cover[category == "other"]) / total,
            .groups = "drop")
stopifnot(all(abs(comp$proportion_am + comp$proportion_em + comp$proportion_nm - 1) < 1e-9))
## total cover (mean of the two quadrats) must match the column already in use
chk <- env |> distinct(node, stand, total_cover) |> left_join(comp, by = c("node", "stand"))
stopifnot(all(abs(chk$total / 2 - chk$total_cover) < 1e-6))
k <- match(paste(env$node, env$stand), paste(comp$node, comp$stand))
for (v in c("proportion_am", "proportion_em", "proportion_nm")) env[[v]] <- comp[[v]][k]
cat("Understory mycorrhizal proportions recomputed for", nrow(comp), "stands.\n")

write_csv(env, env_file)
cat("Saved:", env_file, "\n")

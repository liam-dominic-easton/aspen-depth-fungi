## 01b_build_its_tables.R — build all ITS-derived analysis inputs from the
## DADA2 output of 01_dada2_its.R.
##
## Inputs:  data/derived/its_seqtab_nochim.rds, data/derived/its_taxonomy.rds
##          data/reference/funguild_db_2026-09-26.json
##          data/legacy_inputs/ps_fungi_all.rds (sample_data only)
##          data/legacy_inputs/depth_analysis_data.csv
## Outputs: data/taxonomy/otu_table_for_FUNGuild.guilds.txt
##          data/phyloseq/ps_fungi_all.rds, ps_fungi_clean.rds
##          data/environmental/depth_analysis_data.csv (ITS columns recomputed)
## If data/legacy_inputs/ does not exist, it is created on first run and the
## starting files in data/taxonomy, data/phyloseq and data/environmental are
## copied into it before they are overwritten (see README).
##
## CHECK_MODE=1 recomputes the ITS columns from the FUNGuild table in
## data/legacy_inputs/ and verifies they reproduce data/legacy_inputs/depth_analysis_data.csv.

source("scripts/helpers.R")
source("scripts/funguild_helpers.R")
suppressPackageStartupMessages(library(tidyr))

ARCH <- file.path(DATA_DIR, "legacy_inputs")
if (!dir.exists(ARCH)) {
  dir.create(ARCH)
  file.copy(file.path(TAX_DIR, "otu_table_for_FUNGuild.guilds.txt"), ARCH)
  file.copy(file.path(PHY_DIR, c("ps_fungi_all.rds", "ps_fungi_clean.rds")), ARCH)
  file.copy(file.path(ENV_DIR, "depth_analysis_data.csv"), ARCH)
}
env_legacy <- read_csv(file.path(ARCH, "depth_analysis_data.csv"), show_col_types = FALSE)
ps_legacy  <- readRDS(file.path(ARCH, "ps_fungi_all.rds"))
key <- get_sample_key(ps_legacy)

## per-sample ITS summaries under the guild priority rule (guild_category)
its_summaries <- function(counts, guild) {
  cat_asv <- guild_category(guild)
  cats <- c("ECM", "AMF", "Saprotroph", "Plant pathogen", "Endophyte", "Other", "Unassigned")
  bind_rows(lapply(rownames(counts), function(s) {
    x <- counts[s, ]; tot <- sum(x)
    reads <- sapply(cats, function(k) sum(x[cat_asv == k]))
    rich  <- sapply(cats, function(k) sum(x[cat_asv == k] > 0))
    tibble(sample_name = s, its_total_reads = tot, its_asv_richness = sum(x > 0),
           prop_ecm = reads[["ECM"]] / tot, prop_saprotroph = reads[["Saprotroph"]] / tot,
           prop_plant_pathogen = reads[["Plant pathogen"]] / tot, prop_amf = reads[["AMF"]] / tot,
           prop_endophyte = reads[["Endophyte"]] / tot, prop_unassigned = reads[["Unassigned"]] / tot,
           richness_ECM = rich[["ECM"]], richness_AMF = rich[["AMF"]], richness_Saprotroph = rich[["Saprotroph"]],
           `richness_Plant pathogen` = rich[["Plant pathogen"]], richness_Endophyte = rich[["Endophyte"]],
           richness_Other = rich[["Other"]], richness_Unassigned = rich[["Unassigned"]])
  }))
}
its_cols <- c("its_total_reads", "its_asv_richness", "prop_ecm", "prop_saprotroph", "prop_plant_pathogen",
              "prop_amf", "prop_endophyte", "prop_unassigned", "richness_ECM", "richness_AMF",
              "richness_Saprotroph", "richness_Plant pathogen", "richness_Endophyte", "richness_Other",
              "richness_Unassigned")
update_env <- function(env, summ) {
  summ <- summ |> left_join(key |> select(sample_name, sample_id_clean), by = "sample_name")
  stopifnot(!anyNA(summ$sample_id_clean))
  env <- env |> select(-any_of(its_cols)) |> left_join(summ |> select(sample_id_clean, all_of(its_cols)),
                                                       by = "sample_id_clean")
  env
}

if (Sys.getenv("CHECK_MODE") == "1") {
  fg <- read.delim(file.path(ARCH, "otu_table_for_FUNGuild.guilds.txt"), check.names = FALSE, stringsAsFactors = FALSE)
  counts <- t(as.matrix(fg[, key$sample_name])); colnames(counts) <- fg$OTU_ID
  chk <- update_env(env_legacy, its_summaries(counts, fg$Guild))
  diffs <- sapply(its_cols, function(cc) max(abs(chk[[cc]] - env_legacy[[cc]])))
  print(signif(diffs, 3))
  stopifnot(all(diffs < 1e-9))
  cat("CHECK_MODE: legacy ITS columns reproduced exactly.\n")
  quit(save = "no")
}

## ---- rebuild from the DADA2 output ---------------------------------------------
st  <- readRDS(file.path(DATA_DIR, "derived", "its_seqtab_nochim.rds"))
tax <- readRDS(file.path(DATA_DIR, "derived", "its_taxonomy.rds"))$tax
rownames(st) <- paste0(rownames(st), "_R1_001.fastq.gz")          # sample naming of the phyloseq sample_data
st <- st[rownames(st) %in% key$sample_name, , drop = FALSE]       # 52 biological samples
st <- st[match(key$sample_name, rownames(st)), , drop = FALSE]
stopifnot(identical(rownames(st), key$sample_name))
st <- st[, colSums(st) > 0, drop = FALSE]
tax <- tax[colnames(st), , drop = FALSE]
cat("ITS table:", nrow(st), "samples x", ncol(st), "ASVs;", sum(st), "reads\n")

db <- read_funguild_db(file.path(DATA_DIR, "reference", "funguild_db_2026-09-26.json"))
fgn <- funguild_assign(tax, db)
taxonomy_string <- apply(tax, 1, function(r) paste(ifelse(is.na(r), "", r), collapse = ";"))
fg_out <- data.frame(OTU_ID = colnames(st), t(st), taxonomy = taxonomy_string, fgn,
                     check.names = FALSE, stringsAsFactors = FALSE)
write.table(fg_out, file.path(TAX_DIR, "otu_table_for_FUNGuild.guilds.txt"), sep = "\t", quote = FALSE,
            row.names = FALSE)

sd <- sample_data(ps_legacy)[key$sample_name, ]
ps_new <- phyloseq(otu_table(st, taxa_are_rows = FALSE), tax_table(tax), sd)
saveRDS(ps_new, file.path(PHY_DIR, "ps_fungi_all.rds"))
saveRDS(ps_new, file.path(PHY_DIR, "ps_fungi_clean.rds"))

## replace only the ITS columns of the current table (01c updates the 18S columns)
env_cur <- read_csv(file.path(ENV_DIR, "depth_analysis_data.csv"), show_col_types = FALSE)
env_new <- update_env(env_cur, its_summaries(st, fgn$Guild))[, names(env_cur)]
write_csv(env_new, file.path(ENV_DIR, "depth_analysis_data.csv"))

cmp <- env_legacy |> select(sample_id_clean, its_total_reads, prop_ecm) |>
  left_join(env_new |> select(sample_id_clean, its_total_reads, prop_ecm), by = "sample_id_clean",
            suffix = c("_legacy", "_rebuilt"))
print(summary(cmp$its_total_reads_rebuilt / cmp$its_total_reads_legacy))
print(cmp |> filter(sample_id_clean == "MOO3B"))
cat("Guild categories (ASVs):\n"); print(table(guild_category(fgn$Guild)))
cat("01b complete.\n")

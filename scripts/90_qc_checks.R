## 90_qc_checks.R — data-derived numbers quoted in Methods/Results.
## Output: outputs/tables/qc_reported_numbers.csv

source("scripts/helpers.R")
ps_its <- load_ps("ITS"); ps_18s <- load_ps("18S")
its <- as(otu_table(ps_its), "matrix"); if (taxa_are_rows(ps_its)) its <- t(its)
amf <- as(otu_table(ps_18s), "matrix"); if (taxa_are_rows(ps_18s)) amf <- t(amf)
fg <- read.delim(file.path(TAX_DIR, "otu_table_for_FUNGuild.guilds.txt"), stringsAsFactors = FALSE)
fg <- fg[match(colnames(its), fg$OTU_ID), ]
stopifnot(identical(fg$OTU_ID, colnames(its)))
has_guild <- !(fg$Guild %in% c("-", "", NA))
prob <- fg$Confidence.Ranking %in% c("Probable", "Highly Probable")
key <- get_sample_key(ps_its)
low <- key$sample_id_clean[match(names(sort(rowSums(its)))[1:3], key$sample_name)]
env <- load_env()
qc <- tibble::tribble(
  ~quantity, ~value,
  "ITS samples", nrow(its),
  "ITS ASVs", ncol(its),
  "ITS reads per sample, median", median(rowSums(its)),
  "ITS reads per sample, min", min(rowSums(its)),
  "ITS reads per sample, max", max(rowSums(its)),
  "18S samples", nrow(amf),
  "18S reads per sample, median", median(rowSums(amf)),
  "18S reads per sample, mean", mean(rowSums(amf)),
  "ITS reads with any FUNGuild guild (%)", 100 * sum(its[, has_guild]) / sum(its),
  "ITS reads at probable/highly probable (%)", 100 * sum(its[, has_guild & prob]) / sum(its),
  "Fine roots surface:subsoil ratio of means", mean(env$fine_roots[env$depth == "A"]) / mean(env$fine_roots[env$depth == "B"]),
  "Latitude min", min(env$latitude), "Latitude max", max(env$latitude),
  "Aspen age (stand mean) min", min(env$tree_age), "Aspen age (stand mean) max", max(env$tree_age),
  "Aspen age mean", mean(unique(env[, c("node", "stand", "tree_age")])$tree_age)
)
print(qc, n = Inf)
cat("Lowest-read ITS samples:", paste(low, collapse = ", "), "\n")
write_csv(qc, file.path(TABLES_DIR, "qc_reported_numbers.csv"))

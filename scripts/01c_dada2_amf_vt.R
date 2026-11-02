## 01c_dada2_amf_vt.R — 18S (AMF-targeted) processing from raw FASTQs and
## MaarjAM virtual-taxon (VT) assignment.
##
## Inputs (not shipped): AMF_READS = folder with the 53 18S read pairs
##   (*-amf_*_R1_001.fastq.gz; forward reads only are used; deposited in
##   [SRA accession to be added]).
## Also reads data/legacy_inputs/ps_amf_all.rds (sample_data only) and
##   data/legacy_inputs/depth_analysis_data.csv; missing files are copied there
##   from data/taxonomy, data/phyloseq and data/environmental on first run.
## Requires vsearch 2.28.1 on the PATH (VSEARCH env var to override).
## Settings: DADA2 1.34.0; forward reads; NS31 removed (trimLeft = 21); no
## truncation; maxEE = 2; truncQ = 2; maxN = 0; rm.phix; minLen 50; error model
## from all filtered files; per-sample inference; consensus chimera removal.
## VT: 22-nt conserved leader (AGCAGCCGCGGTAATTCCAGCT) removed when present;
## vsearch --usearch_global --id 0.97 --query_cov 0.90 --maxaccepts 5 against the
## 2021 MaarjAM SSU database; best-identity hits kept; ties resolved to the
## alphabetically first VT (flagged ambiguous).
## Outputs: data/derived/amf_seqtab_nochim.rds, amf_dada2_track.csv,
##          data/taxonomy/easton_18s_asv_map.csv, maarjam_vsearch_mapping_final.csv,
##          data/phyloseq/ps_amf_all.rds, ps_amf_euk_clean.rds,
##          data/environmental/depth_analysis_data.csv (18S columns recomputed).
## CHECK_MODE=1 verifies the 18S column computation on the inputs in data/legacy_inputs/.

source("scripts/helpers.R")
suppressPackageStartupMessages(library(tidyr))

ARCH <- file.path(DATA_DIR, "legacy_inputs"); dir.create(ARCH, showWarnings = FALSE)
for (f in c(file.path(TAX_DIR, c("easton_18s_asv_map.csv", "maarjam_vsearch_mapping_final.csv")),
            file.path(PHY_DIR, c("ps_amf_all.rds", "ps_amf_euk_clean.rds"))))
  if (!file.exists(file.path(ARCH, basename(f)))) file.copy(f, ARCH)
if (!file.exists(file.path(ARCH, "depth_analysis_data.csv")))
  file.copy(file.path(ENV_DIR, "depth_analysis_data.csv"), ARCH)
env_legacy <- read_csv(file.path(ARCH, "depth_analysis_data.csv"), show_col_types = FALSE)
ps_legacy  <- readRDS(file.path(ARCH, "ps_amf_all.rds"))
key <- get_sample_key(ps_legacy)

amf_cols <- c("amf_total_reads", "amf_asv_richness", "amf_reads", "non_amf_reads", "b_amf", "vt_richness")
amf_summaries <- function(counts, vt_of_seq) {
  vt <- vt_of_seq[colnames(counts)]; is_amf <- !is.na(vt) & vt != ""
  tibble(sample_name = rownames(counts),
         amf_total_reads = rowSums(counts), amf_asv_richness = rowSums(counts > 0),
         amf_reads = rowSums(counts[, is_amf, drop = FALSE]),
         non_amf_reads = rowSums(counts[, !is_amf, drop = FALSE]),
         b_amf = log(amf_reads / non_amf_reads),
         vt_richness = apply(counts[, is_amf, drop = FALSE], 1, function(x) length(unique(vt[is_amf][x > 0]))))
}
update_env <- function(env, summ) {
  summ <- summ |> left_join(key |> select(sample_name, sample_id_clean), by = "sample_name")
  stopifnot(!anyNA(summ$sample_id_clean))
  env |> select(-any_of(amf_cols)) |> left_join(summ |> select(sample_id_clean, all_of(amf_cols)), by = "sample_id_clean")
}
counts_of <- function(ps) { m <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) m <- t(m); m }

if (Sys.getenv("CHECK_MODE") == "1") {
  map <- read_csv(file.path(ARCH, "easton_18s_asv_map.csv"), show_col_types = FALSE)
  vtm <- read_csv(file.path(ARCH, "maarjam_vsearch_mapping_final.csv"), show_col_types = FALSE)
  vt_of_seq <- setNames(vtm$vtx_vsearch[match(map$asv_id, vtm$asv_id)], map$sequence)
  chk <- update_env(env_legacy, amf_summaries(counts_of(ps_legacy)[key$sample_name, ], vt_of_seq))
  d <- sapply(amf_cols, function(cc) max(abs(chk[[cc]] - env_legacy[[cc]]), na.rm = TRUE)); print(signif(d, 3))
  stopifnot(all(d < 1e-9)); cat("CHECK_MODE: legacy 18S columns reproduced exactly.\n"); quit(save = "no")
}

## ---- DADA2 -----------------------------------------------------------------------
suppressPackageStartupMessages(library(dada2))
READS <- Sys.getenv("AMF_READS"); stopifnot(dir.exists(READS))
OUT <- file.path(DATA_DIR, "derived"); dir.create(file.path(OUT, "filtered_amf"), recursive = TRUE, showWarnings = FALSE)
fnFs <- sort(list.files(READS, pattern = "amf.*_R1_001\\.fastq\\.gz$", full.names = TRUE))
sn <- sub("_R1_001\\.fastq\\.gz$", "", basename(fnFs))
filtFs <- file.path(OUT, "filtered_amf", basename(fnFs))
ft <- filterAndTrim(fnFs, filtFs, trimLeft = 21, truncLen = 0, maxN = 0, maxEE = 2, truncQ = 2,
                    rm.phix = TRUE, minLen = 50, compress = TRUE, multithread = TRUE)
keep <- file.exists(filtFs)
errF <- learnErrors(filtFs[keep], multithread = TRUE)
derep <- derepFastq(filtFs[keep]); names(derep) <- sn[keep]
dd <- dada(derep, err = errF, multithread = TRUE)
st <- removeBimeraDenovo(makeSequenceTable(dd), method = "consensus", multithread = TRUE)
saveRDS(st, file.path(OUT, "amf_seqtab_nochim.rds"))
getN <- function(x) sum(getUniques(x))
write.csv(data.frame(sample = sn[keep], ft[keep, , drop = FALSE], denoised = sapply(dd, getN), nonchim = rowSums(st)),
          file.path(OUT, "amf_dada2_track.csv"), row.names = FALSE)

## ---- biological samples, ASV map ---------------------------------------------------
st <- st[rownames(st) %in% key$sample_name, , drop = FALSE]
st <- st[match(key$sample_name, rownames(st)), , drop = FALSE]
stopifnot(identical(rownames(st), key$sample_name))
st <- st[, colSums(st) > 0, drop = FALSE]
st <- st[, order(-colSums(st)), drop = FALSE]
map <- tibble(asv_id = paste0("ASV_", seq_len(ncol(st))), sequence = colnames(st),
              total_reads = colSums(st), length = nchar(colnames(st)))
write_csv(map, file.path(TAX_DIR, "easton_18s_asv_map.csv"))

## ---- vsearch against MaarjAM ------------------------------------------------------------
leader <- "AGCAGCCGCGGTAATTCCAGCT"
q <- ifelse(startsWith(map$sequence, leader), substring(map$sequence, nchar(leader) + 1), map$sequence)
qfa <- file.path(OUT, "amf_asvs_trimmed.fa"); hits_file <- file.path(OUT, "vsearch_maarjam_hits.txt")
writeLines(paste0(">", map$asv_id, "\n", q), qfa)
ref <- file.path(DATA_DIR, "reference", "maarjam_SSU_type_2021", "maarjam_database_SSU.fasta")
vs <- Sys.getenv("VSEARCH", "vsearch")
cat(system2(vs, "--version", stdout = TRUE, stderr = TRUE)[1], "\n")
status <- system2(vs, c("--usearch_global", shQuote(qfa), "--db", shQuote(ref), "--id", "0.97",
                        "--blast6out", shQuote(hits_file), "--query_cov", "0.90", "--maxaccepts", "5",
                        "--threads", "8"), stdout = FALSE, stderr = FALSE)
stopifnot(status == 0)
hdr <- grep("^>", readLines(ref), value = TRUE)
id_to_vtx <- setNames(regmatches(hdr, regexpr("VTX[0-9]+", hdr)), sub("^>(\\S+).*", "\\1", hdr)[grepl("VTX[0-9]+", hdr)])
h <- read.delim(hits_file, header = FALSE, stringsAsFactors = FALSE)[, 1:3]
names(h) <- c("asv_id", "target", "pident")
h$vtx <- unname(id_to_vtx[h$target])
best <- h |> group_by(asv_id) |> filter(pident == max(pident)) |>
  summarise(pident = first(pident), all_vtx = paste(sort(unique(vtx)), collapse = ";"),
            vtx_vsearch = sort(unique(vtx))[1], ambiguous = n_distinct(vtx) > 1, .groups = "drop")
vtm <- map |> select(asv_id, total_reads) |> left_join(best, by = "asv_id") |>
  mutate(is_amf = !is.na(vtx_vsearch), vtx_kivlin = NA_character_) |>
  select(asv_id, total_reads, vtx_vsearch, pident, ambiguous, all_vtx, vtx_kivlin, is_amf)
write_csv(vtm, file.path(TAX_DIR, "maarjam_vsearch_mapping_final.csv"))
cat("AMF ASVs:", sum(vtm$is_amf), "| VT:", n_distinct(vtm$vtx_vsearch[vtm$is_amf]),
    "| AMF share of reads:", round(sum(vtm$total_reads[vtm$is_amf]) / sum(vtm$total_reads), 3), "\n")

## ---- phyloseq + analysis table ------------------------------------------------------------
sd <- sample_data(ps_legacy)[key$sample_name, ]
ps_new <- phyloseq(otu_table(st, taxa_are_rows = FALSE), sd)
saveRDS(ps_new, file.path(PHY_DIR, "ps_amf_all.rds")); saveRDS(ps_new, file.path(PHY_DIR, "ps_amf_euk_clean.rds"))
vt_of_seq <- setNames(vtm$vtx_vsearch, map$sequence)
env_cur <- read_csv(file.path(ENV_DIR, "depth_analysis_data.csv"), show_col_types = FALSE)
env_new <- update_env(env_cur, amf_summaries(st, vt_of_seq))[, names(env_cur)]
write_csv(env_new, file.path(ENV_DIR, "depth_analysis_data.csv"))
cat("01c complete.\n")

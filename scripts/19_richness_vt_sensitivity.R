## 19_richness_vt_sensitivity.R — sensitivity of richness to sequencing depth, and of
## the AM fungal results to tied VT assignments (Tables S13 and S14, Figure S2).
##
## A. Rarefaction. ITS counts are rarefied to the smallest ITS library and the
##    AM-fungal (VT-assigned) 18S reads to the smallest AM-fungal read total, 100
##    draws each; richness per sample is the mean over draws. ECM and saprotroph
##    richness use the priority guild rule. Depth contrasts are refitted with the
##    paired Bayesian model of Figure 5 (06), and the AM-fungal latitude
##    model of 16 is refitted. Rarefaction curves are saved as a figure.
## B. VT ties. Some AM fungal ASVs match two or more VT equally well; 01c assigns
##    the alphabetically first. Variants: (i) exclude tied ASVs; (ii) merge every set
##    of VT tied for any ASV into one taxon. For each variant (and the original):
##    VT richness depth contrast, within-stand depth test of composition (pH as
##    covariate), latitude test of composition (nodes permuted), paired overlap and
##    its same-node null.
## C. vsearch rerun with no accept limit (--maxaccepts 0 --maxrejects 0) to check
##    whether the 5-hit limit in 01c hid further ties. Needs vsearch (env VSEARCH or
##    PATH); skipped with a message if unavailable.
##
## Outputs: outputs/tables/richness_rarefaction_sensitivity.csv
##          outputs/tables/vt_tie_sensitivity.csv
##          outputs/tables/vsearch_full_search_check.csv
##          outputs/figures/rarefaction_curves.png

source("scripts/helpers.R")
source("scripts/bayes_helpers.R")
source("scripts/funguild_helpers.R")
suppressPackageStartupMessages({ library(tidyr); library(ggplot2); library(patchwork) })
nsamples <- phyloseq::nsamples
N_DRAWS <- as.integer(Sys.getenv("N_DRAWS", "100"))

env <- load_env() |> mutate(depth_num = ifelse(depth == "A", 0, 1), stand_id = paste0(node, stand))
env$lat_z <- standardize(env$latitude)

## ---- ITS counts and guilds -------------------------------------------------------
fg <- read.delim(file.path(TAX_DIR, "otu_table_for_FUNGuild.guilds.txt"), check.names = FALSE,
                 stringsAsFactors = FALSE)
ps_its <- load_ps("ITS"); key_its <- get_sample_key(ps_its)
samp_cols <- intersect(names(fg), key_its$sample_name)
stopifnot(length(samp_cols) == nrow(key_its))
its_counts <- t(as.matrix(fg[, samp_cols])); colnames(its_counts) <- fg$OTU_ID
its_cat <- guild_category(fg$Guild)
rownames(its_counts) <- key_its$sample_id_clean[match(rownames(its_counts), key_its$sample_name)]

## ---- AM fungal (VT-assigned) counts ---------------------------------------------
ps18 <- load_ps("18S"); key18 <- get_sample_key(ps18)
amf_otu <- as(otu_table(ps18), "matrix"); if (taxa_are_rows(ps18)) amf_otu <- t(amf_otu)
asv_map <- read_csv(file.path(TAX_DIR, "easton_18s_asv_map.csv"), show_col_types = FALSE)
vtx_map <- read_csv(file.path(TAX_DIR, "maarjam_vsearch_mapping_final.csv"), show_col_types = FALSE)
info <- left_join(asv_map, vtx_map[, c("asv_id", "is_amf", "vtx_vsearch", "ambiguous", "all_vtx")], by = "asv_id")
seqs <- colnames(amf_otu)
if (!length(intersect(seqs, info$sequence))) {          # same leader fallback as build_vt_matrix()
  leader <- "AGCAGCCGCGGTAATTCCAGCT"
  seqs <- ifelse(startsWith(seqs, leader), substr(seqs, nchar(leader) + 1, nchar(seqs)), seqs)
}
idx <- match(seqs, info$sequence)
vt <- info$vtx_vsearch[idx]; amb <- info$ambiguous[idx]; allv <- info$all_vtx[idx]
is_amf <- !is.na(vt) & vt != ""
amf_counts <- amf_otu[, is_amf, drop = FALSE]
vt <- vt[is_amf]; amb <- as.logical(amb[is_amf]); allv <- allv[is_amf]
rownames(amf_counts) <- key18$sample_id_clean[match(rownames(amf_counts), key18$sample_name)]
cat("AM fungal ASVs:", ncol(amf_counts), "| tied:", sum(amb), "| tied share of AM reads:",
    round(sum(amf_counts[, amb]) / sum(amf_counts), 3), "\n")

if (Sys.getenv("SKIP_A") != "1") {
## ---- A. rarefaction ---------------------------------------------------------------
set.seed(42)
its_min <- min(rowSums(its_counts)); amf_min <- min(rowSums(amf_counts))
cat("Rarefying ITS to", its_min, "reads and AM fungal reads to", amf_min, "\n")
rar <- replicate(N_DRAWS, {
  r_its <- vegan::rrarefy(its_counts, its_min)
  r_amf <- vegan::rrarefy(amf_counts, amf_min)
  cbind(ecm = rowSums(r_its[, its_cat == "ECM", drop = FALSE] > 0),
        sap = rowSums(r_its[, its_cat == "Saprotroph", drop = FALSE] > 0),
        vt  = apply(r_amf, 1, function(x) length(unique(vt[x > 0]))))
}, simplify = "array")
rar_mean <- apply(rar, c(1, 2), mean)
rd <- env |> select(sample_id_clean, depth, depth_num, node, stand_id, lat_z, richness_ECM,
                    richness_Saprotroph, vt_richness) |>
  inner_join(tibble(sample_id_clean = rownames(rar_mean), ecm_r = rar_mean[, "ecm"],
                    sap_r = rar_mean[, "sap"], vt_r = rar_mean[, "vt"]), by = "sample_id_clean")
stopifnot(nrow(rd) == 52)
amf_reads <- tibble(sample_id_clean = rownames(amf_counts), amf_reads = rowSums(amf_counts)) |>
  left_join(env |> select(sample_id_clean, depth), by = "sample_id_clean")
its_reads <- tibble(sample_id_clean = rownames(its_counts), its_reads = rowSums(its_counts),
                    ecm_reads = rowSums(its_counts[, its_cat == "ECM"])) |>
  left_join(env |> select(sample_id_clean, depth), by = "sample_id_clean")

fit_depth <- function(y, label, version) {
  r <- paired_depth_model(y, rd$depth_num, rd$stand_id, "log05")
  tibble(outcome = label, version = version, surface_mean = mean(y[rd$depth == "A"]),
         subsoil_mean = mean(y[rd$depth == "B"]), !!!as.list(r))
}
rows_a <- bind_rows(
  fit_depth(rd$richness_ECM, "ECM richness", "observed"),
  fit_depth(rd$ecm_r, "ECM richness", paste0("rarefied to ", its_min, " ITS reads")),
  fit_depth(rd$richness_Saprotroph, "Saprotroph richness", "observed"),
  fit_depth(rd$sap_r, "Saprotroph richness", paste0("rarefied to ", its_min, " ITS reads")),
  fit_depth(rd$vt_richness, "AM fungal richness", "observed"),
  fit_depth(rd$vt_r, "AM fungal richness", paste0("rarefied to ", amf_min, " AM fungal reads")))
## latitude model for AM fungal richness (as in 16_latitude_effects.R)
lat_fit <- function(y, version) {
  d <- data.frame(yz = standardize(transform_y(y, "log05")), depth_num = rd$depth_num,
                  lat_z = rd$lat_z, node = rd$node, stand = rd$stand_id)
  fit <- bayes_fit(yz ~ depth_num + lat_z + (1 | node) + (1 | stand), d)
  tibble(outcome = "AM fungal richness: latitude (std.)", version = version,
         !!!as.list(summ(as_draws_df(fit)$b_lat_z, "effect")), !!!as.list(fit_diagnostics(fit)))
}
rows_a <- bind_rows(rows_a, lat_fit(rd$vt_richness, "observed"),
                    lat_fit(rd$vt_r, paste0("rarefied to ", amf_min, " AM fungal reads")))
depth_reads <- bind_rows(
  amf_reads |> group_by(depth) |> summarise(value = mean(amf_reads), .groups = "drop") |> mutate(quantity = "mean AM fungal reads"),
  its_reads |> group_by(depth) |> summarise(value = mean(its_reads), .groups = "drop") |> mutate(quantity = "mean ITS reads"),
  its_reads |> group_by(depth) |> summarise(value = mean(ecm_reads), .groups = "drop") |> mutate(quantity = "mean ECM reads"))
print(depth_reads)
write_csv(rows_a |> mutate(across(where(is.numeric), ~ signif(.x, 4))),
          file.path(TABLES_DIR, "richness_rarefaction_sensitivity.csv"))
write_csv(depth_reads |> mutate(value = round(value)), file.path(TABLES_DIR, "reads_by_depth.csv"))

## rarefaction curves (ECM ASVs; AM fungal VT)
curve_df <- function(mat, label, step) {
  bind_rows(lapply(rownames(mat), function(s) {
    x <- mat[s, ]; x <- x[x > 0]; n <- sum(x)
    depths <- unique(c(seq(step, n, by = step), n))
    tibble(sample = s, reads = depths, richness = as.numeric(vegan::rarefy(x, depths)), panel = label)
  }))
}
ecm_mat <- its_counts[, its_cat == "ECM", drop = FALSE]
vt_mat <- t(rowsum(t(amf_counts), vt))
cd <- bind_rows(curve_df(ecm_mat, "ECM fungi (ASVs)", 250), curve_df(vt_mat, "AM fungi (VT)", 250)) |>
  mutate(depth = ifelse(substr(sample, nchar(sample), nchar(sample)) == "A", "Surface (0–20 cm)", "Subsoil (20–40 cm)"),
         depth = factor(depth, levels = c("Surface (0–20 cm)", "Subsoil (20–40 cm)")))  # surface first in key, as Fig 2
pc <- ggplot(cd, aes(reads, richness, group = sample, colour = depth)) +
  geom_line(alpha = 0.6, linewidth = 0.4) +
  facet_wrap(~panel, scales = "free") +
  scale_colour_manual(values = c("Surface (0–20 cm)" = "#8C510A", "Subsoil (20–40 cm)" = "#01665E"), name = NULL) +
  labs(x = "Reads sampled", y = "Expected richness") + theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(FIGS_DIR, "rarefaction_curves.png"), pc, width = 8, height = 3.8, dpi = 300)

}  # end part A

full_sets <- NULL
## ---- C. vsearch with no accept limit ----------------------------------------------
vs <- Sys.getenv("VSEARCH", "vsearch")
qfa <- file.path(DATA_DIR, "derived", "amf_asvs_trimmed.fa")
ref <- file.path(DATA_DIR, "reference", "maarjam_SSU_type_2021", "maarjam_database_SSU.fasta")
if (nzchar(Sys.which(vs)) || file.exists(vs)) {
  ## only the AM fungal ASVs matter here; searching all 18S ASVs exhaustively is slow
  fa <- readLines(qfa); hdr_i <- grep("^>", fa)
  ids <- sub("^>", "", fa[hdr_i]); keep_ids <- vtx_map$asv_id[vtx_map$is_amf]
  qsub <- tempfile(fileext = ".fa")
  writeLines(unlist(lapply(which(ids %in% keep_ids), function(k) fa[hdr_i[k]:(c(hdr_i, length(fa) + 1)[k + 1] - 1)])), qsub)
  qfa <- qsub
  hits_file <- tempfile(fileext = ".tsv")
  status <- system2(vs, c("--usearch_global", shQuote(qfa), "--db", shQuote(ref), "--id", "0.97",
                          "--blast6out", shQuote(hits_file), "--query_cov", "0.90", "--maxaccepts", "0",
                          "--maxrejects", "0", "--threads", "8"), stdout = FALSE, stderr = FALSE)
  stopifnot(status == 0)
  hdr <- grep("^>", readLines(ref), value = TRUE)
  id_to_vtx <- setNames(regmatches(hdr, regexpr("VTX[0-9]+", hdr)), sub("^>(\\S+).*", "\\1", hdr)[grepl("VTX[0-9]+", hdr)])
  h <- read.delim(hits_file, header = FALSE, stringsAsFactors = FALSE)[, 1:3]
  names(h) <- c("asv_id", "target", "pident"); h$vtx <- unname(id_to_vtx[h$target])
  full <- h |> group_by(asv_id) |> filter(pident == max(pident)) |>
    summarise(all_vtx_full = paste(sort(unique(vtx)), collapse = ";"), .groups = "drop")
  cmp <- vtx_map |> filter(is_amf) |> select(asv_id, all_vtx, total_reads) |> left_join(full, by = "asv_id")
  chk <- tibble(amf_asvs = nrow(cmp),
                tied_with_5_accepts = sum(grepl(";", cmp$all_vtx)),
                tied_with_no_limit = sum(grepl(";", cmp$all_vtx_full)),
                asvs_whose_best_vt_set_changed = sum(cmp$all_vtx != cmp$all_vtx_full, na.rm = TRUE),
                reads_in_changed_asvs_share = sum(cmp$total_reads[cmp$all_vtx != cmp$all_vtx_full], na.rm = TRUE) / sum(cmp$total_reads),
                first_vt_changed = sum(sub(";.*", "", cmp$all_vtx) != sub(";.*", "", cmp$all_vtx_full), na.rm = TRUE))
  print(chk)
  full_sets <- setNames(cmp$all_vtx_full, cmp$asv_id)
  write_csv(chk, file.path(TABLES_DIR, "vsearch_full_search_check.csv"))
} else cat("vsearch not found; part C skipped\n")

## ---- B. VT ties ------------------------------------------------------------------
merge_groups <- function(sets) {                       # union-find over tied VT sets
  parent <- list()
  find <- function(x) { while (!is.null(parent[[x]]) && parent[[x]] != x) x <- parent[[x]]; x }
  for (s in sets) {
    v <- strsplit(s, ";")[[1]]
    for (x in v) if (is.null(parent[[x]])) parent[[x]] <- x
    r <- find(v[1]); for (x in v[-1]) { rx <- find(x); if (rx != r) parent[[rx]] <- r }
  }
  vs <- names(parent); roots <- sapply(vs, find)
  lab <- tapply(vs, roots, function(z) paste(sort(z), collapse = "+"))
  setNames(lab[roots], vs)
}
grp <- merge_groups(unique(allv[amb]))
vt_merged <- ifelse(vt %in% names(grp), grp[vt], vt)
variants <- list(original = list(keep = rep(TRUE, length(vt)), lab = vt),
                 tied_ASVs_excluded = list(keep = !amb, lab = vt),
                 tied_VT_merged = list(keep = rep(TRUE, length(vt)), lab = vt_merged))
if (!is.null(full_sets)) {                              # exhaustive-search best matches (part C)
  asv_ids <- info$asv_id[idx][is_amf]
  fs <- unname(full_sets[asv_ids]); fs[is.na(fs)] <- vt[is.na(fs)]
  fgrp <- merge_groups(unique(fs[grepl(";", fs)]))
  first_full <- sub(";.*", "", fs)
  variants$full_search_ties_merged <- list(keep = rep(TRUE, length(vt)),
                                           lab = ifelse(first_full %in% names(fgrp), fgrp[first_full], first_full))
  variants$full_search_ties_excluded <- list(keep = !grepl(";", fs), lab = first_full)
}
md18 <- match_ps_to_env(ps18, env) |> mutate(stand_id = paste0(node, stand))
md18 <- md18[match(key18$sample_name[match(rownames(amf_counts), key18$sample_id_clean)], md18$sample_name), ]
rows_b <- bind_rows(lapply(names(variants), function(nm) {
  v <- variants[[nm]]
  mat <- t(rowsum(t(amf_counts[, v$keep, drop = FALSE]), v$lab[v$keep]))
  rownames(mat) <- md18$sample_name
  rich <- rowSums(mat > 0)
  r <- paired_depth_model(rich, md18$depth_num, md18$stand_id, "log05")
  d <- aitchison_distance_matrix(mat)
  dep <- test_within_stand(d, md18, "depth_factor", covars = "p_h")
  G <- gower_centre(d)
  lat <- unit_perm_test(G, align_md(md18, d), c("depth_factor", "latitude", "p_h"), "latitude", unit = "node")
  pa <- mat > 0
  sub <- md18 |> filter(depth == "B"); sur <- md18 |> filter(depth == "A")
  paired <- intersect(sub$stand_id, sur$stand_id)
  ov <- function(sb, sf) { s <- pa[sb, ]; if (!any(s)) NA_real_ else sum(s & pa[sf, ]) / sum(s) }
  obs <- sapply(paired, function(sid) ov(sub$sample_name[sub$stand_id == sid], sur$sample_name[sur$stand_id == sid]))
  nul <- sapply(paired, function(sid) {
    nd <- sub$node[sub$stand_id == sid]
    others <- sur$sample_name[sur$node == nd & sur$stand_id != sid]
    if (!length(others)) NA_real_ else mean(sapply(others, ov, sb = sub$sample_name[sub$stand_id == sid]))
  })
  keep <- !is.na(nul)
  tibble(variant = nm, n_taxa = ncol(mat), share_of_am_reads = sum(mat) / sum(amf_counts),
         surface_mean_richness = mean(rich[md18$depth == "A"]), subsoil_mean_richness = mean(rich[md18$depth == "B"]),
         richness_contrast_mean = r[["contrast_mean"]], richness_contrast_lo = r[["contrast_hdi_lo"]],
         richness_contrast_hi = r[["contrast_hdi_hi"]], richness_pr = r[["effect_pr_direction"]],
         depth_F = dep[["F"]], depth_p = dep[["p"]], latitude_F = lat[["F"]], latitude_p = lat[["p"]],
         overlap_same_stand = mean(obs[keep]), overlap_same_node_null = mean(nul[keep]), n_stands_overlap = sum(keep))
}))
print(rows_b)
write_csv(rows_b |> mutate(across(where(is.numeric), ~ signif(.x, 4))), file.path(TABLES_DIR, "vt_tie_sensitivity.csv"))

cat("19_richness_vt_sensitivity.R complete.\n")

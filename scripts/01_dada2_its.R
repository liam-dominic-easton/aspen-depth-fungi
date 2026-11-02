## 01_dada2_its.R — ITS amplicon processing from raw FASTQs.
##
## Inputs (not shipped; set paths with environment variables):
##   ITS_READS    folder with the 53 ITS read pairs (*-fungal_*_R1/R2_001.fastq.gz;
##                deposited in [SRA accession to be added])
##   CUTADAPT     path to cutadapt (v5.2 used; see README)
##   UNITE_FASTA  UNITE general FASTA release for Fungi, v19.02.2025, dynamic file
##                sh_general_release_dynamic_19.02.2025.fasta
##                (doi:10.15156/BIO/3301229; SHA-256 of the .tgz archive:
##                a340864f947c517b671f1054a261e613a6619ef88dd8000de1699f5e1f21413b)
## Outputs: data/derived/its_seqtab_nochim.rds, its_taxonomy.rds, its_dada2_track.csv
## Settings: DADA2 1.34.0; reads with Ns removed, then primers removed with
## cutadapt following the DADA2 ITS workflow (the degenerate 5.8S-Fun/ITS4-Fun
## primers, if left on the reads, split identical sequences into several ASVs);
## no truncation; maxEE = c(2, 3); truncQ = 2; rm.phix; minLen 50;
## error models learned from filtered files up to DADA2's default 1e8-base limit
## (19 files per direction for this dataset); per-sample inference; merging
## with minOverlap = 12; consensus chimera removal; assignTaxonomy (tryRC = TRUE,
## seed 100) against UNITE.
## Run time on 8 cores: ~20 min (DADA2) + ~60-80 min (taxonomy).

suppressPackageStartupMessages(library(dada2))
source("scripts/helpers.R")
READS <- Sys.getenv("ITS_READS"); UNITE <- Sys.getenv("UNITE_FASTA")
CUTADAPT <- Sys.getenv("CUTADAPT", "cutadapt")
SKIP_TAX <- Sys.getenv("SKIP_TAXONOMY") == "1"   # verification runs only
stopifnot(dir.exists(READS), SKIP_TAX || file.exists(UNITE))
OUT <- Sys.getenv("DERIVED_DIR", file.path(DATA_DIR, "derived")); dir.create(file.path(OUT, "filtered_its"), recursive = TRUE, showWarnings = FALSE)

fnFs <- sort(list.files(READS, pattern = "fungal.*_R1_001\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(READS, pattern = "fungal.*_R2_001\\.fastq\\.gz$", full.names = TRUE))
sn <- sub("_R1_001\\.fastq\\.gz$", "", basename(fnFs))
stopifnot(length(fnFs) == length(fnRs), identical(sn, sub("_R2_001\\.fastq\\.gz$", "", basename(fnRs))))

## ---- primer removal ---------------------------------------------------------
## Reads are in reverse orientation: R1 begins with ITS4-Fun and R2 with 5.8S-Fun
## (Taylor et al. 2016). Each read is trimmed of its 5' primer and, where the read
## runs through the amplicon, of the reverse complement of the other primer.
FWD <- "AACTTTYRRCAAYGGATCWCT"        # 5.8S-Fun
REV <- "AGCCTCCGCTTATTGATATGCTTAART"  # ITS4-Fun
rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
dir.create(file.path(OUT, "filtN_its"), showWarnings = FALSE)
dir.create(file.path(OUT, "cutadapt_its"), showWarnings = FALSE)
nFs <- file.path(OUT, "filtN_its", basename(fnFs)); nRs <- file.path(OUT, "filtN_its", basename(fnRs))
cutFs <- file.path(OUT, "cutadapt_its", basename(fnFs)); cutRs <- file.path(OUT, "cutadapt_its", basename(fnRs))
ftN <- filterAndTrim(fnFs, nFs, fnRs, nRs, maxN = 0, multithread = TRUE)
for (i in seq_along(fnFs)) {
  status <- system2(CUTADAPT, c("-g", REV, "-a", rc(FWD), "-G", FWD, "-A", rc(REV), "-n", 2,
                                "--cores", 0, "--quiet",
                                "--json", shQuote(sub("\\.fastq\\.gz$", ".cutadapt.json", cutFs[i])),
                                "-o", shQuote(cutFs[i]), "-p", shQuote(cutRs[i]),
                                shQuote(nFs[i]), shQuote(nRs[i])))
  if (!identical(status, 0L)) stop("cutadapt failed for ", sn[i])
}
## residual primer check on the first three samples (should be ~0)
count_hits <- function(primer, fn) {
  seqs <- ShortRead::sread(ShortRead::readFastq(fn))
  sum(Biostrings::vcountPattern(primer, seqs, fixed = FALSE) > 0)
}
for (i in 1:3) cat(sn[i], "residual 5' primers: R1", count_hits(REV, cutFs[i]),
                   "R2", count_hits(FWD, cutRs[i]), "\n")

filtFs <- file.path(OUT, "filtered_its", basename(fnFs)); filtRs <- file.path(OUT, "filtered_its", basename(fnRs))
ft <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen = c(0, 0), maxN = 0,
                    maxEE = c(2, 3), truncQ = 2, rm.phix = TRUE, minLen = 50, compress = TRUE,
                    multithread = TRUE)
ft[, "reads.in"] <- ftN[, "reads.in"]  # track from the raw read count
keep <- file.exists(filtFs) & file.exists(filtRs)
errF <- learnErrors(filtFs[keep], multithread = TRUE)
errR <- learnErrors(filtRs[keep], multithread = TRUE)
derepFs <- derepFastq(filtFs[keep]); names(derepFs) <- sn[keep]
derepRs <- derepFastq(filtRs[keep]); names(derepRs) <- sn[keep]
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 12)
st <- removeBimeraDenovo(makeSequenceTable(mergers), method = "consensus", multithread = TRUE)
saveRDS(st, file.path(OUT, "its_seqtab_nochim.rds"))

getN <- function(x) sum(getUniques(x))
write.csv(data.frame(sample = sn[keep], ft[keep, , drop = FALSE], denoisedF = sapply(dadaFs, getN),
                     denoisedR = sapply(dadaRs, getN), merged = sapply(mergers, getN),
                     nonchim = rowSums(st)),
          file.path(OUT, "its_dada2_track.csv"), row.names = FALSE)

if (SKIP_TAX) quit(save = "no")
set.seed(100)
tax <- assignTaxonomy(st, UNITE, multithread = TRUE, tryRC = TRUE, outputBootstraps = TRUE)
saveRDS(tax, file.path(OUT, "its_taxonomy.rds"))
cat("01_dada2_its.R complete:", nrow(st), "samples x", ncol(st), "ASVs\n")

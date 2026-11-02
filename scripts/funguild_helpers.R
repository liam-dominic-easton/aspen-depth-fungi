## funguild_helpers.R — FUNGuild annotation in R.
##
## Re-implements the matching rule of FUNGuild's Guilds_v1.1.py (Nguyen et al.
## 2016): an ASV's taxonomy is split into rank names, each database taxon is
## matched against those names, and the most specific match (highest
## taxonomicLevel) is kept. Species epithets from UNITE ("s__roseus") do not
## form "Genus_species" names, so matches are at genus level or above.
## Database: data/reference/funguild_db_<date>.json, a snapshot of
## http://www.stbates.org/funguild_db_2.php (downloaded 2026-09-26).

read_funguild_db <- function(path) {
  db <- jsonlite::fromJSON(path)
  db$taxon_key <- gsub(" ", "_", trimws(db$taxon))
  db$taxonomicLevel <- as.integer(db$taxonomicLevel)
  db[order(-db$taxonomicLevel), ]
}

## tax: character matrix/data.frame of rank columns (values like "g__Ilyonectria")
funguild_assign <- function(tax, db) {
  tax <- as.data.frame(tax, stringsAsFactors = FALSE)
  tokens <- apply(tax, 1, function(r) {
    r <- r[!is.na(r)]
    unique(sub("^([a-z]__)+", "", r))
  })
  if (!is.list(tokens)) tokens <- split(tokens, col(tokens))
  lut <- split(seq_len(nrow(db)), db$taxon_key)
  out <- lapply(tokens, function(tk) {
    hits <- unlist(lut[intersect(tk, names(lut))], use.names = FALSE)
    if (!length(hits)) return(NA_integer_)
    ## highest taxonomicLevel wins; ties (e.g. the snapshot gives the family
    ## Plectosphaerellaceae level 13, the same as its genera) go to the more
    ## specific rank in the ASV's taxonomy, as in the official script
    pos <- match(db$taxon_key[hits], tk)
    hits[order(-db$taxonomicLevel[hits], -pos)][1]
  })
  idx <- unlist(out)
  g <- function(col) ifelse(is.na(idx), "-", db[[col]][idx])
  data.frame(Taxon = g("taxon"), `Taxon Level` = ifelse(is.na(idx), "-", as.character(db$taxonomicLevel[idx])),
             `Trophic Mode` = g("trophicMode"), Guild = g("guild"), `Growth Morphology` = g("growthForm"),
             Trait = g("trait"), `Confidence Ranking` = g("confidenceRanking"), Notes = g("notes"),
             `Citation/Source` = g("citationSource"), check.names = FALSE, stringsAsFactors = FALSE)
}

## priority rule for ASVs with one or more guilds (unassigned first, then
## ECM > AMF > saprotroph > pathogen > endophyte > other)
guild_category <- function(guild) {
  g <- stringr::str_squish(ifelse(is.na(guild), "", guild)); gl <- tolower(g)
  unassigned <- g %in% c("", "-") | grepl("^\\s*(unassigned|unknown|undefined|none)\\s*$", gl)
  dplyr::case_when(unassigned ~ "Unassigned",
                   grepl("ectomycorrhizal", gl) ~ "ECM",
                   grepl("arbuscular", gl) ~ "AMF",
                   grepl("saprotroph", gl) ~ "Saprotroph",
                   ## any FUNGuild pathogen guild; 5 ASVs (26 reads) carry only animal,
                   ## insect or undefined pathogen labels
                   grepl("pathogen", gl) ~ "Plant pathogen",
                   grepl("endophyte", gl) ~ "Endophyte",
                   nzchar(g) ~ "Other",
                   TRUE ~ "Unassigned")
}

# Vertical zonation of soil fungi beneath trembling aspen

Analysis code for:

> Easton L, Kivlin SN, Phillips RP, Bennett J, Karst J. *Vertical zonation of soil fungi in trembling aspen forests persists across latitudes.* [Journal, year, DOI to be added]

## The study

We sampled mineral soil at 0–20 cm and 20–40 cm in 27 trembling aspen (*Populus tremuloides*) stands, grouped into 9 nodes across the parkland–boreal transition in western Canada (52 samples). Fungi were characterised by ITS2 amplicon sequencing, with guilds assigned through UNITE and FUNGuild. Arbuscular mycorrhizal (AM) fungi were characterised by 18S (NS31/AML2) sequencing and assigned to MaarjAM virtual taxa (VT) with vsearch. The scripts test how guild composition differs between depths and along the latitudinal gradient. They use CLR/Aitchison composition analyses, dbRDA and PERMANOVA with design-based permutations, Bayesian mixed models (brms), indicator species and nestedness.

**The data are not included in this repository** (see [Data](#data)).

## Repository layout

```
scripts/           R scripts; run them from the repository root
data/              input data (not included; create the folders listed under Data)
outputs/tables/    CSV results (created and written by the scripts)
outputs/figures/   PNG figures (created and written by the scripts)
R_packages.txt     R packages and the versions used
```

## Software

| Component | Version used | Used by |
|---|---|---|
| R | 4.4.3 | all scripts |
| DADA2 (Bioconductor) | 1.34.0 | `01_dada2_its.R`, `01c_dada2_amf_vt.R` |
| cutadapt | 5.2 | ITS primer removal (`01_dada2_its.R`) |
| vsearch | 2.28.1 | VT assignment (`01c_dada2_amf_vt.R`; part C of `19_richness_vt_sensitivity.R`) |
| phyloseq, vegan, brms, rstan, indicspecies, … | see `R_packages.txt` | analyses |

Install the R packages from CRAN and Bioconductor:

```r
install.packages(c("vegan", "permute", "indicspecies", "ape", "brms", "rstan", "loo",
                   "bayestestR", "dplyr", "tidyr", "readr", "tibble", "stringr", "glue",
                   "ggplot2", "patchwork", "ggrepel", "maps", "scales", "readxl", "jsonlite"))
install.packages("BiocManager")
BiocManager::install(c("phyloseq", "dada2", "Biostrings", "ShortRead"))
```

You can also use conda (conda-forge and bioconda channels). This is an example; the environment name is arbitrary:

```bash
conda create -n aspen-depth-fungi -c conda-forge -c bioconda \
  r-base=4.4 bioconductor-phyloseq=1.50.0 r-vegan=2.7_2 r-brms=2.23.0 r-rstan=2.32.7 \
  r-permute r-indicspecies r-ape r-loo r-bayestestr r-dplyr r-tidyr r-readr r-tibble \
  r-stringr r-glue r-ggplot2 r-patchwork r-ggrepel r-maps r-scales r-readxl r-jsonlite \
  vsearch=2.28.1 cxx-compiler make
```

Notes:

- **RStan from conda-forge** ships oneTBB. Before running the Bayesian scripts, point R at the supplied compiler flags:
  `export R_MAKEVARS_USER="$PWD/scripts/Makevars.stan"`. CRAN installations normally do not need this.
- **DADA2 on macOS (Apple Silicon) with conda:** the bioconda binary can link an incompatible TBB. Installing DADA2 1.34.0 from the Bioconductor 3.20 source package avoids this.
- **cutadapt** can conflict with other Python packages in the same environment. If so, install it in its own environment and pass its path with `CUTADAPT=`.

## Data

The sequence data, processed tables and environmental data are archived separately:

- Raw reads: NCBI SRA [accession to be added].
- Processed analysis inputs (the files in `data/` below): [repository and DOI to be added].

Place the files as follows (paths are relative to the repository root):

| Location | Contents | Source |
|---|---|---|
| any folder; set `ITS_READS` | 53 ITS read pairs (52 samples and a PCR blank), `*-fungal_*_R1_001.fastq.gz` / `*_R2_001.fastq.gz` | NCBI SRA [accession to be added] |
| any folder; set `AMF_READS` | 53 18S read pairs, `*-amf_*_R1_001.fastq.gz` (forward reads are used) | NCBI SRA [accession to be added] |
| any file; set `UNITE_FASTA` | UNITE general FASTA release for Fungi, version 19.02.2025, dynamic file `sh_general_release_dynamic_19.02.2025.fasta`. SHA-256 of the `.tgz` archive: `a340864f947c517b671f1054a261e613a6619ef88dd8000de1699f5e1f21413b` | Abarenkov et al. 2025, [doi:10.15156/BIO/3301229](https://doi.org/10.15156/BIO/3301229) |
| `data/reference/maarjam_SSU_type_2021/maarjam_database_SSU.fasta` | MaarjAM SSU database, 2021 release; VTX codes are read from the FASTA headers | MaarjAM ([maarjam.ut.ee](https://maarjam.ut.ee); Öpik et al. 2010) |
| `data/reference/natural_earth/ne_50m_admin_1_states_provinces_lines.{shp,shx,dbf,prj}` | Natural Earth 1:50m admin-1 boundary lines, version 5.1.0 (public domain); used only for the provincial boundaries in Figure 1A. Without them, `09_figures.R` skips Figure 1 with a warning | [naturalearthdata.com](https://www.naturalearthdata.com/downloads/50m-cultural-vectors/) (Admin 1 – States, Provinces: boundary lines) |
| `data/reference/funguild_db_2026-09-26.json` | FUNGuild database snapshot (JSON) downloaded on 26 September 2026 | `http://www.stbates.org/funguild_db_2.php` (Nguyen et al. 2016) |
| `data/environmental/depth_analysis_data.csv` | One row per sample (52): site, climate, soil, fine-root, ECM colonization and understory variables, plus the ITS and 18S summary columns that `01b` and `01c` recompute | Processed data [DOI to be added]. Site, soil, root, colonization and understory variables are from Karst et al. (2021; data on Borealis, [DOI to be added]). Climate (2008–2018 means) is from the Canadian Forest Service CIPHA network |
| `data/raw_environmental/Soil texture.xls` | Laser-diffraction particle-size output, sheet `Data Table` (one column per sample; rows `Sand`, `Silt`, `Clay`); read by `01d` | Processed data [DOI to be added] |
| `data/raw_environmental/Plant understorey.xlsx` | Percent cover of each understory taxon in two 1-m² quadrats per stand, sheet `raw data`; read by `01d` | Processed data [DOI to be added]; Karst et al. (2021) |
| `data/reference/meng2023_5_11770.csv` | Plant mycorrhizal types for 11,770 species (file `5_11770.csv`, CC0); read by `01d` | Meng et al. 2023, Dryad [doi:10.5061/dryad.x0k6djhq2](https://doi.org/10.5061/dryad.x0k6djhq2) |
| `data/reference/understory_taxon_names.csv` | Name under which each understory taxon is looked up in the Meng et al. table (species, or genus/family for taxa not identified to species or absent from it) | This repository |
| `data/raw_environmental/Master_with pH added.xlsx` | Stand × depth table (sheet `Sheet1`) of site, climate, fine-root and soil-C variables; used by `16_latitude_effects.R` | Processed data [DOI to be added] |
| `data/phyloseq/` | `ps_fungi_all.rds`, `ps_fungi_clean.rds` (written by `01b`); `ps_amf_all.rds`, `ps_amf_euk_clean.rds` (written by `01c`) | Processed data, or rebuild with the `01` scripts |
| `data/taxonomy/` | `otu_table_for_FUNGuild.guilds.txt` (`01b`); `easton_18s_asv_map.csv`, `maarjam_vsearch_mapping_final.csv` (`01c`) | Processed data, or rebuild |
| `data/derived/` | `its_seqtab_nochim.rds`, `its_taxonomy.rds`, `its_dada2_track.csv` (`01`); `amf_seqtab_nochim.rds`, `amf_dada2_track.csv`, `amf_asvs_trimmed.fa`, `vsearch_maarjam_hits.txt` (`01c`) | Processed data, or rebuild |
| `data/legacy_inputs/` | Starting copies of the phyloseq objects, taxonomy tables and analysis table. `01b`/`01c` take the sample metadata from them | Created automatically; see below |

**`data/legacy_inputs/`.** `01b` and `01c` take the sample metadata (phyloseq `sample_data`) and the starting analysis table from this folder. If the folder does not exist, `01b` creates it. It then copies four files there before overwriting them: `ps_fungi_all.rds`, `ps_fungi_clean.rds`, `otu_table_for_FUNGuild.guilds.txt` and `depth_analysis_data.csv`. `01c` adds the 18S files: `ps_amf_all.rds`, `ps_amf_euk_clean.rds`, `easton_18s_asv_map.csv` and `maarjam_vsearch_mapping_final.csv`. To rebuild from the raw reads:

1. Place the distributed `data/` files.
2. Run `01` → `01b` → `01c` in that order.

Do not create an empty `data/legacy_inputs/` folder by hand: `01b` copies the starting files only when the folder is absent.

The analysis scripts (`02`–`90`) only need these inputs:

- `data/phyloseq/`, `data/taxonomy/` and `data/environmental/`;
- the MaarjAM FASTA (for `04`);
- `data/raw_environmental/` (for `01d` and `16`);
- `data/reference/meng2023_5_11770.csv` and `data/reference/understory_taxon_names.csv` (for `01d`);
- `data/derived/amf_asvs_trimmed.fa` (only for part C of `19`, which runs when vsearch is available).

## How to run

Run everything from the repository root.

**1. Amplicon processing (optional; slow).** Skip this step if you use the distributed processed tables.

```bash
ITS_READS=/path/to/ITS_fastq \
UNITE_FASTA=/path/to/sh_general_release_dynamic_19.02.2025.fasta \
CUTADAPT=/path/to/cutadapt \
  Rscript --vanilla scripts/01_dada2_its.R          # ~20 min DADA2 + 60-80 min taxonomy (8 cores)

Rscript --vanilla scripts/01b_build_its_tables.R    # FUNGuild, ITS phyloseq, ITS summary columns

AMF_READS=/path/to/18S_fastq \
  Rscript --vanilla scripts/01c_dada2_amf_vt.R      # DADA2 + vsearch against MaarjAM
```

**2. All analyses, tables and figures.**

```bash
Rscript --vanilla scripts/00_run_all_analyses.R
```

This runs the scripts in dependency order and writes `outputs/session_info.txt`. A full run takes about 30–60 min. The slowest steps are:

- the Bayesian models (`06`, `07`, `10`, `12`, `13`, `14`, `16`, `17`, `19`);
- the depth indicator tests with 99,999 permutations (`03`, `03b`, `04`).

You can run a single script with `Rscript --vanilla scripts/<script>.R` once its inputs exist:

- `10` must run before `06`, `12`, `14`, `16`, `17` and `18`;
- `02` must run before `08`;
- `03` must run before `15`;
- `06` and `07` must run before `09`.

**Settings (environment variables).**

| Variable | Default | Effect |
|---|---|---|
| `NPERM` | `9999` | permutations for permutation tests (`NPERM=999` for a quick check) |
| `NPERM_IND` | `99999` | permutations for depth indicator species (`NPERM_IND=9999` for a quick check) |
| `EXCLUDE_SAMPLES` | none | comma-separated sample IDs to drop everywhere, e.g. `MOO3B` |
| `OUTPUT_TAG` | none | write to `outputs/tables_<tag>/` and `outputs/figures_<tag>/` |
| `VSEARCH` | `vsearch` | path to vsearch (`01c`, `19`) |
| `CUTADAPT` | `cutadapt` | path to cutadapt (`01`) |
| `N_DRAWS` | `100` | rarefaction draws (`19`) |
| `SKIP_TAXONOMY=1` | off | `01`: stop after the sequence table |
| `SKIP_A=1` | off | `19`: skip the rarefaction part |
| `CHECK_MODE=1` | off | `01b`, `01c`: only check that the summary columns computed from the tables in `data/legacy_inputs/` reproduce the legacy analysis table |
| `DERIVED_DIR` | `data/derived` | output folder of `01` |

**Reproducibility notes.**

- Given the raw reads, `01` and `01c` produce the same sequence tables on every run.
- DADA2's `assignTaxonomy()` breaks ties with a random number generator that `set.seed()` does not control. Re-running the taxonomy step can therefore change a few ITS labels. The published results use the taxonomy file distributed as `data/derived/its_taxonomy.rds`.
- Permutation tests, rarefaction and the Bayesian models use fixed seeds (42).

## Scripts and manuscript items

Table and figure numbers refer to the manuscript and its supplementary material. Outputs are in `outputs/tables/` (CSV) and `outputs/figures/` (PNG).

| Item | Script | Output |
|---|---|---|
| Table 1 | `02_permanova_depth_environment.R` | `table1_permanova.csv` |
| Table 2 | `03b_indicators_its_genus.R` | rows of `indicator_genera_by_depth.csv` with BH-adjusted p < 0.05 |
| Table 3 | `12_guild_balance_roots_varpart.R` | `table2_varpart.csv` (fractions), `varpart_four_set_unique_tests.csv` (P values) |
| Figure 1 | `09_figures.R` | `study_area_climate.png` |
| Figure 2 | `18_latitude_figure.R` | `latitude_by_depth.png` |
| Figure 3 | `09_figures.R` | `pcoa_with_envfit.png` |
| Figure 4 | `20_paired_stand_figure.R` | `paired_depth_stands.png` |
| Figure 5 | `06_bayesian_depth_comparisons.R` → `09_figures.R` | `bayesian_depth_comparisons.csv` → `depth_comparison_forest.png` |
| Figure 6 | `07_root_biomass_bayesian.R` → `09_figures.R` | `root_biomass_bayesian_models.csv` → `root_biomass_model_comparison.png` |
| Figure S1 | `12_guild_balance_roots_varpart.R` | `variance_partition_venn.png` |
| Figure S2 | `19_richness_vt_sensitivity.R` | `rarefaction_curves.png` |
| Table S1 | `17_latitude_som_test.R` | `latitude_som_models.csv` |
| Table S2 | `13_additional_sensitivities.R` | `colonization_vs_soil_dna.csv` |
| Table S3 | `02_permanova_depth_environment.R`, `14_basal_area_sensitivity.R` | `depth_x_env_paired_shift.csv`, `depth_x_mat_given_biome.csv`, `depth_x_temperature_metrics.csv`, `basal_area_community_sensitivity.csv` |
| Table S4 | `16_latitude_effects.R` | `latitude_bayes_models.csv` |
| Table S5 | `16_latitude_effects.R`, `13_additional_sensitivities.R` | `latitude_premise_checks.csv`, `understory_vs_latitude.csv` |
| Table S6 | `05_nestedness_betadisper.R` | `nestedness_results.csv` |
| Table S7 | `12_guild_balance_roots_varpart.R` | `sap_ecm_paired_predictor_tests.csv` |
| Table S8 | `07_root_biomass_bayesian.R`, `12_guild_balance_roots_varpart.R` | `root_biomass_bayesian_models.csv`, `sap_ecm_root_som_models.csv` |
| Table S9 | `12_guild_balance_roots_varpart.R` | `varpart_four_set_all_fractions.csv` |
| Table S10 | `10_guild_ambiguity_sensitivity.R` | `guild_ambiguity_depth_inference.csv` |
| Table S11 | `13_additional_sensitivities.R` | `guild_confidence_sensitivity.csv`, `guild_confidence_ratio_models.csv` |
| Table S12 | `13_additional_sensitivities.R` | `guild_ratio_depth_models.csv` |
| Table S13 | `19_richness_vt_sensitivity.R` | `richness_rarefaction_sensitivity.csv` |
| Table S14 | `19_richness_vt_sensitivity.R` | `vt_tie_sensitivity.csv` |
| Table S15 | `15_pathogen_taxon_audit.R` | `pathogen_genera_depth_audit.csv` |
| Numbers quoted in the text | `90_qc_checks.R` | `qc_reported_numbers.csv` |

The scripts write some further outputs that do not appear as separate items, such as `root_biomass_loo.csv`, `guild_component_root_som_models.csv`, `vsearch_full_search_check.csv`, the `*_full.csv` indicator tables, the AM fungal VT indicator tables (`indicator_vtx_by_depth.csv`, `indicator_vtx_by_biome.csv`), the dispersion tests (`betadisper_by_depth.csv` and `.png`), the per-stand VT overlap (`vt_depth_overlap.csv`) and `within_stand_mat.png` from `08_mat_robustness.R`. Results from some of these are quoted in the text. The shared functions are in `helpers.R` (data loading, distances, design-based permutation tests), `bayes_helpers.R` (brms settings and summaries) and `funguild_helpers.R` (FUNGuild matching and guild categories).

## Citation

If you use this code, please cite the article:

> Easton L, Kivlin SN, Phillips RP, Bennett J, Karst J. [Year]. Vertical zonation of soil fungi in trembling aspen forests persists across latitudes. [Journal], [volume], [pages]. [DOI to be added]

and the code archive: [Zenodo DOI to be added].

References for the data sources and databases:

- Abarenkov K, Zirk A, Piirmann T, Pöhönen R, Ivanov F, Nilsson RH, Kõljalg U. 2025. UNITE general FASTA release for Fungi. Version 19.02.2025. UNITE Community. https://doi.org/10.15156/BIO/3301229
- Karst J, Franklin J, Simeon A, Light A, Bennett JA, Erbilgin N. 2021. Assessing the dual-mycorrhizal status of a widespread tree species as a model for studies on stand biogeochemistry. *Mycorrhiza* 31: 313–324. https://doi.org/10.1007/s00572-021-01029-2
- Nguyen NH, Song Z, Bates ST, Branco S, Tedersoo L, Menke J, Schilling JS, Kennedy PG. 2016. FUNGuild: an open annotation tool for parsing fungal community datasets by ecological guild. *Fungal Ecology* 20: 241–248.
- Öpik M, Vanatoa A, Vanatoa E, Moora M, Davison J, Kalwij JM, Reier Ü, Zobel M. 2010. The online database MaarjAM reveals global and ecosystemic distribution patterns in arbuscular mycorrhizal fungi (Glomeromycota). *New Phytologist* 188: 223–241.

## Licence

The code is released under the MIT Licence (see `LICENSE`). The licence covers the code only. The data are not part of this repository, and their own terms of use apply.

## Run all analyses, tables and figures, in order, from the processed inputs
## in data/ (see README.md). Run from the repository root:
##   Rscript --vanilla scripts/00_run_all_analyses.R
##
## Bayesian scripts need a working brms/RStan toolchain. If RStan comes from
## conda-forge, set R_MAKEVARS_USER to scripts/Makevars.stan (see README.md).
## Full run: ~30-60 min (NPERM = 9999 permutations, NPERM_IND = 99999 for
## indicator species). Set NPERM / NPERM_IND lower for a quick check.
## Sensitivity runs: EXCLUDE_SAMPLES=MOO3B OUTPUT_TAG=no_MOO3B.

source("scripts/helpers.R")

active_scripts <- c(
  "scripts/01d_texture_understory.R",          # texture and understory columns from raw files
  "scripts/10_guild_ambiguity_sensitivity.R",  # guild tables used downstream; Table S10
  "scripts/02_permanova_depth_environment.R",  # Table 1; Table S3
  "scripts/03_indicators_its.R",               # ASV-level indicators (input to 15)
  "scripts/03b_indicators_its_genus.R",        # Table 2
  "scripts/04_indicators_amf_vt.R",            # AM fungal VT indicators (quoted in the text)
  "scripts/05_nestedness_betadisper.R",        # Table S6; dispersion tests quoted in the text
  "scripts/06_bayesian_depth_comparisons.R",   # Figure 5 (model estimates)
  "scripts/07_root_biomass_bayesian.R",        # Figure 6 (model estimates); Table S8
  "scripts/12_guild_balance_roots_varpart.R",  # Table 3, Figure S1; Tables S7, S8, S9
  "scripts/13_additional_sensitivities.R",     # Tables S2, S11, S12; understory rows of S5
  "scripts/14_basal_area_sensitivity.R",       # Table S3
  "scripts/15_pathogen_taxon_audit.R",         # Table S15
  "scripts/16_latitude_effects.R",             # latitude models (Tables S4, S5)
  "scripts/17_latitude_som_test.R",            # Table S1: latitude vs organic matter / sand
  "scripts/08_mat_robustness.R",               # within-stand distance vs MAT (figure not in the manuscript)
  "scripts/09_figures.R",                      # Figures 1, 3, 5, 6
  "scripts/18_latitude_figure.R",              # Figure 2
  "scripts/20_paired_stand_figure.R",          # Figure 4
  "scripts/19_richness_vt_sensitivity.R",      # Tables S13, S14; Figure S2 (VSEARCH=path for part C)
  "scripts/90_qc_checks.R"                     # numbers quoted in the text
)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
for (script in active_scripts) {
  cat("\n", strrep("=", 78), "\nRunning ", script, "\n", sep = "")
  status <- system2(rscript, c("--vanilla", "-e",
                               shQuote(sprintf("options(mc.cores = 4); source('%s')", script))))
  if (!identical(status, 0L)) stop("Workflow failed in ", script, " (status ", status, ")")
}

sink(file.path(TABLES_DIR, "..", paste0("session_info", .tag, ".txt")))
print(sessionInfo())
sink()
cat("\nActive workflow complete.\n")

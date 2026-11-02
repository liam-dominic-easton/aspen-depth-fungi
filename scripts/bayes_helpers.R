## bayes_helpers.R — shared Bayesian model settings.
##
## All univariate inference in the package uses brms/RStan with:
##   * the response standardized (mean 0, SD 1) after its transform, so the
##     Normal(0, 1) priors below are weakly informative for every outcome
##     (on a raw scale they would dominate outcomes measured in large units,
##     such as ECM tips per metre of root);
##   * Normal(0, 1) on fixed effects and the intercept; half-Normal(0, 1) on
##     group-level SDs and the residual SD;
##   * 4 chains x 4,000 iterations (1,000 warmup), adapt_delta = 0.99, seed 42;
##   * 90% highest-density intervals (bayestestR::hdi) and Pr(direction), the
##     posterior probability that an effect has its estimated sign.
## Compiled models are cached per formula and reused with update(newdata = ).
##
## RStan from conda-forge needs R_MAKEVARS_USER pointing at Makevars.stan (see README).

suppressPackageStartupMessages({
  library(brms)
  library(bayestestR)
})
nsamples <- phyloseq::nsamples   # brms masks phyloseq::nsamples

BAYES_PRIORS_RE <- c(prior(normal(0, 1), class = "b"),
                     prior(normal(0, 1), class = "Intercept"),
                     prior(normal(0, 1), class = "sd"),
                     prior(normal(0, 1), class = "sigma"))
BAYES_PRIORS_FIXED <- c(prior(normal(0, 1), class = "b"),
                        prior(normal(0, 1), class = "Intercept"),
                        prior(normal(0, 1), class = "sigma"))

.bayes_cache <- new.env()

## fit (or refit a cached compiled model) for formula `f` on `data`
bayes_fit <- function(f, data, seed = 42) {
  key <- paste(deparse(f), collapse = " ")
  has_re <- grepl("|", key, fixed = TRUE)
  cores <- getOption("mc.cores", 4)
  if (is.null(.bayes_cache[[key]])) {
    fit <- brm(f, data = data, prior = if (has_re) BAYES_PRIORS_RE else BAYES_PRIORS_FIXED,
               chains = 4, iter = 4000, warmup = 1000, seed = seed, cores = cores,
               backend = "rstan", silent = 2, refresh = 0,
               control = list(adapt_delta = 0.99))
    assign(key, fit, envir = .bayes_cache)
  } else {
    fit <- suppressMessages(update(.bayes_cache[[key]], newdata = data, seed = seed,
                                   cores = cores, refresh = 0, silent = 2))
  }
  fit
}

standardize <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

## transforms used throughout (inverse needed for response-scale contrasts)
transform_y <- function(x, how, eps = 0.001) {
  switch(how,
         logit  = qlogis(pmin(pmax(x, eps), 1 - eps)),
         log05  = log(x + 0.5),     # counts that can be zero (richness)
         log001 = log(x + 0.01),    # fine-root biomass (g); matches Methods
         none   = x)
}
inverse_y <- function(z, how) {
  switch(how,
         logit  = plogis(z),
         log05  = exp(z) - 0.5,
         log001 = exp(z) - 0.01,
         none   = z)
}

summ <- function(x, prefix = "") {
  h <- bayestestR::hdi(x, ci = 0.90)
  out <- c(mean = mean(x), hdi_lo = h$CI_low, hdi_hi = h$CI_high,
           pr_direction = max(mean(x > 0), mean(x < 0)))
  if (nzchar(prefix)) names(out) <- paste0(prefix, "_", names(out))
  out
}

fit_diagnostics <- function(fit) {
  np <- nuts_params(fit)
  s  <- summary(fit)
  ess <- c(s$fixed$Bulk_ESS, unlist(lapply(s$random, function(r) r$Bulk_ESS)), s$spec_pars$Bulk_ESS)
  c(max_rhat = max(rhat(fit), na.rm = TRUE),
    min_bulk_ess = min(ess, na.rm = TRUE),
    divergences = sum(np$Value[np$Parameter == "divergent__"]))
}

## paired depth model on one outcome: y ~ depth + (1 | stand)
## returns the depth effect (subsoil - surface) on the transformed scale and,
## per posterior draw, the population-level response-scale contrast
## surface - subsoil (the sign convention of Figure 4)
paired_depth_model <- function(y_raw, depth_num, stand, how, scale_factor = 1, eps = 0.001) {
  y <- transform_y(y_raw, how, eps = eps)
  ok <- !is.na(y)
  m <- mean(y[ok]); s <- sd(y[ok])
  d <- data.frame(yz = (y[ok] - m) / s, depth_num = depth_num[ok], stand = stand[ok])
  fit <- bayes_fit(yz ~ depth_num + (1 | stand), d)
  dr <- as_draws_df(fit)
  b  <- dr$b_depth_num * s                       # transformed scale
  surf <- inverse_y(m + s * dr$b_Intercept, how)
  sub  <- inverse_y(m + s * (dr$b_Intercept + dr$b_depth_num), how)
  c(n = sum(ok), summ(b, "effect"), summ((surf - sub) * scale_factor, "contrast"),
    surface_pred = mean(surf) * scale_factor, subsoil_pred = mean(sub) * scale_factor,
    fit_diagnostics(fit))
}

## descriptive paired summary (25 paired stands) plus the Bayesian depth effect
## (all samples) for an already-transformed outcome `y_col`; `raw_col` gives
## the scale used for descriptive means/differences (multiplied by `scale`)
paired_summary_bayes <- function(df, y_col, raw_col = y_col, scale = 1) {
  paired <- df |>
    select(node_stand, depth, v = all_of(raw_col)) |>
    tidyr::pivot_wider(names_from = depth, values_from = v) |>
    filter(!is.na(A), !is.na(B))
  dlt <- (paired$B - paired$A) * scale
  y <- df[[y_col]]; m <- mean(y); s <- sd(y)
  d <- data.frame(yz = (y - m) / s, depth_num = df$depth_num, stand = df$node_stand)
  fit <- bayes_fit(yz ~ depth_num + (1 | stand), d)
  b <- as_draws_df(fit)$b_depth_num * s
  tibble(n_samples = nrow(df), n_stands = nrow(paired),
         surface_mean = mean(paired$A) * scale, subsoil_mean = mean(paired$B) * scale,
         mean_difference = mean(dlt), median_difference = median(dlt),
         n_stands_increase = sum(dlt > 0),
         !!!as.list(summ(b, "bayes_effect")), !!!as.list(fit_diagnostics(fit)))
}

here::i_am('R/fit_single_cohort.R')
# Fit the single coherent cohort benchmark using the B I N N E D sample data.
# Requires cmdstanr and a working CmdStan installation.
# Installation guide: https://mc-stan.org/cmdstanr/articles/cmdstanr.html
library(cmdstanr)
library(posterior)
library(here)
# Set this to the directory containing the accompanying CSV and Stan files.
stopifnot(file.exists(here('R/single_cohort_binned.stan')),
          file.exists(here('data/single_cohort_size_counts.csv')),
          file.exists(here('data/single_cohort_truth.csv')))
counts <- read.csv(here('data/single_cohort_size_counts.csv'),
                   stringsAsFactors = FALSE)
truth <- read.csv(here('data/single_cohort_truth.csv'))

stopifnot(all(c('taxonID','siteID','sampleID','day','area_m2',
                'mass_bin','mass_lower_mg','mass_upper_mg','count') %in% names(counts)))
stopifnot(length(unique(counts$taxonID)) == 1L,
          length(unique(counts$siteID)) == 1L)
stopifnot(!anyNA(counts), all(counts$count >= 0), all(counts$area_m2 > 0))

# Sampling unit is sampleID; preserve zero-count bins and zero-count samples.
samples <- unique(counts[c('sampleID','day','area_m2')])
samples <- samples[order(samples$day, samples$sampleID), ]
stopifnot(!anyDuplicated(samples$sampleID))
dates <- sort(unique(samples$day))
bins <- unique(counts[c('mass_bin','mass_lower_mg','mass_upper_mg')])
bins <- bins[order(bins$mass_bin), ]
S <- nrow(samples); J <- length(dates); K <- nrow(bins)
stopifnot(K >= 2, all(bins$mass_lower_mg > 0),
          all(bins$mass_upper_mg > bins$mass_lower_mg))

# Preserve all sample-bin combinations, and verify all expected bins exist.
lookup <- paste(counts$sampleID, counts$mass_bin, sep = '::')
expected <- as.vector(outer(samples$sampleID, bins$mass_bin, paste, sep = '::'))
stopifnot(!anyDuplicated(lookup), length(lookup) == S*K,
          setequal(expected, lookup))
idx <- match(expected, lookup)
y <- matrix(counts$count[idx], nrow = S, ncol = K)
# outer() flattened column-by-column: sample index moves fastest.
stopifnot(all(rowSums(y) >= 0),
          isTRUE(all.equal(as.numeric(y), as.numeric(round(y)))))
# Make sure the observed size domain is contiguous. The model conditions the
# size likelihood on this interval; bins outside it are assumed unobserved.
stopifnot(isTRUE(all.equal(bins$mass_upper_mg[-K],
                           bins$mass_lower_mg[-1], tolerance = 1e-9)))

stan_data <- list(
  J = J, S = S, K = K,
  day = as.numeric(dates),
  date_id = as.integer(match(samples$day, dates)),
  area_m2 = as.numeric(samples$area_m2),
  n_sample = as.integer(rowSums(y)),
  y = unname(y),
  lower_mg = bins$mass_lower_mg,
  upper_mg = bins$mass_upper_mg,
  horizon_day = 365
)


model <- cmdstan_model(here('R/single_cohort_binned.stan'))
fit <- model$sample(
  data = stan_data,
  seed = 1312,
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.95,
  refresh = 200
)

pars <- c('N0','m0_median','sdlog_m0','r','mu','sigma_patch',
          'P_observed_mg_m2','P_total_mg_m2','B_end_mg_m2','N_end_ind_m2')
print(fit$summary(variables = pars)[, c('variable','mean','sd','q5','median','q95','rhat','ess_bulk')],
      row.names = FALSE)

# Compare known simulation truths. These are NEVER used in the likelihood.
stopifnot(any(truth$day == max(dates)), any(truth$day == 365))
truth_parameter <- data.frame(
  variable = c('N0','m0_median','sdlog_m0','r','mu','sigma_patch',
               'P_observed_mg_m2','P_total_mg_m2'),
  true_value = c(1500,0.010,0.20,0.010,0.005,0.25,
                 truth$cumulative_production_mg_m2[truth$day == max(dates)],
                 truth$cumulative_production_mg_m2[truth$day == 365])
)
recovery <- merge(fit$summary(variables = pars), truth_parameter, by = 'variable')
recovery$covered_90pct <- recovery$true_value >= recovery$q5 &
                          recovery$true_value <= recovery$q95
print(recovery[c('variable','true_value','mean','q5','q95','covered_90pct')],
      row.names = FALSE)

# Check sampler diagnostics; do not ignore warnings about divergent transitions.
print(fit$diagnostic_summary())

# Posterior distributions for production and parameters, for subsequent plots.
write.csv(fit$summary(variables = pars),
          file.path(project_dir, 'single_cohort_posterior_summary.csv'), row.names = FALSE)
write.csv(recovery, file.path(project_dir, 'single_cohort_parameter_recovery.csv'),
          row.names = FALSE)
write.csv(as_draws_df(fit$draws(variables = c('N0','m0_median','sdlog_m0',
                    'r','mu','sigma_patch','P_observed_mg_m2','P_total_mg_m2'))),
          file.path(project_dir, 'single_cohort_posterior_draws.csv'), row.names = FALSE)
fit$save_object(file = file.path(project_dir, 'single_cohort_fit.rds'))

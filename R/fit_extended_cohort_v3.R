# Bayesian extended-recruitment cohort model (Model 2).
# Requires cmdstanr, posterior, and a configured CmdStan installation.
# Run from directory containing these files, or change project_dir.
library(cmdstanr)
library(posterior)
library(stanli)
library(here)
here::i_am('R/fit_extended_cohort_v3.R')
stan_file <- file.path(here('R/extended_cohort_binned_v3.stan'))
data_file <- file.path(here('data/extended_cohort_size_counts.csv'))
truth_file <- file.path(here('data/extended_cohort_parameters.csv'))
stopifnot(file.exists(stan_file), file.exists(data_file))

counts <- read.csv(data_file, stringsAsFactors = FALSE)
needed <- c('taxonID', 'siteID', 'sampleID', 'day', 'area_m2',
            'mass_bin', 'mass_lower_mg', 'mass_upper_mg', 'count')
stopifnot(all(needed %in% names(counts)), !anyNA(counts),
          length(unique(counts$taxonID)) == 1L,
          length(unique(counts$siteID)) == 1L,
          all(counts$area_m2 > 0), all(counts$count >= 0),
          all(counts$mass_lower_mg > 0),
          all(counts$mass_upper_mg > counts$mass_lower_mg))

samples <- unique(counts[c('sampleID', 'day', 'area_m2')])
samples <- samples[order(samples$day, samples$sampleID), ]
stopifnot(!anyDuplicated(samples$sampleID), min(samples$day) == 0)
dates <- sort(unique(samples$day))
bins <- unique(counts[c('mass_bin', 'mass_lower_mg', 'mass_upper_mg')])
bins <- bins[order(bins$mass_bin), ]
S <- nrow(samples); J <- length(dates); K <- nrow(bins)
stopifnot(K >= 2, all(diff(dates) > 0),
          isTRUE(all.equal(bins$mass_upper_mg[-K],
                           bins$mass_lower_mg[-1], tolerance = 1e-9)))

# Require an explicit entry for every size class in every spatial sample.
lookup <- paste(counts$sampleID, counts$mass_bin, sep = '::')
expected <- as.vector(outer(samples$sampleID, bins$mass_bin,
                            paste, sep = '::'))
stopifnot(!anyDuplicated(lookup), length(lookup) == S*K,
          setequal(lookup, expected))
y <- matrix(counts$count[match(expected, lookup)], nrow = S, ncol = K)
stopifnot(all(y == round(y)))

stan_data <- list(
  J = J, S = S, K = K,
  day = as.numeric(dates),
  date_id = as.integer(match(samples$day, dates)),
  area_m2 = as.numeric(samples$area_m2),
  n_sample = as.integer(rowSums(y)),
  y = unname(y),
  lower_mg = as.numeric(bins$mass_lower_mg),
  upper_mg = as.numeric(bins$mass_upper_mg),
  horizon_day = 365
)

# These are plausible generic starting values, NOT the simulation truths.
init_extended <- function(chain_id) {
  list(log_N_R = log(1200), log_m0_median = log(0.012),
       sdlog_m0 = 0.25, r = 0.009, mu = 0.006,
       sigma_patch = 0.25, z_patch = rep(0, S),
       rho = -20, D = 40)
}

m <- stanli::stanli_model(stan_file, data = stan_data, fast_math = TRUE)
fit <- stanli::sample_model(
  m,
  seed = 1312,
  chains = 4,
  # init = init_extended,
  warmup = 100, 
  samples = 100,
  delta = 0.90,
  max_depth = 10,
  refresh = 10
)

vars <- c('N_R', 'rho', 'D', 'recruitment_end_day', 'm0_median',
          'sdlog_m0', 'r', 'mu', 'sigma_patch',
          'P_observed_mg_m2', 'P_total_mg_m2')
summary_tbl <- fit$summary(variables = vars)
print(summary_tbl[, c('variable', 'mean', 'sd', 'q5', 'median', 'q95',
                      'rhat', 'ess_bulk')], row.names = FALSE)
print(fit$diagnostic_summary())

# Ground truth is used for validation only, never for fitting or priors.
if (file.exists(truth_file)) {
  truth <- read.csv(truth_file, stringsAsFactors = FALSE)
  lookup_truth <- c(N_R = 'total_recruits_ind_m2',
                    rho = 'recruitment_onset_day',
                    D = 'recruitment_duration_day',
                    recruitment_end_day = 'recruitment_end_day',
                    m0_median = 'm0_median_mg',
                    sdlog_m0 = 'sdlog_m0', r = 'r_day_inverse',
                    mu = 'mu_day_inverse', sigma_patch = 'sigma_patch',
                    P_observed_mg_m2 = 'P_observed_mg_m2',
                    P_total_mg_m2 = 'P_total_mg_m2')
  true_values <- setNames(truth$true_value, truth$parameter)
  recovery <- merge(summary_tbl,
                    data.frame(variable = names(lookup_truth),
                               true_value = unname(true_values[lookup_truth])),
                    by = 'variable')
  recovery$covered_90pct <- with(recovery,
                                 q5 <= true_value & true_value <= q95)
  print(recovery[c('variable', 'true_value', 'mean', 'q5', 'q95',
                   'covered_90pct')], row.names = FALSE)
  write.csv(recovery, file.path(project_dir, 'extended_cohort_recovery_v3.csv'),
            row.names = FALSE)
}

write.csv(summary_tbl, file.path(project_dir, 'extended_cohort_summary_v3.csv'),
          row.names = FALSE)
write.csv(as.data.frame(fit$draws(variables = vars, format = 'draws_df')),
          file.path(project_dir, 'extended_cohort_draws_v3.csv'), row.names = FALSE)
fit$save_object(file.path(project_dir, 'extended_cohort_fit_v3.rds'))

// Single cohort benchmark, revised 2026-11-11.
// Single cohort, deterministic exponential growth, constant mortality.
// Observed counts are Poisson by sample, multinomial across observed mass bins.
// Robust 8-point Gauss-Legendre quadrature integrates the lognormal density
// on each LOG-MASS bin; avoids subtraction of saturated normal CDFs.
// Units: day, mg dry mass, individuals/m^2.
functions {
  real log_interval_prob(real lo, real hi, real meanlog, real sdlog) {
    // Precondition: 0 < lo < hi; sdlog >= 0.03 (see parameters).
    // Points and weights for 8-point Gauss-Legendre integration on [-1,1].
    vector[8] node = [-0.9602898564975363, -0.7966664774136267,
                      -0.5255324099163290, -0.1834346424956498,
                       0.1834346424956498,  0.5255324099163290,
                       0.7966664774136267,  0.9602898564975363]';
    vector[8] weight = [0.1012285362903763, 0.2223810344533745,
                        0.3137066458778873, 0.3626837833783620,
                        0.3626837833783620, 0.3137066458778873,
                        0.2223810344533745, 0.1012285362903763]';
    real a = log(lo);
    real b = log(hi);
    real mid = 0.5 * (a + b);
    real half = 0.5 * (b - a);
    vector[8] log_integrand;
    for (q in 1:8) {
      log_integrand[q] = log(weight[q])
                       + normal_lpdf(mid + half * node[q] | meanlog, sdlog);
    }
    return log(half) + log_sum_exp(log_integrand);
  }
}

data {
  int<lower=1> J;
  int<lower=1> S;
  int<lower=2> K;
  vector<lower=0>[J] day;
  array[S] int<lower=1, upper=J> date_id;
  vector<lower=0>[S] area_m2;
  array[S] int<lower=0> n_sample;
  array[S, K] int<lower=0> y;
  vector<lower=0>[K] lower_mg;
  vector<lower=0>[K] upper_mg;
  real<lower=0> horizon_day;
}

transformed data {
  for (k in 1:K) {
    if (lower_mg[k] <= 0 || upper_mg[k] <= lower_mg[k])
      reject("Mass-bin limits must satisfy 0 < lower < upper.");
  }
  for (s in 1:S) {
    int total = 0;
    if (area_m2[s] <= 0) reject("Sample area must be positive.");
    for (k in 1:K) total += y[s,k];
    if (total != n_sample[s])
      reject("n_sample must equal sum of mass-bin counts.");
  }
}

parameters {
  real log_N0;
  real log_m0_median;
  real<lower=0.03, upper=1.5> sdlog_m0;
  real<lower=0, upper=0.05> r;
  real<lower=0, upper=0.05> mu;
  real<lower=0, upper=1.5> sigma_patch;
  vector[S] z_patch;
}

transformed parameters {
  real<lower=0> N0 = exp(log_N0);
  real<lower=0> m0_median = exp(log_m0_median);
  array[J] vector[K] log_bin_mass;
  vector[J] log_range_prob;
  vector[S] log_expected_sample_count;

  for (j in 1:J) {
    real meanlog_t = log_m0_median + r * day[j];
    for (k in 1:K)
      log_bin_mass[j][k] = log_interval_prob(
        lower_mg[k], upper_mg[k], meanlog_t, sdlog_m0);
    // Probability of an individual falling inside ANY observed bin.
    log_range_prob[j] = log_sum_exp(log_bin_mass[j]);
  }
  for (s in 1:S) {
    int j = date_id[s];
    log_expected_sample_count[s] = log(area_m2[s]) + log_N0
      - mu * day[j] + sigma_patch * z_patch[s]
      - 0.5 * square(sigma_patch) + log_range_prob[j];
  }
}

model {
  log_N0 ~ normal(log(1000), 1.5);
  log_m0_median ~ normal(log(0.02), 1.2);
  sdlog_m0 ~ normal(0, 0.5);
  r ~ lognormal(log(0.01), 0.8);
  mu ~ lognormal(log(0.005), 1.0);
  sigma_patch ~ normal(0, 0.5);
  z_patch ~ std_normal();

  for (s in 1:S) {
    n_sample[s] ~ poisson_log(log_expected_sample_count[s]);
    // Multinomial-logit normalizes internally, without requiring a simplex.
    y[s] ~ multinomial_logit(log_bin_mass[date_id[s]]);
  }
}

generated quantities {
  real mean_m0_mg = exp(log_m0_median + 0.5 * square(sdlog_m0));
  real B0_mg_m2 = N0 * mean_m0_mg;
  real k_net = r - mu;
  real P_total_mg_m2;
  real P_observed_mg_m2;
  real B_end_mg_m2 = B0_mg_m2 * exp(k_net * horizon_day);
  real N_end_ind_m2 = N0 * exp(-mu * horizon_day);
  vector[J] predicted_density;
  vector[J] predicted_mean_mass;
  vector[J] predicted_production_rate;

  if (fabs(k_net) < 1e-8) {
    P_total_mg_m2 = r * B0_mg_m2 * horizon_day;
    P_observed_mg_m2 = r * B0_mg_m2 * max(day);
  } else {
    P_total_mg_m2 = r * B0_mg_m2 * expm1(k_net * horizon_day) / k_net;
    P_observed_mg_m2 = r * B0_mg_m2 * expm1(k_net * max(day)) / k_net;
  }
  for (j in 1:J) {
    predicted_density[j] = N0 * exp(-mu * day[j]);
    predicted_mean_mass[j] = mean_m0_mg * exp(r * day[j]);
    predicted_production_rate[j] = r * predicted_density[j]
                                 * predicted_mean_mass[j];
  }
}

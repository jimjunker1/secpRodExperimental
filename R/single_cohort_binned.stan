// Single cohort: exact interval-censored lognormal sizes + sample totals.
// Units: days, mg dry mass, individuals / m^2.
// No recruitment after t=0, emergence, or size-dependent mortality.
functions {
  real log_bin_prob(real lower_mass, real upper_mass,
                    real meanlog_mass, real sdlog_mass) {
    // Switch CDF tails to avoid catastrophic cancellation near probability 1.
    if (log(lower_mass) > meanlog_mass) {
      return log_diff_exp(
        normal_lccdf(log(lower_mass) | meanlog_mass, sdlog_mass),
        normal_lccdf(log(upper_mass) | meanlog_mass, sdlog_mass));
    } else {
      return log_diff_exp(
        normal_lcdf(log(upper_mass) | meanlog_mass, sdlog_mass),
        normal_lcdf(log(lower_mass) | meanlog_mass, sdlog_mass));
    }
  }
}
data {
  int<lower=1> J;                            // number of sampling dates
  int<lower=1> S;                            // spatial samples
  int<lower=2> K;                            // mass intervals
  vector<lower=0>[J] day;
  array[S] int<lower=1, upper=J> date_id;
  vector<lower=0>[S] area_m2;
  array[S] int<lower=0> n_sample;
  array[S, K] int<lower=0> y;
  vector<lower=0>[K] lower_mg;
  vector<lower=0>[K] upper_mg;
  real<lower=0> horizon_day;
}
parameters {
  real log_N0;
  real log_m0_median;
  real<lower=0> sdlog_m0;
  real<lower=0> r;                      // relative growth per day
  real<lower=0> mu;                     // instantaneous mortality per day
  real<lower=0> sigma_patch;            // spatial variation in log density
  vector[S] z_patch;                    // non-centered spatial effects
}
transformed parameters {
  real<lower=0> N0 = exp(log_N0);
  real<lower=0> m0_median = exp(log_m0_median);
  array[J] vector[K] p_bin;
  vector[S] expected_sample_count;
  for (j in 1:J) {
    vector[K] log_p;
    real meanlog_t = log_m0_median + r * day[j];
    for (k in 1:K)
      log_p[k] = log_bin_prob(lower_mg[k], upper_mg[k],
                            meanlog_t, sdlog_m0);
    // Condition observed bins on being inside the instrument's mass range.
    p_bin[j] = softmax(log_p);
  }
  for (s in 1:S) {
    expected_sample_count[s] = area_m2[s] * N0 * exp(
      -mu * day[date_id[s]] + sigma_patch * z_patch[s]
      - 0.5 * square(sigma_patch));
  }
}
model {
  // Weakly informative starting priors for this benchmark's units/ranges.
  log_N0 ~ normal(log(1000), 1.5);
  log_m0_median ~ normal(log(0.02), 1.2);
  sdlog_m0 ~ normal(0, 0.5);           // half-normal
  r ~ lognormal(log(0.01), 0.8);
  mu ~ lognormal(log(0.005), 1.0);
  sigma_patch ~ normal(0, 0.5);        // half-normal
  z_patch ~ std_normal();

  for (s in 1:S) {
    n_sample[s] ~ poisson(expected_sample_count[s]);
    y[s] ~ multinomial(p_bin[date_id[s]]);
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
  if (abs(k_net) < 1e-8) {
    P_total_mg_m2 = r * B0_mg_m2 * horizon_day;
    P_observed_mg_m2 = r * B0_mg_m2 * max(day);
  } else {
    P_total_mg_m2 = r * B0_mg_m2 * expm1(k_net * horizon_day) / k_net;
    P_observed_mg_m2 = r * B0_mg_m2 * expm1(k_net * max(day)) / k_net;
  }
  for (j in 1:J) {
    predicted_density[j] = N0 * exp(-mu * day[j]);
    predicted_mean_mass[j] = mean_m0_mg * exp(r * day[j]);
    predicted_production_rate[j] = r * predicted_density[j] * predicted_mean_mass[j];
  }
}

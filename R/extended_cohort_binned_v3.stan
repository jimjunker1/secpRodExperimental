// Model 2: single cohort, uniform recruitment over [rho, rho + D].
// Extension of single_cohort_binned_v2.stan.
// Growth dm/dt = r*m; loss rate mu; lognormal mass at recruitment.
// rho is recruitment ONSET (days relative to the first sampling date).
// Each spatial sample has a lognormal multiplicative abundance effect.
// Binned mass likelihood uses nested Gauss-Legendre quadrature in log-mass
// and recruitment time, evaluated on a log scale for numerical stability.
// Units: days, mg dry mass, individuals/m^2.
functions {
  real log_interval_prob(real lo, real hi, real meanlog, real sdlog) {
    vector[8] node = [-0.9602898564975363, -0.7966664774136267,
                      -0.5255324099163290, -0.1834346424956498,
                       0.1834346424956498, 0.5255324099163290,
                       0.7966664774136267, 0.9602898564975363]';
    vector[8] weight = [0.1012285362903763, 0.2223810344533745,
                        0.3137066458778873, 0.3626837833783620,
                        0.3626837833783620, 0.3137066458778873,
                        0.2223810344533745, 0.1012285362903763]';
    real a = log(lo);
    real b = log(hi);
    real mid = 0.5 * (a + b);
    real half = 0.5 * (b - a);
    vector[8] li;
    for (q in 1:8)
      li[q] = log(weight[q])
              + normal_lpdf(mid + half * node[q] | meanlog, sdlog);
    return log(half) + log_sum_exp(li);
  }

  // Integral from calendar day 0 through T of r*population biomass(t).
  // Integrate over recruitment times, splitting at day zero so that recruits
  // entering before zero contribute ONLY growth after day zero.
  real production_through(real T, real rho, real D, real N_R,
                          real mean_m0, real r, real mu) {
    vector[16] node = [-0.989400934991650, -0.944575023073233,
                       -0.865631202387832, -0.755404408355003,
                       -0.617876244402644, -0.458016777657227,
                       -0.281603550779259, -0.095012509837637,
                        0.095012509837637, 0.281603550779259,
                        0.458016777657227, 0.617876244402644,
                        0.755404408355003, 0.865631202387832,
                        0.944575023073233, 0.989400934991650]';
    vector[16] weight = [0.027152459411754, 0.062253523938648,
                         0.095158511682493, 0.124628971255534,
                         0.149595988816577, 0.169156519395003,
                         0.182603415044924, 0.189450610455068,
                         0.189450610455068, 0.182603415044924,
                         0.169156519395003, 0.149595988816577,
                         0.124628971255534, 0.095158511682493,
                         0.062253523938648, 0.027152459411754]';
    real last_s = fmin(rho + D, T);
    real k = r - mu;
    real accum = 0;
    if (T <= 0 || last_s <= rho) return 0;
    // Segment 1: recruitment before observation starts.
    // Segment 2: recruitment after observation starts.
    for (seg in 1:2) {
      real a;
      real b;
      if (seg == 1) {
        a = rho;
        b = fmin(last_s, 0);
      } else {
        a = fmax(rho, 0);
        b = last_s;
      }
      if (b > a) {
        real half = 0.5 * (b - a);
        real mid = 0.5 * (b + a);
        for (q in 1:16) {
          real s = mid + half * node[q];
          real t0 = fmax(0, s);
          real dt = T - t0;
          real growth_integral;
          if (abs(k) < 1e-8)
            growth_integral = dt;
          else
            growth_integral = expm1(k * dt) / k;
          accum += half * weight[q]
            * exp(k * (t0 - s)) * growth_integral;
        }
      }
    }
    return N_R / D * r * mean_m0 * accum;
  }
}

data {
  int<lower=1> J;  // unique sampling dates
  int<lower=1> S;  // distinct spatial samples
  int<lower=2> K;  // size bins
  vector[J] day;
  array[S] int<lower=1, upper=J> date_id;
  vector<lower=0>[S] area_m2;
  array[S] int<lower=0> n_sample;
  array[S, K] int<lower=0> y;
  vector<lower=0>[K] lower_mg;
  vector<lower=0>[K] upper_mg;
  real<lower=0> horizon_day;
}

transformed data {
  for (k in 1:K)
    if (lower_mg[k] <= 0 || upper_mg[k] <= lower_mg[k])
      reject("Mass-bin limits must satisfy 0 < lower < upper.");
  for (s in 1:S) {
    int total = 0;
    if (area_m2[s] <= 0) reject("Sample area must be positive.");
    for (k in 1:K) total += y[s, k];
    if (total != n_sample[s])
      reject("n_sample must equal sum of binned counts.");
  }
  if (min(day) < 0)
    reject("This benchmark uses nonnegative observation days.");
}

parameters {
  real log_N_R;
  real log_m0_median;
  real<lower=0.03, upper=1.5> sdlog_m0;
  real<lower=0, upper=0.05> r;
  real<lower=0, upper=0.05> mu;
  real<lower=0, upper=1.5> sigma_patch;
  vector[S] z_patch;
  // Benchmark-specific bounds: recruitment began before the first observation.
  // These are NOT generic constraints for later datasets.
  real<lower=-120, upper=-0.001> rho;
  real<lower=1, upper=180> D;
}

transformed parameters {
  real<lower=0> N_R = exp(log_N_R);
  real<lower=0> m0_median = exp(log_m0_median);
  array[J] vector[K] log_bin_mass;
  vector[J] log_range_density;
  vector[S] log_expected_sample_count;
  {
    vector[16] node = [-0.989400934991650, -0.944575023073233,
                       -0.865631202387832, -0.755404408355003,
                       -0.617876244402644, -0.458016777657227,
                       -0.281603550779259, -0.095012509837637,
                        0.095012509837637, 0.281603550779259,
                        0.458016777657227, 0.617876244402644,
                        0.755404408355003, 0.865631202387832,
                        0.944575023073233, 0.989400934991650]';
    vector[16] weight = [0.027152459411754, 0.062253523938648,
                         0.095158511682493, 0.124628971255534,
                         0.149595988816577, 0.169156519395003,
                         0.182603415044924, 0.189450610455068,
                         0.189450610455068, 0.182603415044924,
                         0.169156519395003, 0.149595988816577,
                         0.124628971255534, 0.095158511682493,
                         0.062253523938648, 0.027152459411754]';
    for (j in 1:J) {
      real s_end = fmin(day[j], rho + D);
      real half = 0.5 * (s_end - rho);
      real mid = 0.5 * (s_end + rho);
      vector[16] s_quad;
      vector[16] log_weighted_survival;
      // rho < 0 and day >= 0 ensure the integration interval is nonempty.
      for (q in 1:16) {
        s_quad[q] = mid + half * node[q];
        log_weighted_survival[q] = log(weight[q])
          - mu * (day[j] - s_quad[q]);
      }
      for (k in 1:K) {
        vector[16] log_contrib;
        for (q in 1:16) {
          real meanlog_age = log_m0_median + r * (day[j] - s_quad[q]);
          log_contrib[q] = log_weighted_survival[q]
            + log_interval_prob(lower_mg[k], upper_mg[k],
                                meanlog_age, sdlog_m0);
        }
        // Unnormalized density (individuals m^-2) in bin k.
        log_bin_mass[j][k] = log_N_R - log(D) + log(half)
                             + log_sum_exp(log_contrib);
      }
      log_range_density[j] = log_sum_exp(log_bin_mass[j]);
    }
  }
  for (s in 1:S) {
    int j = date_id[s];
    log_expected_sample_count[s] = log(area_m2[s])
      + log_range_density[j] + sigma_patch * z_patch[s]
      - 0.5 * square(sigma_patch);
  }
}

model {
  log_N_R ~ normal(log(1000), 1.5);
  log_m0_median ~ normal(log(0.02), 1.2);
  sdlog_m0 ~ normal(0, 0.5);
  r ~ lognormal(log(0.01), 0.8);
  mu ~ lognormal(log(0.005), 1.0);
  sigma_patch ~ normal(0, 0.5);
  z_patch ~ std_normal();
  rho ~ normal(-15, 30);
  D ~ lognormal(log(30), 0.8);

  for (s in 1:S) {
    n_sample[s] ~ poisson_log(log_expected_sample_count[s]);
    y[s] ~ multinomial_logit(log_bin_mass[date_id[s]]);
  }
}

generated quantities {
  real mean_m0_mg = exp(log_m0_median + 0.5 * square(sdlog_m0));
  real recruitment_end_day = rho + D;
  real P_observed_mg_m2 = production_through(max(day), rho, D, N_R,
                                             mean_m0_mg, r, mu);
  real P_total_mg_m2 = production_through(horizon_day, rho, D, N_R,
                                          mean_m0_mg, r, mu);
  vector[J] predicted_density;
  vector[J] predicted_mean_mass;
  vector[J] predicted_production_rate;
  for (j in 1:J) {
    // One dimensional, mass-integrated recruitment/survival formulas.
    real s_end = fmin(day[j], rho + D);
    real width_s = s_end - rho;
    real age_min = day[j] - s_end;
    real x_count = mu * width_s;
    real k = r - mu;
    real ratio_count;
    real ratio_mass;
    if (abs(x_count) < 1e-8)
      ratio_count = 1 - 0.5 * x_count + square(x_count) / 6;
    else
      ratio_count = -expm1(-x_count) / x_count;
    if (abs(k * width_s) < 1e-8)
      ratio_mass = 1 + 0.5 * k * width_s + square(k * width_s) / 6;
    else
      ratio_mass = expm1(k * width_s) / (k * width_s);
    predicted_density[j] = N_R / D * width_s
      * exp(-mu * age_min) * ratio_count;
    predicted_production_rate[j] = r * N_R / D * mean_m0_mg * width_s
      * exp(k * age_min) * ratio_mass;
    predicted_mean_mass[j] = predicted_production_rate[j]
      / (r * predicted_density[j]);
  }
}

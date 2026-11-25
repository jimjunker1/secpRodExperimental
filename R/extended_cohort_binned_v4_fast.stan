// Model 2: single cohort, uniform recruitment over [rho, rho + D].
// Extension of single_cohort_binned_v2.stan.
// Growth dm/dt = r*m; loss rate mu; lognormal mass at recruitment.
// rho is recruitment ONSET (days relative to the first sampling date).
// Each spatial sample has a lognormal multiplicative abundance effect.
// Binned likelihood analytically integrates recruitment age, followed by
// one 8-point log-mass quadrature. No nested quadrature in usual sampling.
// Units: days, mg dry mass, individuals/m^2.
functions {
  // Log of the integral of the log-mass density over recruitment age.
  // The age integration is analytic for exponential growth, uniform recruitment,
  // and constant mortality. x is LOG(mass), not mass.
  // Output is density per unit log-mass per m2 (unnormalized in x).
  real log_age_mixture_density(real x, real t, real rho, real D,
                               real log_N_R, real alpha, real sigma,
                               real r, real mu) {
    real s_last = fmin(t, rho + D);
    real age_lo = t - s_last;
    real age_hi = t - rho;
    real y = x - alpha;
    if (s_last <= rho) return negative_infinity();

    // For very small r, avoid the 1/r terms in the exact expression.
    // This is only a numerical fallback, not a different biological model.
    if (r < 1e-3) {
      vector[8] nodes = [-0.9602898564975363, -0.7966664774136267,
                         -0.5255324099163290, -0.1834346424956498,
                          0.1834346424956498, 0.5255324099163290,
                          0.7966664774136267, 0.9602898564975363]';
      vector[8] weights = [0.1012285362903763, 0.2223810344533745,
                           0.3137066458778873, 0.3626837833783620,
                           0.3626837833783620, 0.3137066458778873,
                           0.2223810344533745, 0.1012285362903763]';
      real mid = 0.5 * (age_lo + age_hi);
      real half = 0.5 * (age_hi - age_lo);
      vector[8] log_terms;
      for (q in 1:8) {
        real age = mid + half * nodes[q];
        log_terms[q] = log(weights[q]) - mu * age
                      + normal_lpdf(x | alpha + r * age, sigma);
      }
      return log_N_R - log(D) + log(half) + log_sum_exp(log_terms);
    }

    // Complete the square in the recruitment-age integral.
    // Integral exp(-mu*a) Normal(x|alpha+r*a,sigma) da
    // = exp(-mu*y/r+0.5*(mu*sigma/r)^2)/r
    //   * [Phi(z_hi)-Phi(z_lo)], with z=(r*a-y+mu*sigma^2/r)/sigma.
    {
      real z_lo = (r * age_lo - y + mu * square(sigma) / r) / sigma;
      real z_hi = (r * age_hi - y + mu * square(sigma) / r) / sigma;
      real log_cdf_gap;
      // Upper tails: calculate using survival probabilities to prevent
      // subtraction of CDF values very close to 1.
      if (z_lo > 0) {
        log_cdf_gap = log_diff_exp(normal_lccdf(z_lo | 0, 1),
                                   normal_lccdf(z_hi | 0, 1));
      } else {
        log_cdf_gap = log_diff_exp(normal_lcdf(z_hi | 0, 1),
                                   normal_lcdf(z_lo | 0, 1));
      }
      return log_N_R - log(D) - log(r) - mu * y / r
             + 0.5 * square(mu * sigma / r) + log_cdf_gap;
    }
  }

  // Integrate only over the log-mass bin: recruitment age is already
  // integrated analytically in log_age_mixture_density().
  real log_bin_density(real lo, real hi, real t, real rho, real D,
                       real log_N_R, real alpha, real sigma,
                       real r, real mu) {
    vector[8] nodes = [-0.9602898564975363, -0.7966664774136267,
                       -0.5255324099163290, -0.1834346424956498,
                        0.1834346424956498, 0.5255324099163290,
                        0.7966664774136267, 0.9602898564975363]';
    vector[8] weights = [0.1012285362903763, 0.2223810344533745,
                         0.3137066458778873, 0.3626837833783620,
                         0.3626837833783620, 0.3137066458778873,
                         0.2223810344533745, 0.1012285362903763]';
    real mid = 0.5 * (log(lo) + log(hi));
    real half = 0.5 * (log(hi) - log(lo));
    vector[8] li;
    for (q in 1:8) {
      li[q] = log(weights[q])
              + log_age_mixture_density(mid + half * nodes[q], t,
                                        rho, D, log_N_R, alpha,
                                        sigma, r, mu);
    }
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
  for (j in 1:J) {
    for (k in 1:K) {
      log_bin_mass[j][k] = log_bin_density(
        lower_mg[k], upper_mg[k], day[j], rho, D,
        log_N_R, log_m0_median, sdlog_m0, r, mu);
    }
    log_range_density[j] = log_sum_exp(log_bin_mass[j]);
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

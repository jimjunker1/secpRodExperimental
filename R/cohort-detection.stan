// simple_mass_cohort_quadratic_fixed_sigma.stan

data {
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> K;

  vector<lower=0>[I] m_obs;
  array[I] int<lower=1, upper=J> time_id;
  vector<lower=0>[J] s_unique;

  real<lower=0> sigma_x;
  real<lower=0> order_sd;
  real<lower=0> slope_sd;
}

transformed data {
  vector[I] x_obs;
  real x_min;
  real x_max;
  real x_mid;
  real s_max;

  for (i in 1:I) {
    x_obs[i] = log(m_obs[i]);
  }

  x_min = min(x_obs);
  x_max = max(x_obs);
  x_mid = 0.5 * (x_min + x_max);
  s_max = max(s_unique);
}

parameters {
  // ordered ascending, reversed below so cohort 1 is largest/oldest
  ordered[K] alpha_asc;

  vector<lower=0>[K] g;
  vector<lower=0>[K] h;

  simplex[K] theta;
}

transformed parameters {
  vector[K] alpha;
  matrix[J, K] mu;

  for (k in 1:K) {
    alpha[k] = alpha_asc[K - k + 1];
  }

  for (j in 1:J) {
    for (k in 1:K) {
      mu[j, k] =
        alpha[k] +
        g[k] * s_unique[j] -
        h[k] * square(s_unique[j]);
    }
  }
}

model {
  alpha_asc ~ normal(x_mid, 3);

  // Priors tuned for time in days.
  // g is log-mass growth per day.
  g ~ normal(0, 0.03);

  // h controls deceleration.
  h ~ normal(0, 0.0002);

  theta ~ dirichlet(rep_vector(1.0, K));

  // Soft ordering after computing mu.
  if (K > 1) {
    for (j in 1:J) {
      for (k in 1:(K - 1)) {
        target += normal_lcdf(mu[j, k] - mu[j, k + 1] | 0, order_sd);
      }
    }
  }

  // Soft constraint: slope at the final date should not be strongly negative.
  // slope = g - 2 h s
  for (k in 1:K) {
    target += normal_lcdf(g[k] - 2 * h[k] * s_max | 0, slope_sd);
  }

  for (i in 1:I) {
    vector[K] lp;

    for (k in 1:K) {
      lp[k] =
        log(theta[k]) +
        normal_lpdf(x_obs[i] | mu[time_id[i], k], sigma_x);
    }

    target += log_sum_exp(lp);
  }
}

generated quantities {
  vector[I] log_lik;
  matrix[I, K] p_cohort;
  matrix[J, K] m_center;

  for (j in 1:J) {
    for (k in 1:K) {
      m_center[j, k] = exp(mu[j, k]);
    }
  }

  for (i in 1:I) {
    vector[K] lp;
    vector[K] prob;

    for (k in 1:K) {
      lp[k] =
        log(theta[k]) +
        normal_lpdf(log(m_obs[i]) | mu[time_id[i], k], sigma_x);
    }

    log_lik[i] = log_sum_exp(lp);
    prob = softmax(lp);

    for (k in 1:K) {
      p_cohort[i, k] = prob[k];
    }
  }
}

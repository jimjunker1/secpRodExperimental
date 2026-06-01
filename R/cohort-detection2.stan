data {
  int<lower=1> I;                         // observations
  int<lower=1> J;                         // unique sampling dates
  int<lower=1> K;                         // candidate cohorts

  vector<lower=0>[I] m_obs;               // mass or size-class midpoint
  array[I] int<lower=1, upper=J> time_id;
  vector<lower=0>[J] s_unique;            // elapsed time since first sample

  real<lower=0> prior_sigma_x_sd;         // prior SD for within-cohort log-mass spread
  real<lower=0> order_sd;                 // soft ordering strength
  real<lower=0> slope_sd;                 // soft nonnegative final-slope strength
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
  // Stan orders ascending; reversed below so cohort 1 is largest/oldest.
  ordered[K] alpha_asc;

  // Quadratic decelerating log-growth:
  // mu[j,k] = alpha[k] + g[k] * s[j] - h[k] * s[j]^2
  vector<lower=0>[K] g;
  vector<lower=0>[K] h;

  real<lower=0> sigma_x;

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

  // Priors assume time is in days.
  // Adjust if time units change.
  g ~ normal(0, 0.05);
  h ~ normal(0, 0.0002);

  sigma_x ~ normal(0, prior_sigma_x_sd);

  theta ~ dirichlet(rep_vector(1.0, K));

  // Soft ordering after computing trajectories.
  // Cohort 1 should generally be larger than cohort 2, etc.
  if (K > 1) {
    for (j in 1:J) {
      for (k in 1:(K - 1)) {
        target += normal_lcdf(mu[j, k] - mu[j, k + 1] | 0, order_sd);
      }
    }
  }

  // Soft constraint: final slope should not be strongly negative.
  // slope = g - 2 h s
  for (k in 1:K) {
    target += normal_lcdf(g[k] - 2 * h[k] * s_max | 0, slope_sd);
  }

  // Unweighted mass-only mixture likelihood.
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
  vector[I] m_pred;
  matrix[J, K] m_center;
  matrix[J, K] mu_out;

  for (j in 1:J) {
    for (k in 1:K) {
      mu_out[j, k] = mu[j, k];
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

    m_pred[i] = exp(dot_product(prob, to_vector(mu[time_id[i], ])));
  }
}

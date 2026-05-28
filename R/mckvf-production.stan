functions {
  real rbf_eval(real q, vector beta, vector knots, real ell) {
    real out = 0;
    int K = num_elements(beta);

    for (k in 1:K) {
      out += beta[k] * exp(-0.5 * square((q - knots[k]) / ell));
    }

    return out;
  }
}

data {
  int<lower=2> n;
  vector[n] s;                 // elapsed sampling time since first sample
  vector<lower=0>[n] m_obs;
  vector<lower=0>[n] N_obs;
  
  real<lower=0> m_min_obs;
  real<lower=0> m_max_obs;

  int<lower=2> n_basis;
  vector[n_basis] knots;
  real<lower=0> ell;

  int<lower=2> n_quad;
  vector[n_quad] x_quad;           // Gauss-Legendre nodes on [-1, 1]
  vector[n_quad] w_quad;           // Gauss-Legendre weights

  int<lower=2> n_grid;
  vector<lower=0, upper=1>[n_grid] q_grid;
}

parameters {
  real<lower=0, upper=1> W0_frac;
  real<lower=0> Wexit_extra;

  real<lower=0> N0;

  real alpha_g;
  vector[n_basis] beta_g_raw;
  real<lower=0> sigma_beta_g;

  real alpha_h;
  vector[n_basis] beta_h_raw;
  real<lower=0> sigma_beta_h;

  real<lower=0> sigma_s;
  real<lower=0> sigma_N;
}

transformed parameters {
  real<lower=0> W0;
  real<lower=0>Wexit;
  
  vector[n_basis] beta_g;
  vector[n_basis] beta_h;

  vector[n] q_obs;
  vector[n] t_abs;
  vector[n] H_abs;
  vector[n] s_hat;
  vector[n] log_N_hat;
  
  W0 = W0_frac * m_min_obs;
  Wexit = m_max_obs + Wexit_extra;

  beta_g = sigma_beta_g * beta_g_raw;
  beta_h = sigma_beta_h * beta_h_raw;

  for (i in 1:n) {
    q_obs[i] = (m_obs[i] - W0) / (Wexit - W0);
  }

  for (i in 1:n) {
    real tq = 0;
    real Hq = 0;

    for (j in 1:n_quad) {
      real qj;
      real log_gj;
      real log_hj;
      real gj;
      real hj;

      // transform quadrature node from [-1,1] to [0, q_obs[i]]
      qj = 0.5 * q_obs[i] * (x_quad[j] + 1);

      // log_gj = alpha_g + rbf_eval(qj, beta_g, knots, ell);
      // log_hj = alpha_h + rbf_eval(qj, beta_h, knots, ell);
      // 
      // gj = exp(log_gj);
      // hj = exp(log_hj);
      log_gj = fmin(fmax(log_gj, -10), 5);
      log_hj = fmin(fmax(log_hj, -10), 5);

      gj = exp(log_gj);
      hj = exp(log_hj);

      tq += w_quad[j] / gj;
      Hq += w_quad[j] * hj;
    }

    // dm = (Wexit - W0) dq
    t_abs[i] = 0.5 * q_obs[i] * (Wexit - W0) * tq;
    H_abs[i] = 0.5 * q_obs[i] * (Wexit - W0) * Hq;

    s_hat[i] = t_abs[i] - t_abs[1];
    log_N_hat[i] = log(N0) - H_abs[i];
    // N_hat[i] = N0 * exp(-H_abs[i]);
  }
}

model {
  // Ensure observed masses lie inside cohort domain
  // if (W0 >= min(m_obs)) {
  //   target += negative_infinity();
  // }
  // 
  // if (Wexit <= max(m_obs)) {
  //   target += negative_infinity();
  // }

  // Boundary priors
  target += lognormal_lpdf(W0 | log(0.006), 0.5);
  target += lognormal_lpdf(Wexit | log(5), 0.2);
  Wexit_extra ~ exponential(2);

  // Starting density prior
  N0 ~ lognormal(log(500), 0.5);

  // Smooth growth prior
  alpha_g ~ normal(log(0.03), 1.0);
  beta_g_raw ~ normal(0, 1);
  sigma_beta_g ~ exponential(2);

  // Smooth mortality-per-growth prior
  alpha_h ~ normal(log(0.01), 1.0);
  beta_h_raw ~ normal(0, 1);
  sigma_beta_h ~ exponential(2);

  // Observation error
  sigma_s ~ exponential(5);
  sigma_N ~ exponential(5);

  // Likelihoods
  for(i in 2:n){
  s[i] ~ lognormal(log(s_hat[i]), sigma_s);
  }
  // N_obs ~ lognormal(log(N_hat), sigma_N);
  N_obs ~ lognormal(log_N_hat, sigma_N);
}

generated quantities {
  vector[n_grid] m_grid;
  vector[n_grid] t_grid;
  vector[n_grid] log_N_grid;
  vector[n_grid] N_grid;
  vector[n_grid] g_grid;
  vector[n_grid] h_grid;
  vector[n_grid] mu_grid;

  vector[n_grid - 1] P_interval;
  real P_total;

  real t_exit;
  real N_exit;
  real B_exit;

  // Evaluate smooth functions on prediction grid
  for (i in 1:n_grid) {
    real log_gi;
    real log_hi;

    m_grid[i] = W0 + q_grid[i] * (Wexit - W0);

    log_gi = alpha_g + rbf_eval(q_grid[i], beta_g, knots, ell);
    log_hi = alpha_h + rbf_eval(q_grid[i], beta_h, knots, ell);

    // Numerical safety bounds
    log_gi = fmin(fmax(log_gi, -10), 5);
    log_hi = fmin(fmax(log_hi, -10), 5);

    g_grid[i] = exp(log_gi);
    h_grid[i] = exp(log_hi);

    // mu(m) = g(m) * h(m), where h(m) = mu(m) / g(m)
    mu_grid[i] = g_grid[i] * h_grid[i];
  }

  // Initial conditions at cohort start
  t_grid[1] = 0;
  log_N_grid[1] = log(N0);
  N_grid[1] = N0;

  // Integrate through mass space
  for (i in 2:n_grid) {
    real dm;
    real dt_inc;
    real H_inc;

    dm = m_grid[i] - m_grid[i - 1];

    // dt = dm / g(m)
    dt_inc =
      0.5 *
      (1 / g_grid[i - 1] + 1 / g_grid[i]) *
      dm;

    // cumulative mortality hazard in mass space:
    // dH = h(m) dm = mu(m) / g(m) dm
    H_inc =
      0.5 *
      (h_grid[i - 1] + h_grid[i]) *
      dm;

    t_grid[i] = t_grid[i - 1] + dt_inc;

    log_N_grid[i] = log_N_grid[i - 1] - H_inc;
    N_grid[i] = exp(log_N_grid[i]);
  }

  // Production in mass space:
  // P = integral N(m) dm
  for (i in 1:(n_grid - 1)) {
    P_interval[i] =
      0.5 *
      (N_grid[i] + N_grid[i + 1]) *
      (m_grid[i + 1] - m_grid[i]);
  }

  P_total = sum(P_interval);

  // Terminal quantities
  t_exit = t_grid[n_grid];
  N_exit = N_grid[n_grid];
  B_exit = N_exit * Wexit;
}

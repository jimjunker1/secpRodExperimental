data {
  int<lower=1> n;
  vector[n] s;                         // elapsed sampling time since first sample
  vector<lower=0>[n] m_obs;
  vector<lower=0>[n] N_obs;

  int<lower=2> n_grid;
  vector<lower=0, upper=1>[n_grid] q_grid;  // standardized 0-1 prediction grid

  real<lower=0, upper=1> m_exit_frac;  // e.g. 0.99
}

parameters {
  real<lower=0> Winf;
  real<lower=0> K;

  real<lower=0, upper=1> W1_frac;      // W1 / Winf
  real<lower=0, upper=1> W0_rel;       // W0 / W1, ensures W0 < W1

  real<lower=0> N1;
  real<lower=0> Z;

  real<lower=0> sigma_m;
  real<lower=0> sigma_N;
}

transformed parameters {
  real<lower=0> W1;
  real<lower=0> W0;

  real<lower=0, upper=1> B0;
  real<lower=0, upper=1> B1;

  real<lower=0> age1;
  real<lower=0> N0;

  vector[n] tau_obs;
  vector[n] m_hat;
  vector[n] N_hat;

  W1 = W1_frac * Winf;
  W0 = W0_rel * W1;

  B0 = 1 - pow(W0 / Winf, 1.0 / 3.0);
  B1 = 1 - pow(W1 / Winf, 1.0 / 3.0);

  // Implied cohort age at first sample
  age1 = -log(B1 / B0) / K;

  // Implied density at cohort start
  N0 = N1 * exp(Z * age1);

  for (i in 1:n) {
    tau_obs[i] = age1 + s[i];

    m_hat[i] =
      Winf * pow(
        1 - B0 * exp(-K * tau_obs[i]),
        3
      );

    N_hat[i] =
      N0 * exp(-Z * tau_obs[i]);
  }
}

model {
  // Priors
  Winf ~ lognormal(log(5), 0.15);
  K ~ lognormal(log(0.055), 0.25);

  // Prior on first observed mass
  target += lognormal_lpdf(W1 | log(m_obs[1]), 0.15);

  // Prior on true initial mass
  target += lognormal_lpdf(W0 | log(0.006), 0.50);

  N1 ~ lognormal(log(N_obs[1]), 0.20);
  Z ~ lognormal(log(0.02), 0.30);

  sigma_m ~ exponential(20);
  sigma_N ~ exponential(20);

  // Likelihood
  m_obs ~ lognormal(log(m_hat), sigma_m);
  N_obs ~ lognormal(log(N_hat), sigma_N);
}

generated quantities {
  real<lower=0> t_exit;
  real<lower=0> W_exit;
  real<lower=0> N_exit;
  real<lower=0> B_exit;

  vector[n_grid] tau_pred;
  vector[n_grid] m_pred;
  vector[n_grid] N_pred;
  vector[n_grid] g_pred;

  vector[n_grid - 1] P_interval;
  real<lower=0> P_total;

  // Cohort endpoint: mass reaches m_exit_frac * Winf
  t_exit =
    -log(
      (1 - pow(m_exit_frac, 1.0 / 3.0)) / B0
    ) / K;

  W_exit = m_exit_frac * Winf;
  N_exit = N0 * exp(-Z * t_exit);
  B_exit = W_exit * N_exit;

  for (i in 1:n_grid) {
    tau_pred[i] = q_grid[i] * t_exit;

    m_pred[i] =
      Winf * pow(
        1 - B0 * exp(-K * tau_pred[i]),
        3
      );

    N_pred[i] =
      N0 * exp(-Z * tau_pred[i]);

    // dm / dt over cohort age
    g_pred[i] =
      3 * K * Winf *
      B0 * exp(-K * tau_pred[i]) *
      square(1 - B0 * exp(-K * tau_pred[i]));
  }

  for (i in 1:(n_grid - 1)) {
    P_interval[i] =
      0.5 *
      (N_pred[i] * g_pred[i] + N_pred[i + 1] * g_pred[i + 1]) *
      (tau_pred[i + 1] - tau_pred[i]);
  }

  P_total = sum(P_interval);
}

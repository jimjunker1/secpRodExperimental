extract_indexed_draws <- function(fit, variable, index_name = "i") {
  
  as_draws_df(fit$draws(variable)) |>
    select(.draw, starts_with(variable)) |>
    pivot_longer(
      cols = starts_with(variable),
      names_to = "parameter",
      values_to = "value"
    ) |>
    mutate(
      !!index_name := as.integer(stringr::str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))
    ) |>
    select(.draw, !!index_name, value)
}


prep_mass_unweighted <- function(dat) {
  dat_mass <- dat |>
    group_by(time, m) |>
    summarise(
      N = sum(N),
      .groups = "drop"
    ) |>
    arrange(time, m) |>
    mutate(
      obs_id = row_number(),
      time_id = as.integer(factor(time, levels = sort(unique(time))))
    )
  
  time_lookup <- dat_mass |>
    distinct(time, time_id) |>
    arrange(time_id)
  
  list(
    dat_mass = dat_mass,
    time_lookup = time_lookup,
    s_unique = time_lookup$time - min(time_lookup$time)
  )
}

make_stan_data_simple <- function(
    dat,
    K,
    prior_sigma_x_sd = 0.20,
    order_sd = 0.10,
    slope_sd = 0.02
) {
  prep <- prep_mass_unweighted(dat)
  dat_mass <- prep$dat_mass
  
  list(
    I = nrow(dat_mass),
    J = length(unique(dat_mass$time_id)),
    K = K,
    
    m_obs = dat_mass$m,
    time_id = dat_mass$time_id,
    s_unique = prep$s_unique,
    
    prior_sigma_x_sd = prior_sigma_x_sd,
    order_sd = order_sd,
    slope_sd = slope_sd
  )
}

fit_quad_one <- function(
    dat,
    K,
    prior_sigma_x_sd = 0.20,
    order_sd = 0.10,
    slope_sd = 0.02,
    seed = 1,
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    use_optim_init = TRUE
) {
  stan_data <- make_stan_data_simple(
    dat = dat,
    K = K,
    prior_sigma_x_sd = prior_sigma_x_sd,
    order_sd = order_sd,
    slope_sd = slope_sd
  )
  
  init_list <- NULL
  
  if (use_optim_init) {
    opt <- try(
      mod_quad$optimize(
        data = stan_data,
        seed = seed,
        jacobian = FALSE
      ),
      silent = TRUE
    )
    
    if (!inherits(opt, "try-error")) {
      opt_draw <- posterior::as_draws_df(opt$draws())
      
      # Convert one optimized draw into a named list.
      # Remove CmdStan metadata columns if present.
      opt_vals <- opt_draw |>
        dplyr::select(
          -dplyr::any_of(c(".chain", ".iteration", ".draw"))
        )
      
      init_base <- as.list(opt_vals[1, ])
      
      # Same init for each chain, with small jitter added by Stan internally if needed.
      init_list <- replicate(chains, init_base, simplify = FALSE)
    } else {
      message("Optimization failed for K = ", K, "; using random initialization.")
    }
  }
  
  fit <- mod$sample(
    data = stan_data,
    init = init_list,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed + 1000,
    adapt_delta = 0.97,
    max_treedepth = 12
  )
  
  list(
    K = K,
    fit = fit,
    stan_data = stan_data,
    opt = if (exists("opt")) opt else NULL
  )
}



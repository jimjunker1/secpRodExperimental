# Model 2: a single cohort recruited uniformly between day -10 and +10.
# Base R; new random realization (not byte-for-byte equal to supplied CSVs).
set.seed(20261125)
out_dir <- '.'
N_R <- 1500; tau <- -10; D <- 20; recruit_end <- tau + D
m0_median <- 0.010; sdlog_m0 <- 0.20
r <- 0.010; mu <- 0.005; sigma_patch <- 0.25
area_m2 <- 0.1; n_rep <- 5L
sample_days <- seq(0, 330, 30); grid_days <- 0:365
edges <- exp(seq(log(0.001), log(2), length.out = 65))
q <- N_R / D
mean_m0 <- m0_median * exp(sdlog_m0^2 / 2)
# Stable integral of exp(rate*(t-s)) from tau to min(t, tau + D).
exp_age_integral <- function(t, rate) {
  end <- min(t, recruit_end)
  if (end <= tau) return(0)
  width <- end - tau
  x <- rate * width
  scale <- if (abs(x) < 1e-7) 1 + x/2 + x^2/6 else expm1(x)/x
  width * exp(rate * (t-end)) * scale
}
N_fun <- function(t) q * exp_age_integral(t, -mu)
B_fun <- function(t) q * mean_m0 * exp_age_integral(t, r-mu)
P_rate <- function(t) r * B_fun(t)
P_cum <- function(T) if (T == 0) 0 else integrate(P_rate, lower=0, upper=T,
  subdivisions=200L, rel.tol=1e-10)$value
# Population median mass of survivor-weighted lognormal mixture.
median_fun <- function(t) {
  end <- min(t, recruit_end)
  f <- function(m) {
    weighted_cdf <- integrate(function(s) {
      q * exp(-mu*(t-s)) * pnorm((log(m)-log(m0_median)-r*(t-s))/sdlog_m0)
    }, lower=tau, upper=end)$value
    weighted_cdf / N_fun(t) - 0.5
  }
  uniroot(f, interval=c(1e-6, 10))$root
}
truth <- do.call(rbind, lapply(grid_days, function(t) {
  N <- N_fun(t); B <- B_fun(t)
  data.frame(day=t, density_ind_m2=N, median_mass_mg=median_fun(t),
    mean_mass_mg=B/N, biomass_mg_m2=B,
    production_rate_mg_m2_day=P_rate(t),
    cumulative_production_mg_m2=P_cum(t))
}))
individuals <- list(); binned <- list(); sample_summary <- list()
i <- 0L
for (j in seq_along(sample_days)) for (rep in seq_len(n_rep)) {
  t <- sample_days[j]
  id <- sprintf('S%02d_R%02d', j, rep)
  h <- rlnorm(1L, -sigma_patch^2/2, sigma_patch)
  count <- rpois(1L, area_m2*N_fun(t)*h)
  if (count > 0L) {
    # Recruitment time among SURVIVORS: density proportional to exp(mu*s).
    u <- runif(count)
    s <- tau + log1p(u*expm1(mu*(min(t, recruit_end)-tau)))/mu
    mass <- rlnorm(count, log(m0_median) + r*(t-s), sdlog_m0)
    bin_id <- findInterval(mass, edges, rightmost.closed=FALSE)
    stopifnot(all(bin_id >= 1L & bin_id <= 64L))
  } else {mass <- numeric(); bin_id <- integer()}
  i <- i + 1L
  individuals[[i]] <- data.frame(taxonID='SYN_COHORT_EXT', siteID='SITE_A',
    sampleID=id, day=t, area_m2=area_m2, mass_mg=mass)
  sample_summary[[i]] <- data.frame(taxonID='SYN_COHORT_EXT', siteID='SITE_A',
    sampleID=id, day=t, area_m2=area_m2, observed_count=count,
    observed_density_ind_m2=count/area_m2)
  tab <- tabulate(bin_id, nbins=64L)
  binned[[i]] <- data.frame(taxonID='SYN_COHORT_EXT', siteID='SITE_A',
    sampleID=id, day=t, area_m2=area_m2, mass_bin=seq_len(64L),
    mass_lower_mg=head(edges, -1), mass_upper_mg=tail(edges, -1), count=tab)
}
write.csv(do.call(rbind, individuals), file.path(out_dir,'extended_cohort_individuals.csv'),row.names=FALSE)
write.csv(do.call(rbind, binned), file.path(out_dir,'extended_cohort_size_counts.csv'),row.names=FALSE)
write.csv(do.call(rbind, sample_summary), file.path(out_dir,'extended_cohort_sample_summary.csv'),row.names=FALSE)
write.csv(truth, file.path(out_dir,'extended_cohort_truth.csv'),row.names=FALSE)
parameter <- c('total_recruits_ind_m2','recruitment_onset_day','recruitment_duration_day',
  'recruitment_end_day','recruitment_rate_ind_m2_day','m0_median_mg','sdlog_m0',
  'r_day_inverse','mu_day_inverse','sigma_patch','area_m2','n_replicates',
  'last_sample_day','P_observed_mg_m2','P_total_mg_m2')
true_value <- c(N_R,tau,D,recruit_end,q,m0_median,sdlog_m0,r,mu,
  sigma_patch,area_m2,n_rep,max(sample_days),P_cum(330),P_cum(365))
write.csv(data.frame(parameter,true_value),file.path(out_dir,'extended_cohort_parameters.csv'),row.names=FALSE)
cat(sprintf('P(0,330)=%.3f; P(0,365)=%.3f mg m^-2\n',P_cum(330),P_cum(365)))

# Single-cohort benchmark simulator; base R only.
# Same demographic scenario as the accompanying example CSV files.
# A fresh random realization is produced by R (not identical to the Python-seeded CSV).
set.seed(1312)

# Units: mass mg dry weight, time days, density individuals / m^2.
N0 <- 1500              # density at day 0, individuals / m^2
m0_median <- 0.010      # initial median individual mass, mg
m0_sdlog <- 0.20        # initial spread in log mass (not growth-rate variation)
r <- 0.010              # constant relative growth, day^-1
mu <- 0.005             # constant mortality/loss hazard, day^-1
area_m2 <- 0.10         # area sampled in each independent spatial sample
n_replicates <- 5L
sample_days <- seq(0, 330, by = 30)
truth_days <- 0:365

# Lognormal initial mass distribution, deterministic growth dm/dt = r*m.
mean_m0 <- m0_median * exp(m0_sdlog^2 / 2)
density <- function(t) N0 * exp(-mu * t)
mean_mass <- function(t) mean_m0 * exp(r * t)
median_mass <- function(t) m0_median * exp(r * t)
biomass <- function(t) density(t) * mean_mass(t)  # mg / m^2
production_rate <- function(t) r * biomass(t)      # mg / m^2 / day
cum_production <- function(t) {
  if (abs(r - mu) < 1e-12) return(N0 * mean_m0 * r * t)
  N0 * mean_m0 * r * expm1((r - mu) * t) / (r - mu)
}

truth <- data.frame(
  day = truth_days,
  density_ind_m2 = density(truth_days),
  median_mass_mg = median_mass(truth_days),
  mean_mass_mg = mean_mass(truth_days),
  biomass_mg_m2 = biomass(truth_days),
  production_rate_mg_m2_day = production_rate(truth_days),
  cumulative_production_mg_m2 = vapply(truth_days, cum_production, numeric(1))
)

# Individual sizes are generated conditionally on the sampled count;
# all animals at a given initial mass grow at the SAME r.
# Sample-to-sample lognormal variation in density represents spatial patchiness.
spatial_sdlog <- 0.25
individual_records <- vector('list', length(sample_days) * n_replicates)
index <- 0L
for (j in seq_along(sample_days)) {
  t <- sample_days[j]
  for (s in seq_len(n_replicates)) {
    sample_id <- sprintf('S%02d_R%02d', j, s)
    patch_factor <- rlnorm(1, -spatial_sdlog^2 / 2, spatial_sdlog)
    n_caught <- rpois(1, area_m2 * density(t) * patch_factor)
    index <- index + 1L
    individual_records[[index]] <- data.frame(
      taxonID = rep('SYN_COHORT', n_caught),
      siteID = rep('SITE_A', n_caught),
      sampleID = rep(sample_id, n_caught),
      day = rep(t, n_caught),
      area_m2 = rep(area_m2, n_caught),
      mass_mg = rlnorm(n_caught, log(m0_median) + r * t, m0_sdlog)
    )
  }
}
individuals <- do.call(rbind, individual_records)

# Common size classes across all samples: ideal for a later finite-volume model.
bin_edges <- exp(seq(log(0.001), log(2), length.out = 65))
ids <- cut(individuals$mass_mg, breaks = bin_edges, include.lowest = TRUE,
           labels = FALSE, right = FALSE)
stopifnot(!anyNA(ids))
individuals$bin_id <- ids

# Include ALL size classes, including zeros, for every spatial sample.
sample_meta <- unique(individuals[c('taxonID','siteID','sampleID','day','area_m2')])
# A sample with zero caught animals must still appear in sample_meta; rebuild explicitly.
sample_meta <- expand.grid(day = sample_days, replicate = seq_len(n_replicates))
sample_meta$sampleID <- sprintf('S%02d_R%02d',
                                match(sample_meta$day, sample_days), sample_meta$replicate)
sample_meta$taxonID <- 'SYN_COHORT'
sample_meta$siteID <- 'SITE_A'
sample_meta$area_m2 <- area_m2

table_rows <- vector('list', nrow(sample_meta))
for (s in seq_len(nrow(sample_meta))) {
  meta <- sample_meta[s, ]
  tab <- tabulate(individuals$bin_id[individuals$sampleID == meta$sampleID],
                  nbins = length(bin_edges) - 1L)
  table_rows[[s]] <- data.frame(
    taxonID = meta$taxonID, siteID = meta$siteID,
    sampleID = meta$sampleID, day = meta$day, area_m2 = meta$area_m2,
    mass_bin = seq_along(tab),
    mass_lower_mg = head(bin_edges, -1),
    mass_upper_mg = tail(bin_edges, -1),
    count = as.integer(tab)
  )
}
samples <- do.call(rbind, table_rows)

write.csv(individuals[c('taxonID','siteID','sampleID','day','area_m2','mass_mg')],
          'data/single_cohort_individuals.csv', row.names = FALSE)
write.csv(samples, 'data/single_cohort_size_counts.csv', row.names = FALSE)
write.csv(truth, 'data/single_cohort_truth.csv', row.names = FALSE)
cat(sprintf('True production over 365 days: %.3f mg / m^2\n', cum_production(365)))
cat(sprintf('True production through final sample (day 330): %.3f mg / m^2\n', cum_production(330)))

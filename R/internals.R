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


library(tidyverse)
library(lubridate)
library(MMWRweek)

# -------------------------------------------------------------------
# 1. Load data
# -------------------------------------------------------------------

flu_data <- read_csv(
  "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/refs/heads/main/target-data/target-hospital-admissions.csv"
)

flu_data |>
  mutate(
    year = MMWRweek(date)$MMWRyear,
    week = MMWRweek(date)$MMWRweek,
    resp_season = ifelse(week >= 40, year, year - 1)
  ) |>
  group_by(location, location_name, resp_season) |>
  arrange(date, .by_group = TRUE) |>
  mutate(resp_season_week = seq_along(week)) |>
  ungroup() -> recent_flu_cleaned

forecast_date <- max(recent_flu_cleaned$date, na.rm = TRUE)

# -------------------------------------------------------------------
# 2. Setup
# -------------------------------------------------------------------

source("R/scenariocast-function.R")
load("data/locations-data.rda")
load("data/clean-historic-flu-spline-combined-one-each.rda")

curr_resp_season <- 2025
quantiles_needed <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

# -------------------------------------------------------------------
# 3. MULTI-LOCATION LOOP
# -------------------------------------------------------------------

loc_results <- vector("list", length = nrow(locations))
names(loc_results) <- locations$location_name

for (loc_name in locations$location_name) {
  
  message("\n========================================")
  message("Running scenariocast for: ", loc_name)
  message("========================================")
  
  # Extract location code
  loc_code <- recent_flu_cleaned %>%
    distinct(location, location_name) %>%
    filter(location_name == loc_name) %>%
    pull(location) %>%
    unique()
  
  # -------------------------------------------------------------------
  # A. RECENT DATA → TRAJECTORIES
  # -------------------------------------------------------------------
  
  recent_flu_cleaned |>
    filter(location == loc_code,
           location_name == loc_name,
           resp_season == curr_resp_season,
           date <= forecast_date) |>
    mutate(
      value = value + 1,
      curr_weekly_change = log(lead(value) / value)
    ) |>
    select(resp_season_week, value, curr_weekly_change) |>
    get_traj_forecast(
      db = traj_db_comb,
      nsamps = 1000,
      recent_weeks_touse = 20,
      resp_week_range    = 2,
      forecast_horizon   = 20
    ) |>
    mutate(
      forecast = forecast - 1,
      forecast = pmax(forecast, 0)
    ) -> forecast_trajectories
  
  # -------------------------------------------------------------------
  # B. FILTER EXTREME OUTLIER TRAJECTORIES
  # -------------------------------------------------------------------
  
  prev_season <- curr_resp_season - 1
  
  max_prev_value <- recent_flu_cleaned %>%
    filter(location_name == loc_name,
           resp_season == prev_season) %>%
    summarise(max_value = max(value, na.rm = TRUE)) %>%
    pull(max_value)
  
  id_summary <- forecast_trajectories %>%
    group_by(id) %>%
    summarise(max_forecast = max(forecast, na.rm = TRUE), .groups = "drop")
  
  extreme_ids <- id_summary %>%
    filter(max_forecast > 3 * max_prev_value) %>%
    pull(id)
  
  forecast_trajectories_clean <- forecast_trajectories %>%
    filter(!id %in% extreme_ids)
  
  forecast_trajectories <- forecast_trajectories_clean
  
  # -------------------------------------------------------------------
  # C. OBS + TRAJECTORY MERGED
  # -------------------------------------------------------------------
  
  obs_loc <- recent_flu_cleaned %>%
    filter(location == loc_code,
           resp_season == curr_resp_season,
           date <= forecast_date) %>%
    select(location, location_name, resp_season, resp_season_week, value)
  
  combined_full <- forecast_trajectories %>%
    group_by(id) %>%
    group_modify(~ {
      bind_rows(
        obs_loc,
        .x %>%
          mutate(location = loc_code,
                 location_name = loc_name,
                 resp_season = curr_resp_season) %>%
          select(location, location_name, resp_season,
                 resp_season_week, forecast) %>%
          rename(value = forecast)
      ) %>%
        arrange(resp_season_week)
    }) %>%
    ungroup()
  
  # -------------------------------------------------------------------
  # D. PEAK WEEK & MAGNITUDE
  # -------------------------------------------------------------------
  
  peak_summary <- combined_full %>%
    group_by(id) %>%
    summarise(
      peak_week = resp_season_week[which.max(value)],
      peak_magnitude = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      location       = loc_code,
      forecast_date  = forecast_date,
      target         = "seasonal_peak",
      peak_magnitude = round(peak_magnitude)
    )
  
  # -------------------------------------------------------------------
  # E. WEEKLY QUANTILES
  # -------------------------------------------------------------------
  
  weekly_q <- forecast_trajectories %>%
    group_by(resp_season_week) %>%
    summarise(
      qs = list(quantile(forecast, probs = quantiles_needed)),
      .groups = "drop"
    ) %>%
    mutate(horizon = seq_along(resp_season_week) - 1) %>%
    filter(horizon <= 4) %>%
    unnest_wider(qs) %>%
    pivot_longer(-c(resp_season_week, horizon), names_to = "quantile") %>%
    mutate(
      quantile        = as.numeric(gsub("[\\%,]", "", quantile)) / 100,
      target          = "wk inc flu hosp",
      location        = loc_code,
      reference_date  = forecast_date + 7,
      target_end_date = forecast_date + 7 + horizon * 7,
      output_type     = "quantile",
      output_type_id  = quantile,
      value           = round(value)
    ) %>%
    select(reference_date, target, horizon, target_end_date,
           location, output_type, output_type_id, value)
  
  # -------------------------------------------------------------------
  # F. PEAK WEEK QUANTILES
  # -------------------------------------------------------------------
  
  peak_week_q <- peak_summary %>%
    summarise(qs = list(quantile(peak_week, probs = quantiles_needed))) %>%
    unnest_wider(qs) %>%
    pivot_longer(cols = everything(), names_to = "quantile") %>%
    mutate(
      quantile       = as.numeric(gsub("[\\%,]", "", quantile)) / 100,
      target         = "seasonal_peak_week",
      output_type    = "quantile",
      output_type_id = quantile,
      location       = loc_code,
      reference_date = forecast_date + 7
    )
  
  # -------------------------------------------------------------------
  # G. PEAK MAGNITUDE QUANTILES
  # -------------------------------------------------------------------
  
  peak_mag_q <- peak_summary %>%
    summarise(qs = list(quantile(peak_magnitude, probs = quantiles_needed))) %>%
    unnest_wider(qs) %>%
    pivot_longer(cols = everything(), names_to = "quantile") %>%
    mutate(
      quantile       = as.numeric(gsub("[\\%,]", "", quantile)) / 100,
      target         = "seasonal_peak_inc flu hosp",
      output_type    = "quantile",
      output_type_id = quantile,
      location       = loc_code,
      reference_date = forecast_date + 7
    )
  
  # -------------------------------------------------------------------
  # H. FINAL FORECAST FOR THIS LOCATION
  # -------------------------------------------------------------------
  
  loc_results[[loc_name]] <- bind_rows(
    weekly_q,
    peak_week_q,
    peak_mag_q
  )
}

# -------------------------------------------------------------------
# 4. COMBINE ALL LOCATIONS
# -------------------------------------------------------------------

final_submission <- bind_rows(loc_results, .id = "location_name")

dir.create("processed-data/rt-forecasts", recursive = TRUE, showWarnings = FALSE)

write_csv(
  final_submission,
  paste0(
    "processed-data/rt-forecasts/",
    forecast_date + 7,
    "-UGA_flucast-scenariocast.csv"
  )
)


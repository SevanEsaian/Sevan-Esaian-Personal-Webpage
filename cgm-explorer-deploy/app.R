# =============================================================================
# CGM Explorer: continuous glucose monitoring in children (ages 6 - 16)
#
# An interactive dashboard that turns 14 days of glucose sensor readings
# (one reading every 5 minutes) into clear, usable insights.
#
# TABS
#   1. Start Here       Overview of the data, a guided tour, and the color key
#   2. Glucose Overview How groups of children compare, hour by hour
#   3. Child Explorer   One child's full 14 days, with low episodes marked
#   4. False Lows       Sensor-pressure false lows, and what cleaning them changes
#   5. Models & Forecasts  Group differences over time, and short-term prediction
#
# DATA NOTE FOR CODE READERS
#   The data are generated in SECTION 2 so the app runs anywhere without
#   patient data. To use real data, replace SECTION 2 with your own objects:
#     participants     one row per child (participant_id, group, age, sex, bmi)
#     period_summary   one row per child x period (Day / Night)
#     all_traces       named list of tibbles (day = 0 to 14, glucose = mg/dL)
#
# PACKAGES (mgcv and nlme ship with R):
#   install.packages(c("shiny", "bslib", "dplyr", "ggplot2", "ggiraph"))
#
# DEPLOY
#   rsconnect::deployApp(appDir = "cgm-explorer-deploy", appName = "cgm-explorer")
# =============================================================================

library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(ggiraph)
library(mgcv)


# =============================================================================
# SECTION 1. SETTINGS: colors, glucose ranges, fonts
# =============================================================================

# ---- Color logic: BLUE = good, RED = bad -------------------------------------
# Glucose outside the target range is bad in BOTH directions, so every
# out-of-range value uses a red. A dark red marks lows; lighter and stronger
# reds mark highs. Groups use neutral colors so they never read as good or bad.
range_levels <- c("Low", "In range", "High", "Very high")
range_colors <- c("Low"       = "#7A1230",   # dark red
                  "In range"  = "#2B7BBA",   # blue
                  "High"      = "#E8736C",   # light red
                  "Very high" = "#D62828")   # strong red
range_details <- c("Low"       = "70 mg/dL or below",
                   "In range"  = "71 - 139 mg/dL",
                   "High"      = "140 - 179 mg/dL",
                   "Very high" = "180 mg/dL and above")

group_levels <- c("Type 2 Diabetes", "Prediabetes", "Healthy (No Prediabetes)")
group_colors <- c("Type 2 Diabetes"          = "#6A3D9A",   # purple
                  "Prediabetes"              = "#C98B00",   # gold
                  "Healthy (No Prediabetes)" = "#5C6B7A")   # slate
period_colors <- c("Day" = "#C98B00", "Night" = "#3D3B8E")  # gold / indigo

# ---- Fonts and sizes ---------------------------------------------------------
plot_font    <- "Arial"
base_font_px <- 15
bar_label_size <- base_font_px / ggplot2::.pt * 0.85

tooltip_css <- paste0("background:white;color:#1F2937;border:1px solid #9AA5B1;",
                      "padding:6px 10px;font-family:Arial,Helvetica,sans-serif;",
                      "font-size:14px;border-radius:4px;")

# Classify a glucose value (mg/dL) into one of the four ranges
classify_glucose <- function(x) {
  factor(
    case_when(x <= 70  ~ "Low",
              x < 140  ~ "In range",
              x < 180  ~ "High",
              TRUE     ~ "Very high"),
    levels = range_levels
  )
}


# =============================================================================
# SECTION 2. DATA: 100 children, 14 days of readings each
# =============================================================================

set.seed(2026)

# 100 children: 34 / 33 / 33 per group
group_sizes <- c(34, 33, 33)
age_ranges  <- list(c(6, 16), c(6, 13), c(8, 14))

participants <- bind_rows(lapply(1:3, function(i) {
  tibble(
    group = group_levels[i],
    age   = sample(age_ranges[[i]][1]:age_ranges[[i]][2], group_sizes[i], replace = TRUE),
    sex   = sample(c("Female", "Male"), group_sizes[i], replace = TRUE)
  )
})) %>%
  mutate(
    participant_id = sprintf("P%03d", row_number()),
    group          = factor(group, levels = group_levels),
    # Child-level random effects. Day and Night share them, so each child's
    # two periods stay realistically paired.
    re_events   = rnorm(n()),
    re_duration = rnorm(n()),
    re_nadir    = rnorm(n()),
    re_low_pct  = rnorm(n()),
    re_glucose  = rnorm(n()),
    re_cv       = rnorm(n())
  )

# Day vs. night averages for each group (events are counts per 14 days)
sim_params <- tibble(
  group     = factor(group_levels, levels = group_levels),
  ev_day    = c(3.0, 2.0, 1.5),   ev_night  = c(14, 9, 6),       # false lows
  dur_day   = c(18, 16, 14),      dur_night = c(45, 38, 30),     # minutes
  nad_day   = c(66, 67, 68),      nad_night = c(48, 53, 57),     # lowest mg/dL
  tbr_day   = c(2.0, 1.5, 1.0),   tbr_night = c(8.0, 5.5, 3.5),  # % time low
  glu_day   = c(170, 118, 99),    glu_night = c(148, 105, 89),   # mean mg/dL
  cv_day    = c(34, 27, 20),      cv_night  = c(27, 21, 16)      # % variation
)

# One row per child for one period ("Day" or "Night")
simulate_period <- function(df, period) {
  is_night <- period == "Night"
  df %>%
    left_join(sim_params, by = "group") %>%
    mutate(
      period        = period,
      comp_events   = pmax(1, rpois(n(), (if (is_night) ev_night else ev_day) * exp(0.25 * re_events))),
      comp_duration = pmax(5, (if (is_night) dur_night else dur_day) + 3 * re_duration +
                             rnorm(n(), 0, if (is_night) 7 else 3.5)),
      comp_nadir    = pmax(40, (if (is_night) nad_night else nad_day) + 3 * re_nadir +
                             rnorm(n(), 0, if (is_night) 4 else 2.5)),
      tbr           = pmax(0.1, (if (is_night) tbr_night else tbr_day) * exp(0.15 * re_low_pct) +
                             rnorm(n(), 0, if (is_night) 1.0 else 0.4)),
      mean_glu      = pmax(60, (if (is_night) glu_night else glu_day) + 8 * re_glucose +
                             rnorm(n(), 0, if (is_night) 5 else 4)),
      cv            = pmax(8, (if (is_night) cv_night else cv_day) + 2.5 * re_cv +
                             rnorm(n(), 0, if (is_night) 2.5 else 2))
    ) %>%
    select(participant_id, group, age, sex, period,
           comp_events, comp_duration, comp_nadir, tbr, mean_glu, cv)
}

# BMI: one value per child (separate seed keeps the data above unchanged)
set.seed(2027)
bmi_table <- participants %>%
  select(participant_id, group) %>%
  mutate(
    g   = as.integer(group),
    bmi = rnorm(n(), mean = c(34, 28, 19.5)[g], sd = c(7, 5, 2.5)[g]),
    bmi = pmin(c(50, 42, 24.9)[g], pmax(c(19, 16, 15)[g], bmi)),
    bmi = round(bmi, 1)
  ) %>%
  select(participant_id, bmi)

participants <- participants %>% left_join(bmi_table, by = "participant_id")

# One row per child x period, with a fixed horizontal offset per child so a
# child's Day and Night points line up and can be joined by a line
period_summary <- bind_rows(
  simulate_period(participants, "Day"),
  simulate_period(participants, "Night")
) %>%
  mutate(
    period = factor(period, levels = c("Day", "Night")),
    across(c(comp_duration, comp_nadir, tbr, mean_glu, cv), ~ round(.x, 1)),
    jitter_x = (((as.integer(sub("P", "", participant_id)) * 0.61803398875) %% 1) - 0.5) * 0.28
  ) %>%
  left_join(bmi_table, by = "participant_id")

# Drop-down labels, e.g. "P001 · Type 2 Diabetes · age 12"
participant_choices <- setNames(
  participants$participant_id,
  paste0(participants$participant_id, " \u00b7 ", participants$group,
         " \u00b7 age ", participants$age)
)

# ---- One child's 14-day glucose trace -----------------------------------------
# One reading every 5 minutes for 14 days (4,032 readings). Each trace has a
# day/night rhythm, meal spikes, and false lows that match period_summary.
# Groups differ in how volatile the traces are and how often readings run high.
simulate_trace <- function(pid, n_days = 14, step_min = 5) {
  rows      <- period_summary %>% filter(participant_id == pid)
  day_row   <- rows %>% filter(period == "Day")
  night_row <- rows %>% filter(period == "Night")
  grp       <- as.character(day_row$group)
  
  set.seed(5000 + as.integer(sub("P", "", pid)))   # same child, same trace
  
  sp <- switch(grp,
               "Type 2 Diabetes" = list(amp_mult = 1.6, p_big = 0.25, big_mult = 1.4,
                                        noise_sd = 3.5, noise_ar = 0.94,
                                        nad_lo = 40, nad_hi = 58, extra_hypo = 6,
                                        night_dur_scale = 3.0, extra_dur = 150),
               "Prediabetes"     = list(amp_mult = 1.6, p_big = 0.15, big_mult = 1.5,
                                        noise_sd = 2.5, noise_ar = 0.92,
                                        nad_lo = 42, nad_hi = 58, extra_hypo = 4,
                                        night_dur_scale = 2.5, extra_dur = 120),
               list(amp_mult = 1.3, p_big = 0.05, big_mult = 2.5,            # Healthy
                    noise_sd = 1.2, noise_ar = 0.90,
                    nad_lo = 52, nad_hi = 66, extra_hypo = 0,
                    night_dur_scale = 1.0, extra_dur = 0)
  )
  
  vol <- exp(rnorm(1, 0, 0.25))     # child-to-child variation
  
  t_min <- seq(0, n_days * 1440 - step_min, by = step_min)
  n     <- length(t_min)
  hour  <- (t_min %% 1440) / 60
  
  # Smooth day/night baseline: lowest near 3 AM, highest near 3 PM
  w        <- (1 + sin(2 * pi * (hour - 9) / 24)) / 2
  baseline <- night_row$mean_glu + (day_row$mean_glu - night_row$mean_glu) * w
  
  # Meal spikes: breakfast, lunch, dinner every day
  amp_scale <- day_row$cv / 100 * day_row$mean_glu * sp$amp_mult * vol
  meals     <- c(7.5, 12.5, 18.5) * 60
  bumps     <- numeric(n)
  for (d in 0:(n_days - 1)) {
    for (m in meals) {
      t0  <- d * 1440 + m + rnorm(1, 0, 20)
      amp <- amp_scale * runif(1, 0.5, 1.2)
      if (runif(1) < sp$p_big) amp <- amp * sp$big_mult
      lag   <- pmax(0, t_min - t0)
      bumps <- bumps + amp * (lag / 60) * exp(1 - lag / 60) * (t_min >= t0)
    }
  }
  
  # Slowly wandering sensor noise
  noise   <- as.numeric(stats::filter(rnorm(n, 0, sp$noise_sd * vol), sp$noise_ar, method = "recursive"))
  glucose <- baseline + bumps + noise
  glucose <- ifelse(glucose > 200, 200 + 50 * tanh((glucose - 200) / 50), glucose)  # soft cap near 250
  
  # False lows: sharp drop, short plateau, quick recovery. Each lasts 80+ minutes.
  add_events <- function(glucose, k, hour_range, dur_mean, nadir_mean,
                         nad_lo = sp$nad_lo, nad_hi = sp$nad_hi,
                         dur_sd = 6, min_dur = 80) {
    for (i in seq_len(k)) {
      day_i <- sample(0:(n_days - 1), 1)
      h     <- runif(1, hour_range[1], hour_range[2]) %% 24
      t0    <- day_i * 1440 + h * 60
      dur   <- max(min_dur, rnorm(1, dur_mean, dur_sd))
      nadir <- min(nad_hi, max(nad_lo, rnorm(1, min(nadir_mean, nad_hi - 3), 5)))
      shape <- pmax(0, pmin(1, (t_min - t0) / 15, (t0 + dur - t_min) / 15))
      ev    <- glucose - (glucose - nadir) * shape + rnorm(n, 0, 2) * (shape > 0)
      glucose <- ifelse(shape > 0, pmin(glucose, ev), glucose)
    }
    glucose
  }
  glucose <- add_events(glucose, night_row$comp_events, c(22, 30),          # 10 PM - 6 AM
                        night_row$comp_duration * sp$night_dur_scale,
                        night_row$comp_nadir,
                        dur_sd = 6 * sp$night_dur_scale)
  glucose <- add_events(glucose, day_row$comp_events, c(6, 22),             # 6 AM - 10 PM
                        day_row$comp_duration, day_row$comp_nadir)
  
  # Extra long overnight lows (61 - 69 mg/dL) for prediabetes and type 2
  if (sp$extra_hypo > 0) {
    glucose <- add_events(glucose, sp$extra_hypo, c(22, 30), sp$extra_dur, 65,
                          nad_lo = 61, nad_hi = 69, dur_sd = 20)
  }
  
  tibble(day     = t_min / 1440,                              # 0 to 14
         glucose = round(pmin(250, pmax(40, glucose))))       # whole mg/dL
}

all_traces <- setNames(lapply(participants$participant_id, simulate_trace),
                       participants$participant_id)

# ---- One long table: every reading for every child -----------------------------
weekday_names <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")   # day 1 = Monday

readings <- bind_rows(lapply(names(all_traces), function(pid) {
  tr <- all_traces[[pid]]
  tibble(participant_id = pid, t_min = round(tr$day * 1440), glucose = tr$glucose)
})) %>%
  left_join(participants %>% select(participant_id, group, age, sex, bmi), by = "participant_id") %>%
  mutate(
    day_num  = t_min %/% 1440 + 1,
    hour_dec = (t_min %% 1440) / 60,
    hour     = floor(hour_dec),
    bin      = floor(hour_dec * 4) / 4,                       # 15-minute bins
    weekday  = factor(weekday_names[(day_num - 1) %% 7 + 1], levels = weekday_names),
    period   = factor(ifelse(hour >= 22 | hour < 6, "Night", "Day"), levels = c("Day", "Night"))
  ) %>%
  select(participant_id, group, age, sex, bmi, day_num, hour_dec, hour, bin,
         weekday, period, glucose)

# ---- Child x day x period averages (used by the statistical model) ------------
# 100 children x 14 days x 2 periods = 2,800 rows, so the model fits in seconds.
model_table <- readings %>%
  group_by(participant_id, group, age, sex, bmi, day_num, period) %>%
  summarise(mean_glu = mean(glucose),
            sd_glu   = sd(glucose),
            pct_low  = 100 * mean(glucose <= 70),
            .groups = "drop") %>%
  mutate(
    night    = period,                                        # factor: Day, Night
    time_idx = 2 * (day_num - 1) + as.integer(period == "Night"),
    weekend  = factor(ifelse((day_num - 1) %% 7 >= 5, "Weekend", "Weekday"),
                      levels = c("Weekday", "Weekend")),
    sex            = factor(sex, levels = c("Female", "Male")),
    participant_id = factor(participant_id)
  ) %>%
  arrange(participant_id, time_idx) %>%
  as.data.frame()


# =============================================================================
# SECTION 3. HELPER FUNCTIONS: formatting, summaries, glucose measures
# =============================================================================

# 7.75 -> "7:45 AM"
format_clock <- function(h) {
  hh  <- floor(h)
  mm  <- round((h - hh) * 60)
  h12 <- ifelse(hh %% 12 == 0, 12, hh %% 12)
  paste0(h12, ":", sprintf("%02d", mm), " ", ifelse(hh < 12, "AM", "PM"))
}

clock_breaks <- list(breaks = c(0, 6, 12, 18, 24),
                     labels = c("12 AM", "6 AM", "12 PM", "6 PM", "12 AM"))

# One row per child from any slice of `readings`
child_summary <- function(d) {
  d %>%
    group_by(participant_id, group) %>%
    summarise(mean_glu  = mean(glucose),
              sd_glu    = sd(glucose),
              pct_low   = 100 * mean(glucose <= 70),
              pct_low54 = 100 * mean(glucose < 54),
              pct_inr   = 100 * mean(glucose > 70 & glucose < 140),
              pct_high  = 100 * mean(glucose >= 140 & glucose < 180),
              pct_vhigh = 100 * mean(glucose >= 180),
              .groups = "drop") %>%
    mutate(cv        = 100 * sd_glu / mean_glu,
           gmi       = 3.31 + 0.02392 * mean_glu,
           pct_above = pct_high + pct_vhigh)
}

# Mean amplitude of glycemic excursions (MAGE): the average size of every
# glucose swing larger than 1 standard deviation of the trace
mage <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  amps <- numeric(0)
  dir <- 0                      # 0 = undecided, 1 = rising, -1 = falling
  lo <- x[1]; hi <- x[1]; pivot <- x[1]; cand <- x[1]
  for (v in x[-1]) {
    if (dir == 0) {
      lo <- min(lo, v); hi <- max(hi, v)
      if (v - lo >= s)      { dir <- 1;  pivot <- lo; cand <- v }
      else if (hi - v >= s) { dir <- -1; pivot <- hi; cand <- v }
    } else if (dir == 1) {
      if (v > cand) cand <- v
      else if (cand - v >= s) { amps <- c(amps, cand - pivot); pivot <- cand; dir <- -1; cand <- v }
    } else {
      if (v < cand) cand <- v
      else if (v - cand >= s) { amps <- c(amps, pivot - cand); pivot <- cand; dir <- 1; cand <- v }
    }
  }
  if (length(amps) == 0) NA_real_ else mean(amps)
}

# The four summary measures used in the False Lows tab, for one glucose trace
#   GMI     = 3.31 + 0.02392 x mean glucose
#   J-index = 0.001 x (mean + SD)^2
#   ADRR    = average over days of (largest daily low risk + largest daily high risk)
summary_measure_names <- c("HbA1c estimate (GMI, %)", "Big swings (MAGE, mg/dL)",
                           "Level + swings (J-index)", "Risk of lows and highs (ADRR)")

calc_measures <- function(g, day_index) {
  m  <- mean(g)
  s  <- sd(g)
  f  <- 1.509 * (log(g)^1.084 - 5.381)
  r  <- 10 * f^2
  lr <- tapply(ifelse(f < 0, r, 0), day_index, max)
  hr <- tapply(ifelse(f > 0, r, 0), day_index, max)
  c(3.31 + 0.02392 * m, mage(g), 0.001 * (m + s)^2, mean(lr + hr))
}

# Find low episodes: 15+ minutes (3 readings) in a row at or below `thr`
detect_episodes <- function(g, thr = 70, min_len = 3) {
  r    <- rle(g <= thr)
  ends <- cumsum(r$lengths)
  strt <- ends - r$lengths + 1
  keep <- r$values & r$lengths >= min_len
  if (!any(keep)) return(NULL)
  s <- strt[keep]
  e <- ends[keep]
  tibble(start = s, end = e, minutes = (e - s + 1) * 5,
         lowest = vapply(seq_along(s), function(i) min(g[s[i]:e[i]]), numeric(1)))
}

# Measures a user can pick in the Glucose Overview tab
overview_measure_choices <- c("Time in range, % (71 - 139 mg/dL)" = "pct_inr",
                              "Time low, % (70 or below)"         = "pct_low",
                              "Time below 54 mg/dL, %"            = "pct_low54",
                              "Time high, % (140 and above)"      = "pct_above",
                              "Average glucose (mg/dL)"           = "mean_glu",
                              "Glucose variation, CV (%)"         = "cv",
                              "Estimated HbA1c (GMI, %)"          = "gmi")

# Raise every reading below `floor_val` up to `floor_val` ("clean" false lows)
clean_readings <- function(g, floor_val = NULL) {
  if (is.null(floor_val)) g else pmax(g, floor_val)
}

# ---- Measures before vs. after cleaning (cached per cutoff) --------------------
measures_cache <- new.env()
measures_for <- function(floor_val = NULL) {
  key <- if (is.null(floor_val)) "raw" else as.character(floor_val)
  if (is.null(measures_cache[[key]])) {
    measures_cache[[key]] <- bind_rows(lapply(names(all_traces), function(pid) {
      tr <- all_traces[[pid]]
      tibble(participant_id = pid,
             measure = summary_measure_names,
             value   = unname(calc_measures(clean_readings(tr$glucose, floor_val), floor(tr$day))))
    }))
  }
  measures_cache[[key]]
}
invisible(measures_for(NULL))   # ready when the app starts

# Average % change per measure, with a 95% interval and a plain verdict.
# A measure counts as "unchanged" when its whole interval sits inside +/- margin.
summarise_change <- function(d, margin, by_group) {
  d <- if (by_group) mutate(d, grp = group) else mutate(d, grp = factor("All selected groups"))
  d %>%
    group_by(measure, grp) %>%
    summarise(n = n(), pct_mean = mean(pct), pct_sd = sd(pct), .groups = "drop") %>%
    filter(n >= 2) %>%
    mutate(se    = pct_sd / sqrt(n),
           tcrit = qt(0.975, n - 1),
           lo    = pct_mean - tcrit * se,
           hi    = pct_mean + tcrit * se,
           unchanged = lo > -margin & hi < margin)
}

# ---- Forecasting helpers --------------------------------------------------------
forecast_methods <- c(persist = "Repeat latest reading",
                      trend   = "Continue recent slope",
                      ar      = "Regression model")
forecast_colors  <- c("Actual" = "black",
                      "Repeat latest reading" = "#C98B00",
                      "Continue recent slope" = "#5C6B7A",
                      "Regression model"      = "#2B7BBA")

# Forecast `h` readings ahead (h x 5 minutes) for one child. Each model trains
# only on days BEFORE its test block, so it never sees the future:
#   train days 1-8 -> test 9-10 | train 1-10 -> test 11-12 | train 1-12 -> test 13-14
run_forecast <- function(tr, h, floor_val = NULL, test_starts = c(9, 11, 13)) {
  g <- clean_readings(tr$glucose, floor_val)
  n <- length(g)
  hour <- (((seq_len(n) - 1) * 5) %% 1440) / 60
  lag_by  <- function(k) c(rep(NA_real_, k), g[seq_len(n - k)])
  lead_by <- function(k) c(g[(k + 1):n], rep(NA_real_, k))
  
  d <- data.frame(
    t = seq_len(n), day = floor(tr$day) + 1, actual = lead_by(h),
    g0 = g, l1 = lag_by(1), l2 = lag_by(2), l3 = lag_by(3), l6 = lag_by(6), l12 = lag_by(12),
    s1 = sin(2 * pi * hour / 24), c1 = cos(2 * pi * hour / 24),
    s2 = sin(4 * pi * hour / 24), c2 = cos(4 * pi * hour / 24)
  )
  d$persist <- d$g0
  d$trend   <- d$g0 + 0.5 * (d$g0 - d$l3) / 3 * h
  d$xh      <- ((d$t + h - 1) * 5) / 60 - (d$day - 1) * 24    # clock time being predicted
  
  out <- vector("list", length(test_starts))
  for (k in seq_along(test_starts)) {
    ts           <- test_starts[k]
    first_test_t <- min(d$t[d$day == ts])
    train <- d[d$t + h < first_test_t & complete.cases(d), ]
    test  <- d[d$day %in% c(ts, ts + 1) & !is.na(d$actual) & !is.na(d$l12), ]
    fit   <- lm(actual ~ g0 + l1 + l2 + l3 + l6 + l12 + s1 + c1 + s2 + c2, data = train)
    test$ar    <- as.numeric(predict(fit, test))
    test$sigma <- summary(fit)$sigma
    test$fold  <- paste0("Test days ", ts, "-", ts + 1)
    out[[k]] <- test
  }
  bind_rows(out)
}

# One row per forecast x method
forecasts_long <- function(p) {
  bind_rows(lapply(names(forecast_methods), function(m) {
    tibble(fold = p$fold, day = p$day, xh = p$xh, actual = p$actual,
           method = forecast_methods[[m]], pred = p[[m]])
  }))
}

forecast_accuracy <- function(lp) {
  lp %>%
    group_by(Method = method) %>%
    summarise(`Typical error (mg/dL)` = sqrt(mean((pred - actual)^2)),
              `Average error (mg/dL)` = mean(abs(pred - actual)),
              `Average error (%)`     = 100 * mean(abs(pred - actual) / actual),
              `Forecasts close enough (%)` =
                100 * mean(ifelse(actual <= 100, abs(pred - actual) <= 15,
                                  abs(pred - actual) / actual <= 0.15)),
              .groups = "drop") %>%
    mutate(Method = factor(Method, levels = unname(forecast_methods))) %>%
    arrange(Method) %>%
    mutate(Method = as.character(Method))
}

# ROC curve and area under it (higher score = more likely a low)
roc_curve <- function(score, y) {
  thr <- sort(unique(score), decreasing = TRUE)
  tpr <- c(0, vapply(thr, function(t) mean(score[y] >= t), numeric(1)), 1)
  fpr <- c(0, vapply(thr, function(t) mean(score[!y] >= t), numeric(1)), 1)
  list(df = tibble(fpr = fpr, tpr = tpr),
       auc = sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2))
}

# ---- Model helpers -----------------------------------------------------------
model_outcomes <- c("Average glucose (mg/dL)"     = "mean_glu",
                    "Glucose variation, SD (mg/dL)" = "sd_glu",
                    "Time below 70 mg/dL (%)"     = "pct_low")

model_predictors <- c("Group"                    = "group",
                      "Night vs. day"            = "night",
                      "Weekend vs. weekday"      = "weekend",
                      "Day of monitoring"        = "day_num",
                      "Age"                      = "age",
                      "Body mass index"          = "bmi",
                      "Sex"                      = "sex",
                      "Group \u00d7 Night"       = "grp_night")

clean_term_label <- function(x) {
  map <- c("groupPrediabetes"                         = "Prediabetes (vs. Type 2)",
           "groupHealthy (No Prediabetes)"            = "Healthy (vs. Type 2)",
           "nightNight"                               = "Night (vs. day)",
           "weekendWeekend"                           = "Weekend (vs. weekday)",
           "day_num"                                  = "Each extra day of monitoring",
           "age"                                      = "Each extra year of age",
           "bmi"                                      = "Each extra BMI point",
           "sexMale"                                  = "Male (vs. female)",
           "groupPrediabetes:nightNight"              = "Prediabetes \u00d7 Night",
           "groupHealthy (No Prediabetes):nightNight" = "Healthy \u00d7 Night")
  ifelse(x %in% names(map), unname(map[x]), x)
}


# =============================================================================
# SECTION 4. PLOT HELPERS
# =============================================================================

theme_dashboard <- function() {
  theme_minimal(base_size = base_font_px, base_family = plot_font) +
    theme(
      legend.position    = "none",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90"),
      axis.title         = element_text(color = "grey30"),
      strip.text         = element_text(face = "bold"),
      plot.margin        = margin(8, 12, 8, 8)
    )
}

# Legend on top, no title (used by most charts)
legend_top <- function(rows = 1) {
  list(theme(legend.position = "top", legend.title = element_blank(),
             legend.text = element_text(size = base_font_px)),
       guides(color = guide_legend(nrow = rows), fill = guide_legend(nrow = rows)))
}

# ggiraph wrapper: hover highlights one item and fades the rest
make_girafe <- function(p, w, h, fade = 0.2) {
  girafe(
    ggobj = p, width_svg = w, height_svg = h,
    options = list(
      opts_hover(css = "stroke:black;stroke-width:2px;fill-opacity:1;stroke-opacity:1;"),
      opts_hover_inv(css = paste0("opacity:", fade, ";")),
      opts_tooltip(css = tooltip_css),
      opts_toolbar(saveaspng = FALSE),
      opts_sizing(rescale = TRUE)
    )
  )
}

# Fill scale for mean glucose: dark red (low) -> blue (in range) -> reds (high)
glucose_fill_scale <- function(name) {
  scale_fill_gradientn(
    name = name,
    colours = unname(range_colors[c("Low", "In range", "High", "Very high")]),
    values  = c(0, 0.25, 0.5, 1),                # 50, 100, 150, 250 mg/dL
    limits  = c(50, 250), oob = scales::squish
  )
}

# Stacked bars: average % of time in each range, one bar per group
plot_range_bars <- function(cs) {
  r <- cs %>%
    group_by(group) %>%
    summarise(Low = mean(pct_low), `In range` = mean(pct_inr),
              High = mean(pct_high), `Very high` = mean(pct_vhigh),
              .groups = "drop")
  long <- bind_rows(lapply(range_levels, function(l) {
    tibble(group = r$group, range = l, pct = r[[l]])
  })) %>%
    mutate(range = factor(range, levels = rev(range_levels)),
           tip   = paste0("<b>", group, "</b><br/>", range, " (", range_details[as.character(range)],
                          "): ", sprintf("%.1f", pct), "% of the time"),
           lab   = ifelse(pct >= 4, paste0(round(pct), "%"), ""))
  
  p <- ggplot(long, aes(x = pct, y = group, fill = range)) +
    geom_col_interactive(aes(tooltip = tip, data_id = paste(group, range)),
                         width = 0.65, color = "white") +
    geom_text(aes(label = lab, group = range), position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = bar_label_size, family = plot_font) +
    scale_fill_manual(values = range_colors, breaks = range_levels, drop = FALSE) +
    scale_y_discrete(labels = scales::label_wrap(16)) +
    scale_x_continuous(limits = c(0, 100)) +
    labs(x = "Average % of time", y = NULL) +
    theme_dashboard() + legend_top(2)
  make_girafe(p, 7, 5.5)
}

# Box plot per group with one hoverable dot per child
plot_box_by_group <- function(df, values, ylab, width = 7, height = 5.5) {
  df$value <- values
  df$tip <- paste0("<b>", df$participant_id, "</b><br/>", df$group, "<br/>",
                   ylab, ": ", sprintf("%.1f", df$value))
  p <- ggplot(df, aes(x = group, y = value)) +
    geom_boxplot(aes(fill = group, color = group), alpha = 0.25,
                 outlier.shape = NA, width = 0.55, linewidth = 0.8) +
    geom_point_interactive(aes(color = group, tooltip = tip, data_id = participant_id),
                           position = position_jitter(width = 0.15, height = 0, seed = 7),
                           size = 3, alpha = 0.85) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    scale_x_discrete(labels = scales::label_wrap(14)) +
    labs(x = NULL, y = ylab) +
    theme_dashboard()
  make_girafe(p, width, height, fade = 0.15)
}

# Small colored chip used in the color key
range_chip <- function(level) {
  span(class = "range-chip", style = paste0("background:", range_colors[[level]], ";"),
       level, span(class = "range-chip-detail", range_details[[level]]))
}

# Two-line intro shown at the top of each tab
tab_intro <- function(shows, notice) {
  div(class = "tab-intro",
      p(tags$b("What this shows: "), shows),
      p(tags$b("What to notice: "), notice))
}

# Checkbox labels colored like the chart lines
range_choice_names <- lapply(range_levels, function(l) {
  span(span(l, style = paste0("color:", range_colors[[l]], ";font-weight:bold;")),
       " ", range_details[[l]])
})


# =============================================================================
# SECTION 5. STYLE
# =============================================================================

navy <- "#12344D"

app_css <- paste0("
  body { font-family: Arial, Helvetica, sans-serif; font-size: 16px; color: #1F2937; background: #F6F8FA; }

  /* Top bar */
  .navbar { background: ", navy, " !important; }
  .navbar .navbar-brand, .navbar .nav-link { color: #FFFFFF !important; }
  .navbar .nav-link { font-size: 17px; opacity: 0.8; }
  .navbar .nav-link.active, .navbar .nav-link:hover { opacity: 1; border-bottom: 3px solid #2B7BBA; }
  .brand-main { font-size: 22px; font-weight: bold; line-height: 1.1; }

  /* Cards and headers */
  .card { border: 1px solid #DDE3EA; border-radius: 8px; box-shadow: none; }
  .card-header { font-size: 18px; font-weight: bold; background: #FFFFFF; }
  .card-footer { font-size: 14px; color: #4B5563; background: #FFFFFF; }
  .card-footer p { margin: 0; }
  .sidebar h5 { font-size: 18px; font-weight: bold; }
  .sidebar h6 { font-size: 15px; font-weight: bold; margin-top: 4px; }

  /* Two-line tab intro */
  .tab-intro { border-left: 4px solid #2B7BBA; background: #FFFFFF; padding: 10px 16px;
               margin-bottom: 12px; font-size: 17px; }
  .tab-intro p { margin: 4px 0; }

  /* Landing page */
  .hero { background: ", navy, "; color: #FFFFFF; padding: 36px 40px; border-radius: 8px; margin-bottom: 16px; }
  .hero h1 { font-size: 34px; font-weight: bold; margin: 0 0 10px 0; }
  .hero p  { font-size: 19px; margin: 0 0 8px 0; text-align: justify; }
  .hero .color-key { margin-top: 16px; }
  .color-key-title { font-size: 16px; margin-bottom: 8px; }
  .range-chip { display: inline-block; color: #FFFFFF; font-weight: bold; padding: 6px 12px;
                border-radius: 6px; margin: 0 8px 8px 0; font-size: 15px; }
  .range-chip-detail { font-weight: normal; margin-left: 8px; opacity: 0.95; }
  .tour-card h5 { font-size: 19px; font-weight: bold; }
  .tour-card p  { margin-bottom: 8px; }
  .tour-try { color: #4B5563; font-size: 15px; }
  .about-line { color: #4B5563; font-size: 15px; }

  /* Verdict badges */
  .badge-row { display: flex; flex-wrap: wrap; gap: 10px; margin: 12px 0; }
  .badge-pill { display: inline-block; padding: 6px 12px; border-radius: 6px; font-weight: bold; }
  .badge-good { background: #2B7BBA; color: #FFFFFF; }
  .badge-bad  { background: #D62828; color: #FFFFFF; }
  .badge-info { background: ", navy, "; color: #FFFFFF; }
  .verdict-headline { font-size: 22px; font-weight: bold; margin-bottom: 8px; }

  @media (max-width: 768px) {
    .hero { padding: 24px 20px; }
    .hero h1 { font-size: 26px; }
  }
")


# =============================================================================
# SECTION 6. USER INTERFACE
# =============================================================================

n_children  <- n_distinct(participants$participant_id)
n_readings  <- format(nrow(readings), big.mark = ",")
overall_tir <- sprintf("%.0f%%", 100 * mean(readings$glucose > 70 & readings$glucose < 140))

ui <- page_navbar(
  id           = "main_nav",
  title        = HTML("<span class='brand-main'>CGM Explorer</span>"),
  window_title = "CGM Explorer: continuous glucose monitoring in children",
  theme = bs_theme(version = 5, bootswatch = "flatly",
                   primary = "#1F4E79",
                   base_font = "Arial", heading_font = "Arial", code_font = "Arial"),
  header   = tags$style(HTML(app_css)),
  fillable = FALSE,
  
  # ===========================================================================
  # TAB 1. START HERE
  # ===========================================================================
  nav_panel(
    title = "Start Here", value = "home",
    
    div(class = "hero",
        h1("What can 14 days of glucose data tell us?"),
        p("A continuous glucose monitor (CGM) is a small sensor that reads a child's blood sugar every 5 minutes. This dashboard turns those readings into clear answers about how children with diabetes, prediabetes, and healthy glucose differ."),
        p("Every chart is interactive. Hover over dots for details, and use the filters on the left of each tab."),
        div(class = "color-key",
            div(class = "color-key-title", tags$b("Color key: "), "blue means in range (good). Red means out of range (bad)."),
            lapply(range_levels, range_chip))
    ),
    
    layout_columns(
      fill = FALSE, col_widths = c(3, 3, 3, 3),
      value_box("Children", n_children, theme = "primary"),
      value_box("Days of monitoring", "14", theme = "primary"),
      value_box("Glucose readings", n_readings, theme = "primary"),
      value_box("Time in range", overall_tir, theme = "primary")
    ),
    
    h4("Take the 5-minute tour", class = "mt-3"),
    layout_columns(
      fill = FALSE, col_widths = c(6, 6, 6, 6),
      card(class = "tour-card", card_body(
        h5("Glucose Overview"),
        p("Compare groups of children by time in range, daily glucose curve, and key measures."),
        p(class = "tour-try", "Try: change \"Time of day\" to Night."),
        actionButton("go_overview", "Open Glucose Overview", class = "btn-outline-primary btn-sm"))),
      card(class = "tour-card", card_body(
        h5("Child Explorer"),
        p("Follow one child through all 14 days, with low episodes marked."),
        p(class = "tour-try", "Try: pick a child with type 2 diabetes and turn on \"Mark low episodes\"."),
        actionButton("go_child", "Open Child Explorer", class = "btn-outline-primary btn-sm"))),
      card(class = "tour-card", card_body(
        h5("False Lows"),
        p("See how sensor pressure creates fake low readings, and test whether removing them changes the results."),
        p(class = "tour-try", "Try: turn on \"Track each child\" to connect day and night."),
        actionButton("go_false", "Open False Lows", class = "btn-outline-primary btn-sm"))),
      card(class = "tour-card", card_body(
        h5("Models & Forecasts"),
        p("Test group differences with a statistical model, then predict glucose up to an hour ahead."),
        p(class = "tour-try", "Try: switch off the repeated-readings correction and watch the intervals shrink."),
        actionButton("go_models", "Open Models & Forecasts", class = "btn-outline-primary btn-sm")))
    ),
    
    card(
      card_header("How much time does each group spend in each glucose range?"),
      girafeOutput("home_ranges", width = "100%"),
      card_footer(p("Blue is the target. Every shade of red is a problem: dark red for lows, lighter and stronger reds for highs."))
    ),
    
    p(class = "about-line mt-3",
      "Data: de-identified readings from 100 children aged 6 to 16 (34 with type 2 diabetes, 33 with prediabetes, 33 healthy), one reading every 5 minutes for 14 days. ",
      "Built with R, Shiny, bslib, ggiraph, and mgcv.")
  ),
  
  # ===========================================================================
  # TAB 2. GLUCOSE OVERVIEW
  # ===========================================================================
  nav_panel(
    title = "Glucose Overview", value = "overview",
    
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Filters"),
        checkboxGroupInput("ov_groups", "Group", choices = group_levels, selected = group_levels),
        sliderInput("ov_age", "Age (years)", min = 6, max = 16, value = c(6, 16), step = 1),
        sliderInput("ov_bmi", "Body mass index (kg/m\u00b2)", min = 15, max = 50, value = c(15, 50), step = 1),
        radioButtons("ov_period", "Time of day",
                     choices = c("All day" = "all", "Day (6 AM - 10 PM)" = "Day", "Night (10 PM - 6 AM)" = "Night")),
        sliderInput("ov_days", "Days of monitoring", min = 1, max = 14, value = c(1, 14), step = 1),
        hr(),
        selectInput("ov_measure", "Measure shown for each child",
                    choices = overview_measure_choices)
      ),
      
      tab_intro(
        "How groups of children compare on average glucose, time in range, and the shape of a typical day.",
        "Children with type 2 diabetes spend the least time in range, and their glucose swings the most."
      ),
      
      layout_columns(
        fill = FALSE, col_widths = c(2, 3, 2, 2, 3),
        value_box("Children", textOutput("ov_n"), theme = "primary"),
        value_box("Average glucose", textOutput("ov_mean"), theme = "primary"),
        value_box("Estimated HbA1c", textOutput("ov_gmi"), theme = "primary"),
        value_box("Glucose variation", textOutput("ov_cv"), theme = "primary"),
        value_box("Time in range", textOutput("ov_tir"), theme = "primary")
      ),
      
      card(
        card_header("A typical day, by group"),
        girafeOutput("ov_day", width = "100%"),
        card_footer(p("Line = median glucose at each time of day. Dark band = middle half of readings. Light band = 90% of readings. Dashed lines = 70 and 140 mg/dL."))
      ),
      
      layout_columns(
        fill = FALSE, col_widths = c(6, 6),
        card(
          card_header("Time in each glucose range"),
          girafeOutput("ov_ranges", width = "100%"),
          card_footer(p("Average share of readings in each range."))
        ),
        card(
          card_header("How much do children differ?"),
          girafeOutput("ov_dist", width = "100%"),
          card_footer(p("One dot per child. Hover to see who it is."))
        )
      ),
      
      card(
        card_header("Summary by group"),
        tableOutput("ov_table"),
        card_footer(p("Group averages for the children and days you selected."))
      )
    )
  ),
  
  # ===========================================================================
  # TAB 3. CHILD EXPLORER
  # ===========================================================================
  nav_panel(
    title = "Child Explorer", value = "child",
    
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Choose a child"),
        selectInput("ce_pid", "Participant", choices = participant_choices,
                    selected = participants$participant_id[1]),
        sliderInput("ce_days", "Days shown", min = 1, max = 14, value = c(1, 14), step = 1),
        hr(),
        h6("Glucose ranges to show"),
        checkboxGroupInput("ce_ranges", NULL,
                           choiceNames = range_choice_names, choiceValues = range_levels,
                           selected = range_levels),
        input_switch("ce_episodes", "Mark low episodes", value = FALSE),
        helpText("An episode is 15 minutes or more at 70 mg/dL or below."),
        hr(),
        input_switch("ce_clean", "Clean false lows", value = FALSE),
        conditionalPanel(
          condition = "input.ce_clean",
          sliderInput("ce_floor", "Raise readings below (mg/dL)", min = 50, max = 70, value = 60, step = 1)
        ),
        helpText("Raises every reading below the cutoff up to the cutoff. Sharp dips flatten; the rest of the trace stays the same."),
        hr(),
        selectInput("ce_heat", "Day-by-hour map shows",
                    choices = c("Average glucose (mg/dL)" = "mean_glu",
                                "Time low (%)"            = "pct_low",
                                "Time in range (%)"       = "pct_inr"))
      ),
      
      tab_intro(
        "One child's glucose for up to 14 days, colored by range.",
        "Glucose follows a daily rhythm. Sharp overnight dips are often false lows from sensor pressure."
      ),
      
      layout_columns(
        fill = FALSE, col_widths = c(3, 3, 3, 3),
        value_box("Average glucose", textOutput("ce_mean"), theme = "primary"),
        value_box("Time in range", textOutput("ce_tir"), theme = "primary"),
        value_box("Low episodes", textOutput("ce_eps"), theme = "primary"),
        value_box("Glucose variation", textOutput("ce_cv"), theme = "primary")
      ),
      
      card(
        card_header("Glucose over time"),
        plotOutput("ce_trace", height = "430px"),
        card_footer(p("Each line segment takes the color of its range. Turn ranges on or off in the sidebar."))
      ),
      
      card(
        card_header("Every day, every hour"),
        girafeOutput("ce_heatmap", width = "100%"),
        card_footer(p("Each tile is one hour of one day. Hover for the value. Days run left to right; hours run top to bottom."))
      )
    )
  ),
  
  # ===========================================================================
  # TAB 4. FALSE LOWS
  # ===========================================================================
  nav_panel(
    title = "False Lows", value = "false",
    
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Filters"),
        checkboxGroupInput("fl_groups", "Group", choices = group_levels, selected = group_levels),
        sliderInput("fl_age", "Age (years)", min = 6, max = 16, value = c(6, 16), step = 1),
        sliderInput("fl_bmi", "Body mass index (kg/m\u00b2)", min = 15, max = 50, value = c(15, 50), step = 1),
        checkboxInput("fl_track", "Track each child", value = FALSE),
        helpText("Connects each child's day and night dots. Hover over a dot to highlight that child."),
        hr(),
        h5("Cleaning test"),
        selectInput("fl_floor", "Raise readings below (mg/dL)",
                    choices = c(50, 55, 60, 65, 70), selected = 60),
        sliderInput("fl_margin", "Largest change still called negligible (\u00b1 %)",
                    min = 1, max = 10, value = 5, step = 0.5),
        helpText("A smaller number is a stricter test."),
        input_switch("fl_split", "Split results by group", value = TRUE)
      ),
      
      tab_intro(
        "A false low is a sudden, fake drop in the sensor reading. It happens when a sleeping child lies on the sensor. False lows trigger needless alarms and wrong treatment decisions.",
        "False lows happen far more often at night. Cleaning them barely moves average glucose but can shift the risk measures."
      ),
      
      layout_columns(
        fill = FALSE, col_widths = c(3, 3, 3, 3),
        value_box("Children", textOutput("fl_n"), theme = "primary"),
        value_box("False lows by day", textOutput("fl_day"), theme = "primary"),
        value_box("False lows by night", textOutput("fl_night"), theme = "primary"),
        value_box("Share at night", textOutput("fl_pct"), theme = "primary")
      ),
      
      card(
        card_header("False lows per child: day vs. night"),
        girafeOutput("fl_box", width = "100%"),
        card_footer(p("Each dot is one child's false lows over 14 days. Boxes show the middle half of children."))
      ),
      
      card(
        card_header("Does cleaning false lows change the results?"),
        card_body(
          uiOutput("fl_verdict"),
          p("We recalculate four summary measures after raising low readings to the cutoff:"),
          tags$ul(
            tags$li(tags$b("HbA1c estimate (GMI): "), "estimated average blood sugar over several months."),
            tags$li(tags$b("Big swings (MAGE): "), "typical size of the large glucose swings."),
            tags$li(tags$b("Level + swings (J-index): "), "combines how high and how variable glucose is."),
            tags$li(tags$b("Risk of lows and highs (ADRR): "), "daily risk from dangerous lows and highs.")
          ),
          accordion(open = FALSE, class = "mt-2",
                    accordion_panel("Technical details",
                                    tags$ul(
                                      tags$li("GMI = 3.31 + 0.02392 \u00d7 mean glucose (mg/dL)."),
                                      tags$li("MAGE = mean of the glucose swings larger than 1 standard deviation."),
                                      tags$li("J-index = 0.001 \u00d7 (mean + SD)\u00b2."),
                                      tags$li("ADRR = largest daily low risk plus largest daily high risk, averaged over 14 days."),
                                      tags$li("% change = 100 \u00d7 (cleaned \u2212 original) / original, per child."),
                                      tags$li("\"Unchanged\" means the full 95% confidence interval of the average % change lies inside the \u00b1 margin (equivalence test).")
                                    )))
        )
      ),
      
      layout_columns(
        fill = FALSE, col_widths = c(6, 6),
        card(
          card_header("Average change per measure"),
          girafeOutput("fl_forest", width = "100%"),
          card_footer(p("Dot = average % change, line = 95% interval. Inside the gray band means unchanged."))
        ),
        card(
          card_header("Does any child change a lot?"),
          girafeOutput("fl_scatter", width = "100%"),
          card_footer(p("One dot per child. A dot on the diagonal did not change."))
        )
      )
    )
  ),
  
  # ===========================================================================
  # TAB 5. MODELS & FORECASTS
  # ===========================================================================
  nav_panel(
    title = "Models & Forecasts", value = "models",
    
    navset_card_underline(
      # ---- 5a. Compare groups over time ------------------------------------
      nav_panel(
        "Compare groups",
        layout_sidebar(
          sidebar = sidebar(
            width = 320,
            h5("Build the model"),
            selectInput("md_outcome", "What to compare", choices = model_outcomes),
            checkboxGroupInput("md_terms", "What might explain it",
                               choices = model_predictors,
                               selected = c("group", "night", "weekend")),
            helpText("Every child gets their own baseline. Choosing Group \u00d7 Night adds Group and Night automatically."),
            input_switch("md_ar", "Account for repeated readings from each child", value = TRUE),
            helpText("Readings from the same child resemble each other. Turn this off to see how ignoring that makes results look more certain than they are."),
            actionButton("md_fit", "Fit model", class = "btn-primary w-100")
          ),
          
          tab_intro(
            "Which factors change glucose, after accounting for the fact that each child is measured many times.",
            "A dark dot with a line that avoids zero marks a real difference. Compare the results with and without the repeated-readings correction."
          ),
          
          card(card_header("What was fit?"), card_body(uiOutput("md_summary"))),
          
          layout_columns(
            fill = FALSE, col_widths = c(6, 6),
            card(
              card_header("Effects with 95% confidence intervals"),
              girafeOutput("md_forest", width = "100%"),
              card_footer(p("Dot right of zero = higher than the comparison. Group effects compare with type 2 diabetes; night compares with day."))
            ),
            card(
              card_header("The 24-hour curve, by group"),
              plotOutput("md_gam", height = "400px"),
              card_footer(p("Smooth curves of average glucose over the day, with 95% bands."))
            )
          ),
          
          card(accordion(open = FALSE,
                         accordion_panel("Full model output (technical)", verbatimTextOutput("md_text"))))
        )
      ),
      
      # ---- 5b. Forecast glucose ---------------------------------------------
      nav_panel(
        "Forecast glucose",
        layout_sidebar(
          sidebar = sidebar(
            width = 320,
            h5("Settings"),
            selectInput("fc_pid", "Participant", choices = participant_choices,
                        selected = participants$participant_id[1]),
            selectInput("fc_h", "How far ahead",
                        choices = c("15 minutes" = 3, "30 minutes" = 6, "45 minutes" = 9, "60 minutes" = 12),
                        selected = 6),
            selectInput("fc_day", "Day shown in the chart", choices = 9:14, selected = 11),
            helpText("Days 9 - 14 hold the forecasts. Days 1 - 8 train the models."),
            hr(),
            input_switch("fc_clean", "Clean false lows first", value = FALSE),
            conditionalPanel(
              condition = "input.fc_clean",
              sliderInput("fc_floor", "Raise readings below (mg/dL)", min = 50, max = 70, value = 60, step = 1)
            ),
            helpText("Forecasts are scored against the same cleaned readings.")
          ),
          
          tab_intro(
            "Three methods predict glucose a short time ahead using only the readings available at that moment. Each method trains on earlier days and predicts later days.",
            "The regression model beats the simple methods. Sudden false lows are hard to predict for any method."
          ),
          
          card(
            card_header("Forecast vs. what happened"),
            plotOutput("fc_trace", height = "420px"),
            card_footer(p("Each forecast appears at the time it predicts. Blue band = 95% prediction range for the regression model. Dashed lines = 70 and 140 mg/dL."))
          ),
          
          layout_columns(
            fill = FALSE, col_widths = c(6, 6),
            card(
              card_header("How accurate are the forecasts?"),
              tableOutput("fc_table"),
              card_footer(p("Lower error is better. \"Close enough\" means within 15 mg/dL (at 100 or below) or within 15% (above 100)."))
            ),
            card(
              card_header("Can a forecast warn of a coming low?"),
              plotOutput("fc_roc", height = "380px"),
              card_footer(p("A curve closer to the top-left corner warns better. AUC of 0.5 equals a coin flip; 1.0 is perfect."))
            )
          ),
          
          card(
            card_header("All 100 children: regression model vs. repeating the latest reading"),
            card_body(
              actionButton("fc_all", "Run all 100 children", class = "btn-primary"),
              helpText("Uses the settings above. Takes a few seconds."),
              textOutput("fc_all_text"),
              girafeOutput("fc_all_plot", width = "100%")
            ),
            card_footer(p("One dot per child. Dots below the diagonal mean the regression model made smaller errors."))
          )
        )
      )
    )
  )
)


# =============================================================================
# SECTION 7. SERVER
# =============================================================================

server <- function(input, output, session) {
  
  # ---- Tab 1: navigation buttons and overview chart --------------------------
  observeEvent(input$go_overview, nav_select("main_nav", "overview"))
  observeEvent(input$go_child,    nav_select("main_nav", "child"))
  observeEvent(input$go_false,    nav_select("main_nav", "false"))
  observeEvent(input$go_models,   nav_select("main_nav", "models"))
  
  output$home_ranges <- renderGirafe(plot_range_bars(child_summary(readings)))
  
  # ===========================================================================
  # Tab 2: Glucose Overview
  # ===========================================================================
  
  # Readings filtered by group, age, BMI, and days (used for the daily curve)
  ov_base <- reactive({
    req(input$ov_groups, input$ov_days)
    readings %>%
      filter(group %in% input$ov_groups,
             age >= input$ov_age[1], age <= input$ov_age[2],
             bmi >= input$ov_bmi[1], bmi <= input$ov_bmi[2],
             day_num >= input$ov_days[1], day_num <= input$ov_days[2])
  })
  
  # Same, plus the time-of-day choice (used for everything else)
  ov_data <- reactive({
    d <- ov_base()
    if (input$ov_period != "all") d <- filter(d, period == input$ov_period)
    validate(need(nrow(d) > 0, "No children match these filters. Widen the filters to see results."))
    d
  })
  ov_child <- reactive(child_summary(ov_data()))
  
  output$ov_n    <- renderText(n_distinct(ov_data()$participant_id))
  output$ov_mean <- renderText(sprintf("%.0f mg/dL", mean(ov_child()$mean_glu)))
  output$ov_gmi  <- renderText(sprintf("%.1f%%", mean(ov_child()$gmi)))
  output$ov_cv   <- renderText(sprintf("%.1f%%", mean(ov_child()$cv)))
  output$ov_tir  <- renderText(sprintf("%.0f%%", mean(ov_child()$pct_inr)))
  
  output$ov_day <- renderGirafe({
    d <- ov_base()
    validate(need(nrow(d) > 0, "No children match these filters. Widen the filters to see results."))
    a <- d %>%
      group_by(group, bin) %>%
      summarise(q05 = quantile(glucose, 0.05), q25 = quantile(glucose, 0.25),
                q50 = median(glucose),
                q75 = quantile(glucose, 0.75), q95 = quantile(glucose, 0.95),
                .groups = "drop") %>%
      mutate(tip = paste0("<b>", group, "</b><br/>", format_clock(bin), "<br/>",
                          "Median: ", sprintf("%.0f", q50), " mg/dL<br/>",
                          "Middle half: ", sprintf("%.0f", q25), " to ", sprintf("%.0f", q75), "<br/>",
                          "90% of readings: ", sprintf("%.0f", q05), " to ", sprintf("%.0f", q95)))
    
    p <- ggplot(a, aes(x = bin)) +
      geom_hline(yintercept = c(70, 140), linetype = "dashed", color = "grey60") +
      geom_ribbon(aes(ymin = q05, ymax = q95, fill = group), alpha = 0.18) +
      geom_ribbon(aes(ymin = q25, ymax = q75, fill = group), alpha = 0.35) +
      geom_line(aes(y = q50, color = group), linewidth = 1.1) +
      geom_point_interactive(aes(y = q50, color = group, tooltip = tip,
                                 data_id = paste(group, bin)), size = 1.8) +
      facet_wrap(~ group, nrow = 1) +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors) +
      scale_x_continuous(breaks = clock_breaks$breaks, labels = clock_breaks$labels) +
      coord_cartesian(ylim = c(40, 250)) +
      labs(x = "Time of day", y = "Sensor glucose (mg/dL)") +
      theme_dashboard()
    make_girafe(p, 12, 5, fade = 0.35)
  })
  
  output$ov_ranges <- renderGirafe(plot_range_bars(ov_child()))
  
  output$ov_dist <- renderGirafe({
    cs  <- ov_child()
    met <- input$ov_measure
    lab <- names(overview_measure_choices)[overview_measure_choices == met]
    plot_box_by_group(cs, cs[[met]], lab)
  })
  
  output$ov_table <- renderTable({
    ov_child() %>%
      group_by(Group = as.character(group)) %>%
      summarise(Children = n(),
                `Average glucose (mg/dL)` = mean(mean_glu),
                `HbA1c estimate (%)`      = mean(gmi),
                `Variation, CV (%)`       = mean(cv),
                `Time in range (%)`       = mean(pct_inr),
                `Time low (%)`            = mean(pct_low),
                `Time below 54 (%)`       = mean(pct_low54),
                `Time high (%)`           = mean(pct_above),
                .groups = "drop")
  }, digits = 1, striped = TRUE, bordered = TRUE, width = "100%")
  
  # ===========================================================================
  # Tab 3: Child Explorer
  # ===========================================================================
  
  # The selected child's readings, cleaned if the switch is on
  ce_glucose <- reactive({
    req(input$ce_pid)
    g <- all_traces[[input$ce_pid]]$glucose
    clean_readings(g, if (isTRUE(input$ce_clean)) input$ce_floor else NULL)
  })
  
  output$ce_mean <- renderText(sprintf("%.0f mg/dL", mean(ce_glucose())))
  output$ce_tir  <- renderText(sprintf("%.0f%%", 100 * mean(ce_glucose() > 70 & ce_glucose() < 140)))
  output$ce_cv   <- renderText(sprintf("%.0f%%", 100 * sd(ce_glucose()) / mean(ce_glucose())))
  output$ce_eps  <- renderText({
    e <- detect_episodes(ce_glucose())
    if (is.null(e)) "0" else nrow(e)
  })
  
  output$ce_trace <- renderPlot({
    req(input$ce_days)
    validate(need(length(input$ce_ranges) > 0, "Turn on at least one glucose range to see the trace."))
    tr <- all_traces[[input$ce_pid]]
    g  <- ce_glucose()
    d0 <- input$ce_days[1] - 1
    d1 <- input$ce_days[2]
    
    df <- tibble(x = tr$day, glucose = g) %>%
      filter(x >= d0, x < d1) %>%
      mutate(range = classify_glucose(glucose),
             # Hidden ranges become NA so the line breaks instead of bridging them
             glucose_shown = ifelse(range %in% input$ce_ranges, glucose, NA_real_))
    
    p <- ggplot(df, aes(x = x, y = glucose_shown))
    if (isTRUE(input$ce_episodes)) {
      e <- detect_episodes(g)
      if (!is.null(e)) {
        e <- e %>% mutate(xmin = (start - 1) * 5 / 1440, xmax = end * 5 / 1440) %>%
          filter(xmax > d0, xmin < d1)
        if (nrow(e) > 0) {
          p <- p + geom_rect(data = e, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                             inherit.aes = FALSE, fill = range_colors[["Low"]], alpha = 0.22)
        }
      }
    }
    p +
      geom_vline(xintercept = d0:d1, color = "grey88") +
      geom_hline(yintercept = c(70, 140), linetype = "dashed", color = "grey60") +
      geom_line(aes(color = range, group = 1), linewidth = 0.7) +
      scale_color_manual(values = range_colors, breaks = range_levels, drop = FALSE) +
      scale_x_continuous(breaks = d0:(d1 - 1) + 0.5, labels = (d0 + 1):d1) +
      coord_cartesian(xlim = c(d0, d1), ylim = c(40, 250), expand = FALSE) +
      guides(color = guide_legend(override.aes = list(linewidth = 2))) +
      labs(x = "Day of monitoring", y = "Sensor glucose (mg/dL)") +
      theme_dashboard() +
      theme(legend.position = "top", legend.title = element_blank(),
            legend.text = element_text(size = base_font_px))
  })
  
  output$ce_heatmap <- renderGirafe({
    met <- input$ce_heat
    lab <- c(mean_glu = "Average glucose (mg/dL)", pct_low = "Time low (%)",
             pct_inr = "Time in range (%)")[[met]]
    g <- ce_glucose()
    h <- tibble(glucose = g,
                day_num = floor(all_traces[[input$ce_pid]]$day) + 1,
                hour    = (((seq_along(g) - 1) * 5) %% 1440) %/% 60) %>%
      group_by(day_num, hour) %>%
      summarise(mean_glu = mean(glucose),
                pct_low  = 100 * mean(glucose <= 70),
                pct_inr  = 100 * mean(glucose > 70 & glucose < 140),
                .groups = "drop")
    h$value <- h[[met]]
    h$tip <- paste0("<b>Day ", h$day_num, " (", weekday_names[(h$day_num - 1) %% 7 + 1], ")</b><br/>",
                    format_clock(h$hour), "<br/>", lab, ": ", sprintf("%.1f", h$value))
    
    p <- ggplot(h, aes(x = day_num, y = hour, fill = value)) +
      geom_tile_interactive(aes(tooltip = tip, data_id = paste(day_num, hour)),
                            color = "white", linewidth = 0.2) +
      scale_y_reverse(breaks = c(0, 6, 12, 18), labels = c("12 AM", "6 AM", "12 PM", "6 PM")) +
      scale_x_continuous(breaks = 1:14,
                         labels = paste0(1:14, "\n", weekday_names[((1:14) - 1) %% 7 + 1]),
                         expand = c(0, 0)) +
      labs(x = "Day of monitoring", y = "Time of day") +
      theme_dashboard() + theme(legend.position = "right")
    
    p <- switch(met,
                mean_glu = p + glucose_fill_scale(lab),
                pct_low  = p + scale_fill_gradient(name = lab, low = "white", high = range_colors[["Low"]]),
                pct_inr  = p + scale_fill_gradient(name = lab, low = "white", high = range_colors[["In range"]]))
    make_girafe(p, 12, 5.5, fade = 0.4)
  })
  
  # ===========================================================================
  # Tab 4: False Lows
  # ===========================================================================
  
  fl_filtered <- reactive({
    req(input$fl_groups)
    d <- period_summary %>%
      filter(group %in% input$fl_groups,
             age >= input$fl_age[1], age <= input$fl_age[2],
             bmi >= input$fl_bmi[1], bmi <= input$fl_bmi[2])
    validate(need(nrow(d) > 0, "No children match these filters. Widen the filters to see results."))
    d
  })
  
  output$fl_n     <- renderText(n_distinct(fl_filtered()$participant_id))
  output$fl_day   <- renderText(sum(fl_filtered()$comp_events[fl_filtered()$period == "Day"]))
  output$fl_night <- renderText(sum(fl_filtered()$comp_events[fl_filtered()$period == "Night"]))
  output$fl_pct   <- renderText({
    d <- fl_filtered()
    paste0(round(100 * sum(d$comp_events[d$period == "Night"]) / sum(d$comp_events), 1), "%")
  })
  
  output$fl_box <- renderGirafe({
    df <- fl_filtered() %>%
      mutate(value    = comp_events,
             x_center = as.numeric(period),
             x_point  = x_center + jitter_x,
             tip      = paste0("<b>", participant_id, "</b><br/>Age: ", age,
                               " years<br/>Body mass index: ", bmi, "<br/>",
                               period, " false lows: ", comp_events))
    
    p <- ggplot(df) +
      geom_boxplot(aes(x = x_center, y = value, group = period, fill = period, color = period),
                   alpha = 0.25, outlier.shape = NA, width = 0.55, linewidth = 0.8)
    if (isTRUE(input$fl_track)) {
      p <- p + geom_line_interactive(
        aes(x = x_point, y = value, group = participant_id, data_id = participant_id, tooltip = tip),
        color = "grey70", linewidth = 0.5)
    }
    p <- p +
      geom_point_interactive(aes(x = x_point, y = value, color = period,
                                 data_id = participant_id, tooltip = tip),
                             alpha = 0.85, size = 3) +
      facet_wrap(~ group) +
      scale_fill_manual(values = period_colors) +
      scale_color_manual(values = period_colors) +
      scale_x_continuous(breaks = c(1, 2), labels = c("Day", "Night"), limits = c(0.5, 2.5)) +
      labs(x = NULL, y = "False lows per 14 days") +
      theme_dashboard()
    make_girafe(p, 12, 5.2, fade = 0.12)
  })
  
  # Original vs. cleaned measures for the selected groups
  fl_changes <- reactive({
    req(input$fl_groups, input$fl_floor)
    cleaned <- measures_for(as.numeric(input$fl_floor)) %>% rename(cleaned = value)
    measures_for(NULL) %>% rename(original = value) %>%
      left_join(cleaned, by = c("participant_id", "measure")) %>%
      left_join(participants %>% select(participant_id, group), by = "participant_id") %>%
      filter(group %in% input$fl_groups, !is.na(original), !is.na(cleaned), original > 0) %>%
      mutate(measure = factor(measure, levels = summary_measure_names),
             pct     = 100 * (cleaned - original) / original)
  })
  
  output$fl_verdict <- renderUI({
    s <- summarise_change(fl_changes(), input$fl_margin, by_group = FALSE)
    validate(need(nrow(s) > 0, "Select at least one group to see the answer."))
    n_ok  <- sum(s$unchanged)
    cut   <- input$fl_floor
    headline <- if (n_ok == nrow(s)) {
      sprintf("Cleaning at %s mg/dL leaves all %d measures essentially unchanged.", cut, nrow(s))
    } else if (n_ok == 0) {
      sprintf("Cleaning at %s mg/dL changes all %d measures.", cut, nrow(s))
    } else {
      sprintf("Cleaning at %s mg/dL leaves %d of %d measures essentially unchanged.", cut, n_ok, nrow(s))
    }
    badges <- lapply(seq_len(nrow(s)), function(i) {
      span(class = paste("badge-pill", if (s$unchanged[i]) "badge-good" else "badge-bad"),
           paste0(if (s$unchanged[i]) "\u2714 " else "\u26A0 ", as.character(s$measure[i]),
                  ": ", sprintf("%+.1f%%", s$pct_mean[i])))
    })
    tagList(
      div(class = "verdict-headline", headline),
      div(class = "badge-row", badges),
      div(class = "text-muted",
          sprintf("Blue = unchanged. Red = changed. All selected children combined; margin \u00b1%s%%.", input$fl_margin))
    )
  })
  
  output$fl_forest <- renderGirafe({
    m     <- input$fl_margin
    split <- isTRUE(input$fl_split)
    s <- summarise_change(fl_changes(), m, by_group = split)
    req(nrow(s) > 0)
    s <- s %>%
      mutate(tip = paste0("<b>", measure, "</b><br/>",
                          if (split) paste0(grp, "<br/>") else "",
                          "Average change: ", sprintf("%.2f", pct_mean), "%<br/>",
                          "95% interval: ", sprintf("%.2f", lo), "% to ", sprintf("%.2f", hi), "%<br/>",
                          "Children: ", n, "<br/>",
                          ifelse(unchanged, "Unchanged (inside the margin)", "Changed")),
             measure = factor(measure, levels = rev(summary_measure_names)))
    x_lim <- max(m, abs(s$lo), abs(s$hi), na.rm = TRUE) * 1.15
    
    p <- ggplot(s, aes(x = pct_mean, y = measure, xmin = lo, xmax = hi, color = grp)) +
      annotate("rect", xmin = -m, xmax = m, ymin = -Inf, ymax = Inf, fill = "#2B7BBA", alpha = 0.12) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      geom_pointrange_interactive(aes(tooltip = tip, data_id = paste(measure, grp)),
                                  position = position_dodge(width = 0.6),
                                  linewidth = 0.9, size = 0.8) +
      scale_color_manual(values = c(group_colors, "All selected groups" = "grey20")) +
      scale_y_discrete(labels = scales::label_wrap(18)) +
      coord_cartesian(xlim = c(-x_lim, x_lim)) +
      labs(x = "Average % change after cleaning (95% interval)", y = NULL,
           caption = paste0("Blue band = negligible change (\u00b1", m, "%)")) +
      theme_dashboard() +
      (if (split) legend_top(2) else theme())
    make_girafe(p, 7, 6, fade = 0.3)
  })
  
  output$fl_scatter <- renderGirafe({
    d <- fl_changes()
    req(nrow(d) > 0)
    m <- input$fl_margin
    
    st <- d %>%
      group_by(measure) %>%
      summarise(pct_over = 100 * mean(abs(pct) > m),
                lo = min(original, cleaned), hi = max(original, cleaned), .groups = "drop") %>%
      arrange(measure) %>%
      mutate(panel = paste0(measure, "\n", sprintf("%.0f", pct_over), "% of children change by more than ", m, "%"))
    panel_levels <- st$panel
    
    d <- d %>%
      left_join(st %>% select(measure, panel), by = "measure") %>%
      mutate(panel = factor(panel, levels = panel_levels),
             tip = paste0("<b>", participant_id, "</b><br/>", group, "<br/>", measure, "<br/>",
                          "Original: ", sprintf("%.2f", original), "<br/>",
                          "Cleaned: ", sprintf("%.2f", cleaned), "<br/>",
                          "Change: ", sprintf("%+.1f", pct), "%"))
    
    # Same range on both axes in every panel, so the diagonal means "no change"
    lim <- bind_rows(
      st %>% transmute(panel = factor(panel, levels = panel_levels), x = lo, y = lo),
      st %>% transmute(panel = factor(panel, levels = panel_levels), x = hi, y = hi))
    
    p <- ggplot(d, aes(x = original, y = cleaned)) +
      geom_blank(data = lim, aes(x = x, y = y), inherit.aes = FALSE) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      geom_point_interactive(aes(color = group, tooltip = tip, data_id = participant_id),
                             alpha = 0.85, size = 3) +
      facet_wrap(~ panel, scales = "free", ncol = 2) +
      scale_color_manual(values = group_colors) +
      labs(x = "Original value", y = "Value after cleaning") +
      theme_dashboard() + legend_top(2) +
      theme(strip.text = element_text(size = base_font_px - 2, face = "bold"))
    make_girafe(p, 7, 6, fade = 0.12)
  })
  
  # ===========================================================================
  # Tab 5a: Compare groups (mixed-effects model)
  # ===========================================================================
  
  # Fits on click (and once on first view). Child x day x period averages keep it fast.
  md_model <- eventReactive(input$md_fit, {
    validate(need(length(input$md_terms) > 0, "Choose at least one factor, then select Fit model."))
    chosen <- input$md_terms
    if ("grp_night" %in% chosen) chosen <- union(chosen, c("group", "night"))
    
    term_map <- c(group = "group", night = "night", weekend = "weekend",
                  day_num = "day_num", age = "age", bmi = "bmi", sex = "sex")
    rhs <- unname(term_map[intersect(names(term_map), chosen)])
    if ("grp_night" %in% input$md_terms) rhs <- c(rhs, "group:night")
    f <- as.formula(paste(input$md_outcome, "~", paste(rhs, collapse = " + ")))
    use_ar <- isTRUE(input$md_ar)
    
    fit <- withProgress(message = "Fitting model...", value = 0.3, {
      tryCatch(
        if (use_ar) {
          # Random intercept per child + correlation between neighboring periods
          nlme::lme(f, random = ~ 1 | participant_id,
                    correlation = nlme::corAR1(form = ~ time_idx | participant_id),
                    data = model_table, method = "REML")
        } else {
          nlme::lme(f, random = ~ 1 | participant_id, data = model_table, method = "REML")
        },
        error = function(e) e)
    })
    if (inherits(fit, "error")) {
      validate(paste("The model could not be fit:", conditionMessage(fit)))
    }
    fit$call$fixed <- f
    list(fit = fit, ar = use_ar, outcome = input$md_outcome)
  }, ignoreNULL = FALSE)
  
  output$md_summary <- renderUI({
    m   <- md_model()
    vc  <- nlme::VarCorr(m$fit)
    sds <- suppressWarnings(as.numeric(vc[, "StdDev"]))
    phi <- if (m$ar) as.numeric(coef(m$fit$modelStruct$corStruct, unconstrained = FALSE)) else NA_real_
    out_lab <- names(model_outcomes)[model_outcomes == m$outcome]
    badge <- function(txt) span(class = "badge-pill badge-info", txt)
    
    tagList(
      p(sprintf("The model compares %s using %s averages from %d children (each child has one day average and one night average for 14 days). Each child gets a personal baseline%s.",
                tolower(out_lab), format(nrow(model_table), big.mark = ","),
                n_distinct(model_table$participant_id),
                if (m$ar) ", and the model accounts for neighboring periods resembling each other" else "")),
      if (!m$ar) div(class = "text-danger",
                     "Warning: this version ignores that readings from one child resemble each other. Its intervals are too narrow and its results look more certain than they are."),
      div(class = "badge-row",
          badge(sprintf("Differences between children (SD): %.1f", sds[1])),
          badge(sprintf("Differences within a child (SD): %.1f", sds[2])),
          if (m$ar) badge(sprintf("Similarity of neighboring periods: %.2f", phi)))
    )
  })
  
  output$md_forest <- renderGirafe({
    m  <- md_model()
    tt <- summary(m$fit)$tTable
    d <- tibble(term = rownames(tt), est = tt[, "Value"], se = tt[, "Std.Error"],
                dfr = tt[, "DF"], pval = tt[, "p-value"]) %>%
      filter(term != "(Intercept)") %>%
      mutate(label = clean_term_label(term),
             crit  = qt(0.975, dfr),
             lo    = est - crit * se,
             hi    = est + crit * se,
             real  = lo > 0 | hi < 0,
             tip   = paste0("<b>", label, "</b><br/>Estimate: ", sprintf("%.2f", est), "<br/>",
                            "95% interval: ", sprintf("%.2f", lo), " to ", sprintf("%.2f", hi), "<br/>",
                            "p = ", format.pval(pval, digits = 2, eps = 0.001)))
    d$label <- factor(d$label, levels = rev(d$label))
    
    p <- ggplot(d, aes(x = est, y = label, xmin = lo, xmax = hi, color = real)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      geom_pointrange_interactive(aes(tooltip = tip, data_id = label), linewidth = 0.9, size = 0.8) +
      scale_color_manual(values = c("TRUE" = navy, "FALSE" = "grey60")) +
      scale_y_discrete(labels = scales::label_wrap(22)) +
      labs(x = "Estimated difference (95% interval)", y = NULL) +
      theme_dashboard()
    make_girafe(p, 7, 5.5, fade = 0.3)
  })
  
  # 24-hour curve per group: generalized additive model with a cyclic smooth
  # (hour 0 joins hour 24) and a random intercept per child. Fit once, on demand.
  gam_curves <- reactive({
    withProgress(message = "Fitting 24-hour curves...", {
      hourly <- readings %>%
        group_by(participant_id, group, hour) %>%
        summarise(glu = mean(glucose), .groups = "drop") %>%
        mutate(hour = hour + 0.5, participant_id = factor(participant_id)) %>%
        as.data.frame()
      fit <- mgcv::gam(glu ~ group + s(hour, by = group, bs = "cc", k = 8) + s(participant_id, bs = "re"),
                       data = hourly, method = "REML", knots = list(hour = c(0, 24)))
      new <- expand.grid(hour = seq(0, 24, length.out = 97),
                         group = factor(group_levels, levels = group_levels))
      new$participant_id <- hourly$participant_id[1]
      pr <- predict(fit, new, se.fit = TRUE, exclude = "s(participant_id)")
      new$fit <- as.numeric(pr$fit)
      new$lo  <- new$fit - 1.96 * as.numeric(pr$se.fit)
      new$hi  <- new$fit + 1.96 * as.numeric(pr$se.fit)
      new
    })
  })
  
  output$md_gam <- renderPlot({
    ggplot(gam_curves(), aes(x = hour, y = fit, color = group, fill = group)) +
      geom_hline(yintercept = c(70, 140), linetype = "dashed", color = "grey60") +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors) +
      scale_x_continuous(breaks = clock_breaks$breaks, labels = clock_breaks$labels) +
      labs(x = "Time of day", y = "Average glucose (mg/dL)") +
      theme_dashboard() + legend_top(2)
  })
  
  output$md_text <- renderPrint(summary(md_model()$fit))
  
  # ===========================================================================
  # Tab 5b: Forecast glucose
  # ===========================================================================
  
  fc_floor_val <- reactive(if (isTRUE(input$fc_clean)) input$fc_floor else NULL)
  fc_horizon   <- reactive(as.integer(input$fc_h))
  
  fc_preds <- reactive({
    req(input$fc_pid, input$fc_h)
    run_forecast(all_traces[[input$fc_pid]], fc_horizon(), fc_floor_val())
  })
  
  output$fc_trace <- renderPlot({
    p  <- fc_preds()
    dd <- as.integer(input$fc_day)
    p  <- p[p$day == dd, ]
    validate(need(nrow(p) > 0, "No forecasts exist for this day."))
    
    all_names <- c("Actual", unname(forecast_methods))
    lines <- bind_rows(
      tibble(xh = p$xh, pred = p$actual, method = "Actual"),
      forecasts_long(p) %>% select(xh, pred, method)
    ) %>% mutate(method = factor(method, levels = all_names))
    
    ggplot() +
      geom_hline(yintercept = c(70, 140), linetype = "dashed", color = "grey60") +
      geom_ribbon(data = p, aes(x = xh, ymin = ar - 1.96 * sigma, ymax = ar + 1.96 * sigma),
                  fill = forecast_colors[["Regression model"]], alpha = 0.15) +
      geom_line(data = lines, aes(x = xh, y = pred, color = method, linewidth = method == "Actual")) +
      scale_color_manual(values = forecast_colors, breaks = all_names) +
      scale_linewidth_manual(values = c("TRUE" = 1.2, "FALSE" = 0.8), guide = "none") +
      scale_x_continuous(breaks = clock_breaks$breaks, labels = clock_breaks$labels) +
      coord_cartesian(xlim = c(0, 24), ylim = c(40, 250)) +
      labs(x = paste0("Time of day (day ", dd, ")"), y = "Sensor glucose (mg/dL)") +
      theme_dashboard() + legend_top(1) +
      guides(color = guide_legend(override.aes = list(linewidth = 2)))
  })
  
  output$fc_table <- renderTable(forecast_accuracy(forecasts_long(fc_preds())),
                                 digits = 1, striped = TRUE, bordered = TRUE, width = "100%")
  
  output$fc_roc <- renderPlot({
    p <- fc_preds()
    y <- p$actual <= 70
    validate(need(sum(y) >= 5 && sum(!y) >= 5,
                  "This child has too few low readings in the test days to draw a warning curve."))
    curves <- lapply(names(forecast_methods), function(m) roc_curve(-p[[m]], y))
    labels <- vapply(seq_along(curves), function(i)
      sprintf("%s (AUC %.2f)", forecast_methods[[i]], curves[[i]]$auc), character(1))
    df <- bind_rows(lapply(seq_along(curves), function(i) mutate(curves[[i]]$df, method = labels[i]))) %>%
      mutate(method = factor(method, levels = labels))
    
    ggplot(df, aes(x = fpr, y = tpr, color = method)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = setNames(unname(forecast_colors[unname(forecast_methods)]), labels)) +
      coord_equal() +
      labs(x = "False alarm rate", y = "Lows caught") +
      theme_dashboard() + legend_top(3)
  })
  
  # All 100 children (runs only when the button is clicked)
  fc_all_results <- eventReactive(input$fc_all, {
    h   <- fc_horizon()
    fl  <- fc_floor_val()
    ids <- names(all_traces)
    res <- withProgress(message = "Forecasting all 100 children...", value = 0, {
      bind_rows(lapply(ids, function(pid) {
        incProgress(1 / length(ids))
        p <- run_forecast(all_traces[[pid]], h, fl)
        tibble(participant_id = pid,
               rmse_persist = sqrt(mean((p$persist - p$actual)^2)),
               rmse_ar      = sqrt(mean((p$ar - p$actual)^2)))
      }))
    })
    res %>%
      left_join(participants %>% select(participant_id, group), by = "participant_id") %>%
      mutate(minutes = h * 5)
  })
  
  output$fc_all_text <- renderText({
    r <- fc_all_results()
    sprintf("At %d minutes ahead, the regression model made smaller errors than repeating the latest reading for %d of %d children (median error reduction: %.0f%%).",
            r$minutes[1], sum(r$rmse_ar < r$rmse_persist), nrow(r),
            median(100 * (1 - r$rmse_ar / r$rmse_persist)))
  })
  
  output$fc_all_plot <- renderGirafe({
    r <- fc_all_results() %>%
      mutate(tip = paste0("<b>", participant_id, "</b><br/>", group, "<br/>",
                          "Repeat-latest error: ", sprintf("%.1f", rmse_persist), " mg/dL<br/>",
                          "Regression error: ", sprintf("%.1f", rmse_ar), " mg/dL"))
    lim <- range(c(r$rmse_persist, r$rmse_ar))
    
    p <- ggplot(r, aes(x = rmse_persist, y = rmse_ar)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      geom_point_interactive(aes(color = group, tooltip = tip, data_id = participant_id),
                             alpha = 0.85, size = 3) +
      scale_color_manual(values = group_colors) +
      coord_equal(xlim = lim, ylim = lim) +
      labs(x = "Typical error, repeat latest reading (mg/dL)",
           y = "Typical error, regression model (mg/dL)") +
      theme_dashboard() + legend_top(2)
    make_girafe(p, 7, 6, fade = 0.15)
  })
}

shinyApp(ui, server)
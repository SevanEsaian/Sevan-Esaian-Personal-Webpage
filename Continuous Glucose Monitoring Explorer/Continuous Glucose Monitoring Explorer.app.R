# =============================================================================
# Continuous Glucose Monitoring - Pediatric Population (6 - 16 years of age)
# Tab 1: Compression Lows (Day vs. Night)
#
# NOTE: All data below are SIMULATED placeholders. When you have real data,
# replace the "SIMULATED DATA" section with your own long-format data frame
# (one row per participant x period) using the same column names
# (including a participant-level "bmi" column), and replace simulate_cgm()
# with your real 14-day sensor readings (one row per reading).
#
# PACKAGE NEEDED:  install.packages("ggiraph")
# =============================================================================

library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(ggiraph)

# -----------------------------------------------------------------------------
# SIMULATED DATA
# -----------------------------------------------------------------------------
set.seed(2026)

group_levels <- c("Type 2 Diabetes", "Prediabetes", "Healthy (No Prediabetes)")

# 100 participants total. 33 x 3 = 99, so the extra participant goes to T2D.
n_per_group <- c(34, 33, 33)

# Age ranges (years) by group
age_ranges <- list(c(6, 16), c(6, 13), c(8, 14))

participants <- bind_rows(lapply(1:3, function(i) {
  tibble(
    group = group_levels[i],
    age   = sample(age_ranges[[i]][1]:age_ranges[[i]][2], n_per_group[i], replace = TRUE),
    sex   = sample(c("Female", "Male"), n_per_group[i], replace = TRUE)
  )
})) %>%
  mutate(
    participant_id = sprintf("P%03d", row_number()),
    group = factor(group, levels = group_levels),
    # Participant-level random effects (shared by Day and Night so that each
    # child's Day and Night values are realistically paired)
    u_rate = rnorm(n()),
    u_dur  = rnorm(n()),
    u_nad  = rnorm(n()),
    u_tbr  = rnorm(n()),
    u_glu  = rnorm(n()),
    u_cv   = rnorm(n())
  )

# Group-level simulation parameters (Day vs. Night means for each metric).
# Night values differ strongly from Day values in every group, and the
# Day/Night distributions also differ between groups.
params <- tibble(
  group    = factor(group_levels, levels = group_levels),
  # compression lows per 14 days
  ev_day   = c(3.0, 2.0, 1.5),    ev_night  = c(14, 9, 6),
  # mean compression low duration (min)
  dur_day  = c(18, 16, 14),       dur_night = c(45, 38, 30),
  # nadir glucose during compression lows (mg/dL)
  nad_day  = c(66, 67, 68),       nad_night = c(48, 53, 57),
  # time below 70 mg/dL (%)
  tbr_day  = c(2.0, 1.5, 1.0),    tbr_night = c(8.0, 5.5, 3.5),
  # mean glucose (mg/dL)
  glu_day  = c(170, 118, 99),     glu_night = c(148, 105, 89),
  # glucose CV (%)
  cv_day   = c(34, 27, 20),       cv_night  = c(27, 21, 16)
)

simulate_period <- function(df, period) {
  is_night <- period == "Night"
  df %>%
    left_join(params, by = "group") %>%
    mutate(
      period        = period,
      comp_events   = pmax(1, rpois(n(), (if (is_night) ev_night else ev_day) * exp(0.25 * u_rate))),
      comp_duration = pmax(5,  (if (is_night) dur_night else dur_day) + 3 * u_dur +
                             rnorm(n(), 0, if (is_night) 7 else 3.5)),
      comp_nadir    = pmax(40, (if (is_night) nad_night else nad_day) + 3 * u_nad +
                             rnorm(n(), 0, if (is_night) 4 else 2.5)),
      tbr           = pmax(0.1, (if (is_night) tbr_night else tbr_day) * exp(0.15 * u_tbr) +
                             rnorm(n(), 0, if (is_night) 1.0 else 0.4)),
      mean_glu      = pmax(60, (if (is_night) glu_night else glu_day) + 8 * u_glu +
                             rnorm(n(), 0, if (is_night) 5 else 4)),
      cv            = pmax(8,  (if (is_night) cv_night else cv_day) + 2.5 * u_cv +
                             rnorm(n(), 0, if (is_night) 2.5 else 2))
    ) %>%
    select(participant_id, group, age, sex, period,
           comp_events, comp_duration, comp_nadir, tbr, mean_glu, cv)
}

cgm_data <- bind_rows(
  simulate_period(participants, "Day"),
  simulate_period(participants, "Night")
) %>%
  mutate(
    period = factor(period, levels = c("Day", "Night")),
    across(c(comp_duration, comp_nadir, tbr, mean_glu, cv), ~ round(.x, 1)),
    # Fixed horizontal jitter per participant (the SAME offset for Day and Night)
    # so that a child's Day and Night points can be connected by a line.
    jit_x = (((as.integer(sub("P", "", participant_id)) * 0.61803398875) %% 1) - 0.5) * 0.28
  )

# ---- Simulated BMI (one value per participant) -------------------------------
# Separate seed so the CGM data above are unaffected.
set.seed(2027)
bmi_tbl <- participants %>%
  select(participant_id, group) %>%
  mutate(
    g   = as.integer(group),
    bmi = rnorm(n(), mean = c(34, 28, 19.5)[g], sd = c(7, 5, 2.5)[g]),
    bmi = pmin(c(50, 42, 24.9)[g], pmax(c(19, 16, 15)[g], bmi)),
    bmi = round(bmi, 1)
  ) %>%
  select(participant_id, bmi)

cgm_data <- cgm_data %>% left_join(bmi_tbl, by = "participant_id")

# ---- Dropdown choices for the longitudinal CGM plot --------------------------
participant_choices <- setNames(
  participants$participant_id,
  paste0(participants$participant_id, " (", participants$group,
         ", age ", participants$age, ")")
)

# ---- Simulated 14-day CGM trace for ONE participant --------------------------
# One reading every 5 minutes for 14 days (4,032 readings). Each child's trace
# has a day/night glucose pattern, meal spikes, and compression lows that match
# the number of daytime and nighttime compression lows in cgm_data. Compression
# lows drop sharply to roughly 40 - 58 mg/dL, mostly at night.
#
# Group differences (readings range from about 40 to 250 mg/dL for most children):
#   Healthy      - stable, few severe hyperglycemic readings, fewer lows < 60
#   Prediabetes  - moderately hyperglycemic, more frequently severe, more lows
#   Type 2       - most volatile, very frequently severely hyperglycemic, most lows
simulate_cgm <- function(pid, n_days = 14, step_min = 5) {
  rows      <- cgm_data %>% filter(participant_id == pid)
  day_row   <- rows %>% filter(period == "Day")
  night_row <- rows %>% filter(period == "Night")
  grp       <- as.character(day_row$group)
  
  # Fixed seed per participant so the same child always gets the same trace
  set.seed(5000 + as.integer(sub("P", "", pid)))
  
  # Group-specific volatility and hypoglycemia settings
  sp <- switch(grp,
               "Type 2 Diabetes" = list(amp_mult = 1.6, p_big = 0.25, big_mult = 1.4,
                                        noise_sd = 3.5, noise_ar = 0.94,
                                        nad_lo = 40, nad_hi = 58, extra_hypo = 6,
                                        night_dur_scale = 3.0, extra_dur = 150),
               "Prediabetes"     = list(amp_mult = 1.6, p_big = 0.15, big_mult = 1.5,
                                        noise_sd = 2.5, noise_ar = 0.92,
                                        nad_lo = 42, nad_hi = 58, extra_hypo = 4,
                                        night_dur_scale = 2.5, extra_dur = 120),
               list(amp_mult = 1.3, p_big = 0.05, big_mult = 2.5,                   # Healthy
                    noise_sd = 1.2, noise_ar = 0.90,
                    nad_lo = 52, nad_hi = 66, extra_hypo = 0,
                    night_dur_scale = 1.0, extra_dur = 0)
  )
  
  # Child-to-child variation: some children have less variation, some more
  vol <- exp(rnorm(1, 0, 0.25))
  
  t_min <- seq(0, n_days * 1440 - step_min, by = step_min)
  n     <- length(t_min)
  hour  <- (t_min %% 1440) / 60
  
  # Smooth day/night baseline: lowest around 3 AM, highest around 3 PM
  w        <- (1 + sin(2 * pi * (hour - 9) / 24)) / 2
  baseline <- night_row$mean_glu + (day_row$mean_glu - night_row$mean_glu) * w
  
  # Meal spikes (breakfast, lunch, dinner) every day
  amp_scale <- day_row$cv / 100 * day_row$mean_glu * sp$amp_mult * vol
  meals     <- c(7.5, 12.5, 18.5) * 60
  bumps     <- numeric(n)
  for (d in 0:(n_days - 1)) {
    for (m in meals) {
      t0    <- d * 1440 + m + rnorm(1, 0, 20)
      amp   <- amp_scale * runif(1, 0.5, 1.2)
      if (runif(1) < sp$p_big) amp <- amp * sp$big_mult   # occasional large meal
      lag   <- pmax(0, t_min - t0)
      bumps <- bumps + amp * (lag / 60) * exp(1 - lag / 60) * (t_min >= t0)
    }
  }
  
  # Random sensor noise (slow wandering; larger for prediabetes and type 2)
  noise   <- as.numeric(stats::filter(rnorm(n, 0, sp$noise_sd * vol), sp$noise_ar, method = "recursive"))
  glucose <- baseline + bumps + noise
  
  # Softly cap very high values so readings level off near 250 mg/dL
  glucose <- ifelse(glucose > 200, 200 + 50 * tanh((glucose - 200) / 50), glucose)
  
  # Compression lows: sharp drop, short plateau, quick recovery
  # Every low lasts at least 80 minutes (so glucose stays <= 70 mg/dL for 1 hour or more)
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
  glucose <- add_events(glucose, night_row$comp_events, c(22, 30),
                        night_row$comp_duration * sp$night_dur_scale,
                        night_row$comp_nadir,
                        dur_sd = 6 * sp$night_dur_scale)                 # 10 PM - 6 AM (prolonged for prediabetes / type 2)
  glucose <- add_events(glucose, day_row$comp_events, c(6, 22),
                        day_row$comp_duration, day_row$comp_nadir)       # 6 AM - 10 PM
  
  # Additional prolonged nighttime hypoglycemia (about 61 - 69 mg/dL) for prediabetes and type 2
  if (sp$extra_hypo > 0) {
    glucose <- add_events(glucose, sp$extra_hypo, c(22, 30), sp$extra_dur, 65,
                          nad_lo = 61, nad_hi = 69, dur_sd = 20)
  }
  
  tibble(
    day     = t_min / 1440,                                # 0 to 14
    glucose = round(pmin(250, pmax(40, glucose)))          # mg/dL, whole numbers
  )
}

# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------
app_title <- "Continuous Glucose Monitoring - Pediatric Population (6 - 16 years of age)"

# The plot shows compression lows (events per 14 days)
metric_label <- "Compression lows (events per 14 days)"

# Colors used by the Day vs. Night distribution (box) plot
box_fill_colors  <- c("Day" = "#FFF6C2", "Night" = "#E3D9F2")   # very light yellow / very light purple
box_point_colors <- c("Day" = "#B8860B", "Night" = "#3B1A78")   # dark yellow / dark purple
box_point_size   <- 3.5                                         # size of the points
track_line_color <- "grey70"                                    # faint line connecting Day to Night

# ---- LONGITUDINAL CGM PLOT SETTINGS ------------------------------------------
range_levels <- c("Hypoglycemia", "In range", "Hyperglycemia", "Severe hyperglycemia")
range_colors <- c("Hypoglycemia"        = "#FF1493",   # bright pink
                  "In range"            = "#2E9E5B",   # green
                  "Hyperglycemia"       = "#F28C28",   # orange
                  "Severe hyperglycemia" = "#B22222")  # red
winsor_floor <- 60   # winsorizing replaces every reading below this value (mg/dL)

# Hypoglycemia: 70 or below | In range: 71-139 | Hyperglycemia: 140-179 | Severe: 180+
classify_glucose <- function(x) {
  factor(
    case_when(
      x <= 70  ~ "Hypoglycemia",
      x < 140  ~ "In range",
      x < 180  ~ "Hyperglycemia",
      TRUE     ~ "Severe hyperglycemia"
    ),
    levels = range_levels
  )
}

# Toggle switch whose label is color-coded to match the line colors
range_switch <- function(id, name, detail, col) {
  input_switch(
    id,
    label = span(span(name, style = paste0("color:", col, "; font-weight:bold;")),
                 " ", detail),
    value = TRUE
  )
}

# ---- CGM-DERIVED METRICS: RAW vs. WINSORIZED (plot 3) -------------------------
group_colors  <- c("Type 2 Diabetes"          = "#0072B2",
                   "Prediabetes"              = "#E69F00",
                   "Healthy (No Prediabetes)" = "#009E73")
metric_levels <- c("HbA1c (GMI, %)", "MAGE (mg/dL)", "J-index", "ADRR")

# MAGE: mean amplitude of glycemic excursions. Averages every glucose swing
# (peak-to-nadir or nadir-to-peak) that is larger than 1 SD of the trace.
mage <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  amps  <- numeric(0)
  dir   <- 0                       # 0 = undecided, 1 = rising, -1 = falling
  lo    <- x[1]
  hi    <- x[1]
  pivot <- x[1]
  cand  <- x[1]
  for (v in x[-1]) {
    if (dir == 0) {
      lo <- min(lo, v)
      hi <- max(hi, v)
      if (v - lo >= s) {
        dir <- 1;  pivot <- lo; cand <- v
      } else if (hi - v >= s) {
        dir <- -1; pivot <- hi; cand <- v
      }
    } else if (dir == 1) {
      if (v > cand) {
        cand <- v
      } else if (cand - v >= s) {
        amps <- c(amps, cand - pivot); pivot <- cand; dir <- -1; cand <- v
      }
    } else {
      if (v < cand) {
        cand <- v
      } else if (v - cand >= s) {
        amps <- c(amps, pivot - cand); pivot <- cand; dir <- 1; cand <- v
      }
    }
  }
  if (length(amps) == 0) NA_real_ else mean(amps)
}

# The four CGM-derived metrics for one glucose trace (mg/dL)
#   HbA1c (GMI) = 3.31 + 0.02392 x mean glucose
#   J-index     = 0.001 x (mean + SD)^2
#   ADRR        = average over days of (max daily low-risk + max daily high-risk)
calc_metrics <- function(g, day_index) {
  m  <- mean(g)
  s  <- sd(g)
  f  <- 1.509 * (log(g)^1.084 - 5.381)
  r  <- 10 * f^2
  lr <- tapply(ifelse(f < 0, r, 0), day_index, max)
  hr <- tapply(ifelse(f > 0, r, 0), day_index, max)
  c(3.31 + 0.02392 * m,
    mage(g),
    0.001 * (m + s)^2,
    mean(lr + hr))
}

# Lin's concordance correlation coefficient
lin_ccc <- function(x, y) {
  mx <- mean(x); my <- mean(y)
  2 * mean((x - mx) * (y - my)) / (mean((x - mx)^2) + mean((y - my)^2) + (mx - my)^2)
}

# All 14-day traces (computed once when the app starts)
all_traces <- setNames(lapply(participants$participant_id, simulate_cgm),
                       participants$participant_id)

# One row per participant x metric; floor_val = NULL gives the raw (not winsorized) values
metric_tbl <- function(floor_val = NULL) {
  bind_rows(lapply(names(all_traces), function(pid) {
    tr <- all_traces[[pid]]
    g  <- if (is.null(floor_val)) tr$glucose else pmax(tr$glucose, floor_val)
    tibble(participant_id = pid,
           metric = metric_levels,
           value  = unname(calc_metrics(g, floor(tr$day))))
  }))
}
raw_metrics <- metric_tbl() %>% rename(raw = value)

# Paired comparison of winsorized vs. raw values (percent change per child),
# with a 95% CI, a TOST p-value, and Lin's CCC. A metric is called equivalent
# when its 95% CI lies entirely inside the +/- margin (a conservative TOST).
eq_summary <- function(d, margin, by_group) {
  d <- if (by_group) mutate(d, grp = group) else mutate(d, grp = factor("All selected groups"))
  d %>%
    group_by(metric, grp) %>%
    summarise(n = n(),
              raw_mean  = mean(raw),
              wins_mean = mean(wins),
              diff_mean = mean(wins - raw),
              pct_mean  = mean(pct),
              pct_sd    = sd(pct),
              ccc       = lin_ccc(raw, wins),
              .groups = "drop") %>%
    filter(n >= 2) %>%
    mutate(se     = pct_sd / sqrt(n),
           tcrit  = qt(0.975, n - 1),
           lo     = pct_mean - tcrit * se,
           hi     = pct_mean + tcrit * se,
           p_tost = pmax(1 - pt((pct_mean + margin) / se, n - 1),
                         pt((pct_mean - margin) / se, n - 1)),
           equivalent = lo > -margin & hi < margin)
}

# ---- BMI SCALE SETTINGS ------------------------------------------------------
# NOTE: these are adult-style BMI cutoffs (kg/m^2). Pediatric BMI classification
# is normally based on age- and sex-specific percentiles (CDC growth charts).
bmi_min  <- 15
bmi_max  <- 50
bmi_cuts <- c(25, 30, 40)    # Normal < 25 | Overweight 25-29.9 | Obese 30-39.9 | Morbidly obese 40+
bmi_cols <- c("#2E8B57", "#F2B134", "#E8743B", "#B22222")

bmi_pos <- function(x) 100 * (x - bmi_min) / (bmi_max - bmi_min)

# Color bar with annotations that sits directly under the BMI slider
bmi_scale <- function() {
  edges <- c(bmi_min, bmi_cuts, bmi_max)
  mid   <- function(i) (bmi_pos(edges[i]) + bmi_pos(edges[i + 1])) / 2
  
  segs <- lapply(1:4, function(i) {
    div(class = "bmi-seg",
        style = paste0("left:", round(bmi_pos(edges[i]), 2), "%;",
                       "width:", round(bmi_pos(edges[i + 1]) - bmi_pos(edges[i]), 2), "%;",
                       "background:", bmi_cols[i], ";"))
  })
  
  cuts <- lapply(bmi_cuts, function(x) {
    div(class = "bmi-cut", style = paste0("left:", round(bmi_pos(x), 2), "%;"), x)
  })
  
  labels <- list(
    div(class = "bmi-label",
        style = paste0("top:0px; left:0px; color:", bmi_cols[1], ";"),
        "Normal"),
    div(class = "bmi-label",
        style = paste0("top:18px; left:", round(mid(2), 2), "%; transform:translateX(-50%); color:", bmi_cols[2], ";"),
        "Overweight"),
    div(class = "bmi-label",
        style = paste0("top:0px; left:", round(mid(3), 2), "%; transform:translateX(-50%); color:", bmi_cols[3], ";"),
        "Obese"),
    div(class = "bmi-label",
        style = paste0("top:18px; right:0px; color:", bmi_cols[4], ";"),
        "Morbidly obese")
  )
  
  div(class = "bmi-scale",
      div(class = "bmi-cuts", cuts),
      div(class = "bmi-bar", segs),
      div(class = "bmi-labels", labels))
}

bmi_css <- "
  .bmi-scale  { position: relative; height: 72px; margin: -8px 10px 12px 10px; }
  .bmi-cuts   { position: absolute; top: 0; left: 0; right: 0; height: 14px; }
  .bmi-cut    { position: absolute; transform: translateX(-50%); font-size: 13px; color: #555; }
  .bmi-bar    { position: absolute; top: 18px; left: 0; right: 0; height: 10px;
                border-radius: 5px; overflow: hidden; }
  .bmi-seg    { position: absolute; top: 0; bottom: 0; }
  .bmi-labels { position: absolute; top: 32px; left: 0; right: 0; height: 40px; }
  .bmi-label  { position: absolute; font-size: 13px; font-weight: bold;
                color: #333; white-space: nowrap; }
"

# ---- FONT SETTINGS (change these in one place) -------------------------------
plot_font      <- "Arial"   # use "sans" if Arial doesn't show up in plots (e.g., Linux)
body_size_px   <- 16        # most text: sidebar, plot labels
intro_size_px  <- 18        # "Why compression lows matter" paragraphs
title_size_px  <- 20        # major headers: plot/card titles
tab_size_px    <- 20        # navbar tab names
brand_size_px  <- 22        # app title in the top bar

app_css <- sprintf("
  body, .form-control, .form-select, .form-check-label, .control-label,
  .help-block, .form-text, .card-body, .sidebar, .bslib-value-box {
    font-family: Arial, Helvetica, sans-serif;
    font-size: %dpx;
  }
  .intro-text {
    font-family: Arial, Helvetica, sans-serif;
    font-size: %dpx;
  }
  .navbar-brand {
    font-family: Arial, Helvetica, sans-serif;
    font-size: %dpx;
    font-weight: bold;
    line-height: 1.25;
    white-space: normal;
  }
  .navbar .nav-link {
    font-family: Arial, Helvetica, sans-serif;
    font-size: %dpx;
  }
  .card-header {
    font-family: Arial, Helvetica, sans-serif;
    font-size: %dpx;
    font-weight: bold;
  }
  .sidebar h5 {
    font-family: Arial, Helvetica, sans-serif;
    font-size: %dpx;
    font-weight: bold;
  }
", body_size_px, intro_size_px, brand_size_px, tab_size_px, title_size_px, title_size_px)

# Black-and-white theme with no gridlines (used for both plots)
plot_theme_bw <- function() {
  theme_bw(base_size = body_size_px, base_family = plot_font) +
    theme(
      legend.position = "none",
      panel.grid      = element_blank(),
      axis.text       = element_text(size = body_size_px),
      axis.title      = element_text(size = body_size_px),
      strip.text      = element_text(size = body_size_px, face = "bold")
    )
}

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
ui <- page_navbar(
  # Two-line title (the browser tab still uses the one-line app_title)
  title = HTML("Continuous Glucose Monitoring <br>Pediatric Population (aged 6 - 16)"),
  window_title = app_title,
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = "Arial",
    heading_font = "Arial",
    code_font = "Arial",
    `font-size-base` = "1rem"   # 16px
  ),
  header = tags$style(HTML(paste(app_css, bmi_css))),
  fillable = FALSE,
  
  nav_panel(
    title = "Compression Lows",
    
    card(
      card_header("Why do compression lows matter for endocrinologists?"),
      card_body(
        p(class = "intro-text",
          "Continuous glucose monitors (CGMs) are small wearable sensors that track a child's blood sugar around the clock. In children, a common problem is the \"compression low.\" When a child sleeps with their body weight pressing on the sensor, usually while lying on their side, the pressure can make the sensor report a sudden drop in blood sugar that isn't real. These false lows can trigger unnecessary alarms, lead to wrong treatment decisions, and make the data harder to trust. By comparing daytime and nighttime readings below, we can see how often this happens.")
      )
    ),
    
    # ---- Plot 1: Day vs. Night distribution (with its own filters) ----------
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Filters"),
        checkboxGroupInput(
          "groups", "Group",
          choices  = group_levels,
          selected = group_levels
        ),
        sliderInput("age", "Age (years)", min = 6, max = 16, value = c(6, 16), step = 1),
        sliderInput("bmi", "Body Mass Index (kg/m\u00b2)",
                    min = bmi_min, max = bmi_max,
                    value = c(bmi_min, bmi_max), step = 1, ticks = FALSE),
        bmi_scale(),
        checkboxInput("track", "Track participants", value = FALSE),
        helpText(
          "Adds a line connecting each child's Day and Night points.",
          "Hover over any point to see the participant and highlight their Day and Night values."
        ),
        hr(),
        helpText(
          "Compression lows are sensor-related, transient false drops in glucose,",
          "typically occurring when pressure is applied to the sensor (e.g., while sleeping)."
        )
      ),
      
      # Summary value boxes
      layout_columns(
        fill = FALSE,
        col_widths = c(3, 3, 3, 3),
        value_box("Participants", textOutput("vb_n"), theme = "primary"),
        value_box("Daytime compression lows", textOutput("vb_day"), theme = "warning"),
        value_box("Nighttime compression lows", textOutput("vb_night"), theme = "secondary"),
        value_box("% occurring at night", textOutput("vb_pct"), theme = "info")
      ),
      
      # Plot
      card(
        card_header("Distribution of Daytime vs Nighttime Compression Lows & Paired Individual Change"),
        girafeOutput("box_plot", width = "100%")
      )
    ),
    
    # ---- Plot 2: Longitudinal CGM trace (with its own filters) --------------
    layout_sidebar(
      class = "mt-4",
      sidebar = sidebar(
        width = 320,
        h5("Filters"),
        selectInput("cgm_pid", "Participant",
                    choices = participant_choices,
                    selected = participants$participant_id[1]),
        hr(),
        h6("Glucose ranges"),
        range_switch("show_hypo",   "Hypoglycemia",         "(70 mg/dL or below)", range_colors[["Hypoglycemia"]]),
        range_switch("show_inrange", "In range",            "(71 - 139 mg/dL)",    range_colors[["In range"]]),
        range_switch("show_hyper",  "Hyperglycemia",        "(140 - 179 mg/dL)",   range_colors[["Hyperglycemia"]]),
        range_switch("show_severe", "Severe hyperglycemia", "(180 mg/dL and above)", range_colors[["Severe hyperglycemia"]]),
        hr(),
        input_switch("winsorize", "Winsorize line plot", value = FALSE),
        conditionalPanel(
          condition = "input.winsorize",
          sliderInput("winsor_level", "Winsorizing cutoff (mg/dL)",
                      min = 50, max = 70, value = winsor_floor, step = 1)
        ),
        helpText("Replaces every reading below the cutoff with the cutoff value (60 mg/dL by default).")
      ),
      
      card(
        card_header("Visual effects of winsorizing compression lows"),
        card_body(
          p(class = "intro-text",
            "Winsorizing is a simple way to limit the influence of extreme values. Instead of deleting unusual readings, each one is replaced with a chosen cutoff value. Here, compression lows can push sensor readings as low as 40 mg/dL, so any reading below the selected cutoff (60 mg/dL by default) is replaced with that cutoff. Turn on the winsorize switch, then move the slider to choose the cutoff and see how the sharp dips flatten out while the rest of the child's glucose pattern stays the same."),
          plotOutput("cgm_plot", height = "450px")
        )
      )
    ),
    
    # ---- Plot 3: Effect of winsorizing on CGM-derived metrics ---------------
    layout_sidebar(
      class = "mt-4",
      sidebar = sidebar(
        width = 320,
        h5("Filters"),
        checkboxGroupInput(
          "eq_groups", "Group",
          choices  = group_levels,
          selected = group_levels
        ),
        helpText("Which groups of children to include."),
        sliderInput("eq_margin", "Equivalence margin (\u00b1 % change)",
                    min = 1, max = 10, value = 5, step = 0.5),
        helpText(
          "The largest average change we would call negligible. Smaller is stricter.",
          "It sets the width of the gray band in the chart."
        ),
        selectInput("w_floor", "Winsorizing floor (mg/dL)",
                    choices = c(50, 55, 60, 65, 70), selected = 60),
        helpText(
          "Readings below this value are raised to it before the measures are recalculated.",
          "Separate from the winsorize switch and slider in the line plot above."
        ),
        input_switch("eq_split", "Split by group", value = TRUE),
        helpText("On: one row per group. Off: all selected children combined.")
      ),
      
      card(
        card_header("Does winsorizing change the summary measures?"),
        card_body(
          p(class = "intro-text",
            "Winsorizing raises any reading below a cutoff up to that cutoff. Here we ask whether doing so changes the summary measures that clinicians and researchers rely on. For each child, four measures are calculated from the 14-day trace twice, once from the original readings and once after winsorizing, so every child is compared with themselves. The measures are:"),
          tags$ul(class = "intro-text",
                  tags$li(strong("HbA1c (GMI, %): "),
                          "An estimate of the lab test that reflects a child's average blood sugar over the past few months. Higher means a higher average glucose. ",
                          em("Technical: Glucose Management Indicator = 3.31 + 0.02392 \u00d7 mean glucose (mg/dL).")),
                  tags$li(strong("MAGE (mg/dL): "),
                          "The typical size of the big swings in glucose, up and down, ignoring small wiggles. Higher means bigger swings. ",
                          em("Technical: mean amplitude of glycemic excursions, the average of peak-to-nadir (or nadir-to-peak) swings larger than 1 standard deviation of the trace.")),
                  tags$li(strong("J-index: "),
                          "A single score that combines how high glucose runs on average and how much it varies. Higher means higher and/or less stable glucose. ",
                          em("Technical: 0.001 \u00d7 (mean + SD)\u00b2.")),
                  tags$li(strong("ADRR: "),
                          "A score for how risky a child's day-to-day glucose swings are, counting both dangerous lows and dangerous highs. Higher means higher risk. ",
                          em("Technical: average daily risk range. For each day, the largest low-glucose risk plus the largest high-glucose risk from a standard risk function, averaged over the 14 days. The risk function rises steeply at low glucose, so ADRR can be more sensitive to winsorizing than the other measures."))
          )
        )
      ),
      
      card(
        girafeOutput("eq_forest", width = "100%"),
        card_footer(
          p("Each dot is the average percent change in a measure after winsorizing, with its 95% confidence interval. Hover over a dot for the exact values. The dashed line is no change, and the gray band is the negligible-change zone set by the margin. A measure counts as unchanged (equivalent) when its whole interval sits inside the band, a conservative version of the TOST equivalence test.")
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
# SERVER
# -----------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Filtered data (participant-level filters applied to both periods)
  filtered <- reactive({
    req(input$groups)
    cgm_data %>%
      filter(
        group %in% input$groups,
        age >= input$age[1], age <= input$age[2],
        bmi >= input$bmi[1], bmi <= input$bmi[2]
      )
  })
  
  # Plot data: compression lows plus plotting helpers
  plot_data <- reactive({
    filtered() %>%
      mutate(
        value    = comp_events,
        x_center = as.numeric(period),           # box position (1 = Day, 2 = Night)
        x_pt     = x_center + jit_x,             # point position (same offset Day & Night)
        # Hover text: participant ID, age, and Body Mass Index
        tip      = paste0("<b>", participant_id, "</b><br/>",
                          "Age: ", age, " years<br/>",
                          "Body Mass Index: ", bmi)
      )
  })
  
  # ---- Value boxes ----------------------------------------------------------
  output$vb_n <- renderText({
    n_distinct(filtered()$participant_id)
  })
  
  output$vb_day <- renderText({
    sum(filtered()$comp_events[filtered()$period == "Day"])
  })
  
  output$vb_night <- renderText({
    sum(filtered()$comp_events[filtered()$period == "Night"])
  })
  
  output$vb_pct <- renderText({
    d <- filtered()
    total <- sum(d$comp_events)
    if (total == 0) return("NA")
    paste0(round(100 * sum(d$comp_events[d$period == "Night"]) / total, 1), "%")
  })
  
  # ---- Plot 1: Day vs. Night distribution -----------------------------------
  output$box_plot <- renderGirafe({
    df <- plot_data()
    req(nrow(df) > 0)
    track <- isTRUE(input$track)
    
    # Boxes: very light fill, dark outline
    p <- ggplot(df) +
      geom_boxplot(aes(x = x_center, y = value, group = period,
                       fill = period, color = period),
                   alpha = 1, outlier.shape = NA, width = 0.55, linewidth = 0.8)
    
    # Faint gray line connecting each child's Day and Night points (only when tracking)
    if (track) {
      p <- p +
        geom_line_interactive(aes(x = x_pt, y = value, group = participant_id,
                                  data_id = participant_id, tooltip = tip),
                              color = track_line_color, linewidth = 0.5)
    }
    
    # Points are ALWAYS interactive: hovering shows the participant ID and
    # highlights that child's Day and Night points while the rest fade
    p <- p +
      geom_point_interactive(aes(x = x_pt, y = value, color = period,
                                 data_id = participant_id, tooltip = tip),
                             alpha = 0.85, size = box_point_size)
    
    p <- p +
      facet_wrap(~ group) +
      scale_fill_manual(values = box_fill_colors) +
      scale_color_manual(values = box_point_colors) +
      scale_x_continuous(breaks = c(1, 2), labels = c("Day", "Night"),
                         limits = c(0.5, 2.5)) +
      labs(x = NULL, y = metric_label) +
      plot_theme_bw()
    
    girafe(
      ggobj = p,
      width_svg = 12,
      height_svg = 5.5,
      options = list(
        opts_hover(css = "stroke:black;stroke-width:2px;fill-opacity:1;stroke-opacity:1;"),
        opts_hover_inv(css = "opacity:0.12;"),
        opts_tooltip(css = "background:white;color:black;border:1px solid #666;padding:6px 10px;font-family:Arial,Helvetica,sans-serif;font-size:14px;border-radius:4px;"),
        opts_toolbar(saveaspng = FALSE),
        opts_sizing(rescale = TRUE)
      )
    )
  })
  
  # ---- Plot 2: Longitudinal CGM trace ---------------------------------------
  # Raw 14-day trace for the selected participant (re-simulated only when the
  # participant changes, not when toggles change)
  cgm_series <- reactive({
    req(input$cgm_pid)
    simulate_cgm(input$cgm_pid)
  })
  
  output$cgm_plot <- renderPlot({
    raw <- cgm_series()
    
    # Which glucose ranges are switched on
    visible <- c(
      if (isTRUE(input$show_hypo))    "Hypoglycemia",
      if (isTRUE(input$show_inrange)) "In range",
      if (isTRUE(input$show_hyper))   "Hyperglycemia",
      if (isTRUE(input$show_severe))  "Severe hyperglycemia"
    )
    validate(need(length(visible) > 0,
                  "Turn on at least one glucose range to display the data."))
    
    # Winsorize: replace every reading below 60 mg/dL with 60 mg/dL
    df <- raw %>%
      mutate(
        glucose_w = if (isTRUE(input$winsorize)) pmax(glucose, input$winsor_level) else glucose,
        range     = classify_glucose(glucose_w),
        # Hidden ranges become NA so the line breaks instead of connecting across them
        glucose_plot = ifelse(range %in% visible, glucose_w, NA_real_)
      )
    
    # Fixed y-axis for every participant: always 40 - 250 mg/dL
    y_min <- 40
    y_max <- 250
    
    ggplot(df, aes(x = day, y = glucose_plot, color = range, group = 1)) +
      geom_line(linewidth = 0.7) +
      scale_color_manual(values = range_colors, breaks = range_levels, drop = FALSE) +
      scale_x_continuous(breaks = seq(0.5, 13.5, by = 1), labels = 1:14,
                         limits = c(0, 14), expand = c(0.005, 0)) +
      coord_cartesian(ylim = c(y_min, y_max)) +
      guides(color = guide_legend(override.aes = list(linewidth = 2))) +
      labs(x = "Day of monitoring", y = "Sensor glucose (mg/dL)") +
      plot_theme_bw() +
      theme(legend.position = "top",
            legend.title    = element_blank(),
            legend.text     = element_text(size = body_size_px))
  })
  
  # ---- Plot 3: Effect of winsorizing on CGM-derived metrics -----------------
  # Winsorized metrics (recomputed only when the floor changes)
  wins_tbl <- reactive({
    req(input$w_floor)
    metric_tbl(as.numeric(input$w_floor)) %>% rename(wins = value)
  })
  
  # Paired raw / winsorized values for the selected groups
  eq_data <- reactive({
    req(input$eq_groups)
    raw_metrics %>%
      left_join(wins_tbl(), by = c("participant_id", "metric")) %>%
      left_join(participants %>% select(participant_id, group), by = "participant_id") %>%
      filter(group %in% input$eq_groups, !is.na(raw), !is.na(wins), raw > 0) %>%
      mutate(metric = factor(metric, levels = metric_levels),
             pct    = 100 * (wins - raw) / raw)
  })
  
  eq_stats <- reactive({
    eq_summary(eq_data(), margin = input$eq_margin, by_group = isTRUE(input$eq_split))
  })
  
  # Forest (equivalence) plot (hover over a point for its values)
  output$eq_forest <- renderGirafe({
    s <- eq_stats()
    req(nrow(s) > 0)
    m     <- input$eq_margin
    split <- isTRUE(input$eq_split)
    s     <- s %>%
      mutate(
        tip = paste0("<b>", metric, "</b><br/>",
                     if (split) paste0(grp, "<br/>") else "",
                     "Mean change: ", sprintf("%.2f", pct_mean), "%<br/>",
                     "95% CI: ", sprintf("%.2f", lo), "% to ", sprintf("%.2f", hi), "%<br/>",
                     "Children: ", n, "<br/>",
                     ifelse(equivalent, "Within the margin (unchanged)",
                            "Not shown to be unchanged")),
        metric = factor(metric, levels = rev(metric_levels))
      )
    x_lim <- max(m, abs(s$lo), abs(s$hi), na.rm = TRUE) * 1.15
    
    p <- ggplot(s, aes(x = pct_mean, y = metric, xmin = lo, xmax = hi, color = grp)) +
      annotate("rect", xmin = -m, xmax = m, ymin = -Inf, ymax = Inf,
               fill = "grey85", alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      geom_pointrange_interactive(aes(tooltip = tip, data_id = paste(metric, grp)),
                                  position = position_dodge(width = 0.6),
                                  linewidth = 0.9, size = 0.8) +
      scale_color_manual(values = c(group_colors, "All selected groups" = "grey20")) +
      coord_cartesian(xlim = c(-x_lim, x_lim)) +
      labs(x = "Mean % change after winsorizing (95% CI)", y = NULL,
           caption = paste0("Shaded band = equivalence margin (\u00b1", m, "%)")) +
      plot_theme_bw() +
      theme(legend.position = if (split) "top" else "none",
            legend.title    = element_blank(),
            legend.text     = element_text(size = body_size_px))
    
    girafe(
      ggobj = p,
      width_svg = 12,
      height_svg = 5,
      options = list(
        opts_hover(css = "stroke-width:3px;"),
        opts_hover_inv(css = "opacity:0.3;"),
        opts_tooltip(css = "background:white;color:black;border:1px solid #666;padding:6px 10px;font-family:Arial,Helvetica,sans-serif;font-size:14px;border-radius:4px;"),
        opts_toolbar(saveaspng = FALSE),
        opts_sizing(rescale = TRUE)
      )
    )
  })
}

shinyApp(ui, server)
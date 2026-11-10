#dir.create("cgm-explorer-deploy")
#file.copy(src, "cgm-explorer-deploy/app.R")

#rsconnect::deployApp(appDir = "cgm-explorer-deploy",
#                    appName = "cgm-explorer",
#                   account = "sevanesaian",
#                  forceUpdate = TRUE)

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

# =============================================================================
# TAB 2 HELPERS: PCA + CLUSTERING OF CGM-DERIVED MEASURES   (added section)
# =============================================================================

# The ten measures fed into PCA (one value per child)
#   mean_glu = mean glucose            sd_glu = standard deviation
#   cv       = coefficient of variation (SD / mean, %)
#   tir      = % time in range (71 - 139 mg/dL, same definition as the line plot)
#   tbr      = % time below range (70 mg/dL or below)
#   tar      = % time above range (140 mg/dL or above)
#   gmi, mage, jindex, adrr = same four measures as the winsorizing section
pc_feature_cols <- c("mean_glu", "sd_glu", "cv", "tir", "tbr", "tar",
                     "gmi", "mage", "jindex", "adrr")

# Cluster colors (kept separate from the group colors so the two are not confused)
cluster_colors <- c("#56B4E9", "#CC79A7", "#D55E00", "#666666")

# Point SHAPE shows the child's true group on the cluster map
group_shapes <- c("Type 2 Diabetes"          = 16,   # circle
                  "Prediabetes"              = 17,   # triangle
                  "Healthy (No Prediabetes)" = 15)   # square

pc_highlight_color <- "#3B1A78"   # highlight for the selected component / k
pc_tip_css <- "background:white;color:black;border:1px solid #666;padding:6px 10px;font-family:Arial,Helvetica,sans-serif;font-size:14px;border-radius:4px;"

# Participant-level info used for hover text (age, BMI)
part_info <- participants %>%
  select(participant_id, group, age, sex) %>%
  left_join(bmi_tbl, by = "participant_id")

# Traces used ONLY by the Glucotypes tab. Real children do not fall neatly into
# three boxes, so about a fifth of the simulated children are partly blended with a
# child from the NEIGHBORING group (Healthy <-> Prediabetes <-> Type 2), which
# creates overlap between Glucotypes 1 & 2 and 2 & 3. Tab 1 is not affected.
#   cross_prob = share of children who cross over | w = how strongly they borrow
#   from the neighboring child (0 = none, 1 = fully the neighbor)
glucotype_traces <- local({
  set.seed(2028)
  cross_prob_edges <- 0.22     # Healthy and Type 2 children
  cross_prob_mid   <- 0.30     # Prediabetes children (can go either way)
  w_range          <- c(0.45, 0.80)
  ids <- names(all_traces)
  grp <- as.character(part_info$group[match(ids, part_info$participant_id)])
  pool <- function(g) ids[grp == g]
  out <- all_traces
  for (i in seq_along(ids)) {
    partner <- NULL
    r <- runif(1)
    if (grp[i] %in% c("Healthy (No Prediabetes)", "Type 2 Diabetes") && r < cross_prob_edges) {
      partner <- sample(pool("Prediabetes"), 1)
    } else if (grp[i] == "Prediabetes" && r < cross_prob_mid) {
      partner <- sample(pool(if (runif(1) < 0.5) "Healthy (No Prediabetes)" else "Type 2 Diabetes"), 1)
    }
    if (!is.null(partner)) {
      w <- runif(1, w_range[1], w_range[2])
      out[[i]]$glucose <- round(pmin(250, pmax(40,
                                               (1 - w) * all_traces[[i]]$glucose + w * all_traces[[partner]]$glucose)))
    }
  }
  out
})

# Builds, for one set of readings (raw or winsorized at floor_val), the ten
# per-child measures, the % of time in each of the four ranges, and each
# child's average glucose for every hour of the day.
pc_build <- function(floor_val = NULL) {
  ids  <- names(glucotype_traces)
  mets <- vector("list", length(ids))
  hrs  <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    tr <- glucotype_traces[[i]]
    g  <- if (is.null(floor_val)) tr$glucose else pmax(tr$glucose, floor_val)
    m  <- calc_metrics(g, floor(tr$day))
    mg <- mean(g)
    sg <- sd(g)
    p_low  <- 100 * mean(g <= 70)
    p_inr  <- 100 * mean(g > 70 & g < 140)
    p_high <- 100 * mean(g >= 140 & g < 180)
    p_sev  <- 100 * mean(g >= 180)
    mets[[i]] <- tibble(
      participant_id = ids[i],
      mean_glu = mg, sd_glu = sg, cv = 100 * sg / mg,
      tir = p_inr, tbr = p_low, tar = p_high + p_sev,
      gmi = m[1], mage = m[2], jindex = m[3], adrr = m[4],
      pct_low = p_low, pct_inr = p_inr, pct_high = p_high, pct_sev = p_sev
    )
    # Hour of day for each 5-minute reading (0 - 23)
    # (computed BEFORE tibble() so the new "hour" column cannot mask this vector)
    hour_idx   <- (((seq_along(g) - 1) * 5) %% 1440) %/% 60
    hourly_glu <- as.numeric(tapply(g, factor(hour_idx, levels = 0:23), mean))
    hrs[[i]] <- tibble(
      participant_id = ids[i],
      hour = 0:23,
      glu  = hourly_glu
    )
  }
  list(metrics = bind_rows(mets), hourly = bind_rows(hrs))
}

# Cache so each winsorizing cutoff is only calculated once per app session
pc_cache <- new.env()
pc_get <- function(floor_val = NULL) {
  key <- if (is.null(floor_val)) "raw" else as.character(floor_val)
  if (is.null(pc_cache[[key]])) pc_cache[[key]] <- pc_build(floor_val)
  pc_cache[[key]]
}
invisible(pc_get(NULL))   # raw version is ready when the app starts

# Plain-language name and direction for a principal component, based on which
# measure loads most strongly on it
pc_axis <- function(rot, j) {
  top <- rownames(rot)[which.max(abs(rot[, j]))]
  lab <- if (top %in% c("sd_glu", "cv", "mage", "adrr")) {
    "Glucose swings"
  } else if (top == "tbr") {
    "Time spent low"
  } else {
    "Overall glucose level"
  }
  ref <- switch(lab,
                "Glucose swings" = if ("cv" %in% rownames(rot)) "cv" else top,
                "Time spent low" = "tbr",
                if ("mean_glu" %in% rownames(rot)) "mean_glu" else top)
  list(label = lab, sign = if (rot[ref, j] < 0) -1 else 1)
}

# Fits PCA, keeps enough components to reach var_cut (%), runs k-means on them,
# and returns everything the plots need. Clusters are numbered from the lowest
# to the highest average glucose so labels stay consistent as settings change.
pc_fit <- function(feats, var_cut, k) {
  m <- feats$metrics
  X <- as.matrix(m[, pc_feature_cols])
  rownames(X) <- m$participant_id
  # Fill any missing value (e.g., MAGE undefined) with that measure's median
  for (j in seq_len(ncol(X))) {
    if (anyNA(X[, j])) X[is.na(X[, j]), j] <- median(X[, j], na.rm = TRUE)
  }
  X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
  
  pca   <- prcomp(X, center = TRUE, scale. = TRUE)
  ve    <- pca$sdev^2 / sum(pca$sdev^2)
  cum   <- cumsum(ve)
  ncomp <- which(cum >= var_cut / 100 - 1e-9)[1]
  if (is.na(ncomp)) ncomp <- length(cum)
  
  Z <- pca$x[, seq_len(ncomp), drop = FALSE]
  
  set.seed(123)
  km  <- kmeans(Z, centers = k, nstart = 25, iter.max = 50)
  ord <- order(tapply(m$mean_glu, km$cluster, mean))
  cl  <- match(km$cluster, ord)
  
  # First two components (the cluster map), named and oriented in plain language
  rot <- pca$rotation
  ax  <- lapply(1:2, function(j) pc_axis(rot, j))
  sc  <- pca$x[, 1:2, drop = FALSE]
  sc[, 1] <- sc[, 1] * ax[[1]]$sign
  sc[, 2] <- sc[, 2] * ax[[2]]$sign
  lab1 <- ax[[1]]$label
  lab2 <- ax[[2]]$label
  if (lab2 == lab1) lab2 <- paste(lab2, "(second pattern)")
  axis_x <- paste0(lab1, " (higher \u2192), ", round(100 * ve[1]), "% of information")
  axis_y <- paste0(lab2, " (higher \u2192), ", round(100 * ve[2]), "% of information")
  
  d <- m %>%
    left_join(part_info, by = "participant_id") %>%
    mutate(cluster = factor(paste("Glucotype", cl), levels = paste("Glucotype", seq_len(k))),
           pc1 = sc[, 1],
           pc2 = sc[, 2])
  
  # Plain-language fingerprint of each cluster (typical = median child)
  prof <- d %>%
    group_by(cluster) %>%
    summarise(n        = n(),
              mean_glu = median(mean_glu),
              cv       = median(cv),
              tir      = median(tir),
              tbr      = median(tbr),
              .groups = "drop") %>%
    mutate(label = case_when(
      mean_glu >= 140 & cv >= 25 ~ "High and swinging",
      mean_glu >= 140            ~ "High but steady",
      tir >= 85 & cv < 22        ~ "Stable, mostly in range",
      tbr >= 5                   ~ "Mostly in range, frequent lows",
      TRUE                       ~ "In between, some swings"
    ),
    name = paste0(cluster, ": ", label))
  
  list(data = d, prof = prof, ve = ve, cum = cum, ncomp = ncomp,
       k = k, axis_x = axis_x, axis_y = axis_y, lab1 = lab1, lab2 = lab2)
}

# Share of children who landed in the cluster dominated by their own group
pc_purity <- function(d) {
  tab <- table(d$cluster, d$group)
  sum(apply(tab, 1, max)) / nrow(d)
}

# Axis labels: glucotype number, short plain-language name, and size
pc_xlabs <- function(prof) {
  lab <- vapply(prof$label, function(z) paste(strwrap(z, 13), collapse = "\n"), character(1))
  setNames(paste0(prof$cluster, "\n", lab, "\n(n = ", prof$n, ")"), prof$cluster)
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

# =============================================================================
# NEW: SHARED DATA AND HELPERS FOR TABS 3 - 7
# =============================================================================
# The five new tabs use the SAME simulated 14-day traces (all_traces) that
# drive the Compression Lows tab, so every tab can be compared directly.
#
# PACKAGES: nlme and mgcv ship with R (no install needed).
# ASSUMPTION: monitoring starts on a MONDAY (day 1 = Mon ... days 6, 7, 13, 14
# are weekends). The simulated traces have no weekday effect built in, so the
# weekday / weekend results are expected to be flat.
library(mgcv)    # also attaches nlme

wk_names <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

# 0 - 23.99 hours -> "7:45 AM"
fmt_time <- function(h) {
  hh <- floor(h)
  mm <- round((h - hh) * 60)
  h12 <- ifelse(hh %% 12 == 0, 12, hh %% 12)
  paste0(h12, ":", sprintf("%02d", mm), " ", ifelse(hh < 12, "AM", "PM"))
}

lab_size_g <- body_size_px / ggplot2::.pt * 0.85

# One long table: every 5-minute reading for all 100 children
cgm_long <- bind_rows(lapply(names(all_traces), function(pid) {
  tr <- all_traces[[pid]]
  tibble(participant_id = pid, t_min = round(tr$day * 1440), glucose = tr$glucose)
})) %>%
  left_join(part_info, by = "participant_id") %>%
  mutate(
    day_num  = t_min %/% 1440 + 1,
    hour_dec = (t_min %% 1440) / 60,
    hour     = floor(hour_dec),
    bin      = floor(hour_dec * 4) / 4,                      # 15-minute bins
    weekday  = factor(wk_names[(day_num - 1) %% 7 + 1], levels = wk_names),
    period   = factor(ifelse(hour >= 22 | hour < 6, "Night", "Day"),   # same 10 PM - 6 AM as Tab 1
                      levels = c("Day", "Night"))
  ) %>%
  select(participant_id, group, age, sex, bmi, day_num, hour_dec, hour, bin,
         weekday, period, glucose)

# Hourly averages (used by the Modeling tab: 100 children x 14 days x 24 hours)
hourly_tbl <- cgm_long %>%
  group_by(participant_id, group, age, sex, bmi, day_num, hour) %>%
  summarise(mean_glu = mean(glucose), min_glu = min(glucose),
            sd_glu = sd(glucose), .groups = "drop") %>%
  mutate(
    t       = (day_num - 1) * 24 + hour,
    sin_h   = sin(2 * pi * hour / 24),
    cos_h   = cos(2 * pi * hour / 24),
    night   = factor(ifelse(hour >= 22 | hour < 6, "Night", "Day"), levels = c("Day", "Night")),
    weekend = factor(ifelse((day_num - 1) %% 7 >= 5, "Weekend", "Weekday"),
                     levels = c("Weekday", "Weekend")),
    sex     = factor(sex, levels = c("Female", "Male")),
    participant_id = factor(participant_id)
  ) %>%
  arrange(participant_id, t) %>%
  as.data.frame()

# Shared ggiraph options (same look as the earlier tabs)
girafe_std <- function(p, w, h, inv = 0.2) {
  girafe(
    ggobj = p, width_svg = w, height_svg = h,
    options = list(
      opts_hover(css = "stroke:black;stroke-width:2px;fill-opacity:1;stroke-opacity:1;"),
      opts_hover_inv(css = paste0("opacity:", inv, ";")),
      opts_tooltip(css = pc_tip_css),
      opts_toolbar(saveaspng = FALSE),
      opts_sizing(rescale = TRUE)
    )
  )
}

# Glucose color scale that matches the range colors (pink = low, green = in range ...)
glu_fill_scale <- function(name) {
  scale_fill_gradientn(
    name = name,
    colours = unname(range_colors[c("Hypoglycemia", "In range", "Hyperglycemia", "Severe hyperglycemia")]),
    values  = c(0, 0.25, 0.5, 1),            # 50, 100, 150, 250 mg/dL
    limits  = c(50, 250), oob = scales::squish
  )
}

# Box plot per group with one hoverable dot per child
box_by_group <- function(df, v, ylab) {
  df$value <- v
  df$tip <- paste0("<b>", df$participant_id, "</b><br/>", df$group, "<br/>",
                   ylab, ": ", sprintf("%.1f", df$value))
  ggplot(df, aes(x = group, y = value)) +
    geom_boxplot(aes(fill = group, color = group), alpha = 0.3,
                 outlier.shape = NA, width = 0.55, linewidth = 0.8) +
    geom_point_interactive(aes(color = group, tooltip = tip, data_id = participant_id),
                           position = position_jitter(width = 0.15, height = 0, seed = 7),
                           size = box_point_size, alpha = 0.85) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    scale_x_discrete(labels = scales::label_wrap(14)) +
    labs(x = NULL, y = ylab) +
    plot_theme_bw()
}

# ---- TAB 3 HELPERS: glycemic summary -----------------------------------------
ov_metric_choices <- c(
  "Time in range, % (71 - 139 mg/dL)"   = "pct_inr",
  "Time below range, % (70 or below)"    = "pct_low",
  "Time below 54 mg/dL, %"               = "pct_low2",
  "Time above range, % (140 and above)"  = "pct_above",
  "Mean glucose (mg/dL)"                 = "mean_glu",
  "Glucose CV (%)"                       = "cv",
  "GMI (estimated HbA1c, %)"             = "gmi"
)

# One row per child, from any slice of cgm_long
child_summary <- function(d) {
  d %>%
    group_by(participant_id, group) %>%
    summarise(mean_glu = mean(glucose),
              sd_glu   = sd(glucose),
              pct_low  = 100 * mean(glucose <= 70),
              pct_low2 = 100 * mean(glucose < 54),
              pct_inr  = 100 * mean(glucose > 70 & glucose < 140),
              pct_high = 100 * mean(glucose >= 140 & glucose < 180),
              pct_sev  = 100 * mean(glucose >= 180),
              .groups = "drop") %>%
    mutate(cv = 100 * sd_glu / mean_glu,
           gmi = 3.31 + 0.02392 * mean_glu,
           pct_above = pct_high + pct_sev)
}

# ---- TAB 5 HELPERS: variability and low-glucose episodes ---------------------
var_metric_choices <- c(
  "Standard deviation (mg/dL)"        = "sd",
  "Coefficient of variation (%)"      = "cv",
  "MAGE (mg/dL)"                      = "mage",
  "CONGA, 1 hour (mg/dL)"             = "conga1",
  "J-index"                           = "jindex",
  "Low blood glucose index (LBGI)"    = "lbgi",
  "High blood glucose index (HBGI)"   = "hbgi"
)

# Variability measures for one glucose trace
#   CONGA(1h) = SD of the differences between each reading and the one 1 hour earlier
#   LBGI / HBGI = mean low-side / high-side blood glucose risk (Kovatchev)
var_calc <- function(g) {
  f  <- 1.509 * (log(g)^1.084 - 5.381)
  rl <- ifelse(f < 0, 10 * f^2, 0)
  rh <- ifelse(f > 0, 10 * f^2, 0)
  d12 <- g[13:length(g)] - g[1:(length(g) - 12)]
  tibble(sd     = sd(g),
         cv     = 100 * sd(g) / mean(g),
         mage   = mage(g),
         conga1 = sd(d12),
         jindex = 0.001 * (mean(g) + sd(g))^2,
         lbgi   = mean(rl),
         hbgi   = mean(rh))
}

var_cache <- new.env()
var_get <- function(floor_val = NULL) {
  key <- if (is.null(floor_val)) "raw" else as.character(floor_val)
  if (is.null(var_cache[[key]])) {
    var_cache[[key]] <- bind_rows(lapply(names(all_traces), function(pid) {
      g <- all_traces[[pid]]$glucose
      if (!is.null(floor_val)) g <- pmax(g, floor_val)
      bind_cols(tibble(participant_id = pid), var_calc(g))
    })) %>%
      left_join(part_info, by = "participant_id")
  }
  var_cache[[key]]
}

# Episode = at least 15 minutes (3 readings) in a row at 70 mg/dL or below
detect_episodes <- function(g, thr = 70, min_len = 3) {
  r    <- rle(g <= thr)
  ends <- cumsum(r$lengths)
  strt <- ends - r$lengths + 1
  keep <- r$values & r$lengths >= min_len
  if (!any(keep)) return(NULL)
  s <- strt[keep]
  e <- ends[keep]
  tibble(start   = s,
         end     = e,
         dur_min = (e - s + 1) * 5,
         nadir   = vapply(seq_along(s), function(i) min(g[s[i]:e[i]]), numeric(1)))
}

ep_level_cols <- c("Level 1 (54 - 70)" = "#FF69B4", "Level 2 (below 54)" = "#8B0A50")

ep_cache <- new.env()
ep_get <- function(floor_val = NULL) {
  key <- if (is.null(floor_val)) "raw" else as.character(floor_val)
  if (is.null(ep_cache[[key]])) {
    eps <- bind_rows(lapply(names(all_traces), function(pid) {
      g <- all_traces[[pid]]$glucose
      if (!is.null(floor_val)) g <- pmax(g, floor_val)
      e <- detect_episodes(g)
      if (is.null(e)) return(NULL)
      e$participant_id <- pid
      e
    })) %>%
      left_join(part_info, by = "participant_id") %>%
      mutate(day_num    = floor(((start - 1) * 5) / 1440) + 1,
             onset_hour = (((start - 1) * 5) %% 1440) / 60,
             level  = factor(ifelse(nadir < 54, "Level 2 (below 54)", "Level 1 (54 - 70)"),
                             levels = names(ep_level_cols)),
             period = ifelse(onset_hour >= 22 | onset_hour < 6, "Night", "Day"))
    ep_cache[[key]] <- eps
  }
  ep_cache[[key]]
}
invisible(var_get(NULL))   # raw versions are ready when the app starts
invisible(ep_get(NULL))

# ---- TAB 6 HELPERS: statistical modeling -------------------------------------
md_outcome_choices <- c("Hourly mean glucose (mg/dL)"    = "mean_glu",
                        "Hourly minimum glucose (mg/dL)" = "min_glu",
                        "Hourly SD of glucose (mg/dL)"   = "sd_glu")

md_term_choices <- c("Group"                         = "group",
                     "Night (10 PM - 6 AM)"          = "night",
                     "Time of day (24-hour cycle)"   = "circ",
                     "Weekend"                       = "weekend",
                     "Age"                           = "age",
                     "Body Mass Index"               = "bmi",
                     "Sex"                           = "sex",
                     "Group \u00d7 Night interaction" = "grp_night")

clean_term <- function(x) {
  map <- c("groupPrediabetes"                         = "Prediabetes (vs Type 2)",
           "groupHealthy (No Prediabetes)"            = "Healthy (vs Type 2)",
           "nightNight"                               = "Night (vs Day)",
           "sin_h"                                    = "Time of day: sine",
           "cos_h"                                    = "Time of day: cosine",
           "weekendWeekend"                           = "Weekend (vs weekday)",
           "age"                                      = "Age (per year)",
           "bmi"                                      = "BMI (per kg/m\u00b2)",
           "sexMale"                                  = "Male (vs female)",
           "groupPrediabetes:nightNight"              = "Prediabetes \u00d7 Night",
           "groupHealthy (No Prediabetes):nightNight" = "Healthy \u00d7 Night")
  ifelse(x %in% names(map), unname(map[x]), x)
}

# ---- TAB 7 HELPERS: forecasting ----------------------------------------------
fc_models <- c(persist = "Persistence", trend = "Damped trend", ar = "AR + time of day")
fc_cols   <- c("Actual" = "black", "Persistence" = "#E69F00",
               "Damped trend" = "#56B4E9", "AR + time of day" = "#3B1A78")

# Rolling-origin evaluation for ONE child. Models are trained only on days BEFORE
# each test block (days 1-8 -> test 9-10; days 1-10 -> test 11-12; days 1-12 -> test 13-14).
#   Persistence      : forecast = current reading
#   Damped trend     : current reading + half of the last 15-minute slope carried forward
#   AR + time of day : linear regression on current + lagged readings + 24 h / 12 h cycles
fc_run <- function(tr, h, floor_val = NULL, test_starts = c(9, 11, 13)) {
  g <- if (is.null(floor_val)) tr$glucose else pmax(tr$glucose, floor_val)
  n <- length(g)
  hour <- (((seq_len(n) - 1) * 5) %% 1440) / 60
  lg <- function(k) c(rep(NA_real_, k), g[seq_len(n - k)])
  ld <- function(k) c(g[(k + 1):n], rep(NA_real_, k))
  d <- data.frame(
    t = seq_len(n), day = floor(tr$day) + 1, actual = ld(h),
    g0 = g, l1 = lg(1), l2 = lg(2), l3 = lg(3), l6 = lg(6), l12 = lg(12),
    s1 = sin(2 * pi * hour / 24), c1 = cos(2 * pi * hour / 24),
    s2 = sin(4 * pi * hour / 24), c2 = cos(4 * pi * hour / 24)
  )
  d$persist <- d$g0
  d$trend   <- d$g0 + 0.5 * (d$g0 - d$l3) / 3 * h
  # clock time of the TARGET (hours since midnight of the day the forecast was made)
  d$xh <- ((d$t + h - 1) * 5) / 60 - (d$day - 1) * 24
  
  out <- vector("list", length(test_starts))
  for (k in seq_along(test_starts)) {
    ts <- test_starts[k]
    first_test_t <- min(d$t[d$day == ts])
    tr_d <- d[d$t + h < first_test_t & complete.cases(d), ]
    te_d <- d[d$day %in% c(ts, ts + 1) & !is.na(d$actual) & !is.na(d$l12), ]
    fit  <- lm(actual ~ g0 + l1 + l2 + l3 + l6 + l12 + s1 + c1 + s2 + c2, data = tr_d)
    te_d$ar    <- as.numeric(predict(fit, te_d))
    te_d$sigma <- summary(fit)$sigma
    te_d$fold  <- paste0("Test days ", ts, "-", ts + 1)
    out[[k]] <- te_d
  }
  bind_rows(out)
}

# Long format (one row per forecast x model); model stays a character column
fc_long <- function(p) {
  bind_rows(lapply(names(fc_models), function(m) {
    tibble(fold = p$fold, day = p$day, xh = p$xh, actual = p$actual,
           model = fc_models[[m]], pred = p[[m]])
  }))
}

fc_metric_tbl <- function(lp) {
  lp %>%
    group_by(Model = model) %>%
    summarise(`RMSE (mg/dL)` = sqrt(mean((pred - actual)^2)),
              `MAE (mg/dL)`  = mean(abs(pred - actual)),
              `MARD (%)`     = 100 * mean(abs(pred - actual) / actual),
              `Within 15 mg/dL or 15% (%)` =
                100 * mean(ifelse(actual <= 100, abs(pred - actual) <= 15,
                                  abs(pred - actual) / actual <= 0.15)),
              .groups = "drop") %>%
    mutate(Model = factor(Model, levels = unname(fc_models))) %>%
    arrange(Model) %>%
    mutate(Model = as.character(Model))
}

# ROC curve + AUC (score: higher = more likely to be positive)
roc_curve <- function(score, y) {
  thr <- sort(unique(score), decreasing = TRUE)
  tpr <- vapply(thr, function(t) mean(score[y] >= t), numeric(1))
  fpr <- vapply(thr, function(t) mean(score[!y] >= t), numeric(1))
  tpr <- c(0, tpr, 1)
  fpr <- c(0, fpr, 1)
  auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  list(df = tibble(fpr = fpr, tpr = tpr), auc = auc)
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
        sliderInput("eq_margin", "How much change is too much? (\u00b1 %)",
                    min = 1, max = 10, value = 5, step = 0.5),
        helpText(
          "The largest average change we would still call negligible. Smaller is stricter.",
          "It sets the gray band in the charts and the unchanged / changed verdict."
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
      
      # ---- Short answer first -------------------------------------------------
      card(
        card_header("Does winsorizing change the summary measures?"),
        card_body(
          uiOutput("eq_verdict"),
          p(class = "intro-text",
            "Winsorizing raises readings below a cutoff up to that cutoff. Four measures are calculated for each child from the original and winsorized readings:"),
          tags$ul(class = "intro-text",
                  tags$li(strong("HbA1c (GMI, %): "), "estimated average blood sugar over the past few months."),
                  tags$li(strong("MAGE (mg/dL): "), "typical size of the big glucose swings."),
                  tags$li(strong("J-index: "), "combines how high and how variable glucose is."),
                  tags$li(strong("ADRR: "), "risk from dangerous lows and highs; the most sensitive to winsorizing.")
          ),
          accordion(
            open = FALSE,
            class = "mt-2",
            accordion_panel(
              "Statistical details",
              tags$ul(
                tags$li(strong("GMI: "), "3.31 + 0.02392 \u00d7 mean glucose (mg/dL)."),
                tags$li(strong("MAGE: "), "mean of peak-to-nadir swings larger than 1 SD."),
                tags$li(strong("J-index: "), "0.001 \u00d7 (mean + SD)\u00b2."),
                tags$li(strong("ADRR: "), "daily max low-risk + max high-risk, averaged over 14 days."),
                tags$li(strong("% change: "), "100 \u00d7 (winsorized \u2212 original) / original, per child."),
                tags$li(strong("Unchanged (TOST): "), "the whole 95% CI of the mean % change sits inside the \u00b1 margin."),
                tags$li(strong("CCC: "), "Lin's concordance; 1 = winsorized values identical to originals.")
              )
            )
          )
        )
      ),
      
      # ---- Average change and child-by-child, side by side ---------------------
      layout_columns(
        fill = FALSE,
        col_widths = c(6, 6),
        card(
          card_header("Average change per measure"),
          girafeOutput("eq_forest", width = "100%"),
          card_footer(
            p("Dot = average % change with 95% CI. Inside the gray band = unchanged.")
          )
        ),
        card(
          card_header("Does anyone move?"),
          girafeOutput("eq_scatter", width = "100%"),
          card_footer(
            p("One dot per child; on the diagonal = no change. CCC = agreement (1 = identical). % = children changed by more than the margin.")
          )
        )
      )
    )
  ),
  
  # ===========================================================================
  # TAB 2: Glucotypes (PCA + k-means)   (added)
  # ===========================================================================
  nav_panel(
    title = "Glucotypes",
    
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Settings"),
        sliderInput("pc_k", "Number of glucotypes",
                    min = 2, max = 4, value = 3, step = 1),
        helpText("How many subgroups the computer looks for."),
        sliderInput("pc_var", "Information to keep (%)",
                    min = 50, max = 95, value = 85, step = 1),
        helpText("Share of the original information used to form the glucotypes."),
        hr(),
        input_switch("pc_wins", "Use winsorized readings", value = FALSE),
        conditionalPanel(
          condition = "input.pc_wins",
          sliderInput("pc_floor", "Winsorizing cutoff (mg/dL)",
                      min = 50, max = 70, value = winsor_floor, step = 1)
        ),
        helpText("Raises readings below the cutoff to the cutoff (removes false lows), then recalculates.")
      ),
      
      # PCA map, with a caption above it explaining what it shows
      card(
        card_header("Glucotype map (PCA): each dot is one child"),
        card_body(
          uiOutput("pc_caption"),
          girafeOutput("pc_map", width = "100%")
        ),
        card_footer(textOutput("pc_switch_text"))
      ),
      
      # Stacked bars: who is in each glucotype / time in ranges
      layout_columns(
        fill = FALSE,
        col_widths = c(6, 6),
        card(
          card_header("Who is in each glucotype?"),
          girafeOutput("pc_comp_plot", width = "100%"),
          card_footer(textOutput("pc_comp_text"))
        ),
        card(
          card_header("Time in ranges"),
          girafeOutput("pc_range_plot", width = "100%"),
          card_footer(p("Average share of the 14 days spent in each range."))
        )
      ),
      
      # Typical day
      card(
        card_header("A typical day in each glucotype"),
        girafeOutput("pc_day_plot", width = "100%"),
        card_footer(
          p("Line = median glucose by hour. Band = most children (10th to 90th percentile). Dashed lines = 70 and 140 mg/dL.")
        )
      )
    )
  ),
  
  # ===========================================================================
  # NEW TAB 3: Glycemic Summary
  # ===========================================================================
  nav_panel(
    title = "Glycemic Summary",
    
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Filters"),
        checkboxGroupInput("ov_groups", "Group",
                           choices = group_levels, selected = group_levels),
        selectInput("ov_pid", "Participant",
                    choices = c("All participants" = "all", participant_choices)),
        helpText("Choosing one child overrides the group filter."),
        radioButtons("ov_period", "Time of day",
                     choices = c("All day" = "all",
                                 "Day (6 AM - 10 PM)" = "Day",
                                 "Night (10 PM - 6 AM)" = "Night")),
        sliderInput("ov_days", "Days of monitoring",
                    min = 1, max = 14, value = c(1, 14), step = 1),
        hr(),
        selectInput("ov_metric", "Measure shown for each child",
                    choices = ov_metric_choices),
        helpText("Time in range uses the same cutoffs as the rest of the app: 71 - 139 mg/dL.")
      ),
      
      layout_columns(
        fill = FALSE,
        col_widths = c(2, 3, 2, 2, 3),
        value_box("Participants", textOutput("ov_n"), theme = "primary"),
        value_box("Mean glucose", textOutput("ov_mean"), theme = "secondary"),
        value_box("GMI (est. HbA1c)", textOutput("ov_gmi"), theme = "info"),
        value_box("Glucose CV", textOutput("ov_cv"), theme = "warning"),
        value_box("Time in range", textOutput("ov_tir"), theme = "success")
      ),
      
      layout_columns(
        fill = FALSE,
        col_widths = c(6, 6),
        card(
          card_header("Where does the time go?"),
          girafeOutput("ov_ranges", width = "100%"),
          card_footer(p("Average share of readings in each glucose range, by group."))
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
        card_footer(p("Averages of the per-child values for the selected filters."))
      )
    )
  ),
  
  # ===========================================================================
  # NEW TAB 4: Daily Patterns (AGP + heatmaps)
  # ===========================================================================
  nav_panel(
    title = "Daily Patterns",
    
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Filters"),
        checkboxGroupInput("agp_groups", "Group",
                           choices = group_levels, selected = group_levels),
        helpText("Used for the profile and the weekday plot."),
        hr(),
        selectInput("agp_pid", "Participant (day-by-day heatmap)",
                    choices = participant_choices,
                    selected = participants$participant_id[1]),
        selectInput("agp_metric", "Heatmap shows",
                    choices = c("Mean glucose (mg/dL)" = "mean_glu",
                                "% time below range"   = "pct_low",
                                "% time in range"      = "pct_inr"))
      ),
      
      card(
        card_header("Ambulatory glucose profile (AGP) by group"),
        girafeOutput("agp_plot", width = "100%"),
        card_footer(
          p("Line = median glucose by time of day. Dark band = middle half of readings (25th - 75th percentile). Light band = 5th - 95th percentile. Dashed lines = 70 and 140 mg/dL.")
        )
      ),
      
      card(
        card_header("One child, every day and hour"),
        girafeOutput("agp_heat", width = "100%"),
        card_footer(p("Each tile is one hour of one day. Hover for the value."))
      ),
      
      card(
        card_header("Hour of day by day of the week"),
        girafeOutput("agp_week", width = "100%"),
        card_footer(p("Average glucose for the selected groups. Day 1 is assumed to be a Monday. The simulated data contain no weekday effect, so this is expected to look flat across days."))
      )
    )
  ),
  
  # ===========================================================================
  # NEW TAB 5: Variability & Events
  # ===========================================================================
  nav_panel(
    title = "Variability & Events",
    
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Filters"),
        checkboxGroupInput("vr_groups", "Group",
                           choices = group_levels, selected = group_levels),
        selectInput("vr_metric", "Variability measure", choices = var_metric_choices),
        hr(),
        input_switch("vr_wins", "Use winsorized readings", value = FALSE),
        conditionalPanel(
          condition = "input.vr_wins",
          sliderInput("vr_floor", "Winsorizing cutoff (mg/dL)",
                      min = 50, max = 70, value = winsor_floor, step = 1)
        ),
        helpText("Raises readings below the cutoff to the cutoff, then recalculates the measures and re-detects episodes."),
        hr(),
        selectInput("vr_pid", "Participant (trace viewer)",
                    choices = participant_choices,
                    selected = participants$participant_id[1]),
        sliderInput("vr_days", "Days shown", min = 1, max = 14, value = c(1, 3), step = 1)
      ),
      
      layout_columns(
        fill = FALSE,
        col_widths = c(3, 3, 3, 3),
        value_box("Low-glucose episodes", textOutput("vr_n_ep"), theme = "primary"),
        value_box("Per child (14 days)", textOutput("vr_per_child"), theme = "secondary"),
        value_box("Median duration", textOutput("vr_dur"), theme = "info"),
        value_box("Starting at night", textOutput("vr_night"), theme = "warning")
      ),
      
      card(
        card_header("Glucose variability by group"),
        girafeOutput("vr_box", width = "100%"),
        card_footer(p("One dot per child. MAGE = typical size of big swings; CONGA = change over 1 hour; LBGI / HBGI = risk of lows / highs."))
      ),
      
      layout_columns(
        fill = FALSE,
        col_widths = c(6, 6),
        card(
          card_header("When do lows start?"),
          girafeOutput("vr_onset", width = "100%"),
          card_footer(p("Episodes per child by hour of onset. An episode is 15 minutes or more at 70 mg/dL or below; Level 2 means the lowest reading was below 54."))
        ),
        card(
          card_header("How long and how deep?"),
          girafeOutput("vr_scatter", width = "100%"),
          card_footer(p("Each dot is one episode. Hover to highlight all episodes of that child."))
        )
      ),
      
      card(
        card_header("Trace viewer with episodes highlighted"),
        plotOutput("vr_trace", height = "420px"),
        card_footer(p("Pink bands = detected low-glucose episodes. Dashed lines = 70 and 54 mg/dL. Use the \"Days shown\" slider to zoom."))
      )
    )
  ),
  
  # ===========================================================================
  # NEW TAB 6: Modeling
  # ===========================================================================
  nav_panel(
    title = "Modeling",
    
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Model"),
        selectInput("md_outcome", "Outcome", choices = md_outcome_choices),
        checkboxGroupInput("md_terms", "Predictors",
                           choices  = md_term_choices,
                           selected = c("group", "night", "circ", "weekend")),
        helpText("Each child gets their own baseline (random intercept). Picking the interaction adds Group and Night automatically."),
        input_switch("md_ar", "AR(1) errors (hourly values are correlated)", value = TRUE),
        helpText("Turn this off to see what happens to the diagnostics when the correlation is ignored."),
        actionButton("md_fit", "Fit model", class = "btn-primary w-100"),
        helpText("The AR(1) fit can take 10 - 30 seconds.")
      ),
      
      card(
        card_header("What was fit?"),
        card_body(uiOutput("md_summary"))
      ),
      
      layout_columns(
        fill = FALSE,
        col_widths = c(6, 6),
        card(
          card_header("Effects with 95% confidence intervals"),
          girafeOutput("md_forest", width = "100%"),
          card_footer(p("Dark = interval excludes 0. Group effects are compared with Type 2 Diabetes, night with day."))
        ),
        card(
          card_header("Is leftover correlation gone?"),
          plotOutput("md_acf", height = "380px"),
          card_footer(p("Average autocorrelation of the model's residuals within each child. Bars inside the dashed lines are consistent with no remaining correlation."))
        )
      ),
      
      layout_columns(
        fill = FALSE,
        col_widths = c(6, 6),
        card(
          card_header("Residuals vs. fitted values"),
          plotOutput("md_resid", height = "380px"),
          card_footer(p("A random sample of 3,000 hourly values. A flat line near 0 and an even spread are what we hope to see."))
        ),
        card(
          card_header("Circadian curve by group (GAM)"),
          plotOutput("md_gam", height = "380px"),
          card_footer(p("Smooth 24-hour curves of hourly mean glucose with 95% bands, from a GAM with a cyclic smooth per group and a random intercept per child."))
        )
      ),
      
      card(
        accordion(
          open = FALSE,
          accordion_panel("Full model output", verbatimTextOutput("md_text"))
        )
      )
    )
  ),
  
  # ===========================================================================
  # NEW TAB 7: Forecasting
  # ===========================================================================
  nav_panel(
    title = "Forecasting",
    
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Settings"),
        selectInput("fc_pid", "Participant",
                    choices = participant_choices,
                    selected = participants$participant_id[1]),
        selectInput("fc_h", "Forecast horizon",
                    choices = c("15 minutes" = 3, "30 minutes" = 6,
                                "45 minutes" = 9, "60 minutes" = 12),
                    selected = 6),
        selectInput("fc_day", "Day shown in the trace plot",
                    choices = 9:14, selected = 11),
        helpText("Only days 9 - 14 are forecasts; days 1 - 8 are needed for training."),
        hr(),
        input_switch("fc_wins", "Use winsorized readings", value = FALSE),
        conditionalPanel(
          condition = "input.fc_wins",
          sliderInput("fc_floor", "Winsorizing cutoff (mg/dL)",
                      min = 50, max = 70, value = winsor_floor, step = 1)
        ),
        helpText("Raises readings below the cutoff to the cutoff before forecasting; forecasts are scored against the same readings.")
      ),
      
      card(
        card_header("Can we see a low coming?"),
        card_body(
          p(class = "intro-text",
            "Each forecast uses only the readings available at that moment to predict glucose a short time ahead. Three methods are compared: carrying the current reading forward (persistence), extending the recent slope, and a regression on recent readings plus time of day. To avoid peeking at the future, models are trained only on earlier days and tested on later ones (rolling origin: train days 1 - 8, test 9 - 10; train 1 - 10, test 11 - 12; train 1 - 12, test 13 - 14). Sudden compression lows are very hard to forecast, so compare the error with and without winsorizing.")
        )
      ),
      
      card(
        card_header("Forecast vs. what actually happened"),
        plotOutput("fc_trace", height = "420px"),
        card_footer(p("Each forecast is drawn at the time it is predicting. Purple band = 95% prediction interval for the AR + time-of-day model. Dashed lines = 70 and 140 mg/dL."))
      ),
      
      layout_columns(
        fill = FALSE,
        col_widths = c(6, 6),
        card(
          card_header("Accuracy on all test days"),
          tableOutput("fc_table"),
          card_footer(p("RMSE and MAE: average error in mg/dL. MARD: mean absolute relative difference. Last column: share of forecasts within 15 mg/dL (at 100 or below) or 15% (above 100) of the real value."))
        ),
        card(
          card_header("Accuracy in each test block"),
          girafeOutput("fc_folds", width = "100%"),
          card_footer(p("Lower is better. Stable ordering across blocks means the result is not a fluke of one block."))
        )
      ),
      
      card(
        card_header("Can the forecast flag a coming low?"),
        plotOutput("fc_roc", height = "420px"),
        card_footer(p("ROC curve for predicting that glucose will be 70 mg/dL or below at the chosen horizon. Higher and further left is better; AUC = 0.5 is a coin flip."))
      ),
      
      card(
        card_header("All 100 children: AR + time of day vs. persistence"),
        card_body(
          actionButton("fc_all", "Run all 100 children", class = "btn-primary"),
          helpText("Uses the horizon and winsorizing settings at the moment you click. Takes a few seconds."),
          textOutput("fc_all_text"),
          girafeOutput("fc_all_plot", width = "100%")
        ),
        card_footer(p("One dot per child. Below the diagonal = the AR + time-of-day model had lower error than persistence."))
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
      width_svg = 7,
      height_svg = 6,
      options = list(
        opts_hover(css = "stroke-width:3px;"),
        opts_hover_inv(css = "opacity:0.3;"),
        opts_tooltip(css = "background:white;color:black;border:1px solid #666;padding:6px 10px;font-family:Arial,Helvetica,sans-serif;font-size:14px;border-radius:4px;"),
        opts_toolbar(saveaspng = FALSE),
        opts_sizing(rescale = TRUE)
      )
    )
  })
  
  # ---- Short answer (plain-language verdict, all selected children combined) --
  eq_verdict_stats <- reactive({
    eq_summary(eq_data(), margin = input$eq_margin, by_group = FALSE)
  })
  
  output$eq_verdict <- renderUI({
    s <- eq_verdict_stats()
    validate(need(nrow(s) > 0, "Select at least one group to see the answer."))
    m     <- input$eq_margin
    n_ok  <- sum(s$equivalent)
    n_all <- nrow(s)
    
    headline <- if (n_ok == n_all) {
      sprintf("All %d measures are essentially unchanged by winsorizing at %s mg/dL.",
              n_all, input$w_floor)
    } else if (n_ok == 0) {
      sprintf("None of the %d measures can be called unchanged after winsorizing at %s mg/dL.",
              n_all, input$w_floor)
    } else {
      sprintf("%d of %d measures are essentially unchanged by winsorizing at %s mg/dL.",
              n_ok, n_all, input$w_floor)
    }
    
    flagged <- s[!s$equivalent, ]
    detail <- if (nrow(flagged) == 0) NULL else {
      div(class = "intro-text",
          paste0("Not shown unchanged: ",
                 paste0(as.character(flagged$metric), " (", sprintf("%+.1f", flagged$pct_mean), "%)",
                        collapse = "; "), "."))
    }
    
    badges <- lapply(seq_len(nrow(s)), function(i) {
      ok <- s$equivalent[i]
      span(
        paste0(if (ok) "\u2714 " else "\u26A0 ",
               as.character(s$metric[i]), ": ", sprintf("%+.1f%%", s$pct_mean[i])),
        style = paste0("display:inline-block; padding:6px 12px; border-radius:6px; ",
                       "font-weight:bold; background:", if (ok) "#2E9E5B" else "#F2B134", ";",
                       "color:", if (ok) "white" else "#333", ";")
      )
    })
    
    tagList(
      div(style = "font-size:22px; font-weight:bold; margin-bottom:8px;", headline),
      detail,
      div(style = "display:flex; flex-wrap:wrap; gap:10px; margin:12px 0;", badges),
      div(class = "text-muted",
          sprintf("All selected children combined; margin \u00b1%s%%.", m))
    )
  })
  
  # ---- Child by child: original vs. winsorized -------------------------------
  output$eq_scatter <- renderGirafe({
    d <- eq_data()
    req(nrow(d) > 0)
    m <- input$eq_margin
    
    # Per-measure agreement (Lin's CCC) and share of children changing more than the margin
    st <- d %>%
      group_by(metric) %>%
      summarise(ccc      = lin_ccc(raw, wins),
                pct_over = 100 * mean(abs(pct) > m),
                lo       = min(raw, wins),
                hi       = max(raw, wins),
                .groups = "drop") %>%
      arrange(metric) %>%
      mutate(metric_lab = paste0(metric, "\nCCC ", sprintf("%.3f", ccc), " \u00b7 ",
                                 sprintf("%.0f", pct_over), "% > ", m, "%"))
    lab_levels <- st$metric_lab
    
    d <- d %>%
      left_join(st %>% select(metric, metric_lab), by = "metric") %>%
      mutate(metric_lab = factor(metric_lab, levels = lab_levels),
             tip = paste0("<b>", participant_id, "</b><br/>",
                          group, "<br/>",
                          metric, "<br/>",
                          "Original: ", sprintf("%.2f", raw), "<br/>",
                          "After winsorizing: ", sprintf("%.2f", wins), "<br/>",
                          "Change: ", sprintf("%+.1f", pct), "%"))
    
    # Same range on both axes in every panel so the diagonal means "no change"
    lim <- bind_rows(st %>% transmute(metric_lab = factor(metric_lab, levels = lab_levels), x = lo, y = lo),
                     st %>% transmute(metric_lab = factor(metric_lab, levels = lab_levels), x = hi, y = hi))
    
    p <- ggplot(d, aes(x = raw, y = wins)) +
      geom_blank(data = lim, aes(x = x, y = y), inherit.aes = FALSE) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      geom_point_interactive(aes(color = group, tooltip = tip, data_id = participant_id),
                             alpha = 0.85, size = box_point_size) +
      facet_wrap(~ metric_lab, scales = "free", ncol = 2) +
      scale_color_manual(values = group_colors) +
      labs(x = "Original value", y = "After winsorizing") +
      plot_theme_bw() +
      theme(legend.position = "top",
            legend.title    = element_blank(),
            legend.text     = element_text(size = body_size_px)) +
      guides(color = guide_legend(nrow = 2))
    
    girafe(
      ggobj = p,
      width_svg = 7,
      height_svg = 6,
      options = list(
        opts_hover(css = "stroke:black;stroke-width:2px;fill-opacity:1;stroke-opacity:1;"),
        opts_hover_inv(css = "opacity:0.12;"),
        opts_tooltip(css = "background:white;color:black;border:1px solid #666;padding:6px 10px;font-family:Arial,Helvetica,sans-serif;font-size:14px;border-radius:4px;"),
        opts_toolbar(saveaspng = FALSE),
        opts_sizing(rescale = TRUE)
      )
    )
  })
  
  # ===========================================================================
  # TAB 2: Glucotypes (PCA + k-means)   (added)
  # ===========================================================================
  
  # Selected fit (raw or winsorized readings, depending on the switch)
  pc_sel <- reactive({
    req(input$pc_var, input$pc_k)
    if (isTRUE(input$pc_wins)) {
      fl    <- as.numeric(input$pc_floor)
      feats <- withProgress(pc_get(fl), message = "Recalculating measures...")
    } else {
      feats <- pc_get(NULL)
    }
    pc_fit(feats, input$pc_var, input$pc_k)
  })
  
  # Same settings on the raw readings (used to see who changes cluster)
  pc_raw <- reactive({
    req(input$pc_var, input$pc_k)
    pc_fit(pc_get(NULL), input$pc_var, input$pc_k)
  })
  
  # Shared ggiraph options (same look as the earlier plots)
  pc_girafe <- function(p, w, h) {
    girafe(
      ggobj = p,
      width_svg = w,
      height_svg = h,
      options = list(
        opts_hover(css = "stroke:black;stroke-width:2px;fill-opacity:1;stroke-opacity:1;"),
        opts_hover_inv(css = "opacity:0.25;"),
        opts_tooltip(css = pc_tip_css),
        opts_toolbar(saveaspng = FALSE),
        opts_sizing(rescale = TRUE)
      )
    )
  }
  
  # Size of in-plot text labels, matched to the axis text
  lab_size <- body_size_px / ggplot2::.pt * 0.85
  
  # ---- Caption above the PCA map --------------------------------------------
  output$pc_caption <- renderUI({
    f <- pc_sel()
    div(class = "intro-text", style = "margin-bottom:10px;",
        sprintf("PCA squeezes ten glucose measures into two summary patterns: %s (left to right) and %s (bottom to top). Each dot is one child, and dots close together have similar glucose patterns. Color is the glucotype the computer found. Shape is the child's true group. Where colors overlap, children sit between glucotypes. Grouping uses %d component%s; the map shows the first two.",
                tolower(f$lab1), tolower(f$lab2), f$ncomp, ifelse(f$ncomp == 1, "", "s")))
  })
  
  # ---- 1. Glucotype map -------------------------------------------------------
  output$pc_map <- renderGirafe({
    f <- pc_sel()
    d <- f$data %>%
      mutate(tip = paste0("<b>", participant_id, "</b><br/>",
                          "Group: ", group, "<br/>",
                          "Glucotype: ", cluster, "<br/>",
                          "Age: ", age, " years<br/>",
                          "Body Mass Index: ", bmi, "<br/>",
                          "Mean glucose: ", round(mean_glu), " mg/dL"))
    cols <- setNames(cluster_colors[seq_len(f$k)], levels(d$cluster))
    
    p <- ggplot(d, aes(x = pc1, y = pc2)) +
      geom_hline(yintercept = 0, color = "grey85") +
      geom_vline(xintercept = 0, color = "grey85") +
      geom_point_interactive(aes(color = cluster, shape = group,
                                 data_id = participant_id, tooltip = tip),
                             alpha = 0.85, size = box_point_size) +
      scale_color_manual(values = cols,
                         labels = setNames(f$prof$name, as.character(f$prof$cluster))) +
      scale_shape_manual(values = group_shapes) +
      labs(x = f$axis_x, y = f$axis_y) +
      plot_theme_bw() +
      theme(legend.position = "top",
            legend.box      = "vertical",
            legend.title    = element_blank(),
            legend.text     = element_text(size = body_size_px))
    
    girafe(
      ggobj = p,
      width_svg = 12,
      height_svg = 6.5,
      options = list(
        opts_hover(css = "stroke:black;stroke-width:2px;fill-opacity:1;stroke-opacity:1;"),
        opts_hover_inv(css = "opacity:0.12;"),
        opts_tooltip(css = pc_tip_css),
        opts_toolbar(saveaspng = FALSE),
        opts_sizing(rescale = TRUE)
      )
    )
  })
  
  # Does winsorizing change who lands in which glucotype?
  output$pc_switch_text <- renderText({
    if (!isTRUE(input$pc_wins)) {
      return("Turn on \"Use winsorized readings\" in the settings to see whether replacing false lows changes which glucotype each child lands in.")
    }
    a  <- pc_raw()$data
    b  <- pc_sel()$data
    ch <- sum(as.character(a$cluster) != as.character(b$cluster))
    sprintf("After winsorizing at %s mg/dL, %d of %d children (%.0f%%) land in a different glucotype than with the raw readings (matched by average glucose level).",
            input$pc_floor, ch, nrow(a), 100 * ch / nrow(a))
  })
  
  # ---- 2. Who is in each cluster --------------------------------------------
  output$pc_comp_plot <- renderGirafe({
    f    <- pc_sel()
    prof <- f$prof
    comp <- f$data %>%
      count(cluster, group) %>%
      group_by(cluster) %>%
      mutate(pct = 100 * n / sum(n)) %>%
      ungroup() %>%
      mutate(tip = paste0("<b>", cluster, "</b><br/>", group, "<br/>",
                          sprintf("%.0f", pct), "% of the glucotype (", n, " children)"),
             lab = ifelse(pct >= 8, paste0(round(pct), "%"), ""))
    x_labs <- pc_xlabs(prof)
    
    p <- ggplot(comp, aes(x = cluster, y = pct, fill = group)) +
      geom_col_interactive(aes(tooltip = tip, data_id = paste(cluster, group)),
                           width = 0.7, color = "white") +
      geom_text(aes(label = lab), position = position_stack(vjust = 0.5),
                color = "white", fontface = "bold", size = lab_size, family = plot_font) +
      scale_fill_manual(values = group_colors) +
      scale_x_discrete(labels = x_labs) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = NULL, y = "% of children in the glucotype") +
      plot_theme_bw() +
      theme(legend.position = "top",
            legend.title    = element_blank(),
            legend.text     = element_text(size = body_size_px)) +
      guides(fill = guide_legend(nrow = 2))
    
    pc_girafe(p, 7, 6.5)
  })
  
  output$pc_comp_text <- renderText({
    sprintf("%.0f%% of children landed in the glucotype dominated by their own group.",
            100 * pc_purity(pc_sel()$data))
  })
  
  # ---- 4. Time in ranges ----------------------------------------------------
  output$pc_range_plot <- renderGirafe({
    f    <- pc_sel()
    prof <- f$prof
    r <- f$data %>%
      group_by(cluster) %>%
      summarise(`Hypoglycemia`         = mean(pct_low),
                `In range`             = mean(pct_inr),
                `Hyperglycemia`        = mean(pct_high),
                `Severe hyperglycemia` = mean(pct_sev),
                .groups = "drop")
    long <- bind_rows(lapply(range_levels, function(l) {
      tibble(cluster = r$cluster, range = l, pct = r[[l]])
    })) %>%
      mutate(range = factor(range, levels = rev(range_levels)),
             tip   = paste0("<b>", cluster, "</b><br/>", range, ": ", sprintf("%.1f", pct), "%"),
             lab   = ifelse(pct >= 4, paste0(round(pct), "%"), ""))
    x_labs <- pc_xlabs(prof)
    
    p <- ggplot(long, aes(x = cluster, y = pct, fill = range)) +
      geom_col_interactive(aes(tooltip = tip, data_id = paste(cluster, range)),
                           width = 0.7, color = "white") +
      geom_text(aes(label = lab, group = range), position = position_stack(vjust = 0.5),
                color = "white", fontface = "bold", size = lab_size, family = plot_font) +
      scale_fill_manual(values = range_colors, breaks = range_levels, drop = FALSE) +
      scale_x_discrete(labels = x_labs) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = NULL, y = "Average % of time") +
      plot_theme_bw() +
      theme(legend.position = "top",
            legend.title    = element_blank(),
            legend.text     = element_text(size = body_size_px)) +
      guides(fill = guide_legend(nrow = 2))
    
    pc_girafe(p, 7, 6.5)
  })
  
  # ---- 5. Typical day per cluster -------------------------------------------
  output$pc_day_plot <- renderGirafe({
    f    <- pc_sel()
    prof <- f$prof
    feats <- if (isTRUE(input$pc_wins)) pc_get(as.numeric(input$pc_floor)) else pc_get(NULL)
    
    hr <- feats$hourly %>%
      left_join(f$data %>% select(participant_id, cluster), by = "participant_id") %>%
      group_by(cluster, hour) %>%
      summarise(med = median(glu),
                lo  = quantile(glu, 0.10),
                hi  = quantile(glu, 0.90),
                .groups = "drop") %>%
      left_join(prof %>% select(cluster, name), by = "cluster") %>%
      mutate(name = factor(name, levels = prof$name),
             hlab = ifelse(hour == 0, "12 AM",
                           ifelse(hour < 12, paste0(hour, " AM"),
                                  ifelse(hour == 12, "12 PM", paste0(hour - 12, " PM")))),
             tip  = paste0("<b>", name, "</b><br/>", hlab, "<br/>",
                           "Median: ", sprintf("%.0f", med), " mg/dL<br/>",
                           "Most children: ", sprintf("%.0f", lo), " to ",
                           sprintf("%.0f", hi), " mg/dL"))
    cols <- setNames(cluster_colors[seq_len(f$k)], levels(hr$cluster))
    
    p <- ggplot(hr, aes(x = hour)) +
      geom_hline(yintercept = c(70, 140), linetype = "dashed", color = "grey60") +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = cluster, group = cluster), alpha = 0.25) +
      geom_line(aes(y = med, color = cluster, group = cluster), linewidth = 1.1) +
      geom_point_interactive(aes(y = med, color = cluster, tooltip = tip,
                                 data_id = paste(cluster, hour)), size = 2.5) +
      facet_wrap(~ name, nrow = 1, labeller = label_wrap_gen(24)) +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      scale_x_continuous(breaks = c(0, 6, 12, 18),
                         labels = c("12 AM", "6 AM", "12 PM", "6 PM")) +
      labs(x = "Time of day", y = "Sensor glucose (mg/dL)") +
      plot_theme_bw()
    
    pc_girafe(p, 12, 5.5)
  })
  
  # ===========================================================================
  # NEW TAB 3: Glycemic Summary
  # ===========================================================================
  ov_data <- reactive({
    req(input$ov_groups, input$ov_days)
    d <- cgm_long %>%
      filter(day_num >= input$ov_days[1], day_num <= input$ov_days[2])
    d <- if (input$ov_pid == "all") {
      filter(d, group %in% input$ov_groups)
    } else {
      filter(d, participant_id == input$ov_pid)
    }
    if (input$ov_period != "all") d <- filter(d, period == input$ov_period)
    d
  })
  
  ov_child <- reactive({
    d <- ov_data()
    req(nrow(d) > 0)
    child_summary(d)
  })
  
  output$ov_n    <- renderText(n_distinct(ov_data()$participant_id))
  output$ov_mean <- renderText(sprintf("%.0f mg/dL", mean(ov_child()$mean_glu)))
  output$ov_gmi  <- renderText(sprintf("%.1f%%", mean(ov_child()$gmi)))
  output$ov_cv   <- renderText(sprintf("%.1f%%", mean(ov_child()$cv)))
  output$ov_tir  <- renderText(sprintf("%.0f%%", mean(ov_child()$pct_inr)))
  
  # Stacked horizontal bar: time in each range, by group
  output$ov_ranges <- renderGirafe({
    cs <- ov_child()
    r <- cs %>%
      group_by(group) %>%
      summarise(`Hypoglycemia`         = mean(pct_low),
                `In range`             = mean(pct_inr),
                `Hyperglycemia`        = mean(pct_high),
                `Severe hyperglycemia` = mean(pct_sev),
                .groups = "drop")
    long <- bind_rows(lapply(range_levels, function(l) {
      tibble(group = r$group, range = l, pct = r[[l]])
    })) %>%
      mutate(range = factor(range, levels = rev(range_levels)),
             tip   = paste0("<b>", group, "</b><br/>", range, ": ", sprintf("%.1f", pct), "%"),
             lab   = ifelse(pct >= 4, paste0(round(pct), "%"), ""))
    
    p <- ggplot(long, aes(x = pct, y = group, fill = range)) +
      geom_col_interactive(aes(tooltip = tip, data_id = paste(group, range)),
                           width = 0.65, color = "white") +
      geom_text(aes(label = lab, group = range), position = position_stack(vjust = 0.5),
                color = "white", fontface = "bold", size = lab_size_g, family = plot_font) +
      scale_fill_manual(values = range_colors, breaks = range_levels, drop = FALSE) +
      scale_y_discrete(labels = scales::label_wrap(16)) +
      scale_x_continuous(limits = c(0, 100)) +
      labs(x = "Average % of time", y = NULL) +
      plot_theme_bw() +
      theme(legend.position = "top",
            legend.title    = element_blank(),
            legend.text     = element_text(size = body_size_px)) +
      guides(fill = guide_legend(nrow = 2))
    
    girafe_std(p, 7, 5.5)
  })
  
  # Per-child distribution of the chosen measure
  output$ov_dist <- renderGirafe({
    cs  <- ov_child()
    met <- input$ov_metric
    lab <- names(ov_metric_choices)[ov_metric_choices == met]
    girafe_std(box_by_group(cs, cs[[met]], lab), 7, 5.5, inv = 0.15)
  })
  
  output$ov_table <- renderTable({
    ov_child() %>%
      group_by(Group = as.character(group)) %>%
      summarise(Children = n(),
                `Mean glucose (mg/dL)` = mean(mean_glu),
                `GMI (%)`              = mean(gmi),
                `CV (%)`               = mean(cv),
                `Time in range (%)`    = mean(pct_inr),
                `Time below 70 (%)`    = mean(pct_low),
                `Time below 54 (%)`    = mean(pct_low2),
                `Time above range (%)` = mean(pct_above),
                .groups = "drop")
  }, digits = 1, striped = TRUE, bordered = TRUE, width = "100%")
  
  # ===========================================================================
  # NEW TAB 4: Daily Patterns
  # ===========================================================================
  agp_data <- reactive({
    req(input$agp_groups)
    filter(cgm_long, group %in% input$agp_groups)
  })
  
  output$agp_plot <- renderGirafe({
    a <- agp_data() %>%
      group_by(group, bin) %>%
      summarise(q05 = quantile(glucose, 0.05), q25 = quantile(glucose, 0.25),
                q50 = median(glucose),
                q75 = quantile(glucose, 0.75), q95 = quantile(glucose, 0.95),
                .groups = "drop") %>%
      mutate(tip = paste0("<b>", group, "</b><br/>", fmt_time(bin), "<br/>",
                          "Median: ", sprintf("%.0f", q50), " mg/dL<br/>",
                          "Middle half: ", sprintf("%.0f", q25), " to ", sprintf("%.0f", q75), "<br/>",
                          "5th - 95th: ", sprintf("%.0f", q05), " to ", sprintf("%.0f", q95)))
    
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
      scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
                         labels = c("12 AM", "6 AM", "12 PM", "6 PM", "12 AM")) +
      coord_cartesian(ylim = c(40, 250)) +
      labs(x = "Time of day", y = "Sensor glucose (mg/dL)") +
      plot_theme_bw()
    
    girafe_std(p, 12, 5.5, inv = 0.35)
  })
  
  output$agp_heat <- renderGirafe({
    req(input$agp_pid)
    met <- input$agp_metric
    lab <- names(c("Mean glucose (mg/dL)" = "mean_glu", "% time below range" = "pct_low",
                   "% time in range" = "pct_inr"))[c("mean_glu", "pct_low", "pct_inr") == met]
    h <- cgm_long %>%
      filter(participant_id == input$agp_pid) %>%
      group_by(day_num, weekday, hour) %>%
      summarise(mean_glu = mean(glucose),
                pct_low  = 100 * mean(glucose <= 70),
                pct_inr  = 100 * mean(glucose > 70 & glucose < 140),
                .groups = "drop")
    h$value <- h[[met]]
    h$tip <- paste0("<b>Day ", h$day_num, " (", h$weekday, ")</b><br/>",
                    fmt_time(h$hour), "<br/>", lab, ": ", sprintf("%.1f", h$value))
    
    p <- ggplot(h, aes(x = day_num, y = hour, fill = value)) +
      geom_tile_interactive(aes(tooltip = tip, data_id = paste(day_num, hour)),
                            color = "white", linewidth = 0.2) +
      scale_y_reverse(breaks = c(0, 6, 12, 18),
                      labels = c("12 AM", "6 AM", "12 PM", "6 PM")) +
      scale_x_continuous(breaks = 1:14,
                         labels = paste0(1:14, "\n", wk_names[((1:14) - 1) %% 7 + 1]),
                         expand = c(0, 0)) +
      labs(x = "Day of monitoring", y = "Time of day") +
      plot_theme_bw() +
      theme(legend.position = "right")
    
    p <- if (met == "mean_glu") p + glu_fill_scale(lab) else
      p + scale_fill_viridis_c(name = lab, option = "plasma")
    
    girafe_std(p, 12, 5.5, inv = 0.4)
  })
  
  output$agp_week <- renderGirafe({
    w <- agp_data() %>%
      group_by(weekday, hour) %>%
      summarise(mean_glu = mean(glucose), .groups = "drop") %>%
      mutate(tip = paste0("<b>", weekday, "</b><br/>", fmt_time(hour), "<br/>",
                          "Mean glucose: ", sprintf("%.0f", mean_glu), " mg/dL"))
    
    p <- ggplot(w, aes(x = weekday, y = hour, fill = mean_glu)) +
      geom_tile_interactive(aes(tooltip = tip, data_id = paste(weekday, hour)),
                            color = "white", linewidth = 0.3) +
      glu_fill_scale("Mean glucose (mg/dL)") +
      scale_y_reverse(breaks = c(0, 6, 12, 18),
                      labels = c("12 AM", "6 AM", "12 PM", "6 PM")) +
      labs(x = NULL, y = "Time of day") +
      plot_theme_bw() +
      theme(legend.position = "right")
    
    girafe_std(p, 12, 5.5, inv = 0.4)
  })
  
  # ===========================================================================
  # NEW TAB 5: Variability & Events
  # ===========================================================================
  vr_floor_val <- reactive({
    if (isTRUE(input$vr_wins)) as.numeric(input$vr_floor) else NULL
  })
  
  vr_vars <- reactive({
    req(input$vr_groups)
    var_get(vr_floor_val()) %>% filter(group %in% input$vr_groups)
  })
  
  vr_eps <- reactive({
    req(input$vr_groups)
    ep_get(vr_floor_val()) %>% filter(group %in% input$vr_groups)
  })
  
  output$vr_n_ep      <- renderText({ req(nrow(vr_eps()) > 0); format(nrow(vr_eps()), big.mark = ",") })
  output$vr_per_child <- renderText({ req(nrow(vr_eps()) > 0); sprintf("%.1f", nrow(vr_eps()) / nrow(vr_vars())) })
  output$vr_dur       <- renderText({ req(nrow(vr_eps()) > 0); sprintf("%.0f min", median(vr_eps()$dur_min)) })
  output$vr_night     <- renderText({ req(nrow(vr_eps()) > 0); sprintf("%.0f%%", 100 * mean(vr_eps()$period == "Night")) })
  
  output$vr_box <- renderGirafe({
    v   <- vr_vars()
    req(nrow(v) > 0)
    met <- input$vr_metric
    lab <- names(var_metric_choices)[var_metric_choices == met]
    girafe_std(box_by_group(v, v[[met]], lab), 12, 5.5, inv = 0.15)
  })
  
  output$vr_onset <- renderGirafe({
    e <- vr_eps()
    req(nrow(e) > 0)
    n_child <- vr_vars() %>% count(group, name = "n_child")
    o <- e %>%
      mutate(hr = floor(onset_hour)) %>%
      count(group, hr, level) %>%
      left_join(n_child, by = "group") %>%
      mutate(per_child = n / n_child,
             tip = paste0("<b>", group, "</b><br/>", fmt_time(hr), " start<br/>",
                          level, "<br/>", n, " episodes (",
                          sprintf("%.2f", per_child), " per child)"))
    
    p <- ggplot(o, aes(x = hr, y = per_child, fill = level)) +
      geom_col_interactive(aes(tooltip = tip, data_id = paste(group, hr, level)), width = 0.85) +
      facet_wrap(~ group, ncol = 1) +
      scale_fill_manual(values = ep_level_cols, drop = FALSE) +
      scale_x_continuous(breaks = c(0, 6, 12, 18),
                         labels = c("12 AM", "6 AM", "12 PM", "6 PM")) +
      labs(x = "Hour the episode started", y = "Episodes per child (14 days)") +
      plot_theme_bw() +
      theme(legend.position = "top", legend.title = element_blank(),
            legend.text = element_text(size = body_size_px))
    
    girafe_std(p, 7, 7, inv = 0.4)
  })
  
  output$vr_scatter <- renderGirafe({
    e <- vr_eps()
    req(nrow(e) > 0)
    e <- e %>%
      mutate(tip = paste0("<b>", participant_id, "</b><br/>", group, "<br/>",
                          "Day ", day_num, ", ", fmt_time(onset_hour), "<br/>",
                          "Lasted ", dur_min, " min<br/>",
                          "Lowest reading: ", nadir, " mg/dL"))
    
    p <- ggplot(e, aes(x = dur_min, y = nadir, color = group)) +
      geom_hline(yintercept = 54, linetype = "dashed", color = "grey60") +
      geom_point_interactive(aes(tooltip = tip, data_id = participant_id),
                             alpha = 0.55, size = 2.5) +
      scale_color_manual(values = group_colors) +
      labs(x = "Episode duration (minutes)", y = "Lowest reading (mg/dL)") +
      plot_theme_bw() +
      theme(legend.position = "top", legend.title = element_blank(),
            legend.text = element_text(size = body_size_px)) +
      guides(color = guide_legend(nrow = 2, override.aes = list(alpha = 1)))
    
    girafe_std(p, 7, 7, inv = 0.1)
  })
  
  output$vr_trace <- renderPlot({
    req(input$vr_pid, input$vr_days)
    fl <- vr_floor_val()
    tr <- all_traces[[input$vr_pid]]
    g  <- if (is.null(fl)) tr$glucose else pmax(tr$glucose, fl)
    d0 <- input$vr_days[1] - 1
    d1 <- input$vr_days[2]
    df <- tibble(x = tr$day, glucose = g) %>% filter(x >= d0, x < d1)
    
    e <- ep_get(fl) %>%
      filter(participant_id == input$vr_pid) %>%
      mutate(xmin = (start - 1) * 5 / 1440, xmax = end * 5 / 1440) %>%
      filter(xmax > d0, xmin < d1)
    
    p <- ggplot(df, aes(x = x, y = glucose))
    if (nrow(e) > 0) {
      p <- p + geom_rect(data = e, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = level),
                         inherit.aes = FALSE, alpha = 0.35)
    }
    p +
      geom_vline(xintercept = d0:d1, color = "grey88") +
      geom_hline(yintercept = c(54, 70), linetype = "dashed", color = "grey50") +
      geom_line(linewidth = 0.6, color = "grey20") +
      scale_fill_manual(values = ep_level_cols, drop = FALSE, name = NULL) +
      scale_x_continuous(breaks = d0:(d1 - 1) + 0.5, labels = (d0 + 1):d1) +
      coord_cartesian(xlim = c(d0, d1), ylim = c(40, 250), expand = FALSE) +
      labs(x = "Day of monitoring", y = "Sensor glucose (mg/dL)") +
      plot_theme_bw() +
      theme(legend.position = "top", legend.text = element_text(size = body_size_px))
  })
  
  # ===========================================================================
  # NEW TAB 6: Modeling
  # ===========================================================================
  md_model <- eventReactive(input$md_fit, {
    validate(need(length(input$md_terms) > 0, "Choose at least one predictor."))
    sel <- input$md_terms
    if ("grp_night" %in% sel) sel <- union(sel, c("group", "night"))
    term_map <- c(group = "group", night = "night", circ = "sin_h + cos_h",
                  weekend = "weekend", age = "age", bmi = "bmi", sex = "sex")
    rhs <- unname(term_map[intersect(names(term_map), sel)])
    if ("grp_night" %in% input$md_terms) rhs <- c(rhs, "group:night")
    f <- as.formula(paste(input$md_outcome, "~", paste(rhs, collapse = " + ")))
    use_ar <- isTRUE(input$md_ar)
    
    fit <- withProgress(message = "Fitting mixed-effects model...", value = 0.3, {
      tryCatch(
        if (use_ar) {
          nlme::lme(f, random = ~ 1 | participant_id,
                    correlation = nlme::corAR1(form = ~ t | participant_id),
                    data = hourly_tbl, method = "REML")
        } else {
          nlme::lme(f, random = ~ 1 | participant_id,
                    data = hourly_tbl, method = "REML")
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
    fit <- m$fit
    vc  <- nlme::VarCorr(fit)
    sds <- suppressWarnings(as.numeric(vc[, "StdDev"]))
    phi <- if (m$ar) as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE)) else NA_real_
    out_lab <- names(md_outcome_choices)[md_outcome_choices == m$outcome]
    
    badge <- function(txt) {
      span(txt, style = "display:inline-block; padding:6px 12px; border-radius:6px; font-weight:bold; background:#3B1A78; color:white;")
    }
    badges <- list(
      badge(sprintf("Between-child SD: %.1f", sds[1])),
      badge(sprintf("Residual SD: %.1f", sds[2])),
      if (m$ar) badge(sprintf("AR(1) correlation: %.2f", phi)),
      badge(sprintf("AIC: %.0f", AIC(fit)))
    )
    
    tagList(
      p(class = "intro-text",
        sprintf("%s was modeled on %s hourly values from %d children (14 days each), with a separate baseline for every child%s.",
                out_lab, format(nrow(hourly_tbl), big.mark = ","),
                n_distinct(hourly_tbl$participant_id),
                if (m$ar) " and correlated (AR(1)) errors" else " and independent errors")),
      if (!m$ar) div(class = "text-danger", style = "margin-bottom:8px;",
                     "Warning: hourly values from the same child are strongly correlated, so with independent errors the intervals below are too narrow. Compare the ACF plot with and without AR(1)."),
      div(style = "display:flex; flex-wrap:wrap; gap:10px; margin:8px 0;", badges)
    )
  })
  
  output$md_forest <- renderGirafe({
    m  <- md_model()
    tt <- summary(m$fit)$tTable
    d <- tibble(term = rownames(tt),
                est  = tt[, "Value"],
                se   = tt[, "Std.Error"],
                dfr  = tt[, "DF"],
                pval = tt[, "p-value"]) %>%
      filter(term != "(Intercept)") %>%
      mutate(lab  = clean_term(term),
             crit = qt(0.975, dfr),
             lo   = est - crit * se,
             hi   = est + crit * se,
             sig  = lo > 0 | hi < 0,
             tip  = paste0("<b>", lab, "</b><br/>",
                           "Estimate: ", sprintf("%.2f", est), "<br/>",
                           "95% CI: ", sprintf("%.2f", lo), " to ", sprintf("%.2f", hi), "<br/>",
                           "p = ", format.pval(pval, digits = 2, eps = 0.001)))
    d$lab <- factor(d$lab, levels = rev(d$lab))
    
    p <- ggplot(d, aes(x = est, y = lab, xmin = lo, xmax = hi, color = sig)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      geom_pointrange_interactive(aes(tooltip = tip, data_id = lab), linewidth = 0.9, size = 0.8) +
      scale_color_manual(values = c("TRUE" = "#3B1A78", "FALSE" = "grey55")) +
      labs(x = "Estimated difference (mg/dL, 95% CI)", y = NULL) +
      plot_theme_bw()
    
    girafe_std(p, 7, 5.5, inv = 0.3)
  })
  
  output$md_acf <- renderPlot({
    m <- md_model()
    r <- as.numeric(residuals(m$fit, type = "normalized"))
    by_id <- split(r, hourly_tbl$participant_id)
    ac <- sapply(by_id, function(x) acf(x, lag.max = 24, plot = FALSE)$acf[, 1, 1])
    df <- tibble(lag = 1:24, acf = rowMeans(ac)[-1])
    ci <- 1.96 / sqrt(14 * 24)
    
    ggplot(df, aes(x = lag, y = acf)) +
      geom_hline(yintercept = 0, color = "grey40") +
      geom_hline(yintercept = c(-ci, ci), linetype = "dashed", color = "#B22222") +
      geom_col(fill = "#3B1A78", width = 0.6) +
      scale_x_continuous(breaks = c(1, 6, 12, 18, 24)) +
      coord_cartesian(ylim = c(-0.2, 1)) +
      labs(x = "Lag (hours)", y = "Average residual autocorrelation") +
      plot_theme_bw()
  })
  
  output$md_resid <- renderPlot({
    m <- md_model()
    set.seed(1)
    idx <- sample(nrow(hourly_tbl), 3000)
    df <- tibble(fitted = as.numeric(fitted(m$fit))[idx],
                 resid  = as.numeric(residuals(m$fit, type = "response"))[idx])
    ggplot(df, aes(x = fitted, y = resid)) +
      geom_hline(yintercept = 0, color = "grey40") +
      geom_point(alpha = 0.25, size = 1.5, color = "#3B1A78") +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#B22222") +
      labs(x = "Fitted value (mg/dL)", y = "Residual (mg/dL)") +
      plot_theme_bw()
  })
  
  # Circadian GAM (fit once, the first time it is needed)
  gam_res <- reactive({
    withProgress(message = "Fitting circadian GAM...", {
      hh <- hourly_tbl %>%
        group_by(participant_id, group, hour) %>%
        summarise(glu = mean(mean_glu), .groups = "drop") %>%
        mutate(hour = hour + 0.5) %>%
        as.data.frame()
      fit <- mgcv::gam(glu ~ group + s(hour, by = group, bs = "cc", k = 8) +
                         s(participant_id, bs = "re"),
                       data = hh, method = "REML", knots = list(hour = c(0, 24)))
      nd <- expand.grid(hour = seq(0, 24, length.out = 97),
                        group = factor(group_levels, levels = group_levels))
      nd$participant_id <- hh$participant_id[1]
      pr <- predict(fit, nd, se.fit = TRUE, exclude = "s(participant_id)")
      nd$fit <- as.numeric(pr$fit)
      nd$lo  <- nd$fit - 1.96 * as.numeric(pr$se.fit)
      nd$hi  <- nd$fit + 1.96 * as.numeric(pr$se.fit)
      nd
    })
  })
  
  output$md_gam <- renderPlot({
    nd <- gam_res()
    ggplot(nd, aes(x = hour, y = fit, color = group, fill = group)) +
      geom_hline(yintercept = c(70, 140), linetype = "dashed", color = "grey60") +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors) +
      scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
                         labels = c("12 AM", "6 AM", "12 PM", "6 PM", "12 AM")) +
      labs(x = "Time of day", y = "Mean glucose (mg/dL)") +
      plot_theme_bw() +
      theme(legend.position = "top", legend.title = element_blank(),
            legend.text = element_text(size = body_size_px)) +
      guides(color = guide_legend(nrow = 2))
  })
  
  output$md_text <- renderPrint({
    summary(md_model()$fit)
  })
  
  # ===========================================================================
  # NEW TAB 7: Forecasting
  # ===========================================================================
  fc_floor_val <- reactive({
    if (isTRUE(input$fc_wins)) as.numeric(input$fc_floor) else NULL
  })
  fc_hor <- reactive(as.integer(input$fc_h))
  
  fc_preds <- reactive({
    req(input$fc_pid, input$fc_h)
    fc_run(all_traces[[input$fc_pid]], fc_hor(), fc_floor_val())
  })
  
  output$fc_trace <- renderPlot({
    p  <- fc_preds()
    dd <- as.integer(input$fc_day)
    p  <- p[p$day == dd, ]
    validate(need(nrow(p) > 0, "No forecasts for this day."))
    
    all_names <- c("Actual", unname(fc_models))
    lines <- bind_rows(
      tibble(xh = p$xh, pred = p$actual, model = "Actual"),
      fc_long(p) %>% select(xh, pred, model)
    ) %>%
      mutate(model = factor(model, levels = all_names))
    
    ggplot() +
      geom_hline(yintercept = c(70, 140), linetype = "dashed", color = "grey60") +
      geom_ribbon(data = p, aes(x = xh, ymin = ar - 1.96 * sigma, ymax = ar + 1.96 * sigma),
                  fill = fc_cols[["AR + time of day"]], alpha = 0.15) +
      geom_line(data = lines, aes(x = xh, y = pred, color = model,
                                  linewidth = model == "Actual")) +
      scale_color_manual(values = fc_cols, breaks = all_names) +
      scale_linewidth_manual(values = c("TRUE" = 1.2, "FALSE" = 0.8), guide = "none") +
      scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
                         labels = c("12 AM", "6 AM", "12 PM", "6 PM", "12 AM")) +
      coord_cartesian(xlim = c(0, 24), ylim = c(40, 250)) +
      labs(x = paste0("Time of day (day ", dd, ")"), y = "Sensor glucose (mg/dL)") +
      plot_theme_bw() +
      theme(legend.position = "top", legend.title = element_blank(),
            legend.text = element_text(size = body_size_px)) +
      guides(color = guide_legend(override.aes = list(linewidth = 2)))
  })
  
  output$fc_table <- renderTable({
    fc_metric_tbl(fc_long(fc_preds()))
  }, digits = 1, striped = TRUE, bordered = TRUE, width = "100%")
  
  output$fc_folds <- renderGirafe({
    f <- fc_long(fc_preds()) %>%
      group_by(fold, model) %>%
      summarise(rmse = sqrt(mean((pred - actual)^2)), .groups = "drop") %>%
      mutate(model = factor(model, levels = unname(fc_models)),
             tip = paste0("<b>", model, "</b><br/>", fold, "<br/>RMSE: ",
                          sprintf("%.1f", rmse), " mg/dL"))
    
    p <- ggplot(f, aes(x = fold, y = rmse, color = model, group = model)) +
      geom_line(linewidth = 0.8) +
      geom_point_interactive(aes(tooltip = tip, data_id = paste(fold, model)), size = 3.5) +
      scale_color_manual(values = fc_cols) +
      labs(x = NULL, y = "RMSE (mg/dL)") +
      plot_theme_bw() +
      theme(legend.position = "top", legend.title = element_blank(),
            legend.text = element_text(size = body_size_px)) +
      guides(color = guide_legend(nrow = 2))
    
    girafe_std(p, 7, 5, inv = 0.3)
  })
  
  output$fc_roc <- renderPlot({
    p <- fc_preds()
    y <- p$actual <= 70
    validate(need(sum(y) >= 5 && sum(!y) >= 5,
                  "This child has too few low readings in the test days to draw a ROC curve."))
    curves <- lapply(names(fc_models), function(m) roc_curve(-p[[m]], y))
    labs_m <- vapply(seq_along(curves), function(i)
      sprintf("%s (AUC %.2f)", fc_models[[i]], curves[[i]]$auc), character(1))
    df <- bind_rows(lapply(seq_along(curves), function(i)
      mutate(curves[[i]]$df, model = labs_m[i]))) %>%
      mutate(model = factor(model, levels = labs_m))
    
    ggplot(df, aes(x = fpr, y = tpr, color = model)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = setNames(unname(fc_cols[unname(fc_models)]), labs_m)) +
      coord_equal() +
      labs(x = "False alarm rate (1 - specificity)", y = "Lows caught (sensitivity)") +
      plot_theme_bw() +
      theme(legend.position = "right", legend.title = element_blank(),
            legend.text = element_text(size = body_size_px))
  })
  
  # All 100 children (only when the button is clicked)
  fc_all_res <- eventReactive(input$fc_all, {
    h   <- fc_hor()
    fl  <- fc_floor_val()
    ids <- names(all_traces)
    res <- withProgress(message = "Forecasting all 100 children...", value = 0, {
      out <- lapply(ids, function(pid) {
        incProgress(1 / length(ids))
        p <- fc_run(all_traces[[pid]], h, fl)
        tibble(participant_id = pid,
               rmse_persist = sqrt(mean((p$persist - p$actual)^2)),
               rmse_ar      = sqrt(mean((p$ar - p$actual)^2)))
      })
      bind_rows(out)
    })
    res %>%
      left_join(part_info %>% select(participant_id, group), by = "participant_id") %>%
      mutate(horizon = h * 5)
  })
  
  output$fc_all_text <- renderText({
    r <- fc_all_res()
    sprintf("At a %d-minute horizon, the AR + time-of-day model had a lower RMSE than persistence for %d of %d children (median reduction %.0f%%).",
            r$horizon[1], sum(r$rmse_ar < r$rmse_persist), nrow(r),
            median(100 * (1 - r$rmse_ar / r$rmse_persist)))
  })
  
  output$fc_all_plot <- renderGirafe({
    r <- fc_all_res() %>%
      mutate(tip = paste0("<b>", participant_id, "</b><br/>", group, "<br/>",
                          "Persistence RMSE: ", sprintf("%.1f", rmse_persist), " mg/dL<br/>",
                          "AR + time of day RMSE: ", sprintf("%.1f", rmse_ar), " mg/dL"))
    lim <- range(c(r$rmse_persist, r$rmse_ar))
    
    p <- ggplot(r, aes(x = rmse_persist, y = rmse_ar)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      geom_point_interactive(aes(color = group, tooltip = tip, data_id = participant_id),
                             alpha = 0.85, size = box_point_size) +
      scale_color_manual(values = group_colors) +
      coord_equal(xlim = lim, ylim = lim) +
      labs(x = "Persistence RMSE (mg/dL)", y = "AR + time of day RMSE (mg/dL)") +
      plot_theme_bw() +
      theme(legend.position = "top", legend.title = element_blank(),
            legend.text = element_text(size = body_size_px)) +
      guides(color = guide_legend(nrow = 2))
    
    girafe_std(p, 7, 6, inv = 0.15)
  })
}

shinyApp(ui, server)
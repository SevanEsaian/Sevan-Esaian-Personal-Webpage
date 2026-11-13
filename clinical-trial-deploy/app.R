library(shiny); library(bslib); library(plotly); library(dplyr)
library(survival); library(gsDesign); library(DT)

col <- c(ctl = "#6B7C93", trt = "#2A9D8F", navy = "#1F3A5F", amber = "#E9A23B", red = "#C8553D")
ep_lab <- c(tte = "PFS (time-to-event)", bin = "ORR (binary)", cont = "FACT-L change (continuous)")
lay <- function(p, ...) p |> layout(font = list(family = "Inter, system-ui, sans-serif", size = 12),
                                    hoverlabel = list(bgcolor = "white"), legend = list(orientation = "h", y = -0.2), ...) |>
  config(displaylogo = FALSE)

# ---- Core statistical engine -------------------------------------------------
# P(event) under uniform accrual over A months, total study length Tt, exponential hazard l
pev <- function(l, A, Tt) 1 - (exp(-l * (Tt - A)) - exp(-l * Tt)) / (l * A)

# Expected events (time-to-event only), Schoenfeld framework
ev <- function(N, th, q) {
  if (q$ep != "tte") return(NA_real_ * N)
  l0 <- log(2) / q$m0; f <- q$r / (1 + q$r)
  N * (1 - q$drop) * (f * pev(l0 * th, q$A, q$Tt) + (1 - f) * pev(l0, q$A, q$Tt))
}

# Power at total N and true effect th (HR / treatment ORR / mean difference)
pw <- function(N, th, q) {
  f <- q$r / (1 + q$r)
  z <- switch(q$ep,
              cont = { Ne <- N * (1 - q$drop); abs(th) * sqrt(q$r * Ne) / (q$sd * (1 + q$r)) },
              bin  = { Ne <- N * (1 - q$drop); n1 <- Ne * f; n0 <- Ne - n1
              abs(th - q$p0) / sqrt(th * (1 - th) / n1 + q$p0 * (1 - q$p0) / n0) },
              tte  = abs(log(th)) * sqrt(ev(N, th, q) * q$r) / (1 + q$r))
  pnorm(z - q$za)
}
th_of <- function(q, f = 1) switch(q$ep, tte = exp(log(q$hr) * f),
                                   bin = q$p0 + (q$p1 - q$p0) * f, cont = q$delta * f)
nreq <- function(th, q, target) tryCatch(
  uniroot(function(N) pw(N, th, q) - target, c(4, 5e5))$root, error = function(e) NA_real_)

# Bayesian assurance: power averaged over a normal prior on the effect (as fraction of planned)
assur <- function(N, q, m, s) {
  fr <- seq(0, 2, length.out = 61); w <- dnorm(fr, m, s); w <- w / sum(w)
  sapply(N, function(n) sum(w * pw(n, th_of(q, fr), q)))
}

# Single simulated trial -> one-sided z statistic (positive = favours investigational)
simz <- function(q, th, N, mech) {
  n1 <- round(N * q$r / (1 + q$r)); n0 <- N - n1; d <- q$drop
  d0 <- if (mech == "mnar") d * .5 else d; d1 <- if (mech == "mnar") d * 1.5 else d
  if (q$ep == "cont") {
    y0 <- rnorm(n0, 0, q$sd); y1 <- rnorm(n1, th, q$sd)
    if (mech == "mnar") { k0 <- runif(n0) > pmin(.95, d0 * 2 * plogis(-y0 / q$sd))
    k1 <- runif(n1) > pmin(.95, d1 * 2 * plogis(-(y1 - th) / q$sd)) }
    else { k0 <- runif(n0) > d; k1 <- runif(n1) > d }
    a <- y1[k1]; b <- y0[k0]
    (mean(a) - mean(b)) / sqrt(var(a) / length(a) + var(b) / length(b))
  } else if (q$ep == "bin") {
    y0 <- rbinom(n0, 1, q$p0); y1 <- rbinom(n1, 1, th)
    k0 <- runif(n0) > ifelse(y0 == 1, d0 * .5, d0 * 1.5) * (mech == "mnar") - (mech != "mnar") * -d0 |> pmax(0)
    k1 <- runif(n1) > ifelse(y1 == 1, d1 * .5, d1 * 1.5) * (mech == "mnar") - (mech != "mnar") * -d1 |> pmax(0)
    x1 <- sum(y1[k1]); m1 <- sum(k1); x0 <- sum(y0[k0]); m0 <- sum(k0)
    pp <- (x1 + x0) / (m1 + m0)
    (x1 / m1 - x0 / m0) / sqrt(pp * (1 - pp) * (1 / m1 + 1 / m0))
  } else {
    l0 <- log(2) / q$m0; t <- c(rexp(n0, l0), rexp(n1, l0 * th)); arm <- rep(0:1, c(n0, n1))
    ld <- -log(1 - d) / q$Tt; dr <- c(rep(d0, n0), rep(d1, n1)) / max(d, 1e-9)
    cut <- pmin(q$Tt - runif(N, 0, q$A), rexp(N, ld * dr + 1e-12))
    s <- survdiff(Surv(pmin(t, cut), t <= cut) ~ arm)
    sign(s$exp[2] - s$obs[2]) * sqrt(s$chisq)
  }
}

# Economics ($M). z = list of cost/benefit assumptions
econ <- function(N, q, z) {
  mo <- N / (z$sites * z$rate) + z$fu
  cost <- (N * z$cpt + (N / (1 - z$sf) - N) * z$csf + z$sites * z$su + z$mon * mo) / 1000
  a <- assur(N, q, z$m, z$s); ben <- a * z$papp * (z$V - z$delay * mo)
  data.frame(N = N, months = mo, cost = cost, assur = a, ben = ben, enpv = ben - cost,
             pow = pw(N, th_of(q), q))
}

est <- function(d, ep) {
  if (ep == "tte") { s <- summary(coxph(Surv(AVAL, EVENT) ~ ARM, d))$conf.int; c(s[1], s[3], s[4]) }
  else if (ep == "bin") { m <- glm(RESP ~ ARM, binomial, d); exp(c(coef(m)[2], confint.default(m)[2, ])) }
  else { m <- lm(CHG ~ ARM, d); c(coef(m)[2], confint(m)[2, ]) }
}

# ---- UI helpers --------------------------------------------------------------
css <- "
 body{font-family:Inter,system-ui,-apple-system,'Segoe UI',sans-serif;background:#F6F8FA}
 .navbar{background:#1F3A5F!important} .card{border:1px solid #E3E8EE;box-shadow:none;border-radius:8px}
 .card-header{background:#fff;font-weight:600;color:#1F3A5F} .sidebar h6{color:#1F3A5F;margin-top:.8rem;font-weight:600}
 table{font-size:.85rem} details summary{cursor:pointer;color:#1F3A5F;font-weight:600}
 .about-plain,.about-tech{padding:.6rem .9rem;border-radius:4px;margin-bottom:.7rem;font-size:.92rem}
 .about-plain{background:#EAF5F3;border-left:4px solid #2A9D8F}
 .about-tech{background:#EFF2F7;border-left:4px solid #1F3A5F}
 .about-q h5{color:#1F3A5F;font-weight:600} .about-hero h2{color:#1F3A5F;font-weight:700}
 .about-tech ul,.about-plain ul{margin-bottom:.3rem}"

crd <- function(title, ..., h = 520) card(full_screen = TRUE, min_height = h, card_header(title), ...)

simple_tbl <- function(head, rows) {
  tags$table(class = "table table-sm table-striped",
             tags$thead(tags$tr(lapply(head, tags$th))),
             tags$tbody(lapply(rows, function(r) tags$tr(lapply(r, tags$td)))))
}

glance <- simple_tbl(c("Item", "Detail"), list(
  c("Indication", "Previously untreated metastatic non-small cell lung cancer"),
  c("Design", "Phase III, randomized, double-blind, 2-arm parallel, superiority"),
  c("Arms", "Investigational agent + chemotherapy vs. placebo + chemotherapy"),
  c("Stratification", "PD-L1 (<1%, 1\u201349%, \u226550%), ECOG (0/1), histology"),
  c("Hypotheses", "H0: no difference between arms; H1: investigational arm superior (one-sided)"),
  c("Endpoints", "Primary and key secondary endpoints are selectable in the sidebar"),
  c("Testing order", "PFS, then OS, then ORR (fixed-sequence)")))

estimands <- simple_tbl(c("Endpoint", "Variable", "Intercurrent event: strategy", "Summary measure"), list(
  c("PFS", "Time to progression (BICR, RECIST 1.1) or death",
    "Subsequent therapy: treatment policy. Missed assessments: censor at last adequate scan", "Hazard ratio"),
  c("ORR", "Confirmed complete or partial response",
    "Discontinuation without response: composite (non-responder)", "Risk difference"),
  c("FACT-L", "Change from baseline at week 24",
    "Discontinuation: hypothetical. Rescue therapy: treatment policy", "Mean difference")))

# ---- About tab content -------------------------------------------------------
plain <- function(...) div(class = "about-plain", tags$b(), ...)
tech  <- function(...) div(class = "about-tech",  tags$b("Technical detail. "), ...)
q_card <- function(title, body, where)
  card(class = "about-q", card_body(h5(title), p(body), p(class = "text-muted small mb-0", "Where to look: ", where)))

glossary <- simple_tbl(c("Term", "Meaning"), list(
  c("PFS", "Progression-free survival: time from randomization until the cancer progresses or the patient dies, whichever comes first."),
  c("ORR", "Objective response rate: the proportion of patients whose tumours shrink enough to count as a confirmed complete or partial response."),
  c("FACT-L", "Functional Assessment of Cancer Therapy\u2013Lung: a patient-reported quality-of-life score; here, change from baseline at week 24."),
  c("Hazard ratio (HR)", "Relative rate of events on treatment versus control. HR 0.70 means a 30% lower instantaneous risk of progression or death."),
  c("Type I error (alpha)", "Probability of declaring a benefit when there is none. Oncology Phase III trials typically use a one-sided 2.5%."),
  c("Power", "Probability the trial succeeds if the drug truly works exactly as assumed in the plan."),
  c("Assurance", "Bayesian 'probability of success': power averaged over our uncertainty about the true effect. Usually lower than power when the plan is optimistic."),
  c("Events", "For time-to-event endpoints, statistical information comes from the number of events (progressions or deaths), not the number of patients."),
  c("Dropout / MCAR / MNAR", "Patients who leave early. MCAR = at random; MNAR = related to the outcome itself, which can bias results."),
  c("Alpha spending", "A rule for 'spending' the Type I error budget across interim looks so that the overall error stays at the nominal level."),
  c("Multiplicity", "When several hypotheses are tested, the chance of at least one false positive grows unless the testing procedure is adjusted."),
  c("Estimand", "A precise statement of what treatment effect is being estimated, including how events such as discontinuation or rescue therapy are handled (ICH E9(R1))."),
  c("ENPV", "Expected net present value: expected commercial benefit (probability-weighted) minus trial cost, in $ millions."),
  c("Tipping point", "How bad the outcomes of dropouts would have to be, relative to the observed data, before the trial conclusion flips from success to failure.")))

methods_tbl <- simple_tbl(c("Component", "Method implemented", "Key assumptions"), list(
  c("Events / sample size (PFS)", "Schoenfeld events formula; events converted to patients through the exponential event probability under uniform accrual and a fixed study end; N solved numerically (root-finding on the power function)", "Proportional hazards; exponential survival; 24-month uniform accrual; no crossover or non-administrative censoring beyond dropout"),
  c("Sample size (ORR)", "Two-sample z-test for proportions, unpooled variance, normal approximation", "Independent binary outcomes; dropouts removed before analysis"),
  c("Sample size (FACT-L)", "Two-sample z-test for a mean difference with common SD", "Normality; common variance; known SD"),
  c("Assurance", "Power integrated over a normal prior on the effect (as a fraction of the planned effect), 61-point grid on 0\u20132, weights renormalized", "Prior is symmetric and independent of the data; effect scaled on the log-HR (PFS), risk-difference (ORR) or mean-difference scale"),
  c("Operating characteristics", "Monte Carlo: Welch z (continuous), pooled z (binary), signed log-rank z (PFS); exact binomial interval on the rejection rate", "Dropout mechanism as selected; one-sided test at the chosen alpha"),
  c("Interim analyses", "gsDesign: Lan\u2013DeMets spending (O'Brien\u2013Fleming or Pocock); optional non-binding beta-spending futility (Hwang\u2013Shih\u2013DeCani, gamma = \u22122)", "Equally informative looks scaled to the fixed-design information; boundaries are not folded into the Monte Carlo or economics"),
  c("Multiplicity", "Fixed-sequence (running maximum of p-values), Holm, Hochberg, Bonferroni", "User-entered p-values; Hochberg assumes non-negative dependence"),
  c("Tipping point", "Closed-form attenuation of the expected z-statistic by (investigational-arm dropout) \u00d7 (penalty on imputed outcomes)", "A simplification of delta-adjusted multiple imputation, for design-stage intuition"),
  c("Economics", "Cost build-up (patients, screen failures, sites, monitoring) and duration from enrolment rate; ENPV = assurance \u00d7 P(approval | success) \u00d7 (value \u2212 delay cost \u00d7 duration) \u2212 cost", "Undiscounted; single approval scenario; linear delay cost"),
  c("Mock analysis", "Stratified Cox / Kaplan\u2013Meier, risk difference and odds ratio, ANCOVA, subgroup forest plots on a synthetic dataset", "Synthetic data generated from the sidebar assumptions, with PD-L1 and ECOG effects and informative missingness")))

about_ui <- tagList(
  card(class = "about-hero", card_body(
    h2("Phase III Oncology Trial Design Studio"),
    p(class = "lead", "An interactive workbench for planning a randomized Phase III trial in previously untreated metastatic non-small cell lung cancer. Change one assumption in the sidebar and every tab updates together: sample size, probability of success, simulated trial outcomes, interim boundaries, trial economics and a rehearsal of the final analysis."),
    p("It was built to show how statistical design decisions translate into patients, time, risk and money. Each section below gives a short plain-English description and a technical description for statisticians."))),
  
  layout_columns(col_widths = c(3, 3, 3, 3), fill = FALSE,
                 q_card("How many patients do we need?", "Sample size, required events and how sensitive they are to optimistic assumptions.", "Sample Size & Power"),
                 q_card("How likely is success, and how robust is it?", "Power, Bayesian assurance, simulated trials, interim looks, multiplicity and missing data.", "Sample Size & Power; Operating Characteristics"),
                 q_card("Is the design worth the investment?", "Cost, duration and expected value by sample size, with sensitivity analysis.", "Cost / Benefit"),
                 q_card("What will the results look like?", "A mock dataset with baseline table, primary analysis, subgroups and an analysis plan.", "Data & Analysis Plan")),
  
  card(card_header("Quick start"), card_body(
    tags$ol(
      tags$li("Pick the endpoint in the sidebar (PFS, ORR or FACT-L) and set the planned treatment effect."),
      tags$li("Set the design: Type I error, target power, dropout, allocation ratio and whether sample size is inflated for dropout."),
      tags$li("Set the Bayesian prior to say how confident you are that the planned effect is real."),
      tags$li("Read the headline cards on Study Design and Sample Size & Power."),
      tags$li("Stress-test the design under Operating Characteristics (click 'Run simulation')."),
      tags$li("Open Cost / Benefit to see whether a larger or smaller trial maximizes expected value."),
      tags$li("Finish on Data & Analysis Plan to preview the analysis outputs."))),
    p(class = "text-muted small px-3 pb-2", "The sidebar controls apply to every tab. Charts are interactive: hover for exact values, and use the expand icon on a card for a full-screen view.")),
  
  withMathJax(accordion(id = "about_acc", open = "sidebar",
                        accordion_panel("Sidebar: planning assumptions", value = "sidebar",
                                        plain("The sidebar is the control panel. It holds the assumptions behind the whole trial: how well the drug is expected to work, how rigorous the trial must be, and how many patients will drop out."),
                                        tech(
                                          tags$ul(
                                            tags$li(tags$b("Endpoint. "), "PFS uses a hazard ratio, control median PFS and total study duration (accrual is fixed at 24 months). ORR uses control and investigational response rates. FACT-L uses a mean difference and standard deviation."),
                                            tags$li(tags$b("Type I error and power. "), "One-sided alpha of 0.025, 0.05 or 0.01, and target power from 70% to 99%."),
                                            tags$li(tags$b("Dropout. "), "0\u201340%. A switch controls whether sample size is inflated to compensate."),
                                            tags$li(tags$b("Allocation. "), "1:1, 3:2 or 2:1 (investigational : control)."),
                                            tags$li(tags$b("Prior on the true effect. "), "A normal prior expressed as a percentage of the planned effect (mean and SD). It drives assurance, the economics and the robustness heatmap.")),
                                          "Required sample size is found by root-finding on the exact power function rather than a closed-form shortcut, so accrual, follow-up, dropout and allocation are all reflected.")),
                        
                        accordion_panel("Study Design", value = "design",
                                        plain("A one-page summary of the trial: headline numbers, what is being tested and how effects are defined. The four cards at the top are the quick answers: how many patients, how many events, how likely success is, and what effect was assumed."),
                                        tech(
                                          tags$ul(
                                            tags$li("Cards report total N with arm split, required events (PFS, Schoenfeld) or analyzable patients after dropout (ORR, FACT-L), assurance against the target power, and the planned effect."),
                                            tags$li("The estimand table follows ICH E9(R1), giving the variable, the strategy for each intercurrent event and the summary measure. PFS uses a treatment-policy strategy, ORR a composite strategy and FACT-L a hypothetical strategy for discontinuation."),
                                            tags$li("The testing hierarchy is fixed-sequence: PFS, then OS, then ORR, at one-sided 2.5%."))) ),
                        
                        accordion_panel("Sample Size & Power", value = "power",
                                        plain("Shows how many patients are needed, and what happens if the drug works less well than hoped. Power curves show the chance of success as enrolment grows, for several 'what if the effect is weaker' cases. The second chart shows how much the required sample size changes when key assumptions change."),
                                        tech(
                                          tags$ul(
                                            tags$li(tags$b("Power curves. "), "Power against total N for true effects at 100, 85, 70 and 55% of plan. Overlays show assurance and the target-power line. Effect scaling is on the log-HR scale for PFS, and on the risk-difference or mean-difference scale otherwise."),
                                            tags$li(tags$b("Required sample size. "), "Required N against the true effect, dropout rate or target power."),
                                            tags$li(tags$b("Events (PFS). "), "Schoenfeld: \\(d = (1+r)^2 (z_{1-\\alpha} + z_{1-\\beta})^2 / (r \\log^2 HR)\\). Events are converted to patients with the exponential event probability under uniform accrual, \\(P(\\text{event}) = 1 - [e^{-\\lambda(T-A)} - e^{-\\lambda T}] / (\\lambda A)\\)."),
                                            tags$li(tags$b("Assurance. "), "\\(\\int \\text{Power}(\\theta)\\,\\pi(\\theta)\\,d\\theta\\) with a normal prior \\(\\pi\\) on the fraction of the planned effect, evaluated on a 61-point grid over 0\u20132 and renormalized. Mass at 0 represents a null drug.")))),
                        
                        accordion_panel("Operating Characteristics", value = "oc",
                                        plain("A stress test of the design. It runs thousands of virtual trials to check that the false-positive rate is controlled and that the trial wins as often as the maths says. It also shows how interim looks, multiple hypotheses and missing patient data change the picture."),
                                        tech(
                                          tags$ul(
                                            tags$li(tags$b("Monte Carlo. "), "Scenarios: planned effect (power), null (Type I error) or a 70% effect. Sample size can be varied from 50 to 150% of plan, with 200\u20133,000 simulated trials and a fixed seed. Dropout is either MCAR or informative and differential (MNAR: lower in control, higher in investigational, and related to outcome). Test statistics are Welch z (continuous), pooled z (binary) and signed log-rank z (PFS). Results are compared with analytic power, with an exact binomial interval on the simulated rejection rate."),
                                            tags$li(tags$b("Interim analyses. "), "gsDesign with 2\u20135 looks and Lan\u2013DeMets O'Brien\u2013Fleming or Pocock alpha spending, plus an optional non-binding futility boundary (Hwang\u2013Shih\u2013DeCani, gamma = \u22122). The plot and table show information fraction, cumulative events (PFS) or N, efficacy z, nominal p and cumulative alpha spent. Fixed-design requirements are scaled by the group-sequential inflation factor."),
                                            tags$li(tags$b("Multiplicity. "), "Enter p-values for PFS, OS and ORR. Fixed-sequence (gatekeeping), Holm, Hochberg and Bonferroni adjusted p-values and reject/retain decisions are shown side by side."),
                                            tags$li(tags$b("Missing data (tipping point). "), "Shows how the expected z-statistic erodes as outcomes imputed for investigational-arm dropouts are penalized. The penalty axis runs from 0 to 300%; 100% means no benefit is assumed after dropout (jump-to-reference). The tipping point is where the expected z falls to the critical value. It is a closed-form design-stage approximation, not a full delta-adjusted multiple imputation.")))),
                        
                        accordion_panel("Cost / Benefit", value = "cost",
                                        plain("A bigger trial is more likely to succeed but costs more and takes longer, which delays revenue. This tab finds the sample size with the best expected financial return and shows how robust that answer is."),
                                        tech(
                                          tags$ul(
                                            tags$li(tags$b("Duration and cost. "), "Duration = N / (sites \u00d7 enrolment rate) + follow-up. Cost covers per-patient cost, screen failures, site start-up and monitoring."),
                                            tags$li(tags$b("Expected value. "), "ENPV = assurance \u00d7 P(approval | success) \u00d7 (market value \u2212 delay cost \u00d7 duration) \u2212 trial cost, in $ millions and undiscounted. The value-maximizing N is taken from a 100-point grid from 0.3 to 2.5 times the planned N."),
                                            tags$li(tags$b("Probability of success. "), "Conventional power at the planned effect against Bayesian assurance across N."),
                                            tags$li(tags$b("Robustness heatmap. "), "ENPV across N while one input (prior mean, cost per patient, market value or delay cost) varies from 50 to 150% of its base value."),
                                            tags$li(tags$b("Design scenarios. "), "Side-by-side designs at 80, 85, 90 and 95% power: N, events, months, cost, assurance and ENPV.")))),
                        
                        accordion_panel("Data & Analysis Plan", value = "data",
                                        plain("A rehearsal of the final results. It generates a realistic fake dataset from your assumptions so you can see the tables and charts the analysis would produce, and the analysis plan that governs them. It contains no real patient data."),
                                        tech(
                                          tags$ul(
                                            tags$li(tags$b("Synthetic data. "), "Up to 1,500 subjects with ADaM-style variables (USUBJID, ARM, AVAL, EVENT, RESP, CHG). Covariates are age, sex, ECOG and PD-L1. PD-L1 modifies the treatment hazard ratio (\u00d71.12, \u00d71.00, \u00d70.88 across <1%, 1\u201349%, \u226550%) and ECOG 1 raises the hazard. Missingness is informative through baseline score. The seed is user-controlled."),
                                            tags$li(tags$b("Baseline table. "), "Means (SD) and counts (%) by arm."),
                                            tags$li(tags$b("Primary analysis. "), "PFS: Kaplan\u2013Meier curves with optional confidence bands, and a Cox model stratified by PD-L1 with a one-sided p-value. ORR: response rates with exact CIs, risk difference and odds ratio. FACT-L: ANCOVA adjusted for baseline and PD-L1."),
                                            tags$li(tags$b("Subgroups. "), "Forest plot of HR, OR or mean difference for PD-L1, ECOG, sex and age. These are exploratory."),
                                            tags$li(tags$b("Missing data. "), "Missingness of the final assessment by arm and PD-L1 subgroup."),
                                            tags$li(tags$b("Statistical analysis plan. "), "A summary table covering analysis sets, estimand, primary analysis, intercurrent events, multiplicity, interim analysis, missing data, sensitivity analyses and software. It changes with the selected endpoint.")))),
                        
                        accordion_panel("Methods at a glance", value = "methods",
                                        p("A compact reference for reviewers."), methods_tbl),
                        
                        accordion_panel("Glossary", value = "glossary", glossary),
                        
                        accordion_panel("Assumptions, limitations and intended use", value = "limits",
                                        plain("This is a planning and teaching tool. It is designed to build intuition and support early design discussions, not to replace a protocol, a statistical analysis plan or validated software."),
                                        tech(
                                          tags$ul(
                                            tags$li("PFS planning assumes proportional hazards and exponential survival with fixed 24-month uniform accrual. Delayed or non-proportional effects, as seen with some immunotherapies, are not modelled."),
                                            tags$li("Power and sample size use normal approximations. OS is not modelled as an endpoint; the multiplicity tab uses user-entered p-values."),
                                            tags$li("Interim boundaries are not incorporated into the Monte Carlo power or the economics, so those figures are fixed-design figures."),
                                            tags$li("The MNAR mechanism and the tipping-point calculation are illustrative simplifications."),
                                            tags$li("The economics are stylized: undiscounted, one approval scenario, linear delay cost."),
                                            tags$li("The dataset is synthetic. Subgroup results are exploratory and are not adjusted for multiplicity."),
                                            tags$li("Final designs should be confirmed in validated software (for example rpact, gsDesign or East) and agreed with the study team and regulators.")))
                        )))
)

# ---- UI ----------------------------------------------------------------------
ui <- page_navbar(
  title = "Phase III Oncology Trial Design Studio",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE,
  header = tags$style(HTML(css)),
  sidebar = sidebar(width = 320,
                    h6("Endpoint"),
                    selectInput("ep", NULL, c("PFS (time-to-event)" = "tte", "ORR (binary)" = "bin", "FACT-L change (continuous)" = "cont")),
                    conditionalPanel("input.ep == 'tte'",
                                     sliderInput("hr", "Hazard ratio", .4, .95, .70, .01), sliderInput("m0", "Control median PFS (months)", 3, 12, 6, .5),
                                     sliderInput("Tt", "Total study duration (months)", 30, 72, 42, 1)),
                    conditionalPanel("input.ep == 'bin'",
                                     sliderInput("p0", "Control ORR (%)", 10, 60, 30, 1), sliderInput("p1", "Investigational ORR (%)", 15, 80, 45, 1)),
                    conditionalPanel("input.ep == 'cont'",
                                     sliderInput("delta", "Mean difference (points)", 1, 10, 4, .5), sliderInput("sd", "Standard deviation", 6, 20, 12, .5)),
                    h6("Design"),
                    selectInput("alpha", "Type I error", c("One-sided 0.025" = .025, "One-sided 0.05" = .05, "One-sided 0.01" = .01)),
                    sliderInput("pw", "Target power (%)", 70, 99, 90, 1),
                    sliderInput("drop", "Dropout rate (%)", 0, 40, 15, 1),
                    selectInput("r", "Allocation (investigational : control)", c("1:1" = 1, "3:2" = 1.5, "2:1" = 2)),
                    input_switch("adj", "Inflate sample size for dropout", TRUE),
                    h6("Bayesian prior on true effect"),
                    sliderInput("pm", "Prior mean (% of planned effect)", 40, 120, 85, 5),
                    sliderInput("ps", "Prior SD (% of planned effect)", 5, 60, 25, 5)
  ),
  
  # ---- About ----
  nav_panel("About", about_ui),
  
  # ---- Tab 1: Study design ----
  nav_panel("Study Design",
            uiOutput("vb1"),
            layout_columns(col_widths = c(5, 7),
                           crd("Study at a glance", glance, h = 300),
                           crd("Estimands (ICH E9(R1))", estimands, h = 300))
  ),
  
  # ---- Tab 2: Sample size & power ----
  nav_panel("Sample Size & Power",
            uiOutput("vb2"),
            crd("Power curves",
                layout_sidebar(sidebar = sidebar(width = 230,
                                                 checkboxGroupInput("frac", "True effect (% of planned)", c("100%" = 100, "85%" = 85, "70%" = 70, "55%" = 55), c(100, 70)),
                                                 checkboxGroupInput("ovl", "Overlays", c("Assurance (Bayesian prior)" = "as", "Target power line" = "tl"), c("as", "tl"))),
                               plotlyOutput("pw1", height = "420px"))),
            crd("Required sample size",
                layout_sidebar(sidebar = sidebar(width = 230,
                                                 selectInput("x2", "Vary", c("True effect" = "eff", "Dropout rate" = "drop", "Target power" = "pow"))),
                               plotlyOutput("nq", height = "420px"))),
            crd("Methods and assumptions", h = 150,
                withMathJax(tags$details(
                  tags$summary("Show formulas"),
                  tags$p("Continuous: $n = 2\\sigma^2(z_{1-\\alpha}+z_{1-\\beta})^2/\\Delta^2$ per arm. Binary: unpooled normal approximation. Time-to-event (Schoenfeld): $d = (1+r)^2(z_{1-\\alpha}+z_{1-\\beta})^2 / (r\\,\\log^2 HR)$ events, converted to patients with exponential event probabilities under 24-month uniform accrual. Dropout inflates N by $1/(1-\\text{dropout})$. Assurance averages power over a normal prior on the true effect.")
                )))
  ),
  
  # ---- Tab 3: Operating characteristics ----
  nav_panel("Operating Characteristics",
            crd("Monte Carlo", h = 700,
                layout_sidebar(sidebar = sidebar(width = 250,
                                                 radioButtons("mcs", "Scenario", c("Planned effect (power)" = "alt", "Null (Type I error)" = "null", "Reduced effect (70%)" = "red")),
                                                 sliderInput("mcN", "Sample size (% of planned)", 50, 150, 100, 5),
                                                 sliderInput("nsim", "Simulated trials", 200, 3000, 1000, 100),
                                                 selectInput("mech", "Dropout mechanism", c("Random (MCAR)" = "mcar", "Informative, differential (MNAR)" = "mnar")),
                                                 actionButton("run", "Run simulation", class = "btn-primary")),
                               uiOutput("mcv"), plotlyOutput("mcp", height = "400px"))),
            crd("Interim analyses", h = 800,
                layout_sidebar(sidebar = sidebar(width = 250,
                                                 sliderInput("k", "Number of analyses", 2, 5, 3, 1),
                                                 selectInput("sfu", "Alpha-spending", c("O'Brien-Fleming (Lan-DeMets)" = "obf", "Pocock (Lan-DeMets)" = "poc")),
                                                 input_switch("fut", "Non-binding futility boundary", FALSE)),
                               plotlyOutput("gsp", height = "400px"), tableOutput("gst"))),
            crd("Multiplicity", h = 350,
                layout_sidebar(sidebar = sidebar(width = 250,
                                                 selectInput("mm", "Procedure", c("Fixed-sequence (gatekeeping)" = "seq", "Holm" = "holm", "Hochberg" = "hoch", "Bonferroni" = "bonf")),
                                                 numericInput("mp1", "p-value H1 (PFS)", .008, 0, 1, .001), numericInput("mp2", "p-value H2 (OS)", .02, 0, 1, .001),
                                                 numericInput("mp3", "p-value H3 (ORR)", .04, 0, 1, .001)),
                               tableOutput("mult"))),
            crd("Missing data",
                layout_sidebar(sidebar = sidebar(width = 250,
                                                 sliderInput("tf", "Dropout in investigational arm (%)", 5, 40, 20, 1)),
                               plotlyOutput("tip", height = "420px")))
  ),
  
  # ---- Tab 4: Cost / benefit ----
  nav_panel("Cost / Benefit",
            uiOutput("ecv"),
            layout_columns(col_widths = c(7, 5),
                           crd("Expected value by sample size", h = 620,
                               layout_sidebar(sidebar = sidebar(width = 280, open = "desktop",
                                                                checkboxGroupInput("ser", "Show", c("Expected NPV" = "enpv", "Expected benefit" = "ben", "Trial cost" = "cost"), c("enpv", "ben", "cost")),
                                                                sliderInput("cpt", "Cost per patient ($k)", 30, 250, 95, 5), sliderInput("su", "Site start-up ($k)", 10, 200, 70, 5),
                                                                sliderInput("sites", "Number of sites", 20, 300, 120, 5), sliderInput("rate", "Patients per site per month", .1, 2, .4, .05),
                                                                sliderInput("sf", "Screen-failure rate (%)", 0, 60, 25, 1), sliderInput("csf", "Cost per screen failure ($k)", 0, 20, 6, 1),
                                                                sliderInput("mon", "Monitoring ($k per month)", 100, 1500, 600, 25), sliderInput("fu", "Follow-up after last patient (mo)", 0, 36, 14, 1),
                                                                sliderInput("V", "Value if approved ($M)", 200, 6000, 2500, 100), sliderInput("papp", "P(approval | success) (%)", 50, 100, 85, 1),
                                                                sliderInput("delay", "Revenue lost per month of delay ($M)", 0, 40, 10, 1)),
                                              plotlyOutput("enpvp", height = "520px"))),
                           crd("Probability of success", h = 620, plotlyOutput("posp", height = "520px"))),
            layout_columns(col_widths = c(7, 5),
                           crd("Robustness heatmap",
                               layout_sidebar(sidebar = sidebar(width = 230,
                                                                selectInput("hy", "Vary", c("Prior mean of effect" = "m", "Cost per patient" = "cpt", "Market value" = "V", "Delay cost" = "delay"))),
                                              plotlyOutput("heat", height = "420px"))),
                           crd("Design scenarios", h = 300, DTOutput("scn")))
  ),
  
  # ---- Tab 5: Data & analysis plan ----
  nav_panel("Data & Analysis Plan",
            crd("Baseline", h = 380,
                layout_sidebar(sidebar = sidebar(width = 200, numericInput("seed", "Random seed", 2026, 1, 99999)), tableOutput("base"))),
            crd("Primary analysis", h = 700,
                layout_sidebar(sidebar = sidebar(width = 200, conditionalPanel("input.ep == 'tte'", input_switch("ci", "Show 95% CI", TRUE))),
                               plotlyOutput("prim", height = "420px"), tableOutput("primtab"))),
            crd("Subgroups",
                layout_sidebar(sidebar = sidebar(width = 200,
                                                 checkboxGroupInput("sg", "Subgroups", c("PD-L1" = "PDL1", "ECOG" = "ECOG", "Sex" = "SEX", "Age" = "AGEG"), c("PDL1", "ECOG"))),
                               plotlyOutput("forest", height = "420px"))),
            crd("Missing data", plotlyOutput("miss", height = "420px")),
            crd("Statistical analysis plan", h = 450, tableOutput("sap"))
  )
)

# ---- Server ------------------------------------------------------------------
server <- function(input, output, session) {
  
  q <- reactive({
    req(input$p1 > input$p0)
    list(ep = input$ep, za = qnorm(1 - as.numeric(input$alpha)), drop = input$drop / 100,
         r = as.numeric(input$r), hr = input$hr, m0 = input$m0, A = 24, Tt = input$Tt,
         p0 = input$p0 / 100, p1 = input$p1 / 100, sd = input$sd, delta = input$delta)
  })
  prior <- reactive(list(m = input$pm / 100, s = input$ps / 100))
  
  res <- reactive({
    qq <- q(); qp <- if (isTRUE(input$adj)) qq else modifyList(qq, list(drop = 0))
    N <- nreq(th_of(qq), qp, input$pw / 100)
    validate(need(is.finite(N), "Effect too small to reach the target power."))
    N <- ceiling(N); list(N = N, d = ev(N, th_of(qq), qq), qp = qp)
  })
  
  eff_txt <- function(qq) switch(qq$ep, tte = sprintf("HR %.2f, control median %g mo", qq$hr, qq$m0),
                                 bin = sprintf("ORR %.0f%% vs %.0f%%", 100 * qq$p1, 100 * qq$p0),
                                 cont = sprintf("Mean difference %g pts (SD %g)", qq$delta, qq$sd))
  
  vb <- function() {
    qq <- q(); r <- res(); n1 <- round(r$N * qq$r / (1 + qq$r)); pr <- prior()
    th <- function(x) value_box_theme(bg = x, fg = "white")
    layout_columns(fill = FALSE,
                   value_box("Total sample size", format(r$N, big.mark = ","),
                             sprintf("%d investigational / %d control", n1, r$N - n1), theme = th(col["navy"])),
                   value_box(if (qq$ep == "tte") "Required events" else "Analyzable patients",
                             format(round(if (qq$ep == "tte") r$d else r$N * (1 - qq$drop)), big.mark = ","),
                             if (qq$ep == "tte") "Schoenfeld approximation" else "After dropout", theme = th(col["trt"])),
                   value_box("Assurance (Bayesian)", sprintf("%.0f%%", 100 * assur(r$N, qq, pr$m, pr$s)),
                             sprintf("vs %.0f%% target power", input$pw), theme = th(col["ctl"])),
                   value_box("Planned effect", ep_lab[[qq$ep]], eff_txt(qq), theme = th(col["amber"])))
  }
  output$vb1 <- renderUI(vb()); output$vb2 <- renderUI(vb())
  
  # ---- Tab 2: sample size & power ---------------------------------------------
  output$pw1 <- renderPlotly({
    qq <- q(); r <- res(); Ns <- round(seq(r$N * .4, r$N * 1.6, length.out = 60))
    p <- plot_ly(); pal <- c("#1F3A5F", "#2A9D8F", "#E9A23B", "#C8553D")
    for (i in seq_along(input$frac)) {
      f <- as.numeric(input$frac[i]) / 100; th <- th_of(qq, f); pwr <- pw(Ns, th, qq)
      tx <- sprintf("N = %d<br>Power = %.1f%%<br>True effect = %.0f%% of planned%s", Ns, 100 * pwr, 100 * f,
                    if (qq$ep == "tte") sprintf("<br>Events = %.0f", ev(Ns, th, qq)) else "")
      p <- add_lines(p, x = Ns, y = pwr, name = paste0(100 * f, "% of planned effect"),
                     line = list(color = pal[i %% 4 + 1]), text = tx, hoverinfo = "text")
    }
    if ("as" %in% input$ovl) { pr <- prior(); a <- assur(Ns, qq, pr$m, pr$s)
    p <- add_lines(p, x = Ns, y = a, name = "Assurance", line = list(color = col["red"], dash = "dot"),
                   text = sprintf("N = %d<br>Assurance = %.1f%%", Ns, 100 * a), hoverinfo = "text") }
    shp <- list(list(type = "line", x0 = r$N, x1 = r$N, y0 = 0, y1 = 1, line = list(color = "grey", dash = "dash")))
    if ("tl" %in% input$ovl) shp[[2]] <- list(type = "line", x0 = min(Ns), x1 = max(Ns), y0 = input$pw / 100,
                                              y1 = input$pw / 100, line = list(color = "grey", dash = "dash"))
    lay(p, shapes = shp, xaxis = list(title = "Total sample size"),
        yaxis = list(title = "Power", tickformat = ".0%", range = c(0, 1.02)))
  })
  
  output$nq <- renderPlotly({
    qq <- q(); qp <- res()$qp; x2 <- input$x2
    if (x2 == "eff") { xs <- seq(.4, 1.3, .05); gs <- c(.8, .9, .95); xl <- "True effect (% of planned)"
    fn <- function(x, g) nreq(th_of(qq, x), qp, g); gl <- function(g) paste0(100 * g, "% power") }
    else if (x2 == "drop") { xs <- seq(0, .4, .02); gs <- c(.8, .9, .95); xl <- "Dropout rate"
    fn <- function(x, g) nreq(th_of(qq), modifyList(qp, list(drop = x)), g); gl <- function(g) paste0(100 * g, "% power") }
    else { xs <- seq(.7, .99, .01); gs <- c(1, .85, .7); xl <- "Target power"
    fn <- function(x, g) nreq(th_of(qq, g), qp, x); gl <- function(g) paste0(100 * g, "% of planned effect") }
    p <- plot_ly(); pal <- c("#1F3A5F", "#2A9D8F", "#E9A23B")
    for (i in seq_along(gs)) { y <- ceiling(sapply(xs, fn, g = gs[i]))
    p <- add_lines(p, x = xs, y = y, name = gl(gs[i]), line = list(color = pal[i]),
                   text = sprintf("%s = %s<br>Required N = %s<br>%s", xl, xs, format(y, big.mark = ","), gl(gs[i])), hoverinfo = "text") }
    lay(p, xaxis = list(title = xl, tickformat = if (x2 != "eff") ".0%" else ".0%"),
        yaxis = list(title = "Required total N"))
  })
  
  # ---- Tab 3: operating characteristics ----------------------------------------
  mc <- eventReactive(input$run, {
    qq <- q(); N <- round(res()$N * input$mcN / 100)
    th <- th_of(qq, switch(input$mcs, null = 0, alt = 1, red = .7)); set.seed(2026)
    withProgress(message = "Simulating trials", z <- replicate(input$nsim, simz(qq, th, N, input$mech)))
    list(z = z, N = N, an = pw(N, th, qq), qq = qq, null = input$mcs == "null")
  }, ignoreNULL = FALSE)
  
  output$mcv <- renderUI({
    m <- mc(); rej <- sum(m$z > m$qq$za); ci <- binom.test(rej, length(m$z))$conf.int
    th <- function(x) value_box_theme(bg = x, fg = "white")
    layout_columns(fill = FALSE,
                   value_box(if (m$null) "Simulated Type I error" else "Simulated power", sprintf("%.1f%%", 100 * rej / length(m$z)),
                             sprintf("95%% MC interval %.1f\u2013%.1f%%", 100 * ci[1], 100 * ci[2]), theme = th(col["navy"])),
                   value_box(if (m$null) "Nominal alpha" else "Analytic power", sprintf("%.1f%%", 100 * (if (m$null) 1 - pnorm(m$qq$za) else m$an)),
                             "Normal approximation", theme = th(col["trt"])),
                   value_box("Total N per trial", m$N, sprintf("%d simulated trials", length(m$z)), theme = th(col["ctl"])))
  })
  output$mcp <- renderPlotly({
    m <- mc()
    plot_ly(x = m$z, type = "histogram", nbinsx = 45, marker = list(color = col["trt"], line = list(color = "white", width = 1)),
            hovertemplate = "z in %{x}<br>Trials: %{y}<extra></extra>") |>
      lay(shapes = list(list(type = "line", x0 = m$qq$za, x1 = m$qq$za, y0 = 0, y1 = 1, yref = "paper",
                             line = list(color = col["red"], dash = "dash"))),
          xaxis = list(title = "Test statistic (z); dashed line = critical value"), yaxis = list(title = "Simulated trials"))
  })
  
  gsd <- reactive({
    gsDesign(k = input$k, test.type = if (input$fut) 4 else 1, alpha = as.numeric(input$alpha), beta = 1 - input$pw / 100,
             sfu = switch(input$sfu, obf = sfLDOF, poc = sfLDPocock), sfl = sfHSD, sflpar = -2, n.fix = 1)
  })
  gs_tab <- reactive({
    g <- gsd(); r <- res(); base <- if (q()$ep == "tte") r$d else r$N
    data.frame(Look = seq_len(g$k), `Information` = round(g$timing, 2),
               `Cumulative events / N` = ceiling(base * g$n.I), `Efficacy z` = round(g$upper$bound, 3),
               `Nominal p` = signif(1 - pnorm(g$upper$bound), 3), `Cumulative alpha` = round(cumsum(g$upper$spend), 4),
               check.names = FALSE)
  })
  output$gsp <- renderPlotly({
    g <- gsd(); tb <- gs_tab()
    p <- plot_ly() |> add_trace(x = g$timing, y = g$upper$bound, type = "scatter", mode = "lines+markers", name = "Efficacy",
                                line = list(color = col["navy"]), marker = list(color = col["navy"], size = 9),
                                text = sprintf("Look %d<br>Information %.0f%%<br>z = %.3f<br>Nominal p = %.4f<br>Cumulative alpha = %.4f<br>%s = %s",
                                               tb$Look, 100 * g$timing, g$upper$bound, tb$`Nominal p`, tb$`Cumulative alpha`,
                                               if (q()$ep == "tte") "Events" else "N", tb$`Cumulative events / N`), hoverinfo = "text")
    if (input$fut) p <- add_trace(p, x = g$timing, y = g$lower$bound, type = "scatter", mode = "lines+markers", name = "Futility",
                                  line = list(color = col["amber"]), marker = list(color = col["amber"], size = 9),
                                  text = sprintf("Look %d<br>Futility z = %.3f", tb$Look, g$lower$bound), hoverinfo = "text")
    lay(p, xaxis = list(title = "Information fraction", tickformat = ".0%"), yaxis = list(title = "Boundary (z)"))
  })
  output$gst <- renderTable(gs_tab(), striped = TRUE, width = "100%")
  
  output$mult <- renderTable({
    p <- c(input$mp1, input$mp2, input$mp3); a <- as.numeric(input$alpha)
    adj <- switch(input$mm, bonf = pmin(1, p * 3), holm = p.adjust(p, "holm"), hoch = p.adjust(p, "hochberg"), seq = pmin(1, cummax(p)))
    data.frame(Hypothesis = c("H1: PFS", "H2: OS", "H3: ORR"), `Raw p` = p, `Adjusted p` = round(adj, 4),
               Decision = ifelse(adj <= a, "Reject H0", "Retain H0"), check.names = FALSE)
  }, striped = TRUE, width = "100%")
  
  output$tip <- renderPlotly({
    f <- input$tf / 100; qq <- q(); zexp <- qnorm(input$pw / 100) + qq$za
    s <- seq(0, 3, .05); zz <- zexp * (1 - f * s); tp <- (1 - qq$za / zexp) / f
    p <- plot_ly(x = s, y = zz, type = "scatter", mode = "lines", line = list(color = col["navy"]),
                 text = sprintf("Imputed effect retained: %.0f%%<br>Adjusted z = %.2f<br>One-sided p = %.4f", 100 * (1 - s), zz, 1 - pnorm(zz)),
                 hoverinfo = "text", name = "Adjusted z")
    lay(p, shapes = list(list(type = "line", x0 = 0, x1 = 3, y0 = qq$za, y1 = qq$za, line = list(color = col["red"], dash = "dash"))),
        annotations = list(list(x = min(tp, 3), y = qq$za, text = if (tp < 3) sprintf("Tipping point: %.0f%% penalty", 100 * tp) else "No tipping point in range",
                                showarrow = TRUE, ay = -40)),
        xaxis = list(title = "Penalty applied to imputed outcomes in dropouts (1 = no effect, i.e. jump-to-reference)", tickformat = ".0%"),
        yaxis = list(title = "Expected z (red = critical value)"))
  })
  
  # ---- Tab 4: cost / benefit ---------------------------------------------------
  zz <- reactive({ pr <- prior()
  list(cpt = input$cpt, su = input$su, sites = input$sites, rate = input$rate, sf = input$sf / 100, csf = input$csf,
       mon = input$mon, fu = input$fu, V = input$V, papp = input$papp / 100, delay = input$delay, m = pr$m, s = pr$s) })
  curve <- reactive({ Np <- res()$N; Ns <- round(seq(.3 * Np, 2.5 * Np, length.out = 100)); econ(Ns, q(), zz()) })
  
  output$ecv <- renderUI({
    e <- curve(); o <- e[which.max(e$enpv), ]; pl <- econ(res()$N, q(), zz()); th <- function(x) value_box_theme(bg = x, fg = "white")
    layout_columns(fill = FALSE,
                   value_box("Value-maximizing N", format(o$N, big.mark = ","), sprintf("Power %.0f%%, assurance %.0f%%", 100 * o$pow, 100 * o$assur), theme = th(col["navy"])),
                   value_box("Expected NPV at optimum", sprintf("$%.0fM", o$enpv), sprintf("Cost $%.0fM, %.0f months", o$cost, o$months), theme = th(col["trt"])),
                   value_box("Expected NPV at planned N", sprintf("$%.0fM", pl$enpv), sprintf("N = %s (%.0f%% power)", format(res()$N, big.mark = ","), input$pw), theme = th(col["ctl"])),
                   value_box("Gain from re-sizing", sprintf("$%.0fM", o$enpv - pl$enpv), "Optimum vs. planned", theme = th(col["amber"])))
  })
  output$enpvp <- renderPlotly({
    e <- curve(); o <- which.max(e$enpv); p <- plot_ly()
    ser <- list(enpv = c("Expected NPV", col["navy"], "enpv"), ben = c("Expected benefit", col["trt"], "ben"), cost = c("Trial cost", col["red"], "cost"))
    for (k in input$ser) { s <- ser[[k]]
    p <- add_lines(p, x = e$N, y = e[[s[3]]], name = s[1], line = list(color = s[2]),
                   text = sprintf("N = %d<br>%s = $%.0fM<br>Duration = %.1f months<br>Assurance = %.1f%%", e$N, s[1], e[[s[3]]], e$months, 100 * e$assur), hoverinfo = "text") }
    if ("enpv" %in% input$ser) p <- add_markers(p, x = e$N[o], y = e$enpv[o], name = "Optimum", marker = list(color = col["amber"], size = 12),
                                                text = sprintf("Optimum N = %d<br>ENPV = $%.0fM", e$N[o], e$enpv[o]), hoverinfo = "text")
    lay(p, xaxis = list(title = "Total sample size"), yaxis = list(title = "$ millions"))
  })
  output$posp <- renderPlotly({
    e <- curve(); lay(plot_ly(e, x = ~N) |>
                        add_lines(y = ~pow, name = "Conventional power", line = list(color = col["navy"]),
                                  text = ~sprintf("N = %d<br>Power at planned effect = %.1f%%", N, 100 * pow), hoverinfo = "text") |>
                        add_lines(y = ~assur, name = "Assurance", line = list(color = col["red"], dash = "dot"),
                                  text = ~sprintf("N = %d<br>Assurance = %.1f%%", N, 100 * assur), hoverinfo = "text"),
                      xaxis = list(title = "Total sample size"), yaxis = list(title = "Probability of success", tickformat = ".0%", range = c(0, 1)))
  })
  output$heat <- renderPlotly({
    Np <- res()$N; Ns <- round(seq(.5 * Np, 2 * Np, length.out = 24)); key <- input$hy
    base <- zz()[[key]]; mult <- seq(.5, 1.5, length.out = 9); ys <- base * mult
    m <- sapply(mult, function(k) { z2 <- zz(); z2[[key]] <- base * k; econ(Ns, q(), z2)$enpv })
    lbl <- c(m = "Prior mean (fraction of planned effect)", cpt = "Cost per patient ($k)", V = "Market value ($M)", delay = "Delay cost ($M/mo)")[[key]]
    plot_ly(x = Ns, y = ys, z = t(m), type = "heatmap", zmid = 0, colors = colorRamp(c(col["red"], "white", col["trt"])),
            hovertemplate = paste0("N = %{x}<br>", lbl, " = %{y:.2f}<br>ENPV = $%{z:.0f}M<extra></extra>")) |>
      lay(xaxis = list(title = "Total sample size"), yaxis = list(title = lbl))
  })
  output$scn <- renderDT({
    qq <- q(); th <- th_of(qq)
    rows <- lapply(c(80, 85, 90, 95), function(p) { N <- ceiling(nreq(th, res()$qp, p / 100)); e <- econ(N, qq, zz())
    data.frame(`Target power` = paste0(p, "%"), N = N, Events = round(ev(N, th, qq)), Months = round(e$months, 1),
               `Cost ($M)` = round(e$cost, 1), `Assurance (%)` = round(100 * e$assur, 1), `ENPV ($M)` = round(e$enpv, 1), check.names = FALSE) })
    datatable(do.call(rbind, rows), rownames = FALSE, options = list(dom = "t"))
  })
  
  # ---- Tab 5: simulated analysis dataset --------------------------------------
  dat <- reactive({
    qq <- q(); N <- min(res()$N, 1500); n1 <- round(N * qq$r / (1 + qq$r)); set.seed(input$seed)
    d <- data.frame(USUBJID = sprintf("ONC-%04d", 1:N), ARM = sample(rep(c("Control", "Investigational"), c(N - n1, n1))),
                    AGE = round(rnorm(N, 64, 9)), SEX = sample(c("F", "M"), N, TRUE, c(.45, .55)), ECOG = sample(c("0", "1"), N, TRUE, c(.4, .6)),
                    PDL1 = sample(c("<1%", "1-49%", ">=50%"), N, TRUE, c(.3, .35, .35)), BASE = round(rnorm(N, 70, 12), 1))
    d$AGEG <- ifelse(d$AGE >= 65, ">=65", "<65"); trt <- d$ARM == "Investigational"
    hr <- ifelse(trt, qq$hr * unname(c("<1%" = 1.12, "1-49%" = 1, ">=50%" = .88)[d$PDL1]), 1) * exp(.25 * (d$ECOG == "1"))
    t <- rexp(N, log(2) / qq$m0 * hr); cut <- pmin(qq$Tt - runif(N, 0, qq$A), rexp(N, -log(1 - max(qq$drop, .001)) / qq$Tt))
    d$EVENT <- as.integer(t <= cut); d$AVAL <- round(pmin(t, cut), 2)
    d$RESP <- rbinom(N, 1, ifelse(trt, qq$p1, qq$p0))
    d$CHG <- round(rnorm(N, -2 + qq$delta * trt - .25 * (d$BASE - 70), qq$sd), 1)
    d$MISS <- runif(N) < plogis(qlogis(max(qq$drop, .001)) - .04 * (d$BASE - 70))
    d$CHG[d$MISS] <- NA; d$RESP[d$MISS] <- 0
    d$ARM <- factor(d$ARM, c("Control", "Investigational")); d
  })
  
  output$base <- renderTable({
    d <- dat(); ms <- function(x) sprintf("%.1f (%.1f)", mean(x), sd(x)); np <- function(x) sprintf("%d (%.0f%%)", sum(x), 100 * mean(x))
    rows <- list("Age, mean (SD)" = \(x) ms(x$AGE), "Female, n (%)" = \(x) np(x$SEX == "F"), "ECOG 1, n (%)" = \(x) np(x$ECOG == "1"),
                 "PD-L1 <1%, n (%)" = \(x) np(x$PDL1 == "<1%"), "PD-L1 1-49%, n (%)" = \(x) np(x$PDL1 == "1-49%"),
                 "PD-L1 >=50%, n (%)" = \(x) np(x$PDL1 == ">=50%"), "Baseline FACT-L, mean (SD)" = \(x) ms(x$BASE))
    out <- data.frame(Characteristic = names(rows))
    for (a in levels(d$ARM)) out[[sprintf("%s (N=%d)", a, sum(d$ARM == a))]] <- sapply(rows, \(f) f(d[d$ARM == a, ]))
    out
  }, striped = TRUE, width = "100%")
  
  output$prim <- renderPlotly({
    d <- dat(); ep <- q()$ep; cl <- c(Control = col[["ctl"]], Investigational = col[["trt"]])
    if (ep == "tte") {
      sf <- survfit(Surv(AVAL, EVENT) ~ ARM, d); s <- data.frame(t = sf$time, S = sf$surv, lo = sf$lower, hi = sf$upper, nr = sf$n.risk,
                                                                 a = rep(levels(d$ARM), sf$strata)); p <- plot_ly()
                                                                 for (a in levels(d$ARM)) { x <- s[s$a == a, ]
                                                                 if (isTRUE(input$ci)) p <- add_ribbons(p, data = x, x = ~t, ymin = ~lo, ymax = ~hi, fillcolor = cl[[a]], opacity = .15, line = list(width = 0), showlegend = FALSE, hoverinfo = "none")
                                                                 p <- add_lines(p, data = x, x = ~t, y = ~S, name = a, line = list(color = cl[[a]], shape = "hv"),
                                                                                text = ~sprintf("%s<br>Month %.1f<br>PFS = %.1f%%<br>At risk = %d", a, t, 100 * S, nr), hoverinfo = "text") }
                                                                 lay(p, xaxis = list(title = "Months since randomization"), yaxis = list(title = "Progression-free survival", tickformat = ".0%", range = c(0, 1)))
    } else if (ep == "bin") {
      s <- d |> group_by(ARM) |> summarise(x = sum(RESP), n = n(), .groups = "drop"); ci <- t(sapply(1:2, \(i) binom.test(s$x[i], s$n[i])$conf.int))
      plot_ly(s, x = ~ARM, y = ~x / n, type = "bar", marker = list(color = unname(cl[as.character(s$ARM)])),
              error_y = list(type = "data", symmetric = FALSE, array = ci[, 2] - s$x / s$n, arrayminus = s$x / s$n - ci[, 1]),
              text = ~sprintf("%s<br>ORR = %.1f%% (%d/%d)", ARM, 100 * x / n, x, n), hoverinfo = "text") |>
        lay(xaxis = list(title = ""), yaxis = list(title = "Objective response rate (95% CI)", tickformat = ".0%"), showlegend = FALSE)
    } else {
      plot_ly(d, x = ~ARM, y = ~CHG, color = ~ARM, colors = cl, type = "box", boxpoints = "all", jitter = .4, pointpos = 0, marker = list(opacity = .35),
              text = ~sprintf("%s<br>Baseline = %.1f<br>Change = %.1f", USUBJID, BASE, CHG), hoverinfo = "text") |>
        lay(xaxis = list(title = ""), yaxis = list(title = "Change in FACT-L from baseline"), showlegend = FALSE)
    }
  })
  output$primtab <- renderTable({
    d <- dat(); ep <- q()$ep
    if (ep == "tte") { sf <- survfit(Surv(AVAL, EVENT) ~ ARM, d); cx <- summary(coxph(Surv(AVAL, EVENT) ~ ARM + strata(PDL1), d))
    md <- summary(sf)$table[, "median"]
    data.frame(Metric = c("Stratified Cox HR (95% CI)", "Log-rank p (one-sided)", "Median PFS, control (mo)", "Median PFS, investigational (mo)", "Events"),
               Value = c(sprintf("%.2f (%.2f\u2013%.2f)", cx$conf.int[1], cx$conf.int[3], cx$conf.int[4]), sprintf("%.4f", cx$sctest[3] / 2),
                         sprintf("%.1f", md[1]), sprintf("%.1f", md[2]), sum(d$EVENT)))
    } else if (ep == "bin") { s <- d |> group_by(ARM) |> summarise(x = sum(RESP), n = n()); pt <- prop.test(s$x[2:1], s$n[2:1], alternative = "greater")
    data.frame(Metric = c("ORR difference (95% CI)", "One-sided p", "Odds ratio"),
               Value = c(sprintf("%.1f%% (%.1f%%\u2013%.1f%%)", 100 * diff(rev(s$x / s$n)), 100 * pt$conf.int[1], 100 * min(1, diff(rev(s$x / s$n)) + (diff(rev(s$x / s$n)) - pt$conf.int[1]))),
                         sprintf("%.4f", pt$p.value), sprintf("%.2f", exp(coef(glm(RESP ~ ARM, binomial, d))[2]))))
    } else { m <- lm(CHG ~ ARM + BASE + PDL1, d); ci <- confint(m)[2, ]; p <- summary(m)$coef[2, 4]
    data.frame(Metric = c("ANCOVA adjusted difference (95% CI)", "One-sided p", "Completers"),
               Value = c(sprintf("%.2f (%.2f\u2013%.2f)", coef(m)[2], ci[1], ci[2]), sprintf("%.4f", p / 2), sum(!d$MISS)))
    }
  }, striped = TRUE, width = "100%")
  
  output$forest <- renderPlotly({
    d <- dat(); ep <- q()$ep; ov <- est(d, ep); rows <- list(data.frame(g = "Overall", n = nrow(d), e = ov[1], lo = ov[2], hi = ov[3]))
    for (v in input$sg) for (l in sort(unique(d[[v]]))) { s <- d[d[[v]] == l, ]
    x <- tryCatch(est(s, ep), error = \(e) rep(NA, 3)); rows[[length(rows) + 1]] <- data.frame(g = paste(v, l), n = nrow(s), e = x[1], lo = x[2], hi = x[3]) }
    f <- do.call(rbind, rows); f$g <- factor(f$g, rev(f$g)); nm <- c(tte = "HR", bin = "OR", cont = "Mean diff")[[ep]]
    ref <- if (ep == "cont") 0 else 1
    plot_ly(f, x = ~e, y = ~g, type = "scatter", mode = "markers", marker = list(color = col["navy"], size = 10),
            error_x = list(type = "data", symmetric = FALSE, array = f$hi - f$e, arrayminus = f$e - f$lo, color = col["navy"]),
            text = ~sprintf("%s<br>n = %d<br>%s = %.2f (%.2f\u2013%.2f)", g, n, nm, e, lo, hi), hoverinfo = "text") |>
      lay(shapes = list(list(type = "line", x0 = ref, x1 = ref, y0 = 0, y1 = 1, yref = "paper", line = list(color = "grey", dash = "dash"))),
          xaxis = list(title = paste(nm, "(95% CI)"), type = if (ep == "cont") "linear" else "log"), yaxis = list(title = ""), showlegend = FALSE)
  })
  output$miss <- renderPlotly({
    d <- dat(); s <- d |> group_by(ARM, PDL1) |> summarise(p = mean(MISS), n = n(), k = sum(MISS), .groups = "drop")
    plot_ly(s, x = ~PDL1, y = ~p, color = ~ARM, colors = c(col[["ctl"]], col[["trt"]]), type = "bar",
            text = ~sprintf("%s, PD-L1 %s<br>Missing: %d/%d (%.1f%%)", ARM, PDL1, k, n, 100 * p), hoverinfo = "text") |>
      lay(barmode = "group", xaxis = list(title = "PD-L1 subgroup"), yaxis = list(title = "Missing final assessment", tickformat = ".0%"))
  })
  output$sap <- renderTable({
    ep <- q()$ep
    data.frame(Element = c("Analysis sets", "Primary estimand", "Primary analysis", "Intercurrent events", "Multiplicity", "Interim analysis",
                           "Missing data", "Sensitivity analyses", "Software"),
               Specification = c("ITT (efficacy), safety set (all treated)",
                                 switch(ep, tte = "Treatment-policy PFS: time from randomization to progression or death", bin = "Composite strategy: patients without a confirmed response count as non-responders",
                                        cont = "Hypothetical strategy for treatment discontinuation, treatment-policy for rescue therapy"),
                                 switch(ep, tte = "Stratified log-rank and Cox model (PD-L1, ECOG, histology)", bin = "Cochran\u2013Mantel\u2013Haenszel risk difference, stratified",
                                        cont = "MMRM / ANCOVA adjusting for baseline, strata and region"),
                                 "Pre-specified strategies per ICH E9(R1); documented in the estimand table",
                                 "Fixed-sequence gatekeeping: PFS, OS, ORR at one-sided 2.5%",
                                 "Group-sequential, Lan\u2013DeMets alpha-spending, non-binding futility (optional)",
                                 "MAR via MMRM or multiple imputation; delta-adjusted tipping-point analysis",
                                 "Unstratified analysis, per-protocol set, BICR vs investigator assessment, reference-based imputation",
                                 "R (survival, gsDesign, mmrm, rpact)"))
  }, striped = TRUE, width = "100%")
  
}

shinyApp(ui, server)
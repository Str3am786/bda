# ==============================================================================
# PROJECT: THE GEOGRAPHY OF OPPORTUNITY
# A Bayesian Hierarchical Analysis of Education's Impact on Poverty
# ==============================================================================

# 1. SETUP ---------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  openintro,      # county_complete data
  brms,           # Bayesian modeling
  tidyverse,      # dplyr + ggplot2 + tibble + forcats
  bayesplot,      # diagnostics & posterior plots
  loo,            # LOO & WAIC
  tidybayes,      # tidy posterior draws
  gridExtra,
  forcats
)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
set.seed(2025)

# ==============================================================================
# 2. DATA PREPARATION ----------------------------------------------------------
# ==============================================================================

data("county_complete")

# Step 1: select only the columns we need
df_edu <- county_complete %>%
  select(
    state,
    county = name,
    poverty_2010,
    bachelors_2010,
    unemployment_rate_2010,
    density_2010
  ) %>%
  drop_na()  

cat("Rows after drop_na():", nrow(df_edu), "\n")

# Step 2: create log density safely (avoid log(0) = -Inf)
df_edu <- df_edu %>%
  mutate(
    log_density_raw = log(density_2010 + 1)  # +1 handles zero densities
  ) %>%
  filter(is.finite(log_density_raw))

cat("Rows after removing non-finite log density:", nrow(df_edu), "\n")

# Step 3: scale variables and force them to be plain numeric vectors
df_edu <- df_edu %>%
  mutate(
    poverty_scaled = as.numeric(scale(poverty_2010)),
    edu_scaled     = as.numeric(scale(bachelors_2010)),
    unemp_scaled   = as.numeric(scale(unemployment_rate_2010)),
    dens_scaled    = as.numeric(scale(log_density_raw))
  )

cat("\nSUMMARY OF SCALED VARIABLES:\n")
print(summary(df_edu$poverty_scaled))
print(summary(df_edu$edu_scaled))
print(summary(df_edu$unemp_scaled))
print(summary(df_edu$dens_scaled))

# EDA: Does the slope look different by state?
set.seed(123)
sample_states <- sample(unique(df_edu$state), 6)

p_eda <- df_edu %>%
  filter(state %in% sample_states) %>%
  ggplot(aes(x = edu_scaled, y = poverty_scaled, color = state)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title    = "Does Education reduce Poverty equally everywhere?",
    subtitle = "Visual evidence for different slopes by state",
    x        = "Bachelor Degrees (Scaled)",
    y        = "Poverty Rate (Scaled)"
  )

print(p_eda)

# ==============================================================================
# 3. MODEL DEFINITIONS ---------------------------------------------------------
# ==============================================================================

# Model 1: Pooled (education effect is the same everywhere)
f_pooled <- bf(
  poverty_scaled ~ 1 + edu_scaled + unemp_scaled + dens_scaled
)

# Model 2: Hierarchical (education slope varies by state)
f_hier <- bf(
  poverty_scaled ~ 1 + edu_scaled + unemp_scaled + dens_scaled +
    (1 + edu_scaled | state)
)

# ==============================================================================
# 4. PRIORS --------------------------------------------------------------------
# ==============================================================================

priors_hier <- c(
  prior(normal(0, 1),   class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(1), class = "sigma"),
  prior(exponential(1), class = "sd"),
  prior(lkj(2),         class = "cor")
)

# Sensitivity: "Cynical" prior saying education does almost nothing
priors_cynic <- c(
  prior(normal(0, 1),    class = "Intercept"),
  prior(normal(0, 0.05), class = "b", coef = "edu_scaled"), # VERY tight around 0
  prior(normal(0, 0.5),  class = "b"),                      # others still weak
  prior(exponential(1),  class = "sigma"),
  prior(exponential(1),  class = "sd"),
  prior(lkj(2),          class = "cor")
)

cat("\nCHECKING PRIORS FOR HIERARCHICAL MODEL:\n")
print(get_prior(f_hier, data = df_edu))

# ==============================================================================
# 5. MCMC INFERENCE ------------------------------------------------------------
# ==============================================================================

cat("\nFitting hierarchical model (main)...\n")
fit_hier <- brm(
  formula   = f_hier,
  data      = df_edu,
  prior     = priors_hier,
  family    = gaussian(),
  chains    = 4,
  iter      = 2000,
  warmup    = 1000,
  control   = list(adapt_delta = 0.95),
  seed      = 2025
)

cat("Fitting hierarchical model with cynical prior...\n")
fit_cynic <- brm(
  formula   = f_hier,
  data      = df_edu,
  prior     = priors_cynic,
  family    = gaussian(),
  chains    = 4,
  iter      = 2000,
  warmup    = 1000,
  control   = list(adapt_delta = 0.95),
  seed      = 2025
)

cat("\nSUMMARY OF MAIN HIERARCHICAL MODEL:\n")
print(summary(fit_hier))

# ==============================================================================
# 6. RESULTS & DISCOVERY (PARAMETERS & STATE-SLOPES) ---------------------------
# ==============================================================================

# 1. Compare effect of education under different priors
est_main <- fixef(fit_hier)["edu_scaled", "Estimate"]
est_cyn  <- fixef(fit_cynic)["edu_scaled", "Estimate"]

cat("\nEffect of Education (Main Model):  ", round(est_main, 4), "\n")
cat("Effect of Education (Cynic Model): ", round(est_cyn, 4), "\n")

# Also: posterior probability that the effect is negative
post_draws <- as_draws_df(fit_hier)
prob_neg   <- mean(post_draws$b_edu_scaled < 0)
cat("Posterior P(edu_scaled < 0):       ", round(prob_neg, 3), "\n")

# 2. State-specific random slopes (where does a degree "matter" most?)
state_slopes <- ranef(fit_hier)$state[, , "edu_scaled"] %>%
  as.data.frame() %>%
  rownames_to_column("state") %>%
  mutate(
    Total_Effect = Estimate + est_main
  )

ggplot(
  state_slopes %>%
    mutate(state = fct_reorder(state, Total_Effect)),
  aes(x = state, y = Total_Effect)
) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = Q2.5 + est_main,
      ymax = Q97.5 + est_main
    ),
    width = 0.2
  ) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title    = "Where does a Degree matter most?",
    subtitle = "More negative = education reduces poverty more strongly",
    y        = "Total effect of education on poverty (state-specific slope)",
    x        = "State"
  )

pp_check(fit_hier, ndraws = 50) +
  ggtitle("Posterior Predictive Check: Education–Poverty Model")

# ==============================================================================
# 7. PUBLICATION-QUALITY PLOTS -------------------------------------------------
# ==============================================================================

# ---------------------------
# 7A. Fixed effect posteriors
# ---------------------------
p_fixed <- mcmc_areas(
  fit_hier,
  pars = c("b_edu_scaled", "b_unemp_scaled", "b_dens_scaled"),
  prob = 0.95
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = c(
    "b_edu_scaled"   = "Education (scaled)",
    "b_unemp_scaled" = "Unemployment (scaled)",
    "b_dens_scaled"  = "Log density (scaled)"
  )) +
  labs(
    title = "Posterior distributions of key predictors",
    subtitle = "95% central credible intervals (dashed line = 0 effect)",
    x = "Effect size on poverty (in SD units)",
    y = NULL
  ) +
  theme_bw(base_size = 12)

print(p_fixed)
# ggsave("fig_fixed_effects.pdf", p_fixed, width = 6, height = 4)

# ----------------------------------------
# 7B. Marginal effect of education (global)
# ----------------------------------------
ce_edu <- conditional_effects(fit_hier, "edu_scaled")

p_marginal <- plot(ce_edu, plot = FALSE)[[1]] +
  labs(
    title = "Marginal effect of education on poverty",
    x = "Education (bachelor's degree rate, scaled)",
    y = "Predicted poverty (scaled)"
  ) +
  theme_bw(base_size = 12)

print(p_marginal)
# ggsave("fig_marginal_edu.pdf", p_marginal, width = 6, height = 4)

# --------------------------------------------------
# 7C. Random slopes by state (full forest-style plot)
# --------------------------------------------------
p_states <- state_slopes %>%
  mutate(state = fct_reorder(state, Total_Effect)) %>%
  ggplot(aes(x = state, y = Total_Effect)) +
  geom_point(size = 1.8) +
  geom_errorbar(
    aes(
      ymin = Q2.5 + est_main,
      ymax = Q97.5 + est_main
    ),
    width = 0.2
  ) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title    = "State-specific effect of education on poverty",
    subtitle = "More negative = stronger poverty reduction from higher education",
    x        = "State",
    y        = "Total effect of education (state-specific slope)"
  ) +
  theme_bw(base_size = 11)

print(p_states)
# ggsave("fig_state_slopes.pdf", p_states, width = 7, height = 7)

# --------------------------------------------
# 7D. Posterior predictive check (nice version)
# --------------------------------------------
p_ppc <- pp_check(fit_hier, ndraws = 50) +
  labs(
    title = "Posterior predictive check",
    subtitle = "Hierarchical education–poverty model"
  ) +
  theme_bw(base_size = 11)


print(p_ppc)
# ggsave("fig_ppc.pdf", p_ppc, width = 6, height = 4)

# ==============================================================================
# 8. MODEL COMPARISON: POOLED vs HIERARCHICAL ---------------------------------
# ==============================================================================

# --- Priors for the pooled model (no random effects) --------------------------
priors_pooled <- c(
  prior(normal(0, 1),   class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),       # same weakly informative slopes
  prior(exponential(1), class = "sigma")
)

cat("\nFitting pooled (non-hierarchical) model...\n")
fit_pooled <- brm(
  formula   = f_pooled,
  data      = df_edu,
  prior     = priors_pooled,
  family    = gaussian(),
  chains    = 4,
  iter      = 2000,
  warmup    = 1000,
  control   = list(adapt_delta = 0.95),
  seed      = 2025
)

# --- LOO cross-validation -----------------------------------------------------
cat("\nComputing LOO for pooled and hierarchical models...\n")
loo_pooled <- loo(fit_pooled)
loo_hier   <- loo(fit_hier)

cat("\nLOO results: pooled model\n")
print(loo_pooled)

cat("\nLOO results: hierarchical model\n")
print(loo_hier)

cat("\nLOO comparison (higher elpd_loo = better predictive performance):\n")
loo_comp <- loo_compare(loo_pooled, loo_hier)
print(loo_comp)


# --- WAIC (optional, similar idea to LOO) -------------------------------------
cat("\nComputing WAIC for pooled and hierarchical models...\n")
waic_pooled <- waic(fit_pooled)
waic_hier   <- waic(fit_hier)

cat("\nWAIC results: pooled model\n")
print(waic_pooled)

cat("\nWAIC results: hierarchical model\n")
print(waic_hier)

# Robust WAIC comparison table (avoids class issues with loo_compare)
waic_table <- tibble(
  model = c("Pooled", "Hierarchical"),
  waic  = c(
    waic_pooled$estimates["waic", "Estimate"],
    waic_hier$estimates["waic", "Estimate"]
  ),
  se_waic = c(
    waic_pooled$estimates["waic", "SE"],
    waic_hier$estimates["waic", "SE"]
  )
) %>%
  mutate(delta_waic = waic - min(waic))

cat("\nWAIC comparison table (lower WAIC = better):\n")
print(waic_table)

# ==============================================================================
# 9. HIGHEST-RETURN AND LOWEST-RETURN STATES ----------------------------------
# ==============================================================================

# More negative Total_Effect = stronger poverty reduction from education

# Top 10 states (strongest benefit from education)
top_states <- state_slopes %>%
  arrange(Total_Effect) %>%            # ascending (most negative first)
  head(10) %>%
  select(
    state,
    Mean_Random_Effect = Estimate,
    CI_low             = Q2.5,
    CI_high            = Q97.5,
    Total_Effect
  )

cat("\nTop 10 states: strongest reduction in poverty from higher education\n")
print(top_states)

# Bottom 10 states (weakest benefit / possibly near zero)
bottom_states <- state_slopes %>%
  arrange(desc(Total_Effect)) %>%      # descending (least negative / most positive)
  head(10) %>%
  select(
    state,
    Mean_Random_Effect = Estimate,
    CI_low             = Q2.5,
    CI_high            = Q97.5,
    Total_Effect
  )

cat("\nBottom 10 states: weakest reduction (or possibly no reduction) from education\n")
print(bottom_states)



# ==============================================================================
# 10. RESULTS SECTION (KEY NUMBERS COLLECTED) ----------------------------------
# ==============================================================================

cat("\n\n==================== SUMMARY RESULTS ====================\n")
cat("Global education effect (main prior):   ", round(est_main, 4), "\n")
cat("Global education effect (cynic prior):  ", round(est_cyn, 4), "\n")
cat("Posterior P(edu_scaled < 0):            ", round(prob_neg, 3), "\n\n")

cat("LOO comparison (positive elpd_diff -> better model):\n")
print(loo_comp)

cat("\nWAIC comparison (lower WAIC -> better model):\n")
print(waic_table)

cat("\nTop 5 'highest-return' states (education reduces poverty most):\n")
print(top_states %>% head(5))

cat("\nTop 5 'lowest-return' states (weakest education effect):\n")
print(bottom_states %>% head(5))

cat("==========================================================\n")

# Store everything in a single object for interactive use
results_list <- list(
  est_main        = est_main,
  est_cyn         = est_cyn,
  prob_neg        = prob_neg,
  loo_pooled      = loo_pooled,
  loo_hier        = loo_hier,
  loo_comparison  = loo_comp,
  waic_pooled     = waic_pooled,
  waic_hier       = waic_hier,
  waic_table      = waic_table,
  top_states      = top_states,
  bottom_states   = bottom_states,
  fit_hier        = fit_hier,
  fit_cynic       = fit_cynic,
  fit_pooled      = fit_pooled,
  state_slopes    = state_slopes
)


# Current Alabama data, but with HIGHER education
hypothetical_alabama <- df_edu %>%
  filter(state == "Alabama") %>%
  mutate(edu_scaled = edu_scaled + 1) # Add 1 SD of education

preds <- posterior_predict(fit_hier, newdata = hypothetical_alabama)
diff <- preds - df_edu$poverty_scaled[df_edu$state == "Alabama"]
mean(diff) # The expected reduction in poverty



# PROOFS NEEDED:
# Trace Plots ("The Caterpillars")
mcmc_trace(fit_hier, pars = c("b_edu_scaled", "b_Intercept", "sigma"))

# R hat verification:
rhat_values <- brms::rhat(fit_hier)
print(any(rhat_values > 1.05))
neff_ratio(fit_hier) # Gives Ratio (ESS / Total Draws)
# OR
print(summary(fit_hier)) # Look at 'Bulk_ESS' and 'Tail_ESS' columns

# Checking directly in the object
np <- nuts_params(fit_hier)
sum(subset(np, Parameter == "divergent__")$Value)

np <- nuts_params(fit_hier)

# Count how many times the skater flew off the bowl
div_count <- sum(subset(np, Parameter == "divergent__")$Value)
cat("Number of Divergent Transitions:", div_count, "\n")


rhat_values <- brms::rhat(fit_hier)
print(any(rhat_values > 1.01))

rhat_values <- brms::rhat(fit_pooled)
print(any(rhat_values > 1.01))

rhat_values <- brms::rhat(fit_cynic)
print(any(rhat_values > 1.01))

# --------------------------------------------
# Posterior predictive check: POOLED model
# --------------------------------------------

# Quick check (base bayesplot style)
pp_check(fit_pooled, ndraws = 50) +
  ggtitle("Posterior Predictive Check: Pooled Education–Poverty Model")

# Publication-style version (analogous to p_ppc for fit_hier)
p_ppc_pooled <- pp_check(fit_pooled, ndraws = 50) +
  labs(
    title    = "Posterior predictive check",
    subtitle = "Pooled education–poverty model"
  ) +
  theme_bw(base_size = 11)

print(p_ppc_pooled)
# ggsave("fig_ppc_pooled.pdf", p_ppc_pooled, width = 6, height = 4)

pp_check(fit_hier, ndraws = 50) +
  ggtitle("Posterior Predictive Check: Education–Poverty Model hierarchical")


mcmc_trace(
  fit_pooled,
  pars = c("b_edu_scaled", "b_unemp_scaled", "b_dens_scaled", "b_Intercept", "sigma")
)

# ---------- ESS SUMMARY ----------
mcmc_trace(fit_pooled, pars = c("b_edu_scaled", "b_Intercept", "sigma"))

ce_edu <- conditional_effects(fit_hier, "edu_scaled")

p_marginal <- plot(ce_edu, plot = FALSE, points = TRUE)[[1]] +
  labs(
    title = "Marginal effect of education on poverty",
    x = "Education (bachelor's degree rate, scaled)",
    y = "Predicted poverty (scaled)"
  ) +
  theme_bw(base_size = 12)
print(p_marginal)

# Load data (if not already)
library(openintro)
data("county_complete")

# Core helpers
library(dplyr)
library(ggplot2)

#----------------------------------------------------------
# 1. Basic structure and overview
#----------------------------------------------------------

# How big is the dataset?
dim(county_complete)       # rows, columns
names(county_complete)     # column names
str(county_complete)       # structure (types, first values)

# Quick summary stats for all variables
summary(county_complete)

#----------------------------------------------------------
# 2. Number of missing values (NULLs / NAs)
#----------------------------------------------------------

# Total number of NAs in the whole dataset
sum(is.na(county_complete))

# NAs per column
na_per_col <- sapply(county_complete, function(x) sum(is.na(x)))
na_per_col

# NAs per column in a tidy tibble
na_table <- data.frame(
  variable   = names(na_per_col),
  na_count   = as.integer(na_per_col),
  na_percent = round(100 * na_per_col / nrow(county_complete), 2)
) %>%
  arrange(desc(na_count))

na_table

# Example: poverty_2010
ggplot(county_complete, aes(x = poverty_2010)) +
  geom_histogram(bins = 30, na.rm = TRUE) +
  labs(
    title = "Distribution of poverty_2010",
    x     = "Poverty rate (%)",
    y     = "Count of counties"
  )

# Example: bachelors_2010
ggplot(county_complete, aes(x = bachelors_2010)) +
  geom_histogram(bins = 30, na.rm = TRUE) +
  labs(
    title = "Distribution of bachelors_2010",
    x     = "Bachelor's degree rate (%)",
    y     = "Count of counties"
  )

# Example: unemployment_rate_2010
ggplot(county_complete, aes(x = unemployment_rate_2010)) +
  geom_histogram(bins = 30, na.rm = TRUE) +
  labs(
    title = "Distribution of unemployment_rate_2010",
    x     = "Unemployment rate (%)",
    y     = "Count of counties"
  )

# Example: density_2010 (can be heavily skewed, so log scale is nice)
ggplot(county_complete, aes(x = density_2010)) +
  geom_histogram(bins = 30, na.rm = TRUE) +
  scale_x_log10() +
  labs(
    title = "Distribution of density_2010 (log scale)",
    x     = "Population density (log10 scale)",
    y     = "Count of counties"
  )


num_vars <- county_complete %>%
  select(poverty_2010, bachelors_2010, unemployment_rate_2010, density_2010) %>%
  na.omit()

cor(num_vars)



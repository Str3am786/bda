# ==============================================================================
# PROJECT: DRIVING UNDER PRESSURE
# A Bayesian Hierarchical Analysis of Economic Indicators on US Traffic Fatalities
# ==============================================================================

# 1. SETUP & LIBRARIES ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  AER,          # Dataset
  brms,         # Bayesian Modeling
  tidyverse,    # Data Wrangling
  bayesplot,    # Plotting MCMC
  loo,          # Model Comparison
  gridExtra,    # Plot arrangement
  tidybayes     # Handling posterior draws
)

# Optimization for speed (uses CmdStan if available, otherwise falls back)
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
set.seed(2025)

# ==============================================================================
# 2. DATA PREPARATION & FEATURE ENGINEERING ------------------------------------
# Requirement: "Description of the data... Motivation"
# ==============================================================================

data("Fatalities", package = "AER")

# We are pivoting to ECONOMICS. 
# Theory: Higher unemployment -> Less commuting/vacations -> Fewer deaths?
# Theory: Higher Income -> Safer/Newer cars -> Fewer deaths?

df_econ <- Fatalities %>%
  mutate(
    # 1. OUTCOME: Normalized Fatality Rate (per 10,000 people)
    fatality_rate = (fatal / pop) * 10000,
    
    # 2. MAIN PREDICTOR: Log Income (Wealth effects are multiplicative)
    # Scaling it makes priors easier to set (Normal(0,1) works better on scaled data)
    log_income_scaled = scale(log(income)),
    
    # 3. MAIN PREDICTOR: Unemployment
    # Centering allows us to interpret the intercept as "Average Unemployment"
    unemp_centered = scale(unemp, scale = FALSE),
    
    # 4. CONTROLS: 
    # Alcohol laws are now just controls, not the main story.
    jail_binary = factor(ifelse(jail == "yes", 1, 0)),
    
    # We must control for 'miles' driven per person to isolate the *economic* effect
    # from just "people driving less". 
    miles_per_driver = miles
  ) %>%
  select(state, year, fatality_rate, unemp_centered, log_income_scaled, 
         jail_binary, miles_per_driver)

# ==============================================================================
# 3. EXPLORATORY DATA ANALYSIS (EDA) -------------------------------------------
# Requirement: "Illustrative figure"
# ==============================================================================

# Plot 1: The Raw Economic Relationship (Pooled)
p1 <- ggplot(df_econ, aes(x = unemp_centered, y = fatality_rate)) +
  geom_point(alpha = 0.3, color = "darkblue") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Unemployment vs. Fatalities (Pooled)",
       x = "Unemployment Rate (Centered)", y = "Fatality Rate")

# Plot 2: Heterogeneity by State (Justifying Hierarchical Model)
# Does the relationship look different in NY vs WY?
p2 <- df_econ %>%
  filter(state %in% c("ny", "ca", "wy", "tx", "fl", "nd")) %>% # Sample states
  ggplot(aes(x = unemp_centered, y = fatality_rate, color = state)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Economic Impact Varies by State",
       subtitle = "Slopes represent the effect of Unemployment",
       x = "Unemployment", y = "Fatality Rate")

grid.arrange(p1, p2, ncol = 2)

# ==============================================================================
# 4. MODEL DEFINITIONS ---------------------------------------------------------
# Requirement: "At least two models... Hierarchical vs Non-Hierarchical"
# ==============================================================================

# --- MODEL 1: THE "NAIVE" POOLED MODEL ---
# Assumes the effect of Unemployment is identical across the entire USA.
f_pooled <- bf(fatality_rate ~ 1 + unemp_centered + log_income_scaled + 
                 miles_per_driver + jail_binary)

# --- MODEL 2: HIERARCHICAL RANDOM SLOPES ---
# The "Fresh" Angle: We allow the effect of Unemployment to VARY by state.
# (1 + unemp_centered | state) -> Random Intercept + Random Slope for Unemployment
# This tests: "Does economic stress impact safety differently in different states?"
f_hier <- bf(fatality_rate ~ 1 + unemp_centered + log_income_scaled + 
               miles_per_driver + jail_binary + 
               (1 + unemp_centered | state))

# ==============================================================================
# 5. PRIORS --------------------------------------------------------------------
# Requirement: "Informative or weakly informative priors... Justification"
# ==============================================================================

# Justification:
# Intercept: The mean fatality rate is around 2.0. Prior: Normal(2, 0.5).
# Slopes (unemp, income): We expect economic effects to be subtle, not massive.
# A 1% change in unemployment won't double deaths. Normal(0, 0.5) is conservative.
# Sigma: Must be positive. Exponential(1).

priors <- c(
  prior(normal(2, 0.5), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),        # Weakly informative for predictors
  prior(exponential(1), class = "sigma"),    # Residual noise
  prior(exponential(1), class = "sd"),       # Group-level variability
  prior(lkj(2), class = "cor")               # Correlation between Int and Slope
)

# ==============================================================================
# 6. MCMC INFERENCE ------------------------------------------------------------
# Requirement: "How MCMC was run... Convergence diagnostics"
# ==============================================================================

# Fit Pooled
fit_pooled <- brm(
  formula = f_pooled, data = df_econ, prior = priors[1:3],
  family = gaussian(),
  chains = 4, iter = 2000, warmup = 1000, seed = 2025,
  save_pars = save_pars(all = TRUE)
)

# Fit Hierarchical (Random Slopes)
# 
fit_hier <- brm(
  formula = f_hier, data = df_econ, prior = priors,
  family = gaussian(),
  chains = 4, iter = 2000, warmup = 1000, seed = 2025,
  control = list(adapt_delta = 0.95), # Increased for complex random slopes
  save_pars = save_pars(all = TRUE)
)

# Diagnostics
print(summary(fit_hier))
# Look for Rhat < 1.01 and high ESS.
mcmc_trace(fit_hier, pars = c("b_unemp_centered", "b_log_income_scaled"))

# ==============================================================================
# 7. POSTERIOR PREDICTIVE CHECKS (PPC) -----------------------------------------
# Requirement: "PPC... Interpretation"
# ==============================================================================

# Visualizing if the model captures the data distribution
pp_check(fit_hier, ndraws = 50) + 
  ggtitle("Posterior Predictive Check: Hierarchical Model")

# ==============================================================================
# 8. MODEL COMPARISON (LOO-CV) -------------------------------------------------
# Requirement: "Model comparison with LOO-CV"
# ==============================================================================

loo_pooled <- loo(fit_pooled)
loo_hier <- loo(fit_hier)

# A negative 'elpd_diff' for fit_pooled means fit_hier is better.
print(loo_compare(loo_pooled, loo_hier))

# ==============================================================================
# 9. SENSITIVITY ANALYSIS ------------------------------------------------------
# Requirement: "Sensitivity analysis with respect to prior choices"
# ==============================================================================

# HYPOTHESIS TEST: The "Economic Skeptic"
# What if we assume Unemployment has ALMOST ZERO effect? 
# We use a very tight prior: Normal(0, 0.05). 
# If the posterior still moves away from zero, the data signal is strong.

prior_skeptic <- c(
  prior(normal(2, 0.5), class = "Intercept"),
  prior(normal(0, 0.05), class = "b", coef = "unemp_centered"), # TIGHT PRIOR
  prior(exponential(1), class = "sigma"),
  prior(exponential(1), class = "sd"),
  prior(lkj(2), class = "cor") 
)

fit_sensitivity <- update(fit_hier, prior = prior_skeptic)

# Compare Estimates
est_orig <- fixef(fit_hier)["unemp_centered", "Estimate"]
est_sens <- fixef(fit_sensitivity)["unemp_centered", "Estimate"]

cat(paste("Unemployment Effect (Original):", round(est_orig, 4), "\n"))
cat(paste("Unemployment Effect (Skeptic): ", round(est_sens, 4), "\n"))

# ==============================================================================
# 10. RESULTS & CONCLUSION -----------------------------------------------------
# Requirement: "Conclusion... what was learned"
# ==============================================================================

# 1. Did Unemployment decrease fatalities?
# Plot the main fixed effect
mcmc_areas(fit_hier, pars = c("b_unemp_centered", "b_log_income_scaled"), prob = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Effect of Economy on Fatalities (95% CI)")

# 2. The Hierarchical "Reveal" (Random Slopes)
# Which states are MOST sensitive to economic changes?
# We extract the random slopes for unemployment.
state_effects <- ranef(fit_hier)$state[, , "unemp_centered"] %>%
  as.data.frame() %>%
  rownames_to_column("state") %>%
  arrange(Estimate)

# Plotting the variation of the Unemployment Effect by State
ggplot(state_effects, aes(x = reorder(state, Estimate), y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Heterogeneity of Economic Impact",
       subtitle = "Does Unemployment reduce fatalities more in some states than others?",
       y = "Effect of Unemployment on Fatality Rate", x = "State")

# INTERPRETATION GUIDE FOR REPORT:
# 1. If 'b_unemp_centered' is NEGATIVE: Economic stress (unemployment) SAVES lives.
#    (Likely mechanism: Less discretionary driving, fewer road trips).
# 2. If 'b_log_income_scaled' is NEGATIVE: Wealthier states are safer.
#    (Likely mechanism: Better infrastructure, newer cars with airbags/ABS).
# 3. If State Effects vary (the last plot): The economy impacts rural vs urban states differently.
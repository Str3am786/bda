# ==============================================================================
# PROJECT: STEEL BARS VS. STEEL BARRELS
# A Bayesian Hierarchical Comparison of Incarceration vs. Gun Laws on Crime
# ==============================================================================

# 1. SETUP & LIBRARIES ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  AER,          # Dataset ('Guns')
  brms,         # Bayesian Modeling
  tidyverse,    # Data Wrangling
  bayesplot,    # Plotting MCMC
  loo,          # Model Comparison
  gridExtra,    # Plot arrangement
  tidybayes     # Handling posterior draws
)

# Optimization for speed
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
set.seed(2025)

# ==============================================================================
# 2. DATA PREPARATION & FEATURE ENGINEERING ------------------------------------
# Requirement: "Description of the data... Motivation"
# ==============================================================================

data("Guns", package = "AER")

# PIVOT: We are comparing "Incarceration" (prisoners) vs "Gun Law" (law).
# We also want to see if Density changes how well Prison works.

df_crime <- Guns %>%
  mutate(
    # 1. OUTCOME: Violent Crime Rate
    # Crime rates are strictly positive and often skewed. Log-transforming 
    # stabilizes variance and implies multiplicative effects (standard in criminology).
    log_violent = log(violent),
    
    # 2. MAIN PREDICTOR A: Incarceration Rate (Lagged by 1 year in dataset)
    # We log this because the 1000th prisoner likely matters less than the 10th.
    # We scale it to make priors easier (mean=0, sd=1).
    log_prisoners_scaled = scale(log(prisoners)),
    
    # 3. MAIN PREDICTOR B: Shall-Carry Law (The "Gun" variable)
    # We keep this to compare its effect size against prisoners.
    law_binary = factor(law, levels = c("no", "yes")),
    
    # 4. INTERACTION VARIABLE: Population Density
    # Does prison work differently in crowded cities vs. rural areas?
    log_density_scaled = scale(log(density)),
    
    # 5. CONTROLS:
    # Income and Demographics (Young Males are the most crime-prone group)
    log_income_scaled = scale(log(income)),
    male_scaled = scale(male)
  ) %>%
  select(state, year, log_violent, log_prisoners_scaled, law_binary, 
         log_density_scaled, log_income_scaled, male_scaled)

# ==============================================================================
# 3. EXPLORATORY DATA ANALYSIS (EDA) -------------------------------------------
# Requirement: "Illustrative figure"
# ==============================================================================

# Plot 1: The "Simple" View (Prisoners vs Crime)
p1 <- ggplot(df_crime, aes(x = log_prisoners_scaled, y = log_violent)) +
  geom_point(alpha = 0.3, color = "darkred") +
  geom_smooth(method = "lm", color = "black") +
  labs(title = "Incarceration vs. Violence (Pooled)",
       x = "Log Prisoners (Scaled)", y = "Log Violent Crime Rate")

# Plot 2: The Hierarchical Hint (Density Interaction)
# We split states into "High Density" and "Low Density" just for the plot
p2 <- df_crime %>%
  mutate(density_group = ifelse(log_density_scaled > 0, "High Density", "Low Density")) %>%
  ggplot(aes(x = log_prisoners_scaled, y = log_violent, color = density_group)) +
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.1) +
  labs(title = "Does Density Change Effectiveness?",
       subtitle = "Slope comparison: High vs Low Density States",
       x = "Incarceration", y = "Violent Crime")

grid.arrange(p1, p2, ncol = 2)

# ==============================================================================
# 4. MODEL DEFINITIONS ---------------------------------------------------------
# Requirement: "At least two models... Hierarchical vs Non-Hierarchical"
# ==============================================================================

# --- MODEL 1: POOLED INTERACTION MODEL ---
# Assumes the relationship is the same everywhere.
# We include 'log_prisoners_scaled * log_density_scaled' to test if 
# density modifies the effect of prison.
f_pooled <- bf(log_violent ~ 1 + log_prisoners_scaled * log_density_scaled + 
                 law_binary + log_income_scaled + male_scaled)

# --- MODEL 2: HIERARCHICAL RANDOM SLOPES ---
# The "Fresh" Angle: 
# 1. We keep the density interaction (Fixed Effect).
# 2. We allow the "Effect of Incarceration" to VARY by state (Random Slope).
#    Maybe Prison reduces crime in Texas but increases it in California?
# (1 + log_prisoners_scaled | state) -> Random Intercept + Random Slope
f_hier <- bf(log_violent ~ 1 + log_prisoners_scaled * log_density_scaled + 
               law_binary + log_income_scaled + male_scaled + 
               (1 + log_prisoners_scaled | state))

# ==============================================================================
# 5. PRIORS --------------------------------------------------------------------
# Requirement: "Informative or weakly informative priors... Justification"
# ==============================================================================

# Justification:
# Intercept: Log crime rates are usually around 5 to 7. Prior: Normal(6, 1).
# Slopes: We assume incarceration REDUCES crime (negative slope), but we don't 
# want to force it. Normal(0, 0.5) allows it to go positive if data says so.
# Gun Law: The debate is fierce, so we use a Skeptical prior Normal(0, 0.2) 
# to assume the law has a small effect unless data proves otherwise.

priors <- c(
  prior(normal(6, 1), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),        # Generic predictors
  prior(normal(0, 0.2), class = "b", coef = "law_binaryyes"), # Gun Law specific
  prior(exponential(1), class = "sigma"),    # Residual noise
  prior(exponential(1), class = "sd"),       # State variability
  prior(lkj(2), class = "cor")               # Correlation
)

# ==============================================================================
# 6. MCMC INFERENCE ------------------------------------------------------------
# Requirement: "How MCMC was run... Convergence diagnostics"
# ==============================================================================

# Fit Pooled
fit_pooled <- brm(
  formula = f_pooled, data = df_crime, prior = priors[1:3],
  family = gaussian(),
  chains = 4, iter = 2000, warmup = 1000, seed = 2025,
  save_pars = save_pars(all = TRUE)
)

# Fit Hierarchical (Random Slopes)
# [Image of random slope visualization]
fit_hier <- brm(
  formula = f_hier, data = df_crime, prior = priors,
  family = gaussian(),
  chains = 4, iter = 2000, warmup = 1000, seed = 2025,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE)
)

# Diagnostics
print(summary(fit_hier))
# Check Rhat < 1.01 and good ESS

# ==============================================================================
# 7. POSTERIOR PREDICTIVE CHECKS (PPC) -----------------------------------------
# Requirement: "PPC... Interpretation"
# ==============================================================================

pp_check(fit_hier, ndraws = 50) + 
  ggtitle("PPC: Does the model simulate realistic crime rates?")

# ==============================================================================
# 8. MODEL COMPARISON (LOO-CV) -------------------------------------------------
# Requirement: "Model comparison with LOO-CV"
# ==============================================================================

loo_pooled <- loo(fit_pooled)
loo_hier <- loo(fit_hier)

# Comparison
print(loo_compare(loo_pooled, loo_hier))

# ==============================================================================
# 9. SENSITIVITY ANALYSIS ------------------------------------------------------
# Requirement: "Sensitivity analysis with respect to prior choices"
# ==============================================================================

# HYPOTHESIS TEST: The "Abolitionist" Prior
# What if we assume Prisons DO NOT reduce crime? 
# We set a tight prior centered at 0 for 'log_prisoners_scaled'.
prior_skeptic <- c(
  prior(normal(6, 1), class = "Intercept"),
  prior(normal(0, 0.05), class = "b", coef = "log_prisoners_scaled"), # TIGHT
  prior(exponential(1), class = "sigma"),
  prior(exponential(1), class = "sd"),
  prior(lkj(2), class = "cor")
)

fit_sensitivity <- update(fit_hier, prior = prior_skeptic)

# Compare Results
est_main <- fixef(fit_hier)["log_prisoners_scaled", "Estimate"]
est_sens <- fixef(fit_sensitivity)["log_prisoners_scaled", "Estimate"]

cat(paste("Prison Effect (Original):", round(est_main, 4), "\n"))
cat(paste("Prison Effect (Skeptic): ", round(est_sens, 4), "\n"))

# ==============================================================================
# 10. RESULTS & CONCLUSION -----------------------------------------------------
# Requirement: "Conclusion... what was learned"
# ==============================================================================

# 1. Compare Effect Sizes (Prison vs Guns)
# We plot the posterior distributions of the coefficients side-by-side
# Note: 'law_binaryyes' vs 'log_prisoners_scaled'
mcmc_areas(fit_hier, 
           pars = c("b_log_prisoners_scaled", "b_law_binaryyes", "b_log_prisoners_scaled:log_density_scaled"),
           prob = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Relative Impact: Incarceration vs. Gun Laws")

# 2. Visualizing the Interaction (Density Effect)
# If the interaction term (prisoners:density) is significant, we visualize it.
# conditional_effects plots the "marginal effects"
conditional_effects(fit_hier, effects = "log_prisoners_scaled:log_density_scaled")

# 3. State Heterogeneity
# Which states have the strongest "Prison Effect"?
state_slopes <- ranef(fit_hier)$state[, , "log_prisoners_scaled"] %>%
  as.data.frame() %>%
  rownames_to_column("state") %>%
  arrange(Estimate)

ggplot(state_slopes, aes(x = reorder(state, Estimate), y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5)) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Variation in Incarceration Effectiveness",
       subtitle = "Negative values = More prisoners lead to less crime",
       y = "Effect of Prisoners on Crime", x = "State")

# INTERPRETATION FOR REPORT:
# 1. 'b_log_prisoners_scaled': If negative, prisons reduce crime overall.
# 2. 'b_law_binaryyes': Compare the width/location of this posterior to the prison one. 
#    Likely, the prison effect is stronger (further from 0).
# 3. 'b_log_prisoners_scaled:log_density_scaled': 
#    - If POSITIVE, it means the negative effect of prison gets WEAKER (closer to 0) 
#      as density increases. (i.e., Prisons work better in rural areas).
#    - If NEGATIVE, it means prisons work BETTER in dense cities.
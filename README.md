# The Geography of Opportunity: A Bayesian Hierarchical Analysis

**Authors:** Andrei Iliescu & Miguel Arroyo Marquez  
**Date:** December 2025

![Project Status](https://img.shields.io/badge/Status-Completed-success)
![R Version](https://img.shields.io/badge/R-4.0%2B-blue)
![Stan](https://img.shields.io/badge/Stan-brms-orange)

## Project Overview

Does a college degree offer the same economic protection against poverty in the industrial Midwest as it does in the coastal South? 

This project explores the relationship between **educational attainment (bachelor's degrees)** and **poverty rates** across 3,139 U.S. counties using data from the 2010 Decennial Census. 

Using **Bayesian Hierarchical Modeling**, we move beyond a "one-size-fits-all" national regression to determine if the "return on education" varies by state geography.

### Key Research Questions
1. Does the impact of education on poverty reduction remain constant across the country?
2. Is a hierarchical model statistically superior to a pooled model for this data?
3. How robust are these findings against a "cynical" (skeptical) prior?

---

## 📊 Data

The analysis uses the `county_complete` dataset from the **openintro** R package.
* **Source:** United States Census Bureau (2010 Decennial Census).
* **Observation Level:** U.S. Counties ($N=3,139$) nested within 51 States/regions.
* **Outcome Variable:** `poverty_2010` (Percentage of population below federal poverty line).
* **Key Predictors:** * `bachelors_2010` (Education)
    * `unemployment_rate_2010` (Control)
    * `density_2010` (Control, log-transformed)

---

## Methodology

We utilized **Bayesian Inference** via the `brms` interface for Stan.

### Models Compared
1.  **Pooled Model (Baseline):** Assumes a single, fixed slope for education across the entire U.S.
2.  **Hierarchical Model:** Allows for **state-specific intercepts** and **state-specific education slopes**. This enables partial pooling, where states borrow statistical strength from one another.

### Sensitivity Analysis: The "Cynical Prior"
To test the robustness of our findings, we fitted a model with a "Cynical Prior"—a highly skeptical prior distribution ($\mathcal{N}(0, 0.05)$) placed on the education coefficient. This tests the hypothesis: *"If I assume education does nothing, is the signal in the data strong enough to convince me otherwise?"*

### Model Checking
* **Convergence:** Assessed via R-hat ($\hat{R} \leq 1.01$), Trace Plots, and Effective Sample Size (ESS).
* **Comparison:** Leave-One-Out Cross-Validation (LOO-CV) and WAIC.
* **Posterior Predictive Checks:** Validated that simulated data matched observed distributions.

---

## Key Findings

1.  **Geography Matters:** The hierarchical model decisively outperformed the pooled model (LOO-CV difference of ~428 ELPD). The relationship between education and poverty is **not** uniform across the U.S.
2.  **Robust Negative Association:** Even under the "Cynical Prior," the posterior distribution for education's effect remained negative and significant.
3.  **State Variation:**
    * **High Return:** States like *Kentucky, South Dakota, and Georgia* show the strongest reduction in poverty for every increase in bachelor's degrees.
    * **Low Return:** States like *Indiana and Wyoming* show much flatter slopes, suggesting education is less of a differentiating factor for poverty in those economies.
4.  **Counterfactual Simulation:** A simulation of **Alabama** indicated that increasing the county-level bachelor's rate by 1 standard deviation would result in a significant decrease in poverty, holding unemployment and density constant.

*(See the `reports/` folder for the full PDF analysis including detailed plots).*

---

## Installation & Usage

### Prerequisites
* R (version 4.0 or higher)
* A working C++ compiler (for Stan/brms)

### Required Libraries
Run the following in R to install dependencies:
```r
install.packages("pacman")
pacman::p_load(openintro, brms, tidyverse, bayesplot, loo, tidybayes, forcats, gridExtra)
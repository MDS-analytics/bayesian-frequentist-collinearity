#    /\_/\  
#   ( o.o ) 
#    > ^ <    "Prior beliefs about treats: Very high!"


# --- Purpose ---------------------------------------------------------------
# Simulate linear data under varying sample sizes (n) and collinearity levels (rho),
# fit both frequentist (lm) and Bayesian (brms) regressions before and after
# residualizing x2 on x1 (a simple FWL-style deconfounding step).
#
# Naming convention for saved model objects (created via `assign()`):
#   - Pre-residualization:
#       lm_n{n}_rho{rho}        : OLS model with raw x1, x2, x3
#       brm_n{n}_rho{rho}       : Bayesian model with raw x1, x2, x3
#   - Post-residualization (x2_resid replaces x2):
#       lm_resid_n{n}_rho{rho}  : OLS model with x2 residualized on x1
#       brm_resid_n{n}_rho{rho} : Bayesian model with x2 residualized on x1
#
# Example object names: "lm_n200_rho0.2", "brm_resid_n4000_rho0.9"
# This makes it easy to grep or programmatically select models by n and rho.

set.seed(42)  # ensure reproducibility across all simulation settings

# Grid of experimental conditions:
# - n controls sampling variability (200, 1000, 4000)
# - rho controls the induced correlation between x1 and x2 (0.2=low, 0.9=high)
n_vals   <- c(200, 1000, 4000)
rho_vals <- c(0.2, 0.5, 0.9)

for (n in n_vals) {
  for (rho in rho_vals) {
    
    # --- 1) Simulate data ---------------------------------------------------
    # x1 and x3 are standard normals; x2 is constructed to be correlated with x1.
    # Increasing rho increases collinearity between x1 and x2.
    x1 <- rnorm(n)
    x2 <- rho * x1 + rnorm(n, sd=0.1)   # strongly tied to x1 when rho is large
    x3 <- rnorm(n)
    
    # Linear DGP with moderate SNR; opposing signs on x1 and x2 highlight collinearity effects
    y  <- 1 + 1.0*x1 - 1.0*x2 + 0.5*x3 + rnorm(n, sd=1)
    
    dat <- data.frame(y, x1, x2, x3)
    
    # --- 2) Fit pre-residualization models ---------------------------------
    # Save models with names that encode n and rho for easy retrieval later.
    lm_name  <- sprintf("lm_n%d_rho%.1f",  n, rho)
    brm_name <- sprintf("brm_n%d_rho%.1f", n, rho)
    
    # OLS baseline
    assign(lm_name, lm(y ~ x1 + x2 + x3, data=dat))
    
    # Bayesian Gaussian regression (weak default priors).
    # Note: iter/warmup kept modest for speed; increase for real inference.
    assign(
      brm_name,
      brm(
        y ~ x1 + x2 + x3,
        data  =dat,
        family=gaussian(),
        chains=4, iter=6000, warmup=1000,
        refresh=0, seed=42
      )
    )
    
    # --- 3) Residualize x2 on x1 (FWL-style step) ---------------------------
    # This constructs x2_resid=x2 - E[x2 | x1]; i.e., the portion of x2
    # orthogonal to x1. Including x2_resid instead of x2 reduces collinearity
    # between regressors and can stabilize coefficient estimates.
    x2_resid  <- resid(lm(x2 ~ x1, data=dat))
    dat_resid <- transform(dat, x2_resid=x2_resid)
    
    # --- 4) Fit post-residualization models --------------------------------
    lm_resid_name  <- sprintf("lm_resid_n%d_rho%.1f",  n, rho)
    brm_resid_name <- sprintf("brm_resid_n%d_rho%.1f", n, rho)
    
    assign(lm_resid_name, lm(y ~ x1 + x2_resid + x3, data=dat_resid))
    
    assign(
      brm_resid_name,
      brm(
        y ~ x1 + x2_resid + x3,
        data  =dat_resid,
        family=gaussian(),
        chains=4, iter=6000, warmup=1000,
        refresh=0, seed=42
      )
    )
  }
}

# --- Notes ------------------------------------------------------------------
# • The residualization illustrates the Frisch–Waugh–Lovell idea for isolating
#   the unique contribution of a correlated regressor.
# • Because object names encode (n, rho) and whether residualization was applied,
#   you can iterate or compare easily, e.g.:
#     get(sprintf("lm_n%d_rho%.1f", 1000, 0.5))
#     get(sprintf("brm_resid_n%d_rho%.1f", 4000, 0.9))
# • Consider storing models in lists (nested by n and rho) instead of `assign()`
#   if you prefer stricter scope and easier programmatic handling.

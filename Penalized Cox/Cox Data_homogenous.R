rm(list = ls())

library(mvtnorm)
library(survival)
library(survminer)
library(ggplot2)
library(mvtnorm)

sim_cox_no_clusters <- function(N = 100, p = 50, q = 10, rho = 0,
                                   lambda0 = 0.5, censor_rate = 0.01) {
  

  grp <- factor(1:N)

  # AR(1) covariance matrix
  covmat <- matrix(0, p, p)
  for (i in 1:p) for (j in 1:p) covmat[i, j] <- rho^abs(i - j)
  
  # Covariates (with intercept)
  # x <- cbind(1, rmvnorm(n_subjects, mean = rep(1, p), sigma = covmat))
  # colnames(x) <- c("Intercept", paste0("x.", 1:p))
  
  x <- cbind(1, rmvnorm(N, mean = rep(1, p), sigma = covmat))
  colnames(x) <- paste0("x.", 1:(p + 1))  # intercept becomes x.1
  
  # Random effects design matrix
  z <- x[, 1:q, drop = FALSE]
  colnames(z) <- paste0("z.", 1:q)
  
  # Fixed effects
  beta <- c(0.1, 0.2, 0, 0.3, 0.3, 0, rep(0, p - 5))

  # Random effects per grp (q × N)
  bi <- rbind(rnorm(N, 0, 0.3),                 # intercept
               matrix(0, q - 3, N),              # zero middle
               rnorm(N, 0, 0.6),                 # slope q-1
               rnorm(N, 0, 0.4))                 # slope q
  
  
  # Linear predictor and outcome
  linpred <- numeric(N)
  for (i in 1:N) {
    b_i <- bi[, i]
    linpred[i] <- x[i, ] %*% beta + z[i, ] %*% b_i
  }
  
  # Simulate survival times
  U <- runif(N)
  true_time <- pmax(-log(U) / (lambda0 * exp(linpred)), 0.1)
  
  # Apply censoring
  censor_time <- rexp(N, rate = censor_rate)
  observed_time <- pmin(true_time, censor_time)
  event <- as.numeric(true_time <= censor_time)
  
  # Create combined data frame
  df <- data.frame(
    time = observed_time,
    status = event,
    x,
    z,
    grp = grp
  )
  
  return(df)
}


# === Run Simulation ===
set.seed(1613)
data2 <- sim_cox_no_clusters()

# My Cox data
cox_data <- data2


# Save
save(data2, file = "sim_cox_no_clusters.RData")
write.csv(cox_data, "sim_cox_no_clusters.csv", row.names = FALSE)

# === Kaplan-Meier by Cluster ===
surv_obj <- Surv(cox_data$time, cox_data$status)
km_fit <- survfit(surv_obj ~ 1, data = cox_data)

ggsurvplot(
  km_fit,
  data = cox_data,
  risk.table = TRUE,
  pval = TRUE,
  surv.median.line = "hv",
  ggtheme = theme_minimal(),
)

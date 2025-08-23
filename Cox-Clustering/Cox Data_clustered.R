rm(list = ls())

library(mvtnorm)
library(survival)
library(survminer)
library(ggplot2)

library(mvtnorm)

sim_cox_fixed_clusters <- function(N1 = 50, N2 = 50, p = 50, q = 10, rho = 0,
                                   lambda0 = 0.5, censor_rate = 0.01) {
  
  n_subjects <- N1 + N2
  grp <- factor(1:n_subjects)
  cluster_label <- c(rep("One", N1), rep("Two", N2))
  
  # AR(1) covariance matrix
  covmat <- matrix(0, p, p)
  for (i in 1:p) for (j in 1:p) covmat[i, j] <- rho^abs(i - j)
  
  # Covariates (with intercept)
  # x <- cbind(1, rmvnorm(n_subjects, mean = rep(1, p), sigma = covmat))
  # colnames(x) <- c("Intercept", paste0("x.", 1:p))
  
  x <- cbind(1, rmvnorm(n_subjects, mean = rep(1, p), sigma = covmat))
  colnames(x) <- paste0("x.", 1:(p + 1))  # intercept becomes x.1
  
  # Random effects design matrix
  z <- x[, 1:q, drop = FALSE]
  colnames(z) <- paste0("z.", 1:q)
  
  # Fixed effects
  beta1 <- c(0.1, 0.2, 0, 0.3, 0.3, 0, rep(0, p - 5))
  beta2 <- c(-0.1, -0.2, 0, -0.3, -0.3, 0, rep(0, p - 5))
  
  # Random effects per cluster (q × N1/N2)
  bi1 <- rbind(rnorm(N1, 0, 0.3),                 # intercept
               matrix(0, q - 3, N1),              # zero middle
               rnorm(N1, 0, 0.6),                 # slope q-1
               rnorm(N1, 0, 0.4))                 # slope q
  
  bi2 <- rbind(rnorm(N2, 0, 0.3),
               matrix(0, q - 3, N2),
               rnorm(N2, 0, 0.6),
               rnorm(N2, 0, 0.4))
  
  # Linear predictor and outcome
  linpred <- numeric(n_subjects)
  for (i in 1:n_subjects) {
    if (i <= N1) {
      beta <- beta1
      bi <- bi1[, i]
    } else {
      beta <- beta2
      bi <- bi2[, i - N1]
    }
    linpred[i] <- x[i, ] %*% beta + z[i, ] %*% bi
  }
  
  # Simulate survival times
  U <- runif(n_subjects)
  true_time <- pmax(-log(U) / (lambda0 * exp(linpred)), 0.1)
  
  # Apply censoring
  censor_time <- rexp(n_subjects, rate = censor_rate)
  observed_time <- pmin(true_time, censor_time)
  event <- as.numeric(true_time <= censor_time)
  
  # Create combined data frame
  df <- data.frame(
    time = observed_time,
    status = event,
    x,
    z,
    grp = grp,
    cluster = factor(cluster_label)
  )
  
  return(df)
}


# === Run Simulation ===
set.seed(1613)
data3 <- sim_cox_fixed_clusters()

# My Cox data
cox_data <- data3


# Save
save(data3, file = "sim_cox_fixed_clusters.RData")
write.csv(cox_data, "sim_cox_fixed_clusters.csv", row.names = FALSE)

# === Kaplan-Meier by Cluster ===
surv_obj <- Surv(cox_data$time, cox_data$status)
km_fit <- survfit(surv_obj ~ cluster, data = cox_data)

ggsurvplot(
  km_fit,
  data = cox_data,
  risk.table = TRUE,
  pval = TRUE,
  surv.median.line = "hv",
  ggtheme = theme_minimal(),
  legend.title = "Cluster",
  legend.labs = c("One", "Two"),
  palette = c("#E69F00", "#56B4E9")
)

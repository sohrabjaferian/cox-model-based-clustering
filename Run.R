library(survival)
library(glmnet)
library(mclust)

source("helpers_Cox.R", keep.source = TRUE)
source("clusterspCox.R", keep.source = TRUE)
source("spcoxControl.R", keep.source = TRUE)



# Load the RData file
load("sim_cox_fixed_clusters.RData")

# Extract variables from data3
time <- data3$time
event <- data3$status
grp <- data3$grp

length(grp)             
length(unique(grp)) 
length(levels(grp))


# Covariates
x <- as.matrix(data3[, grep("^x\\.", names(data3))])  # already includes x.1 (intercept)
z <- as.matrix(data3[, grep("^z\\.", names(data3))])  # already includes z.1 (intercept)
# Response
y <- Surv(time, event)

# Run the clustering model
fit <- clustercox(
  x = x,
  z = z,
  grp = grp,
  time = time,
  event = event,
  lam1 = 0.1,
  lam2 = 0.1,
  nCluster = 2,
  nonpen.b = 1,
  nonpen.L = 1,
  penalty.b = "lasso",
  penalty.L = "lasso",
  standardize = TRUE,
  control = spCoxControl()
)


# View the results
str(fit)


# 6. Inspect the result
print(fit$coefficients)
print(fit$membership)
print(fit$logLik)

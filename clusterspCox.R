clustercox <- function(x, z, grp, time, event,
                       lam1, lam2, nCluster,
                       nonpen.b = 1, nonpen.L = 1,
                       penalty.b = c("lasso","scad"),
                       penalty.L = c("lasso","scad"),
                       standardize = TRUE,
                       control = spCoxControl()) {
  
  set.seed(123)
  
  penalty.b <- match.arg(penalty.b)
  penalty.L <- match.arg(penalty.L)
  
  # ------ checks ------
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.data.frame(z)) z <- as.matrix(z)
  if (!is.matrix(x)) stop("x has to be a matrix or data frame")
  if (!is.matrix(z)) stop("z has to be a matrix or data frame")
  if (any(is.na(x)) || any(is.na(z))) stop("Missing values in x/z not allowed")
  
  if (!is.numeric(time) || !is.numeric(event)) stop("time and event must be numeric")
  if (length(time) != length(event)) stop("time and event must have equal length")
  y <- Surv(time, event)
  if (any(is.na(y))) stop("Missing values in y not allowed")
  
  if (nrow(x) != nrow(y)) stop("x and y have not correct dimensions")
  if (any(x[,1] != 1)) stop("first column of x must be intercept 1s")
  
  if (!all(nonpen.b %in% seq_len(ncol(x)))) stop("bad nonpen.b")
  if (!all(nonpen.L %in% seq_len(ncol(z)))) stop("bad nonpen.L")
  
  if (any(lam1 < 0) || any(lam2 < 0)) stop("lam1/lam2 must be >= 0")
  if (length(lam1) == 1) lam1 <- rep(lam1, nCluster)
  if (length(lam2) == 1) lam2 <- rep(lam2, nCluster)
  
  # ------ standardize (except intercepts) ------
  if (standardize) {
    xOr <- x
    meanx <- apply(x[, -1, drop = FALSE], 2, mean)
    sdx   <- apply(x[, -1, drop = FALSE], 2, sd)
    x <- cbind(1, scale(x[, -1, drop = FALSE], center = meanx, scale = sdx))
    
    zOr <- z
    meanz <- apply(z[, -1, drop = FALSE], 2, mean)
    sdz   <- apply(z[, -1, drop = FALSE], 2, sd)
    z <- cbind(1, scale(z[, -1, drop = FALSE], center = meanz, scale = sdz))
  }
  
  # ------ allocate ------
  grp <- factor(grp)
  N   <- length(unique(grp))      # subjects
  p   <- ncol(x)
  q   <- ncol(z)
  ntot <- nrow(x)
  
  # initial global scaling of lambdas (will be re-scaled by π_g later)
  lambda1 <- lam1 * (N / nCluster)
  lambda2 <- lam2 * (N / nCluster)
  
  # --- base per-subject lists: NO weighting/scaling here ---
  xGrp <- lapply(seq_len(nrow(x)), function(i) matrix(x[i, ], nrow = 1))
  zGrp <- lapply(seq_len(nrow(z)), function(i) matrix(z[i, ], nrow = 1))
  yGrp <- lapply(seq_len(nrow(x)), function(i) {
    m <- matrix(c(time[i], event[i]), nrow = 1)
    colnames(m) <- c("time","event"); m
  })
  
  # --- initial cluster labels via Mclust on X (no weights applied to data) ---
  # --- initial cluster labels via Mclust on X ---
  x_kmeans <- x[, -1, drop = FALSE]
  ini.fit  <- Mclust(as.matrix(x_kmeans), G = nCluster)
  
  # responsibilities: N x G -> G x N
  post <- pmax(ini.fit$z, 1e-12)                 # N x G
  post <- post / rowSums(post)                   # normalize per subject
  membership <- t(post)                          # G x N
  memb.prob  <- rowMeans(membership)
  
  betaStart  <- LStart <- DStart <- vector("list", nCluster)
  fctStart   <- numeric(nCluster)
  
  pen_obj <- function(beta, L, w, lam1_i, lam2_i) {
    D <- L %*% t(L)
    ll <- cox_laplace_loglik(xGrp, yGrp, zGrp, beta, D, wGroup = w)$loglik
    pen_b <- if (penalty.b == "lasso") lam1_i * sum(abs(beta[-nonpen.b])) else sum(scad_group(beta[-nonpen.b], lam1_i))
    L2 <- sum(apply(L[-nonpen.L, , drop = FALSE], 1, function(r) sqrt(sum(r^2))))
    pen_L <- if (penalty.L == "lasso") lam2_i * L2 else sum(scad_group(L2, lam2_i))
    -ll + pen_b + pen_L
  }
  
  cat("fitting ...\n")
  pb <- txtProgressBar(min = 0, max = nCluster, style = 3)
  for (g in 1:nCluster) {
    w_g <- as.numeric(membership[g, ])  # length N
    
    yi <- Surv(time, event)
    betaStart[[g]] <- tryCatch({
      fit <- glmnet(
        x = x[, -1, drop = FALSE], y = yi, family = "cox",
        lambda = lambda1[g],
        alpha  = if (penalty.b == "lasso") 1 else 0,
        weights = pmax(w_g, 1e-6)
      )
      c(0, as.numeric(fit$beta[, 1]))
    }, error = function(e) {
      fit2 <- tryCatch(coxph(yi ~ x[, -1], weights = pmax(w_g, 1e-6), ties = "breslow"),
                       error = function(e2) NULL)
      if (!is.null(fit2)) {
        b <- coef(fit2); bb <- rep(0, ncol(x) - 1)
        nm <- intersect(names(b), colnames(x)[-1])
        bb[match(nm, colnames(x)[-1])] <- b[nm]
        return(c(0, bb))
      }
      rep(0, ncol(x))
    })
    
    # init D from your updated helper (pass weights, don't pass zId/N if your new signature dropped them)
    covInit <- covStartingValues(xGrp, yGrp, zGrp, b = betaStart[[g]], wGroup = w_g)
    tau <- if (is.finite(covInit$tau)) covInit$tau else 1
    DStart[[g]] <- diag(tau, ncol(z))
    LStart[[g]] <- chol(DStart[[g]])
    
    fctStart[g] <- pen_obj(betaStart[[g]], LStart[[g]], w_g, lambda1[g], lambda2[g])
  }
  
  
  
  betaIter <- betaStart
  LIter    <- LStart
  DIter    <- DStart
  LvecIter <- lapply(LStart, function(L) L[lower.tri(L, TRUE)])
  
  fctIter   <- fctStart
  covIter   <- LvecIter
  converged <- 0
  counter   <- 0
  counterIn <- 0
  doAll     <- FALSE
  convFct2  <- -10
  
  repeat {
    if (!(counter < control$maxIter && (convFct2 < 0 || counter < 1))) break
    counter <- counter + 1
    
    betaOld <- betaIter; LOld <- LIter; fctOld <- fctIter; covOld <- covIter
    
    activeSet <- lapply(betaIter, function(b) which(b != 0))
    if (counterIn == 0 || counterIn > control$number) {
      doAll <- TRUE
      activeSet <- lapply(activeSet, function(.) seq_len(p))
      counterIn <- 1
    } else {
      doAll <- FALSE
      counterIn <- counterIn + 1
    }
    
    # ------------------ M-step ------------------
    for (g in 1:nCluster) {
      w_g <- membership[g, ]
      
      Xst <- x
      Zst <- z
      Yst <- cbind(time, event)
      
      # curvature proxy for β: diagonal info (stable & >0)
      Hjj <- function(j) {
        # event-wise Var_w(x_j) accumulation with weights w_g
        timev <- Yst[,1]; eventv <- Yst[,2]
        ord <- order(timev, -eventv)
        r   <- exp(as.vector(Xst %*% betaIter[[g]]))
        r   <- r[ord]; w  <- w_g[ord]; xj <- Xst[ord, j]
        tme <- timev[ord]; evt <- eventv[ord]
        val <- 0
        for (i in which(evt == 1)) {
          Ri <- which(tme >= tme[i])
          wr <- w[Ri] * r[Ri]; den <- sum(wr); if (!is.finite(den) || den <= 0) next
          mu  <- sum(xj[Ri] * wr) / den
          e2  <- sum((xj[Ri]^2) * wr) / den
          val <- val + w[i] * (e2 - mu^2)
        }
        max(val, control$lower)
      }
      
      activeSet[[g]] <- if (doAll) seq_len(ncol(Xst)) else which(betaIter[[g]] != 0)
      for (j in activeSet[[g]]) {
        score <- cox_partial_score_component(x = Xst, y = Yst, beta = betaIter[[g]], j = j, w = w_g)
        Hj    <- Hjj(j)
        if (j %in% nonpen.b) {
          betaIter[[g]][j] <- score / Hj
        } else if (penalty.b == "lasso") {
          betaIter[[g]][j] <- SoftThreshold(score, lambda1[g]) / Hj
        } else {
          scada <- 3.7
          betaIter[[g]][j] <- SoftThreshold(score, lambda1[g]) / (Hj * (1 - 1 / scada))
        }
      }
      betaIter[[g]][!is.finite(betaIter[[g]])] <- 0
      betaIter[[g]][abs(betaIter[[g]]) < 0.05] <- 0
      
      
      # ---- L update (group penalty) ----
      D.grad <- D_Gradient(xGrp, zGrp, NULL, yGrp, b = betaIter[[g]], wGroup = w_g)
      D.hessian <- D_HessianMatrix(xGrp, zGrp, NULL, yGrp, b = betaIter[[g]], q = ncol(z), wGroup = w_g)
      L.grad <- t(LIter[[g]] %*% (D.grad + t(D.grad)))

      
      # simple diagonal curvature for stability
      for (k in 1:q) {
        for (l in k:q) {
          L.lk.grad  <- L.grad[l, k]
          L.lk.Hess  <- max(control$lower, min(control$upper, D.hessian[l, l]))
          
          linNonpen <- l %in% nonpen.L
          # one-step group update with pseudo-Armijo
          if (linNonpen) {
            dk <- - L.lk.grad / L.lk.Hess
          } else {
            row_norm <- sqrt(sum(LIter[[g]][l, ]^2))
            if (row_norm == 0) row_norm <- 1e-8
            if (penalty.L == "lasso") {
              dk <- (-L.lk.grad - lambda2[g] / row_norm * LIter[[g]][l, k]) /
                (L.lk.Hess + lambda2[g] / row_norm)
            } else {
              # SCAD group — approximate with same form but replace lambda2 with
              # its SCAD derivative scaling at row-norm.
              a <- 3.7
              deriv_scale <- if (row_norm <= lambda2[g]) 1 else if (row_norm <= a*lambda2[g]) (a*lambda2[g]-row_norm)/((a-1)*row_norm) else 0
              dk <- (-L.lk.grad - deriv_scale * LIter[[g]][l, k]) /
                (L.lk.Hess + max(deriv_scale, 1e-8))
            }
          }
          step <- control$a_init
          LIter[[g]][l, k] <- LIter[[g]][l, k] + step * dk
        }
      }
      
      LIter[[g]][abs(LIter[[g]]) < 1e-2] <- 0
      DIter[[g]]  <- LIter[[g]] %*% t(LIter[[g]])
      LvecIter[[g]] <- LIter[[g]][lower.tri(LIter[[g]], TRUE)]
      
      # new objective
      fctIter[g] <- pen_obj(betaIter[[g]], LIter[[g]], w_g, lambda1[g], lambda2[g])
    }
    
    # ------------------ E-step ------------------             # <<< CHANGED
    # p(y_i | X_i, Θ_g) for each i,g using *unweighted* subject i
    # --- E-step (compute responsibilities) ---
    log_ind <- sapply(seq_len(nCluster), function(g) {
      ll_i <- mapply(function(xi, yi, zi) {
        cox_laplace_loglik(list(xi), list(yi), list(zi),
                           betaIter[[g]], DIter[[g]])$loglik
      }, xGrp, yGrp, zGrp)
      log(pmax(memb.prob[g], 1e-12)) + ll_i
    })
    # log_ind is N x G
    
    # stabilize per SUBJECT (row-wise)
    log_ind  <- sweep(log_ind, 1, apply(log_ind, 1, max), "-")
    
    temp <- 1.5  # cool towards 1 over iterations if you like
    log_ind <- log_ind / temp

    
    ind.prob <- exp(log_ind)
    
    # normalize per SUBJECT (rows)
    post_ig     <- sweep(ind.prob, 1, rowSums(ind.prob), "/")  # N x G
    membership  <- t(post_ig)                                  # G x N (what the rest expects)
    memb.prob   <- rowMeans(membership)                        # mixture weights π_g

    
    
    # ---- reinit near-empty components (put this block here) ----
    eps_pi <- 1e-6
    for (g in 1:nCluster) {
      if (memb.prob[g] < eps_pi) {
        ref <- which.max(memb.prob)          # pick a healthy component as template
        betaIter[[g]] <- betaIter[[ref]] + rnorm(length(betaIter[[ref]]), 0, 0.1)
        LIter[[g]]    <- diag(diag(LIter[[ref]])) * runif(1, 0.8, 1.2)
        DIter[[g]]    <- LIter[[g]] %*% t(LIter[[g]])
        
        # seed some soft responsibilities for this component
        membership[g, ] <- eps_pi
        seeds <- sample.int(ncol(membership), min(5, ncol(membership)))
        membership[g, seeds] <- 1
      }
    }
    # renormalize responsibilities and update mixture weights
    membership <- sweep(membership, 2, colSums(membership), "/")
    memb.prob  <- rowMeans(membership)
    
    
    
    # rescale lambdas by cluster mass
    # lambda1 <- rep(lam1, nCluster)  # or lam1 * (N/nCluster)
    # lambda2 <- rep(lam2, nCluster)  # or lam2 * (N/nCluster)
    # lambda1 <- rep_len(lam1, nCluster)
    # lambda2 <- rep_len(lam2, nCluster)
    lambda1_eff <- lambda1 * memb.prob
    lambda2_eff <- lambda2 * memb.prob
    
    
    # then use lambda1_eff[g], lambda2_eff[g] in the M-step,
    # but keep lambda1/lambda2 unchanged for returning in fit$lambda*
    
    # ------------- convergence checks -------------
    convPar <- max(mapply(function(a,b) sqrt(crossprod(a-b))/(1+sqrt(crossprod(a))),
                          betaIter, betaOld))
    convFct  <- abs((sum(fctOld) - sum(fctIter)) / (1 + abs(sum(fctIter))))
    convFct2 <- sum(fctIter) - sum(fctOld)
    convCov  <- max(mapply(function(a,b) sqrt(crossprod(a-b)), LvecIter, covOld))
    
    if (!any(is.na(c(convPar, convFct, convCov))) &&
        convPar <= control$tol && convFct <= control$tol && convCov <= control$tol) {
      counterIn <- 0
    }
    
    if (counter >= control$maxIter) break
  } # end repeat
  
  # ------ unstandardize coefficients back to original scale ------
  if (standardize) {
    betaIter <- lapply(betaIter, function(b) {
      b[-1] <- b[-1] / sdx
      b[1]  <- b[1] - sum(meanx * b[-1])
      b
    })
    x <- xOr; z <- zOr
  }
  
  # ------ final model fit metrics (using weighted objective) ------
  npar <- sum(unlist(lapply(betaIter, function(b) sum(b != 0)))) + length(unlist(LvecIter))
  
  logLik <- sum(mapply(function(beta, L, w) {
    cox_laplace_loglik(
      xGroup = xGrp, yGroup = yGrp, zGroup = zGrp,
      beta = beta, D = L %*% t(L), wGroup = w
    )$loglik
  }, betaIter, LIter, as.data.frame(t(membership)) ))          # membership rows as w
  
  deviance <- -2 * logLik
  aic <- -2 * logLik + 2 * npar
  bic <- -2 * logLik + log(ntot) * npar
  
  p.nz <- sum(unlist(lapply(betaIter, function(b) sum(b != 0))))
  q.nz <- sum(unlist(lapply(DIter, function(D) sum(diag(D) != 0))))
  bbic <- -2 * logLik + max(1, log(log(p.nz + q.nz))) * log(ntot) * npar
  ebic <- -2 * logLik + (log(ntot) + 2 * log(p.nz + q.nz)) * npar
  
  out <- list(
    data = list(x = x, y = y, z = z, grp = grp),
    membership = membership,
    penalty.b = penalty.b, penalty.L = penalty.L,
    nonpen.b = nonpen.b, nonpen.L = nonpen.L,
    lambda1 = lambda1, lambda2 = lambda2,
    Lvec = LvecIter,
    coefficients = betaIter, D = DIter,
    converged = converged,
    logLik = logLik, npar = npar,
    deviance = deviance, aic = aic, bic = bic, bbic = bbic, ebic = ebic,
    counter = counter, control = control, call = match.call()
  )
  structure(out, class = "spcox")
}

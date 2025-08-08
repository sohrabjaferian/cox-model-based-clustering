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
  
  # ---- defaults for new control fields ----
  ctrl <- modifyList(list(
    tol = 1e-4,
    lower = 1e-8, upper = 1e8,
    maxIter = 10,
    inner_maxit = 10,
    armijo_c = 1e-4,          # sufficient decrease
    armijo_rho = 0.5,         # backtracking shrink
    a_init = 1.0,             # initial step
    number = 5                # your active-set cycling trigger
  ), control)
  
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
  N   <- nrow(x)                 # per-row subjects in your current setup
  p   <- ncol(x)
  q   <- ncol(z)

  # initial global scaling of lambdas; later rescaled by pi_g
  lambda1_base <- lam1 * (N / nCluster)
  lambda2_base <- lam2 * (N / nCluster)
  
  # per-"subject" lists used by your cox_laplace_loglik
  xGrp <- lapply(seq_len(N), function(i) matrix(x[i, ], nrow = 1))
  zGrp <- lapply(seq_len(N), function(i) matrix(z[i, ], nrow = 1))
  yGrp <- lapply(seq_len(N), function(i) {
    m <- matrix(c(time[i], event[i]), nrow = 1)
    colnames(m) <- c("time","event"); m
  })
  
  # ---- initialization via Mclust on X (no weights) ----
  ini.fit  <- mclust::Mclust(as.matrix(x[, -1, drop = FALSE]), G = nCluster)
  post     <- pmax(ini.fit$z, 1e-12)
  post     <- post / rowSums(post)
  membership <- t(post)                      # G x N
  memb.prob  <- rowMeans(membership)        # pi_g
  
  # ---- helpers ----
  scad_deriv <- function(t, lam, a = 3.7) {
    if (t <= lam) return(lam)
    if (t <= a * lam) return((a * lam - t) / (a - 1))
    0
  }
  
  pen_obj <- function(beta, L, w, lam1_i, lam2_i) {
    D <- L %*% t(L)
    ll <- cox_laplace_loglik(xGrp, yGrp, zGrp, beta, D, wGroup = w)$loglik
    # beta penalty (intercept non-penalized)
    if (penalty.b == "lasso") {
      pen_b <- lam1_i * sum(abs(beta[-nonpen.b]))
    } else {
      # LLA on SCAD: use derivative*|beta| surrogate (classic local linear approx)
      pen_b <- sum(sapply(setdiff(seq_along(beta), nonpen.b), function(j) {
        wj <- scad_deriv(abs(beta[j]), lam1_i)
        wj * abs(beta[j])
      }))
    }
    # group penalty on rows of L (by row ℓ2 norm), with nonpen rows exempt
    row_norms <- apply(L, 1, function(r) sqrt(sum(r^2)))
    if (penalty.L == "lasso") {
      pen_L <- lam2_i * sum(row_norms[setdiff(seq_len(q), nonpen.L)])
    } else {
      pen_L <- sum(sapply(setdiff(seq_len(q), nonpen.L), function(l) {
        wl <- scad_deriv(row_norms[l], lam2_i)
        wl * row_norms[l]
      }))
    }
    -(ll) + pen_b + pen_L
  }
  
  # curvature proxy for β_j (positive)
  Hjj_fun <- function(j, Xst, Yst, beta, w) {
    timev <- Yst[,1]; eventv <- Yst[,2]
    ord <- order(timev, -eventv)
    r   <- exp(as.vector(Xst %*% beta))
    r   <- r[ord]; w  <- w[ord]; xj <- Xst[ord, j]
    tme <- timev[ord]; evt <- eventv[ord]
    val <- 0
    for (i in which(evt == 1)) {
      Ri <- which(tme >= tme[i])
      wr <- w[Ri] * r[Ri]; den <- sum(wr); if (!is.finite(den) || den <= 0) next
      mu  <- sum(xj[Ri] * wr) / den
      e2  <- sum((xj[Ri]^2) * wr) / den
      val <- val + w[i] * (e2 - mu^2)
    }
    max(val, ctrl$lower)
  }
  
  # ---- initialize beta/L/D per cluster (like before) ----
  betaIter <- LIter <- DIter <- vector("list", nCluster)
  for (g in seq_len(nCluster)) {
    w_g <- as.numeric(membership[g, ])
    yi  <- Surv(time, event)
    betaIter[[g]] <- tryCatch({
      fit <- glmnet::glmnet(
        x = x[, -1, drop = FALSE], y = yi, family = "cox",
        lambda = lambda1_base[g], alpha = if (penalty.b == "lasso") 1 else 0,
        weights = pmax(w_g, 1e-6)
      )
      c(0, as.numeric(fit$beta[, 1]))
    }, error = function(e) {
      fit2 <- tryCatch(survival::coxph(yi ~ x[, -1], weights = pmax(w_g, 1e-6), ties = "breslow"),
                       error = function(e2) NULL)
      if (!is.null(fit2)) {
        b <- stats::coef(fit2); bb <- rep(0, ncol(x) - 1)
        nm <- intersect(names(b), colnames(x)[-1])
        bb[match(nm, colnames(x)[-1])] <- b[nm]
        return(c(0, bb))
      }
      rep(0, ncol(x))
    })
    
    covInit <- covStartingValues(xGrp, yGrp, zGrp, b = betaIter[[g]], wGroup = w_g)
    tau <- if (is.finite(covInit$tau)) covInit$tau else 1
    DIter[[g]] <- diag(tau, ncol(z))
    LIter[[g]] <- chol(DIter[[g]])
  }
  
  # ---------- EM loop ----------
  outer <- 0
  repeat {
    outer <- outer + 1
    
    # scale lambdas by cluster mass (π_g)
    lambda1_eff <- lambda1_base * memb.prob
    lambda2_eff <- lambda2_base * memb.prob
    
    fct_before <- sum(mapply(function(b, L, w, l1, l2)
      pen_obj(b, L, w, l1, l2),
      betaIter, LIter, as.data.frame(t(membership)), lambda1_eff, lambda2_eff))
    
    # ------------------ M-step: for each g, fully CGD until convergence ------------------
    for (g in seq_len(nCluster)) {
      w_g <- as.numeric(membership[g, ])
      Xst <- x; Zst <- z; Yst <- cbind(time, event)
      beta_g <- betaIter[[g]]
      L_g    <- LIter[[g]]
      
      inner <- 0
      repeat {
        inner <- inner + 1
        beta_old <- beta_g
        L_old    <- L_g
        
        # ---- β update (coordinate-wise with Armijo) ----
        active_b <- seq_len(p)  # doAll every time to match spec
        for (j in active_b) {
          # gradient component (score) and curvature proxy
          score_j <- cox_partial_score_component(x = Xst, y = Yst, beta = beta_g, j = j, w = w_g)
          Hj      <- Hjj_fun(j, Xst, Yst, beta_g, w_g)
          
          # proximal "full step" target for this coord
          beta_star_j <- if (j %in% nonpen.b) {
            beta_g[j] + score_j / Hj
          } else if (penalty.b == "lasso") {
            SoftThreshold(beta_g[j] + score_j / Hj, lambda1_eff[g] / Hj)
          } else {
            # SCAD-LLA: weight = scad_deriv(|beta_j|, lambda), prox-l1 with that weight
            wj <- scad_deriv(abs(beta_g[j]), lambda1_eff[g])
            SoftThreshold(beta_g[j] + score_j / Hj, wj / Hj)
          }
          
          d_j <- beta_star_j - beta_g[j]
          if (!is.finite(d_j) || d_j == 0) next
          
          # Armijo backtracking along this coordinate
          a <- ctrl$a_init
          obj0 <- pen_obj(beta_g, L_g, w_g, lambda1_eff[g], lambda2_eff[g])
          repeat {
            beta_try <- beta_g
            beta_try[j] <- beta_g[j] + a * d_j
            obj1 <- pen_obj(beta_try, L_g, w_g, lambda1_eff[g], lambda2_eff[g])
            if (is.finite(obj1) && obj1 <= obj0 - ctrl$armijo_c * a * abs(obj0)) break
            a <- a * ctrl$armijo_rho
            if (a < 1e-8) break
          }
          beta_g[j] <- beta_g[j] + a * d_j
        }
        
        # ---- L update (row-wise group with Armijo) ----
        # gradients/Hessian proxies from your helpers
        D.grad     <- D_Gradient(xGrp, zGrp, NULL, yGrp, b = beta_g, wGroup = w_g)
        D.hessian  <- D_HessianMatrix(xGrp, zGrp, NULL, yGrp, b = beta_g, q = ncol(z), wGroup = w_g)
        L_grad     <- t(L_g %*% (D.grad + t(D.grad)))   # same as before
        
        for (l in seq_len(q)) {
          g_l   <- as.numeric(L_grad[l, ])
          # curvature proxy for row l: positive scalar (use diag block or safe lower bound)
          H_l   <- max(ctrl$lower, min(ctrl$upper, D.hessian[l, l]))
          
          # unpenalized gradient step for the whole row
          l_tilde <- L_g[l, ] - (1 / H_l) * g_l
          
          # group shrink (lasso or SCAD-LLA on row norm)
          if (l %in% nonpen.L) {
            l_star <- l_tilde
          } else {
            rn <- sqrt(sum(l_tilde^2))
            if (penalty.L == "lasso") {
              shrink <- max(0, 1 - (lambda2_eff[g] / (H_l * rn)))
            } else {
              w_row  <- scad_deriv(sqrt(sum(L_g[l, ]^2)), lambda2_eff[g])
              shrink <- max(0, 1 - (w_row / (H_l * rn)))
            }
            l_star <- shrink * l_tilde
          }
          
          d_l <- l_star - L_g[l, ]
          if (!all(is.finite(d_l)) || all(abs(d_l) < 1e-15)) next
          
          # Armijo backtracking along the row direction
          a <- ctrl$a_init
          obj0 <- pen_obj(beta_g, L_g, w_g, lambda1_eff[g], lambda2_eff[g])
          repeat {
            L_try <- L_g
            L_try[l, ] <- L_g[l, ] + a * d_l
            obj1 <- pen_obj(beta_g, L_try, w_g, lambda1_eff[g], lambda2_eff[g])
            if (is.finite(obj1) && obj1 <= obj0 - ctrl$armijo_c * a * abs(obj0)) break
            a <- a * ctrl$armijo_rho
            if (a < 1e-8) break
          }
          L_g[l, ] <- L_g[l, ] + a * d_l
        }
        
        # check inner convergence
        max_beta <- max(abs(beta_g - beta_old))
        max_L    <- max(abs(L_g - L_old))
        if (max(max_beta, max_L) <= ctrl$tol || inner >= ctrl$inner_maxit) break
      }
      
      betaIter[[g]] <- beta_g
      LIter[[g]]    <- L_g
      DIter[[g]]    <- L_g %*% t(L_g)
    } # end M-step
    
    # ------------------ E-step (exact, no temperature, no reseed) ------------------
    log_ind <- sapply(seq_len(nCluster), function(g) {
      ll_i <- mapply(function(xi, yi, zi) {
        cox_laplace_loglik(list(xi), list(yi), list(zi),
                           betaIter[[g]], DIter[[g]])$loglik
      }, xGrp, yGrp, zGrp)
      log(pmax(memb.prob[g], 1e-12)) + ll_i
    })               # N x G
    # stabilize and normalize row-wise
    log_ind <- sweep(log_ind, 1, apply(log_ind, 1, max), "-")
    ind.prob <- exp(log_ind)
    post_ig  <- sweep(ind.prob, 1, rowSums(ind.prob), "/")   # N x G
    membership <- t(post_ig)                                  # G x N
    memb.prob  <- rowMeans(membership)
    
    # convergence check on penalized objective
    fct_after <- sum(mapply(function(b, L, w, l1, l2)
      pen_obj(b, L, w, l1, l2),
      betaIter, LIter, as.data.frame(t(membership)), 
      lambda1_base * memb.prob, lambda2_base * memb.prob))
    
    rel_dec <- abs(fct_before - fct_after) / (1 + abs(fct_after))
    if (!is.finite(rel_dec) || rel_dec < ctrl$tol || outer >= ctrl$maxIter) break
  } # EM repeat
  
  # ------ unstandardize coefficients back to original scale ------
  if (standardize) {
    betaIter <- lapply(betaIter, function(b) {
      b[-1] <- b[-1] / sdx
      b[1]  <- b[1] - sum(meanx * b[-1])
      b
    })
    x <- xOr; z <- zOr
  }
  
  # ------ final metrics ------
  npar <- sum(unlist(lapply(betaIter, function(b) sum(b != 0)))) +
    sum(unlist(lapply(LIter, function(L) length(L))))
  
  logLik <- sum(mapply(function(beta, L, w) {
    cox_laplace_loglik(xGroup = xGrp, yGroup = yGrp, zGroup = zGrp,
                       beta = beta, D = L %*% t(L), wGroup = w)$loglik
  }, betaIter, LIter, as.data.frame(t(membership))))
  
  deviance <- -2 * logLik
  aic <- -2 * logLik + 2 * npar
  bic <- -2 * logLik + log(N) * npar
  
  out <- list(
    data = list(x = x, y = y, z = z, grp = grp),
    membership = membership, pi = memb.prob,
    penalty.b = penalty.b, penalty.L = penalty.L,
    nonpen.b = nonpen.b, nonpen.L = nonpen.L,
    lambda1 = lambda1_base, lambda2 = lambda2_base,
    coefficients = betaIter, L = LIter, D = lapply(LIter, function(L) L %*% t(L)),
    logLik = logLik, npar = npar, deviance = deviance, aic = aic, bic = bic,
    iters = outer, control = ctrl, call = match.call()
  )
  structure(out, class = "spcox")
}

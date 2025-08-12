# multiplication <- function(x,y){
#   x*y
# }

multiplication <- function(x, y) {
  if (inherits(x, "Surv") && is.numeric(y)) {
    x_mat <- as.matrix(x)
    scaled_time <- pmax(x_mat[, 1] * y, 1e-2)  # cap too-small times
    return(Surv(time = scaled_time, event = x_mat[, 2]))
  } else if (inherits(y, "Surv") && is.numeric(x)) {
    y_mat <- as.matrix(y)
    return(Surv(time = y_mat[, 1] * x, event = y_mat[, 2]))
  } else if (is.numeric(y) && is.matrix(x)) {
    return(x * y)
  } else if (is.numeric(y) && is.vector(x)) {
    return(x * y)
  } else if (is.list(x) && is.numeric(y)) {
    return(lapply(x, function(elem) multiplication(elem, y)))
  } else if (is.numeric(x) && is.list(y)) {
    return(lapply(y, function(elem) multiplication(x, elem)))
  } else {
    stop(paste("Unsupported types in multiplication:", class(x), "and", class(y)))
  }
}






# ================================
# 1) frailty_gradient_hessian
# ================================
# Gradient/Hessian of the weighted Cox partial loglik w.r.t. u
# (with linear predictor eta = X b + Z u), Breslow handling.
frailty_gradient_hessian <- function(x, y, z, beta, u,
                                     w = NULL, ridge = 1e-10) {
  if (!is.matrix(x) || !is.matrix(z) || !is.matrix(y))
    stop("x, y, z must be matrices")
  if (ncol(y) < 2) stop("y must have columns: time, event")
  
  n  <- nrow(x)
  q  <- ncol(z)
  if (is.null(w)) w <- rep(1, n)
  
  time  <- y[, 1]
  event <- y[, 2]
  
  # Order: increasing time, events first (Breslow-like)
  ord   <- order(time, -event)
  x <- x[ord, , drop = FALSE]
  z <- z[ord, , drop = FALSE]
  time  <- time[ord]
  event <- event[ord]
  w     <- w[ord]
  
  eta  <- as.vector(x %*% beta + z %*% u)
  r    <- exp(eta)
  
  grad <- numeric(q)
  hess <- matrix(0, q, q)
  
  ev_idx <- which(event == 1)
  for (i in ev_idx) {
    Ri <- which(time >= time[i])
    
    wr    <- w[Ri] * r[Ri]
    denom <- sum(wr)
    if (!is.finite(denom) || denom <= 0) next
    
    Zi <- z[Ri, , drop = FALSE]
    
    # E_w[Z] and E_w[ZZ^T] under weights proportional to wr
    mu  <- colSums(Zi * wr) / denom
    Ezz <- crossprod(Zi, Zi * (wr / denom))
    
    # Score and observed info in u-direction
    grad <- grad + w[i] * (z[i, ] - mu)
    hess <- hess - w[i] * (Ezz - tcrossprod(mu))
  }
  
  # Symmetrize / tiny ridge
  hess <- (hess + t(hess)) / 2
  if (any(!is.finite(hess))) {
    hess[!is.finite(hess)] <- 0
    hess <- hess + diag(ridge, q)
  }
  
  list(grad = grad, hess = hess, eta = eta)
}


# ==========================================
# 2) cox_partial_score_component (for beta_j)
# ==========================================
cox_partial_score_component <- function(x, y, beta, j, w = NULL) {
  if (!is.matrix(x) || !is.matrix(y)) stop("x and y must be matrices")
  if (ncol(y) < 2) stop("y must have columns: time, event")
  
  n <- nrow(x)
  if (is.null(w)) w <- rep(1, n)
  
  time  <- y[, 1]
  event <- y[, 2]
  
  ord   <- order(time, -event)
  x <- x[ord, , drop = FALSE]
  time  <- time[ord]
  event <- event[ord]
  w     <- w[ord]
  
  eta <- as.vector(x %*% beta)
  r   <- exp(eta)
  
  score <- 0
  ev_idx <- which(event == 1)
  for (i in ev_idx) {
    Ri <- which(time >= time[i])
    wr <- w[Ri] * r[Ri]
    
    denom <- sum(wr)
    if (!is.finite(denom) || denom <= 0) next
    
    Ezj <- sum(x[Ri, j] * wr) / denom
    score <- score + w[i] * (x[i, j] - Ezj)
  }
  score
}


# ==================
# 3) matsplitter
# ==================
# Split a (q^2 x q^2) matrix into a list of q^2 blocks, each (q x q).
# This is strict on dimensions to avoid silent misuse.
matsplitter <- function(M, q) {
  nr <- nrow(M); nc <- ncol(M)
  if (nr != q * q || nc != q * q) {
    stop(sprintf("matsplitter: expected a %d x %d matrix; got %d x %d",
                 q*q, q*q, nr, nc))
  }
  blocks <- vector("list", q * q)
  for (i in seq_len(q)) {
    for (j in seq_len(q)) {
      r <- ((i - 1) * q + 1):(i * q)
      c <- ((j - 1) * q + 1):(j * q)
      blocks[[ (i - 1) * q + j ]] <- M[r, c, drop = FALSE]
    }
  }
  blocks
}


# =========================
# 4) covStartingValues (Cox)
# =========================
# Scalar frailty variance init via Laplace-approximated Cox marginal loglik.
# Parameterization: D = tau^2 * I_q with tau = exp(gamma).
covStartingValues <- function(xGroup, yGroup, zGroup, zIdGroup,  # zIdGroup kept for signature compatibility
                              b, N,
                              wGroup = NULL,
                              lower = -6, upper = 2,
                              ridge = 1e-6) {
  if (length(zGroup) == 0L || ncol(zGroup[[1]]) == 0L) {
    return(list(tau = 0, sigma = NA_real_, opt = NA_real_))
  }
  q <- ncol(zGroup[[1]])
  
  # Objective for optimize(): minimize negative Laplace loglik
  nll <- function(gamma) {
    tau2 <- exp(2 * gamma)
    D    <- diag(tau2, q)
    out  <- tryCatch(
      cox_laplace_loglik(
        xGroup = xGroup, yGroup = yGroup, zGroup = zGroup,
        beta = b, D = D, wGroup = wGroup, ridge = ridge
      )$loglik,
      error = function(e) -1e6
    )
    # We minimize
    if (!is.finite(out)) out <- -1e6
    -out
  }
  
  opt <- optimize(nll, interval = c(lower, upper))
  gamma_hat <- opt$minimum
  tau_hat   <- exp(gamma_hat)
  
  list(tau = tau_hat, sigma = NA_real_, opt = opt$objective)
}



ZIdentity <- function(Z) {
  diag(nrow(Z))  # not list()
}




nlogdet_Cox <- function(V_list) {
  sum(sapply(V_list, function(V) {
    -0.5 * determinant(V, logarithm = TRUE)$modulus[1]
  }))
}


VInv <- function(x, y, z, beta, D) {
  time <- y[,1]; event <- y[,2]
  eta  <- as.vector(x %*% beta)
  ord  <- order(time, -event)
  time <- time[ord]; event <- event[ord]; z <- z[ord,,drop=FALSE]; eta <- eta[ord]
  r    <- exp(eta)
  
  H <- matrix(0, ncol(z), ncol(z))
  for (i in which(event == 1)) {
    Ri <- which(time >= time[i])
    w  <- r[Ri]
    Zw <- sweep(z[Ri,,drop=FALSE], 1, w, `*`)
    mu <- colSums(Zw) / sum(w)                 # E_w[z]
    S2 <- crossprod(z[Ri,,drop=FALSE], Zw) / sum(w)   # E_w[zz^T]
    H  <- H + (S2 - tcrossprod(mu))            # Var_w(z)
  }
  H + solve(D)
}






VnotInv <- function(x, y, z, beta, D) {
  # Compute eta
  eta <- x %*% beta
  risk <- exp(eta)
  event <- y[, 2]
  
  # Compute the Hessian of the log-likelihood w.r.t. u
  W <- risk * cumsum(event / rev(cumsum(rev(risk))))  # approximate weights
  H_u <- t(z) %*% (W * z)  # observed Fisher information
  
  # Add the precision term from the prior: D^{-1}
  V_inv <- H_u + solve(D)
  
  # Return the covariance matrix (inverse of the above)
  V <- solve(V_inv)
  
  return(V)
}






nlogdet <- function(LGroup)
{
  nlogdetfun <- function(L)
  {
    -1/2*determinant(L)$modulus[1]
  }
  
  sum(mapply(nlogdetfun,LGroup))
}


# Efron

cox_laplace_loglik <- function(xGroup, yGroup, zGroup, beta, D,
                               wGroup = NULL, tol = 1e-6, maxiter = 50,
                               ridge = 1e-6,
                               ties = c("efron","breslow")) {
  ties <- match.arg(ties)
  
  # Stack per-subject lists
  X <- do.call(rbind, xGroup)
  Z <- do.call(rbind, zGroup)
  Y <- do.call(rbind, yGroup)
  if (!is.matrix(X) || !is.matrix(Z) || !is.matrix(Y))
    stop("xGroup/yGroup/zGroup must be lists of matrices that rbind cleanly.")
  
  n <- nrow(X); q <- ncol(Z)
  if (q == 0L || n == 0L) return(list(loglik = -1e6, uhat = rep(0, 0), Hessian = matrix(,0,0)))
  
  # Weights
  if (is.null(wGroup)) {
    w <- rep(1, n)
  } else if (is.list(wGroup)) {
    w <- as.numeric(unlist(wGroup))
  } else {
    w <- as.numeric(wGroup)
  }
  if (length(w) != n) {
    warning("wGroup length mismatch; using equal weights.")
    w <- rep(1, n)
  }
  
  # Order by time asc, break ties by events first
  time  <- Y[, 1]
  event <- Y[, 2]
  ord   <- order(time, -event)
  X <- X[ord, , drop = FALSE]
  Z <- Z[ord, , drop = FALSE]
  time  <- time[ord]
  event <- event[ord]
  w     <- w[ord]
  
  # Safe inverse for D
  D_inv <- tryCatch(solve(D), error = function(e) solve(D + diag(ridge, nrow(D))))
  
  u <- rep(0, q)
  
  # --- Efron/Breslow partial loglik + score/Hessian w.r.t u ---
  w_cox_stats <- function(u) {
    eta <- as.vector(X %*% beta + Z %*% u)
    r   <- exp(eta)
    wr  <- w * r
    
    loglik <- 0
    grad   <- numeric(q)
    hess   <- matrix(0, q, q)
    
    # Unique event times
    ev_times <- sort(unique(time[event == 1]))
    for (t in ev_times) {
      Et <- which(time == t & event == 1)        # indices of failures at t
      dt <- length(Et)
      Rt <- which(time >= t)                     # risk set just before t
      
      # risk-set sums (weighted by wr)
      R_denom <- sum(wr[Rt])
      if (!is.finite(R_denom) || R_denom <= 0) next
      
      Z_R <- Z[Rt, , drop = FALSE]
      Z_E <- Z[Et, , drop = FALSE]
      
      # first moments (numerators): sum wr * Z
      numZ_R <- colSums(Z_R * wr[Rt])
      numZ_E <- if (dt > 0) colSums(Z_E * wr[Et]) else rep(0, q)
      
      # second moments (numerators): sum wr * Z Z^T
      M2_R <- crossprod(Z_R, Z_R * wr[Rt])
      M2_E <- if (dt > 0) crossprod(Z_E, Z_E * wr[Et]) else matrix(0, q, q)
      
      # Numerator contribution from failures to loglik
      loglik <- loglik + sum(w[Et] * eta[Et])
      
      if (ties == "breslow" || dt == 1) {
        # Breslow (or single failure)
        denom <- R_denom
        if (!is.finite(denom) || denom <= 0) next
        
        EZ   <- numZ_R / denom
        EZZ  <- M2_R   / denom
        CovZ <- EZZ - tcrossprod(EZ)
        
        grad <- grad + colSums(Z[Et, , drop = FALSE] * w[Et]) - EZ
        hess <- hess - CovZ
        
        if (ties == "breslow") {
          # subtract dt * log(denom) for dt failures
          loglik <- loglik - dt * log(denom)
        } else { # dt==1 case already covered: subtract log(denom)
          loglik <- loglik - log(denom)
        }
      } else {
        # Efron: sum over k = 0..dt-1
        F_denom <- sum(wr[Et])  # sum of wr over the tied failures
        numZ_E2 <- numZ_E
        M2_E2   <- M2_E
        
        # add the score piece from failures once
        grad <- grad + colSums(Z[Et, , drop = FALSE] * w[Et])
        
        for (k in 0:(dt - 1)) {
          frac   <- k / dt
          denomk <- R_denom - frac * F_denom
          if (!is.finite(denomk) || denomk <= 0) next
          
          numZk  <- numZ_R - frac * numZ_E2
          M2k    <- M2_R   - frac * M2_E2
          
          EZk    <- numZk / denomk
          EZZk   <- M2k   / denomk
          CovZk  <- EZZk - tcrossprod(EZk)
          
          loglik <- loglik - log(denomk)
          grad   <- grad   - EZk
          hess   <- hess   - CovZk
        }
      }
    }
    
    list(loglik = loglik, grad = grad, hess = hess)
  }
  
  # Newton on log posterior: ℓ_partial(u) - 1/2 uᵀ D^{-1} u
  for (iter in seq_len(maxiter)) {
    s   <- w_cox_stats(u)
    g   <- s$grad - as.vector(D_inv %*% u)        # gradient of log posterior
    H   <- s$hess - D_inv                         # Hessian  of log posterior (negative-definite)
    
    if (any(!is.finite(g)) || any(!is.finite(H))) {
      warning("Non-finite grad/Hess in cox_laplace_loglik; aborting NR.")
      break
    }
    
    H_stable <- H + diag(ridge, q)
    step <- tryCatch(solve(H_stable, g), error = function(e) rep(NA_real_, q))
    if (anyNA(step)) {
      warning("solve() failed in cox_laplace_loglik; aborting NR.")
      break
    }
    
    u_new <- as.vector(u - step)
    if (max(abs(u_new - u)) < tol) { u <- u_new; break }
    u <- u_new
  }
  
  # Final Laplace pieces at û
  s_final <- w_cox_stats(u)
  H_neg <- -(s_final$hess - D_inv)               # -∇² log posterior at û
  
  if (any(!is.finite(H_neg))) H_neg <- H_neg + diag(ridge, q)
  
  logdetHu <- tryCatch(
    determinant(H_neg, logarithm = TRUE)$modulus[1],
    error = function(e) { warning("determinant(H_neg) failed; using fallback."); log(ridge) * q }
  )
  
  quad_pen <- as.numeric(t(u) %*% D_inv %*% u)
  laplace_approx <- as.numeric(s_final$loglik - 0.5 * quad_pen - 0.5 * logdetHu)
  if (!is.finite(laplace_approx)) laplace_approx <- -1e6
  
  list(loglik = laplace_approx, uhat = u, Hessian = H_neg)
}


## Breslow
# cox_laplace_loglik <- function(xGroup, yGroup, zGroup, beta, D,
#                                wGroup = NULL, tol = 1e-6, maxiter = 50,
#                                ridge = 1e-6) {
#   # Stack per-subject lists
#   X <- do.call(rbind, xGroup)
#   Z <- do.call(rbind, zGroup)
#   Y <- do.call(rbind, yGroup)
#   
#   if (!is.matrix(X) || !is.matrix(Z) || !is.matrix(Y))
#     stop("xGroup/yGroup/zGroup must be lists of matrices that rbind cleanly.")
#   
#   n <- nrow(X); q <- ncol(Z)
#   if (q == 0L || n == 0L) return(list(loglik = -1e6, uhat = rep(0, 0), Hessian = matrix(,0,0)))
#   
#   # Weights (soft memberships). Accept vector or list; default 1.
#   if (is.null(wGroup)) {
#     w <- rep(1, n)
#   } else if (is.list(wGroup)) {
#     w <- as.numeric(unlist(wGroup))
#   } else {
#     w <- as.numeric(wGroup)
#   }
#   if (length(w) != n) {
#     warning("wGroup length mismatch; using equal weights.")
#     w <- rep(1, n)
#   }
#   
#   # Order by time asc, break ties by events first (Breslow-like)
#   time  <- Y[, 1]
#   event <- Y[, 2]
#   ord   <- order(time, -event)
#   X <- X[ord, , drop = FALSE]
#   Z <- Z[ord, , drop = FALSE]
#   time  <- time[ord]
#   event <- event[ord]
#   w     <- w[ord]
#   
#   # Safe inverse for D
#   D_inv <- tryCatch(solve(D), error = function(e) solve(D + diag(ridge, nrow(D))))
#   
#   u <- rep(0, q)
#   
#   # Helper: compute weighted Cox loglik, grad_u, hess_u for current u
#   w_cox_stats <- function(u) {
#     eta <- as.vector(X %*% beta + Z %*% u)
#     r   <- exp(eta)
#     
#     loglik <- 0
#     grad   <- numeric(q)
#     hess   <- matrix(0, q, q)
#     
#     # Loop over event times; Breslow weighting with case-weights w
#     ev_idx <- which(event == 1)
#     for (i in ev_idx) {
#       Ri <- which(time >= time[i])
#       
#       wr   <- w[Ri] * r[Ri]
#       denom <- sum(wr)
#       if (!is.finite(denom) || denom <= 0) next
#       
#       Zi <- Z[Ri, , drop = FALSE]
#       mu <- colSums(Zi * wr) / denom                         # E_w[Z]
#       # E_w[ZZ^T]
#       Ezz <- crossprod(Zi, Zi * (wr / denom))
#       
#       loglik <- loglik + w[i] * (eta[i] - log(denom))
#       grad   <- grad + w[i] * (Z[i, ] - mu)
#       hess   <- hess - w[i] * (Ezz - tcrossprod(mu))
#     }
#     
#     list(loglik = loglik, grad = grad, hess = hess)
#   }
#   
#   # Newton on log posterior: ℓ(u) - 1/2 uᵀD^{-1}u
#   for (iter in seq_len(maxiter)) {
#     s   <- w_cox_stats(u)
#     g   <- s$grad - as.vector(D_inv %*% u)        # gradient of log posterior
#     H   <- s$hess - D_inv                         # Hessian  of log posterior (negative-definite)
#     
#     if (any(!is.finite(g)) || any(!is.finite(H))) {
#       warning("Non-finite grad/Hess in cox_laplace_loglik; aborting NR.")
#       break
#     }
#     
#     H_stable <- H + diag(ridge, q)
#     step <- tryCatch(solve(H_stable, g), error = function(e) rep(NA_real_, q))
#     if (anyNA(step)) {
#       warning("solve() failed in cox_laplace_loglik; aborting NR.")
#       break
#     }
#     
#     u_new <- as.vector(u - step)
#     if (max(abs(u_new - u)) < tol) { u <- u_new; break }
#     u <- u_new
#   }
#   
#   # Final stats at û
#   s_final <- w_cox_stats(u)
#   # Negative Hessian of log posterior at û (must be PD for Laplace)
#   H_neg <- -(s_final$hess - D_inv)
#   
#   # Stabilize if necessary
#   if (any(!is.finite(H_neg))) H_neg <- H_neg + diag(ridge, q)
#   
#   logdetHu <- tryCatch(
#     determinant(H_neg, logarithm = TRUE)$modulus[1],
#     error = function(e) { warning("determinant(H_neg) failed; using fallback."); log(ridge) * q }
#   )
#   
#   quad_pen <- as.numeric(t(u) %*% D_inv %*% u)
#   
#   laplace_approx <- as.numeric(s_final$loglik - 0.5 * quad_pen - 0.5 * logdetHu)
#   
#   if (!is.finite(laplace_approx)) laplace_approx <- -1e6
#   
#   list(loglik = laplace_approx, uhat = u, Hessian = H_neg)
# }








  


quad.form.inv <- function(A, x) {
  A_reg <- A + diag(1e-6, nrow(A))  # small ridge for numerical stability
  return(c(crossprod(x, solve(A_reg, x))))
}



cox_partial_loglik <- function(y, eta) {
  event <- y[, 2]
  risk <- exp(eta)
  log_risk_sum <- log(rev(cumsum(rev(risk))))
  return(sum(event * (eta - log_risk_sum)))
}











ObjFunction <- function(xGroup, yGroup, zGroup, beta, L,
                           lambda1, lambda2,
                           nonpen.b, nonpen.L,
                           penalty_b, penalty_L) {
  D <- L %*% t(L)
  
  laplace_res <- cox_laplace_loglik(xGroup, yGroup, zGroup, beta, D)
  
  loglik <- laplace_res$loglik
  
  # Penalties
  pen_b <- if (penalty_b == "lasso") {
    lambda1 * sum(abs(beta[-nonpen.b]))
  } else sum(scad_group(beta[-nonpen.b], lambda1))
  
  L2norm <- sum(apply(L[-nonpen.L, , drop = FALSE], 1, function(x) sqrt(sum(x^2))))
  pen_L <- if (penalty_L == "lasso") {
    lambda2 * L2norm
  } else sum(scad_group(L2norm, lambda2))
  
  return(-loglik + pen_b + pen_L)
}


as1 <- function(xGroup,LGroup,activeSet,N)
{
  fs <- function(x,l,a) {l%*%x[,a]}
  
  SGroup <- mapply(fs,xGroup,LGroup,MoreArgs=list(a=activeSet),SIMPLIFY=FALSE)
  
  return(SGroup)
}

as2 <- function(x,y,b,j,activeSet,group,sGroup)
{
  r <- y-x[,-c(j),drop=FALSE]%*%b[-j]
  rGroup <- split(r,group)
  
  as3 <- function(s,r,j) crossprod(r,s[,j])
  
  ma <- mapply(as3,sGroup,rGroup,MoreArgs=list(j=match(j,activeSet)))
  
  sumMa <- sum(ma)
  
  return(sumMa)
}

SoftThreshold <- function(z,g)
{
  sign(z)*max(abs(z)-g,0)
}












ArmijoRule_b <- function(x, y, z, b, j, cut, HkOldJ, HkJ, JinNonpen,
                         lambda, nonpen, penalty, L, converged, control)
{
  b.new <- b
  bJ <- b[j]
  grad <- -cut + bJ * HkOldJ
  
  # Objective value at current b
  fctOld <- objFunction_L(
    xGroup = list(x), yGroup = list(y), zGroup = list(z),
    zIdGrp = list(zId), b = b, L = L, lambda = lambda,
    nonpen = nonpen, penalty = penalty, ll1 = ll1, ntot = ntot
  )
  
  # Compute descent direction
  if (JinNonpen) {
    dk <- -grad / HkJ
  } else {
    if (penalty == "lasso") {
      dk <- MedianValue(grad, HkJ, lambda, bJ)
    } else if (penalty == "scad") {
      dk <- ScadValue(grad, HkJ, lambda, bJ)
    }
  }
  
  # Check dk validity
  if (is.na(dk) || !is.finite(dk)) {
    warning("dk is NA or not finite — skipping Armijo step for variable j = ", j)
    return(list(b = b, fct = fctOld, converged = converged + 1))
  }
  
  # Recompute objective
  fctOld <- objFunction_b(
    xGroup = list(x), yGroup = list(y), zGroup = list(z),
    beta = b, L = L, lambda = lambda, 
    nonpen = nonpen, penalty = penalty
  )
  
  if (!is.na(dk) && dk != 0) {
    if (JinNonpen) {
      deltak <- dk * grad + control$gamma * dk^2 * HkJ
    } else {
      deltak <- dk * grad + control$gamma * dk^2 * HkJ + lambda * (abs(bJ + dk) - abs(bJ))
    }
    
    for (l in 0:control$maxArmijo) {
      step_size <- control$a_init * control$delta^l
      b.new[j] <- bJ + step_size * dk
      
      fctNew <- objFunction_b(
        xGroup = list(x), yGroup = list(y), zGroup = list(z),
        beta = b.new, L = L, lambda = lambda,
        nonpen = nonpen, penalty = penalty
      )
      
      addDelta <- step_size * control$rho * deltak
      
      # Check for valid objective values
      if (any(is.na(c(fctNew, fctOld, addDelta))) || any(!is.finite(c(fctNew, fctOld, addDelta)))) {
        warning("Non-finite or NA values in Armijo condition — skipping update at j = ", j)
        return(list(b = b, fct = fctOld, converged = converged + 1))
      }
      
      if (fctNew <= fctOld + addDelta) {
        b[j] <- bJ + step_size * dk
        return(list(b = b, fct = fctNew, converged = converged))
      }
    }
    
    # No step size succeeded
    converged <- converged + 2
  }
  
  return(list(b = b, fct = fctOld, converged = converged))
}



objFunction_b <- function(xGroup, yGroup, zGroup, beta, L, lambda, nonpen, penalty) {
  # Construct D = L L^T
  D <- L %*% t(L)
  
  # Log-likelihood via Laplace approximation
  laplace_out <- cox_laplace_loglik(xGroup = xGroup, yGroup = yGroup, zGroup = zGroup, beta = beta, D = D)
  loglik <- laplace_out$loglik
  
  # Penalty term
  pen_term <- if (penalty == "lasso") {
    lambda * sum(abs(beta[-nonpen]))
  } else if (penalty == "scad") {
    sum(scad_group(beta[-nonpen], lambda))
  } else {
    0
  }
  
  return(-loglik + pen_term)
}



scad_group <- function(beta, lam1, scada = 3.7) {
  beta_orig <- beta
  beta <- abs(beta)
  
  s <- rep(0, length(beta))
  
  # Region 1: beta < lam1
  tmp1 <- which(!is.na(beta) & beta < lam1)
  s[tmp1] <- lam1 * beta[tmp1]
  
  # Region 2: lam1 < beta <= scada * lam1
  tmp2 <- which(!is.na(beta) & beta > lam1 & beta <= scada * lam1)
  s[tmp2] <- -(beta[tmp2]^2 - 2 * scada * lam1 * beta[tmp2] + lam1^2) / (2 * (scada - 1))
  
  # Region 3: beta > scada * lam1
  tmp3 <- which(!is.na(beta) & beta > scada * lam1)
  s[tmp3] <- (scada + 1) * lam1^2 / 2
  
  return(s)
}


scad = function(bj, lambda, a=3.7){
  temp = abs(bj)
  if(temp<=lambda){
    return(lambda)
  }else if(temp>lambda&temp<=a*lambda){
    return((a*lambda-bj)/(a-1))
  }else{
    return(0)
  }
}

MedianValue <- function(grad,hessian,lambda,bj)
{
  median(c((lambda-grad)/hessian,-bj,(-lambda-grad)/hessian))
}

ScadValue <- function(grad,hessian,lambda,bj, a=3.7)
{
  median(c((lambda-grad)/(hessian*(1-1/a)),-bj,(-lambda-grad)/hessian*(1-1/a)))
}


ResAsSplit <- function(x,y,b,f,activeset)
{
  r <- y-x[,activeset,drop=FALSE]%*%b[activeset,drop=FALSE]
  resGroup <- split(r,f)
  return(resGroup)
}

# ---- Utility: expand weights from subjects -> rows ----
.expand_weights <- function(wGroup, xGroup) {
  if (is.null(wGroup)) {
    # one row per subject assumed
    return(rep(1, sum(vapply(xGroup, nrow, 1L))))
  }
  if (is.list(wGroup)) {
    # list of scalars or vectors aligned with subjects
    unlist(mapply(function(wi, xi) {
      if (length(wi) == 1L) rep(wi, nrow(xi)) else wi
    }, wGroup, xGroup, SIMPLIFY = FALSE), use.names = FALSE)
  } else {
    # numeric vector; recycle if needed
    rep(wGroup, length.out = sum(vapply(xGroup, nrow, 1L)))
  }
}

# ---- Sum of event-wise Var_w(z) (observed info in u-direction) ----
D_Gradient <- function(xGroup, zGroup, LGroup = NULL, yGroup, b, N = NULL,
                       wGroup = NULL, verbose = FALSE, ridge = 1e-10) {
  # Stack lists
  X <- do.call(rbind, xGroup)
  Z <- do.call(rbind, zGroup)
  Y <- do.call(rbind, yGroup)
  
  if (!is.matrix(X) || !is.matrix(Z) || !is.matrix(Y))
    stop("xGroup/zGroup/yGroup must rbind to matrices.")
  
  q <- ncol(Z)
  if (q == 0L) return(matrix(0, 0, 0))
  
  w <- .expand_weights(wGroup, xGroup)
  
  time  <- Y[, 1]
  event <- Y[, 2]
  ord   <- order(time, -event)     # Breslow-like
  X <- X[ord, , drop = FALSE]
  Z <- Z[ord, , drop = FALSE]
  time  <- time[ord]
  event <- event[ord]
  w     <- w[ord]
  
  eta <- as.vector(X %*% b)
  r   <- exp(eta)
  
  G <- matrix(0, q, q)
  ev_idx <- which(event == 1)
  for (i in ev_idx) {
    Ri <- which(time >= time[i])
    
    wr    <- w[Ri] * r[Ri]
    denom <- sum(wr)
    if (!is.finite(denom) || denom <= 0) next
    
    Zi <- Z[Ri, , drop = FALSE]
    mu <- colSums(Zi * wr) / denom                     # E_w[Z]
    Ezz <- crossprod(Zi, Zi * (wr / denom))            # E_w[ZZ^T]
    VarZ <- Ezz - tcrossprod(mu)
    
    # Weight this event's contribution by the event weight
    G <- G + w[i] * VarZ
  }
  
  # Numerical guard: symmetrize & ridge if needed
  G <- (G + t(G)) / 2
  if (any(!is.finite(G))) {
    if (verbose) warning("Non-finite entries in D_Gradient; applying ridge.")
    G[!is.finite(G)] <- 0
    G <- G + diag(ridge, q)
  }
  G
}

# ---- Same structure for a curvature proxy; PSD and symmetric ----
D_HessianMatrix <- function(xGroup, zGroup, LGroup = NULL, yGroup, b, N = NULL,
                            q = ncol(zGroup[[1]]), wGroup = NULL,
                            verbose = FALSE, ridge = 1e-10) {
  # Use the same event-wise Var_w(z) accumulation as a stable curvature proxy.
  H <- D_Gradient(xGroup = xGroup, zGroup = zGroup, LGroup = LGroup,
                  yGroup = yGroup, b = b, N = N, wGroup = wGroup,
                  verbose = verbose, ridge = ridge)
  
  # Ensure symmetry/PSD
  H <- (H + t(H)) / 2
  eig <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values,
                  error = function(e) NA_real_)
  if (any(is.na(eig)) || min(eig) < ridge)
    H <- H + diag(ridge - min(0, min(eig, na.rm = TRUE)) + ridge, nrow(H))
  
  H
}






















armijoRule_L <- function(xGroup, yGroup, zGroup, L, l, k, grad, hessian,
                         b, zIdGrp, linNonpen, lambda, nonpen, penalty,
                         ll1, converged, control, ntot, fctOld) {
  
  L.new <- L
  Llk <- L[l, k]
  
  grad <- grad / ntot
  hessian <- hessian / ntot
  
  # Compute step size dk
  if (linNonpen) {
    dk <- -grad / hessian
  } else {
    L2norm <- sqrt(sum(L[l, ]^2))
    if (L2norm == 0 || is.na(L2norm)) {
      warning("L2norm is zero or NA, skipping this update")
      return(list(L = L, fct = fctOld, converged = converged + 1))
    }
    
    if (penalty == "lasso") {
      dk <- (-grad - lambda / L2norm * Llk) / (hessian + lambda / L2norm)
    } else if (penalty == "scad") {
      group_scad <- scad_group(L2norm, lambda)
      if (group_scad == 0 || is.na(group_scad)) {
        warning("group_scad is zero or NA, skipping this update")
        return(list(L = L, fct = fctOld, converged = converged + 1))
      }
      dk <- (-grad - lambda / L2norm * Llk) / (hessian + lambda / group_scad)
    }
  }
  
  if (is.na(dk) || !is.finite(dk)) {
    warning("dk is NA or not finite — skipping Armijo update")
    return(list(L = L, fct = fctOld, converged = converged + 1))
  }
  
  # Evaluate function
  if (dk != 0) {
    # Compute delta_k
    if (!is.na(linNonpen) && linNonpen == TRUE) {
      deltak <- dk * grad + control$gamma * dk^2 * hessian
    } else {
      L.tmp <- L
      L.tmp[l, k] <- Llk + dk
      deltak <- dk * grad + control$gamma * dk^2 * hessian +
        lambda * (sqrt(sum(L.tmp[l, ]^2)) - sqrt(sum(L[l, ]^2)))
    }
    
    for (j in 0:control$maxArmijo) {
      L.new[l, k] <- Llk + control$a_init * control$delta^j * dk
      
      fctNew <- objFunction_L(
        xGroup = xGroup,
        yGroup = yGroup,
        zGroup = zGroup,
        zIdGrp = zIdGrp,
        b = b,
        L = L.new,
        lambda = lambda,
        nonpen = nonpen,
        penalty = penalty,
        ll1 = ll1,
        ntot = ntot
      )
      
      addDelta <- control$a_init * control$delta^j * control$rho * deltak * ntot
      
      # Skip update if invalid
      if (any(!is.finite(c(fctNew, fctOld, addDelta)))) {
        warning("Skipping update due to NA/NaN/Inf in Armijo condition.")
        return(list(L = L, fct = fctOld, converged = converged + 1))
      }
      
      if (fctNew <= fctOld + addDelta) {
        L[l, k] <- Llk + control$a_init * control$delta^j * dk
        fct <- fctNew
        break
      }
      
      if (j == control$maxArmijo) {
        converged <- converged + 2
        fct <- fctOld
      }
    }
  } else {
    fct <- objFunction_L(
      xGroup = xGroup,
      yGroup = yGroup,
      zGroup = zGroup,
      zIdGrp = zIdGrp,
      b = b,
      L = L,
      lambda = lambda,
      nonpen = nonpen,
      penalty = penalty,
      ll1 = ll1,
      ntot = ntot
    )
  }
  
  return(list(L = L, fct = fct, converged = converged))
}

  



objFunction_L <- function(xGroup, yGroup, zGroup, zIdGrp, b, L, lambda, nonpen, penalty, ll1 = NULL, ntot = NULL) {
  
  # Construct D = LL^T and initialize
  D <- L %*% t(L)
  
  # Compute marginal (Laplace-approximated) log-likelihood
  laplace_terms <- mapply(function(x, y, z) {
    
    # Initial frailty u = 0
    u0 <- rep(0, ncol(z))
    
    # Gradient and Hessian of frailty at u = 0
    gh <- frailty_gradient_hessian(x, y, z, b, u = u0)
    
    # Laplace approx = -eta'y + log sum(exp(eta)) + 1/2 log|H|
    loglik <- -sum(y[, 2] * gh$eta) + sum(log(cumsum(rev(exp(rev(gh$eta))))))
    
    # Add Gaussian prior penalty (u ~ N(0, D))
    quad_pen <- sum((u0)^2 / diag(D))
    
    # Hessian contribution
    logdetH <- determinant(gh$hess, logarithm = TRUE)$modulus
    
    return(loglik + 0.5 * quad_pen + 0.5 * logdetH)
    
  }, xGroup, yGroup, zGroup, SIMPLIFY = TRUE)
  
  laplace_total <- sum(unlist(laplace_terms))

  
  # Penalty on L (group L2 or SCAD)
  L2norm <- sum(apply(L[-nonpen, , drop = FALSE], 1, function(row) sqrt(sum(row^2))))
  
  pen.term.L <- if (penalty == "lasso") {
    lambda * L2norm
  } else if (penalty == "scad") {
    sum(scad_group(L2norm, lambda))
  } else {
    0
  }
  
  # Final penalized objective
  obj <- laplace_total + pen.term.L
  return(obj)
}



densityfunc <- function(xGrp, yGrp, V, b,ll1.tmp){
  
  detv <- det(V)
  #if(detv==0){
  #  detv <- 0.0000001
  #}
  logv <- log(detv)
  
  ri <- yGrp-xGrp%*%b
  Vinv <- solve(V)
  rtvr <- t(ri)%*%Vinv%*%ri
  
  -0.5*(logv+rtvr)-ll1.tmp*length(yGrp)
}


em_Q <- function(betaList, LList, membership, pi_vec, lam1_base, lam2_base) {
  # membership: G x N; pi_vec: length G
  # Effective lambdas scaled by cluster mass (your current scheme)
  lam1_eff <- lam1_base * pi_vec
  lam2_eff <- lam2_base * pi_vec
  
  # Parameter part: sum_g [ penalized (negative) Laplace objective with weights w_g ]
  param_term <- sum(mapply(function(b, L, w, l1, l2) {
    pen_obj(b, L, as.numeric(w), l1, l2)
  }, betaList, LList, as.data.frame(t(membership)), lam1_eff, lam2_eff))
  
  # Mixing proportions part: - sum_{i,g} w_{ig} log pi_g
  # (Add a small floor for numerical safety)
  pi_safe <- pmax(pi_vec, 1e-12)
  mix_term <- - sum(membership * log(pi_safe))  # membership is G x N
  
  param_term + mix_term
}



obj_with_mix <- function(betaList, LList, membership, pi_vec, lam1_eff, lam2_eff) {
  # param part: sum_g [ - Laplace loglik_g (given w_g) + penalties_g ]
  param_term <- sum(mapply(function(b, L, w, l1, l2) {
    pen_obj(b, L, as.numeric(w), l1, l2)  # your pen_obj already = -ℓ_Laplace + penalties
  }, betaList, LList, as.data.frame(t(membership)), lam1_eff, lam2_eff))
  # mixing proportions part: - sum_{i,g} w_{ig} log pi_g
  pi_safe <- pmax(pi_vec, 1e-12)
  mix_term <- - sum(membership * log(pi_safe))  # membership is G x N
  param_term + mix_term
}


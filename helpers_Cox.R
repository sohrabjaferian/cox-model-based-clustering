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








# ZIdentity <- function(Z)
# {
#   ZId <- diag(dim(Z)[[1]])
#   return(list(ZId))
# }


ZIdentity <- function(Z) {
  diag(nrow(Z))  # not list()
}


covStartingValues <- function(xGroup,yGroup,zGroup,zIdGroup,b,N,lower=-10,upper=10)
{
  optimize1 <- function(x,y,b) y-x%*%b
  
  optimize2 <- function(gamma,zId,ZtZ)
  {
    H <- zId + exp(2*gamma)*ZtZ
    return(list(H=H))
  }
  
  optimize3 <- function(res,zId,ZtZ,gamma)
  {
    lambda <- optimize2(gamma,zId,ZtZ)
    logdetH <- determinant(lambda$H)$modulus
    quadH <- quad.form.inv(lambda$H,res)
    return(c(logdetH,quadH))
  }
  
  optimize4 <- function(gamma) {
    optH <- mapply(optimize3, resGroup, zIdGroup, ZtZ = ztzGroup, MoreArgs = list(gamma = gamma))
    H1 <- optH[1, ]
    H2 <- optH[2, ]
    
    if (any(!is.finite(H2)) || sum(H2) <= 0) {
      return(1e10)  # Large penalty instead of -Inf
    }
    
    fn <- N * log(sum(H2)) + sum(H1)
    return(fn)
  }
  
  optimize5 <- function(z) tcrossprod(z)
  
  resGroup <- mapply(optimize1,x=xGroup,y=yGroup,MoreArgs=list(b=b),SIMPLIFY=FALSE)
  ztzGroup <- mapply(optimize5,z=zGroup,SIMPLIFY=FALSE)
  
  optRes <- optimize(f=optimize4,interval=c(lower,upper))
  
  gamma <- optRes$minimum
  
  quadH <- mapply(optimize3,resGroup,zIdGroup,ztzGroup,MoreArgs=list(gamma=gamma))[2,]
  
  sig <- sqrt(1/N*sum(quadH))
  tau <- exp(gamma)*sig
  objfct <- 1/2*(optRes$objective + N*(1-log(N)))
  
  return(list(tau=tau,sigma=sig,opt=objfct))
}


nlogdet_Cox <- function(V_list) {
  sum(sapply(V_list, function(V) {
    -0.5 * determinant(V, logarithm = TRUE)$modulus[1]
  }))
}


VInv <- function(x, y, z, beta, D) {
  # Safety: Ensure all inputs are matrices
  if (!is.matrix(x) || !is.matrix(y) || !is.matrix(z)) {
    stop("Inputs x, y, and z must all be matrices.")
  }
  
  if (nrow(x) != nrow(y) || nrow(x) != nrow(z)) {
    stop("x, y, and z must have the same number of rows.")
  }
  
  # Compute linear predictor and risk
  
  
  if (!is.matrix(x)) stop("x is not a matrix")
  if (!is.numeric(x)) stop("x is not numeric")
  if (!is.numeric(beta)) stop("beta is not numeric")
  if (!is.vector(beta)) stop("beta is not a vector")
  

  
  eta <- x %*% beta
  risk <- exp(eta)
  
  # Order data by increasing survival time
  ord <- order(y[, 1])
  risk <- risk[ord]
  z_ord <- z[ord, , drop = FALSE]
  
  # Compute risk weights
  risk_cumsum <- rev(cumsum(rev(risk)))
  risk_weights <- risk / risk_cumsum
  risk_weights[!is.finite(risk_weights)] <- 0
  
  # Ensure risk_weights length matches z_ord rows
  if (length(risk_weights) != nrow(z_ord)) {
    warning("Mismatch between risk weights and z_ord. Returning identity matrix.")
    return(diag(ncol(z)))
  }
  
  # Construct diagonal weight matrix
  W <- diag(risk_weights)
  
  # Compute Hessian-like term with error handling
  H_u <- tryCatch({
    t(z_ord) %*% W %*% z_ord + solve(D)
  }, error = function(e) {
    warning("Matrix multiplication or inversion failed in VInv: ", conditionMessage(e))
    diag(ncol(z))  # Fallback: return identity
  })
  
  return(H_u)
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




frailty_gradient_hessian <- function(x, y, z, beta, u) {
  eta <- x %*% beta + z %*% u
  risk <- exp(eta)
  time <- y[, 1]
  event <- y[, 2]
  
  n <- nrow(x)
  grad <- matrix(0, ncol = ncol(z), nrow = 1)
  hess <- matrix(0, ncol = ncol(z), nrow = ncol(z))
  
  for (i in which(event == 1)) {
    # Risk set: those at risk at time[i]
    Ri <- which(time >= time[i])
    w <- risk[Ri]
    
    # Weighted means over risk set
    zw <- sweep(z[Ri,,drop=FALSE], 1, w, `*`)
    Ewz <- colSums(zw) / sum(w)
    
    grad <- grad + (z[i,,drop=FALSE] - Ewz)
    
    # Hessian approximation
    outer_Ewz <- tcrossprod(Ewz)
    Ezz <- crossprod(sqrt(w) * z[Ri,,drop=FALSE]) / sum(w)
    hess <- hess - (Ezz - outer_Ewz)
  }
  
  return(list(grad = grad, hess = hess, eta = eta))
}


nlogdet <- function(LGroup)
{
  nlogdetfun <- function(L)
  {
    -1/2*determinant(L)$modulus[1]
  }
  
  sum(mapply(nlogdetfun,LGroup))
}



cox_laplace_loglik <- function(xGroup, yGroup, zGroup, beta, D, tol = 1e-6, maxiter = 25) {
  if (length(zGroup) == 0 || is.null(zGroup[[1]]) || !is.matrix(zGroup[[1]]) || ncol(zGroup[[1]]) == 0) {
    warning("zGroup is malformed or empty in cox_laplace_loglik.")
    return(list(loglik = -1e6, uhat = NA, Hessian = diag(1e-2, ncol(D))))
  }
  
  q <- ncol(zGroup[[1]])
  u <- rep(0, q)
  
  for (iter in 1:maxiter) {
    grad_sum <- matrix(0, nrow = 1, ncol = q)
    hess_sum <- matrix(0, nrow = q, ncol = q)
    
    for (i in seq_along(xGroup)) {
      out <- frailty_gradient_hessian(xGroup[[i]], yGroup[[i]], zGroup[[i]], beta, u)
      grad_sum <- grad_sum + out$grad
      hess_sum <- hess_sum + out$hess
    }
    
    D_inv <- tryCatch({
      solve(D)
    }, error = function(e) {
      warning("D is near-singular; applying ridge regularization")
      solve(D + diag(1e-6, nrow(D)))
    })
    
    penalized_grad <- grad_sum - t(u) %*% D_inv
    penalized_hess <- hess_sum - D_inv
    
    if (anyNA(penalized_hess) || any(!is.finite(penalized_hess))) {
      warning("Non-finite penalized Hessian — skipping Laplace update")
      return(list(loglik = -1e6, uhat = rep(0, q), Hessian = diag(1e-2, q)))
    }
    
    if (anyNA(penalized_grad) || any(!is.finite(penalized_grad))) {
      stop("Non-finite penalized gradient")
    }
    
    penalized_hess_stable <- penalized_hess + diag(1e-6, q)
    
    step <- tryCatch({
      solve(penalized_hess_stable, t(penalized_grad))
    }, error = function(e) {
      warning("solve() failed in Laplace update; returning NA step")
      return(rep(NA, ncol(penalized_grad)))
    })
    
    if (anyNA(step)) {
      warning("Skipping update due to invalid step")
      break
    }
    
    u_new <- u - as.vector(step)
    
    if (max(abs(u_new - u)) < tol) {
      u <- u_new
      break
    }
    
    u <- u_new
  }
  
  # Final Laplace correction
  H_u <- hess_sum + D_inv
  
  logdetHu <- tryCatch({
    determinant(H_u, logarithm = TRUE)$modulus[1]
  }, error = function(e) {
    warning("Determinant failed; using fallback")
    1e-6
  })
  
  quad_penalty <- t(u) %*% D_inv %*% u
  
  loglik <- 0
  for (i in seq_along(xGroup)) {
    eta <- xGroup[[i]] %*% beta + zGroup[[i]] %*% u
    loglik <- loglik + cox_partial_loglik(yGroup[[i]], eta)
  }
  
  laplace_approx <- loglik - 0.5 * quad_penalty - 0.5 * logdetHu
  
  return(list(loglik = laplace_approx, uhat = u, Hessian = H_u))
}


  


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


cox_partial_score_component <- function(x, y, beta, j) {
  eta <- x %*% beta
  event <- y[, 2]
  risk <- exp(eta)
  risk_set_sum <- rev(cumsum(rev(risk)))
  weight <- risk / risk_set_sum
  
  weighted_means <- colSums(sweep(x, 1, weight, `*`))
  score <- sum(event * (x[, j] - weighted_means[j]))
  return(score)
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

D_Gradient <- function(xGroup, zGroup, LGroup, yGroup, b, N, verbose = FALSE) {
  q <- ncol(zGroup[[1]])
  mat <- matrix(0, nrow = q, ncol = q)
  
  for (i in seq_along(yGroup)) {
    x_i <- xGroup[[i]]
    z_i <- zGroup[[i]]
    y_i <- yGroup[[i]]
    
    # Skip empty groups
    if (nrow(z_i) == 0 || nrow(x_i) == 0) next
    
    # Linear predictor and risk
    eta_i <- x_i %*% b
    risk <- as.numeric(exp(eta_i))
    
    if (!all(is.finite(risk)) || length(risk) == 0) {
      if (verbose) warning(paste("Skipping group", i, "- invalid risk"))
      next
    }
    
    # Cumulative risk with stability adjustment
    risk_sum <- rev(cumsum(rev(risk)))
    risk_sum[risk_sum == 0] <- 1e-10
    w <- risk / risk_sum
    
    # Final check for weight vector
    if (!all(is.finite(w)) || length(w) == 0 || any(is.na(w))) {
      if (verbose) warning(paste("Skipping group", i, "- invalid weights"))
      next
    }
    
    # Weighted Z^T W Z
    W <- diag(w, nrow = length(w), ncol = length(w))
    ztwz <- t(z_i) %*% W %*% z_i
    mat <- mat + ztwz
  }
  
  return(mat)
}


  



D_HessianMatrix <- function(xGroup, zGroup, LGroup, yGroup, b, N, q, verbose = FALSE) {
  hessian <- matrix(0, nrow = q, ncol = q)
  
  for (i in seq_along(yGroup)) {
    x_i <- xGroup[[i]]
    z_i <- zGroup[[i]]
    y_i <- yGroup[[i]]
    
    # Skip if any input is empty
    if (nrow(z_i) == 0 || nrow(x_i) == 0 || length(y_i) == 0) next
    
    time <- y_i[, 1]
    event <- y_i[, 2]
    
    eta_i <- x_i %*% b
    risk <- as.numeric(exp(eta_i))
    
    if (!all(is.finite(risk)) || length(risk) == 0) {
      if (verbose) warning(paste("Skipping group", i, "- invalid risk"))
      next
    }
    
    risk_sum <- rev(cumsum(rev(risk)))
    risk_sum[risk_sum == 0] <- 1e-10
    w <- risk / risk_sum
    
    # Validate w before forming diag
    if (anyNA(w) || any(!is.finite(w)) || length(w) != length(risk)) {
      if (verbose) warning(paste("Skipping group", i, "- invalid weights"))
      next
    }
    
    # Safe diagonal matrix creation
    W <- diag(w, nrow = length(w), ncol = length(w))
    
    ztwz <- t(z_i) %*% W %*% z_i
    hessian <- hessian + ztwz
  }
  
  return(hessian)
}













matsplitter <- function(M, q) {
  n_blocks <- q * q
  block_list <- vector("list", n_blocks)
  
  for (i in 1:q) {
    for (j in 1:q) {
      idx <- (i - 1) * q + j
      row_idx <- ((i - 1) * q + 1):(i * q)
      col_idx <- ((j - 1) * q + 1):(j * q)
      block_list[[idx]] <- M[row_idx, col_idx, drop = FALSE]
    }
  }
  return(block_list)
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



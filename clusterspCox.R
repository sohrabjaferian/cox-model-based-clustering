clustercox <- function(x,z,grp,time,event,lam1,lam2,nCluster,nonpen.b=1,nonpen.L=1,penalty.b=c("lasso","scad"),
                         penalty.L=c("lasso","scad"),
                       standardize = TRUE,
                       control=spCoxControl()){

  

  set.seed(123)
  
  penalty.b <- match.arg(penalty.b)
  penalty.L <- match.arg(penalty.L)
  
  # transform data
  if(is.data.frame(x)) x = as.matrix(x)
  if(is.data.frame(z)) z = as.matrix(z)

  
  # some checks
  if (!is.matrix(x)) stop("x has to be a matrix or data frame")
  if (any(is.na(x))) stop("Missing values in x not allowed")
  
  if (!is.numeric(time) || !is.numeric(event)) stop("time and event must be numeric")
  if (length(time) != length(event)) stop("time and event must have equal length")
  
  y <- Surv(time, event)
  

  if (any(is.na(y))) stop("Missing values in y not allowed")
  
  
  if (nrow(x) != nrow(y)) stop("x and y have not correct dimensions")
  
  
  if (!is.matrix(z)) stop("z has to be a matrix")
  if (any(is.na(z))) stop("Missing values in z not allowed")
  
  if (any(x[,1]!=rep(1,dim(x)[[1]]))) stop("first column is not the intercept")
  
  if (!all(nonpen.b%in%1:dim(x)[[2]])) stop("Error with the argument nonpen for beta")
  if (!all(nonpen.L%in%1:dim(z)[[2]])) stop("Error with the argument nonpen for L")
  
  if (length(which(lam1<0))>0) stop("lam1 must be positive")
  if (length(which(lam2<0))>0) stop("lam2 must be positive")
  
  if (length(lam1)==1) lam1 <- rep(lam1, nCluster)
  if (length(lam2)==1) lam2 <- rep(lam2, nCluster)
  
  
  # grp <- factor(grp)  # Make it a factor
  # grp <- droplevels(grp)
  # 
  # timeGrp  <- split(time, grp)
  # eventGrp <- split(event, grp)
  
  
  ##### Standardize covariates
  
  print("typeof(standardize):")
  print(typeof(standardize))
  print("standardize value:")
  print(standardize)
  
  # standardize <- isTRUE(standardize)
  
  if (standardize)
  {
    xOr <- x
    meanx <- apply(x[,-1,drop = FALSE],2,mean)
    sdx <- apply(x[,-1,drop = FALSE],2,sd)
    x <- cbind(1,scale(x[,-1],center=meanx,scale=sdx))

    zOr <- z
    meanz <- apply(z[,-1,drop = FALSE],2,mean)
    sdz <- apply(z[,-1,drop = FALSE],2,sd)
    z <- cbind(1,scale(z[,-1,drop = FALSE],center=meanz,scale=sdz))
  }

  ##### allocate variables
  
  grp <- factor(grp)
  N <- length(unique(grp))  # or length(levels(grp)) if grp is a factor
  p <- dim(x)[[2]]         # p is the number of covariates
  q <- dim(z)[[2]]         # q is the number of random effects variables
  ntot <- length(time)
  Q <- q*(q+1)/2           # maximum number of variance components parameters
  lambda1 <- lam1*(N/nCluster)
  lambda2 <- lam2*(N/nCluster)

  
  ###### initiation

  x_kmeans <- x[, -1, drop = FALSE]
  # ini.clusters <- kmeans(x_kmeans, centers = nCluster, nstart = 10)$cluster
  
  ini.fit <- Mclust(as.matrix(x_kmeans), G = nCluster)
  
  yGrp.mixture <- xGrp.mixture <- zGrp.mixture <- ziGrp.mixture <- zIdGrp.mixture <- list()
  memb.prob <- vector()
  betaStart <- covStart <- parsStart <- LStart <- DStart <- VInvGrp <- list()
  membership <- matrix(ncol = N, nrow = nCluster)
  

  
  xGrp <- lapply(seq_len(nrow(x)), function(i) matrix(x[i, ], nrow = 1))
  zGrp <- lapply(seq_len(nrow(z)), function(i) matrix(z[i, ], nrow = 1))
  yGrp <- lapply(seq_len(nrow(x)), function(i) {
    mat <- matrix(c(time[i], event[i]), nrow = 1)
    colnames(mat) <- c("time", "event")
    mat
  })
  zIdGrp <- lapply(zGrp, ZIdentity)
  

  
  ##ll1 <- 1/2*ntot*log(2*pi)
  ll1 <- list()
  ntot.mix <- rep(1,ntot)
  ntot.Grp <- split(ntot.mix,grp)
  ntot.mixture <-list()
  N.mixture <- list()
  
  fctStart <- vector()
  
  
  for (i in 1:nCluster) {
    
    # Membership from soft clustering (initial classification is hard)
    zi <- as.integer(ini.fit$classification == i)
    memb.prob[i] <- sum(zi) / N
    membership[i, ] <- zi
    N.mixture[[i]] <- sum(zi)
    ll1[[i]] <- NA
    
    # ----- Soft assignment across all N subjects -----
    
    # Scale survival times only (not events)
    safe_multiply_surv <- function(weight, surv_obj) {
      Surv(time = pmax(weight * surv_obj[, 1], 1e-4), event = surv_obj[, 2])
    }
    
    yGrp.mixture[[i]]   <- mapply(safe_multiply_surv, membership[i, ], yGrp, SIMPLIFY = FALSE)
    xGrp.mixture[[i]]   <- mapply(multiplication, membership[i, ], xGrp, SIMPLIFY = FALSE)
    zGrp.mixture[[i]]   <- mapply(multiplication, membership[i, ], zGrp, SIMPLIFY = FALSE)
    zIdGrp.mixture[[i]] <- mapply(multiplication, membership[i, ], zIdGrp, SIMPLIFY = FALSE)
    
    y_i <- yGrp.mixture[[i]]
    x_i <- xGrp.mixture[[i]]
    z_i <- zGrp.mixture[[i]]
    
    
    # ----- Build full stacked matrices -----
    
    xi <- do.call(rbind, xGrp.mixture[[i]])
    times_i  <- unlist(lapply(yGrp.mixture[[i]], function(x) x[, 1]))
    status_i <- unlist(lapply(yGrp.mixture[[i]], function(x) x[, 2]))
    yi <- Surv(times_i, status_i)

    
    
    cat("Cluster", i, ":\n")
    cat("xi dim:", dim(xi), "\n")
    cat("Any NA:", any(is.na(xi)), "\n")
    cat("Any Inf:", any(is.infinite(xi)), "\n")
    cat("Col SDs:\n")
    print(apply(xi, 2, sd))
    
    cat("First few survival times:\n")
    print(head(yi))
    
    
    coxfit <- glmnet(x = xi[,-1], y = yi, family = "cox", lambda = lambda1[i],
                     alpha = if (penalty.b == "lasso") 1 else 0)
    betaStart[[i]] <- c(0, as.numeric(coxfit$beta[, 1]))  # prepend intercept as 0
    
    # Frailty variance initialization using Cox-specific Laplace function
    y_time_only <- lapply(y_i, function(y) y[, 1])  # extract time
    
    N_i <- N.mixture[[i]]
    
    covStart[[i]] <- covStartingValues(
      xGroup = x_i,
      yGroup = y_time_only,
      zGroup = z_i,
      zIdGroup = zIdGrp.mixture[[i]], 
      b = betaStart[[i]],
      N = N_i
    )
    
    
    tau <- covStart[[i]]$tau
    
    DStart[[i]] <- diag(tau, q)
    LStart[[i]] <- chol(DStart[[i]])
    VInvGrp[[i]] <- mapply(VInv,
                           x = x_i,
                           y = y_i,
                           z = z_i,
                           MoreArgs = list(beta = betaStart[[i]], D = DStart[[i]]),
                           SIMPLIFY = FALSE)
    
    
    # --- Calculate objective function for the starting values ---
    # ------------------------------------------------------------
    
    fctStart[i] <- ObjFunction(
      xGroup = x_i,
      yGroup = y_i,
      zGroup = z_i,
      beta = betaStart[[i]],
      L = LStart[[i]],
      lambda1 = lambda1[i],
      lambda2 = lambda2[i],
      nonpen.b = nonpen.b,
      nonpen.L = nonpen.L,
      penalty_b = penalty.b,
      penalty_L = penalty.L
    )
  }
  
  
  
  # some necessary allocations:
  betaIter <- betaStart
  LIter <- LStart
  DIter <- DStart
  
  LvecIter <- lapply(LStart, function(x) x[lower.tri(x,diag = TRUE)])
  convPar <- max(unlist(lapply(betaIter, function(x) crossprod(x))))
  convCov <- max(sapply(LvecIter, function(x) crossprod(x)))
  
  
  
  fctIter <- convFct <- fctStart
  hessian0 <- rep(0,p)
  mat0 <- matrix(0,ncol=p,nrow=N)
  covIter <- LvecIter
  
  ##### algorithm parameters
  
  stopped <- FALSE
  doAll <- FALSE
  converged <- 0
  counterIn <- 0
  counter <- 0       # counts the number of outer iterations
  
  convFct2 <- -10
  
  #while (counter<2) {
  
  
  while((counter<control$maxIter)&(convFct2<0|counter<1)&((convPar>control$tol|convFct[1]>control$tol|convCov>control$tol|!doAll ))) {
    #while((counter<control$maxIter)&((convPar>control$tol|convFct>control$tol|convCov>control$tol|!doAll ))) {
    counter <- counter + 1 
    
    
    betaIterOld <- betaIter
    LIterOld <- LIter
    fctIterOld <- fctIter
    covIterOld <- covIter
    
    
    
    activeSet <- lapply(betaIter, function(x) which(x!=0))
    #activeSet <- lapply(activeSet, function(x) 1:p)
    #if ((length(activeSet)>min(p,ntot))&(lambda1>0)&(counter>2)) {stopped <- TRUE ; break}
    
    if (counterIn==0 | counterIn>control$number)
    {
      doAll <- TRUE
      activeSet <- lapply(activeSet, function(x) 1:p)
      counterIn <- 1    
    } else
    {
      doAll <- FALSE
      counterIn <- counterIn+1
    }
    
    
    for (i in 1:nCluster) {
      # --- optimization w.r.t the fixed effects vector beta ---
      # --------------------------------------------------------
      
      x <- do.call(rbind, xGrp.mixture[[i]])
      z <- do.call(rbind, zGrp.mixture[[i]])
      y <- do.call(rbind, yGrp.mixture[[i]])
      
      if (is.null(x) || is.null(z) || is.null(y)) {
        stop("One of x, y, or z is NULL — mixture component ", i)
      }
      
      x <- as.matrix(x)
      z <- as.matrix(z)
      
      HessOut <- frailty_gradient_hessian(
        x,
        y,
        z,
        beta = betaIter[[i]],
        u = rep(0, ncol(zGrp.mixture[[i]][[1]]))
      )
      HessIter <- HessIterTrunc <- HessOut$hess
      
      HessIter[activeSet[[i]]] <- pmin(pmax(HessIter[activeSet[[i]]],control$lower),control$upper)
      # LxGrp <- as1(xGrp.mixture[[i]],VInvGrp[[i]],activeSet[[i]],N=N)
      ll2 <- nlogdet(LGroup=VInvGrp[[i]])
      
      for (j in activeSet[[i]])
      {
        cut1 <- cox_partial_score_component(
          x = do.call(rbind, xGrp.mixture[[i]]),
          y = do.call(rbind, yGrp.mixture[[i]]),
          beta = betaIter[[i]],
          j = j
        )
        
        JinNonpen <- j%in%nonpen.b
        
        # optimum can be calculated analytically
        if (!is.na(HessIterTrunc[j]) && !is.na(HessIter[j]) && HessIterTrunc[j]==HessIter[j])
        {
          if (JinNonpen) {betaIter[[i]][j] <- cut1/HessIter[j]} else {
            if(penalty.b=="scad"){
              scada = 3.7
              betaIter[[i]][j] <- SoftThreshold(cut1,lambda1[i])/(HessIter[j]*(1-1/scada))
              
            }else if(penalty.b=="lasso"){
              betaIter[[i]][j] <- SoftThreshold(cut1,lambda1[i])/HessIter[j]
            }
            
          }
        }else
          
          # optimimum is determined by the armijo rule
        {
          armijo <- ArmijoRule_b(
            x = do.call(rbind, xGrp.mixture[[i]]),
            y = do.call(rbind, yGrp.mixture[[i]]),
            z = do.call(rbind, zGrp.mixture[[i]]),
            b = betaIter[[i]],
            j = j,
            cut = cut1,
            HkOldJ = HessIterTrunc[j],
            HkJ = HessIter[j],
            JinNonpen = JinNonpen,
            lambda = lambda1[i],
            nonpen = nonpen.b,
            penalty = penalty.b,
            L = LIter[[i]],
            converged = converged,
            control = control
          )
          
          
          
          
          betaIter[[i]] <- armijo$b
          converged <- armijo$converged
          fctIter[i] <- armijo$fct
        }
        
      }
      betaIter[[i]][abs(betaIter[[i]])<0.05]=0
      
      # --- optimization w.r.t the variance components parameters ---
      # -------------------------------------------------------------
      
      # calculations before the covariance optimization
      activeSet[[i]] <- which(betaIter[[i]]!=0)
      ll4 <- lambda1[i]*sum(abs(betaIter[[i]][-nonpen.b]))
      
      # optimization of L
      
      activeSet.L = which(rowSums(abs(LIter[[i]]))!=0)
      
      
      # calculate the hessian matrices for k in the activeSet
      
      D.grad = D_Gradient(xGroup=xGrp.mixture[[i]],zGroup=zGrp.mixture[[i]],LGroup=VInvGrp[[i]],yGroup=yGrp.mixture[[i]],b=betaIter[[i]],N=N)
      L.grad = t(LIter[[i]]%*%(D.grad+t(D.grad)))
      
      D.hessian = D_HessianMatrix(xGroup=xGrp.mixture[[i]],zGroup=zGrp.mixture[[i]],LGroup=VInvGrp[[i]],yGroup=yGrp.mixture[[i]],b=betaIter[[i]],N=N,q=q)
      # D.hessian.submatrix = matsplitter(D.hessian,q)
      D.hessian.submatrix <- replicate(q^2, D.hessian, simplify = FALSE)
      
      
      for (k in 1:q) {
        
        L.hessian.sub <- sapply(D.hessian.submatrix[((k - 1) * q + 1):(k * q)], function(x) {
          if (is.matrix(x) && ncol(x) >= k) {
            return(x[, k])
          } else {
            warning("Invalid or undersized matrix in D.hessian.submatrix")
            return(rep(NA, nrow(x)))  # or rep(0, q) or NA depending on how you want to handle it
          }
        })
        
        
        submat_vector <- unlist(D.hessian.submatrix[[(k - 1) * q + k]])
        
        if (length(submat_vector) == q * q) {
          Dkk_mat <- matrix(submat_vector, ncol = q, byrow = TRUE)
        } else {
          warning(sprintf("D.hessian.submatrix[[%d]] is malformed: expected length %d, got %d", 
                          (k - 1) * q + k, q*q, length(submat_vector)))
          Dkk_mat <- matrix(0, ncol = q, nrow = q)  # safe default
        }
        
        L.hessian = diag(2*D.grad[k,k],q,q)+2*LIter[[i]]%*%(matrix(unlist(D.hessian.submatrix[[(k-1)*q+k]]),ncol = q,byrow = TRUE)+L.hessian.sub)%*%t(LIter[[i]])
        
        for (l in intersect(k:q,activeSet.L)) {
          
          L.lk.grad <- L.grad[l,k]
          L.lk.Hess <- L.hessian[l,l]
          L.lk.Hess <- min(max(L.lk.Hess,control$lower),control$upper)
          
          linNonpen <- l%in%nonpen.L
          
          ##armijo <- ArmijoRule_L(xGroup=xGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], yGroup=yGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], zGroup=zGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], L=LIter[[i]], l=l-1,k=k-1,grad=L.lk.grad,hessian=L.lk.Hess, 
          ##                       b=betaIter[[i]],sigma=sigmaIter[[i]],zIdGrp=zIdGrp.mixture[[i]], linNonpen=linNonpen, nonpen=nonpen.L-1, lambda=lambda2[i], penalty=penalty.L, 
          ##                       ll1=ll1[[i]], gamma=control$gamma, maxArmijo=control$maxArmijo, a_init=control$a_init, delta=control$delta, rho=control$rho, converged=converged, ntot=1)
          
          valid_idx <- which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x)) != 0)))
          
          
          fctOld <- objFunction_L(
            xGroup = xGrp.mixture[[i]][valid_idx],
            yGroup = yGrp.mixture[[i]][valid_idx],
            zGroup = zGrp.mixture[[i]][valid_idx],
            zIdGrp = zIdGrp[valid_idx],
            b = betaIter[[i]],
            L = LIter[[i]],
            lambda = lambda2[i],
            nonpen = nonpen.L,
            penalty = penalty.L,
            ll1 = ll1[[i]],
            ntot = N.mixture[[i]]
          )
          
          
          
          armijo <- armijoRule_L(xGroup=xGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], yGroup=yGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], zGroup=zGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], L=LIter[[i]], l=l,k=k,grad=L.lk.grad,hessian=L.lk.Hess, 
                                 b=betaIter[[i]], zIdGrp=zIdGrp[which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], linNonpen=linNonpen,nonpen=nonpen.L, lambda=lambda2[i], penalty=penalty.L, 
                                 ll1=ll1[[i]],converged=converged,control=control, ntot=N.mixture[[i]], fctOld = fctOld)
          
          LIter[[i]] <- armijo$L
          converged <- armijo$converged
          fctIter[i] <- armijo$fct
          
        }
        
        
      }
      
      LIter[[i]][abs(LIter[[i]]) < 1e-2] <- 0
      
      DIter[[i]] = LIter[[i]]%*%t(LIter[[i]])
      LvecIter[[i]] = LIter[[i]][lower.tri(LIter[[i]],diag = TRUE)]
      
  
      VInvGrp[[i]] <- mapply(VInv,
                             x = xGrp.mixture[[i]],
                             y = yGrp.mixture[[i]],
                             z = zGrp.mixture[[i]],
                             MoreArgs = list(beta = betaIter[[i]], D = DStart[[i]]),
                             SIMPLIFY = FALSE)


      
      covIter[[i]] <- LvecIter[[i]]
      
      
    }

    
    ##### clustering
    
    ##### update mixture membership
    
    ind.prob <- matrix(nrow = nCluster, ncol = N)
    
    
    for (i in 1:nCluster) {
      
      # Use negative Laplace-approximated log-likelihood per subject (unpenalized)
      if (length(xGrp.mixture[[i]]) == 0 || length(yGrp.mixture[[i]]) == 0) {
        warning(sprintf("Cluster %d is empty. Skipping update.", i))
        next
      }
    
      
      
      ll_per_subject <- mapply(function(xj, yj, zj) {
        cox_laplace_loglik(
          xGroup = list(xj),
          yGroup = list(yj),
          zGroup = list(zj),
          beta = betaIter[[i]],
          D = LIter[[i]] %*% t(LIter[[i]])
        )$loglik
      }, xGrp, yGrp, zGrp)
      
      
      
      # Convert to unnormalized densities
      cluster_den <- exp(ll_per_subject)
      ind.prob[i, ] <- cluster_den * memb.prob[i]
    }
    
    # Handle potential zero columns in ind.prob (all clusters 0 for some subjects)
    zero_cols <- which(colSums(ind.prob) == 0)
    if (length(zero_cols) > 0) {
      ind.prob[, zero_cols] <- 1 / nCluster  # Avoid division by zero
    }
    
    

    # Soft assignment: assign each subject to the cluster with max posterior
    ind.prob[, which(colSums(ind.prob) == 0)] <- 1 / nCluster
    # Enforce minimum membership probability
    eps <- 1e-4
    membership[membership < eps] <- eps
    membership <- membership / matrix(colSums(membership), nrow = nCluster, ncol = N, byrow = TRUE)
    
    
    
    # Update mixture weights
    memb.prob <- rowSums(membership) / N
    
    
    # Update per-cluster data splits using soft membership
    for (i in 1:nCluster) {
      yGrp.mixture[[i]]    <- mapply(multiplication, membership[i, ], yGrp, SIMPLIFY = FALSE)
      xGrp.mixture[[i]]    <- mapply(multiplication, membership[i, ], xGrp, SIMPLIFY = FALSE)
      zGrp.mixture[[i]]    <- mapply(multiplication, membership[i, ], zGrp, SIMPLIFY = FALSE)
      zIdGrp.mixture[[i]]  <- mapply(multiplication, membership[i, ], zIdGrp, SIMPLIFY = FALSE)
      
      N.mixture[[i]]       <- sum(membership[i, ])
      ntot.mixture[[i]]    <- sum(unlist(mapply(multiplication, membership[i, ], ntot.Grp, SIMPLIFY = FALSE)))
    }
    
    
    
      
      # Optional: update ll1 (not strictly needed for Cox)
      ll1[[i]] <- 1/2 * ntot.mixture[[i]] * log(2 * pi)
    }
    
    
    #lambda1 <- lam1*(rowSums(membership))
    #lambda2 <- lam2*(rowSums(membership))
    
    ##### reorder cluster
    
    
    
    
    # --- check convergence ---
    
    convPar <- max(mapply(function(x,y) sqrt(crossprod(x-y))/(1+sqrt(crossprod(x))), x=betaIter, y=betaIterOld))
    #convPar <- max(mapply(function(x,y) sqrt(crossprod(x-y)),x=betaIter, y=betaIterOld))
    convFct <- abs((sum(fctIterOld)-sum(fctIter))/(1+abs(sum(fctIter))))
    #convFct <- abs((sum(fctIterOld)-sum(fctIter)))
    convFct2 <- sum(fctIter) - sum(fctIterOld)
    #convCov <- max(mapply(function(x,y) sqrt(crossprod(x-y))/(1+sqrt(crossprod(x))), x=covIter, y=covIterOld))
    convCov <- max(mapply(function(x,y) sqrt(crossprod(x-y)), x=covIter, y=covIterOld))
    
    
    if (!any(is.na(c(convPar, convFct, convCov))) &&
        convPar <= control$tol &&
        convFct <= control$tol &&
        convCov <= control$tol) {
      counterIn <- 0
    }
    
    if (
      !any(is.na(c(convPar, convFct, convCov))) &&
      convPar <= control$tol &&
      convFct <= control$tol &&
      convCov <= control$tol
    ) {
      counterIn <- 0
    }
    

    
  if (standardize) {
    betaIter <- lapply(betaIter, function(x) {
      x[-1] <- x[-1]/sdx
      x[1] <- x[1] - sum(meanx * x[-1])
      return(x)
    })
    
    x <- xOr
    xGrp <- lapply(seq_len(nrow(x)), function(i) matrix(x[i, ], nrow = 1))
    
    z <- zOr
    zGrp <- lapply(seq_len(nrow(z)), function(i) matrix(z[i, ], nrow = 1))
  }
    
  
  ntot <- nrow(x)


  
  
  # --- summary information ---
  # ---------------------------
  ntot.mix <- rep(1,ntot)
  ntot.Grp <- split(ntot.mix,grp)
  ntot.mixture <-list()
  
  for (i in 1:nCluster) {

    xGrp.mixture[[i]] <- mapply(multiplication, membership[i,], xGrp, SIMPLIFY = FALSE)
    
    yGrp.mixture[[i]] <- mapply(multiplication, membership[i,], yGrp, SIMPLIFY = FALSE)
    
    zGrp.mixture[[i]] <- mapply(multiplication, membership[i,], zGrp, SIMPLIFY = FALSE)
    
    zIdGrp.mixture[[i]] <- mapply(multiplication, membership[i,], zIdGrp, SIMPLIFY = FALSE)
    
    
    
    
    ntot.mixture[[i]] <- mapply(multiplication,membership[i,],ntot.Grp, SIMPLIFY = FALSE)
    VInvGrp[[i]] <- mapply(
      VInv,
      x = xGrp.mixture[[i]],
      y = yGrp.mixture[[i]],
      z = zGrp.mixture[[i]],
      MoreArgs = list(beta = betaIter[[i]], D = DIter[[i]]),
      SIMPLIFY = FALSE
    )
  }
  
  # Sum of group membership weights (i.e., number of subjects per cluster)
  N.mixture <- apply(membership, 1, sum)
  
  # Just record it as vector (we no longer need nested structure)
  ntot.mixture <- as.numeric(N.mixture)
  
  # Optional: Ensure no zero-cluster (if you allow zero-weighted clusters)
  ntot.mixture[ntot.mixture < 1] <- 1
  
  
  N.mixture[N.mixture<1] <- 1
  D <- DIter
  #npar <- sum(unlist(lapply(betaIter, function(x) sum(x!=0)))) + sum(unlist(lapply(D, function(x) sum(diag(x)!=0)))) + 1
  npar <- sum(unlist(lapply(betaIter, function(x) sum(x != 0)))) + length(unlist(LvecIter))
  
  logLik <- sum(mapply(function(x, y, z, beta, L) {
    D <- L %*% t(L)
    cox_laplace_loglik(xGroup = list(x),
                       yGroup = list(y),
                       zGroup = list(z),
                       beta = beta,
                       D = D)$loglik
  },
  xGrp.mixture, yGrp.mixture, zGrp.mixture, betaIter, LIter))
  
  deviance <- -2*logLik
  aic <- -2* logLik + 2*npar
  bic <- -2* logLik + log(ntot)*npar
  
  p <- sum(unlist(lapply(betaIter, function(x) sum(x!=0))))
  q <- sum(unlist(lapply(D, function(x) sum(diag(x)!=0))))
  bbic <- -2*logLik + max(1,log(log(p+q)))*log(ntot)*npar
  ebic <- -2*logLik + (log(ntot)+2*log(p+q))*npar
  
  if (converged>0) cat("Algorithm does not properly converge.","\n")
  if (stopped) {cat("|activeSet|>=min(p,ntot): Increase lambda or set stopSat=FALSE.","\n")
    ;  LvecIter <- nlogLik <- aic <- bic <- NA ; betaIter <- rep(NA,p) ; bi <- fitted <- residuals <- NULL}
  
  out <- list(
    data = list(x = x, y = y, z = z, grp = grp),
    membership = membership,
    coefInit = list(betaStart = betaStart, parsStart = parsStart),
    penalty.b = penalty.b,
    penalty.L = penalty.L,
    nonpen.b = nonpen.b,
    nonpen.L = nonpen.L,
    lambda1 = lambda1,
    lambda2 = lambda2,
    Lvec = LvecIter,
    coefficients = betaIter,
    D = DIter,
    converged = converged,
    logLik = logLik,
    npar = npar,
    deviance = deviance,
    aic = aic,
    bic = bic,
    bbic = bbic,
    ebic = ebic,
    counter = counter,
    control = control,
    call = match.call(),
    stopped = stopped,
    objective = fctIter
  )
  structure(out, class = "spcox")  # Change class name to reflect Cox model

}

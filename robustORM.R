## Robust Ordinal Response Model via Divergence Approach

## INPUT
# y (factor): an ordered categorical response vector
# X (matrix): (n,p)-matrix of covariates without intercept term (n: total sample size; p: the number of covariates)
# method (character): Name of link function to be used
# divergence (character): Name of divergence to be used; "KL": Kullback-Liebler, "DP": Density-Power, "gamma": $\gamma$-divergence
# tnp (numeric): the value of tuning parameter for DP and $\gamma$-divergences; defalt is 0.5
# init_parameters (numeric): (p+M-1) vector of initial value of parameters; Input initial values for cutpoints for the first M-1 elements, and initial values for coefficients for the rest.

## OUTPUT
# (list): Estimates of coefficient parameters and cutoffs

# source
source("links.R")

robustORM = function(y, X, method=c("probit", "logistic", "cloglog", "loglog"), 
                     divergence=c("KL", "DP", "gamma"), tnp, init_parameters) {
  # Preparation -------------------------------------------------------
  method <- match.arg(method)
  divergence <- match.arg(divergence)
  if (missing(init_parameters)) init_parameters <- NULL
  if(missing(tnp)) tnp = 0.5
  if(tnp==0.0) divergence <- "KL"
  ## Parameters setting -------------------------------------------------------
  n = dim(X)[1] # sample size
  p = dim(X)[2] # the number of coefficient parameters
  M = max( as.numeric(y) ) # the number of categories
  q = p + M-1 # the number of parameters
  ## Setting the initial values -------------------------------------------------------
  if (is.null(init_parameters)) {
    ### Setting the start value for parameters "start <- clm.start(y.levels = rorm.struct$y.levels, threshold = threshold, X = rorm.struct$X, NOM = rorm.struct$NOM, has.intercept = TRUE)" ----------
    #### Setting the start value for cutpoints "st <- start.threshold(y.levels, threshold)" ----------
    init_delta <- as.vector( qlogis((1:(M-1))/((M-1) + 1)) )
    init_delta_tilde <- c(init_delta[1])
    for (m in 2:(M-1)) {
      init_delta_tilde[m] <- sqrt(init_delta[m] - init_delta[m-1])
    }
    #### Setting the start value for coefficients "start.beta(X, has.intercept = TRUE)" ----------
    init_beta <- rep(0, p)
    # init_parameters = c(init_beta, init_delta_tilde)
  }else if(all(init_parameters=="MLE")) {
    d = data.frame(y=as.ordered(y),X)
    fit = ordinal::clm(y~., data=d, link=method)
    init_beta = as.vector(fit$beta)
    init_delta = as.vector(fit$Theta)
    init_delta_tilde <- c(init_delta[1])
    for (m in 2:(M-1)) {
      init_delta_tilde[m] <- sqrt(init_delta[m] - init_delta[m-1])
    }
  }else{
    init_delta = init_parameters[1:(M-1)]
    init_delta_tilde <- c(init_delta[1])
    for (m in 2:(M-1)) {
      init_delta_tilde[m] <- sqrt(init_delta[m] - init_delta[m-1])
    }
    init_beta = init_parameters[-(1:(M-1))]
  }
  ## Setting initial values of optimaization "control <- do.call(clm.control, c(control, list(...)))" ----------------------------------------------
  maxIter = 100L
  gradTol = 1e-06
  relTol = 1e-06
  tol = sqrt(.Machine$double.eps)
  maxLineIter = 15L
  maxModIter = 5L
  ## Setting objective function -------------------------------------------------------
  obj_rorm = switch(divergence,
                    "KL" = function(beta,delta_tilde,y,X,link) {
                      n = dim(X)[1]
                      p = dim(X)[2]
                      M = max( as.numeric(y) )
                      delta = c(delta_tilde[1], delta_tilde[1]+cumsum(delta_tilde[-1]^2), Inf)
                      delta_0_M = c(-Inf, delta[-M])
                      # Delta = matrix( rep(delta,n),n,M,T )
                      # Delta_0_M = matrix( rep(delta_0_M,n),n,M,T )
                      Xb = X%*%beta
                      # XB = matrix( rep(X%*%beta,M),n,M )
                      -sum(log( drop(pfunc(delta[y]-Xb,method=link)-pfunc(delta_0_M[y]-Xb,method=link)) ))
                    },
                    "DP" = function(beta,delta_tilde,y,X,link,tnp) {
                      n = dim(X)[1]
                      p = dim(X)[2]
                      M = max( as.numeric(y) )
                      delta = c(delta_tilde[1], delta_tilde[1]+cumsum(delta_tilde[-1]^2), Inf)
                      delta_0_M = c(-Inf, delta[-M])
                      Delta = matrix( rep(delta,n),n,M,T )
                      Delta_0_M = matrix( rep(delta_0_M,n),n,M,T )
                      Xb = X%*%beta
                      XB = matrix( rep(X%*%beta,M),n,M )
                      -(1/tnp)*{(1/n)*sum( drop(pfunc(delta[y]-Xb,method=link)-pfunc(delta_0_M[y]-Xb,method=link))^tnp )} +
                        (1/(1+tnp))*{(1/n)*sum( (pfunc(Delta-XB,method=link)-pfunc(Delta_0_M-XB,method=link))^(1+tnp) )}
                    },
                    "gamma" = function(beta,delta_tilde,y,X,link,tnp) {
                      n = dim(X)[1]
                      p = dim(X)[2]
                      M = max( as.numeric(y) )
                      delta = c(delta_tilde[1], delta_tilde[1]+cumsum(delta_tilde[-1]^2), Inf)
                      delta_0_M = c(-Inf, delta[-M])
                      Delta = matrix( rep(delta,n),n,M,T )
                      Delta_0_M = matrix( rep(delta_0_M,n),n,M,T )
                      Xb = X%*%beta
                      XB = matrix( rep(X%*%beta,M),n,M )
                      -(1/tnp)*log( (1/n)*sum( drop(pfunc(delta[y]-Xb,method=link)-pfunc(delta_0_M[y]-Xb,method=link))^tnp ) ) +
                        (1/(1+tnp))*log( (1/n)*sum( (pfunc(Delta-XB,method=link)-pfunc(Delta_0_M-XB,method=link))^(1+tnp) ) )
                    } )
  
  # Optimization -------------------------------------------------------
  time = proc.time()[3]
  beta = init_beta
  delta_tilde = init_delta_tilde
  for (r in 1:maxIter) {
    old_beta = beta
    old_delta_tilde = delta_tilde
    
    
    opt_beta = switch(divergence,
                      "KL" = optim(par=beta, fn=obj_rorm, 
                                   delta_tilde=delta_tilde,
                                   y=y, X=X, link=method),
                      "DP" = optim(par=beta, fn=obj_rorm, 
                                   delta_tilde=delta_tilde,
                                   y=y, X=X, link=method, tnp=tnp),
                      "gamma" = optim(par=beta, fn=obj_rorm, 
                                      delta_tilde=delta_tilde,
                                      y=y, X=X, link=method, tnp=tnp))
    beta = opt_beta$par
    
    opt_delta_tilde = switch(divergence,
                      "KL" = optim(par=delta_tilde, fn=obj_rorm,
                                   beta=beta,
                                   y=y, X=X, link=method),
                      "DP" = optim(par=delta_tilde, fn=obj_rorm,
                                   beta=beta,
                                   y=y, X=X, link=method, tnp=tnp),
                      "gamma" = optim(par=delta_tilde, fn=obj_rorm,
                                      beta=beta,
                                      y=y, X=X, link=method, tnp=tnp))
    delta_tilde = opt_delta_tilde$par
    
    b_diff = sqrt(crossprod(beta-old_beta)/p)
    d_diff = sqrt(crossprod(delta_tilde-old_delta_tilde)/{M-1})
    if( b_diff<tol & d_diff<tol ) {
      break
    }else if(r%%20==0){
      cat(r,"th iteration", "\n",
          proc.time()[3]-time, "seconds gone.", "\n",
          "convergence criterion values: \n", 
          " for beta: ", b_diff, "\n", 
          " for delta: ", d_diff, "\n", "\n")
    } 
  }
  delta = c(delta_tilde[1], delta_tilde[1]+cumsum(delta_tilde[-1]^2))
  
  # results -------------------------------------------------------
  list(
    Coefficients = beta,
    Cutpoints = delta
  )
}

# system.time(
#   robustORM(y,X,divergence="KL")
# )
# system.time(
#   ordinal::clm(formula=formula, data=data, link="probit")
# )
# robustORM(y,X,divergence="DP",tnp=0.5)
# robustORM(y,X,divergence="gamma",tnp=0.5)


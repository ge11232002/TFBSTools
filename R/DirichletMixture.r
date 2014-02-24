### -----------------------------------------------------------------
### Train the Dirichlet mixture model from matrice
### Exported!

setMethod("dmmEM", signature(x="matrix"),
          function(x, K, alpha0=NULL, pmix=NULL){
            dirichletMixtureEMEstimation(t(x), K, alpha0, pmix)
          }
          )

setMethod("dmmEM", signature(x="PFMatrixList"),
          function(x, K, alpha0=NULL, pmix=NULL){
            allMatrix <- do.call(cbind, Matrix(x))
            dirichletMixtureEMEstimation(t(allMatrix), K, alpha0, pmix)
          }
          )

setMethod("dmmEM", signature(x="ANY"),
          function(x, K, alpha0=NULL, pmix=NULL){
            allMatrix <- getMatrixSet(x, opts=list(all=TRUE))
            dmmEM(allMatrix, K, alpha0, pmix)
          }
          )


dirichletMixtureEMEstimation <- function(inputMatrix, K, 
                                         alpha0=NULL, pmix=NULL){
## inputMatrix: the samples summarized as counts for each of A letters, N x A
  ## A <- 4 for DNA. This matrix should be concatenated by all matrices 
  ## from Jaspar.
## K: the number of sought component.
## alpha0: the estimated Dirichlet parameters A x K
## pmix: mixing proportions 1 x K
  K = as.integer(K)
  N <- nrow(inputMatrix)
  A <- ncol(inputMatrix)
  rowSumsInputMatrix <- rowSums(inputMatrix)
  oN <- rep(1, N)
  oK <- rep(1, K)
  oA <- rep(1, A)

  gam <- matrix(0, nrow=N, ncol=K)
  contrib_n <- matrix(0, nrow=A, ncol=K)
  contrib_d <- rep(0, K)

  ## random initialization
  if(is.null(alpha0)){
    epsilon <- 0.1
    alpha0 <- epsilon * (2 * matrix(runif(A * K), nrow=A, ncol=K) - 1) + 2
  }else{
    stopifnot(all(dim(alpha0) == c(nrow(inputMatrix), K)))
  }
  Alpha0 <- colSums(alpha0)
  if(is.null(pmix)){
    pmix <- rep(1, K) / K
  }else{
    stopifnot(sum(pmix) == 1)
  }

  ## minimum alpha0 value
  alpha0_min <- 1e-8
  ftol <- 1e-10 
  iteouter_max <- 10000 
  iteinner_max <- 1
  dll_min <- 1e-6 

  ite <- 0; 
  dll <- Inf;
  ll <- numeric(iteouter_max)
  while(ite < iteouter_max && (dll > dll_min || ite == 1)){
    ite <- ite + 1;

    ## E-step
    contrib <- log(pmix) + lgamma(Alpha0) - colSums(lgamma(alpha0))
    for(k in 1:K){
      alpha0_k <- alpha0[ , k]
      gam[ ,k] <- rowSums(lgamma(inputMatrix + 
                                 matrix(rep(alpha0_k, N), nrow=N, ncol=A, 
                                        byrow=TRUE))) -
                     lgamma(rowSumsInputMatrix + rep(Alpha0[k], N))
    }
    gam <- gam + matrix(rep(contrib, N), nrow=N, ncol=K, byrow=TRUE)
    maxgam <- apply(gam, 1, max)
    gam <- exp(gam - matrix(rep(maxgam, K), nrow=N, ncol=K, byrow=FALSE))

    ll[ite] <- sum(maxgam + log(rowSums(gam)))

    dll <- ll[ite] - ll[max(ite-1,1)]
    sumgam <- 1 / rowSums(gam)
    gam <- gam * matrix(rep(sumgam, K), nrow=N, ncol=K, byrow=FALSE)

    ## M-step
    pmix <- colSums(gam) / N
    dalpha0 <- Inf  
    iteinner <- 0 
    while(sum(abs(dalpha0)) > ftol && iteinner < iteinner_max){
      iteinner = iteinner + 1
      for(k in 1:K){
        alpha0_k <- alpha0[ , k]
        contrib_n[ , k] <- psigamma(t(inputMatrix) + 
                             matrix(rep(alpha0_k, N), nrow=A, 
                                    ncol=N, byrow=FALSE)) %*% gam[ , k]
        contrib_d[k] <- sum(psigamma(rowSumsInputMatrix + Alpha0[k]) * 
          gam[ ,k]) - pmix[k] * N * psigamma(Alpha0[k])
      }
      contrib_n <- contrib_n - N * psigamma(alpha0) * 
        matrix(rep(pmix, A), nrow=A, ncol=length(pmix), byrow=TRUE)
      dalpha0 <- pmax(alpha0 * (contrib_n / 
                              matrix(rep(contrib_d, A), nrow=A, 
                                     ncol=length(contrib_d), byrow=TRUE)), 
                    alpha0_min) - alpha0
      alpha0 <- alpha0 + dalpha0 
      Alpha0 <- colSums(alpha0)
    }
    if(isTRUE(all.equal(ite %% iteouter_max / 10, 0))){
      print(dll)
      print(sum(abs(dalpha0)))
      print(pmix)
      print(alpha0)
    }
  }
  return(list(alpha0=alpha0, pmix=pmix, ll=ll))
}

repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}

## Sampling from Dirichlet distribution is implemented in gtools::rdirichlet.
#dirichletSample <- function(a, n=1){
  ## Sample from Dirichlet distribution.
  ## DIRICHLET_SAMPLE(a) returns a probability vector sampled from a
  ## Dirichlet distribution with parameter vector A.
  ## DIRICHLET_SAMPLE(a,n) returns N samples, collected into a matrix, each
  ## vector having the same orientation as A.
  ## References:
  ##   [1]  L. Devroye, "Non-Uniform Random Variate Generation",
  ##   Springer-Verlag, 1986
  ## This is essentially a generalization of the method for Beta rv's.
  ## Theorem 4.1, p.594
  #rgamma(n*length(a), rep(a, n))
#}

### -----------------------------------------------------------------
### sample the matrix from Dirichlet mixture model
### Exported!

setMethod("rPWMDmm", signature(x="matrix"),
          function(x, alpha0, pmix, N=1, W=6){
            PWMrandomizeBayes(x, alpha0, pmix, N, W)
          }
          )

setMethod("rPWMDmm", signature(x="PFMatrix"),
          function(x, alpha0, pmix, N=1, W=6){
            PWMrandomizeBayes(Matrix(x), alpha0, pmix, N, W)
          }
          )

setMethod("rPWMDmm", signature(x="PFMatrixList"),
          function(x, alpha0, pmix, N=1, W=6){
            allMatrix <- do.call(cbind, Matrix(x))
            PWMrandomizeBayes(allMatrix, alpha0, pmix, N, W)
          }
          )

#setMethod("rPWMDmm", signature(x="ANY"),
#          function(x, alpha0, pmix, N=1, W=6){
#            allMatrix <- getMatrixSet(x, opts=list(all=TRUE))
#            rPWMDmm(allMatrix, alpha0, pmix, N, W)
#          }
#          )


PWMrandomizeBayes <- function(PCM, alpha0, pmix, N=1, W=6){
  ## generates N (default 1) random PWM of width drawn from 
  ## the posterior distribution
  ## of PWMs. The posterior is propertional to the 
  ## Dirichlet prior distribution Dprior (which
  ## might be a mixture) times the mulitnomial likelihood with count matrix PCM.
  A <- nrow(PCM)
  Win <- ncol(PCM)
  ### parameters of prior
  #alpha0 <- Dprior$alpha0 
  ##alpha0 is the estimated Dirichlet parameters A x K
  #pmix <- Dprior$pmix
  ## pmix mixing proportions 1 x K
  K <- length(pmix)

  ### accumulative distribution for mixing proportions
  apmix <- cumsum(pmix)

  PWM = list() 
  for(n in 1:N){
    PWM[[n]] <- matrix(NA, ncol=W, nrow=nrow(alpha0))
    if(W == Win){
      for(w in 1:W){
        ## draw from the mixture component
        k = which(runif(1) < apmix)[1]
        ## draw from component k of dirichlet posterior
        PWM[[n]][ ,w] <- 
          rdirichlet(n=1, alpha=(alpha0[ ,k] + PCM[ ,w]))
      }
    }else{
      for(w in 1:W){
        ## draw from the mixture component
        k = which(runif(1) < apmix)[1]
        ## draw from component k of dirichlet posterior and use random
        ## column in PCM
        PWM[[n]][, w] <- 
          rdirichlet(n=1, alpha=(alpha0[ ,k] + 
                                         PCM[ ,ceiling(Win * runif(1))]))
      }
    }
  }
  stopifnot(all(colSums(PWM) == 1))
  return(PWM)
}



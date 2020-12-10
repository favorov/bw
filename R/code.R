#' @export
#---------- Compute Wasserstein barycenter -----
Wbarycenter <- function(population, weights=NULL, err=1e-6) {
  ##Compute the Wasserstein barycenter of a population GPs.
  #population is a list with m Gaussian distributions given by list-type 
  #consisting of 'mean' and 'covariance' variables
  ##weights contains weights
  ##err is the error margin
  
  m <- length(population) #num of Gaussians
  #If weights are not specified, uniform weights are chosen.
  if (is.null(weights)) {
    weights <- rep(1, m) / m
  }
  ##Barycenter of means
  ##
  d        <- length(population[[1]]$mean) #Dimension of the Gaussians
  mean.bar <- rep(0, d) #barycenter of means
  for (i in 1:m) {
    mean.bar <- mean.bar + weights[i] * population[[i]]$mean
  }
  
  ## Barycenter of covariances
  covs <- list() #list of covariances
  for (i in 1:m) {
    covs[[i]] <- population[[i]]$covariance
  }
  
  #The barycenter is a solution of a fixed-point equation
  uplimit <- 1e3 #num of iterations
  Q       <- diag(d)
  count   <- 0;
  repeat {
    sqrt.Q    <- Sqrtm(Q)
    sqrtinv.Q <- SqrtmInv(Q) 
    mean.T <- matrix(rep(0, d^2), d, d)
    for (i in 1:m) {
      mean.T <- mean.T + weights[i] * Sqrtm(sqrt.Q %*% covs[[i]] %*% sqrt.Q)
    }
    mean.T <- sqrtinv.Q %*% mean.T %*% sqrtinv.Q
    Q <- mean.T %*% Q %*% mean.T
    count  <- count + 1;
    if (norm(mean.T - diag(d), "f") < err || count == uplimit) {
      break
    }
  }
  if (count == uplimit) {
    print('Iterations number exceeded!')
  } else {
    #print(paste('Iterations number:', count))
  }
  #population barycenter
  mean.gp <- list('mean' = mean.bar, 'covariance' = Q)
  return(mean.gp)
}

#' @export
Wbarycenter_2 <- function(population, weights=NULL, err=1e-6) {
  ##Compute the Wasserstein barycenter of a population GPs.
  #population is a list with m Gaussian distributions given by list-type 
  #consisting of 'mean' and 'covariance' variables
  ##weights contains weights
  ##err is the error margin
  
  m <- length(population) #num of Gaussians
  #If weights are not specified, uniform weights are chosen.
  if (is.null(weights)) {
    weights <- rep(1, m) / m
  }
  ##Barycenter of means
  ##
  d <- length(population[[1]]$mean) #Dimension of the Gaussians
  g <- d * (d + 1) / 2
  mean.bar <- rep(0, d) #barycenter of means
  for (i in 1:m) {
    mean.bar <- mean.bar + weights[i] * population[[i]]$mean
  }
  
  ## Barycenter of covariances
  covs <- list() #list of covariances
  for (i in 1:m) {
    covs[[i]] <- population[[i]]$covariance
  }
  
  #The barycenter is a solution of a fixed-point equation
  uplimit <- 1e3 #num of iterations
  Q <- diag(d)
  next.Q  <- Fpm(Q, covs, weights)
  count   <- 0;
  while (W(list(0, Q), list(0, next.Q)) > err && count < uplimit) {
    Q <- next.Q
    D <- Sqrtm(Q)
    T.mean <- matrix(rep(0, d^2), d, d)
    dT.mean <- matrix(rep(0, g^2), g, g)
    for (i in 1:m) {
      T.mean <- T.mean + weights[i] * Sqrtm(D %*% covs[[i]] %*% D)
      dT.mean <- dT.mean + weights[i] * ComputeReprdT(Q, covs[[i]])
    }
    T.mean <- MInv(D) %*% T.mean %*% MInv(D)
    T.diff.repr <- GetRepr(T.mean - diag(d))
    Delta.repr <- - solve(dT.mean, T.diff.repr)
    Delta <- Reconstruct(Delta.repr)
    b <- min(eigen(MInv(D) %*% Delta %*% MInv(D), symmetric = TRUE)$values)
    a <- (4 - b) / (2 - b)^2
    next.Q <- Q + a * Delta
    
    count  <- count + 1;
  }
  if (count == uplimit) {
    print('Iterations number exceeded!')
  } else {
    print(paste('Iterations number:', count))
  }
  #population barycenter
  mean.gp <- list('mean' = mean.bar, 'covariance' = next.Q)
  return(mean.gp)
}

#---------- Compute Wasserstein distance -----
W <- function(f1,f2) {
  #Compute the 2-Wasserstein distance between Gaussian distributions
  #f1 = N(m1,K1) and f2 = N(m2,K2); f1, f2 are lists
  
  m1 <- f1[[1]] 
  m2 <- f2[[1]]
  K1 <- f1[[2]] 
  K2 <- f2[[2]]
  
  #2-W distance
  cov.dist <- BW(K1, K2)
  l2norm   <- norm(as.matrix(m1 - m2), 'F')
  d        <- sqrt(cov.dist^2 + l2norm^2)
  return(d)
}

#' @export
BW <- function(K1, K2) {
  #Compute the Bures-Wasserstein distance between covariance matrices
  return(sqrt(abs(Tr(K1) + Tr(K2) - 2 * Tr(Sqrtm(Sqrtm(K1) %*% K2 %*% Sqrtm(K1))))))
} 

bootstrap <- function(population, alpha, M) {
  #population is an observed sample
  # alpha is a confidence level, by default set to 5%
  # M is num of iterations in resampling, by default set to 1000
  
  if (is.null(alpha)) {
    alpha <- .05
  }
  if (is.null(M)) {
    M <- 1000
  }
  stat  <- rep(0, M)
  for (m in 1 : M) {
    weights <- rpois(n, 1)
    weights <- weights / sum(weights)
    bw_boot <- Wbarycenter(population, weights = weights)
    stat[m] <- BW(bw$covariance, bw_boot$covariance)
  }
 return(quantile(stat, 1-alpha))
}

#' @export
genRanGauss <- function(d, k) {
  #generats a test population of Gaussian distributions
  #n is the dimensionality of GPs
  #k is the number of GPs
  out <- list()
  for (i in 1 : k) {
    A        <- matrix(runif(d^2) + 100, ncol = d) 
    sigma    <- crossprod(A) # A' A
    E        <- eigen(sigma, symmetric = TRUE)
    values   <- abs(rnorm(d, 0, 4))^2
    #for (j in 1 : d) {
    #    values[j] <- values[j] + j
    #}
    Sigma    <- t(E$vectors) %*% diag(abs(values)) %*% E$vectors #diag(rnorm(n, 20, 3)) #
    #out[[i]] <- list('mean' = rnorm(d), 'covariance' = Sigma)
    out[[i]] <- list('mean' = rep(0, d), 'covariance' = Sigma)
  }
  return(out)
}

# Example: population <- genRanGauss(10, 5)

#' @export
GenONbasis <- function(d) {
  #Generate orthonormal basis in space of d \tims d symmetric matrices
  out     <- list()
  counter <- 0
  for (i in 1 : d) {
    ei    <- rep(0, d)
    ei[i] <- 1
    for (j in 1 : i) {
      ej      <- rep(0, d)
      ej[j]   <- 1
      counter <- counter + 1
      
      if (i == j) {
        out[[counter]] <- list('count' = counter, 'u' = ei %o% ei)
      } else {
        out[[counter]] <- list('count' = counter, 'u' = (ei %o% ej + ej %o% ei) / sqrt(2))
      }
    }
  }
  return(out)
}


#---------- Auxiliary functions -----
Fpm <- function(K, covs, weights) {
  #step in the solution of fixed-point equation
  sqrt.K    <- Sqrtm(K)
  sqrtinv.K <- SqrtmInv(K) 
  d         <- dim(covs[[1]])[1]
  m         <- length(weights)
  M         <- matrix(rep(0, d^2), d, d)
  for (i in 1:m) {
    M <- M + weights[i] * Sqrtm(sqrt.K %*% covs[[i]] %*% sqrt.K)
  }
  M <- sqrtinv.K %*% Sqrm(M) %*% sqrtinv.K
  return(M)
}

L2Norm <- function(X, Y, B) {
  #Frobenius norm of matrix X w.r.t matrix B
  sqrti.B <- SqrtmInv(B)
  sqrt.B  <- Sqrtm(B)
  R       <- Sqrtm(sqrt.B %*% X %*% sqrt.B) - Sqrtm(sqrt.B %*% Y %*% sqrt.B)
  Z       <- sqrti.B %*% R %*% sqrti.B
  out     <- Tr(Z)
  return(out)
}

Tr <- function(X) {
  #computes trace
  out <- sum(diag(X))
  return(out)
}

Sqrtm <- function(X) {
  #computes square root of a symmetric matrix
  E <- eigen(X, symmetric = TRUE) 
  V <- E$values
  Q <- E$vectors 
  Y <- Q %*% diag(sqrt(V)) %*% t(Q) 
  return(Y)
}

Sqrm <- function(X) {
  #computes a squared symmetric matrix
  E <- eigen(X, symmetric = TRUE)
  V <- E$values 
  Q <- E$vectors 
  Y <- Q %*% diag((V)^2) %*% t(Q) 
  return(Y)
}

SqrtmInv <- function(X) {
  #Inversion of a square root of X
  E <- eigen(X, symmetric = TRUE) 
  V <- E$values 
  Q <- E$vectors 
  Y <- Q %*% diag(1 / sqrt(V)) %*% t(Q) 
  return(Y)
}

MInv <- function(X) {
  #Inversion of X
  E <- eigen(X, symmetric = TRUE) 
  V <- E$values 
  Q <- E$vectors 
  Y <- Q %*% diag(1 / V) %*% t(Q) 
  return(Y)  
}

GetOTmap <- function(Q, S) {
  #Generate OT map from N(0, Q) to N(0, S)
  D  <- SqrtmInv(Q)
  K  <- Sqrtm(Q)
  G  <- Sqrtm(K %*% S %*% K)
  OT <- D %*% G %*% D
  return(OT)
}




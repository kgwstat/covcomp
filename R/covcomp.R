# COMPLETION
# Calculate completion using truncation

#' Pivot
#'
#' Finds the first pivot of Omega.
#' @param Omega A matrix of 0s and 1s.
#' @return Vector of dimension 2: entries of the pivot.
#' @examples
#' Omega <- matrix(c(1,1,0,1,1,1,0,1,1), ncol = 3)
#' Pivot(Omega)
Pivot = function(Omega) {
  # Calculates pivot of Omega.
  i <- 1
  j <- 1
  n <- nrow(Omega)
  
  while ((i <= n) && (Omega[i,j] == 1)) i <- i + 1
  if(i > n) stop('no pivots found')
  while ((j <= n) && (Omega[i,j] == 0)) j <- j + 1
  c(i, j-1)
}

#' C Index
#'
#' Finds the C Index of Omega.
#' @param Omega A matrix of 0s and 1s.
#' @return Vector; range of C Index
#' @examples
#' Omega <- matrix(c(1,1,0,1,1,1,0,1,1), ncol = 3)
#' IndexCof(Omega)
IndexCof = function(Omega) {
  # Calculates the IndexC of Omega
  n <- nrow(Omega)
  p <- Pivot(Omega)
  i <- p[1]
  j <- p[2]+1
  
  while((i < n) && (Omega[i,j] == 1)) i <- i + 1
  p[1]: if(i == n) n else i-1
} 

#' Completion
#'
#' Calculates the canonical completion of partial covariance assuming the domain is given by Omega.
#' @param Cov The partial covariance matrix with same dimensions as Omega.
#' @param Omega A matrix of 0s and 1s; indicating the domain of partial covariance.
#' @param VectorOfTuning vector of truncation parameters with dimension at least m-1 where Omega is an m-serrated domain
#' @return A matrix of the same dimensions: the covariance completion of Cov
Compln = function(Cov, Omega, VectorOfTuning) {
  # Returns the Completion of Cov for the given Omega and VectorOfTuningParameters
  
  # Return Cov if Omega is trivial.
  if (all(Omega == 1)) {
    return(Cov)
  }
  
  # Calculating indices for matrix inversion and multiplication.
  p <- Pivot(Omega)
  i <- p[1]
  j <- p[2]
  
  IndexA <- 1:j
  IndexB <- (j+1):(i-1)
  IndexC <- IndexCof(Omega)
  
  # Calculating truncated inverse of Cov[IndexB, IndexB] of rank r
  r <- VectorOfTuning[1]
  E <- eigen(Cov[IndexB,IndexB])
  D <- diag(E$values[1:r]^{-1},r,r)
  U <- E$vectors[,1:r]
  Binv <- U %*% D %*% t(U)
  
  # Calculating part of completion and 
  # Updating Cov
  R <- Cov[IndexA, IndexB] %*% Binv %*% Cov[IndexB, IndexC]
  Cov[IndexA,IndexC] <- R
  Cov[IndexC,IndexA] <- t(R)
  
  # Updating Omega
  Omega[IndexC,1:j] <- 1
  Omega[1:j,IndexC] <- 1
  
  # Continue with updated Cov and Omega
  Compln(Cov, Omega, VectorOfTuning[-1])
}


#' Completion
#' 
#' Calculates the canonical completion of partial covariance assuming the domain is given by Omega.
#' @param Cov The partial covariance matrix with same dimensions as Omega.
#' @param Omega A matrix of 0s and 1s; indicating the domain of partial covariance.
#' @param FVE Fraction of variance explained; 0.95 by default.
#' @return A matrix of the same dimensions: the covariance completion of Cov

ComplnFVE = function(Cov, Omega, FVE = 0.95) {
  # Returns the Completion of Cov for the given Omega and tunes according to FVE criterion.
  
  # Return Cov if Omega is trivial.
  if (all(Omega == 1)) {
    return(Cov)
  }
  
  # Calculating indices for matrix inversion and multiplication.
  p <- Pivot(Omega)
  i <- p[1]
  j <- p[2]
  IndexA <- 1:j
  IndexB <- (j+1):(i-1)
  IndexC <- IndexCof(Omega)
  
  # Calculating truncated inverse of Cov[IndexB, IndexB] of rank r
  M <- Cov[IndexB,IndexB]
  E <- eigen(M)
  
  trace <- sum(diag(M)) / nrow(M)
  r <- 1
  while(sum(E$values[1:r]) < FVE * trace) {r <- r + 1}
  
  D <- diag(E$values[1:r]^{-1},r,r)
  U <- E$vectors[,1:r]
  Binv <- U %*% D %*% t(U)
  
  # Calculating part of completion and 
  # Updating Cov
  R <- Cov[IndexA, IndexB] %*% Binv %*% Cov[IndexB, IndexC]
  Cov[IndexA,IndexC] <- R
  Cov[IndexC,IndexA] <- t(R)
  
  # Updating Omega
  Omega[IndexC,1:j] <- 1
  Omega[1:j,IndexC] <- 1
  
  # Continue with updated Cov and Omega
  ComplnFVE(Cov, Omega, FVE)
}
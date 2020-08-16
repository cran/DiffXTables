# Method to determine second order (interaction) differential patterns
# Input : a list of sampled matrices, must be greater or equal to 2
# Output: Statistics and p-value, significance of differentiality
# Created by: Ruby Sharma and Dr. Joe Song
# Date Created: 28 December 2018
# Modified: 
#   12 June 2019 by Dr. Joe Song
# 
#   October 29, 2019:
#     - Added ifelse to check if U is a matrix, if, yes, multiply by 
#       nrow(U) else multiply by 1
#     - Added a test case for 2x2 tables
#     - Renamed matrix U to Q to be consistent with the manuscript.
#     - Simplified the calculation of Df.
# 
#   March 18, 2020: The rank of the covariance matrix is always
#     K-1. No need to compute the rank numerically.
#
#   May 14, 2020: Corrected the normalized frequence table A inside 
#     function get.e.mat() when an entry has both zero row and zero 
#     column sums.

#######################################
#### Imports for the whole package ####
# No longer need this import 
#XXXXXXXXXX' @importFrom Matrix rankMatrix
#######################################

#' @export
sharma.song.test <- function(tables)
{
  if (mode(tables) != "list" || length(tables) < 2) {
    stop("Input must be a list of 2 or more matrices!")
  }

  for(k in 1:(length(tables)-1)){
    if(!identical(dim(tables[[k]]), dim(tables[[k+1]])))
    stop("All matrices must have identical dimensions!")
  }

  # Get independent standard normal variables (E matrix) using Helmert tranfrom 
  EN <- get.e.mat(tables)
  n <- EN$SamS

  if(all(n >= 0) && any(n > 0)) {

    K <- length(tables)

    EList <- EN$EList

    EP <- Reduce('+', EList) # Pooled E vectors of all conditions

    # Get the scaling coffecients:
    b <- sqrt(n) / sum(sqrt(n))

    # Difference of individual E vectors from E pooled:
    Q <- sapply(seq(K), function(k){
      return( EList[[k]] - b[k] * EP )
    })

    # Calculate covariance matrix:
    C <- cpc.cov.matrix(b)

    # Find the rank of the cov matrix
    R <- K - 1 # rankMatrix(C) 

    # Perform eigen-decomposition of covariance matrix:
    eig <- eigen(C)

    nonzero.eigenvalues <- eig$values[seq(R)]

    # Eigenvectors corresponding to non-zero eigenvalues:
    S <- eig$vectors[, seq(R), drop=FALSE]

    Z <- diag( 1 / sqrt( nonzero.eigenvalues ), nrow = R, ncol = R)

    # Calculating Mahalanobis distance squared from eigenvalue
    #   decomposition since rank is not full
    Stat <- sum( ( Q %*% S %*% Z ) ^ 2 )

    Df <- R * length(EP)  # Df <- R * ifelse(is.matrix(Q), nrow(Q), 1)
      
    P.val <- pchisq(Stat, Df, lower.tail = FALSE)

  } else {

    Stat = 0
    Df = 0
    P.val = 1

  }

  names(Stat) <- "X-squared"
  names(Df) <- "df"
  DNAME <- deparse(substitute(tables))

  return(structure(list(
    statistic = Stat,
    parameter = Df,
    p.value = P.val,
    data.name = DNAME,
    method = "Sharma-Song Test for Second-Order Differential Contingency Tables"),
    class = "htest"))
}


get.e.mat <- function(tables)
{ 
  K <- length(tables)
  EList <- vector("list", length = K)
  n <- vector("numeric", length = K)

  for(k in seq(K)) {

    # Obtain row and column Helmert matrices
    rowSum <- rowSums(tables[[k]])  # (expec)
    colSum <- colSums(tables[[k]])  # (expec)
    totalSum <- sum(rowSum)

    prdot <- rowSum / totalSum
    pcdot <- colSum / totalSum

    V <- helmert.matrix(prdot)
    W <- helmert.matrix(pcdot)

    # Calculate expected table
    expec <- expected(tables[[k]])
    
    # Get normalized frequency table of sampled data
    A <- (tables[[k]] - expec) / sqrt(ifelse(expec == 0, 1, expec))
    
    # Correct the frequency for zero entries with both 
    #   a zero row sum and a zero column sum.
    for(r in seq_along(rowSum)) {
       if(rowSum[r] != 0) next 
      for(c in seq_along(colSum)) {
         if(colSum[c] == 0) A[r, c] <- sqrt(totalSum)
      }
     }

    # Transform normalized frequency matrix A to independent normal 
    #   variables by multiplying row and column helmert matrices
    E <- V %*% A %*% t(W)

    # Remove the first row and first column from E, which
    #   are alwarys zero. The other entries in E are random normal
    #   variables.
    E <- E[-1, -1]

    # Vectorize E in column major 
    EList[[k]] <- as.vector(E)

    n[k] <- totalSum
  }
  return(list(SamS = n, EList = EList))
}


expected <- function(table)
{
  rowSum = rowSums(table)
  colSum = colSums(table)
  totalSum = sum(rowSum)
  prod = outer(rowSum, colSum, "*")

  if(totalSum != 0) {
    Exp <- prod / totalSum
  } else {
    Exp <- prod
  }
  
  return (Exp)
}


helmert.matrix <- function(p)
{
  n <- length(p)
  L <- matrix(0, nrow=n, ncol=n)
  L[1, ] <- sqrt(p)

  sumpi.minus.1 <- p[1]

  if(n >= 2) {
    for(i in 2:n) {
      sumpi <- sumpi.minus.1 + p[i]

      # L[i,i] <- - sqrt(sumpi.minus.1  / sumpi) 
      L[i,i] <- - sqrt(
       ifelse(sumpi.minus.1 == 0,
              0,
              sumpi.minus.1  / sumpi) )

      # L[i,1:(i-1)] <- sqrt(p[i] * p[1:(i-1)] / (sumpi.minus.1 * sumpi))
      
      L[i,1:(i-1)] <- sqrt(
        ifelse( p[i] * p[1:(i-1)] == 0,
                0,
                p[i] * p[1:(i-1)] / (sumpi.minus.1 * sumpi)
              ) )

      sumpi.minus.1 <- sumpi
    }
  }
  return(L)
}


cpc.cov.matrix <- function(b)
{
  # bj > 0 and sum bj=1
  K <- length(b)
  JK <- matrix(1, nrow=K, ncol=K)
  covmat <- diag(K) - JK %*% diag(b) - diag(b) %*% JK + K * b %*% t(b)
  return(covmat)
}


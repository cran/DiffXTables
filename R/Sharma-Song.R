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
    stop("All matries must have identical dimensions!")
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

    # Calculating expected table and normalized table of sampled data

    expec <- expected(tables[[k]])

    A <- (tables[[k]] - expec) / sqrt(ifelse(expec == 0, 1, expec))

    rowSum <- rowSums(expec)
    colSum <- colSums(expec)
    totalSum <- sum(rowSum)

    prdot <- rowSum / totalSum
    pcdot <- colSum / totalSum

    # Obtaining row and column Helmert matrices
    V <- helmert.matrix(prdot)
    W <- helmert.matrix(pcdot)

    # Product of row and column helmert matrix with normalized sampled data
    E <- V %*% A %*% t(W)

    E <- E[-1, -1]

    # List of E vector
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

  # t <-  ifelse(prod==0, 0, prod/totalSum)
  # return (t)
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

      L[i,i] <- - sqrt(ifelse (sumpi.minus.1 == 0, 0, sumpi.minus.1  / sumpi) )

      L[i,1:(i-1)] <- sqrt( ifelse( (p[i]*p[1:(i-1)]) == 0, 0, (p[i]*p[1:(i-1)])  / (sumpi.minus.1 * sumpi) ) )

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


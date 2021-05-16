# marginal-change-test.R
#
# Method to determine first-order (marginal) differential patterns.
# Input : a list of sampled matrices, must be greater or equal to 2
# Output: Statistics and p-value, significance of marginal differentiality
# Created by: Ruby Sharma and Dr. Joe Song
# Date Created: May 17 2020
# Updated: August 2, 2020. MS. Turned off correction on 2x2 tables.

#' @importFrom  Rdpack reprompt

#' @export
marginal.change.test <- function(tables)
{
  K <- length(tables)
  if (mode(tables) != "list" || K < 2) {
    stop("Input must be a list of 2 or more matrices!")
  }
  
  for (k in seq(K-1)) {
    if (!identical(dim(tables[[k]]), dim(tables[[k + 1]]))) 
      stop("All matries should have same number of rows and columns")
  }
  
  ## Constructing row marginal table
  r = nrow(tables[[1]])
  
  Xc = matrix(nrow = K, ncol = r)
  for (i in seq(K)) {
    Xc[i, ] = rowSums(tables[[i]])
  }
  
  ## Constructing column marginal table
  s = ncol(tables[[1]])
  Yc = matrix(nrow = K, ncol = s) 
  for (i in seq(K)) {
    Yc[i, ] = colSums(tables[[i]])
  }
  
  ## Evaluating Chisq of all the tables 
  Xchisq = 0
  Ychisq = 0
  if(sum(Xc) > 0 || sum(Yc) > 0) {
    if(sum(Xc)!=0){
      Xchisq = chisq.test.stat(Xc)
    }
    if(sum(Yc)!=0){
      Ychisq = chisq.test.stat(Yc)
    }
    
    ## Determining the strength of marginal difference between the tables
    Stat = Xchisq+Ychisq	
    Df = (K-1)*((r+s)-2)
    P.val = pchisq(Stat , Df, lower.tail = FALSE)
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
    method = "Test for Marginal Change Across Contingency Tables"),
    class = "htest"))
  
}

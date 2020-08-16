# strength-test.R
#
# Method to determine strength of association (zero-order) in each table.
# Input : a list of sampled matrices, must be greater or equal to 2
# Output: Statistics and p-value, significance of association strength.
# Created by: Ruby Sharma 
# Date Created: May 17 2020
# Updated: August 2, 2020. MS

#' @export
strength.test <- function(tables)
{
  K = length(tables)
  
  if (mode(tables) != "list" || K < 2) {
    stop("Input tables must be a list of two or more matrices!")
  }
  
  KStats <- vector("numeric", K)
  df <- vector("integer", K)
  
  for(i in seq(K)){
    if(sum(tables[[i]]>0)){
      KStats[i] <- chisq.test.stat(tables[[i]])
    }else{
      KStats[i] <- 0
    }
    df[i] <- prod(dim(tables[[i]]) - 1)
  }
  
  Stat <- sum(KStats)
  Df <- sum(df)
  
  P.val <- pchisq(Stat, Df, lower.tail = FALSE)
  
  names(Stat) <- "X-squared"
  names(Df) <- "df"
  DNAME <- deparse(substitute(tables))
  
  return(structure(list(
    statistic = Stat,
    parameter = Df,
    p.value = P.val,
    data.name = DNAME,
    method = "Strength Test for Association in Multiple Contingency Tables"),
    class = "htest"))
}

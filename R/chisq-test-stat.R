# chisq-test-stat.R
#
# Extracted from heterogeneity.R and marginal-change-test.R 
#   to compute the chi-square statistic only 

# #' @importFrom stats chisq.test

# Calculate chi-square test after removing zero rows or columns
#chisq.test.stat <- function(table)
#{
#  non.zero.rows <- apply(table, 1, function(row) { 0 != sum(abs(row)) } )
#  non.zero.cols <- apply(table, 2, function(col) { 0 != sum(abs(col)) } )

#  perform Pearson chi-square test
#  table = table[non.zero.rows, non.zero.cols]
#  if(all(dim(table) > 1)) { # tables needs to have more than one row and one column
#    if(length(table) > 1) {
#    chisq <- suppressWarnings(
#      chisq.test(table, correct = FALSE)$statistic)
#  } else {
#    # stop("'Table' must have at least 2 elements")
#    chisq <- 0
#  }
  
#  return(chisq)
#}

# Calculate chi-square statistics
chisq.test.stat <- function(table){
  
  if(all(dim(table) > 1)) {
    nr <- nrow(table)
    nc <- ncol(table)
    sr <- rowSums(table)
    sc <- colSums(table)
    n <- sum(sr)
    E <- outer(sr, sc, "*")/n
    chisq <- sum((table-E)^2/E, na.rm=TRUE)
  } else {
    stop("'table' must be least 2x2")
    chisq <- 0
  }
  return(chisq)
}

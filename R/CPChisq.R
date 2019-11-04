# CPChisq.R -- comparative chi-square test for heterogeneous associations
#
# MS
# Created: July 3, 2015 generalized from two contingency tables in
#                       "hetero-chisq-steel.R"
# Modified: Feb 5, 2016. MS
#   Modified the method name
#   Added log.p argument to specify log p-value


#' @importFrom stats pchisq pnorm
#' @export
cp.chisq.test <- function(tables, method = c("chisq", "nchisq", "default", "normalized"), log.p=FALSE)
{
  DNAME <- deparse(substitute(tables))
  if(mode(tables)!="list" || length(tables)<2 )
  {
    stop("only accept list of 2 or more matrices as input!")
  }

  for(k in 1:(length(tables)-1)){
    if(!identical(dim(tables[[k]]), dim(tables[[k+1]])))
      stop("All matries should have same number of rows and columns")
  }

  method <- match.arg(method)
  if(method == "default") {
    warning(paste0("method=\"", method, "\" is deprecated. Use \"chisq\" instead."))
    method <- "chisq"
  } else if(method == "normalized") {
    warning(paste0("method=\"", method, "\" is deprecated. Use \"nchisq\" instead."))
    method <- "nchisq"
  }

  pooled <- tables[[1]]
  for(i in 2:length(tables)) {
    pooled <- pooled + tables[[i]]
  }

  # Row sum of matrix pooled
  pooled.row.sum <- apply(pooled, 1, sum)

  # Column sum of matrix pooled
  pooled.col.sum <- apply(pooled, 2, sum)

  #Sum of matrix pooled
  pooled.sum <- sum(pooled.col.sum)

  pooled.chisq <- 0

  if(pooled.sum > 0) {
    #Expected value in pooled
    pooled.expected <- pooled.row.sum %*% t(pooled.col.sum) / pooled.sum

    #Chi-square value of pooled
    pooled.chisq <- sum((pooled - pooled.expected)^2 / pooled.expected, na.rm = TRUE)

    total.chisq <- 0

    for(i in seq_along(tables)) {
      #Expected values in matrix i
      expected <- pooled.expected / pooled.sum * sum(tables[[i]])

      #Chi-square value of matrix i
      chisq <- sum((tables[[i]] - expected)^2 / expected, na.rm = TRUE)

      total.chisq <- total.chisq + chisq
    }

    hetero.chisq  <- total.chisq - pooled.chisq
    df <- (length(tables) - 1) * (sum(pooled.row.sum != 0) - 1) *
      (sum(pooled.col.sum != 0) - 1)

  } else {

    hetero.chisq <- 0
    df <- 0

  }

  if(method=="chisq") {

    finalStat <- hetero.chisq
    names(finalStat) <- "X-squared"

    finalDf <- df
    names(finalDf) <- "df"

    p.val <- pchisq(finalStat, df = finalDf, lower.tail=FALSE, log.p=log.p)
    names(p.val) <- "p.value"

    return(structure(list( statistic=finalStat, parameter=finalDf, p.value=p.val,
                           method = "Comparative chi-squared test",
                           data.name= DNAME),
                     class = "htest"))

  } else if(method=="nchisq") {

    finalStat <- ifelse(df > 0, ( hetero.chisq - df) / sqrt( 2 * df ), 0)
    names(finalStat) <- "Normalized X-squared"

    finalDf <- df
    names(finalDf) <- "df"

    p.val <- ifelse(df > 0, pnorm( finalStat, lower.tail=FALSE, log.p=log.p), 1)
    names(p.val) <- "p.value"

    return(structure(list(statistic = finalStat, parameter = finalDf, p.value = p.val,
                          method = "Nomalized comparative chi-squared test",
                          data.name= DNAME),
                     class = "htest"))
  }
}

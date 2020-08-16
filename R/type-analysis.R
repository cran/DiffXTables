# type-analysis.R
#
# Method to determine Type-0, Type-1 and Type-2 patterns in the set of 
# contingency tables.
# Input : a list of sampled matrices, must be greater or equal to 2
# Output: Statistics and p-value, significance of zero-order, first-order and
#         second-order test.
# Created by: Ruby Sharma and Dr. Joe Song
# Date Created: July 20 2020
# Updated: August 2, 2020 MS

#' @importFrom  pander pandoc.emphasis
#' @importFrom  pander pander

#' @export
type.analysis <- function(tables, alpha = 0.05)
{
  K <- length(tables)
  if (mode(tables) != "list" || K < 2) {
    stop("Must input a list of two or more matrices!")
  }
  
  for (k in seq(K-1)) {
    if (!identical(dim(tables[[k]]), dim(tables[[k + 1]]))) 
      stop("All matries must have the same dimension!")
  }
  
  SecondOrd <- sharma.song.test(tables)
  MarginalChange <- marginal.change.test(tables)
  StrengthCond <- strength.test(tables)
  
  Type <- decide.type(
    StrengthCond$p.value, MarginalChange$p.value, SecondOrd$p.value, 
    alpha)
  
  return(structure(list(
    Zeroth.order = StrengthCond,
    First.order = MarginalChange,
    Second.order = SecondOrd,
    Type = Type,
    Method = "Type Analysis Across Multiple Contingency Tables"),
    class = "DiffXTableTypeAnalysis"
  ))
}

#' @export
print.DiffXTableTypeAnalysis <- function(x, ...)
{
  # Printing differential table type analysis 
  
  with(x, {
    # TypeNull <- "X and Y are statistically independent in all K conditions
            # and same X marginal distribution across conditions so does Y."
    
    # Constructing Result table
    Result = data.frame(matrix(nrow = 3, ncol = 3))
    Result[1,] = c(
      '0th', AddStars(Zeroth.order$p.value), 
      "association present in tables" )
    Result[2,] = c(
      '1st', AddStars(First.order$p.value), 
      "differential in marginal distribution across tables" )
    Result[3,] = c(
      '2nd', AddStars(Second.order$p.value), 
      "differential in deviation from joint distribution to product of marginals" )
    
    colnames(Result) = c("Order", "P value", "Description")
    
    pander("\nComprehensive Type Analysis of Multiple Contingency Tables:\n\n")
    
    if(is.null(Type)) {
      pander("\tType null: association absent")
    } else if(Type == 0) {
      pander("\tType 0: association present; patterns conserved")
    } else if(Type == 1) {
      pander("\tType 1: association present; difference up to first order")
    } else if(Type == 2) {
      pander("\tType 2: association present; difference up to second order")
    }
    
    pander(Result, style = "grid")
    pander("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1" )
    # cat("\n\n")
    # pandoc.emphasis(Message)
  } )
  invisible(x)
}

# Deciding the differential type based on order 0, 1, and 2 p.value
decide.type <- function(
  Active.P.value, Marginal.P.value, SecondOrd.P.value, alpha
)
{
  if(Active.P.value <= alpha && 
     Marginal.P.value > alpha && 
     SecondOrd.P.value > alpha){

    Type <- 0
  } else if(Active.P.value <= alpha && 
            Marginal.P.value <= alpha && 
            SecondOrd.P.value > alpha){
    Type <- 1
  } else if(Active.P.value <= alpha && 
            SecondOrd.P.value <= alpha){
    Type <- 2

  } else {
    Type <- NULL
  }
  
  return(Type)
}

## Assigning significance stars to p.value
AddStars <- function(P.val)
{
  if(P.val > 0.1){
    stars = " "
  } else if(P.val > 0.05 && P.val <= 0.1 ){
    stars = "."
  } else if(P.val > 0.01 && P.val <= 0.05){
    stars = "*"
  } else if(P.val > 0.001  && P.val <= 0.01){
    stars = "**"
  } else {
    stars = "***"
  }
  P.val <- signif(P.val, digits = 5)
  P.val <- paste0(P.val, " ", stars) 
  return(P.val)
}


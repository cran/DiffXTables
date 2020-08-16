# heterogeneity.R -- heterogeneity chi square test
# Created by: Ruby Sharma
# Created: January 7, 2019
# Modified: Feb 5, 2016. MS

#' @export
heterogeneity.test <- function(tables)
{

   if(mode(tables)!="list" || length(tables)<2 )
   {
      stop("only accept list of 2 or more matrices as input!")
   }

   for(k in 1:(length(tables)-1)){
      if(!identical(dim(tables[[k]]), dim(tables[[k+1]])))
        stop("All matries should have same number of rows and columns")
   }

   Hp <- Reduce('+', tables)
   if(all(Hp == 0))
   {
     Stat <- 0
      df <- 0
   }else{
      xall <-  unlist(lapply(seq_along(tables),  function(k){

               return(chisq.test.stat(tables[[k]]))
       }))
       xHp <- chisq.test.stat(Hp)

       # heterogeneity chi square, summation of chisq - pooled chisq
       Stat <- sum(xall) -xHp
       df <- (length(tables)* prod(dim(Hp) - 1)) - prod(dim(Hp) - 1)
  }

  p.val = pchisq(abs( Stat), df, lower.tail=F)
  DNAME <- deparse(substitute(tables))
  
  names( Stat) <- "X-squared"
  names(df) <- "df"
  
  return(structure(list(
    statistic = abs(Stat), 
    parameter = df, 
    p.value = p.val,
    data.name = DNAME,
    method = "Heterogeneity Test for Pattern Difference"),
    class = "htest"))

}



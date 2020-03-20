\name{simulate_diff_tables}
\alias{simulate_diff_tables}
\title{
  Simulating Contingency Tables that are Differential in Joint Distribution
}
\description{
  Generate contingency tables that are first-order, second-order or full-order differential in the joint distribution of row and column variables.
}

\usage{
simulate_diff_tables(
  K = 2, nrow = 3, ncol = 3, n = 100, B = 100, 
  type = c("second-order", "first-order", "full-order")
)
}

\arguments{
  \item{K}{
    the number of tables that are differential. It must be an integer greater than one.
  }
  
  \item{nrow}{
    the number of rows for all tables to be generated. It must be an integer greater than one.
  }
  
  \item{ncol}{
    the number of columns for all tables to be generated. It must be an integer greater than one.
  }
  
  \item{n}{
    the sample size for each table to be generated. It must be a positive integer.
  }
  
  \item{B}{
    the number of iterations indicating the level of differentiality. It must be a positive integer. The greater the value, the stronger the differentiality across tables.
  }
  
  \item{type}{
    the type of differential tables to be generated. Options are \code{"first-order"}, \code{"second-order"} (default), and \code{"full-order"}. See Details.}
 
}

\details{
The function randomly generates contingency tables differential in the joint distribution of the row and column variables. Specifically, three types of differential contingency tables can be simulated:

First-order differential contingency tables only differ in row or column marginal distribution. Such tables are differential in joint distribution, but different from second-order differential tables.

Second-order differential contingency tables differ in joint distribution The difference is not attributed to row or column marginal distributions.

Full-order differential contingency tables are tables that are both first-order and second-order differential.

The simulation starts with randomly generated probability tables where row and column variables are independent. The probability tables are modified to \code{K} tables such that they represent specific distributions that strictly satisfy the type requirement. Finally, contingency tables are generated using multinomial distribution using these probability tables parameters and the required sample size. 
}

\value{
  A list containing the following components:
    \item{contingency.tables}{a list of \code{K} contingency tables that are differential in joint distribution according to the \code{type} argument. They contain non-negative integers following the multinomial distribution with probability parameters from \code{probability.tables}.}

  \item{probability.tables}{a list of \code{K} tables representing randomly generated differential joint probabilities that reflect the specified \code{type}.}

   \item{method}{a string that specifies the type of the differential tables.}

}

\author{
  Ruby Sharma and Joe Song
}

\seealso{
  Differential tables are simulated to evaluate the following tests for comparing contingency tables: 
  
  {The Sharma-Song test}{\code{\link[DiffXTables]{sharma.song.test}}}
  
  {The comparative chi-squared test}{\code{\link[DiffXTables]{cp.chisq.test}}}
  
  {The heterogeneity test}{\code{\link[DiffXTables]{heterogeneity.test}}}
  
  The \code{\link[FunChisq]{simulate_tables}} function in package \pkg{FunChisq} can generate a variety of tables.
}

\examples{
# Three first-order differential tables:
  res <- simulate_diff_tables(K = 3, nrow = 4, ncol = 3, n = 150, B = 200, type = "first-order")
  print(res)
  
# Two second-order differential tables:
  res <- simulate_diff_tables(K = 2, nrow = 2, ncol = 5, n = 100, B = 100, type = "second-order")
  print(res)
  
# Four full-order differential tables:
  res <- simulate_diff_tables(K = 4, nrow = 3, ncol = 4, n = 250, B = 200, type = "full-order")
  print(res)
}

\keyword{datagen}


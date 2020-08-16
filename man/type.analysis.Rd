\name{type.analysis}
\alias{type.analysis}
\title{Comprehensive Type Analysis for Difference Across Contingency Tables}

\description{
Four types (0, 1, 2, null) are assigned to a collection of contingency tables to categorize their differences.
}

\usage{
  type.analysis(tables, alpha = 0.05)
}

\arguments{
  \item{tables}{a list of at least two non-negative matrices or data frames representing contingency tables of the same dimensions.}
  
  \item{alpha}{a numerical value representing the significance level for all involved hypothesis tests. Default is 0.05.}
}

\details{
      The function determines whether differences across patterns underlying multiple input contingency tables are type 0, 1, 2, or null. The function calls \code{\link[DiffXTables]{strength.test}}, \code{\link[DiffXTables]{marginal.change.test}}, and \code{\link[DiffXTables]{sharma.song.test}} to obtain three p-values and use them to decide on the type of difference among input contingency tables \insertCite{sharma2020COPTYP}{DiffXTables}.

Type null: association absent from the input tables. These tables can still differ in joint distribution, but a strong relationship is lacking between row and column variables. No mechanisms are implied.

Type 0: association present and patterns are conserved. 
These tables show strong row and column association but have no difference in distribution. Conserved mechanisms with conserved trajectories are implied.

Type 1: association present; the difference is up to the first order. These tables show strong row and column association, differ in marginal distribution, and do not differ deviation from joint distribution to product of marginals. Conserved mechanisms with differential trajectories are implied. Differences in trajecotry can be due to changed stimuli.

Type 2: association present; the difference is up to the second order. These tables show strong row and column association and differ in deviation from joint distribution to the product of marginals. Differential mechanisms are implied.
}

\value{
  A list with class "\code{DiffXTableTypeAnalysis}" containing the following components:

\item{Zeroth.order}{a list with class "\code{htest}" containing the test statistic, degrees of freedom and p-value of the \code{strength.test}.}

\item{First.order}{a list with class "\code{htest}" containing the test statistic, degrees of freedom and p-value of the \code{marginal.change.test}.}

\item{Second.order}{a list with class "\code{htest}" containing the test statistic, degrees of freedom and p-value of the \code{sharma.song.test}.}

\item{Type}{the type of differences across the input contingency tables. Possible values are 0, 1, 2, and \code{NULL}.}
}

\references{
\insertAllCited{} 
}

\seealso{
\code{\link[DiffXTables]{strength.test}}, \code{\link[DiffXTables]{marginal.change.test}}, and \code{\link[DiffXTables]{sharma.song.test}}.

A related function is the comparative chi-squared test \code{\link[DiffXTables]{cp.chisq.test}}.
}

\examples{
  # Type-null tables:
  tables <- list(
   matrix(c(7, 4, 1,
            3, 6, 9,
            2, 4, 4), nrow=3),
   matrix(c(2, 1, 2,
            2, 0, 8,
            6, 2, 7), nrow=3)
  )
  type.analysis(tables, alpha = 0.05)

  # Type-0 tables:
  tables <- list(
   matrix(c(30,0,0,
            0,10,0,
            0,0,20), nrow=3),
   matrix(c(30,0,0,
            0,10,0,
            0,0,20), nrow=3)
  )
  type.analysis(tables, alpha = 0.05)
  
  # Type-1 differential tables:
  tables <- list(
   matrix(c(30,0,0,
            0,10,0,
            0,0,20), nrow=3),
   matrix(c(10,0,0,
            0,20,0,
            0,0,30), nrow=3)
  )
  type.analysis(tables, alpha = 0.05)
  
  # Type-2 differential tables:
  tables <- list(
    matrix(c(4,0,0,
             0,4,0,
             0,0,4), nrow=3),
    matrix(c(0,4,4,
             4,0,4,
             4,4,0), nrow=3)
  )
  type.analysis(tables, alpha = 0.05)

}

\author{
Ruby Sharma and Joe Song
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\keyword{htest}
\keyword{nonparametric}
\keyword{distribution}

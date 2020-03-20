\name{heterogeneity.test}
\alias{heterogeneity.test}

\title{Heterogeneity Test for Distributional Differences Across Contingency Tables}

\description{
Across given contingency tables, the test admits any type of differences in either the joint or marginal distributions of the tables.}

\usage{
  heterogeneity.test(tables)
}

\arguments{
  \item{tables}{  a list of at least two non-negative matrices or data frames representing contingency tables of the same dimensions.}
}

\details{
    The heterogeneity test determines whether the patterns underlying multiple contingency tables are heterogeneous or differential. The chi-squared distribution is used for the null distribution of its test statistic (Zar, 2010).
}

\value{
  A list with class "\code{htest}" containing the following components:

\item{statistic}{heterogeneity test statistic.}
\item{parameter}{degrees of freedom of used for the null distribution of the heterogeneity test statistic.}
\item{p.value}{p-value of the heterogeneity test, computed using the chi-squared distribution.}
}

\references{
Zar, J. H. (2010) \emph{Biostatistical Analysis}, 5th Ed., New Jersey: Prentice Hall.
}

\seealso{
The comparative chi-squared test \code{\link[DiffXTables]{cp.chisq.test}}.

The Sharma-Song test \code{\link[DiffXTables]{sharma.song.test}}.
}

\examples{
  # Two second-order differential tables:
  tables <- list(
    matrix(c(4,0,0,
             0,4,0,
             0,0,4), nrow=3),
    matrix(c(0,4,4,
             4,0,4,
             4,4,0), nrow=3)
  )
  heterogeneity.test(tables)
  
  # Three tables differ in the first-order but not second-order:
  tables <- list(
    matrix(c(2, 4,  6,  8, 
             3, 6,  9, 12, 
             4, 8, 12, 16), nrow=4),
    matrix(c( 2, 1,  3,  7,
              2, 1,  3,  7,
             10, 5, 15, 35), nrow=4),
    matrix(c(40, 16, 72, 16, 
             45, 18, 81, 18,
             25, 10, 45, 10), nrow=4)
  )
  heterogeneity.test(tables)
}

\keyword{htest}
\keyword{nonparametric}

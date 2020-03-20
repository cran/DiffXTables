\name{sharma.song.test}
\alias{sharma.song.test}
\title{Sharma-Song Test for Second-Order Distributional Differences Across Contingency Tables}

\description{
Across given contingency tables, the test admits only second-order difference in the joint distributions underlying the tables.
}

\usage{
  sharma.song.test(tables)
}

\arguments{
  \item{tables}{a list of at least two non-negative matrices or data frames representing contingency tables of the same dimensions.}
}

\details{
      The Sharma-Song test determines whether the patterns underlying multiple input contingency tables are second-order differential. Its null test statistic is proved to asymptotically follow the chi-squared distribution.
}

\value{
  A list with class "\code{htest}" containing the following components:

\item{statistic}{the Sharma-Song chi-squared test statistic.}
\item{parameter}{degrees of freedom of the chi-squared test statistic.}
\item{p.value}{p-value of the Sharma-Song test, computed using the chi-squared distribution.}
}

\seealso{
The comparative chi-squared test \code{\link[DiffXTables]{cp.chisq.test}}.

The heterogeneity test \code{\link[DiffXTables]{heterogeneity.test}}.
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
  sharma.song.test(tables)
  
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
  sharma.song.test(tables)
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

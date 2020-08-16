\name{marginal.change.test}
\alias{marginal.change.test}
\title{Test for Marginal Change Across Contingency Tables}

\description{
The test detects change in either row or column marginal distributions across given contingency tables.
}

\usage{
  marginal.change.test(tables)
}

\arguments{
  \item{tables}{a list of at least two non-negative matrices or data frames representing contingency tables of the same dimensions.}
}

\details{
      The marginal change test determines whether the patterns underlying multiple input contingency tables have changed row or column marginal distributions. Its test statistic is proved to asymptotically follow the chi-squared distribution under the null hypothesis of same row and marginal distributions across tables \insertCite{sharma2020COPTYP}{DiffXTables}.
      
      The test statistic is minimized to zero if and only if observed row marginal distributions are the same across tables and so do the column marginal distributions.
}

\value{
  A list with class "\code{htest}" containing the following components:

\item{statistic}{the chi-squared test statistic.}
\item{parameter}{the degrees of freedom of the null chi-squared distribution.}
\item{p.value}{the p-value for the test, computed using the chi-squared distribution.}
}


\references{
\insertAllCited{} 
}


\seealso{
\code{\link[DiffXTables]{sharma.song.test}}, \code{\link[DiffXTables]{strength.test}}, and
\code{\link[DiffXTables]{type.analysis}}.
}

\examples{
  # Two first-order differential tables:
  tables <- list(
   matrix(c(30,0,0,
            0,10,0,
            0,0,20), nrow=3),
   matrix(c(10,0,0,
            0,20,0,
            0,0,30), nrow=3)
  )
  marginal.change.test(tables)
  
  # Tables differ in the second-order but not first-order:
  tables <- list(
    matrix(c(4,0,0,
             0,4,0,
             0,0,4), nrow=3),
    matrix(c(0,0,4,
             0,4,0,
             4,0,0), nrow=3)
  )
   marginal.change.test(tables)
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

\name{strength.test}
\alias{strength.test}
\title{Strength Test for Association in Multiple Contingency Tables}

\description{
The test determines the total strength of association in multiple contingency tables.
}

\usage{
  strength.test(tables)
}

\arguments{
  \item{tables}{a list of at least two non-negative matrices or data frames representing contingency tables.}
}

\details{
      The strength test determines total amount of association in multiple input contingency tables. Its test statistic asymptotically follows the chi-squared distribution under the null hypothesis of each table having independent row and column variables \insertCite{sharma2020COPTYP}{DiffXTables}.
      
      The test statistic is minimized to zero if and only if row and column variables are empirically independent of each other in every table.
      
      This test is considered a zeroth-order test in the function \code{type.analysis} that characterizes the difference across multiple contingency tables.
}

\value{
  A list with class "\code{htest}" containing the following components:

\item{statistic}{the strength test statistic.}
\item{parameter}{the degrees of freedom of null chi-squared distribution.}
\item{p.value}{the p-value for the test, computed using the null chi-squared distribution.}
}

\references{
\insertAllCited{} 
}


\seealso{
A second-order different test \code{\link[DiffXTables]{sharma.song.test}}. A first-order difference test \code{\link[DiffXTables]{marginal.change.test}}. The comparative type of differences across contingency tables \code{\link[DiffXTables]{type.analysis}}.
}

\examples{
  # Both tables have strong association:
  tables <- list(
   matrix(c(30,0,0,
            0,10,0,
            0,0,20), nrow=3),
   matrix(c(10,0,0,
            0,20,0,
            0,0,30), nrow=3)
  )
  strength.test(tables)
  
  # One table has strong association:
  tables <- list(
    matrix(c(4,0,0,
             0,4,0,
             0,0,4), nrow=3),
    matrix(c(4,0,4,
             8,4,8,
             4,0,4), nrow=3)
  )
  strength.test(tables)
  
  # Both tables has no association:
  tables <- list(
    matrix(c(4,0,4,
             8,4,8,
             4,0,4), nrow=3),
    matrix(c(4,0,4,
             8,4,8,
             4,0,4), nrow=3)
  )
  strength.test(tables)
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

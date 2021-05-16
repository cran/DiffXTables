\name{sharma.song.test}
\alias{sharma.song.test}
\title{Sharma-Song Test for Second-Order Difference Across Contingency Tables}

\description{
The test detects differential departure from independence via second-order difference in joint distributions underlying two or more contingency tables.
}

\usage{
sharma.song.test(
  tables, null.table.marginal = c("observed", "uniform"), 
  compensated = FALSE
)
}

\arguments{
  \item{tables}{a list of at least two non-negative matrices or data frames representing contingency tables of the same dimensions.}
  \item{null.table.marginal}{a character string to specify marginal distributions of null tables. The options are \code{"observed"}(default) and \code{"uniform"}.}
  \item{compensated}{a logical value to compensate for the Cochran's condition. It is only used if \code{null.table.marginal="observed"}. Default is FALSE.}
}

\details{
      The Sharma-Song test determines whether the patterns underlying multiple input contingency tables are second-order differential. The test statistic measures differential departure form independence. Its null test statistic is proved to asymptotically follow the chi-squared distribution. For full detail of the test, see \insertCite{ruby2021/BIOINF-2020-0795}{DiffXTables}.

      If \code{null.table.marginal} is set to \code{"observed"}, the null hypothesis uses the observed marginals. The \code{compensated} parameter, if set to \code{TRUE}, adds a small constant to each entry of the tables to address the Cochran's condition that the expected count in any table entry is 5 or less. %It can detect all differential patterns.

      If the \code{null.table.marginal} is set to \code{"uniform"}, the null tables are set to have uniform marginals. No longer suffering from the Cochran's condition, it detects second-order differential patterns with additional robustness. 
}

\value{
  A list with class "\code{htest}" containing the following components:

\item{statistic}{the Sharma-Song chi-squared test statistic.}
\item{parameter}{degrees of freedom of the chi-squared test statistic.}
\item{p.value}{p-value of the Sharma-Song test, computed using the chi-squared distribution.}
}


\references{
\insertAllCited{} 
}


\seealso{
\code{\link[DiffXTables]{cp.chisq.test}}, \code{\link[DiffXTables]{heterogeneity.test}},  \code{\link[DiffXTables]{strength.test}}, \code{\link[DiffXTables]{marginal.change.test}}, and
\code{\link[DiffXTables]{type.analysis}}.
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
  sharma.song.test(tables, null.table.marginal = "uniform")
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

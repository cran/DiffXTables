\name{cp.chisq.test}
\alias{cp.chisq.test}
\title{
Comparative Chi-Squared Test for Distributional Differences Across Contingency Tables
}
\description{
Across given contingency tables, the test admits any type of differences in either the joint or marginal distributions underlying the tables.
}

\usage{
cp.chisq.test(
  tables, method=c("chisq", "nchisq", "default", "normalized"),
  log.p = FALSE
)
}

\arguments{
  \item{tables}{
  a list of at least two nonnegative matrices or data frames representing contingency tables of the same dimensions.
}

  \item{method}{
  a character string to specify the method to compute the chi-squared statistic and its p-value. The default is \code{"chisq"}. See Details.

    Note: \code{"default"} and \code{"normalized"} are deprecated.
}
  \item{log.p}{
  logical; if \code{TRUE}, the p-value is given as \code{log(p)}. Taking the log improves the accuracy when p-value is close to zero. The default is \code{FALSE}.
  }
}

\details{
The comparative chi-squared test determines whether the patterns underlying multiple contingency tables are heterogeneous. Its null test statistic is proved to asymptotically follow the chi-squared distribution (Song et al., 2014; Zhang et al., 2015). This test is different from the heterogeneity test (Zar, 2010).

Two methods are provided to compute the chi-squared statistic and its p-value. When \code{method = "chisq"} (or \code{"default"}), the p-value is computed using the chi-squared distribution; when \code{method =} \code{"nchisq"} (or \code{"normalized"}) a normalized statistic is obtained by shifting and scaling the original chi-squared test statistic and a p-value is computed using the standard normal distribution (Box et al., 2005). The normalized test is more conservative on the degrees of freedom.
}

\value{
A list with class "\code{htest}" containing the following components:

\item{statistic}{chi-squared test statistic if \code{method = "chisq"} (equivalent to \code{"default"}), or normalized test statistic if \code{method = "nchisq"} (equivalent to \code{"normalized"}).}
\item{parameter}{degrees of freedom of the chi-squared statistic.}
\item{p.value}{p-value of the comparative chi-squared test. By default, it is computed by the chi-squared distribution (\code{method = "chisq"} or \code{"default"}). If \code{method = "nchisq"} (or \code{"normalized"}), it is the p-value of the normalized chi-squared statistic using the standard normal distribution.}
}

\references{
Box, G. E., Hunter, J. S. and Hunter, W. G. (2005) \emph{Statistics for Experimenters: Design, Innovation and Discovery}, 2nd Ed., New York: Wiley-Interscience.

Song M., Zhang Y., Katzaroff A. J., Edgar B. A. and Buttitta L. (2014) Hunting complex differential gene interaction patterns across molecular contexts. \emph{Nucleic Acids Research} \bold{42}(7), e57. Retrieved from \url{https://doi.org/10.1093/nar/gku086}

Zar, J. H. (2010) \emph{Biostatistical Analysis}, 5th Ed., New Jersey: Prentice Hall.

Zhang, Y., Liu, Z. L. and Song, M. (2015) ChiNet uncovers rewired transcription subnetworks in tolerant yeast for advanced biofuels conversion. \emph{Nucleic Acids Research} \bold{43}(9), 4393--4407. Retrieved from \url{https://doi.org/10.1093/nar/gkv358}
}

\author{
Joe Song
}

\seealso{
The Sharma-Song test \code{\link[DiffXTables]{sharma.song.test}}.

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
  cp.chisq.test(tables)
  
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
  cp.chisq.test(tables)
}

\keyword{htest}
\keyword{nonparametric}

\name{pneumo}
\alias{pneumo}
\docType{data}
\title{Pneumoconiosis in Coalminers Data}
\description{
The \code{pneumo} data frame has 8 rows and 4 columns.
Exposure time is explanatory, and there are 3 ordinal response variables.

}
\usage{data(pneumo)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{exposure.time}{a numeric vector, in years}
    \item{normal}{a numeric vector, counts}
    \item{mild}{a numeric vector, counts}
    \item{severe}{a numeric vector, counts}
  }
}
\details{
These were collected from coalface workers. In the original
data set, the two most severe categories were combined.


}
\source{

 Ashford, J.R., 1959. An approach to the analysis of data for
 semi-quantal responses in biological assay.
 \emph{Biometrics}, \bold{15}, 573--581.


}
\seealso{
  \code{\link{cumulative}}.


}

\references{

  McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}
\examples{
# Fit the proportional odds model, p.179, in McCullagh and Nelder (1989)
pneumo <- transform(pneumo, let = log(exposure.time))
vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo)
}
\keyword{datasets}

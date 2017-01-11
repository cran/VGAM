\name{flourbeetle}
\alias{flourbeetle}
\docType{data}
\title{Mortality of Flour Beetles from Carbon Disulphide}
\description{
The \code{flourbeetle} data frame has 8 rows and 4 columns.
Two columns are explanatory, the other two are responses.


}
\usage{data(flourbeetle)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{logdose}{\code{\link[base]{log10}} applied to \code{CS2mgL}. }
    \item{CS2mgL}{a numeric vector, the concentration of gaseous
                  carbon disulphide in mg per litre. }
    \item{exposed}{a numeric vector, counts; the number of
                   beetles exposed to the poison. }
    \item{killed}{a numeric vector, counts; the numbers killed. }
  }
}
\details{
These data were originally given in Table IV of Bliss (1935) and
are the combination of
two series of toxicological experiments involving
\emph{Tribolium confusum}, also known as the flour beetle.
Groups of such adult beetles were exposed for 5 hours of
gaseous carbon disulphide at different concentrations,
and their mortality measured.



}
\source{

 Bliss, C.I., 1935.
 The calculation of the dosage-mortality curve.
 \emph{Annals of Applied Biology}, \bold{22}, 134--167.


}
\seealso{
  \code{\link{binomialff}},
  \code{\link{probit}}.


}

%\references{
%
%
%
%
%
%}

\examples{
fit1 <- vglm(cbind(killed, exposed - killed) ~ logdose, binomialff(link = probit),
             data = flourbeetle, trace = TRUE)
summary(fit1)
}
\keyword{datasets}

\name{coalminers}
\alias{coalminers}
\docType{data}
\title{ Breathlessness and Wheeze Amongst Coalminers Data}
\description{
  Coalminers who are smokers without radiological pneumoconiosis,
  classified by age, breathlessness and wheeze.
}
\usage{data(coalminers)}
\format{
  A data frame with 9 age groups with the following 5 columns.
  \describe{
    \item{BW}{Counts with breathlessness and wheeze. }
    \item{BnW}{Counts with breathlessness but no wheeze. }
    \item{nBW}{Counts with no breathlessness but wheeze. }
    \item{nBnW}{Counts with neither breathlessness or wheeze. }
    \item{age}{Age of the coal miners (actually, the
      midpoints of the 5-year category ranges). }
  }
}
\details{
  The data were published in Ashford and Sowden (1970).
  A more recent analysis is McCullagh and Nelder (1989, Section 6.6).



}
\source{
  Ashford, J. R. and Sowden, R. R. (1970)
  Multi-variate probit analysis.
  \emph{Biometrics}, \bold{26}, 535--546.

}
\references{

  McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}. 2nd ed. London: Chapman & Hall.



}
\examples{
str(coalminers)
}
\keyword{datasets}

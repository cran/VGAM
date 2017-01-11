\name{cfibrosis}
\alias{cfibrosis}
\docType{data}
\title{ Cystic Fibrosis Data
%%   ~~ data name/kind ... ~~

}
\description{
  This data frame concerns families data
  and cystic fibrosis.


}
\usage{
data(cfibrosis)
}

\format{
  A data frame with 24 rows on the following 4 variables.

\describe{
  \item{siblings, affected, ascertained, families}{
    Over ascertained families, the \eqn{k}th ascertained family
    has \eqn{s_k} siblings of whom \eqn{r_k}
    are affected and \eqn{a_k} are ascertained.

  }
}

}
\details{

  The data set allows a classical segregation analysis
  to be peformed. In particular,
  to test Mendelian segregation ratios in nuclear family data.
The likelihood has similarities with \code{\link{seq2binomial}}.



%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{

  The data is originally from Crow (1965) and
  appears as Table 2.3 of Lange (2002).


Crow, J. F. (1965)
Problems of ascertainment in the analysis of family data.
Epidemiology and Genetics of Chronic Disease.
Public Health Service Publication 1163,
Neel J. V., Shaw M. W., Schull W. J., editors,
Department of Health, Education, and Welfare, Washington, DC,
USA.


Lange, K. (2002)
Mathematical and Statistical Methods for Genetic Analysis.
Second Edition.
Springer-Verlag: New York, USA.


}
\examples{
cfibrosis
summary(cfibrosis)
}
\keyword{datasets}




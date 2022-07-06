\name{lpossums}
\alias{lpossums}
\docType{data}
\title{ Leadbeater's Possums }

\description{
  Abundance of Leadbeater's Possums observed in the field.


}
\usage{
data(lpossums)
}
\format{
  A data frame with the following variables.

  \describe{

    \item{number}{
      Values between 0 and 10 excluding 6.


    }
    \item{ofreq}{
      Observed frequency, i.e., the number of sites.

    }

  }
}
\details{
  A small data set recording the abundance of Leadbeater's Possums
  \emph{Gymnobelideus leadbeateri} observed in
  the montane ash forests of 
  the Central Highlands of Victoria,
  in south-eastern Australia.
  There are 151 3-hectare sites.
  The data has more 0s than usual relative to the Poisson,
  as well as exhibiting overdispersion too.
  

}
\source{

Welsh, A. H., Cunningham, R. B., Donnelly, C. F. and
Lindenmayer, D. B. (1996).
Modelling the abundances of rare species: statistical models
for counts with extra zeros.
\emph{Ecological Modelling},
\bold{88},
297--308.


}
\seealso{
  \code{\link[VGAM]{zipoissonff}}.


}
\examples{
lpossums
(samplemean <- with(lpossums, weighted.mean(number, ofreq)))
with(lpossums,  var(rep(number, times = ofreq)) / samplemean)
sum(with(lpossums, ofreq))
\dontrun{ spikeplot(with(lpossums, rep(number, times = ofreq)),
  main = "Leadbeater's possums", col = "blue", xlab = "Number") }
}
\keyword{datasets}


%
%

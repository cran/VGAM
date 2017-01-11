\name{corbet}
\alias{corbet}
\docType{data}
\title{
  Corbet's Butterfly Data

%%   ~~ data name/kind ... ~~
}
\description{
  About 3300 individual butterflies were caught in Malaya
  by naturalist Corbet trapping butterflies.
  They were classified to about 500 species.

%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(corbet)}
\format{
  A data frame with 24 observations on the following 2 variables.
\describe{
    \item{\code{species}}{Number of species. }
    \item{\code{ofreq}}{Observed frequency of individual
      butterflies of that species. }
  }
}
%%\format{
%%  The format is:
%% chr "corbet"
%%}
\details{
  In the early 1940s Corbet spent two years trapping butterflies
  in Malaya. Of interest was the total number of species.
  Some species were so rare (e.g., 118 species had only
  one specimen) that it was thought likely that there were
  many unknown species.


%%  ~~ If necessary, more details than the __description__ above
}


%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
\references{

  Fisher, R. A., Corbet, A. S. and Williams, C. B. (1943)
  The Relation Between the Number of Species and
     the Number of Individuals in a Random Sample of an Animal
     Population.
\emph{Journal of Animal Ecology},
\bold{12}, 42--58.

}
\examples{
summary(corbet)
}
\keyword{datasets}

\name{budworm}
\alias{budworm}
\docType{data}
\title{ Western Spuce Budworm }

\description{
Counts of western spuce budworm 
(Choristoneura freemani) across seven developmental
stages
(five larval instars, pupae, and adults)
on 12 sampling occasions.


}
\usage{
data(budworm)
}
\format{
  A data frame with the following variables.

  \describe{

    \item{ddeg}{
      Degree days.


    }
    \item{total}{
      Sum of stages 1--7.

    }
    \item{stage1, stage2, stage3, stage4}{
      Successive stages.

    }
    \item{stage5, stage6, stage7 }{
      Successive stages.

    }

  }
}
\details{
  This data concerns the development of a defoliating
  moth widespread in western North America
  (i.e., north of Mexico).
  According to Boersch-Supan (2021),
  the insect passes hrough successive stages
  \eqn{j=1,\ldots,r}, delimited by \eqn{r-1} moults.
  The data was originally used in a 1986 publication
  but has been corrected for two sampling occasions;
  the data appears in Candy (1990) and
  was analyzed in Boersch-Supan (2021).
  See the latter for more references.



}
\source{


Candy, S. G. (1990).
\emph{Biology of the mountain pinhole borer,
Platypus subgranosus Scheld, in Tasmania}.
MA thesis, University of Tasmania, Australia.
\code{https://eprints.utas.edu.au/18864/}.


}
\references{


Boersch-Supan, P. H. (2021).
Modeling insect phenology using ordinal
regression and continuation ratio models.
\emph{ReScience C},
\bold{7.1}, 1--14.



%\bold{7.1}(#5), 1--14.



}
%\seealso{
%  \code{\link[VGAM]{zipoissonff}}.


%}
\examples{
budworm
summary(budworm)
}
\keyword{datasets}


%
%

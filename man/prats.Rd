\name{prats}
\alias{prats}
\docType{data}
\title{ Pregnant Rats Toxological Experiment Data }

\description{
  A small toxological experiment data.
  The subjects are fetuses from
  two randomized groups of pregnant rats,
  and they were given a placebo or chemical treatment.
  The number with birth defects were recorded, as well
  as each litter size.


}
\usage{
data(prats)
}
\format{
  A data frame with the following variables.

  \describe{

    \item{treatment}{
    A \code{0} means control;
    a \code{1} means the chemical treatment.


    }
    \item{alive, litter.size}{
    The number of fetuses alive at 21 days, out of
    the number of fetuses alive at  4 days (the litter size).


    }

  }
}
\details{
  The data concerns a toxological experiment where
  the subjects are fetuses from
  two randomized groups of 16 pregnant rats each,
  and they were given a placebo or chemical treatment.
  The number with birth defects andn the litter size were recorded.
  Half the rats were fed a control diet during pregnancy and
  lactation, and the diet of the other half was treated with a
  chemical. For each litter the number of pups alive at 4 days
  and the number of pups that survived the 21 day lactation period,
  were recorded.


}
\source{

  Weil, C. S. (1970)
  Selection of the valid number of sampling units and a consideration
  of their combination in toxicological studies involving
  reproduction, teratogenesis or carcinogenesis.
  \emph{Food and Cosmetics Toxicology},
  \bold{8}(2), 177--182.


%Food and Cosmetics Toxicology
%Fd. Cosmet. Toxicol.



}
\references{

  Williams, D. A. (1975)
  The Analysis of Binary Responses From Toxicological
            Experiments Involving Reproduction and Teratogenicity.
  \emph{Biometrics},
  \bold{31}(4), 949--952.


}
\seealso{
  \code{\link[VGAM]{betabinomial}},
  \code{\link[VGAM]{betabinomialff}}.


}
\examples{
prats
colSums(subset(prats, treatment == 0))
colSums(subset(prats, treatment == 1))
summary(prats)
}
\keyword{datasets}


%
%

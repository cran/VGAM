\name{finney44}
\alias{finney44}
\docType{data}
\title{ Toxicity trial for insects
%%   ~~ data name/kind ... ~~
}
\description{
  A data frame of a toxicity trial.

%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(finney44)}
\format{
  A data frame with 6 observations on the following 3 variables.
  \describe{
    \item{\code{pconc}}{a numeric vector, percent concentration
of pyrethrins. }
    \item{\code{hatched}}{number of eggs that hatched. }
    \item{\code{unhatched}}{number of eggs that did not hatch. }
  }
}
\details{

Finney (1944) describes a toxicity trial of five different
concentrations of pyrethrins (percent) plus a control that were
administered to eggs of \emph{Ephestia kuhniella}.
The natural mortality rate is large, and a common adjustment
is to use Abbott's formula.

}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{

  Finney, D. J., 1944.
  The application of the probit method to toxicity test
  data adjusted for mortality in the controls.
  \emph{Annals of Applied Biology}, \bold{31}, 68--74.


  Abbott, W. S. (1925).
  A method of computing the effectiveness of an insecticide.
  \emph{Journal of Economic Entomology}, 18, 265--7.


%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(finney44)
transform(finney44, mortality = unhatched / (hatched + unhatched))
}
\keyword{datasets}

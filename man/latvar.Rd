\name{latvar}
\alias{lv}
\alias{latvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Latent Variables }
\description{
  Generic function for the \emph{latent variables} of a model.


}
\usage{
latvar(object, ...)
    lv(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An object for which the extraction of latent
  variables is meaningful.


  }
  \item{\dots}{
  Other arguments fed into the specific
  methods function of the model. Sometimes they are fed
  into the methods function for \code{\link{Coef}}.


  }
}
\details{
  Latent variables occur in reduced-rank regression models,
  as well as in quadratic and additive ordination models.
  For the latter two,
  latent variable values are often called \emph{site scores}
  by ecologists.
  Latent variables are linear combinations of the explanatory
  variables.


}
\value{
  The value returned depends specifically on the methods
  function invoked.


}
\references{
Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.


}
\author{ Thomas W. Yee }

\section{Warning}{
  \code{\link{latvar}} and \code{\link{lv}} are identical,
  but the latter will be deprecated soon.


  Latent variables are not really applicable to
  \code{\link{vglm}}/\code{\link{vgam}} models.


}


\seealso{
  \code{latvar.qrrvglm},
  \code{latvar.rrvglm},
  \code{latvar.cao},
  \code{\link{lvplot}}.


}

\examples{
\dontrun{
hspider[, 1:6] <- scale(hspider[, 1:6])  # Standardized environmental vars
set.seed(123)
p1 <- cao(cbind(Pardlugu, Pardmont, Pardnigr, Pardpull, Zoraspin) ~
          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
          family = poissonff, data = hspider, Rank = 1, df1.nl =
          c(Zoraspin = 2.5, 3), Bestof = 3, Crow1positive = TRUE)

var(latvar(p1))  # Scaled to unit variance  # Scaled to unit variance
c(latvar(p1))    # Estimated site scores
}
}
\keyword{models}
\keyword{regression}


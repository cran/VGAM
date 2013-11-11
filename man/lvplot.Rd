\name{lvplot}
\alias{lvplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Latent Variable Plot }
\description{
  Generic function for a \emph{latent variable plot}
  (also known as an \emph{ordination diagram} by ecologists).
}
\usage{
lvplot(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object for a latent
    variable plot is meaningful.
  }
  \item{\dots}{ Other arguments fed into the specific
    methods function of the model. They usually are graphical
    parameters, and sometimes they are fed
    into the methods function for \code{\link{Coef}}.
  }
}
\details{
  Latent variables occur in reduced-rank regression models,
  as well as in
  quadratic and additive ordination. For the latter,
  latent variables are often called the \emph{site scores}.
  Latent variable plots were coined by Yee (2004), and have
  the latent variable as at least one of its axes.


}
\value{
  The value returned depends specifically on the methods
  function invoked.


}
\references{
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

\note{
  Latent variables are not really applicable to
  \code{\link{vglm}}/\code{\link{vgam}} models.
}


\seealso{
  \code{\link{lvplot.qrrvglm}},
  \code{lvplot.cao},
  \code{\link{latvar}},
  \code{\link{trplot}}.


}

\examples{
\dontrun{
hspider[,1:6] <- scale(hspider[,1:6])  # Standardized environmental vars
set.seed(123)
p1 <- cao(cbind(Pardlugu, Pardmont, Pardnigr, Pardpull, Zoraspin) ~
          WaterCon + BareSand + FallTwig +
          CoveMoss + CoveHerb + ReflLux,
          family = poissonff, data = hspider, Bestof = 3,
          df1.nl = c(Zoraspin = 2.5, 3), Crow1positive = TRUE)
index <- 1:ncol(depvar(p1))
lvplot(p1, lcol = index, pcol = index, y = TRUE, las = 1)
}
}
\keyword{models}
\keyword{regression}


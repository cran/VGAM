\name{simplex}
\alias{simplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simplex distribution }
\description{
  The two parameters of the univariate simplex distribution are estimated.
}
\usage{
simplex(lmu = "logit", lsigma = "loge", emu=list(), esigma=list(),
        imu = NULL, isigma = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu}{
  Link function for \code{mu}. 
  See \code{\link{Links}} for more choices.
  The parameter lies in the unit interval.

  }
  \item{lsigma}{
  Link function for \code{sigma}. 
  See \code{\link{Links}} for more choices.
  The parameter is positive, therefore the log link is the default.

  }
  \item{emu, esigma}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{imu, isigma}{
  Optional initial values for \code{mu} and \code{sigma}.
  A \code{NULL} means a value is obtained internally.

  }
}
\details{
  See Jorgensen (1997) for details. 

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
\references{
Jorgensen, B. (1997)
\emph{The Theory of Dispersion Models}.
London: Chapman & Hall
}
\author{ T. W. Yee }
\note{ 
  This distribution is potentially useful for dispersion modelling.
  This family function only works for intercept-only models, i.e.,
  the formula should have \code{~ 1}.

}

\seealso{ 
  \code{\link{rig}}.
}
\examples{
x = runif(n <- 100)
y = rbeta(n, shape1=3+x, shape2=4-x)
(fit = vglm(y ~ 1, simplex, trace=TRUE, cri="coef"))
coef(fit, matrix=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}



% zz fitted values and formulas needed here





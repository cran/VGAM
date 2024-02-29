\name{gensh}
\alias{gensh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Secant Hyperbolic Regression
  Family Function }
\description{
  Estimation of the parameters of the
  generalized secant hyperbolic
  distribution.


}
\usage{
gensh(shape, llocation = "identitylink",
      lscale = "loglink", zero = "scale",
      ilocation = NULL, iscale = NULL, imethod = 1,
      glocation.mux = exp((-4:4)/2),
      gscale.mux = exp((-4:4)/2),
      probs.y = 0.3, tol0 = 1e-4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shape}{Numeric of length 1.
    Shape parameter, called \eqn{t} in
    Vaughan (2002).
    Valid values are
    \eqn{-\pi/2 < t}{-0.5 * pi < t}.
  
 
  }
  \item{llocation, lscale}{
    Parameter link functions applied to the
    two parameters.
    See \code{\link{Links}} for more choices.
    See \code{\link{CommonVGAMffArguments}}
    for more information.

  }
  \item{zero, imethod}{
    See \code{\link{CommonVGAMffArguments}}
    for information.

  }
  \item{ilocation, iscale}{
    See \code{\link{CommonVGAMffArguments}}
    for information.

  }
  \item{glocation.mux, gscale.mux}{
    See \code{\link{CommonVGAMffArguments}}
    for information.

  }

  \item{probs.y, tol0}{
    See \code{\link{CommonVGAMffArguments}}
    for information.

  }

}
\details{
  The probability density function of the
  hyperbolic secant distribution
  is given by
  \deqn{f(y; a, b, s) =
[(c_1 / b) \; \exp(c_2 z)] / [
\exp(2 c_2 z) + 2 C_3 \exp(c_2 z) + 1]}{%
((c1 / b) exp(c2 z)) / (
exp(2 c2 z) + 2 C3 exp(c2 z) + 1)}
  for shape
  parameter \eqn{-\pi < s}{-pi < s}
  and all real \eqn{y}.
  The scalars \eqn{c_1}{c1}, \eqn{c_2}{c2},
  \eqn{C_3}{C3} are functions of \eqn{s}.
  The mean of \eqn{Y} is
  the location parameter \eqn{a}
  (returned as the fitted values).
  All moments of the distribution are finite.


  
  Further details about
  the parameterization can be found
  in Vaughan (2002).
  Fisher scoring is implemented and it has
  a diagonal EIM.
  More details are at
  \code{\link{Gensh}}.



}

\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

  Vaughan, D. C. (2002).
  The generalized secant hyperbolic
  distribution and its properties.
  \emph{Communications in Statistics---Theory
    and Methods},
  \bold{31}(2): 219--238.


}

\author{ T. W. Yee }
%\note{

%}
\seealso{
  \code{\link{hypersecant}},
  \code{\link{logistic}}.



% \code{\link{nefghs}},


}
\examples{sh <- -pi / 2; loc <- 2
hdata <- data.frame(x2 = rnorm(nn <- 200))
hdata <- transform(hdata, y = rgensh(nn, sh, loc))
fit <- vglm(y ~ x2, gensh(sh), hdata, trace = TRUE)
coef(fit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}

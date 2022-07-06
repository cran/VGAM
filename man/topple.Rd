\name{topple}
\alias{topple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Topp-Leone Distribution Family Function }
\description{
 Estimating the parameter of the Topp-Leone distribution by
 maximum likelihood estimation.

}
\usage{
topple(lshape = "logitlink", zero = NULL, gshape = ppoints(8),
       parallel = FALSE,
       type.fitted = c("mean", "percentiles", "Qlink"),
       percentiles = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, gshape}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{zero, parallel}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{type.fitted, percentiles}{
  See \code{\link{CommonVGAMffArguments}} for information.
  Using \code{"Qlink"} is for quantile-links in \pkg{VGAMextra}.


  }
}
\details{
  The Topple distribution
  has a probability density function that can be written
  \deqn{f(y;s) = 2 s (1 - y) [y (2-y)]^{s-1}}{%
        f(y;s) = 2 * s * (1 - y) * (y * (2-y))^(s-1)}
  for \eqn{0<y<1} and shape parameter \eqn{0<s<1}.
  The mean of \eqn{Y} is
  \eqn{1 - 4^s [\Gamma(1+s)]^2 / \Gamma(2 + 2s)}{1 - 4^s
    [Gamma(1+s)]^2 / Gamma(2 + 2s)}
  (returned as the fitted values).


}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Topp, C. W. and F. C. Leone (1955).
  A family of J-shaped frequency functions.
  \emph{Journal of the American Statistical Association},
  \bold{50}, 209--219.



}
\author{ T. W. Yee }
\note{
  Fisher-scoring and Newton-Raphson are the same here.
  A related distribution is the triangle distribution.
  This \pkg{VGAM} family function handles multiple responses.



}

\seealso{
  \code{\link{Topple}},
  \code{\link{Triangle}}.



}
\examples{
tdata <- data.frame(y = rtopple(1000, logitlink(1, inverse = TRUE)))
tfit <- vglm(y ~ 1, topple, tdata, trace = TRUE, crit = "coef")
coef(tfit, matrix = TRUE)
Coef(tfit)
}
\keyword{models}
\keyword{regression}

\name{erlang}
\alias{erlang}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Erlang Distribution }
\description{
  Estimates the scale parameter of the Erlang distribution
  by maximum likelihood estimation.


}
\usage{
erlang(shape.arg, lscale = "loge", imethod = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shape.arg}{
  The shape parameters.
  The user must specify a positive integer, or integers for multiple responses.
  They are recycled \code{by.row = TRUE} according to \code{\link[base]{matrix}}.


  }
  \item{lscale}{
  Link function applied to the (positive) \eqn{scale} parameter.
  See \code{\link{Links}} for more choices.


  }
  \item{imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for more details.


  }

}
\details{
  The Erlang distribution is a special case of the gamma distribution
  with \emph{shape} that is a positive integer. If \code{shape.arg = 1}
  then it simplifies to the exponential distribution. As illustrated
  in the example below, the Erlang distribution is the distribution of
  the sum of \code{shape.arg} independent and identically distributed
  exponential random variates.


  The probability density function of the Erlang
  distribution is given by
  \deqn{f(y) = \exp(-y/scale) y^{shape-1} scale^{-shape} / \Gamma(shape)}{%
        f(y) =  exp(-y/scale) y^(shape-1) scale^(-shape) / gamma(shape)}
  for known positive integer \eqn{shape},
  unknown \eqn{scale > 0}  and \eqn{y > 0}.
  Here,
  \eqn{\Gamma(shape)}{gamma(shape)} is the gamma
  function, as in \code{\link[base:Special]{gamma}}.
  The mean of \emph{Y}
  is \eqn{\mu=shape \times scale}{mu=shape*scale} and
  its variance is \eqn{shape \times scale^2}{shape*scale^2}.
  The linear/additive predictor, by default, is
  \eqn{\eta=\log(scale)}{eta=log(scale)}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{
  Most standard texts on statistical distributions describe
  this distribution, e.g.,


Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}
\author{ T. W. Yee }

\note{
  Multiple responses are permitted.
  The \code{rate} parameter found in \code{\link{gammaR}}
  is \code{1/scale} here---see also \code{\link[stats]{rgamma}}.


}

\seealso{
  \code{\link{gammaR}},
  \code{\link{exponential}},
  \code{\link{simulate.vlm}}.


}
\examples{
rate <- exp(2); myshape <- 3
edata <- data.frame(y = rep(0, nn <- 1000))
for (ii in 1:myshape)
  edata <- transform(edata, y = y + rexp(nn, rate = rate))
fit <- vglm(y ~ 1, erlang(shape = myshape), data = edata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)  # Answer = 1/rate
1/rate
summary(fit)
}
\keyword{models}
\keyword{regression}


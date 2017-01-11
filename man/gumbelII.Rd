\name{gumbelII}
\alias{gumbelII}
%\alias{gumbelIIff}
%\alias{gumbelII.lsh}
%\alias{gumbelII3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gumbel-II Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter Gumbel-II distribution.

}
\usage{
gumbelII(lscale = "loge", lshape = "loge", iscale = NULL, ishape = NULL,
         probs.y = c(0.2, 0.5, 0.8), perc.out = NULL, imethod = 1,
         zero = "shape", nowarning = FALSE)
}
%- maybe also 'usage' for other objects documented here.
%        zero = "scale", nowarning = FALSE 20151128
\arguments{

  \item{nowarning}{ Logical. Suppress a warning? }


  \item{lshape, lscale}{
  Parameter link functions applied to the
  (positive) shape parameter (called \eqn{s} below) and
  (positive) scale parameter (called \eqn{b} below).
  See \code{\link{Links}} for more choices.


  }

% \item{eshape, escale}{
%         eshape = list(), escale = list(),
% Extra argument for the respective links.
% See \code{earg} in \code{\link{Links}} for general information.
% }

  Parameter link functions applied to the
  \item{ishape, iscale}{
  Optional initial values for the shape and scale parameters.


  }
  \item{imethod}{
  See \code{\link{weibullR}}.


  }
  \item{zero, probs.y}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{perc.out}{
  If the fitted values are to be quantiles then set this
  argument to be the percentiles of these, e.g., 50 for median.

  }

}
\details{
  The Gumbel-II density for a response \eqn{Y} is
  \deqn{f(y;b,s) = s y^{s-1} \exp[-(y/b)^s] / (b^s)}{%
        f(y;b,s) = s y^(s-1) * exp(-(y/b)^s) / [b^s]}
  for \eqn{b > 0}, \eqn{s > 0}, \eqn{y > 0}.
  The cumulative distribution function is
  \deqn{F(y;b,s) = \exp[-(y/b)^{-s}].}{%
        F(y;b,s) = exp(-(y/b)^(-s)).}
  The mean of \eqn{Y} is \eqn{b \, \Gamma(1 - 1/s)}{b * gamma(1 - 1/s)}
  (returned as the fitted values)
  when \eqn{s>1},
  and the variance is \eqn{b^2\,\Gamma(1-2/s)}{b^2 * Gamma(1-2/s)} when
  \eqn{s>2}.
  This distribution looks similar to \code{\link{weibullR}}, and is
  due to Gumbel (1954).


  This \pkg{VGAM} family function currently does not handle censored data.
  Fisher scoring is used to estimate the two parameters.
  Probably similar regularity conditions hold for this distribution
  compared to the Weibull distribution.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{


Gumbel, E. J. (1954).
Statistical theory of extreme values and some practical applications.
\emph{Applied Mathematics Series}, volume 33,
U.S. Department of Commerce, National Bureau of Standards, USA.



}
\author{ T. W. Yee }
\note{
  See \code{\link{weibullR}}.
  This \pkg{VGAM} family function handles multiple responses.



}
%\section{Warning}{
%  This function is under development to handle other censoring situations.
%  The version of this function which will handle censored data will be
%  called \code{cengumbelII()}. It is currently being written and will use
%  \code{\link{SurvS4}} as input.
%  It should be released in later versions of \pkg{VGAM}.
%
%
%  If the shape parameter is less than two then misleading inference may
%  result, e.g., in the \code{summary} and \code{vcov} of the object.
%
%
%}

\seealso{
    \code{\link{dgumbelII}},
    \code{\link{gumbel}},
    \code{\link{gev}}.


}
\examples{
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, heta1  = +1,
                          heta2  = -1 + 0.1 * x2,
                          ceta1 =  0,
                          ceta2 =  1)
gdata <- transform(gdata, shape1 = exp(heta1),
                          shape2 = exp(heta2),
                          scale1 = exp(ceta1),
                          scale2 = exp(ceta2))
gdata <- transform(gdata,
                   y1 = rgumbelII(nn, scale = scale1, shape = shape1),
                   y2 = rgumbelII(nn, scale = scale2, shape = shape2))

fit <- vglm(cbind(y1, y2) ~ x2,
            gumbelII(zero = c(1, 2, 3)), data = gdata, trace = TRUE)
coef(fit, matrix = TRUE)
vcov(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

\name{gompertz}
\alias{gompertz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gompertz Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter
  Gompertz distribution.

}
\usage{
gompertz(lscale = "loge", lshape = "loge",
         iscale = NULL,   ishape = NULL,
         nsimEIM = 500, zero = NULL, nowarning = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nowarning}{ Logical. Suppress a warning?
  Ignored for \pkg{VGAM} 0.9-7 and higher.


  }


  \item{lshape, lscale}{
  Parameter link functions applied to the
  shape parameter \code{a},
  scale parameter \code{scale}.
  All parameters are positive.
  See \code{\link{Links}} for more choices.


  }

% \item{eshape, escale}{
% List. Extra argument for each of the links.
%        eshape = list(), escale = list(),
% See \code{earg} in \code{\link{Links}} for general information.
% }

  \item{ishape, iscale}{
  Optional initial values.
  A \code{NULL} means a value is computed internally.


  }
  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}}.

  }
}
\details{
The Gompertz distribution has a cumulative distribution function
  \deqn{F(x;\alpha, \beta) = 1 - \exp[-(\alpha/\beta) \times (\exp(\beta x) - 1) ]}{%
        F(x;alpha, beta) = 1 - exp(-(alpha/beta) * (exp(beta * x) - 1) )}
which leads to a probability density function
  \deqn{f(x; \alpha, \beta) = \alpha \exp(\beta x)
                     \exp [-(\alpha/\beta) \times (\exp(\beta x) - 1) ]}{%
        f(x; alpha, beta) = alpha * exp[-beta * x] * exp[-(alpha/beta) * (exp(beta * x) - 1) ]}
  for \eqn{\alpha > 0}{a > 0},
      \eqn{\beta > 0}{b > 0},
      \eqn{x > 0}.
Here, \eqn{\beta} is called the scale parameter \code{scale},
and \eqn{\alpha} is called the shape parameter
(one could refer to \eqn{\alpha}{a} as a location parameter and \eqn{\beta}{b} as
a shape parameter---see Lenart (2012)).
The mean is involves an exponential integral function.
Simulated Fisher scoring is used and multiple responses are handled.


The Makeham distibution has an additional parameter compared to
the Gompertz distribution.
If \eqn{X} is defined to be the result of sampling from a Gumbel
distribution until a negative value \eqn{Z} is produced,
then \eqn{X = -Z} has a Gompertz distribution.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

Lenart, A. (2012)
The moments of the Gompertz distribution
and maximum likelihood estimation of its parameters.
\emph{Scandinavian Actuarial Journal}, in press.


}

\author{ T. W. Yee }
\section{Warning }{
The same warnings in \code{\link{makeham}} apply here too.


}

\seealso{
  \code{\link{dgompertz}},
  \code{\link{makeham}},
  \code{\link{simulate.vlm}}.


}

\examples{
\dontrun{
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, eta1  = -1,
                          eta2  = -1 + 0.2 * x2,
                          ceta1 =  1,
                          ceta2 = -1 + 0.2 * x2)
gdata <- transform(gdata, shape1 = exp(eta1),
                          shape2 = exp(eta2),
                          scale1 = exp(ceta1),
                          scale2 = exp(ceta2))
gdata <- transform(gdata, y1 = rgompertz(nn, scale = scale1, shape = shape1),
                          y2 = rgompertz(nn, scale = scale2, shape = shape2))

fit1 <- vglm(y1 ~ 1,  gompertz, data = gdata, trace = TRUE)
fit2 <- vglm(y2 ~ x2, gompertz, data = gdata, trace = TRUE)
coef(fit1, matrix = TRUE)
Coef(fit1)
summary(fit1)
coef(fit2, matrix = TRUE)
summary(fit2)
}
}
\keyword{models}
\keyword{regression}

% probs.y = c(0.20, 0.50, 0.80)





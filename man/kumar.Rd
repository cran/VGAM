\name{kumar}
\alias{kumar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kumaraswamy Distribution Family Function}
\description{
  Estimates the two parameters of the Kumaraswamy distribution
  by maximum likelihood estimation.

}
\usage{
kumar(lshape1 = "loge", lshape2 = "loge",
      ishape1 = NULL,   ishape2 = NULL, gshape1 = exp(2*ppoints(5) - 1),
      tol12 = 1.0e-4, zero = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape1, lshape2}{
  Link function for the two positive shape parameters,
  respectively, called \eqn{a} and \eqn{b} below.
  See \code{\link{Links}} for more choices.

  }

% \item{eshape1, eshape2}{
% List. Extra argument for each of the links.
% See \code{earg} in \code{\link{Links}} for general information.
%     eshape1 = list(), eshape2 = list(),
% }

  \item{ishape1, ishape2}{
  Numeric.
  Optional initial values for the two positive shape parameters.

  }
  \item{tol12}{
  Numeric and positive.
  Tolerance for testing whether the second shape parameter is either 1 or 2.
  If so then the working weights need to handle these singularities.

  }
  \item{gshape1}{
  Values for a grid search for the first shape parameter.
  See \code{\link{CommonVGAMffArguments}} for more information.


% Lower and upper limits for a grid search for the first shape parameter.



  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}}.

  }
}
\details{
  The Kumaraswamy distribution has density function
  \deqn{f(y;a = shape1,b = shape2)  =
  a b y^{a-1} (1-y^{a})^{b-1}}{%
  a*b*y^(a-1)*(1-y^a)^(b-1)}
  where \eqn{0 < y < 1} and the two shape parameters,
  \eqn{a} and \eqn{b}, are positive.
  The mean is \eqn{b \times Beta(1+1/a,b)}{b * Beta(1+1/a,b)}
  (returned as the fitted values) and the variance is
  \eqn{b \times Beta(1+2/a,b) -
      (b \times Beta(1+1/a,b))^2}{b * Beta(1+2/a,b) - (b * Beta(1+1/a,b))^2}.
  Applications of the Kumaraswamy distribution include the storage
  volume of a water reservoir.
  Fisher scoring is implemented.
  Handles multiple responses (matrix input).



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{
  Kumaraswamy, P. (1980).
  A generalized probability density function
  for double-bounded random processes.
  \emph{Journal of Hydrology},
  \bold{46}, 79--88.


  Jones, M. C. (2009).
  Kumaraswamy's distribution: A beta-type distribution with some
  tractability advantages.
  \emph{Statistical Methodology},
  \bold{6}, 70--81.


}
\author{ T. W. Yee }
%\note{
%
%}

\seealso{
  \code{\link{dkumar}},
  \code{\link{betaff}},
  \code{\link{simulate.vlm}}.


}
\examples{
shape1 <- exp(1); shape2 <- exp(2)
kdata <- data.frame(y = rkumar(n = 1000, shape1, shape2))
fit <- vglm(y ~ 1, kumar, data = kdata, trace = TRUE)
c(with(kdata, mean(y)), head(fitted(fit), 1))
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}


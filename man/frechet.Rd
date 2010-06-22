\name{frechet}
% \alias{frechet}
\alias{frechet2}
\alias{frechet3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Frechet Distribution Family Function }
\description{
   Maximum likelihood estimation of the 2-parameter and 3-parameter
   Frechet distribution.

}
\usage{
frechet2(location=0, lscale="loge", lshape="loglog", escale = list(),
         eshape = list(), iscale=NULL, ishape=3, zero=NULL)
frechet3(anchor=NULL, ldifference="loge", lscale="loge", lshape="loglog",
         edifference=list(), escale = list(), eshape = list(),
         ilocation=NULL, iscale=NULL, ishape=3, zero=NULL,
         effpos = .Machine$double.eps^0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{location}{
  Numeric. Location parameter.
  It is called \eqn{a} below.

  }
  \item{lscale, lshape}{
  Link functions for the parameters.
  See \code{\link{Links}} for more choices.

  }
  \item{iscale, ishape}{
  Initial value for the scale and shape parameters respectively. 
  Both parameters must be positive. 
  A \code{NULL} means it is chosen internally.

  }
  \item{edifference, escale, eshape}{
  Extra argument for the respective links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.  The value
  must be from the set \{1,2\}, corresponding respectively to \eqn{scale}
  and \eqn{shape}.  By default all linear/additive predictors are modelled
  as a linear combination of the explanatory variables.

  }
  \item{anchor}{
  An ``anchor'' point for estimating the location parameter. This must
  be a value no greater than \code{min(y)} where \code{y} is the response.
  The location parameter is \eqn{A - D} where
  \eqn{A} is the anchor,
  \eqn{D} is the ``difference'' (default is to make this positive).
  The default value of \code{anchor} means \code{min(y)} is chosen.

  }
  \item{ldifference}{
  Parameter link function for the difference \eqn{D} between the anchor
  point and the location parameter estimate.
  The default keeps this difference positive so that numerical
  problems are less likely to occur.

  }
  \item{ilocation}{
  Optional initial value for the location parameter.
  A good choice can speed up the convergence rate markedly.
  A \code{NULL} means it is chosen internally.

  }
  \item{effpos}{
  A small number, usually positive, to pass into the quasi-Newton
  updating function. A value greater than \code{effpos} is considered
  effectively positive.

  }
}
\details{
  The Frechet distribution has a density function that can be written
  \deqn{f(y) = \frac{sb}{(y-a)^2} [b/(y-a)]^{s-1} \, \exp[-(b/(y-a))^s] }{%
        f(y) = ((s*b) / (y-a)^2) * exp[-(b/(y-a))^s] * [b/(y-a)]^(s-1)}
  for \eqn{y>a} and scale parameter \eqn{b>0}.
  The positive shape parameter is \eqn{s}.
  The cumulative distribution function is 
  \deqn{F(y) = \exp[-(b/(y-a))^s]. }{%
        F(y) =  exp[-(b/(y-a))^s].}
  The mean of \eqn{Y} is \eqn{a + b \Gamma(1-1/s)}{a + b*gamma(1-1/s)}
  for \eqn{s>1} (these are returned as the fitted values).
  Note that the \code{\link{loglog}} link ensures \eqn{s>1}.
  The variance of \eqn{Y} is
  \eqn{b^2 [ \Gamma(1-2/s) - \Gamma^2(1-1/s)]}{b^2 * [gamma(1-2/s) - gamma(1-1/s)^2]}
  for \eqn{s>2}.


  \code{frechet2} has \eqn{a} known whereas \code{frechet3}
  estimates it.  Estimating \eqn{a} well requires a lot of data and
  a good choice of \code{ilocation} will help speed up convergence.
  For \code{frechet3} the default linear/additive predictors are
  \eqn{\log(D)}{log(D)}, 
  \eqn{\log(b)}{log(b)}, 
  \eqn{\log\log(s)}{log(log(s))}, respectively.
  It would be great if the first linear/additive predictor was a direct
  function of the location parameter, but this can run the risk that
  the estimate is out of range (i.e., greater than \code{min(y)}).

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\references{
Castillo, E., Hadi, A. S., Balakrishnan, N. Sarabia, J. S. (2005)
\emph{Extreme Value and Related Models with Applications in Engineering and Science},
Hoboken, N.J.: Wiley-Interscience.

}
\author{ T. W. Yee }
\section{Warning}{
  Convergence for \code{frechet3} can be very slow, especially if the
  initial value for the location parameter is poor. Setting something
  like \code{maxit=200, trace=TRUE} is a good idea.

}

\note{ 
  These family functions use the BFGS quasi-Newton update formula for the
  working weight matrices.  Consequently the estimated variance-covariance
  matrix may be inaccurate or simply wrong! The standard errors must be
  therefore treated with caution; these are computed in functions such
  as \code{vcov()} and \code{summary()}.


  If \code{fit} is a \code{frechet3} fit then \code{fit@extra$location}
  is the final estimate of the location parameter, and
  \code{fit@extra$LHSanchor} is the anchor point.


}
\seealso{
  \code{\link{rfrechet}}.

}
\examples{
y = rfrechet(n <- 1000, shape=exp(exp(0)))
\dontrun{hist(y)}
fit2 = vglm(y ~ 1, frechet2, trace=TRUE, maxit=155)
fit3 = vglm(y ~ 1, frechet3(ilocation=0), trace=TRUE, maxit=155)
coef(fit3, matrix=TRUE)
Coef(fit3)
head(fitted(fit3))
mean(y)
head(weights(fit3, type="w"))
vcov(fit3)   # caution needed!
head(fit3@extra$location)  # Estimate of the location parameter
fit3@extra$LHSanchor       # Anchor point
min(y)
}
\keyword{models}
\keyword{regression}


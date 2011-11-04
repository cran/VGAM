\name{poissonff}
%\alias{poisson}
\alias{poissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Poisson Family Function }
\description{
  Family function for a generalized linear model fitted to
  Poisson responses.
  The dispersion parameters may be known or unknown.

}
\usage{
poissonff(link = "loge", earg=list(), dispersion = 1, onedpar = FALSE,
          imu = NULL, imethod = 1, parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, earg}{
  Link function and extra argument applied to the mean or means.
  See \code{\link{Links}} for more choices
  and information.

  }
  \item{dispersion}{
  Dispersion parameter. By default, maximum
  likelihood is used to estimate the model because it is known.
  However, the user can specify
  \code{dispersion = 0} to have it estimated, or
  else specify a known positive value (or values if the response
  is a matrix---one value per column).

  }
  \item{onedpar}{
  One dispersion parameter? If the response is a matrix,
  then a separate
  dispersion parameter will be computed for each response (column), by default.
  Setting \code{onedpar=TRUE} will pool them so that there is only
  one dispersion parameter to be estimated.

  }
  \item{parallel}{
  A logical or formula. Used only if the response is a matrix.

  }
  \item{imu, imethod}{
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{zero}{
  An integer-valued vector specifying which linear/additive predictors
  are modelled as intercepts only.  The values must be from the set
  \{1,2,\ldots,\eqn{M}\}, where \eqn{M} is the number of columns of the
  matrix response.

  }
}
\details{
  \eqn{M} defined above is the number of linear/additive predictors.


  If the dispersion parameter is unknown, then the resulting estimate
  is not fully a maximum likelihood estimate.


  A dispersion parameter that is less/greater than unity corresponds to
  under-/over-dispersion relative to the Poisson model.  Over-dispersion
  is more common in practice.


  When fitting a Quadratic RR-VGLM (see \code{\link{cqo}}), the
  response is a matrix of \eqn{M}, say, columns (e.g., one column
  per species). Then there will be \eqn{M} dispersion parameters
  (one per column of the response matrix) if \code{dispersion = 0} and
  \code{onedpar = FALSE}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as
  \code{\link{vglm}},
  \code{\link{vgam}},
  \code{\link{rrvglm}},
  \code{\link{cqo}},
  and \code{\link{cao}}.


}
\references{

 McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}

\author{ Thomas W. Yee }

\note{

  This function will handle a matrix response automatically.


  The call \code{poissonff(dispersion=0, ...)} is equivalent to
  \code{quasipoissonff(...)}.  The latter was written so that R users
  of \code{quasipoisson()} would only need to add a  ``\code{ff}''
  to the end of the family function name.


  Regardless of whether the dispersion parameter is to be estimated or
  not, its value can be seen from the output from the \code{summary()}
  of the object.


% With the introduction of name spaces for the \pkg{VGAM} package,
% \code{"ff"} can be dropped for this family function.


}
\section{Warning }{
    With a multivariate response, assigning a known dispersion parameter
    for \emph{each} response is not handled well yet.  Currently, only
    a single known dispersion parameter is handled well.


}

\seealso{
    \code{\link{Links}},
    \code{\link{quasipoissonff}},
    \code{\link{genpoisson}},
    \code{\link{zipoisson}},
    \code{\link{skellam}},
    \code{\link{mix2poisson}},
    \code{\link{cenpoisson}},
    \code{\link{ordpoisson}},
    \code{\link{amlpoisson}},
    \code{\link{invbinomial}},
    \code{\link{loge}},
    \code{\link{polf}},
    \code{\link{rrvglm}},
    \code{\link{cqo}},
    \code{\link{cao}},
    \code{\link{binomialff}},
    \code{\link{quasibinomialff}},
    \code{\link[stats]{poisson}},
    \code{\link{poissonp}}.


}
\examples{
poissonff()

pdat = data.frame(x = rnorm(nn <- 100))
pdat = transform(pdat, y = rpois(nn, exp(1+x)))
(fit = vglm(y ~ x, family = poissonff, pdat))
coef(fit, matrix = TRUE)

nn = 200
cdat = data.frame(x2 = rnorm(nn), x3 = rnorm(nn), x4 = rnorm(nn))
cdat = transform(cdat, lv1 = 0 + x3 - 2*x4)
cdat = transform(cdat, lambda1 = exp(3 - 0.5 *  (lv1-0)^2),
                       lambda2 = exp(2 - 0.5 *  (lv1-1)^2),
                       lambda3 = exp(2 - 0.5 * ((lv1+4)/2)^2))
cdat = transform(cdat, y1 = rpois(nn, lambda1),
                       y2 = rpois(nn, lambda2),
                       y3 = rpois(nn, lambda3))
# vvv p1 = cqo(cbind(y1,y2,y3) ~ x2 + x3 + x4, poissonff, cdat,
# vvv          EqualTol = FALSE, ITol = FALSE)
# vvv summary(p1)  # # Three dispersion parameters are all unity
\dontrun{ lvplot(p1, y = TRUE, lcol = 2:4, pch = 2:4, pcol = 2:4, rug = FALSE) }
}
\keyword{models}
\keyword{regression}


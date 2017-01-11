\name{bifgmcop}
\alias{bifgmcop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Farlie-Gumbel-Morgenstern's Bivariate Distribution Family Function }
\description{
  Estimate the association parameter of
  Farlie-Gumbel-Morgenstern's bivariate
  distribution by maximum likelihood estimation.

}
\usage{
bifgmcop(lapar = "rhobit", iapar = NULL, imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lapar, iapar, imethod}{
  Details at \code{\link{CommonVGAMffArguments}}.
  See \code{\link{Links}} for more link function choices.


  }
}
\details{
  The cumulative distribution function is
  \deqn{P(Y_1 \leq y_1, Y_2 \leq y_2) = y_1 y_2
             ( 1 + \alpha (1 - y_1) (1 - y_2) ) }{%
        P(Y1 <= y1, Y2 <= y2) =
        y1 * y2 * ( 1 + alpha * (1 - y1) * (1 - y2) ) }
  for \eqn{-1 < \alpha < 1}{-1 < alpha < 1}.
  The support of the function is the unit square.
  The marginal distributions are the standard uniform distributions.
  When \eqn{\alpha = 0}{alpha=0} the random variables are
  independent.


% A variant of Newton-Raphson is used, which only seems to work for an
% intercept model.
% It is a very good idea to set \code{trace=TRUE}.
% This \pkg{VGAM} family function is prone to numerical difficulties.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}

\references{

Castillo, E., Hadi, A. S., Balakrishnan, N. Sarabia, J. S. (2005)
\emph{Extreme Value and Related Models with Applications in Engineering and Science},
Hoboken, NJ, USA: Wiley-Interscience.


Smith, M. D. (2007)
Invariance theorems for Fisher information.
\emph{Communications in Statistics---Theory and Methods},
\bold{36}(12), 2213--2222.



}
\author{ T. W. Yee }
\note{
  The response must be a two-column matrix.  Currently, the fitted
  value is a matrix with two columns and values equal to 0.5.
  This is because each marginal distribution corresponds to a standard
  uniform distribution.


% This \pkg{VGAM} family function should be used with caution.


}

\seealso{
  \code{\link{rbifgmcop}},
  \code{\link{bifrankcop}},
  \code{\link{bifgmexp}},
  \code{\link{simulate.vlm}}.


}
\examples{
ymat <- rbifgmcop(n = 1000, apar = rhobit(3, inverse = TRUE))
\dontrun{plot(ymat, col = "blue")}
fit <- vglm(ymat ~ 1, fam = bifgmcop, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
head(fitted(fit))
}
\keyword{models}
\keyword{regression}

% for real \eqn{\alpha}{alpha} (the range is data-dependent).

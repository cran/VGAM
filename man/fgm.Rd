\name{fgm}
\alias{fgm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Farlie-Gumbel-Morgenstern's Bivariate Distribution Family Function }
\description{
  Estimate the association parameter of 
  Farlie-Gumbel-Morgenstern's bivariate
  distribution by maximum likelihood estimation.

}
\usage{
fgm(lapar="rhobit", earg=list(), iapar=NULL, method.init=1, nsimEIM=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lapar}{
  Link function applied to the association parameter
  \eqn{\alpha}{alpha}, which is real.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{iapar}{
  Numeric. Optional initial value for \eqn{\alpha}{alpha}.
  By default, an initial value is chosen internally.
  If a convergence failure occurs try assigning a different value.
  Assigning a value will override the argument \code{method.init}.

  }
  \item{method.init}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method. If failure to converge occurs
  try the other value, or else specify a value for \code{ia}.

  }
  \item{nsimEIM}{
  See \code{\link{CommonVGAMffArguments}} for more information.

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
Hoboken, N.J.: Wiley-Interscience.

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
  \code{\link{rfgm}},
  \code{\link{frank}},
  \code{\link{morgenstern}}.

}
\examples{
ymat = rfgm(n = 1000, alpha=rhobit(3, inverse=TRUE))
\dontrun{plot(ymat, col="blue")}
fit = vglm(ymat ~ 1, fam=fgm, trace=TRUE)
coef(fit, matrix=TRUE)
Coef(fit)
head(fitted(fit))
}
\keyword{models}
\keyword{regression}

% for real \eqn{\alpha}{alpha} (the range is data-dependent).

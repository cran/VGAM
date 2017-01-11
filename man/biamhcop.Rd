\name{biamhcop}
\alias{biamhcop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Ali-Mikhail-Haq Distribution Family Function }
\description{
  Estimate the association parameter of
  Ali-Mikhail-Haq's bivariate
  distribution by maximum likelihood estimation.

}
\usage{
biamhcop(lapar = "rhobit", iapar = NULL, imethod = 1, nsimEIM = 250)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lapar}{
  Link function applied to the association parameter
  \eqn{\alpha}{alpha}, which is real
  and \eqn{-1 < \alpha < 1}{-1 < alpha < 1}.
  See \code{\link{Links}} for more choices.


  }
  \item{iapar}{
  Numeric. Optional initial value for \eqn{\alpha}{alpha}.
  By default, an initial value is chosen internally.
  If a convergence failure occurs try assigning a different value.
  Assigning a value will override the argument \code{imethod}.


  }
  \item{imethod}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method. If failure to converge occurs
  try the other value, or else specify a value for \code{iapar}.


  }
  \item{nsimEIM}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  The cumulative distribution function is
  \deqn{P(Y_1 \leq y_1, Y_2 \leq y_2) = y_1 y_2
          / ( 1 - \alpha (1 - y_1) (1 - y_2) ) }{%
        P(Y1 < =  y1, Y2 < =  y2)  =
        y1 * y2 / ( 1 - alpha * (1 - y1) * (1 - y2) ) }
  for \eqn{-1 < \alpha < 1}{-1 < alpha < 1}.
  The support of the function is the unit square.
  The marginal distributions are the standard uniform distributions.
  When \eqn{\alpha = 0}{alpha = 0} the random variables are
  independent.
  This is an Archimedean copula.



% A variant of Newton-Raphson is used, which only seems to work for an
% intercept model.
% It is a very good idea to set \code{trace = TRUE}.
% This \pkg{VGAM} family function is prone to numerical difficulties.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}

\references{

%Hutchinson, T. P. and Lai, C. D. (1990)
%\emph{Continuous Bivariate Distributions, Emphasising Applications},
%Adelaide, South Australia: Rumsby Scientific Publishing.

Balakrishnan, N. and Lai, C.-D. (2009)
\emph{Continuous Bivariate Distributions},
2nd ed.
New York: Springer.


}
\author{ T. W. Yee and C. S. Chee }
\note{
  The response must be a two-column matrix.  Currently, the fitted
  value is a matrix with two columns and values equal to 0.5.
  This is because each marginal distribution corresponds to a standard
  uniform distribution.


}

\seealso{
  \code{\link{rbiamhcop}},
  \code{\link{bifgmcop}},
  \code{\link{bigumbelIexp}},
  \code{\link{rbilogis}},
  \code{\link{simulate.vlm}}.


}
\examples{
ymat <- rbiamhcop(1000, apar = rhobit(2, inverse = TRUE))
fit <- vglm(ymat ~ 1, biamhcop, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{models}
\keyword{regression}





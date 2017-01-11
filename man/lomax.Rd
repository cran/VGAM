\name{lomax}
\alias{lomax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lomax Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter
  Lomax distribution.

}
\usage{
lomax(lscale = "loge", lshape3.q = "loge", iscale = NULL,
      ishape3.q = NULL, imethod = 1, gscale = exp(-5:5),
      gshape3.q = seq(0.75, 4, by = 0.25),
      probs.y = c(0.25, 0.5, 0.75), zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale, lshape3.q}{
  Parameter link function applied to the
  (positive) parameters \code{scale} and \code{q}.
  See \code{\link{Links}} for more choices.


  }
  \item{iscale, ishape3.q, imethod}{
  See \code{\link{CommonVGAMffArguments}} for information.
  For \code{imethod = 2} a good initial value for
  \code{iscale} is needed to obtain a good estimate for
  the other parameter.


  }

  \item{gscale, gshape3.q, zero, probs.y}{
  See
  \code{\link{CommonVGAMffArguments}}.


  }

% \item{zero}{
% An integer-valued vector specifying which
% linear/additive predictors are modelled as intercepts only.
% Here, the values must be from the set \{1,2\} which correspond to
% \code{scale}, \code{q}, respectively.
% }




}
\details{
  The 2-parameter Lomax distribution is the 4-parameter
  generalized beta II distribution with shape parameters \eqn{a=p=1}.
It is probably more widely known as the Pareto (II) distribution.
  It is also the 3-parameter Singh-Maddala distribution
  with shape parameter \eqn{a=1}, as well as the
  beta distribution of the second kind with \eqn{p=1}.
  More details can be found in Kleiber and Kotz (2003).


The Lomax distribution has density
  \deqn{f(y) = q / [b \{1 + y/b\}^{1+q}]}{%
        f(y) = q / [b (1 + y/b)^(1+q)]}
  for \eqn{b > 0}, \eqn{q > 0}, \eqn{y \geq 0}{y >= 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and \code{q} is a shape parameter.
The cumulative distribution function is
  \deqn{F(y) = 1 - [1 + (y/b)]^{-q}.}{%
        F(y) = 1 - [1 + (y/b)]^(-q).}
The mean is
  \deqn{E(Y) = b/(q-1)}{%
        E(Y) = b/(q-1)}
provided \eqn{q > 1}; these are returned as the fitted values.
This family function handles multiple responses.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.


}

\author{ T. W. Yee }
\note{
  See the notes in \code{\link{genbetaII}}.


}

\seealso{
    \code{\link{Lomax}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{inv.lomax}},
    \code{\link{paralogistic}},
    \code{\link{inv.paralogistic}},
    \code{\link{simulate.vlm}}.

}

\examples{
ldata <- data.frame(y = rlomax(n = 1000, scale =  exp(1), exp(2)))
fit <- vglm(y ~ 1, lomax, data = ldata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

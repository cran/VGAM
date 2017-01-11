\name{quasipoissonff}
%\alias{quasipoisson}
\alias{quasipoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quasi-Poisson Family Function }
\description{
  Fits a generalized linear model to a Poisson response, where
  the dispersion parameter is unknown.
}
\usage{
quasipoissonff(link = "loge", onedpar = FALSE,
               parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function. See \code{\link{Links}} for more choices.

  }
  \item{onedpar}{
  One dispersion parameter? If the response is a matrix,
  then a separate
  dispersion parameter will be computed for each response (column),
  by default.
  Setting \code{onedpar=TRUE} will pool them so that there is only
  one dispersion parameter to be estimated.

  }
  \item{parallel}{
  A logical or formula. Used only if the response is a matrix.

  }
  \item{zero}{
  Can be an integer-valued vector specifying which linear/additive predictors
  are modelled as intercepts only.  The values must be from the set
  \{1,2,\ldots,\eqn{M}\}, where \eqn{M} is the number of columns of the
  matrix response.
  See \code{\link{CommonVGAMffArguments}} for more information.



  }
}
\details{
  \eqn{M} defined above is the number of linear/additive predictors.


  If the dispersion parameter is unknown, then the resulting estimate
  is not fully a maximum likelihood estimate.


  A dispersion parameter that is less/greater than unity corresponds to
  under-/over-dispersion relative to the Poisson model.  Over-dispersion
  is more common in practice.


  When fitting a Quadratic RR-VGLM, the response is a matrix of \eqn{M},
  say, columns (e.g., one column per species). Then there will be \eqn{M}
  dispersion parameters (one per column of the response matrix).


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


  The call \code{poissonff(dispersion = 0, ...)} is equivalent to
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
  See the warning in \code{\link{quasibinomialff}}.


}

\seealso{
    \code{\link{poissonff}},
    \code{\link{negbinomial}},
    \code{\link{loge}},
    \code{\link{rrvglm}},
    \code{\link{cqo}},
    \code{\link{cao}},
    \code{\link{binomialff}},
    \code{\link{quasibinomialff}},
    \code{\link[stats]{quasipoisson}}.


}
\examples{
quasipoissonff()

\dontrun{n <- 200; p <- 5; S <- 5
mydata <- rcqo(n, p, S, fam = "poisson", eq.tol = FALSE)
myform <- attr(mydata, "formula")
p1 <- cqo(myform, fam = quasipoissonff, eq.tol = FALSE, data = mydata)
sort(deviance(p1, history = TRUE))  # A history of all the iterations
lvplot(p1, y = TRUE, lcol = 1:S, pch = 1:S, pcol = 1:S)
summary(p1)  # The dispersion parameters are estimated
}}
\keyword{models}
\keyword{regression}

%n = 100
%x2 = rnorm(n)
%x3 = rnorm(n)
%x4 = rnorm(n)
%lv1 = 0 + x3 - 2*x4
%lambda1 = exp(3 - 0.5 * (lv1-0)^2)
%lambda2 = exp(2 - 0.5 * (lv1-1)^2)
%lambda3 = exp(2 - 0.5 * ((lv1+4)/2)^2)
%y1 = rpois(n, lambda1)
%y2 = rpois(n, lambda2)
%y3 = rpois(n, lambda3)
%p1 = cqo(cbind(y1,y2,y3) ~ x2 + x3 + x4, quasipoissonff)
%lvplot(p1, y=TRUE, lcol=1:3, pch=1:3, pcol=1:3)
%summary(p1)  # Three dispersion parameters are estimated

\name{quasibinomialff}
%\alias{quasibinomial}
\alias{quasibinomialff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quasi-Binomial Family Function }
\description{
  Family function for fitting generalized linear models to binomial
  responses, where the dispersion parameters are unknown.

}
\usage{
quasibinomialff(link = "logit", multiple.responses = FALSE,
                onedpar = !multiple.responses, parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{ Link function. See \code{\link{Links}} for more choices.

  }
  \item{multiple.responses}{
   Multiple responses? If \code{TRUE}, then the response is interpreted
   as \eqn{M} binary responses, where \eqn{M} is the number of columns
   of the response matrix. In this case, the response matrix should have
   zero/one values only.


  If \code{FALSE} and the response is a (2-column) matrix, then the
  number of successes is given in the first column and the second column
  is the number of failures.


  }
  \item{onedpar}{
  One dispersion parameter? If \code{multiple.responses}, then
  a separate dispersion
  parameter will be computed for each response (column), by default.
  Setting \code{onedpar=TRUE} will pool them so that there is only one
  dispersion parameter to be estimated.


  }
  \item{parallel}{
  A logical or formula. Used only if \code{multiple.responses}
  is \code{TRUE}.  This
  argument allows for the parallelism assumption whereby the regression
  coefficients for a variable is constrained to be equal over the \eqn{M}
  linear/additive predictors.


  }
  \item{zero}{
  Can be an integer-valued vector specifying which linear/additive predictors
  are modelled as intercepts only.  The values must be from the set
  \{1,2,\ldots,\eqn{M}\}, where \eqn{M} is the number of columns of
  the matrix response.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  The final model is not fully estimated by maximum likelihood since the
  dispersion parameter is unknown (see pp.124--8 of McCullagh and Nelder
  (1989) for more details).


  A dispersion parameter that is less/greater than unity corresponds to
  under-/over-dispersion relative to the binomial model.  Over-dispersion
  is more common in practice.


  Setting \code{multiple.responses=TRUE} is necessary
  when fitting a Quadratic RR-VGLM (see \code{\link{cqo}})
  because the response will be a matrix of \eqn{M} columns
  (e.g., one column per species). Then there will be \eqn{M}
  dispersion parameters (one per column of the response).


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
  If \code{multiple.responses} is \code{FALSE} (the default),
  then the response can be of one of three formats: a factor
  (first level taken as success), a vector of proportions of
  success, or a 2-column matrix (first column = successes)
  of counts.  The argument \code{weights} in the modelling
  function can also be specified. In particular, for a
  general vector of proportions, you will need to specify
  \code{weights} because the number of trials is needed.


  If \code{multiple.responses} is \code{TRUE}, then the
  matrix response can only be of one format: a matrix of
  1s and 0s (1=success).


  This function is only a front-end to the \pkg{VGAM} family function
  \code{binomialff()}; indeed, \code{quasibinomialff(...)} is equivalent
  to \code{binomialff(..., dispersion=0)}.  Here,  the argument
  \code{dispersion=0} signifies that the dispersion parameter is to
  be estimated.


  Regardless of whether the dispersion parameter is to be estimated or
  not, its value can be seen from the output from the \code{summary()}
  of the object.


% With the introduction of name spaces for the \pkg{VGAM} package,
% \code{"ff"} can be dropped for this family function.


}

\section{Warning }{
The log-likelihood pertaining to the ordinary family
is used to test for convergence during estimation,
and is printed out in the summary.


}


\seealso{
    \code{\link{binomialff}},
    \code{\link{rrvglm}},
    \code{\link{cqo}},
    \code{\link{cao}},
    \code{\link{logit}},
    \code{\link{probit}},
    \code{\link{cloglog}},
    \code{\link{cauchit}},
    \code{\link{poissonff}},
    \code{\link{quasipoissonff}},
    \code{\link[stats]{quasibinomial}}.


}
\examples{
quasibinomialff()
quasibinomialff(link = "probit")

# Nonparametric logistic regression
hunua <- transform(hunua, a.5 = sqrt(altitude))  # Transformation of altitude
fit1 <- vglm(agaaus ~ poly(a.5, 2), quasibinomialff, hunua)
fit2 <- vgam(agaaus ~ s(a.5, df = 2), quasibinomialff, hunua)
\dontrun{
plot(fit2, se = TRUE, llwd = 2, lcol = "orange", scol = "orange",
     xlab = "sqrt(altitude)", ylim = c(-3, 1),
     main = "GAM and quadratic GLM fitted to species data")
plotvgam(fit1, se = TRUE, lcol = "blue", scol = "blue", add = TRUE, llwd = 2)
}
fit1@misc$dispersion # dispersion parameter
logLik(fit1)

# Here, the dispersion parameter defaults to 1
fit0 <- vglm(agaaus ~ poly(a.5, 2), binomialff, hunua)
fit0@misc$dispersion # dispersion parameter
}
\keyword{models}
\keyword{regression}


%AIC(fit1)
%AIC(fit0)

\name{binomialff}
%\alias{binomial}
\alias{binomialff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Binomial Family Function }
\description{
  Family function for fitting generalized linear models to binomial
  responses, where the dispersion parameter may be known or unknown.

}
\usage{
binomialff(link = "logit", dispersion = 1, multiple.responses = FALSE,
           onedpar = !multiple.responses, parallel = FALSE,
           zero = NULL, bred = FALSE, earg.link = FALSE)

}
%- maybe also 'usage' for other objects documented here.
% apply.parint = FALSE,
\arguments{

  \item{link}{
  Link function;
  see \code{\link{Links}} and
  \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{dispersion}{
  Dispersion parameter. By default, maximum likelihood is used to
  estimate the model because it is known.  However, the user can specify
  \code{dispersion = 0} to have it estimated, or else specify a known
  positive value (or values if \code{multiple.responses} is \code{TRUE}).


  }
  \item{multiple.responses}{
  Multivariate response? If \code{TRUE}, then the response is interpreted
  as \eqn{M} independent binary responses, where \eqn{M} is the number
  of columns of the response matrix. In this case, the response matrix
  should have \eqn{Q} columns consisting of counts (successes),
  and the \code{weights} argument should have \eqn{Q} columns
  consisting of the number of trials (successes plus failures).


% zero/one values only.



  If \code{FALSE} and the response is a (2-column) matrix, then the number
  of successes is given in the first column, and the second column is
  the number of failures.


  }
  \item{onedpar}{
  One dispersion parameter? If \code{multiple.responses}, then a separate dispersion
  parameter will be computed for each response (column), by default.
  Setting \code{onedpar = TRUE} will pool them so that there is only one
  dispersion parameter to be estimated.


  }
  \item{parallel}{
  A logical or formula. Used only if \code{multiple.responses} is \code{TRUE}.  This
  argument allows for the parallelism assumption whereby the regression
  coefficients for a variable is constrained to be equal over the \eqn{M}
  linear/additive predictors.
  If \code{parallel = TRUE} then the constraint is not applied to the
  intercepts.


  }
  \item{zero}{
  An integer-valued vector specifying which linear/additive predictors
  are modelled as intercepts only.  The values must be from the set
  \{1,2,\ldots,\eqn{M}\}, where \eqn{M} is the number of columns of the
  matrix response.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{earg.link}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{bred}{
  Details at \code{\link{CommonVGAMffArguments}}.
  Setting \code{bred = TRUE} should work for
  multiple responses (\code{multiple.responses = TRUE}) and
  all \pkg{VGAM} link functions;
  it has been tested for
  \code{\link{logit}} only (and it gives similar
  results to \pkg{brglm} but not identical),
  and further testing is required.
  One result from fitting bias reduced binary regression
  is that finite regression coefficients occur when
  the data is separable (see example below).


  }
}
\details{
  This function is largely to mimic \code{\link[stats:Binomial]{binomial}},
  however there are some differences.


  If the dispersion parameter is unknown, then the resulting estimate
  is not fully a maximum likelihood estimate (see pp.124--8 of McCullagh
  and Nelder, 1989).


  A dispersion parameter that is less/greater than unity corresponds to
  under-/over-dispersion relative to the binomial model.  Over-dispersion
  is more common in practice.


  Setting \code{multiple.responses = TRUE} is necessary when fitting a Quadratic RR-VGLM
  (see \code{\link{cqo}}) because the response is a matrix of \eqn{M}
  columns (e.g., one column per species). Then there will be \eqn{M}
  dispersion parameters (one per column of the response matrix).


  When used with \code{\link{cqo}} and \code{\link{cao}}, it may be
  preferable to use the \code{\link{cloglog}} link.


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


  Altman, M. and Gill, J. and McDonald, M. P. (2004)
  \emph{Numerical Issues in Statistical Computing for the Social
        Scientist}, Hoboken, NJ, USA: Wiley-Interscience.


  Ridout, M. S. (1990)
  Non-convergence of Fisher's method of scoring---a simple example.
  \emph{GLIM Newsletter}, 20(6).


}
\author{ Thomas W. Yee }

\note{
  If \code{multiple.responses} is \code{FALSE} (default) then
  the response can be of one
  of two formats:
  a factor (first level taken as failure),
  or a 2-column matrix (first column = successes) of counts.
  The argument \code{weights} in the modelling function can
  also be specified as any vector of positive values.
  In general, 1 means success and 0 means failure
  (to check, see the \code{y} slot of the fitted object).
  Note that a general vector of proportions of success is no
  longer accepted.


  The notation \eqn{M} is used to denote the number of linear/additive
  predictors.


  If \code{multiple.responses} is \code{TRUE}, then the matrix response
  can only be of one format: a matrix of 1's and 0's (1 = success).


  The call \code{binomialff(dispersion = 0, ...)} is equivalent to
  \code{quasibinomialff(...)}.  The latter was written so that R users
  of \code{quasibinomial()} would only need to add a  ``\code{ff}''
  to the end of the family function name.


  Regardless of whether the dispersion parameter is to be estimated or
  not, its value can be seen from the output from the \code{summary()}
  of the object.


% With the introduction of name spaces for the \pkg{VGAM} package,
% \code{"ff"} can be dropped for this family function.


  Fisher scoring is used. This can sometimes fail to converge by
  oscillating between successive iterations (Ridout, 1990).
  See the example below.


}
\seealso{
    \code{\link{quasibinomialff}},
    \code{\link{Links}},
    \code{\link{rrvglm}},
    \code{\link{cqo}},
    \code{\link{cao}},
    \code{\link{betabinomial}},
    \code{\link{posbinomial}},
    \code{\link{zibinomial}},
    \code{\link{double.expbinomial}},
    \code{\link{seq2binomial}},
    \code{\link{amlbinomial}},
    \code{\link{simplex}},
    \code{\link[stats:Binomial]{binomial}},
    \code{\link{simulate.vlm}},
    \pkg{safeBinaryRegression}.



%   \code{\link{matched.binomial}},



}
\section{Warning }{
    With a multivariate response, assigning a known dispersion parameter
    for \emph{each} response is not handled well yet.
    Currently, only a single known dispersion parameter is handled well.


    See the above note regarding \code{bred}.



  The maximum likelihood estimate will not exist if the data is
  \emph{completely separable} or \emph{quasi-completely separable}.
  See Chapter 10 of Altman et al. (2004) for more details,
  and \pkg{safeBinaryRegression}.
  Yet to do: add a \code{sepcheck = TRUE}, say, argument to detect this
  problem and give an appropriate warning.


}

\examples{
quasibinomialff()
quasibinomialff(link = "probit")

shunua <- hunua[sort.list(with(hunua, altitude)), ]  # Sort by altitude
fit <- vglm(agaaus ~ poly(altitude, 2), binomialff(link = cloglog),
            data = shunua)
\dontrun{
plot(agaaus ~ jitter(altitude), shunua, col = "blue", ylab = "P(Agaaus = 1)",
     main = "Presence/absence of Agathis australis", las = 1)
with(shunua, lines(altitude, fitted(fit), col = "orange", lwd = 2)) }


# Fit two species simultaneously
fit2 <- vgam(cbind(agaaus, kniexc) ~ s(altitude),
             binomialff(multiple.responses = TRUE), data = shunua)
\dontrun{
with(shunua, matplot(altitude, fitted(fit2), type = "l",
     main = "Two species response curves", las = 1)) }


# Shows that Fisher scoring can sometime fail. See Ridout (1990).
ridout <- data.frame(v = c(1000, 100, 10), r = c(4, 3, 3), n = c(5, 5, 5))
(ridout <- transform(ridout, logv = log(v)))
# The iterations oscillates between two local solutions:
glm.fail <- glm(r / n ~ offset(logv) + 1, weight = n,
               binomial(link = 'cloglog'), ridout, trace = TRUE)
coef(glm.fail)
# vglm()'s half-stepping ensures the MLE of -5.4007 is obtained:
vglm.ok <- vglm(cbind(r, n-r) ~ offset(logv) + 1,
               binomialff(link = cloglog), ridout, trace = TRUE)
coef(vglm.ok)


# Separable data
set.seed(123)
threshold <- 0
bdata <- data.frame(x2 = sort(rnorm(nn <- 100)))
bdata <- transform(bdata, y1 = ifelse(x2 < threshold, 0, 1))
fit <- vglm(y1 ~ x2, binomialff(bred = TRUE),
            data = bdata, criter = "coef", trace = TRUE)
coef(fit, matrix = TRUE)  # Finite!!
summary(fit)
\dontrun{ plot(depvar(fit) ~ x2, data = bdata, col = "blue", las = 1)
lines(fitted(fit) ~ x2, data = bdata, col = "orange")
abline(v = threshold, col = "gray", lty = "dashed") }
}
\keyword{models}
\keyword{regression}


% a vector of proportions of success,
% In particular, for a general vector of proportions,
% you will need to specify \code{weights} because the number of trials
% is needed.
% To input general positive values into the \code{weights} argument of
% \code{\link{vglm}}/\code{\link{vgam}} one needs to input a 2-column
% response.




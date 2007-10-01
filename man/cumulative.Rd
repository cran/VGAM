\name{cumulative}
\alias{cumulative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Ordinal Regression with Cumulative Probabilities }
\description{
  Fits a cumulative logit/probit/cloglog/cauchit/...
  regression model to an ordered (preferably) factor response. 

}
\usage{
cumulative(link = "logit", earg = list(),
           parallel = FALSE, reverse = FALSE,
           mv = FALSE, intercept.apply = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  In the following, the response \eqn{Y} is assumed to be a factor
  with ordered values \eqn{1,2,\dots,M+1}, so that
  \eqn{M} is the number of linear/additive predictors
  \eqn{\eta_j}{eta_j}.

  \item{link}{
  Link function applied to the \eqn{M} cumulative probabilities. 
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link function.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{parallel}{
  A logical, or formula specifying which terms have
  equal/unequal coefficients.

  }
  \item{reverse}{
  Logical.
  By default, the cumulative probabilities used are
  \eqn{P(Y\leq 1)}{P(Y<=1)}, \eqn{P(Y\leq 2)}{P(Y<=2)},
  \dots, \eqn{P(Y\leq M)}{P(Y<=M)}.
  If \code{reverse} is \code{TRUE}, then 
  \eqn{P(Y\geq 2)}{P(Y>=2)}, \eqn{P(Y\geq 3)}{P(Y>=3)}, \dots,
  \eqn{P(Y\geq M+1)}{P(Y>=M+1)} will be used.

  This should be set to \code{TRUE} for \code{link=}
  \code{\link{golf}},
  \code{\link{polf}},
  \code{\link{nbolf}}.
  For these links the cutpoints must be an increasing sequence;
  if \code{reverse=FALSE} for then the cutpoints must be an decreasing sequence.

  }
  \item{mv}{
  Logical.
  Multivariate response? If \code{TRUE} then the input should be
  a matrix with values \eqn{1,2,\dots,L}, where \eqn{L} is the
  number of levels.
  Each column of the matrix is a response, i.e., multivariate response.
  A suitable matrix can be obtained from \code{Cut}.

  }
  \item{intercept.apply}{
  Logical.
  Whether the \code{parallel} argument should be applied to the intercept term.
  This should be set to \code{TRUE} for \code{link=}
  \code{\link{golf}},
  \code{\link{polf}},
  \code{\link{nbolf}}.

  }
}
\details{
  By default, the non-parallel cumulative logit model is fitted, i.e.,
  \deqn{\eta_j = logit(P[Y \leq j])}{%
         eta_j = logit(P[Y<=j])}
       where \eqn{j=1,2,\dots,M} and
       the \eqn{\eta_j}{eta_j} are not constrained to be parallel.
  This is also known as the \emph{non-proportional odds model}.
  If the logit link is replaced by a complementary log-log link 
  (\code{\link{cloglog}}) then
  this is known as the \emph{proportional-hazards model}.

  In almost all the literature, the constraint matrices associated
  with this family of models are known. For example, setting
  \code{parallel=TRUE} will make all constraint matrices (except for
  the intercept) equal to a vector of \eqn{M} 1's.
  If the constraint matrices are equal, unknown and to be estimated, then
  this can be achieved by fitting the model as a
  reduced-rank vector generalized
  linear model (RR-VGLM; see \code{\link{rrvglm}}).
  Currently, reduced-rank vector generalized additive models
  (RR-VGAMs) have not been implemented here.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.

}
\references{
Agresti, A. (2002)
\emph{Categorical Data Analysis},
2nd ed. New York: Wiley.

Dobson, A. J. (2001)
\emph{An Introduction to Generalized Linear Models},
2nd ed. Boca Raton: Chapman & Hall/CRC Press.

McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.

Simonoff, J. S. (2003)
\emph{Analyzing Categorical Data},
New York: Springer-Verlag.

Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.

}
\author{ Thomas W. Yee }
\note{
  The response should be either a matrix of counts (with row sums that
  are all positive), or a factor. In both cases, the \code{y} slot
  returned by \code{vglm}/\code{vgam}/\code{rrvglm} is the matrix
  of counts.

  For a nominal (unordered) factor response, the multinomial
  logit model (\code{\link{multinomial}}) is more appropriate.

  With the logit link, setting \code{parallel=TRUE} will fit a
  proportional odds model. Note that the \code{TRUE} here does
  not apply to the intercept term. 
  In practice, the validity of the proportional odds
  assumption needs to be checked, e.g., by a likelihood ratio test.
  If acceptable on the data,
  then numerical problems are less likely to occur during the fitting,
  and there are less parameters. Numerical problems occur when
  the linear/additive predictors cross, which results in probabilities
  outside of (0,1); setting \code{parallel=TRUE} will help avoid
  this problem.

  Here is an example of the usage of the \code{parallel} argument.
  If there are covariates \code{x1}, \code{x2} and \code{x3}, then
  \code{parallel = TRUE ~ x1 + x2 -1} and
  \code{parallel = FALSE ~ x3} are equivalent. This would constrain
  the regression coefficients for \code{x1} and \code{x2} to be
  equal; those of the intercepts and \code{x3} would be different.

  In the future, this family function may be renamed to
  ``\code{cups}'' (for \bold{cu}mulative \bold{p}robabilitie\bold{s})
  or ``\code{cute}'' (for \bold{cu}mulative probabili\bold{t}i\bold{e}s).

% Please let me know if you strongly agree or disagree about this.

}
\section{Warning }{
  No check is made to verify that the response is ordinal;
  see \code{\link[base:factor]{ordered}}.

}

\seealso{
  \code{\link{acat}},
  \code{\link{cratio}},
  \code{\link{sratio}},
  \code{\link{multinomial}},
  \code{\link{pneumo}},
  \code{\link{logit}},
  \code{\link{probit}},
  \code{\link{cloglog}},
  \code{\link{cauchit}},
  \code{\link{golf}},
  \code{\link{polf}},
  \code{\link{nbolf}}.

}
\examples{
# Fit the proportional odds model, p.179, in McCullagh and Nelder (1989)
data(pneumo)
pneumo = transform(pneumo, let=log(exposure.time))
(fit = vglm(cbind(normal, mild, severe) ~ let,
            cumulative(parallel=TRUE, reverse=TRUE), pneumo))
fit@y   # Sample proportions
weights(fit, type="prior")   # Number of observations
coef(fit, matrix=TRUE)
constraints(fit)   # Constraint matrices

# Check that the model is linear in let
fit2 = vgam(cbind(normal, mild, severe) ~ s(let, df=2),
            cumulative(reverse=TRUE), pneumo)
\dontrun{
plot(fit2, se=TRUE, overlay=TRUE, lcol=1:2, scol=1:2)
}

# Check the proportional odds assumption with a likelihood ratio test
(fit3 = vglm(cbind(normal, mild, severe) ~ let,
             cumulative(parallel=FALSE, reverse=TRUE), pneumo))
1 - pchisq(2*(logLik(fit3)-logLik(fit)),
           df=length(coef(fit3))-length(coef(fit)))
}
\keyword{models}
\keyword{regression}

% pneumo$let = log(pneumo$exposure.time)

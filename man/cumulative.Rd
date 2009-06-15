\name{cumulative}
\alias{cumulative}
%\alias{scumulative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Ordinal Regression with Cumulative Probabilities }
\description{
  Fits a cumulative link
  regression model to a (preferably ordered) factor response. 

}
\usage{
cumulative(link = "logit", earg = list(),
           parallel = FALSE, reverse = FALSE,
           mv = FALSE, intercept.apply = FALSE)
}
%scumulative(link="logit", earg = list(),
%            lscale="loge", escale = list(),
%            parallel=FALSE, sparallel=TRUE, reverse=FALSE, iscale = 1)
%- maybe also 'usage' for other objects documented here.
\arguments{
  In the following, the response \eqn{Y} is assumed to be a factor
  with ordered values \eqn{1,2,\dots,J+1}.
  Hence \eqn{M} is the number of linear/additive predictors
  \eqn{\eta_j}{eta_j};
  for \code{cumulative()} \eqn{M=J}.
% and for \code{scumulative()} \eqn{M=2J}.

  \item{link}{
  Link function applied to the \eqn{J} cumulative probabilities. 
  See \code{\link{Links}} for more choices,
  e.g., for the cumulative
  \code{\link{probit}}/\code{\link{cloglog}}/\code{\link{cauchit}}/\ldots
  models.

  }
% \item{lscale}{
% Link function applied to the \eqn{J} scaling parameters.
% See \code{\link{Links}} for more choices.
%
% }
  \item{earg}{
  List. Extra argument for the link function.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{parallel}{
  A logical  or formula specifying which terms have
  equal/unequal coefficients.
  See below for more information about the parallelism assumption.

  }
% \item{sparallel}{
% For the scaling parameters.
% A logical, or formula specifying which terms have
% equal/unequal coefficients.
% This argument is not applied to the intercept.
% The \code{scumulative()} function requires covariates; for
% intercept models use \code{cumulative()}.

% }
  \item{reverse}{
  Logical.
  By default, the cumulative probabilities used are
  \eqn{P(Y\leq 1)}{P(Y<=1)}, \eqn{P(Y\leq 2)}{P(Y<=2)},
  \dots, \eqn{P(Y\leq J)}{P(Y<=J)}.
  If \code{reverse} is \code{TRUE}, then 
  \eqn{P(Y\geq 2)}{P(Y>=2)}, \eqn{P(Y\geq 3)}{P(Y>=3)}, \dots,
  \eqn{P(Y\geq J+1)}{P(Y>=J+1)} will be used.

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
  a matrix with values \eqn{1,2,\dots,L}, where \eqn{L=J+1} is the
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
% \item{iscale}{
% Numeric. Initial values for the scale parameters.

% }
}
\details{
  This \pkg{VGAM} family function fits the class of
  cumulative link models to (hopefully) an ordinal response.
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
  \code{parallel=TRUE} will make all constraint matrices
  (except for the intercept) equal to a vector of \eqn{M} 1's.
  If the constraint matrices are equal, unknown and to be estimated, then
  this can be achieved by fitting the model as a
  reduced-rank vector generalized
  linear model (RR-VGLM; see \code{\link{rrvglm}}).
  Currently, reduced-rank vector generalized additive models
  (RR-VGAMs) have not been implemented here.

% The scaled version of \code{cumulative()}, called \code{scumulative()},
% has \eqn{J} positive scaling factors.
% They are described in pages 154 and 177 of McCullagh and Nelder (1989);
% see their equation (5.4) in particular,
% which they call the \emph{generalized rational model}.

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
  returned by \code{\link{vglm}}/\code{\link{vgam}}/\code{\link{rrvglm}}
  is the matrix
  of counts.
  The formula must contain an intercept term.
  Other \pkg{VGAM} family functions for an ordinal response include
  \code{\link{acat}},
  \code{\link{cratio}},
  \code{\link{sratio}}.

  For a nominal (unordered) factor response, the multinomial
  logit model (\code{\link{multinomial}}) is more appropriate.

  With the logit link, setting \code{parallel=TRUE} will fit a
  proportional odds model. Note that the \code{TRUE} here does
  not apply to the intercept term. 
  In practice, the validity of the proportional odds assumption
  needs to be checked, e.g., by a likelihood ratio test (LRT).
  If acceptable on the data,
  then numerical problems are less likely to occur during the fitting,
  and there are less parameters. Numerical problems occur when
  the linear/additive predictors cross, which results in probabilities
  outside of \eqn{(0,1)}; setting \code{parallel=TRUE} will help avoid
  this problem.

  Here is an example of the usage of the \code{parallel} argument.
  If there are covariates \code{x2}, \code{x3} and \code{x4}, then
  \code{parallel = TRUE ~ x2 + x3 -1} and
  \code{parallel = FALSE ~ x4} are equivalent. This would constrain
  the regression coefficients for \code{x2} and \code{x3} to be
  equal; those of the intercepts and \code{x4} would be different.

  If the data is inputted in \emph{long} format
  (not \emph{wide} format, as in \code{\link{pneumo}} below)
  and the self-starting initial values are not good enough then
  try using
  \code{mustart},
  \code{coefstart} and/or 
  \code{etatstart}.
  See the example below.


% In the future, this family function may be renamed to
% ``\code{cups}'' (for \bold{cu}mulative \bold{p}robabilitie\bold{s})
% or ``\code{cute}'' (for \bold{cu}mulative probabili\bold{t}i\bold{e}s).

% Please let me know if you strongly agree or disagree about this.

}
\section{Warning }{
  No check is made to verify that the response is ordinal;
  see \code{\link[base:factor]{ordered}}.

}

\seealso{
  \code{\link{margeff}},
  \code{\link{acat}},
  \code{\link{cratio}},
  \code{\link{sratio}},
  \code{\link{multinomial}},
  \code{\link{pneumo}},
  \code{\link{Links}},
  \code{\link{logit}},
  \code{\link{probit}},
  \code{\link{cloglog}},
  \code{\link{cauchit}},
  \code{\link{golf}},
  \code{\link{polf}},
  \code{\link{nbolf}},
  \code{\link{logistic1}}.

}
\examples{
# Fit the proportional odds model, p.179, in McCullagh and Nelder (1989)
pneumo = transform(pneumo, let=log(exposure.time))
(fit = vglm(cbind(normal, mild, severe) ~ let,
            cumulative(parallel=TRUE, reverse=TRUE), pneumo))
fit@y   # Sample proportions
weights(fit, type="prior")   # Number of observations
coef(fit, matrix=TRUE)
constraints(fit)   # Constraint matrices

# Check that the model is linear in let ----------------------
fit2 = vgam(cbind(normal, mild, severe) ~ s(let, df=2),
            cumulative(reverse=TRUE), pneumo)
\dontrun{
plot(fit2, se=TRUE, overlay=TRUE, lcol=1:2, scol=1:2)
}

# Check the proportional odds assumption with a LRT ----------
(fit3 = vglm(cbind(normal, mild, severe) ~ let,
             cumulative(parallel=FALSE, reverse=TRUE), pneumo))
1 - pchisq(2*(logLik(fit3)-logLik(fit)),
           df=length(coef(fit3))-length(coef(fit)))

# A factor() version of fit ----------------------------------
# This is in long format (cf. wide format above)
nobs = round(fit@y * c(weights(fit, type="prior")))
sumnobs = colSums(nobs) # apply(nobs, 2, sum)

pneumo.long = data.frame(symptoms=ordered(rep(rep(colnames(nobs),
                                  nrow(nobs)),
                                  times=c(t(nobs))),
                                  levels = colnames(nobs)),
                         let = rep(rep(with(pneumo, let), each=ncol(nobs)),
                                   times=c(t(nobs))))
with(pneumo.long, table(let, symptoms)) # check it; should be same as pneumo


(fit.long1 = vglm(symptoms ~ let, data=pneumo.long,
             cumulative(parallel=TRUE, reverse=TRUE), trace=TRUE))
coef(fit.long1, matrix=TRUE) # Should be same as coef(fit, matrix=TRUE)
# Could try using mustart if fit.long1 failed to converge.
mymustart = matrix(sumnobs / sum(sumnobs),
                   nrow(pneumo.long), ncol(nobs), byrow=TRUE)
fit.long2 = vglm(symptoms ~ let,
                fam = cumulative(parallel=TRUE, reverse=TRUE),
                mustart=mymustart, data = pneumo.long, trace=TRUE)
coef(fit.long2, matrix=TRUE) # Should be same as coef(fit, matrix=TRUE)
}
\keyword{models}
\keyword{regression}

% pneumo$let = log(pneumo$exposure.time)

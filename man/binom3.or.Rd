\name{binom3.or}
\alias{binom3.or}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Trivariate Binary Regression with Three
        Odds Ratios (Family Function) }

\description{
  Fits three Palmgren (bivariate odds-ratio model, or
  bivariate logistic regression) models simultaneously
  to three binary responses.
  The odds ratios are used to measure dependencies
  between the responses.
  Several options for the joint probabilities are
  available.


% Actually, a trivariate logistic/probit/cloglog/cauchit
% model can be fitted.

  
}
\usage{
binom3.or(lmu = "logitlink", lmu1 = lmu, lmu2 = lmu,
  lmu3 = lmu, loratio = "loglink", zero = "oratio",
  exchangeable = FALSE, eq.or = FALSE, jpmethod =
  c("min", "mean", "median", "max", "1", "2", "3"),
  imu1 = NULL, imu2 = NULL, imu3 = NULL,
  ioratio12 = NULL, ioratio13 = NULL, ioratio23 = NULL,
  tol = 0.001, more.robust = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu}{
  Same as \code{\link{binom2.or}}.


  }
  \item{lmu1, lmu2, lmu3}{
  Same as \code{\link{binom2.or}}.


  }
  \item{loratio}{
  Same as \code{\link{binom2.or}}.
  Applied to all three odds ratios,
  called \code{oratio12}, \code{oratio13},
  \code{oratio23}.
  


  }
  \item{imu1, imu2, imu3}{
  Similar to \code{\link{binom2.or}}.


  }
  \item{ioratio12, ioratio13, ioratio23}{
  Similar to \code{\link{binom2.or}}.


  }
  \item{zero, exchangeable}{
  Same as \code{\link{binom2.or}}.


  }
  \item{eq.or}{
    Logical. Constrain all the odds ratios to
    be equal? Setting \code{exchangeable = TRUE}
    implies that this is \code{TRUE} also.
    Setting \code{eq.or = TRUE} sometimes is
    a good way to obtain a more stable model,
    because too many different odds ratios
    can easily create numerical problems,
    especially if \code{zero = NULL}.


  }
  \item{jpmethod}{
    See \code{\link{dbinom3.or}}.


  }
  

  \item{tol, more.robust}{
  Same as \code{\link{binom2.or}}.


  }
}
\details{
  This heuristic model is an extension of
  \code{\link{binom2.or}} for handling \emph{three}
  binary responses.
  Rather than allowing something like
  \code{vglm(cbind(y1,y2, y1,y3, y2,y3) ~ x2, binom2.or)},
  which has three pairs of bivariate responses
  in the usual form of multiple responses
  allowed in \pkg{VGAM},
  I have decided to write \code{\link{binom3.or}}
  which operates on \code{cbind(y1, y2, y3)}
  instead.
  This model thus uses three odds ratios to allow
  for dependencies between the pairs of responses.
  It is heuristic because the joint probability
  \eqn{P(y_1=1,y_2=1,y_3=1)=p_{123}}{P(y1=1,y2=1,y3=1)=p123}
  is computed by a number of
  conditional independence assumptions.


  
  This trivariate logistic model has a
  fully specified likelihood.
  Explicitly,
  the default model is
\deqn{logit\;P(Y_j=1)] = \eta_j,\ \ \ j=1,2}{%
      logit  P(Y_j=1)] =  eta_j,\ \ \ j=1,2}
    for the first two marginals,
\deqn{logit\; P(Y_1=1)  = \eta_4,}{%
      logit   P(Y_1=1)  =  eta_4,}
\deqn{logit\; P(Y_3=1)  = \eta_5,}{%
      logit   P(Y_3=1)  =  eta_5,}
\deqn{logit\; P(Y_2=1)  = \eta_7,}{%
      logit   P(Y_2=1)  =  eta_7,}
\deqn{logit\; P(Y_3=1)  = \eta_8,}{%
      logit   P(Y_3=1)  =  eta_8,}
    and
\deqn{\log \psi_{12} = \eta_3,}{%
       log  psi12    =  eta_3,}
\deqn{\log \psi_{13} = \eta_6,}{%
       log  psi13    =  eta_6,}
\deqn{\log \psi_{23} = \eta_9,}{%
       log  psi23    =  eta_9,}
  specifies the dependency between each possible
  pair of responses.
  Many details on such quantities are similar to
  \code{\link{binom2.or}}.


  By default, all odds ratios are intercept-only.
  The \code{exchangeable} argument should be
  used when the error structure is exchangeable.
  However, there is a difference between
  \emph{full}
  and \emph{partial} exchangeability and
  setting \code{exchangeable = TRUE} results in
  the full version.
  The partial version would require the manual
  input of certain constraint matrices via
  \code{constraints}.


  

}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}}
  and \code{\link{vgam}}.


  When fitted, the \code{fitted.values} slot
  of the object is a matrix with successive
  columns equalling the eight joint
  probabilities, labelled
  as \eqn{(Y_1,Y_2,Y_3)}{(Y1,Y2)} =
  (0,0,0), (0,0,1), (0,1,0), (0,1,1),
  (1,0,0), (1,0,1), (1,1,0), (1,1,1),
  respectively.  These estimated probabilities
  should be extracted
  with the \code{fitted} generic function.


}
%\references{
%
%}
%\author{ Thomas W. Yee }
\note{
  At present we call \code{\link{binom3.or}}
  a \emph{trivariate odds-ratio model} (TOM).
  The response should be either a 8-column
  matrix of counts
  (whose columns correspond to
  \eqn{(Y_1,Y_2,Y_3)}{(Y1,Y2)} ordered
  as above),
  or a three-column matrix where each column
  has two distinct values, or a factor with
  8 levels.
  The function \code{\link{rbinom3.or}} may be
  used to generate such data.



% Successful convergence requires at least one
% case of each of the 8 possible outcomes zz.
% (Note: W S D M == Wald Statistic Derivative Metric)




Because some of the \eqn{\eta_j}{eta_j}
are repeated, the constraint matrices have
a special form in order to provide consistency.


  

  By default, intercept-only odds ratios
  are fitted because \code{zero = "oratio"}.
  Set \code{zero = NULL} for the odds
  ratios to be modelled as a function of the
  explanatory variables; however, numerical
  problems are more likely to occur.




  The argument \code{lmu}, which is actually
  redundant, is used for convenience and for
  upward compatibility: specifying \code{lmu}
  only means the link function will be applied
  to \code{lmu1}, \code{lmu2} and \code{lmu3}.
  Users who want a different link function for
  each of the marginal probabilities should use
  \code{lmu1}, \code{lmu2} and \code{lmu3},
  and the argument \code{lmu} is then
  ignored.  It doesn't make sense to specify
  \code{exchangeable = TRUE} and have different
  link functions for the marginal probabilities.

}
\section{Warning }{
  Because the parameter space of this model is
  restricted (see \code{\link{rbinom3.or}}),
  this family function is more limited than
  \code{\link{loglinb3}}. However, this model
  is probably more interpretable since the
  marginal probabilities and odds ratios are
  modelled by conventional link functions
  directly.


  
  If the data is very sparse then convergence
  problems will occur. It is recommended that the
  sample size is several hundred at least.
  Opinion: anything
  less than \eqn{n=100} is liable for failure.
  Setting \code{trace = TRUE} is urged.



}


\seealso{
  \code{\link{rbinom3.or}},
  \code{\link{binom2.or}},
  \code{\link{loglinb3}}.


}
\examples{  set.seed(1)
\dontrun{
nn <- 1000  # Example 1
ymat <- rbinom3.or(nn, mu1 = logitlink(0.5, inv = TRUE),
                   oratio12 = exp(1), exch = TRUE)
fit1 <- vglm(ymat ~ 1, binom3.or(exc = TRUE), tra = TRUE)
coef(fit1, matrix = TRUE)
constraints(fit1)

bdata <- data.frame(x2 = sort(runif(nn)))  # Example 2
bdata <- transform(bdata,
       mu1 = logitlink(-1 + 1 * x2, inv = TRUE),
       mu2 = logitlink(-1 + 2 * x2, inv = TRUE),
       mu3 = logitlink( 2 - 1 * x2, inv = TRUE))
ymat2 <- with(bdata,
     rbinom3.or(nn, mu1, mu2, mu3, exp(0.25),
                oratio13 = exp(0.25), exp(0.25)))
fit2 <- vglm(ymat2 ~ x2, binom3.or(eq.or = TRUE),
             bdata, trace = TRUE)
coef(fit2, matrix = TRUE)
}}
\keyword{models}
\keyword{regression}

%dmat <- with(bdata,
%     dbinom3.or(mu1, mu2, mu3, oratio12 = exp(0.25),
%                oratio13 = exp(0.25), exp(0.25)))
%
%
%













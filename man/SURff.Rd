\name{SURff}
\alias{SURff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Seemingly Unrelated Regressions Family Function
%%  ~~function to do ... ~~
}
\description{
Fits a system of seemingly unrelated regressions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SURff(mle.normal = FALSE,
      divisor = c("n", "n-max(pj,pk)", "sqrt((n-pj)*(n-pk))"),
      parallel = FALSE, Varcov = NULL, matrix.arg = FALSE)
}
%- maybe also 'usage' for other objects documented here.
%apply.parint = TRUE,
\arguments{
%  \item{estimator}{
%Character.
%What estimator is computed.
%%     ~~Describe \code{estimator} here~~
%}

  \item{mle.normal}{
  Logical.
  If \code{TRUE} then the MLE, assuming multivariate normal errors,
  is computed;
  the effect is just to add a \code{loglikelihood} slot to the
  returned object.
  Then it results in the \emph{maximum likelihood estimator}.


}


  \item{divisor}{
Character, partial matching allowed and the first choice is the default.
The divisor for the estimate of the covariances.
If \code{"n"} then the estimate will be biased.
If the others then the estimate will be unbiased for some elements.
If \code{mle.normal = TRUE} and this argument is not \code{"n"} then
a warning or an error will result.


}
  \item{parallel}{
  See
  \code{\link{CommonVGAMffArguments}}.
  If \code{parallel = TRUE} then the constraint applies to
  the intercept too.


}
  \item{Varcov}{
  Numeric.
  This may be assigned a variance-covariance of the errors.
  If \code{matrix.arg} then this is a \eqn{M \times M}{M x M} matrix.
  If \code{!matrix.arg} then this is a \eqn{M \times M}{M x M} matrix in
  matrix-band format (a vector with at least \eqn{M} and
  at most \code{M*(M+1)/2} elements).


}
  \item{matrix.arg}{
  Logical.
  Of single length.


}
}
\details{
  Proposed by Zellner (1962), the basic
  seemingly unrelated regressions (SUR)
  model is a set of LMs (\eqn{M > 1} of them) tied together
  at the error term level.
  Each LM's model matrix may potentially have its own set
  of predictor variables.


  Zellner's efficient (ZEF) estimator (also known as
  \emph{Zellner's two-stage Aitken estimator})
  can be obtained by setting
  \code{maxit = 1}
  (and possibly \code{divisor = "sqrt"} or
  \code{divisor = "n-max"}).


  The default value of \code{maxit} (in \code{\link{vglm.control}})
  probably means \emph{iterative GLS} (IGLS) estimator is computed because
  IRLS will probably iterate to convergence.
  IGLS means, at each iteration, the residuals are used to estimate
  the error variance-covariance matrix, and then the matrix is used
  in the GLS.
  The IGLS estimator is also known
  as \emph{Zellner's iterative Aitken estimator}, or IZEF.


}
\value{

  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as
  \code{\link{vglm}} and \code{\link{vgam}}.


}
\references{

  Zellner, A. (1962)
  An Efficient Method of Estimating Seemingly Unrelated
  Regressions and Tests for Aggregation Bias.
  \emph{J. Amer. Statist. Assoc.},
  \bold{57}(298), 348--368.


  Kmenta, J. and Gilbert, R. F. (1968)
  Small Sample Properties of Alternative Estimators
  of Seemingly Unrelated Regressions.
  \emph{J. Amer. Statist. Assoc.},
  \bold{63}(324), 1180--1200.


}
\author{

  T. W. Yee.


}
\section{Warning }{
  The default convergence criterion may be a little loose.
  Try setting \code{epsilon = 1e-11}, especially
  with \code{mle.normal =  TRUE}.


}
\note{
  The fitted object has slot \code{@extra$ncols.X.lm} which is
  a \eqn{M} vector with the number of parameters for each LM.
  Also, \code{@misc$values.divisor} is the \eqn{M}-vector of
  \code{divisor} values.


  Constraint matrices are needed in order to specify which response
  variables that each term on the RHS of the formula is a
  regressor for.
  See the \code{constraints} argument of \code{\link{vglm}}
  for more information.


% This \pkg{VGAM} family function is currently experimental.


}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{uninormal}},
  \code{\link{gew}}.


}
\examples{
# Obtain some of the results of p.1199 of Kmenta and Gilbert (1968)
clist <- list("(Intercept)" = diag(2),
              "capital.g"   = rbind(1, 0),
              "value.g"     = rbind(1, 0),
              "capital.w"   = rbind(0, 1),
              "value.w"     = rbind(0, 1))
zef1 <- vglm(cbind(invest.g, invest.w) ~
             capital.g + value.g + capital.w + value.w,
             SURff(divisor = "sqrt"), maxit = 1,
             data = gew, trace = TRUE, constraints = clist)

round(coef(zef1, matrix = TRUE), digits = 4)  # ZEF
zef1@extra$ncols.X.lm
zef1@misc$divisor
zef1@misc$values.divisor
round(sqrt(diag(vcov(zef1))),    digits = 4)  # SEs
nobs(zef1, type = "lm")
df.residual(zef1, type = "lm")


mle1 <- vglm(cbind(invest.g, invest.w) ~
             capital.g + value.g + capital.w + value.w,
             SURff(mle.normal = TRUE),
             epsilon = 1e-11,
             data = gew, trace = TRUE, constraints = clist)
round(coef(mle1, matrix = TRUE), digits = 4)  # MLE
round(sqrt(diag(vcov(mle1))),    digits = 4)  # SEs
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}



% Prior to 20141108:
%            SURff(mle.normal = TRUE, divisor = "n-max"),


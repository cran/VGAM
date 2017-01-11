\name{Qvar}
\alias{Qvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quasi-variances Preprocessing Function

%%  ~~function to do ... ~~
}
\description{
  Takes a \code{\link{vglm}} fit or a variance-covariance matrix,
  and preprocesses it for \code{\link{rcim}} and
  \code{\link{uninormal}} so that quasi-variances can be computed.


%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Qvar(object, factorname = NULL, which.linpred = 1,
     coef.indices = NULL, labels = NULL,
     dispersion = NULL, reference.name = "(reference)", estimates = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A \code{"\link[=vglmff-class]{vglm}"} object or a variance-covariance
  matrix, e.g., \code{vcov(vglm.object)}.
  The former is preferred since it contains all the information needed.
  If a matrix then \code{factorname} and/or \code{coef.indices}
  should be specified to identify the factor.


}
\item{which.linpred}{
  A single integer from the set \code{1:M}.
  Specifies which linear predictor to use.
  Let the value of \code{which.linpred} be called \eqn{j}.
  Then the factor should appear in that linear predictor, hence
  the \eqn{j}th row of the constraint matrix corresponding
  to the factor should have at least one nonzero value.
  Currently the \eqn{j}th row must have exactly one nonzero value
  because programming it for more than one nonzero value is difficult.


}

\item{factorname}{
  Character.
  If the \code{\link{vglm}} object contains more than one
  factor as explanatory variable then this argument should
  be the name of the factor of interest.
  If \code{object} is a variance-covariance matrix then
  this argument should also be specified.


}
\item{labels}{
  Character.
  Optional, for labelling the variance-covariance matrix.


}
\item{dispersion}{
  Numeric.
  Optional, passed into \code{vcov()} with the same argument name.


}
\item{reference.name}{
  Character.
  Label for for the reference level.


}
\item{coef.indices}{
  Optional numeric vector of length at least 3 specifying
  the indices of the factor from the variance-covariance
  matrix.


}
\item{estimates}{
  an optional vector of estimated coefficients
  (redundant if \code{object} is a model).


}
}
\details{


  Suppose a factor with \eqn{L} levels is an explanatory variable in a
  regression model. By default, R treats the first level as baseline so
  that its coefficient is set to zero. It estimates the other \eqn{L-1}
  coefficients, and with its associated standard errors, this is the
  conventional output. From the complete variance-covariance matrix one
  can compute \eqn{L} quasi-variances based on all pairwise difference
  of the coefficients. They are based on an approximation, and can be
  treated as uncorrelated. In minimizing the relative (not absolute)
  errors it is not hard to see that the estimation involves a RCIM
  (\code{\link{rcim}}) with an exponential link function
  (\code{\link{explink}}).


  If \code{object} is a model, then at least one of \code{factorname} or
  \code{coef.indices} must be non-\code{NULL}. The value of
  \code{coef.indices}, if non-\code{NULL}, determines which rows and
  columns of the model's variance-covariance matrix to use. If
  \code{coef.indices} contains a zero, an extra row and column are
  included at the indicated position, to represent the zero variances
  and covariances associated with a reference level. If
  \code{coef.indices} is \code{NULL}, then \code{factorname} should be
  the name of a factor effect in the model, and is used in order to
  extract the necessary variance-covariance estimates.


  Quasi-variances were first implemented in R with \pkg{qvcalc}.
  This implementation draws heavily from that.


}
\value{
  A \eqn{L} by \eqn{L} matrix whose \eqn{i}-\eqn{j} element
  is the logarithm of the variance of the \eqn{i}th coefficient
  minus the \eqn{j}th coefficient, for all values of \eqn{i}
  and \eqn{j}. The diagonal elements are abitrary and are set
  to zero.


  The matrix has an attribute that corresponds to the prior
  weight matrix; it is accessed by \code{\link{uninormal}}
  and replaces the usual \code{weights} argument.
  of \code{\link{vglm}}. This weight matrix has ones on
  the off-diagonals and some small positive number on
  the diagonals.


}
\references{


  Firth, D. (2003)
  Overcoming the reference category problem in the
  presentation of statistical models.
  \emph{Sociological Methodology} \bold{33}, 1--18.


  Firth, D. and de Menezes, R. X. (2004)
  Quasi-variances.
  \emph{Biometrika} \bold{91}, 65--80.


Yee, T. W. and Hadi, A. F. (2014)
Row-column interaction models, with an R implementation.
\emph{Computational Statistics},
\bold{29}, 1427--1445.


}

\author{


  T. W. Yee, based heavily on \code{qvcalc()} in \pkg{qvcalc}
  written by David Firth.


}

\note{
  This is an adaptation of \code{qvcalc()} in \pkg{qvcalc}.
  It should work for all \code{\link{vglm}}
  models with one linear predictor, i.e., \eqn{M = 1}.
  For \eqn{M > 1} the factor should appear only in one of the
  linear predictors.


  It is important to set \code{maxit} to be larger than usual for
  \code{\link{rcim}} since convergence is slow. Upon successful
  convergence the \eqn{i}th row effect and the \eqn{i}th column effect
  should be equal. A simple computation involving the fitted and
  predicted values allows the quasi-variances to be extracted (see
  example below).


  A function to plot \emph{comparison intervals} has not been
  written here.

}

\section{Warning }{
  Negative quasi-variances may occur (one of them and
  only one), though they are rare in practice. If
  so then numerical problems may occur. See
  \code{qvcalc()} for more information.



}


\seealso{
  \code{\link{rcim}},
  \code{\link{vglm}},
  \code{\link{qvar}},
  \code{\link{uninormal}},
  \code{\link{explink}},
  \code{qvcalc()} in \pkg{qvcalc},
  \code{\link[MASS]{ships}}.


%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Example 1
data("ships", package = "MASS")

Shipmodel <- vglm(incidents ~ type + year + period,
                  quasipoissonff, offset = log(service),
#                 trace = TRUE, model = TRUE,
                  data = ships, subset = (service > 0))

# Easiest form of input
fit1 <- rcim(Qvar(Shipmodel, "type"), uninormal("explink"), maxit = 99)
qvar(fit1)              # Easy method to get the quasi-variances
qvar(fit1, se = TRUE)   # Easy method to get the quasi-standard errors

(quasiVar <- exp(diag(fitted(fit1))) / 2)                 # Version 1
(quasiVar <- diag(predict(fit1)[, c(TRUE, FALSE)]) / 2)   # Version 2
(quasiSE  <- sqrt(quasiVar))

# Another form of input
fit2 <- rcim(Qvar(Shipmodel, coef.ind = c(0, 2:5), reference.name = "typeA"),
             uninormal("explink"), maxit = 99)
\dontrun{ qvplot(fit2, col = "green", lwd = 3, scol = "blue", slwd = 2, las = 1) }

# The variance-covariance matrix is another form of input (not recommended)
fit3 <- rcim(Qvar(cbind(0, rbind(0, vcov(Shipmodel)[2:5, 2:5])),
                  labels = c("typeA", "typeB", "typeC", "typeD", "typeE"),
                  estimates = c(typeA = 0, coef(Shipmodel)[2:5])),
             uninormal("explink"), maxit = 99)
(QuasiVar <- exp(diag(fitted(fit3))) / 2)                 # Version 1
(QuasiVar <- diag(predict(fit3)[, c(TRUE, FALSE)]) / 2)   # Version 2
(QuasiSE  <- sqrt(quasiVar))
\dontrun{ qvplot(fit3) }


# Example 2: a model with M > 1 linear predictors
\dontrun{ require("VGAMdata")
xs.nz.f <- subset(xs.nz, sex == "F")
xs.nz.f <- subset(xs.nz.f, !is.na(babies)  & !is.na(age) & !is.na(ethnicity))
xs.nz.f <- subset(xs.nz.f, ethnicity != "Other")

clist <- list("sm.bs(age, df = 4)" = rbind(1, 0),
              "sm.bs(age, df = 3)" = rbind(0, 1),
              "ethnicity"          = diag(2),
              "(Intercept)"        = diag(2))
fit1 <- vglm(babies ~ sm.bs(age, df = 4) + sm.bs(age, df = 3) + ethnicity,
            zipoissonff(zero = NULL), xs.nz.f,
            constraints = clist, trace = TRUE)
Fit1 <- rcim(Qvar(fit1, "ethnicity", which.linpred = 1),
             uninormal("explink", imethod = 1), maxit = 99, trace = TRUE)
Fit2 <- rcim(Qvar(fit1, "ethnicity", which.linpred = 2),
             uninormal("explink", imethod = 1), maxit = 99, trace = TRUE)
}
\dontrun{ par(mfrow = c(1, 2))
qvplot(Fit1, scol = "blue", pch = 16, main = expression(eta[1]),
       slwd = 1.5, las = 1, length.arrows = 0.07)
qvplot(Fit2, scol = "blue", pch = 16, main = expression(eta[2]),
       slwd = 1.5, las = 1, length.arrows = 0.07)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
% \code{\link[qvcalc:qvcalc]{qvcalc}} in \pkg{qvcalc}

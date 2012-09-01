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
  \code{\link{normal1}} so that quasi-variances can be computed.


%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Qvar(object, factorname = NULL, coef.indices = NULL, labels = NULL,
     dispersion = NULL, reference.name = "(reference)", estimates = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \code{"\link[=vglmff-class]{vglm}"} object or a variance-covariance
    matrix, e.g., \code{vcov(vglm.object)}.
    The former is preferred since it contains all the information
    needed.
  If a matrix then \code{factorname} and/or \code{coef.indices}
  should be specified to identify the factor.


%%     ~~Describe \code{object} here~~
}
\item{factorname}{
  Character.
  If the \code{\link{vglm}} object contains more than one
  factor as explanatory variable then this argument should
  be the name of the factor of interest.
  If \code{object} is a variance-covariance matrix then
  this argument should also be specified.


%%     ~~Describe \code{factor.name} here~~
}
\item{labels}{
  Character.
  Optional, for labelling the variance-covariance matrix.

%%     ~~Describe \code{level1.name} here~~
}
\item{dispersion}{
  Numeric.
  Optional, passed into \code{vcov()} with the same argument name.

%%     ~~Describe \code{level1.name} here~~
}
\item{reference.name}{
  Character.
  Label for for the reference level.

%%     ~~Describe \code{level1.name} here~~
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



%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A \eqn{L} by \eqn{L} matrix whose \eqn{i}-\eqn{j} element
  is the logarithm of the variance of the \eqn{i}th coefficient
  minus the \eqn{j}th coefficient, for all values of \eqn{i}
  and \eqn{j}. The diagonal elements are abitrary and are set
  to zero.

  
  The matrix has an attribute that corresponds to the prior
  weight matrix; it is accessed by \code{\link{normal1}}
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
  \code{\link{normal1}},
  \code{\link{explink}},
  \code{qvcalc()} in \pkg{qvcalc},
  \code{\link[MASS]{ships}}.


%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("ships", package = "MASS")

Shipmodel <- vglm(incidents ~ type + year + period,
                  quasipoissonff, offset = log(service),
#                 trace = TRUE, model = TRUE,
                  data = ships, subset = (service > 0))

# Easiest form of input
fit1 <- rcim(Qvar(Shipmodel, "type"), normal1("explink"), maxit = 99)
(quasiVar <- exp(diag(fitted(fit1))) / 2)                 # Version 1
(quasiVar <- diag(predict(fit1)[, c(TRUE, FALSE)]) / 2)   # Version 2
(quasiSE  <- sqrt(quasiVar))

# Another form of input
fit2 <- rcim(Qvar(Shipmodel, coef.ind = c(0,2:5), reference.name = "typeA"),
             normal1("explink"), maxit = 99)
\dontrun{ plotqvar(fit2, col = "green", lwd = 3, scol = "blue", slwd = 2, las = 1) }

# The variance-covariance matrix is another form of input (not recommended)
fit3 <- rcim(Qvar(cbind(0, rbind(0, vcov(Shipmodel)[2:5, 2:5])),
                  labels = c("typeA", "typeB", "typeC", "typeD", "typeE"),
                  estimates = c(typeA = 0, coef(Shipmodel)[2:5])),
             normal1("explink"), maxit = 99)
(QuasiVar <- exp(diag(fitted(fit3))) / 2)                 # Version 1
(QuasiVar <- diag(predict(fit3)[, c(TRUE, FALSE)]) / 2)   # Version 2
(QuasiSE  <- sqrt(quasiVar))
\dontrun{ plotqvar(fit3) }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
% \code{\link[qvcalc:qvcalc]{qvcalc}} in \pkg{qvcalc}

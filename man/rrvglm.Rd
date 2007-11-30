\name{rrvglm}
\alias{rrvglm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Fitting Reduced-Rank Vector Generalized Linear Models (RR-VGLMs) }
\description{
  A \emph{reduced-rank vector generalized linear model} (RR-VGLM) is fitted.
  RR-VGLMs are VGLMs but some of the constraint matrices are estimated.
  In this documentation, \eqn{M} is the number of linear predictors.

}
\usage{
rrvglm(formula, family, data = list(), weights = NULL, subset = NULL,
       na.action = na.fail, etastart = NULL, mustart = NULL,
       coefstart = NULL, control = rrvglm.control(...), offset = NULL,
       method = "rrvglm.fit", model = FALSE, x.arg = TRUE, y.arg = TRUE,
       contrasts = NULL, constraints = NULL, extra = NULL,
       qr.arg = FALSE, smart = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{

  \item{formula}{
  a symbolic description of the model to be fit.
  The RHS of the formula is applied to each linear predictor. Different
  variables in each linear predictor can be chosen by specifying
  constraint matrices. 

  }
  \item{family}{
  a function of class \code{"vglmff"} (see \code{\link{vglmff-class}})
  describing what statistical model is to be fitted. This is called a
  ``\pkg{VGAM} family function''.  See \code{\link{CommonVGAMffArguments}}
  for general information about many types of arguments found in this
  type of function.

  }
  \item{data}{
  an optional data frame containing the variables in the model.
  By default the variables are taken from \code{environment(formula)},
  typically the environment from which \code{rrvglm} is called.

  }
  \item{weights}{
  an optional vector or matrix of (prior) weights 
  to be used in the fitting process.
  If \code{weights} is a matrix, then it must be in
  \emph{matrix-band} form, whereby the first \eqn{M} 
  columns of the matrix are the
  diagonals, followed by the upper-diagonal band, followed by the
  band above that, etc. In this case, there can be up to \eqn{M(M+1)}
  columns, with the last column corresponding to the
  (1,\eqn{M}) elements of the weight matrices.

  }
  \item{subset}{
  an optional logical vector specifying a subset of observations to be
  used in the fitting process.

  }
  \item{na.action}{
  a function which indicates what should happen when the data contain
  \code{NA}s.
  The default is set by the \code{na.action} setting
  of \code{\link[base]{options}}, and is \code{na.fail} if that is unset.
  The ``factory-fresh'' default is \code{na.omit}.

  }
  \item{etastart}{
  starting values for the linear predictors.
  It is a \eqn{M}-column matrix.
  If \eqn{M=1} then it may be a vector.
  }
  \item{mustart}{
  starting values for the fitted values. It can be a vector or a matrix.
  Some family functions do not make use of this argument.

  }
  \item{coefstart}{
  starting values for the coefficient vector.

  }
  \item{control}{
  a list of parameters for controlling the fitting process. 
  See \code{\link{rrvglm.control}} for details.

  }
  \item{offset}{
  a vector or \eqn{M}-column matrix of offset values.
  These are \emph{a priori} known and are
  added to the linear predictors during fitting.

  }
  \item{method}{
  the method to be used in fitting the model.
  The default (and presently only) method \code{rrvglm.fit}
  uses iteratively reweighted least squares (IRLS).

  }
  \item{model}{
  a logical value indicating whether the \emph{model frame}
  should be assigned in the \code{model} slot.

  }
  \item{x.arg, y.arg}{
  logical values indicating whether
  the model matrix and response vector/matrix used in the fitting
  process should be assigned in the \code{x} and \code{y} slots.
  Note the model matrix is the LM model matrix; to get the VGLM
  model matrix type \code{model.matrix(vglmfit)} where
  \code{vglmfit} is a \code{vglm} object. 

  }
  \item{contrasts}{
  an optional list. See the \code{contrasts.arg}
  of \code{\link{model.matrix.default}}.

  }
  \item{constraints}{
  an optional list  of constraint matrices.
  The components of the list must be named with the term it
  corresponds to (and it must match in character format). 
  Each constraint matrix must have \eqn{M} rows, and be of
  full-column rank.
  By default, constraint matrices are the \eqn{M} by \eqn{M}
  identity
  matrix unless arguments in the family function itself override
  these values. 
  If \code{constraints} is used it must contain \emph{all} the
  terms; an incomplete list is not accepted.

  }
  \item{extra}{
  an optional list with any extra information that might be needed
  by the family function.

  }
  \item{qr.arg}{
  logical value indicating whether
  the slot \code{qr}, which returns the QR decomposition of the
  VLM model matrix, is returned on the object.

  }
  \item{smart}{
  logical value indicating whether smart prediction
  (\code{\link{smartpred}}) will be used.

  }
  \item{\dots}{
  further arguments passed into \code{\link{rrvglm.control}}.

  }
}
\details{
  The central formula is given by
  \deqn{\eta = B_1^T x_1 + A \nu}{%
         eta = B_1^T x_1 + A nu}
  where \eqn{x_1}{x1} is a vector (usually just a 1 for an intercept),
  \eqn{x_2}{x2} is another vector of explanatory variables, and
  \eqn{\nu=C^T x_2}{nu=C^T x_2} is an \eqn{R}-vector of
  latent variables.
  Here, \eqn{\eta}{eta} is a vector of linear predictors, e.g., the
  \eqn{m}th element is \eqn{\eta_m = \log(E[Y_m])}{eta_m = log(E[Y_m])}
  for the \eqn{m}th Poisson response.  The matrices \eqn{B_1}, \eqn{A}
  and \eqn{C} are estimated from the data, i.e., contain the regression
  coefficients.  For ecologists, the central formula represents a
  \emph{constrained linear ordination} (CLO) since it is linear in
  the latent variables. It means that the response is a monotonically
  increasing or decreasing function of the latent variables.

  The underlying algorithm of RR-VGLMs is iteratively reweighted least
  squares (IRLS) with an optimizing algorithm applied within each IRLS
  iteration (e.g., alternating algorithm).

  In theory, any \pkg{VGAM} family function that works for
  \code{\link{vglm}} and \code{\link{vgam}} should work for \code{rrvglm}
  too.

  \code{rrvglm.fit} is the function that actually does the work. It is
  \code{vglm.fit} with some extra code.

}
\value{
  An object of class \code{"rrvglm"}, which has the the same slots as
  a \code{"vglm"} object. The only difference is that the some of the
  constraint matrices are estimates rather than known. But \pkg{VGAM}
  stores the models the same internally. The slots of \code{"vglm"}
  objects are described in \code{\link{vglm-class}}.

}
\references{

Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.

Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.

Anderson, J. A. (1984)
Regression and ordered categorical variables.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{46}, 1--30.

  Documentation accompanying the \pkg{VGAM} package at
  \url{http://www.stat.auckland.ac.nz/~yee}
  contains further information and examples.

}

\author{ Thomas W. Yee }
\note{
  The smart prediction (\code{\link{smartpred}}) library is packed with
  the \pkg{VGAM} library.

  The arguments of \code{rrvglm} are the same as those of
  \code{\link{vglm}} but with some extras in \code{\link{rrvglm.control}}.

  In the example below, a rank-1 stereotype model of Anderson (1984)
  is fitted to some car data.  The reduced-rank regression is performed,
  adjusting for two covariates. Setting a trivial constraint matrix for
  the latent variable variables in \eqn{x_2}{x2} avoids a warning message
  when it is overwritten by a (common) estimated constraint matrix.
  It shows that German cars tend to be more expensive than American cars,
  given a car of fixed weight and width.

  If \code{fit <- rrvglm(..., data=mydata)} then \code{summary(fit)}
  requires corner constraints and no missing values in \code{mydata}.
  Often the estimated variance-covariance matrix of the parameters is
  not positive-definite; if this occurs, try refitting the model with
  a different value for \code{Index.corner}.

  For \emph{constrained quadratic ordination} (CQO) see \code{\link{cqo}}
  for more details about QRR-VGLMs.

  With multivariate binary responses, one must use
  \code{binomialff(mv=TRUE)} to indicate that the response (matrix)
  is multivariate. Otherwise, it is interpreted as a single binary
  response variable.

}

% zzz; arguments of \code{\link{vglm}} are definitive. They're copied here.

\seealso{
    \code{\link{rrvglm.control}},
%   \code{\link{qrrvglm.control}},
    \code{\link{lvplot.rrvglm}}
    (same as \code{\link{biplot.rrvglm}}),
%   \code{\link{vcovqrrvglm}},
    \code{\link{rrvglm-class}},
    \code{\link{grc}},
    \code{\link{cqo}},
    \code{\link{vglmff-class}},
    \code{\link{vglm}},
    \code{\link{vglm-class}},
    \code{\link{smartpred}},
    \code{rrvglm.fit}.
    Methods functions include
    \code{\link{Coef.rrvglm}},
    \code{summary.rrvglm},
    etc.
}

\examples{
data(car.all)
attach(car.all)
index = Country == "Germany" | Country == "USA" |
        Country == "Japan" | Country == "Korea"
detach(car.all)
scar = car.all[index, ]  # standardized car data
fcols = c(13,14,18:20,22:26,29:31,33,34,36)  # These are factors
scar[,-fcols] = scale(scar[,-fcols]) # Standardize all numerical vars
ones = matrix(1, 3, 1)
cms = list("(Intercept)"=diag(3), Width=ones, Weight=ones,
           Disp.=diag(3), Tank=diag(3), Price=diag(3), 
           Frt.Leg.Room=diag(3))
set.seed(111)
fit = rrvglm(Country ~ Width + Weight + Disp. + Tank + Price + Frt.Leg.Room,
             multinomial, data =  scar, Rank = 2, trace = TRUE,
             constraints=cms, Norrr = ~ 1 + Width + Weight,
             Uncor=TRUE, Corner=FALSE, Bestof=2)
fit@misc$deviance  # A history of the fits
Coef(fit)
\dontrun{
biplot(fit, chull=TRUE, scores=TRUE, clty=2, ccol="blue", scol="red",
       Ccol="darkgreen", Clwd=2, Ccex=2,
       main="1=Germany, 2=Japan, 3=Korea, 4=USA")
}
}
\keyword{models}
\keyword{regression}


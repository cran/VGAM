\name{rrvglm}
\alias{rrvglm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Fitting Reduced-Rank Vector Generalized Linear Models
  (RR-VGLMs)
  and Doubly Constrained RR-VGLMs (DRR-VGLMs)
}
\description{
  A \emph{reduced-rank vector generalized linear model} (RR-VGLM)
  is fitted.  RR-VGLMs are VGLMs but some of the constraint
  matrices are estimated.
  \emph{Doubly constrained} RR-VGLMs (DRR-VGLMs)
  can also be fitted,
  and these provide structure for the two other 
  outer product matrices.

}
\usage{
rrvglm(formula, family = stop("'family' is unassigned"),
       data = list(), weights = NULL, subset = NULL,
       na.action = na.fail, etastart = NULL, mustart = NULL,
       coefstart = NULL, control = rrvglm.control(...),
       offset = NULL, method = "rrvglm.fit", model = FALSE,
       x.arg = TRUE, y.arg = TRUE, contrasts = NULL,
       constraints = NULL, extra = NULL, qr.arg = FALSE,
       smart = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{

  \item{formula, family}{
  See \code{\link{vglm}}.

  }

  \item{weights, data}{
  See \code{\link{vglm}}.


% an optional data frame containing the
%   variables in the model.
% By default the variables are taken from
% \code{environment(formula)}, typically
% the environment from
% which \code{rrvglm} is called.


  }
  \item{subset, na.action}{
  See \code{\link{vglm}}.

  }
  \item{etastart, mustart, coefstart}{
  See \code{\link{vglm}}.

  }
  \item{control}{
    a list of parameters for controlling
    the fitting process.
  See \code{\link{rrvglm.control}} for details.

  }
  \item{offset, model, contrasts}{
  See \code{\link{vglm}}.

  }
  \item{method}{
  the method to be used in fitting the model.
  The default (and presently only)
  method \code{rrvglm.fit}
  uses iteratively reweighted least squares (IRLS).

  }
  \item{x.arg, y.arg}{
  logical values indicating whether the model matrix
  and response vector/matrix used in the fitting
  process should be assigned in the \code{x} and \code{y} slots.
  Note the model matrix is the LM model matrix; to get the VGLM
  model matrix type \code{model.matrix(vglmfit)} where
  \code{vglmfit} is a \code{vglm} object.

  }
  \item{constraints}{
  See \code{\link{vglm}}.

  }
  \item{extra, smart, qr.arg}{
  See \code{\link{vglm}}.

  }
  \item{\dots}{
  further arguments passed into \code{\link{rrvglm.control}}.

  }
}
\details{
  In this documentation, \eqn{M} is the
  number of linear predictors.
  For RR-VGLMs,
  the central formula is given by
  \deqn{\eta = B_1^T x_1 + A \nu}{%
         eta = B_1^T x_1 + A nu}
  where \eqn{x_1}{x1} is a vector
  (usually just a 1 for an intercept),
  \eqn{x_2}{x2} is another vector of explanatory variables, and
  \eqn{\nu = C^T x_2}{nu = C^T x_2} is an \eqn{R}-vector of
  latent variables.
  Here, \eqn{\eta}{eta} is a vector of linear predictors,
  e.g., the \eqn{m}th element is
  \eqn{\eta_m = \log(E[Y_m])}{eta_m = log(E[Y_m])} for the
  \eqn{m}th Poisson response.
  The dimension of \eqn{\eta} is \eqn{M} by
  definition.
  The matrices \eqn{B_1}, \eqn{A} and
  \eqn{C} are estimated from the data, i.e., contain the
  regression coefficients.  For ecologists, the central
  formula represents a \emph{constrained linear ordination}
  (CLO) since it is linear in the latent variables. It
  means that the response is a monotonically increasing or
  decreasing function of the latent variables.



  For identifiability it is common to enforce
  \emph{corner constraints} on \bold{A}:
  by default, for RR-VGLMs, the top \eqn{R}
  by \eqn{R} submatrix is fixed to be
  the order-\eqn{R} identity matrix and
  the remainder of \bold{A} is estimated.
  And by default, for DRR-VGLMs, there is
  also an order-\eqn{R} identity matrix
  embedded in \bold{A} because the RRR must
  be \emph{separable} (this is so that any
  existing structure in \bold{A} is preserved).


  The underlying algorithm of RR-VGLMs is
  iteratively reweighted least squares (IRLS)
  with an optimizing algorithm applied within
  each IRLS iteration (e.g., alternating
  algorithm).


  In theory, any \pkg{VGAM} family function
  that works for \code{\link{vglm}}
  and \code{\link{vgam}} should work for
  \code{rrvglm} too.  The function that
  actually does the work is \code{rrvglm.fit};
  it is essentially \code{vglm.fit} with some
  extra code.


}
\value{
  For RR-VGLMs,
  an object of class \code{"rrvglm"}, which
  has the the same slots as a \code{"vglm"}
  object. The only difference is that the some
  of the constraint matrices are estimates
  rather than known. But \pkg{VGAM} stores
  the models the same internally. The slots
  of \code{"vglm"} objects are described in
  \code{\link{vglm-class}}.


  For DRR-VGLMs,
  an object of class \code{"drrvglm"}.


}
\references{

Yee, T. W. and Hastie, T. J. (2003).
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


Yee, T. W. (2004).
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


Anderson, J. A. (1984).
Regression and ordered categorical variables.
\emph{Journal of the Royal Statistical Society,
      Series B, Methodological},
\bold{46}, 1--30.


Yee, T. W. (2014).
Reduced-rank vector generalized linear models with
two linear predictors.
\emph{Computational Statistics and Data Analysis},
\bold{71}, 889--902.



Yee, T. W., Frigau, L. and Ma, C. (2024).
Heaping and seeping,
GAITD regression and
doubly constrained reduced rank vector
generalized linear models,
in smoking studies.
\emph{In preparation}.



%  Documentation accompanying the \pkg{VGAM} package at
%  \url{http://www.stat.auckland.ac.nz/~yee}
%  contains further information and examples.

}

\author{ Thomas W. Yee }
\note{
  The arguments of \code{rrvglm} are in general the same
  as those of \code{\link{vglm}} but with some extras in
  \code{\link{rrvglm.control}}.


  The smart prediction (\code{\link{smartpred}}) library
  is packed with the \pkg{VGAM} library.


  In an example below, a rank-1 \emph{stereotype}
  (reduced-rank multinomial logit)
  model of Anderson (1984) is fitted to some car data.
  The reduced-rank regression is performed, adjusting for
  two covariates. Setting a trivial constraint matrix
  (\code{diag(M)})
  for the latent variable variables in \eqn{x_2}{x2} avoids
  a warning message when it is overwritten by a (common)
  estimated constraint matrix.  It shows that German cars
  tend to be more expensive than American cars, given a
  car of fixed weight and width.


  If \code{fit <- rrvglm(..., data = mydata)} then
  \code{summary(fit)} requires corner constraints and no
  missing values in \code{mydata}. 
  Sometimes the estimated variance-covariance
  matrix of the parameters is not
  positive-definite; if this occurs, try
  refitting the model with a different value
  for \code{Index.corner}.


  For \emph{constrained quadratic ordination} (CQO) see
  \code{\link{cqo}} for more details about QRR-VGLMs.


  With multiple binary responses, one must use
  \code{binomialff(multiple.responses = TRUE)} to indicate
  that the response is a matrix with one response per column.
  Otherwise, it is interpreted as a single binary response
  variable.


  To fit DRR-VGLMs see the arguments
  \code{H.A.thy} and
  \code{H.C} in \code{\link{rrvglm.control}}.
  DRR-VGLMs provide structure to the \bold{A} and
  \bold{C} matrices via constraint matrices.
  So instead of them being general unstructured
  matrices, one can make specified elements to
  be identically equal to 0, for example.
  This gives greater control
  over what is modelled as a latent variable,
  e.g., in a health study,
  if one subset of the covariates are physical
  variables and the remainder are psychological
  variables then a rank-2 model might have each
  latent variable a linear combination of each
  of the types of variables separately.


% Called "variant1" internally: 20231228:


  Incidentally
  before I forget, since \code{Corner = TRUE},
  then
  the differences between the \code{@H.A.thy} and
  \code{@H.A.alt} slots
  are due to \code{Index.corner},
  which specifies which rows of \bold{A}
  are not estimated.
  However,
  in the alternating algorithm,
  it is more efficient to estimate the entire
  \bold{A}, bar (effectively) rows \code{str0},
  and then normalize it.
  In contrast, optimizing over the subset of
  \bold{A} to be estimated is slow.


% 20240325:

In the \code{@misc} slot are logical components
\code{is.drrvglm} and
\code{is.rrvglm}.
Only one is \code{TRUE}.
If \code{is.rrvglm} then (full) corner constraints
are used.
If \code{is.drrvglm} then
\emph{restricted corner constraints} (RCCs)
are used and the reduced rank regression
(RRR) must be \emph{separable}.
The case \code{is.rrvglm} means that
\code{H.A.thy} is a \code{vector("list", Rank)}
with \code{H.A.thy[[r]] <- diag(M)} assigned
to all \eqn{r=1,\ldots,R}.
Because DRR-VGLMs are implemented only for
separable problems, this means that
all columns of \code{H.A.thy[[s]]}
are orthogonal to all columns from
\code{H.A.try[[t]]}, for all \eqn{s} and \eqn{t}.
DRR-VGLMs are proposed in
Yee et al. (2024) in the context of GAITD regression
for heaped and seeped survey data.




}

% zzz; arguments of \code{\link{vglm}} are definitive.
% Theyre  copied here.



\seealso{
    \code{\link{rrvglm.control}},
    \code{\link{summary.drrvglm}},
    \code{\link{lvplot.rrvglm}}
    (same as \code{\link{biplot.rrvglm}}),
    \code{\link{rrvglm-class}},
    \code{\link{grc}},
    \code{\link{cqo}},
    \code{\link{vglmff-class}},
    \code{\link{vglm}},
    \code{\link{vglm-class}},
    \code{\link{smartpred}},
    \code{rrvglm.fit}.
    Special family functions include
    \code{\link{negbinomial}}
    \code{\link{zipoisson}}
    and \code{\link{zinegbinomial}}.
    (see Yee (2014)
    and what was formerly in \pkg{COZIGAM}).
    Methods functions include
    \code{\link{Coef.rrvglm}},
    \code{\link{calibrate.rrvglm}},
    etc.
    Data include
    \code{\link{crashi}}.


%   \code{\link{qrrvglm.control}},
%   \code{\link{vcovqrrvglm}},




}

\examples{
\dontrun{
# Example 1: RR NB with Var(Y) = mu + delta1 * mu^delta2
nn <- 1000       # Number of observations
delta1 <- 3.0    # Specify this
delta2 <- 1.5    # Specify this; should be greater than 1
a21 <- 2 - delta2
mydata <- data.frame(x2 = runif(nn), x3 = runif(nn))
mydata <- transform(mydata, mu = exp(2 + 3 * x2 + 0 * x3))
mydata <- transform(mydata,
    y2 = rnbinom(nn, mu = mu, size = (1/delta1)*mu^a21))
plot(y2 ~ x2, mydata, pch = "+", col = 'blue', las = 1,
  main = paste0("Var(Y) = mu + ", delta1, " * mu^", delta2))
rrnb2 <- rrvglm(y2 ~ x2 + x3, negbinomial(zero = NULL),
                data = mydata, trace = TRUE)

a21.hat <- (Coef(rrnb2)@A)["loglink(size)", 1]
beta11.hat <- Coef(rrnb2)@B1["(Intercept)", "loglink(mu)"]
beta21.hat <- Coef(rrnb2)@B1["(Intercept)", "loglink(size)"]
(delta1.hat <- exp(a21.hat * beta11.hat - beta21.hat))
(delta2.hat <- 2 - a21.hat)
# delta1.hat:
# exp(a21.hat * predict(rrnb2)[1,1] - predict(rrnb2)[1,2])
summary(rrnb2)

# Obtain a 95 percent CI for delta2:
se.a21.hat <- sqrt(vcov(rrnb2)["I(latvar.mat)", "I(latvar.mat)"])
ci.a21 <- a21.hat +  c(-1, 1) * 1.96 * se.a21.hat
(ci.delta2 <- 2 - rev(ci.a21))  # The 95 percent CI

Confint.rrnb(rrnb2)  # Quick way to get it

# Plot the abundances and fitted values vs the latent variable
plot(y2 ~ latvar(rrnb2), data = mydata, col = "blue",
     xlab = "Latent variable", las = 1)
ooo <- order(latvar(rrnb2))
lines(fitted(rrnb2)[ooo] ~ latvar(rrnb2)[ooo], col = "red")

# Example 2: stereotype model (RR multinomial logit model)
data(car.all)
scar <- subset(car.all,
    is.element(Country, c("Germany", "USA", "Japan", "Korea")))
fcols <- c(13,14,18:20,22:26,29:31,33,34,36)  # These are factors
scar[, -fcols] <- scale(scar[, -fcols])  # Stdze all numerical vars
ones <- CM.ones(3)  # matrix(1, 3, 1)
clist <- list("(Intercept)" = diag(3), Width = ones, Weight = ones,
              Disp. = diag(3), Tank = diag(3), Price = diag(3),
              Frt.Leg.Room = diag(3))
set.seed(111)
fit <- rrvglm(Country ~ Width + Weight + Disp. + Tank +
              Price + Frt.Leg.Room,
              multinomial, data = scar, Rank = 2, trace = TRUE,
              constraints = clist, noRRR = ~ 1 + Width + Weight,
#             Uncor = TRUE, Corner = FALSE,  # orig.
              Index.corner = c(1, 3),  # Less correlation
              Bestof = 3)
fit@misc$deviance  # A history of the fits
Coef(fit)
biplot(fit, chull = TRUE, scores = TRUE, clty = 2, Ccex = 2,
       ccol = "blue", scol = "orange", Ccol = "darkgreen",
       Clwd = 2, main = "1=Germany, 2=Japan, 3=Korea, 4=USA")
}
}
\keyword{models}
\keyword{regression}
\concept{Reduced-Rank Vector Generalized Linear Model}
\concept{Reduced-rank regression}


%index <- with(car.all, Country == "Germany" | Country == "USA" |
%                       Country == "Japan"   | Country == "Korea")
%scar <- car.all[index, ]  # standardized car data


%scar <- subset(car.all,
%     is.element(Country, c("Germany", "USA", "Japan", "Korea")) |
%     is.na(Country))

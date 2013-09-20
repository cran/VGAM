\name{cqo}
\alias{cqo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting Constrained Quadratic Ordination (CQO)}
\description{
  A \emph{constrained quadratic ordination} (CQO; formerly called
  \emph{canonical Gaussian ordination} or CGO) model is fitted using
  the \emph{quadratic reduced-rank vector generalized linear model}
  (QRR-VGLM) framework.

}
\usage{
cqo(formula, family, data = list(), weights = NULL, subset = NULL,
    na.action = na.fail, etastart = NULL, mustart = NULL,
    coefstart = NULL, control = qrrvglm.control(...), offset = NULL,
    method = "cqo.fit", model = FALSE, x.arg = TRUE, y.arg = TRUE,
    contrasts = NULL, constraints = NULL, extra = NULL,
    smart = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{formula}{ a symbolic description of the model to be fit.
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
  Currently the following families are supported:
  \code{\link{poissonff}},
  \code{\link{binomialff}}
      (\code{\link{logit}} and \code{\link{cloglog}} links available),
  \code{\link{negbinomial}},
  \code{\link{gamma2}},
  \code{\link{gaussianff}}.
  Sometimes special arguments are required for \code{cqo()}, e.g.,
  \code{binomialff(mv = TRUE)}.
  Also, \code{\link{quasipoissonff}} and \code{\link{quasibinomialff}}
  may or may not work.


% \code{negbinomial(deviance = TRUE)},
% \code{gamma2(deviance = TRUE)}.


  }
  \item{data}{
   an optional data frame containing the variables in the model.
   By default the variables are taken from \code{environment(formula)},
   typically the environment from which \code{cqo} is called.


  }
  \item{weights}{ an optional vector or matrix of (prior) weights 
    to be used in the fitting process.
    Currently, this argument should not be used.


  }
  \item{subset}{
  an optional logical vector specifying a subset of
  observations to be used in the fitting process.


  }
  \item{na.action}{ 
  a function which indicates what should happen when the data contain
  \code{NA}s.  The default is set by the \code{na.action} setting of
  \code{\link[base]{options}}, and is \code{na.fail} if that is unset.
  The ``factory-fresh'' default is \code{na.omit}.


  }
  \item{etastart}{
  starting values for the linear predictors.
  It is a \eqn{M}-column matrix.
  If \eqn{M = 1} then it may be a vector.
  Currently, this argument probably should not be used.


  }
  \item{mustart}{
  starting values for the 
  fitted values. It can be a vector or a matrix. 
  Some family functions do not make use of this argument.
  Currently, this argument probably should not be used.


  }
  \item{coefstart}{
  starting values for the
  coefficient vector.
  Currently, this argument probably should not be used.


  }
  \item{control}{
  a list of parameters for controlling the fitting process.
  See \code{\link{qrrvglm.control}} for details.


  }
  \item{offset}{ 
  This argument must not be used.

%   especially when \code{ITolerances = TRUE}.
%   a vector or \eqn{M}-column matrix of offset values.
%   These are \emph{a priori} known and are
%   added to the linear predictors during fitting.

 }
  \item{method}{
  the method to be used in fitting the model.
  The default (and presently only) method \code{cqo.fit}
  uses \emph{iteratively reweighted least squares} (IRLS).


  }
  \item{model}{
  a logical value indicating whether the \emph{model frame}
  should be assigned in the \code{model} slot.


  }
  \item{x.arg, y.arg}{
  logical values indicating whether
  the model matrix and response matrix used in the fitting
  process should be assigned in the \code{x} and \code{y} slots.
  Note the model matrix is the LM model matrix.


%    ; to get the VGLM
%    model matrix type \code{model.matrix(vglmfit)} where
%    \code{vglmfit} is a \code{vglm} object. 


  }
  \item{contrasts}{
  an optional list. See the \code{contrasts.arg}
  of \code{model.matrix.default}.


  }
  \item{constraints}{
  an optional list  of constraint matrices.
  The components of the list must be named with the term it corresponds
  to (and it must match in character format). 
  Each constraint matrix must have \eqn{M} rows, and be of full-column
  rank. By default, constraint matrices are the \eqn{M} by \eqn{M}
  identity
  matrix unless arguments in the family function itself override these values. 
  If \code{constraints} is used it must contain \emph{all} the
  terms; an incomplete list is not accepted.
  Constraint matrices for \eqn{x_2}{x_2} variables are taken as the
  identity matrix.


  }
  \item{extra}{
  an optional list with any extra information that might be needed
  by the family function.


  }
% \item{qr.arg}{ logical value indicating whether
%   the slot \code{qr}, which returns the QR decomposition of the
%   VLM model matrix, is returned on the object.
%   }
  \item{smart}{
  logical value indicating whether smart prediction
  (\code{\link{smartpred}}) will be used.


  }
  \item{\dots}{ 
  further arguments passed into \code{\link{qrrvglm.control}}.


  }
}
\details{
  QRR-VGLMs or \emph{constrained quadratic ordination} (CQO) models
  are estimated here by maximum likelihood estimation. Optimal linear
  combinations of the environmental variables are computed, called
  \emph{latent variables} (these appear as \code{latvar} for \eqn{R=1}
  else \code{latvar1}, \code{latvar2}, etc. in the output).  Here, \eqn{R}
  is the \emph{rank} or the number of ordination axes.  Each species'
  response is then a regression of these latent variables using quadratic
  polynomials on a transformed scale (e.g., log for Poisson counts, logit
  for presence/absence responses).  The solution is obtained iteratively
  in order to maximize the log-likelihood function, or equivalently,
  minimize the deviance.


  The central formula (for Poisson and binomial species data) is
  given by
  \deqn{\eta = B_1^T x_1 + A \nu +
               \sum_{m=1}^M (\nu^T D_m \nu) e_m}{%
         eta = B_1^T x_1 + A nu +
         sum_{m=1}^M (nu^T D_m nu) e_m}
  where \eqn{x_1}{x_1} is a vector (usually just a 1 for an intercept),
  \eqn{x_2}{x_2} is a vector of environmental variables, \eqn{\nu=C^T
  x_2}{nu=C^T x_2} is a \eqn{R}-vector of latent variables, \eqn{e_m} is
  a vector of 0s but with a 1 in the \eqn{m}th position.
  The \eqn{\eta}{eta} are a vector of linear/additive predictors,
  e.g., the \eqn{m}th element is \eqn{\eta_m = \log(E[Y_m])}{eta_m =
  log(E[Y_m])} for the \eqn{m}th species.  The matrices \eqn{B_1},
  \eqn{A}, \eqn{C} and \eqn{D_m} are estimated from the data, i.e.,
  contain the regression coefficients.  The tolerance matrices
  satisfy \eqn{T_s = -\frac12 D_s^{-1}}{T_s = -(0.5 D_s^(-1)}.
  Many important CQO details are directly related to arguments
  in \code{\link{qrrvglm.control}}, e.g., the argument \code{noRRR}
  specifies which variables comprise \eqn{x_1}{x_1}.


  Theoretically, the four most popular \pkg{VGAM} family functions
  to be used with \code{cqo} correspond to the Poisson, binomial,
  normal, and negative binomial distributions. The latter is a
  2-parameter model. All of these are implemented, as well as the
  2-parameter gamma.  The Poisson is or should be catered for by
  \code{\link{quasipoissonff}} and \code{\link{poissonff}}, and the
  binomial by \code{\link{quasibinomialff}} and \code{\link{binomialff}}.
  Those beginning with \code{"quasi"} have dispersion parameters that
  are estimated for each species.


  %the negative binomial by \code{\link{negbinomial}}, and the normal by
  %\code{gaussianff}.

  %For overdispersed Poisson data, using \code{\link{quasipoissonff}} is
  %strongly recommended over \code{\link{negbinomial}}; the latter is
  %\emph{very} sensitive to departures from the model assumptions.


  For initial values, the function \code{.Init.Poisson.QO} should
  work reasonably well if the data is Poisson with species having equal
  tolerances.  It can be quite good on binary data too.  Otherwise the
  \code{Cinit} argument in \code{\link{qrrvglm.control}} can be used.


  %(and negative binomial)


  It is possible to relax the quadratic form to an additive model.  The
  result is a data-driven approach rather than a model-driven approach,
  so that CQO is extended to \emph{constrained additive ordination}
  (CAO) when \eqn{R=1}.  See \code{\link{cao}} for more details.


  In this documentation, \eqn{M} is the number of linear predictors,
  \eqn{S} is the number of responses (species). Then
  \eqn{M=S} for Poisson and binomial species data,
  and \eqn{M=2S} for negative binomial and gamma distributed species data.

}
\value{
  An object of class \code{"qrrvglm"}. 
  Note that the slot \code{misc} has a list component called
  \code{deviance.Bestof} which gives the history of deviances over all
  the iterations.

}
\references{

Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


ter Braak, C. J. F. and Prentice, I. C. (1988)
A theory of gradient analysis.
\emph{Advances in Ecological Research},
\bold{18}, 271--317.


%Yee, T. W. (2005)
%On constrained and unconstrained 
%quadratic ordination.
%\emph{Manuscript in preparation}.


Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.

}
\author{
Thomas W. Yee.
Thanks to Alvin Sou for converting a lot of the
original FORTRAN code into C.


} 

\note{

  By default, a rank-1 equal-tolerances QRR-VGLM model is fitted
  (see \code{\link{qrrvglm.control}} for the default control
  parameters).
  The latent variables are always transformed so that they
  are uncorrelated.
  By default, the argument \code{trace} is \code{TRUE} meaning a running
  log is printed out while the computations are taking place.  This is
  because the algorithm is computationally expensive, therefore users
  might think that their computers have frozen if \code{trace = FALSE}!


  The argument \code{Bestof} in \code{\link{qrrvglm.control}} controls
  the number of models fitted (each uses different starting values) to
  the data. This argument is important because convergence may be to a
  \emph{local} solution rather than the \emph{global} solution. Using
  more starting values increases the chances of finding the global
  solution.  Always plot an ordination diagram (use the generic function
  \code{\link{lvplot}}) and see if it looks sensible.  Local solutions
  arise because the optimization problem is highly nonlinear, and this is
  particularly true for CAO.


  %Convergence of QRR-VGLMs can be difficult, especially for binary
  %data. If this is so, then setting \code{ITolerances = TRUE} or
  %\code{EqualTolerances = TRUE} may help, especially when the number of sites,
  %\eqn{n}, is small.

  %If the negative binomial family function \code{\link{negbinomial}} is 
  %used for \code{cqo} then set \code{negbinomial(deviance = TRUE)}
  %is necessary. This means to minimize the deviance, which the fast
  %algorithm can handle. 


  Many of the arguments applicable to \code{cqo} are common to
  \code{\link{vglm}} and \code{\link{rrvglm.control}}.
  The most important arguments are
  \code{Rank},
  \code{noRRR},
  \code{Bestof}, 
  \code{ITolerances},
  \code{EqualTolerances},
  \code{isd.latvar}, and
  \code{MUXfactor}.


  When fitting a 2-parameter model such as the negative binomial
  or gamma, it pays to have \code{EqualTolerances = TRUE} and
  \code{ITolerances = FALSE}. This is because numerical problems can
  occur when fitting the model far away from the global solution when
  \code{ITolerances = TRUE}. Setting the two arguments as described will
  slow down the computation considerably, however it is numerically
  more stable.


  In Example 1 below, an unequal-tolerances rank-1 QRR-VGLM is fitted to the
  hunting spiders dataset, and
  Example 2 is the equal-tolerances version. The latter is less likely to
  have convergence problems compared to the unequal-tolerances model.
  In Example 3 below, an equal-tolerances rank-2 QRR-VGLM is fitted to the
  hunting spiders dataset.
  The numerical difficulties encountered in fitting the rank-2 model
  suggests a rank-1 model is probably preferable.
  In Example 4 below, constrained binary quadratic ordination (in old
  nomenclature, constrained Gaussian logit ordination) is fitted to some
  simulated data coming from a species packing model.
  With multivariate binary responses, one must use \code{mv = TRUE} to
  indicate that the response (matrix) is multivariate. Otherwise, it is
  interpreted as a single binary response variable.
  In Example 5 below, the deviance residuals are plotted for each species.
  This is useful as a diagnostic plot.
  This is done by (re)regressing each species separately against the latent
  variable.


  Sometime in the future, this function might handle input of the form
  \code{cqo(x, y)}, where \code{x} and \code{y} are matrices containing
  the environmental and species data respectively.

}
\section{Warning }{
  Local solutions are not uncommon when fitting CQO models.  To increase
  the chances of obtaining the global solution, increase the value
  of the argument \code{Bestof} in \code{\link{qrrvglm.control}}.
  For reproducibility of the results, it pays to set a different
  random number seed before calling \code{cqo} (the function
  \code{\link[base:Random]{set.seed}} does this).  The function \code{cqo}
  chooses initial values for \bold{C} using \code{.Init.Poisson.QO()}
  if \code{Use.Init.Poisson.QO = TRUE}, else random numbers.


  Unless \code{ITolerances = TRUE} or \code{EqualTolerances = FALSE},
  CQO is computationally expensive. It pays to keep the rank down to 1
  or 2.  If \code{EqualTolerances = TRUE} and \code{ITolerances = FALSE} then
  the cost grows quickly with the number of species and sites (in terms of
  memory requirements and time).  The data needs to conform quite closely
  to the statistical model, and the environmental range of the data should
  be wide in order for the quadratics to fit the data well (bell-shaped
  response surfaces).  If not, RR-VGLMs will be more appropriate because
  the response is linear on the transformed scale (e.g., log or logit)
  and the ordination is called \emph{constrained linear ordination} or CLO.


  Like many regression models, CQO is sensitive to outliers (in the
  environmental and species data), sparse data, high leverage points,
  multicollinearity etc.  For these reasons, it is necessary to examine
  the data carefully for these features and take corrective action
  (e.g., omitting certain species, sites, environmental variables from
  the analysis, transforming certain environmental variables, etc.).
  Any optimum lying outside the convex hull of the site scores should not
  be trusted.  Fitting a CAO is recommended first, then upon transformations
  etc., possibly a CQO can be fitted.


  For binary data, it is necessary to have `enough' data.  In general,
  the number of sites \eqn{n} ought to be much larger than the number of
  species \emph{S}, e.g., at least 100 sites for two species. Compared
  to count (Poisson) data, numerical problems occur more frequently
  with presence/absence (binary) data.  For example, if \code{Rank = 1}
  and if the response data for each species is a string of all absences,
  then all presences, then all absences (when enumerated along the latent
  variable) then infinite parameter estimates will occur.  In general,
  setting \code{ITolerances = TRUE} may help.


  This function was formerly called \code{cgo}. It has been renamed to
  reinforce a new nomenclature described in Yee (2006).

}

\seealso{
  \code{\link{qrrvglm.control}},
  \code{\link{Coef.qrrvglm}},
  \code{\link{predictqrrvglm}},
  \code{\link{rcqo}},
  \code{\link{cao}},
  \code{\link{uqo}},
  \code{\link{rrvglm}},
%   \code{\link{rrvglm.control}},
  \code{\link{poissonff}},
  \code{\link{binomialff}},
  \code{\link{negbinomial}},
  \code{\link{gamma2}},
  \code{\link{lvplot.qrrvglm}},
  \code{\link{perspqrrvglm}},
  \code{\link{trplot.qrrvglm}},
% \code{\link{vcovqrrvglm}},
  \code{\link{vglm}},
  \code{\link[base:Random]{set.seed}},
  \code{\link{hspider}}.


Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.


}
\examples{
\dontrun{
# Example 1; Fit an unequal tolerances model to the hunting spiders data
hspider[,1:6] <- scale(hspider[,1:6]) # Standardized environmental variables
set.seed(1234) # For reproducibility of the results
p1ut <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
                  Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
                  Trocterr, Zoraspin) ~
            WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
            fam = poissonff, data = hspider, Crow1positive = FALSE,
            EqualTol = FALSE)
sort(p1ut@misc$deviance.Bestof) # A history of all the iterations
if(deviance(p1ut) > 1177) warning("suboptimal fit obtained")

S <- ncol(depvar(p1ut)) # Number of species
clr <- (1:(S+1))[-7] # Omits yellow
lvplot(p1ut, y = TRUE, lcol = clr, pch = 1:S, pcol = clr, las = 1) # ordination diagram
legend("topright", leg = colnames(depvar(p1ut)), col = clr,
       pch = 1:S, merge = TRUE, bty = "n", lty = 1:S, lwd = 2)
(cp <- Coef(p1ut))

(a <- cp@latvar[cp@latvar.order])  # The ordered site scores along the gradient
# Names of the ordered sites along the gradient:
rownames(cp@latvar)[cp@latvar.order]
(aa <- (cp@Optimum)[,cp@Optimum.order]) # The ordered optima along the gradient
aa <- aa[!is.na(aa)] # Delete the species that is not unimodal
names(aa) # Names of the ordered optima along the gradient

trplot(p1ut, whichSpecies = 1:3, log = "xy", type = "b", lty = 1, lwd = 2,
       col = c("blue","red","green"), label = TRUE) -> ii # trajectory plot
legend(0.00005, 0.3, paste(ii$species[, 1], ii$species[, 2], sep = " and "),
       lwd = 2, lty = 1, col = c("blue", "red", "green"))
abline(a = 0, b = 1, lty = "dashed")

S <- ncol(depvar(p1ut)) # Number of species
clr <- (1:(S+1))[-7] # Omits yellow
persp(p1ut, col = clr, label = TRUE, las = 1) # perspective plot


# Example 2; Fit an equal tolerances model. Less numerically fraught.
set.seed(1234)
p1et <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
                  Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
                  Trocterr, Zoraspin) ~
            WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
            poissonff, data = hspider, Crow1positive = FALSE)
sort(p1et@misc$deviance.Bestof) # A history of all the iterations
if (deviance(p1et) > 1586) warning("suboptimal fit obtained")
S <- ncol(depvar(p1et)) # Number of species
clr <- (1:(S+1))[-7] # Omits yellow
persp(p1et, col = clr, label = TRUE, las = 1)


# Example 3: A rank-2 equal tolerances CQO model with Poisson data
# This example is numerically fraught... need IToler = TRUE too.
set.seed(555)
p2 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
                Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
                Trocterr, Zoraspin) ~
          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
          poissonff, data = hspider, Crow1positive = FALSE,
          IToler = TRUE, Rank = 2, Bestof = 3, isd.latvar = c(2.1, 0.9))
sort(p2@misc$deviance.Bestof) # A history of all the iterations
if(deviance(p2) > 1127) warning("suboptimal fit obtained")
lvplot(p2, ellips = FALSE, label = TRUE, xlim = c(-3,4),
       C = TRUE, Ccol = "brown", sites = TRUE, scol = "grey", 
       pcol = "blue", pch = "+", chull = TRUE, ccol = "grey")


# Example 4: species packing model with presence/absence data
set.seed(2345)
n <- 200; p <- 5; S <- 5
mydata <- rcqo(n, p, S, fam = "binomial", hiabundance = 4,
               EqualTol = TRUE, ESOpt = TRUE, EqualMax = TRUE)
myform <- attr(mydata, "formula")
set.seed(1234)
b1et <- cqo(myform, binomialff(mv = TRUE, link = "cloglog"), data = mydata)
sort(b1et@misc$deviance.Bestof) # A history of all the iterations
lvplot(b1et, y = TRUE, lcol = 1:S, pch = 1:S, pcol = 1:S, las = 1)
Coef(b1et)

# Compare the fitted model with the 'truth'
cbind(truth = attr(mydata, "concoefficients"), fitted = concoef(b1et))


# Example 5: Plot the deviance residuals for diagnostic purposes
set.seed(1234)
p1et <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
                  Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
                  Trocterr, Zoraspin) ~
            WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
            poissonff, data = hspider, EqualTol = TRUE, trace = FALSE)
sort(p1et@misc$deviance.Bestof) # A history of all the iterations
if(deviance(p1et) > 1586) warning("suboptimal fit obtained")
S <- ncol(depvar(p1et))
par(mfrow = c(3, 4))
for(ii in 1:S) {
  tempdata <- data.frame(latvar1 = c(latvar(p1et)), sppCounts = depvar(p1et)[, ii])
  tempdata <- transform(tempdata, myOffset = -0.5 * latvar1^2)

# For species ii, refit the model to get the deviance residuals
  fit1 <- vglm(sppCounts ~ offset(myOffset) + latvar1, fam = poissonff,
               data = tempdata, trace = FALSE)

# For checking: this should be 0
  print("max(abs(c(Coef(p1et)@B1[1,ii], Coef(p1et)@A[ii,1]) - coef(fit1)))")
  print( max(abs(c(Coef(p1et)@B1[1,ii], Coef(p1et)@A[ii,1]) - coef(fit1))) )

# Plot the deviance residuals
  devresid <- resid(fit1, type = "deviance")
  predvalues <- predict(fit1) + fit1@offset
  ooo <- with(tempdata, order(latvar1))
  with(tempdata, plot(latvar1, predvalues + devresid, col = "darkgreen",
                      xlab = "latvar1", ylab = "", main = colnames(depvar(p1et))[ii]))
  with(tempdata, lines(latvar1[ooo], predvalues[ooo], col = "blue"))
}
}
}
\keyword{models}
\keyword{regression}

%legend("topright", x=1, y=135, leg = colnames(depvar(p1ut)), col = clr,
%       pch = 1:S, merge = TRUE, bty = "n", lty = 1:S, lwd = 2)

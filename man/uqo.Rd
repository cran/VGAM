\name{uqo}
\alias{uqo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting Unconstrained Quadratic Ordination (UQO)}
\description{
  An \emph{unconstrained quadratic ordination} (UQO)
  (equivalently, noncanonical Gaussian ordination) model
  is fitted using the 
  \emph{quadratic unconstrained vector generalized linear model}
  (QU-VGLM) framework.
  In this documentation, \eqn{M} is the number of linear predictors
  or species.

}
\usage{
uqo(formula, family, data = list(), weights = NULL, subset = NULL,
    na.action = na.fail, etastart = NULL, mustart = NULL,
    coefstart = NULL, control = uqo.control(...), offset = NULL,
    method = "uqo.fit", model = FALSE, x.arg = TRUE, y.arg = TRUE,
    contrasts = NULL, constraints = NULL, extra = NULL,
    qr.arg = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ a symbolic description of the model to be fit.
    Since there is no \eqn{x_2} vector by definition, the RHS of
    the formula has all terms belonging to the \eqn{x_1} vector.

  }
  \item{family}{ a function of class \code{"vglmff"} describing
    what statistical model is to be fitted. Currently two families
    are supported: Poisson and binomial.
  }
  \item{data}{ an optional data frame containing the variables
    in the model. By default the variables are taken from
    \code{environment(formula)}, typically the environment from
    which \code{uqo} is called.
 }
  \item{weights}{ an optional vector or matrix of (prior) weights 
    to be used in the fitting process.
    This argument should not be used.

}
  \item{subset}{ an optional logical vector specifying a subset of
          observations to 
          be used in the fitting process.
    }
    \item{na.action}{
      a function which indicates what should happen when
      the data contain \code{NA}s. 
      The default is set by the \code{na.action} setting
      of \code{\link[base]{options}}, and is \code{na.fail}
      if that is unset.
      The ``factory-fresh'' default is \code{na.omit}.
    }
  \item{etastart}{ starting values for the linear predictors.
    It is a \eqn{M}-column matrix. If \eqn{M=1} then it may be a vector.
    }
  \item{mustart}{ starting values for the 
    fitted values. It can be a vector or a matrix. 
    Some family functions do not make use of this argument.
  }
  \item{coefstart}{ starting values for the
    coefficient vector. }
  \item{control}{ a list of parameters for controlling the fitting process. 
          See \code{\link{uqo.control}} for details.
 }
  \item{offset}{ a vector or \eqn{M}-column matrix of offset values.
   This argument should not be used.
 }
  \item{method}{
    the method to be used in fitting the model.
    The default (and presently only) method \code{uqo.fit}
    uses iteratively reweighted least squares (IRLS).
    }
  \item{model}{ a logical value indicating whether the
    \emph{model frame}
    should be assigned in the \code{model} slot. }
  \item{x.arg, y.arg}{ logical values indicating whether
    the model matrix and response matrix used in the fitting
    process should be assigned in the \code{x} and \code{y} slots.
    Note the model matrix is the LM model matrix.

    }
  \item{contrasts}{ an optional list. See the \code{contrasts.arg}
    of \code{\link{model.matrix.default}}. }
  \item{constraints}{ an optional list  of constraint matrices.
    This argument should not be used.
    }
  \item{extra}{ an optional list with any extra information that  
    might be needed by the family function. 
    }
  \item{qr.arg}{ logical value indicating whether
    the slot \code{qr}, which returns the QR decomposition of the
    VLM model matrix, is returned on the object.
    This argument should not be set \code{TRUE}.
    }
  \item{\dots}{ further arguments passed into \code{\link{uqo.control}}. }

}

\details{
  \emph{Unconstrained quadratic ordination} models fit symmetric bell-shaped
  response curves/surfaces to response data, but the latent variables
  are largely free parameters and are not constrained to be linear
  combinations of the environmental variables.  This poses a
  difficult optimization problem.  The current algorithm is very simple
  and will often fail (even for \code{Rank=1}) but hopefully this will
  be improved in the future.

  The central formula is given by
  \deqn{\eta = B_1^T x_1 + A \nu +
               \sum_{m=1}^M (\nu^T D_m \nu) e_m}{%
         eta = B_1^T x_1 + A nu +
         sum_{m=1}^M (nu^T D_m nu) e_m}
  where \eqn{x_1}{x_1} is a vector (usually just a 1 for an intercept),
  \eqn{\nu}{nu} is a \eqn{R}-vector of latent variables, \eqn{e_m} is
  a vector of 0s but with a 1 in the \eqn{m}th position.
  The \eqn{\eta}{eta} are a vector of linear/additive predictors,
  e.g., the \eqn{m}th element is \eqn{\eta_m = \log(E[Y_m])}{eta_m =
  log(E[Y_m])} for the \eqn{m}th species.  The matrices \eqn{B_1},
  \eqn{A}, and \eqn{D_m} are estimated from the data, i.e.,
  contain the regression coefficients. Also, \eqn{\nu}{nu} is
  estimated.
  The tolerance matrices satisfy \eqn{T_s = -\frac12 D_s^{-1}}{T_s =
  -(0.5 D_s^(-1)}.  Many important UQO details are directly related to
  arguments in \code{\link{uqo.control}};
  see also \code{\link{cqo}} and \code{\link{qrrvglm.control}}.

Currently, only Poisson and binomial \pkg{VGAM} family functions are
implemented for this function, and dispersion parameters for these are
assumed known.  Thus the Poisson is catered for by
\code{\link{poissonff}}, and the binomial by \code{\link{binomialff}}.
Those beginning with \code{"quasi"} have dispersion parameters that are
estimated for each species, hence will give an error message here.

}
\value{
  An object of class \code{"uqo"}
  (this may change to \code{"quvglm"} in the future).
}
\references{

Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.

Yee, T. W. (2005)
On constrained and unconstrained
quadratic ordination.
\emph{Manuscript in preparation}.

Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.

}
\author{Thomas W. Yee} 

\note{

  The site scores are centered.
  When \eqn{R>1}, they are uncorrelated and should be unique up
  to a rotation.

The argument \code{Bestof} in \code{\link{uqo.control}} controls
the number of models fitted (each uses different starting values) to
the data. This argument is important because convergence may be to a
\emph{local} solution rather than the \emph{global} solution. Using more
starting values increases the chances of finding the global solution.
Local solutions arise because the optimization problem is highly
nonlinear.

In the example below, a CQO model is fitted and used for providing
initial values for a UQO model.

}
\section{Warning }{

  Local solutions are not uncommon when fitting UQO models.  To increase
  the chances of obtaining the global solution, set
  \code{ITolerances=TRUE} or \code{EqualTolerances=TRUE} and increase
  the value of the argument \code{Bestof} in \code{\link{uqo.control}}.
  For reproducibility of the results, it pays to set a different random
  number seed before calling \code{uqo} (the function
  \code{\link[base:Random]{set.seed}} does this).

The function \code{uqo} is very sensitive to initial values, and there
is a lot of room for improvement here.

UQO is computationally expensive.  It pays to keep the rank to no more
than 2, and 1 is much preferred over 2.
The data needs to conform closely to the statistical model.

Currently there is a bug with the argument \code{Crow1positive}
in \code{\link{uqo.control}}. This argument might be interpreted
as controlling the sign of the first site score, but currently
this is not done.

}

\seealso{
  \code{\link{uqo.control}},
  \code{\link{cqo}},
  \code{\link{qrrvglm.control}},
  \code{\link{rcqo}},
% \code{\link{cao}},
\code{\link{poissonff}},
\code{\link{binomialff}},
  \code{Coef.uqo},
  \code{lvplot.uqo},
  \code{persp.uqo},
  \code{trplot.uqo},
  \code{vcov.uqo},
  \code{\link[base:Random]{set.seed}},
  \code{\link{hspider}}.
}
\examples{
\dontrun{
data(hspider)
set.seed(123)  # This leads to the global solution
hspider[,1:6] = scale(hspider[,1:6]) # Standardized environmental vars
p1 = cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
               Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
               Trocterr, Zoraspin) ~
         WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
         ITolerances = TRUE, fam = poissonff, data = hspider, 
         Crow1positive=TRUE, Bestof=3, trace=FALSE)
if(deviance(p1) > 1589.0) stop("suboptimal fit obtained")

set.seed(111)
up1 = uqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
                Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
                Trocterr, Zoraspin) ~ 1,
          family = poissonff, data = hspider,
          ITolerances = TRUE,
          Crow1positive = TRUE, lvstart = -lv(p1))
if(deviance(up1) > 1310.0) stop("suboptimal fit obtained")

nos = ncol(up1@y) # Number of species
clr = (1:(nos+1))[-7]  # to omit yellow
lvplot(up1, las=1, y=TRUE, pch=1:nos, scol=clr, lcol=clr, 
       pcol=clr, llty=1:nos, llwd=2)
legend(x=2, y=135, dimnames(up1@y)[[2]], col=clr, lty=1:nos,
       lwd=2, merge=FALSE, ncol=1, x.inter=4.0, bty="l", cex=0.9)

# Compare the site scores between the two models
plot(lv(p1), lv(up1), xlim=c(-3,4), ylim=c(-3,4), las=1)
abline(a=0, b=-1, lty=2, col="blue", xpd=FALSE)
cor(lv(p1, ITol=TRUE), lv(up1))

# Another comparison between the constrained and unconstrained models
# The signs are not right so they are similar when reflected about 0 
par(mfrow=c(2,1))
persp(up1, main="Red/Blue are the constrained/unconstrained models",
      label=TRUE, col="blue", las=1)
persp(p1, add=FALSE, col="red")
1-pchisq(deviance(p1) - deviance(up1), df=52-30)
}
}
\keyword{models}
\keyword{regression}

% 6/10/06; when the bug is fixed:
%persp(p1, add=TRUE, col="red")



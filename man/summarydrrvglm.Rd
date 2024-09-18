% Adapted from summaryvglm.Rd
% 20240103

\name{summary.drrvglm}
\alias{summary.drrvglm}
\alias{summary.rrvglm}
\alias{show.summary.drrvglm}
\alias{show.summary.rrvglm}
\title{Summarizing
  Reduced Rank
  Vector Generalized Linear Model (RR-VGLM) and
  Doubly constrained RR-VGLM Fits

}
\description{
  These functions are all \code{\link{methods}}
  for class \code{"drrvglm"} or
  \code{"summary.drrvglm"} objects, or
  for class \code{"rrvglm"} or
  \code{"summary.rrvglm"} objects.


}
\usage{
\method{summary}{drrvglm}(object, correlation = FALSE, dispersion = NULL,
    digits = NULL, numerical = TRUE, h.step = 0.005, omit123 = FALSE,
     omit13 = FALSE, fixA = FALSE, presid = FALSE,
    signif.stars = getOption("show.signif.stars"),
    nopredictors = FALSE, eval0 = TRUE, ...)
\method{show}{summary.drrvglm}(x, digits = NULL,
    quote = TRUE, prefix = "", signif.stars = NULL)
\method{summary}{rrvglm}(object, correlation = FALSE, dispersion = NULL,
    digits = NULL, numerical = TRUE, h.step = 0.005, omit123 = FALSE,
     omit13 = FALSE, fixA = TRUE, presid = FALSE,
    signif.stars = getOption("show.signif.stars"),
    nopredictors = FALSE, upgrade = FALSE, ...)
\method{show}{summary.rrvglm}(x, digits = NULL,
    quote = TRUE, prefix = "", signif.stars = NULL)
}
\arguments{
  \item{object}{an object of class
    \code{"drrvglm"} or \code{"rrvglm"},
    a result of a call to
    \code{\link{rrvglm}}.}
  \item{x}{an object of class
    \code{"summary.drrvglm"} or
    \code{"summary.rrvglm"},
    a result of a call to
    \code{\link{summary.drrvglm}} or
    \code{\link{summary.rrvglm}}.

  }
  \item{dispersion}{
    used mainly for GLMs.
    Not really implemented in \pkg{VGAM} so
    should not be used.
%    See \code{\link[stats]{summary.glm}}.

  }
  \item{correlation}{
    See \code{\link{summaryvglm}}.
  }
  \item{digits}{
    See \code{\link{summaryvglm}}.
  }
  \item{signif.stars}{
    See \code{\link{summaryvglm}}.
  }
  \item{presid, quote}{
    See \code{\link{summaryvglm}}.
  }
  \item{nopredictors}{
    See \code{\link{summaryvglm}}.
  }
  \item{upgrade}{
    Logical.
    Upgrade \code{object} to
    \code{\link{drrvglm-class}}?
    Treating the object as a DRR-VGLM has
    advantages since the framework is larger.
    The code for ordinary RR-VGLMs was written
    a long time ago so it is a good idea
    to check that both give the same answer.
  }





% Arguments more specific to RR-VGLMs:
\item{numerical}{Logical,
  use a finite difference approximation
  for partial derivatives?
  If \code{FALSE} then theoretical formulas
  are used (however this option may no longer
  be implemented).


  
  }
\item{h.step}{Numeric,
  positive and close to 0.
  If \code{numerical} then
  this is the forward step
  for each finite difference approximation.
  That is, it plays the role of
  \eqn{h} in \eqn{(f(x+h)-f(x))/h} for
  some function \eqn{f}.
  If the overall variance-covariance matrix
  is not positive-definite, varying
  this argument can make a difference,
  e.g., increasing it to
  \code{0.01} is recommended.


  }
  \item{fixA}{Logical,
    if \code{TRUE} then the largest block matrix
    is for \bold{B1} and \bold{C}, else
    it is for \bold{A} and \bold{B1}.
    This should not make any difference because
    both estimates of \bold{B1} should be
    extremely similar, including the SEs.


  }
  \item{omit13}{Logical,
    if \code{TRUE} then the (1,3) block matrix
    is set to \bold{O}. That is,
    \bold{A} and \bold{C} are assumed to
    asymptotically uncorrelated.
    Setting this \code{TRUE} is an option when
    \bold{V} (see below) is not
    positive-definite.
    If this fails,
    another option that is often better
    is to set \code{omit123 = TRUE}.


% 20240226; \code{kill.all} renamed to omit123.


  }
  \item{omit123}{Logical.
    If \code{TRUE} then \emph{two}
    block matrices are set to \bold{O}
  (blocks (1,2) and (1,3), else
  blocks (1,3) and (2,3),
  depending on \code{fixA}),
  This will almost surely result in an
  overall variance-covariance matrix
  that is positive-definite, however, the
  SEs will be biased.
  This argument is more extreme than
  \code{omit13}.


% set to 0 (cov12=cov13 = 0, or cov13=cov23 = 0).

  }
% \item{\dots}{further arguments passed to
%  or from other methods.}
% \item{HDEtest}{logical;
% }
% \item{hde.NA}{logical;
% }
  \item{prefix}{
    See \code{\link{summaryvglm}}.
  }

  \item{eval0}{
    Logical.
    Check if \bold{V} is positive-definite?
    That is, all its eigenvalues are positive.
  }

  \item{\dots}{
    Logical argument \code{check.2} might work here.
    If \code{TRUE} then some quantities
    are printed out, for checking and debugging.


%   Not used currently.
  }


}
\details{
  Most of this document concerns DRR-VGLMs but
  also apply equally well to RR-VGLMs as a special
  case.



  The overall variance-covariance matrix
  The overall variance-covariance matrix
  (called \bold{V} below)
  is computed. Since the parameters
  comprise the elements of
  the matrices \bold{A}, \bold{B1} and \bold{C}
  (called here block matrices 1, 2, 3
  respectively), 
  and an alternating algorithm is used for
  estimation, then there are two overlapping
  submodels that are fitted within an IRLS
  algorithm. These have blocks 1 and 2, and
  2 and 3, so that \bold{B1} is common to both.
  They are combined into one large overall
  variance-covariance matrix.
  Argument \code{fixA} specifies which submodel
  the \bold{B1} block is taken from.
  Block (1,3) is the most difficult to
  compute and numerical approximations based on
  first derivatives are used by default for this.



  Sometimes the computed \bold{V}
  is not positive-definite.
  If so,
  then the standard errors will be \code{NA}.
  To avoid this problem,
  try varying \code{h.step}
  or refitting the model with a different
  \code{Index.corner}.
  Argument \code{omit13} and
  \code{omit123}
  can also be used to
  give approximate answers.
  If \bold{V} is not positive-definite
  then this may indicate
  that the model does not fit the
  data very well, e.g.,
  \code{Rank} is not a good value.
  Potentially, there are many ways why
  the model may be ill-conditioned.
  Try several options and set \code{trace = TRUE}
  to monitor convergence---this is informative
  about how well the model and data agree.




  How can one fit an ordinary RR-VGLM as
  a DRR-VGLM?
  If one uses corner constraints (default) then
  one should input \code{H.A} as a list
  containing \code{Rank} \code{diag(M)}
  matrices---one for each column of \bold{A}.
  Then since \code{Corner = TRUE}
  by default, then
  \code{object@H.A.alt} has certain columns
  deleted due to corner constraints.
  In contrast, 
  \code{object@H.A.thy} is the
  \code{H.A} that was inputted.
  FYI, the
  \code{alt} suffix indicates the alternating
  algorithm, while
  the suffix \code{thy} stands for \emph{theory}.



}
\value{
  \code{summarydrrvglm} returns an object
  of class \code{"summary.drrvglm"}.

% see \code{\link{summary.drrvglm-class}}.


}
\references{


Chapter 5 of:
Yee, T. W. (2015).
Vector Generalized Linear and Additive Models:
With an Implementation in R.
New York, USA: \emph{Springer}.
Sections 5.2.2 and 5.3 are particularly relevant.


}



\author{ T. W. Yee.  }
\note{
  Note that \code{\link[stats:vcov]{vcov}}
  methods exist for \code{\link{rrvglm-class}}
  and \code{\link{drrvglm-class}} objects.


  Sometimes this function can take a long time
  and this is because numerical derivatives
  are computed.



}

\section{Warning }{
  DRR-VGLMs are a recent development so
  it will take some time to get things
  totally ironed out.
  RR-VGLMs were developed a long time ago and
  are more well-established, however they
  have only recently been documented here.


}



\seealso{
  \code{\link{rrvglm}},
  \code{\link{rrvglm.control}},
  \code{\link{vcovdrrvglm}},
  \code{\link{CM.free}},
  \code{\link{summaryvglm}},
  \code{\link{summary.rrvglm-class}},
  \code{\link{summary.drrvglm-class}}.



}
\examples{
\dontrun{  # Fit a rank-1 RR-VGLM as a DRR-VGLM.
set.seed(1); n <- 1000; S <- 6  # S must be even
myrank <- 1
rdata <- data.frame(x1 = runif(n), x2 = runif(n),
           x3 = runif(n), x4 = runif(n))
dval <- ncol(rdata)  # Number of covariates
# Involves x1, x2, ... a rank-1 model:
ymatrix <- with(rdata,
  matrix(rpois(n*S, exp(3 + x1 - 0.5*x2)), n, S))
H.C <- vector("list", dval)  # Ordinary "rrvglm"
for (i in 1:dval) H.C[[i]] <- CM.free(myrank)
names(H.C) <- paste0("x", 1:dval)
H.A <- list(CM.free(S))  # rank-1

rfit1 <- rrvglm(ymatrix ~ x1 + x2 + x3 + x4,
           poissonff, rdata, trace = TRUE)
class(rfit1)
dfit1 <- rrvglm(ymatrix ~ x1 + x2 + x3 + x4,
           poissonff, rdata, trace = TRUE,
           H.A = H.A,    # drrvglm
           H.C = H.C)    # drrvglm
class(dfit1)
Coef(rfit1)  # The RR-VGLM is the same as
Coef(dfit1)  # the DRR-VGLM.
max(abs(predict(rfit1) - predict(dfit1)))  # 0
abs(logLik(rfit1) - logLik(dfit1))  # 0
summary(rfit1)
summary(dfit1)
}
}
\keyword{models}
\keyword{regression}


% yettodo: add argument \code{zz}




%\method{summary}{vglm}(object, correlation = FALSE,
%           dispersion = NULL, digits = NULL,
%           presid = TRUE,
%           signif.stars = getOption("show.signif.stars"))





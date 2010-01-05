\name{vgam.control}
\alias{vgam.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Control function for vgam }
\description{
 Algorithmic constants and parameters for running \code{\link{vgam}}
 are set using this function. 
}
\usage{
vgam.control(all.knots = FALSE, bf.epsilon = 1e-07, bf.maxit = 30, 
             checkwz=TRUE, criterion = names(.min.criterion.VGAM),
             epsilon = 1e-07, maxit = 30, na.action = na.fail,
             nk = NULL, save.weight = FALSE, se.fit = TRUE,
             trace = FALSE, wzepsilon = .Machine$double.eps^0.75,
             ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
% zz na.action differs from vglm 


  \item{all.knots}{
  logical indicating if all distinct points of 
  the smoothing variables are to be used as knots.
  By default, \code{all.knots=TRUE} for
  \eqn{n \leq 40}{n <= 40}, and
  for \eqn{n > 40},
  the number of knots is approximately
  \eqn{40 + (n-40)^{0.25}}{40 + (n-40)^0.25}.
  This increases very slowly with \eqn{n}
  so that the number of knots is approximately between 50 and 60
  for large \eqn{n}.

  }
  \item{bf.epsilon}{
  tolerance used by the modified vector
  backfitting algorithm for testing convergence.
  Must be a positive number.

  }
  \item{bf.maxit}{
  maximum number of iterations allowed in
  the modified vector
  backfitting algorithm. Must be a positive integer.

  }
  \item{checkwz}{
  logical indicating whether the diagonal elements of
  the working weight matrices should be checked whether they are
  sufficiently positive, i.e., greater than \code{wzepsilon}. If not,
  any values less than \code{wzepsilon} are replaced with this value.

  }
  \item{criterion}{
  character variable describing what criterion is to
  be used to test for convergence. 
  The possibilities are listed in \code{.min.criterion.VGAM}, but
  most family functions only implement a few of these.

  }
  \item{epsilon}{
  positive convergence tolerance epsilon. Roughly
  speaking, the Newton-Raphson/Fisher-scoring/local-scoring iterations
  are assumed to have
  converged when two successive \code{criterion} values are within
  \code{epsilon} of each other.

  }
  \item{maxit}{
  maximum number of
  Newton-Raphson/Fisher-scoring/local-scoring iterations allowed.

  }
  \item{na.action}{
  how to handle missing values.
  Unlike the SPLUS \code{gam} function, \code{\link{vgam}} cannot handle
  \code{NA}s when smoothing.

  }
  \item{nk}{
  vector of length \eqn{d} containing positive integers.
  where \eqn{d} be the number of \code{\link{s}} terms
  in the formula. 
  Recycling is used if necessary. 
  The \eqn{i}th value is the number of B-spline coefficients to be
  estimated for each component function of the \eqn{i}th 
  \code{s()} term.
  \code{nk} differs from the number of knots by some constant.
  If specified, \code{nk} overrides the automatic knot selection procedure.

  }
  \item{save.weight}{
  logical indicating whether the \code{weights} slot
  of a \code{"vglm"} object will be saved on the object.
  If not, it will be reconstructed when needed, e.g., \code{summary}.

  }
  \item{se.fit}{
  logical indicating whether approximate
  pointwise standard errors are to be saved on the object.
  If \code{TRUE}, then these can be plotted with \code{plot(..., se=TRUE)}. 

  }
  \item{trace}{
  logical indicating if output should be produced for each iteration.

  }
  \item{wzepsilon}{
  Small positive number used to test whether the diagonals of the working
  weight matrices are sufficiently positive.

  }

% \item{xij}{
% formula giving terms making up a covariate-dependent term.
%
% }

  \item{\dots}{
  other parameters that may be picked up from control
  functions that are specific to the \pkg{VGAM} family function.
% zz see later.

  }
}
\details{

  Most of the control parameters are used within \code{vgam.fit} and
  you will have to look at that to understand the full details.  Many of
  the control parameters are used in a similar manner by \code{vglm.fit}
  (\code{\link{vglm}}) because the algorithm (IRLS) is very similar.

  Setting \code{save.weight=FALSE} is useful for some models because the
  \code{weights} slot of the object is often the largest and so less
  memory is used to store the object. However, for some \pkg{VGAM}
  family function, it is necessary to set \code{save.weight=TRUE} because
  the \code{weights} slot cannot be reconstructed later.

}
\value{
  A list with components matching the input names. A little error
  checking is done, but not much.
  The list is assigned to the \code{control} slot of \code{\link{vgam}} objects.

}
\references{ 
Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.

\url{http://www.stat.auckland.ac.nz/~yee}

}
\author{ Thomas W. Yee}

\note{
  \code{\link{vgam}} does not implement half-stepsizing, therefore parametric
  models should be fitted with \code{\link{vglm}}. Also, \code{\link{vgam}} is
  slower than  \code{\link{vglm}} too.

}

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{vgam}},
  \code{\link{vglm.control}},
  \code{\link{vsmooth.spline}},
  \code{\link{vglm}}. 

}

\examples{
data(pneumo)
pneumo = transform(pneumo, let=log(exposure.time))
vgam(cbind(normal, mild, severe) ~ s(let, df=3), multinomial,
     pneumo, trace=TRUE, eps=1e-4, maxit=10)
}
\keyword{models}
\keyword{regression}
\keyword{smooth}



% xij = NULL,




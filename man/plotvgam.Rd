\name{plotvgam}
\alias{plotvgam}
\alias{plot.vgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Default VGAM Plotting }
\description{
  Component functions of a \code{\link{vgam-class}} object can
  be plotted with \code{plotvgam()}. These are on the scale of
  the linear/additive predictor.

}
\usage{
plotvgam(x, newdata = NULL, y = NULL, residuals = NULL,
         rugplot = TRUE, se = FALSE, scale = 0, raw = TRUE,
         offset.arg = 0, deriv.arg = 0, overlay = FALSE,
         type.residuals = c("deviance", "working", "pearson", "response"),
         plot.arg = TRUE, which.term = NULL, which.cf = NULL,
         control = plotvgam.control(...), varxij = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{ A fitted \pkg{VGAM} object, e.g., produced by
  \code{\link{vgam}}, \code{\link{vglm}}, or \code{\link{rrvglm}}.


  }
  \item{newdata}{ Data frame.
  May be used to reconstruct the original data set.


  }
  \item{y}{ Unused.
  }
  \item{residuals}{
  Logical. If \code{TRUE} then residuals are plotted.
  See \code{type.residuals}


  }
\item{rugplot}{
  Logical. If \code{TRUE} then a rug plot is plotted at the
  foot of each plot. These values are jittered to expose ties.

  }
  \item{se}{
  Logical. If \code{TRUE} then approximate \eqn{\pm 2}{+-2} pointwise
  standard error bands are included in the plot.


  }
  \item{scale}{
  Numerical. By default, each plot will have its own
  y-axis scale. However, by specifying a value, each plot's y-axis
  scale will be at least \code{scale} wide.

    }
  \item{raw}{
  Logical. If \code{TRUE} then the smooth functions are those
  obtained directly by the algorithm, and are plotted without
  having to premultiply with the constraint matrices.
  If \code{FALSE} then the smooth functions have been premultiply by
  the constraint matrices.
  The \code{raw} argument is directly fed into \code{predict.vgam()}.


  }
  \item{offset.arg}{
  Numerical vector of length \eqn{r}.
  These are added to the component functions. Useful for
  separating out the functions when \code{overlay} is \code{TRUE}.
  If \code{overlay} is \code{TRUE} and there is one covariate then
  using the intercept values as the offsets can be a good idea.


  }
  \item{deriv.arg}{
  Numerical. The order of the derivative.
  Should be assigned an small
  integer such as 0, 1, 2. Only applying to \code{s()} terms,
  it plots the derivative.


  }
  \item{overlay}{
  Logical. If \code{TRUE} then component functions of the same
  covariate are overlaid on each other.
  The functions are centered, so \code{offset.arg} can be useful
  when \code{overlay} is \code{TRUE}.


  }
  \item{type.residuals}{
  if \code{residuals} is \code{TRUE} then the first
  possible value
  of this vector, is used to specify the type of residual.


  }

  \item{plot.arg}{
   Logical. If \code{FALSE} then no plot is produced.


  }
  \item{which.term}{
  Character or integer vector containing all terms to be
  plotted, e.g., \code{which.term = c("s(age)", "s(height"))} or
  \code{which.term = c(2, 5, 9)}.
  By default, all are plotted.


  }
  \item{which.cf}{ An integer-valued vector specifying which
  linear/additive predictors are to be plotted.
  The values must be from the set \{1,2,\ldots,\eqn{r}\}.
  By default, all are plotted.


  }
  \item{control}{
  Other control parameters. See \code{\link{plotvgam.control}}.


  }
  \item{\dots}{
  Other arguments that can be fed into
  \code{\link{plotvgam.control}}. This includes line colors,
  line widths, line types, etc.


  }

  \item{varxij}{ Positive integer.
  Used if \code{xij} of \code{\link{vglm.control}} was used,
  this chooses which inner argument the component is plotted against.
  This argument is related to \code{raw = TRUE} and terms such as
  \code{NS(dum1, dum2)} and constraint matrices that have more than
  one column. The default would plot the smooth against \code{dum1}
  but setting \code{varxij = 2} could mean plotting the smooth against
  \code{dum2}.
  See the \pkg{VGAM} website for further information.


  }

}
\details{
  In this help file \eqn{M} is the number of linear/additive
  predictors, and \eqn{r} is the number of columns of the
  constraint matrix of interest.


  Many of \code{plotvgam()}'s options can be found in
  \code{\link{plotvgam.control}}, e.g., line types, line widths,
  colors.


}
\value{
  The original object, but with the \code{preplot} slot of the object
  assigned information regarding the plot.



}
%\references{
%
%
%Yee, T. W. and Wild, C. J. (1996)
%Vector generalized additive models.
%\emph{Journal of the Royal Statistical Society, Series B, Methodological},
%\bold{58}, 481--493.
%
%
%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.
%
%
%}
\author{ Thomas W. Yee }

\note{
  While \code{plot(fit)} will work if \code{class(fit)}
  is \code{"vgam"}, it is necessary to use \code{plotvgam(fit)}
  explicitly otherwise.


  \code{plotvgam()} is quite buggy at the moment.


% \code{plotvgam()} works in a similar
% manner to S-PLUS's \code{plot.gam()}, however, there is no
% options for interactive construction of the plots yet.


}

\seealso{
  \code{\link{vgam}},
  \code{\link{plotvgam.control}},
  \code{predict.vgam},
  \code{\link{plotvglm}},
  \code{\link{vglm}}.


}
\examples{
coalminers <- transform(coalminers, Age = (age - 42) / 5)
fit <- vgam(cbind(nBnW, nBW, BnW, BW) ~ s(Age),
            binom2.or(zero = NULL), data = coalminers)
\dontrun{ par(mfrow = c(1,3))
plot(fit, se = TRUE, ylim = c(-3, 2), las = 1)
plot(fit, se = TRUE, which.cf = 1:2, lcol = "blue", scol = "orange",
     ylim = c(-3, 2))
plot(fit, se = TRUE, which.cf = 1:2, lcol = "blue", scol = "orange",
     overlay = TRUE) }
}
\keyword{models}
\keyword{regression}
\keyword{smooth}
\keyword{graphs}

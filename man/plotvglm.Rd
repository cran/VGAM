\name{plotvglm}
\alias{plotvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots for VGLMs }
\description{
  Currently not working, this function can be used to feed
  the object to the VGAM plotting function.
  In the future some diagnostic plots will be plotted.

}
\usage{
plotvglm(x, type = c("vglm", "vgam"),
         newdata = NULL, y = NULL, residuals = NULL,
         rugplot = TRUE, se = FALSE, scale = 0, raw = TRUE,
         offset.arg = 0, deriv.arg = 0, overlay = FALSE,
         type.residuals = c("deviance", "working", "pearson", "response"),
         plot.arg = TRUE, which.term = NULL, which.cf = NULL,
         control = plotvgam.control(...), varxij = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{
  Same as \code{\link{plotvgam}}.


  }
  \item{type}{
  Default is the first choice.
  Currently the first choice gives an error (not written yet).
  If \code{"vgam"} then all the arguments are fed into
  \code{\link{plotvgam}}.


  }

  \item{newdata, y, residuals, rugplot}{
    Same as \code{\link{plotvgam}}.


  }
  \item{se, scale, raw, offset.arg}{
    Same as \code{\link{plotvgam}}.

  }
  \item{deriv.arg, overlay, type.residuals}{
    Same as \code{\link{plotvgam}}.

  }
  \item{plot.arg, which.term, which.cf, control}{
    Same as \code{\link{plotvgam}}.

  }
  \item{\dots, varxij}{
    Same as \code{\link{plotvgam}}.

  }

}
\details{
  Currently this function has not been written.
  When this is done some diagnostic plots based on
  residuals and hatvalues will be done.
  In the meanwhile, this function can be used to
  call the plotting function for \code{\link{vgam}} objects.


}
\value{
  Same as \code{\link{plotvgam}}.


}
%\references{
%}
%\author{ Thomas W. Yee }

%\note{
% \code{plotvglm()} is quite buggy at the moment.

  
% \code{plotvglm()} works in a similar
% manner to S-PLUS's \code{plot.gam()}, however, there is no
% options for interactive construction of the plots yet. 

  
%}

\seealso{
  \code{\link{plotvgam}},
  \code{\link{plotvgam.control}},
  \code{\link{vglm}}.
  
  
}
\examples{
coalminers <- transform(coalminers, Age = (age - 42) / 5)
fit <- vglm(cbind(nBnW, nBW, BnW, BW) ~ sm.bs(Age),
            binom2.or(zero = NULL), data = coalminers)
\dontrun{ par(mfrow = c(1, 3))
plot(fit, type = "vgam", se = TRUE, ylim = c(-3, 2), las = 1)
plot(fit, type = "vgam", se = TRUE, which.cf = 1:2,
     lcol = "blue", scol = "orange", ylim = c(-3, 2))
plot(fit, type = "vgam", se = TRUE, which.cf = 1:2,
     lcol = "blue", scol = "orange", overlay = TRUE) }
}
\keyword{models}
\keyword{regression}
\keyword{smooth}
\keyword{graphs}

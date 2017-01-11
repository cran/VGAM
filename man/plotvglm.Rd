\name{plotvglm}
\alias{plotvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots for VGLMs }
\description{
  Currently this function plots the Pearson residuals versus
  the linear predictors (\eqn{M} plots) and
  plots the Pearson residuals versus
  the hat values (\eqn{M} plots).

}
\usage{
plotvglm(x, which = "(All)", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{
    An object of class \code{"vglm"} (see \code{\link{vglm-class}})
    or inherits from that class.


% Same as \code{\link{plotvgam}}.


  }
  \item{which}{
    If a subset of the plots is required, specify a subset of the
    numbers \code{1:(2*M)}.
    The default is to plot them all.


  }
  \item{\dots}{
    Arguments fed into the primitive \code{\link[graphics]{plot}}
    functions.

  }

}
\details{
  This function is under development.
  Currently it plots the Pearson residuals
  against the predicted
  values (on the transformed scale) and the hat values.
  There are \eqn{2M} plots in total, therefore
  users should call \code{\link[graphics]{par}}
  to assign, e.g., the \code{mfrow} argument.
  Note: Section 3.7 of Yee (2015) describes the
  Pearson residuals and hat values for VGLMs.


}
\value{
  Returns the object invisibly.


% Same as \code{\link{plotvgam}}.


}
%\references{
%}
\author{ T. W. Yee }

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
\dontrun{
ndata <- data.frame(x2 = runif(nn <- 200))
ndata <- transform(ndata, y1 = rnbinom(nn, mu = exp(3+x2), size = exp(1)))
fit1 <- vglm(y1 ~ x2, negbinomial, data = ndata, trace = TRUE)
coef(fit1, matrix = TRUE)
par(mfrow = c(2, 2))
plot(fit1)

# Manually produce the four plots
plot(fit1, which = 1, col = "blue", las = 1, main = "main1")
abline(h = 0, lty = "dashed", col = "gray50")
plot(fit1, which = 2, col = "blue", las = 1, main = "main2")
abline(h = 0, lty = "dashed", col = "gray50")
plot(fit1, which = 3, col = "blue", las = 1, main = "main3")
plot(fit1, which = 4, col = "blue", las = 1, main = "main4")
}
}
\keyword{models}
\keyword{regression}
\keyword{smooth}
\keyword{graphs}

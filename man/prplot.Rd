\name{prplot}
\alias{prplot}
\alias{prplot.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability Plots for Categorical Data Analysis

}
\description{
Plots the fitted probabilities for some very simplified special
cases of categorical data analysis models.

}
\usage{
prplot(object, control = prplot.control(...), ...)

prplot.control(xlab = NULL, ylab = "Probability", main = NULL, xlim = NULL,
    ylim = NULL, lty = par()$lty, col = par()$col, rcol = par()$col,
    lwd = par()$lwd, rlwd = par()$lwd, las = par()$las, rug.arg = FALSE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Currently only an \code{\link{cumulative}} object.
  This includes a \code{\link{propodds}} object since that
  \pkg{VGAM} family function is a special case of \code{\link{cumulative}}.

}
  \item{control}{
  List containing some basic graphical parameters.

}
  \item{xlab, ylab, main, xlim, ylim, lty }{
  See \code{\link[graphics]{par}} and \code{...} below.

}
  \item{col, rcol, lwd, rlwd, las, rug.arg}{
  See \code{\link[graphics]{par}} and \code{...} below.
  Arguments starting with \code{r} refer to the rug.
  Argument \code{rug.arg} is logical: add a rug for the distinct values of the
  explanatory variable?


}
  \item{\dots}{
  Arguments such as \code{xlab} which are fed into \code{prplot.control()}.
  Only a small selection of graphical arguments from
  \code{\link[graphics]{par}} are offered.


}
}
\details{
  For models involving one term in the RHS of the formula this function
  plots the fitted probabilities against the single explanatory variable.


}
\value{
  The object is returned invisibly with the \code{preplot} slot assigned.
  This is obtained by a call to \code{plotvgam()}.


}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%T. W. Yee
%}
\note{
  This function is rudimentary.


}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{cumulative}}.


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo)
M <- npred(fit)  # Or fit@misc$M
\dontrun{ prplot(fit)
prplot(fit, lty = 1:M, col = (1:M)+2, rug = TRUE, las = 1,
       ylim = c(0, 1), rlwd = 2) }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphs}
\keyword{models}
\keyword{regression}


\name{plotqrrvglm}
\alias{plotqrrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Model Diagnostic Plots for QRR-VGLMs }
\description{
  The residuals of a QRR-VGLM are plotted for model diagnostic purposes.
}
\usage{
plotqrrvglm(object, rtype = c("response", "pearson", "deviance", "working"),
            ask = FALSE,
            main = paste(Rtype, "residuals vs latent variable(s)"),
            xlab = "Latent Variable",
            I.tolerances = object@control$eq.tolerances, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{"qrrvglm"}. }
  \item{rtype}{ Character string giving residual type.
                By default, the first one is chosen. }
  \item{ask}{ Logical. If \code{TRUE}, the user is asked to hit the return
  key for the next plot. }
  \item{main}{ Character string giving the title of the plot. }
  \item{xlab}{ Character string giving the x-axis caption. }
  \item{I.tolerances}{ Logical. This argument is fed into
    \code{Coef(object, I.tolerances = I.tolerances)}.
  }
  \item{\dots}{ Other plotting arguments (see \code{\link[graphics]{par}}). }
}
\details{
  Plotting the residuals can be potentially very useful for checking
  that the model fit is adequate.


}
\value{
  The original object.


}
\references{


Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


}
\author{Thomas W. Yee}

\note{
  An ordination plot of a QRR-VGLM can be obtained
  by \code{\link{lvplot.qrrvglm}}.

}

\seealso{
  \code{\link{lvplot.qrrvglm}},
  \code{\link{cqo}}.


}

\examples{\dontrun{
# QRR-VGLM on the hunting spiders data
# This is computationally expensive
set.seed(111)  # This leads to the global solution
# hspider[, 1:6] <- scale(hspider[, 1:6])  # Standardize the environmental variables
p1 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
                Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
                Trocterr, Zoraspin) ~
          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
          quasipoissonff, data = hspider, Crow1positive = FALSE)
par(mfrow = c(3, 4))
plot(p1, rtype = "response", col = "blue", pch = 4, las = 1, main = "")
}
}
\keyword{dplot}
\keyword{models}
\keyword{regression}



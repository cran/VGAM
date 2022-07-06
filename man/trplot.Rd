\name{trplot}
\alias{trplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Trajectory Plot }
\description{
  Generic function for a trajectory plot.
}
\usage{
trplot(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object for which a trajectory plot
    is meaningful.
  }
  \item{\dots}{ Other arguments fed into the specific
    methods function of the model. They usually are graphical
    parameters, and sometimes they are fed
    into the methods function for \code{\link{Coef}}.


  }
}
\details{
  Trajectory plots can be defined in different ways for different
  models.
  Many models have no such notion or definition.


  For quadratic and additive ordination models they plot the
  fitted values of two species against each other (more than
  two is theoretically possible, but not implemented
  in this software
  yet).


}
\value{
  The value returned depends specifically on the methods
  function invoked.


}
\references{


Yee, T. W. (2020).
On constrained and unconstrained
quadratic ordination.
\emph{Manuscript in preparation}.


}
\author{ Thomas W. Yee }

%\note{
%}


\seealso{
  \code{\link{trplot.qrrvglm}},
  \code{\link{perspqrrvglm}},
  \code{\link{lvplot}}.


}

\examples{
\dontrun{ set.seed(123)
hspider[, 1:6] <- scale(hspider[, 1:6])  # Stdze environ. vars
p1cqo <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute,
                   Arctperi, Auloalbi, Pardlugu, Pardmont,
                   Pardnigr, Pardpull, Trocterr, Zoraspin) ~
            WaterCon + BareSand + FallTwig +
            CoveMoss + CoveHerb + ReflLux,
            poissonff, data = hspider, Crow1positive = FALSE)

nos <- ncol(depvar(p1cqo))
clr <- 1:nos  # OR (1:(nos+1))[-7]  to omit yellow

trplot(p1cqo, which.species = 1:3, log = "xy", lwd = 2,
       col = c("blue", "orange", "green"), label = TRUE) -> ii
legend(0.00005, 0.3, paste(ii$species[, 1], ii$species[, 2],
                           sep = " and "),
       lwd = 2, lty = 1, col = c("blue", "orange", "green"))
abline(a = 0, b = 1, lty = "dashed", col = "grey") }
}
\keyword{models}
\keyword{regression}


\name{calibrate.qrrvglm.control}
\alias{calibrate.qrrvglm.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Control Function for CQO/CAO Calibration }
\description{
  Algorithmic constants and parameters for running
  \code{\link{calibrate.qrrvglm}} are set using this function.

}
\usage{
calibrate.qrrvglm.control(object, trace = FALSE, Method.optim = "BFGS",
    gridSize = if (Rank == 1) 9 else 5, varI.latvar = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    The fitted CQO/CAO model. The user should ignore this argument.


%   The fitted CQO/UQO/CAO model. The user should ignore this argument.


  }
  \item{trace}{
    Logical indicating if output should be produced for each iteration. It
    is a good idea to set this argument to be \code{TRUE} since the
    computations are expensive.


  }
\item{Method.optim}{
  Character. Fed into the \code{method} argument of
  \code{\link[stats]{optim}}.


  }
\item{gridSize}{
  Numeric, recycled to length \code{Rank}.  Controls the resolution
  of the grid used for initial values.  For each latent variable,
  an equally spaced grid of length \code{gridSize} is cast from the
  smallest site score to the largest site score. Then the likelihood
  function is evaluated on the grid, and the best fit is chosen as the
  initial value. Thus increasing the value of \code{gridSize} increases
  the chance of obtaining the global solution, however, the computing
  time increases proportionately.


  }
  \item{varI.latvar}{
  Logical. For CQO objects only, this argument is fed into
  \code{\link{Coef.qrrvglm}}.


  }
  \item{\dots}{
  Avoids an error message for extraneous arguments.


  }
}
\details{
  Most CQO/CAO users will only need to make use of \code{trace}
  and \code{gridSize}. These arguments should be used inside their
  call to \code{\link{calibrate.qrrvglm}}, not this function
  directly.


  }
\value{
  A list which with the following components.
  \item{trace}{Numeric (even though the input can be logical). }
  \item{gridSize}{Positive integer. }
  \item{varI.latvar}{Logical.}


}
\references{

Yee, T. W. (2013)
On constrained and unconstrained quadratic ordination.
\emph{Manuscript in preparation}.


}
\author{T. W. Yee}
\note{
  Despite the name of this function, CAO models are handled
  as well.


% Despite the name of this function, UQO and CAO models are handled

}



\seealso{
  \code{\link{calibrate.qrrvglm}},
  \code{\link{Coef.qrrvglm}}.


}
\examples{
\dontrun{ hspider[, 1:6] <- scale(hspider[, 1:6])  # Needed when I.tol = TRUE
set.seed(123)
p1 <- cqo(cbind(Alopacce, Alopcune, Pardlugu, Pardnigr,
                Pardpull, Trocterr, Zoraspin) ~
          WaterCon + BareSand + FallTwig +
          CoveMoss + CoveHerb + ReflLux,
          family = poissonff, data = hspider, I.tol = TRUE)
sort(deviance(p1, history = TRUE))  # A history of all the iterations

siteNos <- 3:4  # Calibrate these sites
cp1 <- calibrate(p1, new = data.frame(depvar(p1)[siteNos, ]), trace = TRUE)
}

\dontrun{
# Graphically compare the actual site scores with their calibrated values
persp(p1, main = "Site scores: solid=actual, dashed=calibrated",
      label = TRUE, col = "blue", las = 1)
abline(v = latvar(p1)[siteNos], lty = 1,
       col = 1:length(siteNos))  # Actual site scores
abline(v = cp1, lty = 2, col = 1:length(siteNos))  # Calibrated values
}
}
\keyword{models}
\keyword{regression}


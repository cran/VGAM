\name{Max}
\alias{Max}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maximums }
\description{
  Generic function for the \emph{maximums} (maxima) of a model.

}
\usage{
Max(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object for which the computation or
    extraction of
    a maximum (or maximums) is meaningful.

  }
  \item{\dots}{ Other arguments fed into the specific
    methods function of the model. Sometimes they are fed
    into the methods function for \code{\link{Coef}}.

  }
}
\details{
  Different models can define a maximum in different ways.
  Many models have no such notion or definition.


  Maximums occur in quadratic and additive ordination,
  e.g., CQO or CAO.
  For these models the maximum is the fitted value at the
  optimum. For quadratic ordination models there is a formula
  for the optimum but for additive ordination models the
  optimum must be searched for numerically. If it occurs
  on the boundary, then the optimum is undefined. For
  a valid optimum, the fitted value at the optimum is the maximum.



% e.g., CQO or UQO or CAO.


}
\value{
  The value returned depends specifically on the methods
  function invoked.


}
\references{

Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.


}
\author{ Thomas W. Yee }

%\note{
%}
%\section{Warning }{
%}

\seealso{
  \code{Max.qrrvglm},
  \code{\link{Tol}},
  \code{\link{Opt}}.


}

\examples{
\dontrun{
set.seed(111)  # This leads to the global solution
hspider[,1:6] <- scale(hspider[,1:6])  # Standardized environmental vars
p1 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
                Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
                Trocterr, Zoraspin) ~
          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
          quasipoissonff, Bestof = 2, data = hspider, Crow1positive = FALSE)

Max(p1)

index <- 1:ncol(depvar(p1))
persp(p1, col = index, las = 1, llwd = 2)
abline(h = Max(p1), lty = 2, col = index)

}
}
\keyword{models}
\keyword{regression}


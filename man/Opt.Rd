\name{Opt}
\alias{Opt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Optimums }
\description{
  Generic function for the \emph{optimums} (or optima) of a model.


}
\usage{
Opt(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object for which the computation or
    extraction of an optimum (or optimums) is meaningful.


  }
  \item{\dots}{ Other arguments fed into the specific
    methods function of the model. Sometimes they are fed
    into the methods function for \code{\link{Coef}}.


  }
}
\details{
  Different models can define an optimum in different ways.
  Many models have no such notion or definition.


  Optimums occur in quadratic and additive ordination,
  e.g., CQO or CAO.
  For these models the optimum is the value of the latent
  variable where the maximum occurs, i.e., where the fitted value
  achieves its highest value.
  For quadratic ordination models there is a formula
  for the optimum but for additive ordination models the
  optimum must be searched for numerically. If it occurs
  on the boundary, then the optimum is undefined.
  At an optimum, the fitted value of the response is
  called the \emph{maximum}.


}
\value{
  The value returned depends specifically on the methods
  function invoked.


}
\references{

Yee, T. W. (2004).
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


Yee, T. W. (2006).
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.


}
\author{ Thomas W. Yee }

\note{
In ordination, the optimum of a species is sometimes
called the \emph{species score}.


}
%\section{Warning }{
%}

\seealso{
  \code{Opt.qrrvglm},
  \code{\link{Max}},
  \code{\link{Tol}}.


}

\examples{
\dontrun{
set.seed(111)  # This leads to the global solution
hspider[,1:6] <- scale(hspider[,1:6])  # Standardized environmental vars
p1 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
                Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
                Trocterr, Zoraspin) ~
          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
          family = poissonff, data = hspider, Crow1positive = FALSE)
Opt(p1)

clr <- (1:(ncol(depvar(p1))+1))[-7]  # Omits yellow
persp(p1, col = clr, las = 1, main = "Vertical lines at the optimums")
abline(v = Opt(p1), lty = 2, col = clr)
}
}
\keyword{models}
\keyword{regression}




% index <- 1:ncol(depvar(p1))
% persp(p1, col = index, las = 1, main = "Vertical lines at the optimums")
% # abline(v = Opt(p1), lty = 2, col = index)
% rug(Opt(p1), col = clr, side = 3)





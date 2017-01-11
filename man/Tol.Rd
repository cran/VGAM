\name{Tol}
\alias{Tol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Tolerances }
\description{
  Generic function for the \emph{tolerances} of a model.
}
\usage{
Tol(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object for which the computation or
    extraction of a tolerance or tolerances is meaningful.


  }
  \item{\dots}{ Other arguments fed into the specific
    methods function of the model. Sometimes they are fed
    into the methods function for \code{\link{Coef}}.


  }
}
\details{
  Different models can define an optimum in different ways.
  Many models have no such notion or definition.


  Tolerances occur in quadratic ordination, i.e., CQO and UQO.
  They have ecological meaning because a high tolerance
  for a species means the species can survive over a large
  environmental range (stenoecous species), whereas a
  small tolerance means the species' niche is small
  (eurycous species).
  Mathematically, the tolerance is like the variance of
  a normal distribution.


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


\note{
  Tolerances are undefined for `linear' and additive
  ordination models.
  They are well-defined for quadratic ordination models.


}
\section{Warning }{
  There is a direct inverse relationship between the scaling of
  the latent variables (site scores) and the tolerances.
  One normalization is for the latent variables to have unit
  variance.
  Another normalization is for all the tolerances to be unit.
  These two normalization cannot simultaneously hold in general.
  For rank-\emph{R>1} models it becomes more complicated because
  the latent variables are also uncorrelated. An important
  argument when fitting quadratic ordination models is whether
  \code{eq.tolerances} is \code{TRUE} or \code{FALSE}.
  See Yee (2004) for details.



}

\seealso{
  \code{Tol.qrrvglm}.
  \code{\link{Max}},
  \code{\link{Opt}},
  \code{\link{cqo}},
  \code{\link{rcim}} for UQO.


}

\examples{
\dontrun{
set.seed(111)  # This leads to the global solution
hspider[,1:6] <- scale(hspider[, 1:6])  # Standardized environmental vars
p1 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
               Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
               Trocterr, Zoraspin) ~
         WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
         Bestof = 2, quasipoissonff, data = hspider, Crow1positive = FALSE)

Tol(p1)
}
}
\keyword{models}
\keyword{regression}


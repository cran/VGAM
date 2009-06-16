\name{ccoef}
\alias{ccoef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract Model Constrained/Canonical Coefficients }
\description{
  \code{ccoef} is a generic function which extracts the constrained
  (canonical) coefficients from objects returned by certain modelling
  functions.

}
\usage{
ccoef(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object for which the extraction of canonical
    coefficients is meaningful.
  }
  \item{\dots}{ Other arguments fed into the specific
    methods function of the model.
  }
}
\details{
  For constrained quadratic and ordination models, \emph{canonical
  coefficients} are the elements of the \bold{C} matrix used to form
  the latent variables.  They are highly interpretable in ecology,
  and are looked at as weights or loadings.

  They are also applicable for reduced-rank VGLMs.

}
\value{
  The value returned depends specifically on the methods function invoked.

}
\references{
Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.

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

\section{Warning }{
  For QO models, there is a direct inverse relationship between the
  scaling of the latent variables (site scores) and the tolerances.
  One normalization is for the latent variables to have unit variance.
  Another normalization is for all the species' tolerances to be
  unit (provided \code{EqualTolerances} is \code{TRUE}).  These two
  normalizations cannot simultaneously hold in general.  For rank
  \eqn{R} models with \eqn{R>1} it becomes more complicated because
  the latent variables are also uncorrelated. An important argument when
  fitting quadratic ordination models is whether \code{EqualTolerances}
  is \code{TRUE} or \code{FALSE}.  See Yee (2004) for details.

}

\seealso{
   \code{\link{ccoef-method}},
   \code{ccoef.qrrvglm},
   \code{ccoef.cao},
   \code{\link[stats]{coef}}.

}
\examples{
set.seed(111)  # This leads to the global solution
hspider[,1:6] = scale(hspider[,1:6]) # Standardized environmental vars
p1 = cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
               Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
               Trocterr, Zoraspin) ~
         WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
         fam = quasipoissonff, data = hspider, Crow1positive=FALSE)
ccoef(p1)
}
\keyword{models}
\keyword{regression}


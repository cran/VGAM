\name{studentt}
\alias{studentt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Student t Distribution }
\description{
  Estimation of the degrees of freedom for a Student t distribution.
}
\usage{
studentt(link.df = "loglog", earg=list(), idf=NULL, nsimEIM=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.df}{
  Parameter link function for the degrees of freedom \eqn{\nu}{nu}.
  See \code{\link{Links}} for more choices.
  The default ensures the parameter is greater than unity.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{idf}{
  Optional initial value.
  If given, its value must be greater than 1.
  The default is to compute an initial value internally.

  }

  \item{nsimEIM}{
  See \code{\link{CommonVGAMffArguments}}.

  }
}
\details{
  The density function is
  \deqn{f(y) = \frac{\Gamma((\nu+1)/2)}{\sqrt{\nu \pi} \Gamma(\nu/2)}
               \left(1 + \frac{y^2}{\nu} \right)^{-(\nu+1)/2}}{%
        f(y) = (gamma((nu+1)/2) / (sqrt(nu*pi) gamma(nu/2))) *
               (1 + y^2 / nu)^{-(nu+1)/2}  }
  for all real \eqn{y}.
  Then \eqn{E(Y)=0} if \eqn{\nu>1}{nu>1} (returned as the fitted values),
  and \eqn{Var(Y)= \nu/(\nu-2)}{Var(Y)= nu/(nu-2)}
  for \eqn{\nu > 2}{nu > 2}.
  When \eqn{\nu=1}{nu=1} then the Student \eqn{t}-distribution 
  corresponds to the standard Cauchy distribution.
  The degrees of freedom is treated as a parameter to be estimated,
  and as real and not integer.

  Simulation is used to estimate the EIM.
  Consequently the results will be reproducible only if
  a function such as \code{\link[base:Random]{set.seed}} is used.
  Increasing the value of \code{nsimEIM} will give more accurate results.
  In general convergence will be slow, especially when there are
  covariates.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{
Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.

Student (1908)
The probable error of a mean.
\emph{Biometrika}, \bold{6}, 1--25.

}

\author{ T. W. Yee }
\note{
  Practical experience has shown reasonably good initial values are
  required. If convergence failure occurs try using \code{idf}.
  Local solutions are also possible.

  A standard normal distribution corresponds to a \emph{t} distribution
  with infinite degrees of freedom.  Consequently, if the data is close
  to normal, there may be convergence problems.

}
\seealso{
    \code{\link{normal1}},
    \code{\link{loglog}},
    \code{\link[stats]{TDist}}.

}
\examples{
n = 500
y = rt(n, df=exp(exp(0.5)))
fit = vglm(y ~ 1, studentt, trace=TRUE)
coef(fit, matrix=TRUE)
Coef(fit) 
}
\keyword{models}
\keyword{regression}

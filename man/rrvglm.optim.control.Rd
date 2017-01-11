\name{rrvglm.optim.control}
\alias{rrvglm.optim.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Control Function for rrvglm() Calling optim() }
\description{
  Algorithmic constants and parameters for running \code{optim} within
  \code{rrvglm} are set using this function.


}
\usage{
rrvglm.optim.control(Fnscale = 1, Maxit = 100,
                     Switch.optimizer = 3, Abstol = -Inf,
                     Reltol = sqrt(.Machine$double.eps), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fnscale}{ Passed into \code{optim} as \code{fnscale}. }
  \item{Maxit}{ Passed into \code{optim} as \code{maxit}. }
  \item{Switch.optimizer}{ Iteration number when the "Nelder-Mead" method
  of \code{optim} is switched to the quasi-Newton "BFGS" method.
  Assigning \code{Switch.optimizer} a negative number
  means always BFGS, while assigning \code{Switch.optimizer} a value
  greater than \code{maxits} means always use Nelder-Mead.
  }
  \item{Abstol}{ Passed into \code{optim} as \code{abstol}. }
  \item{Reltol}{ Passed into \code{optim} as \code{reltol}. }
  \item{\dots}{ Ignored. }

}
\details{
  See \code{\link[stats]{optim}} for more details.


}
\value{
  A list with components equal to the arguments.


}
%\references{ ~put references to the literature/web site here ~ }
\author{ Thomas W. Yee }
\note{
  The transition between optimization methods may be
  unstable, so users may have to vary the value of
  \code{Switch.optimizer}.


  Practical experience with \code{Switch.optimizer} shows
  that setting it to too large a value may lead to a local
  solution, whereas setting it to a low value will obtain
  the global solution.  It appears that, if BFGS kicks in too
  late when the Nelder-Mead algorithm is starting to converge
  to a local solution, then switching to BFGS will not be
  sufficient to bypass convergence to that local solution.


}

\seealso{
  \code{\link{rrvglm.control}},
  \code{\link[stats]{optim}}.


}
%\examples{
%}
\keyword{models}
\keyword{regression}




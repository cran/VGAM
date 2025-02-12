\name{niters}
\alias{niters}
\alias{niters.vlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Number of Iterations

}
\description{
  Extracts the number of IRLS iterations
  performed for a VGLM object.



}
\usage{
niters(object, ...)
niters.vlm(object, history = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
% 
\arguments{
  \item{object}{
    A \code{\link{vglm}} object.
    Currently a \code{\link{vgam}} object
    is accepted but the correct value
    is not returned.


  }
  \item{history}{
    Logical, if \code{TRUE} it returns the
    convergence history with respect to the
    criterion, e.g.,
    \code{vglm.control()[["criterion"]]}.


  }
  \item{\ldots}{
    Currently unused.


  }
}
\details{The number of
  iteratively reweighted least squares (IRLS)
  iterations needed for convergence
  (or non-convergence) does say something
  about the model.
  Since Fisher scoring has a linear convergence
  rate in general, it should take no more than
  10 iterations, say, for successful
  convergence. Much more indicates potential
  problems, e.g.,
  a large disagreement between data and the
  specified model.


}
\value{
  A non-negative integer by default.
  If \code{history = TRUE} then a matrix.


}
% \references{
% }

%\section{Warning }{
%    See \code{\link{gaitdpoisson}}

%}



%\author{ T. W. Yee}
\note{
  Step-halving may or may not affect the answer.

}

\seealso{
  \code{\link{vglm}},
  \code{\link{vglm.control}}.



% \code{\link[base]{setdiff}},

}

\examples{
fit <- vglm(rpois(9, 2) ~ 1, poissonff, crit = "c")
niters(fit)
niters(fit, history = TRUE)
}
\keyword{models}
\keyword{regression}
\keyword{utilities}




\name{loglink}
%\name{loge}
\alias{loglink}
%\alias{loge}
\alias{negloglink}
%\alias{negloge}
\alias{logneglink}
%\alias{logneg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Log Link Function, and Variants }
\description{
  Computes the log transformation, including its inverse and the
  first two derivatives.


}
\usage{
loglink(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
        short = TRUE, tag = FALSE)
negloglink(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE)
logneglink(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.


  }
  \item{bvalue}{
  See \code{\link{Links}}.


  }

  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The log link function is very commonly used for parameters that
  are positive.
  Here, all logarithms are natural logarithms, i.e., to base
  \eqn{e}.  Numerical values of \code{theta} close to 0 or out of
  range result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.


  The function \code{loglink} computes
  \eqn{\log(\theta)}{log(theta)} whereas \code{negloglink} computes
  \eqn{-\log(\theta)=\log(1/\theta)}{-log(theta)=log(1/theta)}.


  The function \code{logneglink} computes
  \eqn{\log(-\theta)}{log(-theta)}, hence is suitable for parameters
  that are negative, e.g.,
  a trap-shy effect in \code{\link{posbernoulli.b}}.


}
\value{
  The following concerns \code{loglink}.
  For \code{deriv = 0}, the log of \code{theta}, i.e.,
  \code{log(theta)} when \code{inverse = FALSE}, and if
  \code{inverse = TRUE} then \code{exp(theta)}.
  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a function of
  \code{theta} if \code{inverse = FALSE}, else if
  \code{inverse = TRUE} then it returns the reciprocal.


}
\references{
  McCullagh, P. and Nelder, J. A. (1989).
  \emph{Generalized Linear Models},
  2nd ed. London: Chapman & Hall.


}
\author{ Thomas W. Yee }

\note{
  This function was called \code{loge} to avoid conflict with the
  \code{\link[base:Log]{log}} function.
  Numerical instability may occur when \code{theta} is close to
  0 unless \code{bvalue} is used.


}

\seealso{
    \code{\link{Links}},
    \code{\link{explink}},
    \code{\link{logitlink}},
    \code{\link{logclink}},
    \code{\link{logloglink}},
    \code{\link[base:Log]{log}},
    \code{\link{logofflink}},
    \code{\link{lambertW}},
    \code{\link{posbernoulli.b}}.


}
\examples{
\dontrun{ loglink(seq(-0.2, 0.5, by = 0.1))
 loglink(seq(-0.2, 0.5, by = 0.1), bvalue = .Machine$double.xmin)
negloglink(seq(-0.2, 0.5, by = 0.1))
negloglink(seq(-0.2, 0.5, by = 0.1), bvalue = .Machine$double.xmin) }
logneglink(seq(-0.5, -0.2, by = 0.1))
}
\keyword{math}
\keyword{models}
\keyword{regression}

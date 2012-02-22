\name{loge}
\alias{loge}
\alias{nloge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Log link function }
\description{
  Computes the log transformation, including its inverse and the first
  two derivatives.

}
\usage{
loge(theta, earg = list(), inverse = FALSE, deriv = 0,
     short = TRUE, tag = FALSE)
nloge(theta, earg = list(), inverse = FALSE, deriv = 0,
      short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.


  }
  \item{earg}{
  Optional list. Extra argument for passing in additional information.
  Values of \code{theta} which are less than or equal to 0 can be
  replaced by the \code{bvalue} component of the list \code{earg}
  before computing the link function value.
  The component name \code{bvalue} stands for ``boundary value''.
  See \code{\link{Links}} for general information about \code{earg}.


  }
  \item{inverse}{ Logical. If \code{TRUE} the inverse function is computed. }
  \item{deriv}{ Order of the derivative. Integer with value 0, 1 or 2. }
  \item{short}{
  Used for labelling the \code{blurb} slot of a
  \code{\link{vglmff-class}} object.

  }
  \item{tag}{
  Used for labelling the linear/additive predictor in the
  \code{initialize} slot of a \code{\link{vglmff-class}} object.
  Contains a little more information if \code{TRUE}.


  }
}
\details{
  The log link function is very commonly used for parameters that
  are positive.
  Numerical values of \code{theta} close to 0 or out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.
  The arguments \code{short} and \code{tag} are used only if
  \code{theta} is character.


  The function \code{loge} computes
  \eqn{\log(\theta)}{log(theta)} whereas \code{nloge} computes
  \eqn{-\log(\theta)=\log(1/\theta)}{-log(theta)=log(1/theta)}.


}
\value{
  The following concerns \code{loge}.
  For \code{deriv = 0}, the log of \code{theta}, i.e., \code{log(theta)}
  when \code{inverse = FALSE}, and if \code{inverse = TRUE} then
  \code{exp(theta)}.
  For \code{deriv = 1}, then the function returns
  \emph{d} \code{theta} / \emph{d} \code{eta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


  Here, all logarithms are natural logarithms, i.e., to base \eqn{e}.


}
\references{
    McCullagh, P. and Nelder, J. A. (1989)
    \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}
\author{ Thomas W. Yee }

\note{
  This function is called \code{loge} to avoid conflict with the
  \code{\link[base:Log]{log}} function.


  Numerical instability may occur when \code{theta} is close to 0 unless
  \code{earg} is used.


}

\seealso{ 
    \code{\link{Links}},
    \code{\link{explink}},
    \code{\link{logit}},
    \code{\link{logc}},
    \code{\link{loglog}},
    \code{\link[base:Log]{log}},
    \code{\link{logoff}},
    \code{\link{lambertW}}.


}
\examples{
\dontrun{ loge(seq(-0.2, 0.5, by = 0.1))
 loge(seq(-0.2, 0.5, by = 0.1), earg = list(bvalue = .Machine$double.xmin))
nloge(seq(-0.2, 0.5, by = 0.1))
nloge(seq(-0.2, 0.5, by = 0.1), earg = list(bvalue = .Machine$double.xmin)) }
}
\keyword{math}
\keyword{models}
\keyword{regression}

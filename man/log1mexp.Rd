\name{log1mexp}
\alias{log1mexp}
\alias{log1pexp}

\title{
  Logarithms with an Unit Offset and Exponential Term
}
\description{
Computes \code{log(1 + exp(x))} and \code{log(1 - exp(-x))} accurately.

}
\usage{
log1mexp(x)
log1pexp(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A vector of reals (numeric). Complex numbers not allowed since
  \code{\link[base]{expm1}} and \code{\link[base]{log1p}} do not handle these.


}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Computes \code{log(1 + exp(x))} and \code{log(1 - exp(-x))}
accurately. An adjustment is made when \eqn{x} is away from 0
in value.


}
\value{

\code{log1mexp(x)} gives the value of \eqn{\log(1-\exp(-x))}{log(1-exp(-x))}.



\code{log1pexp(x)} gives the value of \eqn{\log(1+\exp(x))}{log(1+exp(x))}.



}
\references{

Maechler, Martin (2012).
Accurately Computing log(1-exp(-|a|)).
Assessed from the \pkg{Rmpfr} package.


}
\author{
This is a direct translation of the function in Martin Maechler's
(2012) paper by Xiangjie Xue
and T. W. Yee.


}
\note{
If \code{NA} or \code{NaN} is present in the input, the
corresponding output will be \code{NA}.


}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[base]{log1p}},
  \code{\link[base]{expm1}},
  \code{\link[base]{exp}},
  \code{\link[base]{log}}



}
\examples{
x <-  c(10, 50, 100, 200, 400, 500, 800, 1000, 1e4, 1e5, 1e20, Inf, NA)
log1pexp(x)
log(1 + exp(x))  # Naive; suffers from overflow
log1mexp(x)
log(1 - exp(-x))
y <- -x
log1pexp(y)
log(1 + exp(y))  # Naive; suffers from inaccuracy
}





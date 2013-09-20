\name{log1pexp}
\alias{log1pexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logarithms with an Unit Offset and Exponential Term


}
\description{
Computes \code{log(1 + exp(x))} accurately.

}
\usage{
log1pexp(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of reals (numeric).
Complex numbers not allowed since \code{\link{log1p}} does
not handle these.


}
}
\details{
  Computes \code{log(1 + exp(x))} accurately.
  An adjustment is made when \code{x} is positive and large in value.


}
\value{
  Returns \code{log(1 + exp(x))}.


}
%\references{
%
%}
%\author{
%T. W. Yee
%
%}
%\note{
%
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[base:log]{log1p}},
  \code{\link[base:log]{exp}}.


}
\examples{
x <-  c(10, 50, 100, 200, 400, 500, 800, 1000, 1e4, 1e5, 1e20, Inf)
log1pexp(x)
log(1 + exp(x))  # Naive; suffers from overflow
x <- -c(10, 50, 100, 200, 400, 500, 800, 1000, 1e4, 1e5, 1e20, Inf)
log1pexp(x)
log(1 + exp(x))  # Naive; suffers from inaccuracy
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}

\name{kendall.tau}
\alias{kendall.tau}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kendall's Tau Statistic

}
\description{
  Computes Kendall's Tau, which is a rank-based correlation measure,
  between two vectors.


}
\usage{
kendall.tau(x, y, exact = FALSE, max.n = 3000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
    Numeric vectors. Must be of equal length.
    Ideally their values are continuous and not too discrete.
    Let \code{length(x)} be \eqn{N}, say.


}
  \item{exact}{
    Logical. If \code{TRUE} then the exact value is computed.


}
\item{max.n}{
  Numeric. If \code{exact = FALSE} and \code{length(x)}
  is more than \code{max.n} then a random sample
  of \code{max.n} pairs are chosen.


}
}
\details{
  Kendall's tau is a measure of dependency in a bivariate distribution.
  Loosely, two random variables are \emph{concordant} if large values
  of one random variable are associated with large values of the
  other random variable.
  Similarly, two random variables are \emph{disconcordant} if large values
  of one random variable are associated with small values of the
  other random variable.
  More formally, if \code{(x[i] - x[j])*(y[i] - y[j]) > 0} then
  that comparison is concordant \eqn{(i \neq j)}.
  And if \code{(x[i] - x[j])*(y[i] - y[j]) < 0} then
  that comparison is disconcordant \eqn{(i \neq j)}.
  Out of \code{choose(N, 2}) comparisons,
  let \eqn{c} and \eqn{d} be the
  number of concordant and disconcordant pairs.
  Then Kendall's tau can be estimated by \eqn{(c-d)/(c+d)}.
  If there are ties then half the ties are deemed concordant and
  half disconcordant so that \eqn{(c-d)/(c+d+t)} is used.


}
\value{
  Kendall's tau, which lies between \eqn{-1} and \eqn{1}.



}
%\references{


%}
%\author{
%  T. W. Yee.


%}
%\note{
%This function has not been tested thoroughly.


%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\section{Warning}{
  If \code{length(x)} is large then
  the cost is \eqn{O(N^2)}, which is expensive!
  Under these circumstances
  it is not advisable to set \code{exact = TRUE} or \code{max.n} to a very
  large number.



}
\seealso{
  \code{\link{binormalcop}},
  \code{\link[stats]{cor}}.


}
\examples{
N <- 5000; x <- 1:N; y <- runif(N)
true.rho <- -0.8
ymat <- rbinorm(N, cov12 =  true.rho)  # Bivariate normal, aka N_2
x <- ymat[, 1]
y <- ymat[, 2]

\dontrun{plot(x, y, col = "blue")}

kendall.tau(x, y)  # A random sample is taken here
kendall.tau(x, y)  # A random sample is taken here

kendall.tau(x, y, exact = TRUE)  # Costly if length(x) is large
kendall.tau(x, y, max.n = N)     # Same as exact = TRUE

(rhohat <- sin(kendall.tau(x, y) * pi / 2))  # This formula holds for N_2 actually
true.rho  # rhohat should be near this value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}

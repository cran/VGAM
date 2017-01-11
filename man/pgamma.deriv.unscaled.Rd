\name{pgamma.deriv.unscaled}
\alias{pgamma.deriv.unscaled}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Derivatives of the Incomplete Gamma Integral (Unscaled Version)
}
\description{
The first two derivatives of the incomplete gamma integral
with scaling.

}
\usage{
pgamma.deriv.unscaled(q, shape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{q, shape}{
    As in \code{\link[stats]{pgamma}}
    and \code{\link{pgamma.deriv}} but
    these must be vectors of positive values only and finite.


}
}
\details{
  Define
  \deqn{G(x, a) = \int_0^x t^{a-1} e^{-t} dt}{G(a,x) =
    integral_0^x t^(a-1) exp(-t) dt}
  so that
  \eqn{G(x, a)} is \code{pgamma(x, a) * gamma(a)}.
  Write \eqn{x = q} and \code{shape =} \eqn{a}.
  The 0th and first and second derivatives with respect to \eqn{a}
  of \eqn{G} are returned. This function is similar in spirit to
  \code{\link{pgamma.deriv}}
  but here there is no gamma function to scale things.
  Currently a 3-column matrix is returned (in the future this
  may change and an argument may be supplied so that only what
  is required by the user is computed.)
  This function is based on Wingo (1989).


}
\value{
  The 3 columns, running from left to right, are the \code{0:2}th derivatives
  with respect to \eqn{a}.


}
\references{

  See \code{\link{truncweibull}}.


}
\author{
  T. W. Yee.


}


%\note{
% If convergence does not occur then try increasing the value of
% \code{tmax}.
%
%}

\section{Warning }{
  These function seems inaccurate for \code{q = 1} and \code{q = 2};
  see the plot below.


}




\seealso{
  \code{\link{pgamma.deriv}},
  \code{\link[stats]{pgamma}}.


}
\examples{
x <- 3; aa <- seq(0.3, 04, by = 0.01)
ans.u <- pgamma.deriv.unscaled(x, aa)
head(ans.u)

\dontrun{ par(mfrow = c(1, 3))
for (jay in 1:3) {
  plot(aa, ans.u[, jay], type = "l", col = "blue", cex.lab = 1.5,
       cex.axis = 1.5, las = 1, main = colnames(ans.u)[jay],
       log = "", xlab = "shape", ylab = "")
  abline(h = 0, v = 1:2, lty = "dashed", col = "gray")  # Inaccurate at 1 and 2
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}



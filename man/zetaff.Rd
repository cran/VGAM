\name{zetaff}
\alias{zetaff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zeta Distribution Family Function }
\description{
  Estimates the parameter of the zeta distribution.

}
\usage{
zetaff(lshape = "loge", ishape = NULL, gshape = exp(-3:4)/4, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, ishape, zero}{
  These arguments apply to the (positive) parameter \eqn{p}.
  See \code{\link{Links}} for more choices.
  Choosing \code{\link{loglog}} constrains \eqn{p>1}, but
  may fail if the maximum likelihood estimate is less than one.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{gshape}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
In this long tailed distribution
the response must be a positive integer.
The probability function for a response \eqn{Y} is
\deqn{P(Y=y) = 1/[y^{p+1} \zeta(p+1)],\ \ \  p>0,\ \ \  y=1,2,...}{%
      P(Y=y) = 1/(y^(p+1) zeta(p+1)), p>0, y=1,2,...}
where \eqn{\zeta}{zeta} is Riemann's zeta function.
The parameter \eqn{p} is positive, therefore a log link
is the default.
The mean of \eqn{Y} is
\eqn{\mu = \zeta(p) / \zeta(p+1)}{mu = zeta(p)/zeta(p+1)}
(provided \eqn{p>1}) and these are the fitted values.
The variance of \eqn{Y} is
\eqn{\zeta(p-1) / \zeta(p+1) - \mu^2}{zeta(p-1) / zeta(p+1) - mu^2}
provided \eqn{p>2}.



It appears that good initial values are needed for successful
convergence. If convergence is not obtained, try several values
ranging from values near 0 to values about 10 or more.



Multiple responses are handled.



}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.



}

%Lindsey, J. K. (1995)
%\emph{Modelling Frequency and Count Data}.
%Oxford: Clarendon Press.

\references{


pp.527-- of Chapter 11 of
    Johnson N. L., Kemp, A. W. and Kotz S. (2005)
    \emph{Univariate Discrete Distributions},
    3rd edition,
    Hoboken, New Jersey: Wiley.


Knight, K. (2000)
\emph{Mathematical Statistics}.
Boca Raton: Chapman & Hall/CRC Press.



}
\author{ T. W. Yee }
\note{
  The \code{\link{zeta}} function may be used to compute values
  of the zeta function.


}

\seealso{
  \code{\link{zeta}},
  \code{\link{oazeta}},
  \code{\link{oizeta}},
  \code{\link{otzeta}},
  \code{\link{diffzeta}},
  \code{\link{dzeta}},
  \code{\link{hzeta}},
  \code{\link{zipf}}.


}
\examples{
zdata <- data.frame(y = 1:5, w =  c(63, 14, 5, 1, 2))  # Knight, p.304
fit <- vglm(y ~ 1, zetaff, data = zdata, trace = TRUE, weight = w, crit = "c")
(phat <- Coef(fit))  # 1.682557
with(zdata, cbind(round(dzeta(y, phat) * sum(w), 1), w))

with(zdata, weighted.mean(y, w))
fitted(fit, matrix = FALSE)
predict(fit)

# The following should be zero at the MLE:
with(zdata, mean(log(rep(y, w))) + zeta(1+phat, deriv = 1) / zeta(1+phat))
}
\keyword{models}
\keyword{regression}
% Also known as the Joos model or discrete Pareto distribution.

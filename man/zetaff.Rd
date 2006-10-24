\name{zetaff}
\alias{zetaff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zeta Distribution Family Function }
\description{
  Estimates the parameter of the zeta distribution.
}
\usage{
zetaff(link = "loge", init.p = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Parameter link function applied to the (positive) parameter \eqn{p}.
  See \code{\link{Links}} for more choices.
  Choosing \code{\link{loglog}} constrains \eqn{p>1}, but
  may fail if the maximum likelihood estimate is less than one.

  }
  \item{init.p}{
  Optional initial value for the parameter \eqn{p}.
  The default is to choose an initial value internally.
  If converge failure occurs use this argument to input a value.

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
provided \eqn{p>1}.
The variance of \eqn{Y} is
\eqn{\zeta(p-1) / \zeta(p+1) - \mu^2}{zeta(p-1) / zeta(p+1) - mu^2}
provided \eqn{p>2}.

It appears that good initial values are needed for successful
convergence. If convergence is not obtained, try several values
ranging from values near 0 to values about 10 or more.
}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.
}

%Lindsey, J. K. (1995)
%\emph{Modelling Frequency and Count Data}.
%Oxford: Clarendon Press. 

\references{ 
pp.465--471, Chapter 11 of
Johnson NL, Kotz S, and Kemp AW (1993)
\emph{Univariate Discrete Distributions},
2nd ed.
New York: Wiley.

Knight, K. (2000)
\emph{Mathematical Statistics}. 
Boca Raton: Chapman & Hall/CRC Press.
}
\author{ T. W. Yee }
\note{ The \code{\link{zeta}} function may be used to
compute values of the zeta function.
}

\seealso{ 
\code{\link{zeta}},
\code{\link{dzeta}},
\code{\link{hzeta}},
\code{\link{zipf}}.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.

}
\examples{
y = 1:5     # Knight, p.304
w =  c(63, 14, 5, 1, 2)
fit = vglm(y ~ 1, zetaff, trace=TRUE, wei=w, crit="c")
(phat = Coef(fit)) # 1.682557
cbind(dzeta(y, phat) * sum(w), w)

weighted.mean(y, w)
fitted(fit, mat=FALSE)
predict(fit)

# MLE should satisfy the following:
mean(log(rep(y, w))) + zeta(1+phat, deriv=1)/zeta(1+phat) # Should be 0
}
\keyword{models}
\keyword{regression}

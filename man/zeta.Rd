\name{zeta}
\alias{zeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Riemann's Zeta Function }
\description{
  Computes Riemann's zeta function and its first two derivatives.

}
\usage{
zeta(x, deriv = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ 
   A complex-valued vector/matrix whose real values must be \eqn{\geq
   1}{>= 1}. Otherwise, if \code{x} may be real.  If \code{deriv} is 1
   or 2 then \code{x} must be real and positive.

  }
  \item{deriv}{ 
  An integer equalling 0 or 1 or 2, which is the order of the derivative.
  The default means it is computed ordinarily.

  }
}
\details{
  While the usual definition involves an infinite series, more efficient
  methods have been devised to compute the value.  In particular,
  this function uses Euler-Maclaurin summation.  Theoretically, the
  zeta function can be computed over the whole complex plane because of
  analytic continuation.

  The formula used here for analytic continuation is
  \deqn{\zeta(s) = 2^s \pi^{s-1} \sin(\pi s/2) \Gamma(1-s) \zeta(1-s).}{%
        zeta(s)  = 2^s * pi^(s-1) * sin(pi*s/2) * gamma(1-s) *  zeta(1-s).}
  This is actually one of several formulas, but this one was discovered
  by Riemann himself and is called the \emph{functional equation}.

}

\section{Warning}{
  This function has not been fully tested, especially the derivatives.
  In particular, analytic continuation does not work here for
  complex \code{x} with \code{Re(x)<1} because currently the
  \code{\link[base:Special]{gamma}} function does not handle complex
  arguments.

}

\value{
  A vector/matrix of computed values.

%  The derivative is attached as an attribute zz.
}
\references{ 

Riemann, B. (1859)
Ueber die Anzahl der Primzahlen unter einer gegebenen Grosse.
\emph{Monatsberichte der Berliner Akademie, November 1859}.

Edwards, H. M. (1974)
\emph{Riemann's Zeta Function}.
Academic Press: New York. 

Markman, B. (1965)
The Riemann zeta function.
\emph{BIT},
\bold{5},
138--141.

Abramowitz, M. and Stegun, I. A. (1972)
\emph{Handbook of Mathematical Functions with Formulas, Graphs, and
Mathematical Tables},
New York: Dover Publications Inc.

}
\author{ T. W. Yee, with the help of G. J. Tee. }
\note{ 
  Estimation of the parameter of the zeta distribution can be achieved
  with \code{\link{zetaff}}.

}

\seealso{ 
  \code{\link{zetaff}},
  \code{\link{lerch}},
  \code{\link[base:Special]{gamma}}.
}
\examples{
zeta(2:10)

\dontrun{
x = seq(1.2, 10, len=100)
plot(x, zeta(x), type="l", las=1, xlim=c(-12,10), ylim=c(-1,4), col="red")
x = seq(-12, 0.8, len=100)
lines(x, zeta(x), col="red")
abline(v=0, h=c(0,1), lty="dashed")

# Close up plot
x = seq(-14, -0.4, len=100)
plot(x, zeta(x), type="l", las=1, col="red")
abline(v=0, h=0, lty="dashed")

# Plot of the first derivatives
x = seq(0.04, 0.8, len=100)
plot(x, zeta(x, deriv=1), type="l", las=1, col="blue", xlim=c(0.04,3), ylim=c(-6,0))
x = seq(1.2, 3, len=100)
lines(x, zeta(x, deriv=1), col="blue")
abline(v=0, h=0, lty="dashed") }

zeta(2) - pi^2 / 6      # Should be zero
zeta(4) - pi^4 / 90     # Should be zero
zeta(6) - pi^6 / 945    # Should be 0
zeta(8) - pi^8 / 9450   # Should be 0
# zeta(0, deriv=1) + 0.5 * log(2*pi) # Should be 0
}
\keyword{math}


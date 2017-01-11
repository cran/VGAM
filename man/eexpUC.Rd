\name{Expectiles-Exponential}
\alias{Expectiles-Exponential}
\alias{eexp}
\alias{deexp}
\alias{peexp}
\alias{qeexp}
\alias{reexp}
\title{ Expectiles of the Exponential Distribution }
\description{
  Density function, distribution function, and
  expectile function and random generation for the distribution
  associated with the expectiles of an exponential distribution.


}
\usage{
deexp(x, rate = 1, log = FALSE)
peexp(q, rate = 1, lower.tail = TRUE, log.p = FALSE)
qeexp(p, rate = 1, Maxit.nr = 10, Tol.nr = 1.0e-6,
      lower.tail = TRUE, log.p = FALSE)
reexp(n, rate = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, p, q}{
  See \code{\link{deunif}}.


  }
  \item{n, rate, log}{
  See \code{\link[stats:Exponential]{rexp}}.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Exponential]{pexp}}
  or \code{\link[stats:Exponential]{qexp}}.


  }
  \item{Maxit.nr, Tol.nr}{
  See \code{\link{deunif}}.


  }
}
\details{


General details are given in \code{\link{deunif}}
including
a note regarding the terminology used.
Here,
\code{exp} corresponds to the distribution of interest, \eqn{F}, and
\code{eexp} corresponds to \eqn{G}.
The addition of ``\code{e}'' is for the `other'
distribution associated with the parent distribution.
Thus
\code{deexp} is for \eqn{g},
\code{peexp} is for \eqn{G},
\code{qeexp} is for the inverse of \eqn{G},
\code{reexp} generates random variates from \eqn{g}.


For \code{qeexp} the Newton-Raphson algorithm is used to solve for
\eqn{y} satisfying \eqn{p = G(y)}.
Numerical problems may occur when values of \code{p} are
very close to 0 or 1.


}
\value{
  \code{deexp(x)} gives the density function \eqn{g(x)}.
  \code{peexp(q)} gives the distribution function \eqn{G(q)}.
  \code{qeexp(p)} gives the expectile function:
  the value \eqn{y} such that \eqn{G(y)=p}.
  \code{reexp(n)} gives \eqn{n} random variates from \eqn{G}.


}

%\references{
%
%Jones, M. C. (1994)
%Expectiles and M-quantiles are quantiles.
%\emph{Statistics and Probability Letters},
%\bold{20}, 149--153.
%
%}
\author{ T. W. Yee and Kai Huang }

%\note{
%The ``\code{q}'', as the first character of ``\code{qeunif}'',
%may be changed to ``\code{e}'' in the future,
%the reason being to emphasize that the expectiles are returned.
%Ditto for the argument ``\code{q}'' in \code{peunif}.
%
%}

\seealso{
  \code{\link{deunif}},
  \code{\link{denorm}},
  \code{\link{dexp}}.


}

\examples{
my.p <- 0.25; y <- rexp(nn <- 1000)
(myexp <- qeexp(my.p))
sum(myexp - y[y <= myexp]) / sum(abs(myexp - y))  # Should be my.p

\dontrun{ par(mfrow = c(2,1))
yy <- seq(-0, 4, len = nn)
plot(yy, deexp(yy),  col = "blue", ylim = 0:1, xlab = "y", ylab = "g(y)",
     type = "l", main = "g(y) for Exp(1); dotted green is f(y) = dexp(y)")
lines(yy, dexp(yy), col = "darkgreen", lty = "dotted", lwd = 2)  # 'original'

plot(yy, peexp(yy), type = "l", col = "blue", ylim = 0:1,
     xlab = "y", ylab = "G(y)", main = "G(y) for Exp(1)")
abline(v = 1, h = 0.5, col = "red", lty = "dashed")
lines(yy, pexp(yy), col = "darkgreen", lty = "dotted", lwd = 2) }
}
\keyword{distribution}

%# Equivalently:
%I1 <- mean(y <= myexp) * mean( myexp - y[y <= myexp])
%I2 <- mean(y >  myexp) * mean(-myexp + y[y >  myexp])
%I1 / (I1 + I2)  # Should be my.p
%# Or:
%I1 <- sum( myexp - y[y <= myexp])
%I2 <- sum(-myexp + y[y >  myexp])


%# Non-standard exponential
%myrate <- 8
%yy <- rexp(nn, rate = myrate)
%(myexp <- qeexp(my.p, rate = myrate))
%sum(myexp - yy[yy <= myexp]) / sum(abs(myexp - yy))  # Should be my.p
%peexp(-Inf, rate = myrate)      #  Should be 0
%peexp( Inf, rate = myrate)      #  Should be 1
%peexp(mean(yy), rate = myrate)  #  Should be 0.5
%abs(qeexp(0.5, rate = myrate) - mean(yy))  #  Should be 0
%abs(peexp(myexp, rate = myrate) - my.p)  #  Should be 0
%integrate(f = deexp, lower = -1, upper = Inf, rate = myrate)  #  Should be 1









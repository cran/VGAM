\name{Expectiles-Normal}
\alias{Expectiles-Normal}
\alias{enorm}
\alias{denorm}
\alias{penorm}
\alias{qenorm}
\alias{renorm}
\title{ Expectiles of the Normal Distribution }
\description{
  Density function, distribution function, and
  expectile function and random generation for the distribution
  associated with the expectiles of a normal distribution.


}
\usage{
denorm(x, mean = 0, sd = 1, log = FALSE)
penorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
qenorm(p, mean = 0, sd = 1, Maxit.nr = 10, Tol.nr = 1.0e-6,
       lower.tail = TRUE, log.p = FALSE)
renorm(n, mean = 0, sd = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, p, q}{
  See \code{\link{deunif}}.


  }
  \item{n, mean, sd, log}{
  See \code{\link[stats:Normal]{rnorm}}.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


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
\code{norm} corresponds to the distribution of interest, \eqn{F}, and
\code{enorm} corresponds to \eqn{G}.
The addition of ``\code{e}'' is for the `other'
distribution associated with the parent distribution.
Thus
\code{denorm} is for \eqn{g},
\code{penorm} is for \eqn{G},
\code{qenorm} is for the inverse of \eqn{G},
\code{renorm} generates random variates from \eqn{g}.



For \code{qenorm} the Newton-Raphson algorithm is used to solve for
\eqn{y} satisfying \eqn{p = G(y)}.
Numerical problems may occur when values of \code{p} are
very close to 0 or 1.



}
\value{
  \code{denorm(x)} gives the density function \eqn{g(x)}.
  \code{penorm(q)} gives the distribution function \eqn{G(q)}.
  \code{qenorm(p)} gives the expectile function:
  the value \eqn{y} such that \eqn{G(y)=p}.
  \code{renorm(n)} gives \eqn{n} random variates from \eqn{G}.


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
  \code{\link{deexp}},
  \code{\link{dnorm}},
  \code{\link{amlnormal}},
  \code{\link{lms.bcn}}.


}

\examples{
my.p <- 0.25; y <- rnorm(nn <- 1000)
(myexp <- qenorm(my.p))
sum(myexp - y[y <= myexp]) / sum(abs(myexp - y))  # Should be my.p

# Non-standard normal
mymean <- 1; mysd <- 2
yy <- rnorm(nn, mymean, mysd)
(myexp <- qenorm(my.p, mymean, mysd))
sum(myexp - yy[yy <= myexp]) / sum(abs(myexp - yy))  # Should be my.p
penorm(-Inf, mymean, mysd)      #  Should be 0
penorm( Inf, mymean, mysd)      #  Should be 1
penorm(mean(yy), mymean, mysd)  #  Should be 0.5
abs(qenorm(0.5, mymean, mysd) - mean(yy))  #  Should be 0
abs(penorm(myexp, mymean, mysd) - my.p)    #  Should be 0
integrate(f = denorm, lower = -Inf, upper = Inf,
          mymean, mysd)  #  Should be 1

\dontrun{
par(mfrow = c(2, 1))
yy <- seq(-3, 3, len = nn)
plot(yy, denorm(yy), type = "l", col="blue", xlab = "y", ylab = "g(y)",
     main = "g(y) for N(0,1); dotted green is f(y) = dnorm(y)")
lines(yy, dnorm(yy), col = "darkgreen", lty = "dotted", lwd = 2)  # 'original'

plot(yy, penorm(yy), type = "l", col = "blue", ylim = 0:1,
     xlab = "y", ylab = "G(y)", main = "G(y) for N(0,1)")
abline(v = 0, h = 0.5, col = "red", lty = "dashed")
lines(yy, pnorm(yy), col = "darkgreen", lty = "dotted", lwd = 2) }
}
\keyword{distribution}

%# Equivalently:
%I1 = mean(y <= myexp) * mean( myexp - y[y <= myexp])
%I2 = mean(y >  myexp) * mean(-myexp + y[y >  myexp])
%I1 / (I1 + I2)  # Should be my.p
%# Or:
%I1 = sum( myexp - y[y <= myexp])
%I2 = sum(-myexp + y[y >  myexp])






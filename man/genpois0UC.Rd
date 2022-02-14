\name{Genpois0}
%\alias{dgenpois}
 \alias{Genpois0}
 \alias{dgenpois0}
 \alias{pgenpois0}
 \alias{qgenpois0}
 \alias{rgenpois0}
\title{Generalized Poisson Distribution
  (Original Parameterization)}
\description{
  Density, distribution function, quantile function and
  random generation
  for the original parameterization of the
  generalized Poisson distribution.

}
% dgenpois(x, lambda = 0, theta, log = FALSE)  % 20200808; withdrawn
\usage{
dgenpois0(x, theta, lambda = 0, log = FALSE)
pgenpois0(q, theta, lambda = 0, lower.tail = TRUE)
qgenpois0(p, theta, lambda = 0)
rgenpois0(n, theta, lambda = 0, algorithm = c("qgenpois0",
          "inv", "bup","chdn", "napp", "bran"))
}
\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities. }
  \item{n}{Similar to \code{\link[stats]{runif}}.}


  \item{theta, lambda}{
    See \code{\link{genpoisson0}}.
    The default value of \code{lambda} corresponds to an
    ordinary Poisson distribution.
    \emph{Nonnegative} values
    of \code{lambda} are currently required.


%    ---but this might change in the future.


  }
  \item{lower.tail, log}{
  Similar to \code{\link[stats]{Poisson}}.


  }
  \item{algorithm}{
  Character.
  Six choices are available, standing for the
  \emph{qgenpois0},
  \emph{inversion}, \emph{build-up}, \emph{chop-down},
  \emph{normal approximation} and \emph{branching} methods.
  The first one is the default and
  calls \code{qgenpois0} with \code{\link[stats]{runif}} as
  its first argument.
  The value inputted may be abbreviated, e.g., \code{alg = "n"}.
  The last 5 algorithms are a direct implementation of Demirtas (2017)
  and the relative performance of the algorithms are
  described there---however, the vectorization here may render
  the comments on relative speed as no longer holding.
  
  

  }
}
\value{
  \code{dgenpois0} gives the density,
  \code{pgenpois0} gives the distribution function,
  \code{qgenpois0} gives the quantile function, and
  \code{rgenpois} generates random deviates.
  For some of these functions such as
  \code{dgenpois0} and \code{pgenpois0}
  the value \code{NaN} is returned for elements not satisfying
  the parameter restrictions, e.g., if \eqn{\lambda > 1}.
  For some of these functions such as
  \code{rgenpois0}
  the input must not contain \code{NA}s or \code{NaN}s, etc. since
  the implemented algorithms are fragile.

  

% for the last 5 algorithms



}
\author{ T. W. Yee.
  For \code{rgenpois0()} the last 5 algorithms are based on
  code written in H. Demirtas (2017) and vectorized by T. W. Yee;
  but the \code{"bran"} algorithm was rewritten from
  Famoye (1997).


  
}
\details{
  Most of the background to these functions are given
  in \code{\link{genpoisson0}}.
  Some warnings relevant to this distribution are given there.
  The complicated range of the
  parameter \code{lambda} when negative is no longer
  supported because the distribution is not normalized.
  For other GPD variants see \code{\link{Genpois1}}.

  

%  Note that numerical round off errors etc. can occur; see
%  below for an example.


  
}
\note{For \code{rgentpois0()}:
(1). \code{"inv"}, \code{"bup"} and \code{"chdn"} appear similar and
seem to work okay.
(2). \code{"napp"} works only when theta is large, away from 0.
It suffers from 0-inflation.
(3). \code{"bran"} has a relatively heavy RHS tail and
requires positive \code{lambda}.
More details can be found in
Famoye (1997) and
Demirtas (2017).



The function \code{dgenpois0}
uses \code{\link[base:Special]{lfactorial}}, which
equals \code{Inf} when \code{x} is approximately \code{1e306}
on many machines.
So the density is returned as \code{0} in very extreme cases;
see \code{\link[base]{.Machine}}.



}
\references{

Demirtas, H. (2017).
On accurate and precise generation of generalized Poisson variates.
\emph{Communications in Statistics---Simulation and Computation},
\bold{46}, 489--499.


Famoye, F. (1997).
Generalized Poisson random variate generation.
\emph{Amer. J. Mathematical and Management Sciences},
\bold{17}, 219--237.



}
\section{Warning }{
  These have not been tested thoroughly.



% The default value of \code{algorithm} might change
% in the future, depending on experience gained by using it
% over a long period of time.

  
For \code{pgentpois0()} \code{\link[base]{mapply}} is called
with \code{0:q} as input, hence will be very slow and
memory-hungry for large values of \code{q}.
Likewise \code{qgentpois0()} and \code{rgentpois0()}
may suffer from the same limitations.




}

\seealso{
  \code{\link{genpoisson0}},
  \code{\link{Genpois1}},
  \code{\link[stats:Poisson]{dpois}}.


}
\examples{
sum(dgenpois0(0:1000, theta = 2, lambda = 0.5))
\dontrun{theta <- 2; lambda <- 0.2; y <- 0:10
proby <- dgenpois0(y, theta = theta, lambda = lambda, log = FALSE)
plot(y, proby, type = "h", col = "blue", lwd = 2, ylab = "Pr(Y=y)",
     main = paste0("Y ~ GP-0(theta=", theta, ", lambda=",
                   lambda, ")"), las = 1, ylim = c(0, 0.3),
     sub = "Orange is the Poisson probability function")
lines(y + 0.1, dpois(y, theta), type = "h", lwd = 2, col = "orange") }
}
\keyword{distribution}


%sum(dgenpois(0:1000, lambda = -0.5, theta = 2))  # Not perfect...

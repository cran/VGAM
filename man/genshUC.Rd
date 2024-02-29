\name{Gensh}
\alias{Gensh}
\alias{dgensh}
\alias{pgensh}
\alias{qgensh}
\alias{rgensh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Secant Hyperbolic
  Distribution



}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generalized secant hyperbolic
  distribution.


}
\usage{
dgensh(x, shape, location = 0, scale = 1, tol0 = 1e-4,
       log = FALSE)
pgensh(q, shape, location = 0, scale = 1, tol0 = 1e-4,
       lower.tail = TRUE)
qgensh(p, shape, location = 0, scale = 1, tol0 = 1e-4)
rgensh(n, shape, location = 0, scale = 1, tol0 = 1e-4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log, lower.tail}{
  Similar meaning as in \code{\link{Normal}}.

  }

  \item{shape}{Numeric.
    Shape parameter, called \eqn{t} in
    Vaughan (2002).
    Valid values are
    \eqn{-\pi/2 < t}{-0.5 * pi < t}.
  
 
  }
  \item{location, scale}{Numeric.
    The location and (positive) scale
    parameters.


  }
  \item{tol0}{Numeric.
    Used to test whether the shape parameter
    is close enough to be treated as 0.

 

  }

  
}
\details{
  This is an implementation of the family of
  symmetric densities described
  by Vaughan (2002).
  By default, the mean and variance are 0
  and 1, for all \eqn{t}.
  Some special (default) cases are:
  \eqn{t=0}: logistic
  (which is similar to
  \code{\link[stats:TDist]{stats:dt}}
  with 9 degrees of freedom);
  \eqn{t=-\pi/2}{t=pi/2}: the standard secant
  hyperbolic (whence the name);
  \eqn{t=\infty}{t=Inf}:
  uniform(-sqrt(3), sqrt(3)).


}
\section{Warning }{
  Numerical problems may occur when some
  argument values are extreme.




}



\value{
  \code{dgensh} gives the density,
  \code{pgensh} gives the distribution function,
  \code{qgensh} gives the quantile function, and
  \code{rgensh} generates random deviates.
  


}
%\references{


% Vaughan, D. C. (2002).
% The generalized secant hyperbolic
% distribution and its properties.
% \emph{Communications in Statistics---Theory
%   and Methods},
% \bold{31}(2): 219--238.



%}

\author{ T. W. Yee. }
%\note{
%  See \code{\link{Gaitdpois}} for
%  general information also relevant
%  to this parent distribution.
  

%}




% \code{\link{gaitpoisson.mlm}},


\seealso{
  \code{\link{gensh}},
  \code{\link{logistic}},
  \code{\link{hypersecant}},
  \code{\link[stats:Logistic]{Logistic}}.




}
\examples{
x <- seq(-2, 4, by = 0.01)
loc <- 1; shape <- -pi /2
\dontrun{plot(x, dgensh(x, shape, loc), type = "l",
     main = "Blue is density, orange is the CDF",
     ylim = 0:1, las = 1, ylab = "", 
     sub = "Purple are 5, 10, ..., 95 percentiles",
     col = "blue")
abline(h = 0, col = "blue", lty = 2)
lines(qgensh((1:19) / 20, shape, loc), type = "h",
      dgensh(qgensh((1:19) / 20, shape, loc),
             shape, loc), col = "purple", lty = 3)
lines(x, pgensh(x, shape, loc), col = "orange")
abline(h = 0, lty = 2) }

pp <- (1:19) / 20  # Test two functions
max(abs(pgensh(qgensh(pp, shape, loc),
               shape,loc) - pp))  # Should be 0
}
\keyword{distribution}







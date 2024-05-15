\name{N1pois}
\alias{dN1pois}
\alias{rN1pois}
\title{ Linear Model and Poisson
        Mixed Data Type
        Distribution}
\description{
  Density,
  and random generation
  for the (four parameter bivariate)
  Linear Model--Poisson copula distribution.


% distribution function,

}
\usage{
dN1pois(x1, x2, mean = 0, sd = 1, lambda, apar = 0,
        doff = 5, copula = "gaussian", log = FALSE)
rN1pois(n, mean = 0, sd = 1, lambda, apar = 0,
        doff = 5, copula = "gaussian")
}
\arguments{
  \item{x1, x2}{vector of quantiles.
  The valid values of \code{x2} are nonnegative
  integers.


  }
  \item{n}{number of observations.
    Same as \code{\link[stats]{rnorm}}.

  }
  \item{copula}{
   See \code{\link{N1poisson}}.

  }
  \item{mean, sd, lambda, apar}{
   See \code{\link{N1poisson}}.


  }
  \item{doff}{
   See \code{\link{N1poisson}}.


  }
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm is returned.

%   Same as \code{\link[stats]{rnorm}}.

  }
}
\value{
  \code{dN1pois} gives the probability density/mass function,
  \code{rN1pois} generates random deviate and returns
  a two-column matrix.


% \code{pN1pois} gives the distribution function, and

}
%\references{
%
%}

\author{ T. W. Yee }
\details{
  See \code{\link{N1poisson}}, the \pkg{VGAM}
  family functions for estimating the
  parameter by maximum likelihood estimation,
  for details.


}
%\note{


%}
\seealso{
  \code{\link{N1poisson}},
  \code{\link[stats]{rnorm}},
  \code{\link[stats]{rpois}}.


}
\examples{
\dontrun{ 
nn <- 1000; mymu <- 1; sdev <- exp(1)
apar <- rhobitlink(0.4, inverse = TRUE)
lambda <- loglink(1, inverse = TRUE)
mat <- rN1pois(nn, mymu, sdev, lambda, apar)
pndata <- data.frame(y1 = mat[, 1], y2 = mat[, 2])
with(pndata, plot(jitter(y1), jitter(y2), col = 4))
}
}
\keyword{distribution}


%plot(r <- rN1pois(n = 3000, rho = Rho), col = "blue")
%par(mfrow = c(1, 2))
%hist(r[, 1])  # Should be uniform
%hist(r[, 2])  # Should be uniform

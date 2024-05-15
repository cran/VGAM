\name{N1binom}
\alias{dN1binom}
\alias{rN1binom}
\title{ Linear Model and Binomial
        Mixed Data Type
        Distribution}
\description{
  Density,
  and random generation
  for the (four parameter bivariate)
  Linear Model--Bernoulli copula distribution.


% distribution function,

}
\usage{
dN1binom(x1, x2, mean = 0, sd = 1, prob, apar = 0,
         copula = "gaussian", log = FALSE)
rN1binom(n, mean = 0, sd = 1, prob,
         apar = 0, copula = "gaussian")
}
\arguments{
  \item{x1, x2}{vector of quantiles.
  The valid values of \code{x2} are \eqn{0}
  and \eqn{1}.


  }
  \item{n}{number of observations.
    Same as \code{\link[stats]{rnorm}}.

  }
  \item{copula}{
   See \code{\link{N1binomial}}.

  }
  \item{mean, sd, prob, apar}{
   See \code{\link{N1binomial}}.


  }
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm is returned.

%   Same as \code{\link[stats]{rnorm}}.

  }
}
\value{
  \code{dN1binom} gives the probability density/mass function,
  \code{rN1binom} generates random deviate and returns
  a two-column matrix.


% \code{pN1binom} gives the distribution function, and

}
%\references{
%
%}

\author{ T. W. Yee }
\details{
  See \code{\link{N1binomial}}, the \pkg{VGAM}
  family functions for estimating the
  parameter by maximum likelihood estimation,
  for details.


}
%\note{


%}
\seealso{
  \code{\link{N1binomial}},
  \code{\link[stats]{rnorm}},
  \code{\link[stats]{rbinom}}.


}
\examples{
\dontrun{ 
nn <- 1000;  apar <- rhobitlink(1.5, inverse = TRUE)
prob <- logitlink(0.5, inverse = TRUE)
mymu <- 1; sdev <- exp(1)
mat <- rN1binom(nn, mymu, sdev, prob, apar)
bndata <- data.frame(y1 = mat[, 1], y2 = mat[, 2])
with(bndata, plot(jitter(y1), jitter(y2), col = "blue"))
}
}
\keyword{distribution}


%plot(r <- rN1binom(n = 3000, rho = Rho), col = "blue")
%par(mfrow = c(1, 2))
%hist(r[, 1])  # Should be uniform
%hist(r[, 2])  # Should be uniform

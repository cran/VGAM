\name{Extbetabinom}
\alias{Extbetabinom}
\alias{dextbetabinom}
\alias{pextbetabinom}
\alias{qextbetabinom}
\alias{rextbetabinom}
\title{The Beta-Binomial Distribution}
\description{
  Density, distribution function,
  quantile function and random generation
  for the
  extended beta-binomial distribution.


}
\usage{
dextbetabinom(x, size, prob, rho = 0,
     log = FALSE, forbycol = TRUE)
pextbetabinom(q, size, prob, rho = 0,
     lower.tail = TRUE, forbycol = TRUE)
qextbetabinom(p, size, prob, rho = 0,
     forbycol = TRUE)
rextbetabinom(n, size, prob, rho = 0)
}



\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{size}{number of trials.}
  \item{n}{number of observations.
  Same as \code{\link[stats]{runif}}.


  }
  \item{prob}{
  the probability of success \eqn{\mu}{mu}.
  Must be in the unit closed interval \eqn{[0,1]}.


  }
  \item{rho}{
    the correlation parameter \eqn{\rho}{rho},
    which
    may be negative for underdispersion or
    else be in the interval \eqn{[0, 1)}
    for overdispersion.
  The default value of 0 corresponds to the
  usual binomial distribution with
  probability \code{prob}.



  }
  \item{log, lower.tail}{
  Same meaning as \code{\link[stats]{runif}}.

% log.p, 



  }
  \item{forbycol}{
    Logical.
    A \code{for} loop cycles over either the rows or
    columns and this argument determines which.
    The rows are \code{1:length(x)} and the
    columns are \code{0:max(size)}.
    The best choice is data set dependent.



  }

}
\value{
  \code{dextbetabinom} gives the density,
  \code{pextbetabinom} gives the
  distribution function,
  \code{qextbetabinom} gives the quantile function
  and
  \code{rextbetabinom} generates random
  deviates.


}
%\author{ T. W. Yee}
\details{
  The \emph{extended} beta-binomial
  distribution allows for a slightly negative
  correlation parameter between binary
  responses within a cluster (e.g., a litter).
  An exchangeable error structure with
  correlation \eqn{\rho} is assumed.
  

}
\note{
  Currently most of the code is quite slow.
  Speed improvements are a future project.
  Use \code{forbycol} optimally.


}
\section{Warning }{

  Setting \code{rho = 1} is not recommended
  as \code{NaN} is returned,
  however the code may be
  modified in the future to handle this
  special case.


}


\seealso{
  \code{\link{extbetabinomial}},
  \code{\link{Betabinom}},
  \code{\link[stats:Binomial]{Binomial}}.


}
\examples{
set.seed(1); rextbetabinom(10, 100, 0.5)
set.seed(1);        rbinom(10, 100, 0.5)  # Same

\dontrun{N <- 9; xx <- 0:N; prob <- 0.5; rho <- -0.02
dy <- dextbetabinom(xx, N, prob, rho)
barplot(rbind(dy, dbinom(xx, size = N, prob)),
  beside = TRUE, col = c("blue","green"), las = 1,
  main = paste0("Beta-binom(size=", N, 
  ", prob=", prob, ", rho=", rho, ") (blue) vs\n",
  " Binom(size=", N, ", prob=", prob, ") (green)"),
  names.arg = as.character(xx), cex.main = 0.8)
sum(dy * xx)  # Check expected values are equal
sum(dbinom(xx, size = N, prob = prob) * xx)
cumsum(dy) - pextbetabinom(xx, N, prob, rho)  # 0?
}
}
\keyword{distribution}





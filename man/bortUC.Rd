\name{Bort}
\alias{Bort}
\alias{dbort}
%\alias{pbort}
%\alias{qbort}
\alias{rbort}
\title{The Borel-Tanner Distribution}
\description{
  Density
  and random generation for the Borel-Tanner distribution.


% distribution function, quantile function


}
\usage{
dbort(x, Qsize = 1, a = 0.5, log = FALSE)
rbort(n, Qsize = 1, a = 0.5)
}

%pbort(q, Qsize = 1, a = 0.5)
%qbort(p, Qsize = 1, a = 0.5)

\arguments{
  \item{x}{vector of quantiles.}
% \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{Qsize, a}{
  See \code{\link{borel.tanner}}.


  }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.



  }

}
\value{
  \code{dbort} gives the density,
  \code{rbort} generates random deviates.


% \code{pbort} gives the distribution function,
% \code{qbort} gives the quantile function, and


}
\author{ T. W. Yee }
\details{
  See \code{\link{borel.tanner}}, the \pkg{VGAM} family function
  for estimating the parameter,
  for the formula of the probability density function and other details.


}
\section{Warning }{
  Looping is used for \code{\link{rbort}}, therefore
  values of \code{a} close to 1 will result in long (or infinite!)
  computational times.
  The default value of \code{a} is subjective.


}
\seealso{
  \code{\link{borel.tanner}}.


}
\examples{
\dontrun{ qsize <- 1; a <- 0.5; x <- qsize:(qsize+10)
plot(x, dbort(x, qsize, a), type = "h", las = 1, col = "blue",
     ylab = paste("fbort(qsize=", qsize, ", a=", a, ")"), log = "y",
     main = "Borel-Tanner density function") }
}
\keyword{distribution}



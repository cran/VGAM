\name{rdiric}
\alias{rdiric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Dirichlet distribution }
\description{
  Generates Dirichlet random variates.

}
\usage{
rdiric(n, shape, dimension = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of observations. }
  \item{shape}{
  the shape parameters. These must be positive.
  If \code{dimension} is specifed, values
  are recycled if necessary to length \code{dimension}.

  }
  \item{dimension}{
  the dimension of the distribution. 
  If \code{dimension} is not numeric then it is taken to be
  \code{length(shape)}.

  }
}
\details{
  This function is based on a relationship between the gamma and
  Dirichlet distribution. Random gamma variates are generated, and
  then Dirichlet random variates are formed from these. 

}
\value{
  A \code{n} by \code{dimension} matrix of Dirichlet random variates.
  Each element is positive, and each row will sum to unity.

} 
 
\references{ 
Lange, K. (2002)
\emph{Mathematical and Statistical Methods for Genetic Analysis},
2nd ed.
New York: Springer-Verlag.

}
\author{ Thomas W. Yee }
\seealso{ 
  \code{\link{dirichlet}} is a \pkg{VGAM} family function for
  fitting a Dirichlet distribution to data.

} 

\examples{
y <- rdiric(n = 1000, shape = c(3, 1, 4))
fit <- vglm(y ~ 1, dirichlet, trace = TRUE, crit = "c")
Coef(fit)
coef(fit, matrix = TRUE)
}
\keyword{distribution}


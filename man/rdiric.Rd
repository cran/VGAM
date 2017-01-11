\name{rdiric}
\alias{rdiric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Dirichlet distribution }
\description{
  Generates Dirichlet random variates.

}
\usage{
rdiric(n, shape, dimension = NULL, is.matrix.shape = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of observations.
  Note it has two meanings, see \code{is.matrix.shape} below.



  }
  \item{shape}{
  the shape parameters. These must be positive.
  If \code{dimension} is specifed, values
  are recycled if necessary to length \code{dimension}.


  }
  \item{dimension}{
  the dimension of the distribution.
  If \code{dimension} is not numeric then it is taken to be
  \code{length(shape)}
  (or \code{ncol(shape)} if \code{is.matrix.shape == TRUE}).


  }
  \item{is.matrix.shape}{
  Logical.
  If \code{TRUE} then \code{shape} must be a matrix, and then
  \code{n} is no longer the number of rows of the answer but the
  answer has \code{n * nrow(shape)} rows.
  If \code{FALSE} (the default) then \code{shape} is a vector and each
  of the \code{n} rows of the answer have \code{shape} as
  its shape parameters.


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
  If \code{shape} has names then these will become the column names
  of the answer.


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
ddata <- data.frame(rdiric(n = 1000, shape = c(y1 = 3, y2 = 1, y3 = 4)))
fit <- vglm(cbind(y1, y2, y3) ~ 1, dirichlet, data = ddata, trace = TRUE)
Coef(fit)
coef(fit, matrix = TRUE)
}
\keyword{distribution}


\name{iam}
\alias{iam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Index from Array to Matrix }
\description{
  Maps the elements of an array containing symmetric positive-definite
  matrices to a matrix with sufficient columns to hold them
  (called matrix-band format.)

}
\usage{
iam(j, k, M, both = FALSE, diag = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{
  An integer from the set \{\code{1:M}\} giving the row number
  of an element.

  }
  \item{k}{
  An integer from the set \{\code{1:M}\} giving the column number
  of an element.

  }
  \item{M}{
  The number of linear/additive predictors. This is the
  dimension of each positive-definite symmetric matrix.

  }
  \item{both}{
  Logical. Return both the row and column indices?
  See below for more details.

  }
  \item{diag}{
  Logical. Return the indices for the diagonal elements?
  If \code{FALSE} then only the strictly upper triangular part of the matrix
  elements are used.

  }
}
\details{
  Suppose we have \eqn{n} symmetric positive-definite square matrices,
  each \eqn{M} by \eqn{M}, and
  these are stored in an \code{array} of dimension \code{c(n,M,M)}.
  Then these can be more compactly represented by a
  \code{matrix} of dimension \code{c(n,K)} where \code{K} is an integer
  between \code{M} and \code{M*(M+1)/2} inclusive.
  The mapping between these two representations is given by this function.
  It firstly enumerates by the diagonal elements, followed by the band
  immediately above the diagonal, then the band above that one, etc.
  The last element is \code{(1,M)}.
  This function performs the mapping from elements \code{(j,k)}
  of symmetric positive-definite square matrices to the columns of
  another matrix representing such.
  This is called the \emph{matrix-band} format and is used by
  the \pkg{VGAM} package.


}
\value{
  This function has a dual purpose depending on the value of \code{both}.
  If \code{both=FALSE} then the column number corresponding
  to the \code{j}-\code{k} element of the matrix is returned.
  If \code{both = TRUE} then \code{j} and \code{k} are ignored and a list
  with the following components are returned.


  \item{row.index}{
  The row indices of the upper triangular part of the
  matrix (This may or may not include the diagonal elements, depending
  on the argument \code{diagonal}).


  }
  \item{col.index}{
  The column indices of the upper triangular part of the
  matrix (This may or may not include the diagonal elements, depending
  on the argument \code{diagonal}).


  }
}
%\references{
%  The website \url{http://www.stat.auckland.ac.nz/~yee} contains
%  some additional information.
%
%
%}
\author{ T. W. Yee }
\note{
  This function is used in the \code{weight} slot of many \pkg{VGAM}
  family functions (see \code{\link{vglmff-class}}), especially those
  whose \eqn{M} is determined by the data, e.g., \code{\link{dirichlet}},
  \code{\link{multinomial}}.


}

\seealso{
  \code{\link{vglmff-class}}.


%\code{ima}.


}
\examples{
iam(1, 2, M = 3)  # The 4th column represents element (1,2) of a 3x3 matrix
iam(NULL, NULL, M = 3, both = TRUE)  # Return the row and column indices

dirichlet()@weight

M <- 4
temp1 <- iam(NA, NA, M = M, both = TRUE)
mat1 <- matrix(NA, M, M)
mat1[cbind(temp1$row, temp1$col)] = 1:length(temp1$row)
mat1  # More commonly used

temp2 <- iam(NA, NA, M = M, both = TRUE, diag = FALSE)
mat2 <- matrix(NA, M, M)
mat2[cbind(temp2$row, temp2$col)] = 1:length(temp2$row)
mat2  # Rarely used
}
\keyword{manip}
\keyword{programming}

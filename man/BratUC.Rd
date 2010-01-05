\name{Brat}
\alias{Brat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Inputting Data to fit a Bradley Terry Model }
\description{
  Takes in a square matrix of counts and outputs
  them in a form that is accessible to the \code{\link{brat}}
  and \code{\link{bratt}} family functions.
}
\usage{
Brat(mat, ties=0*mat, string=c(" > "," == "))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{ Matrix of counts,
which is considered \eqn{M} by \eqn{M} in dimension when
there are ties, and \eqn{M+1} by \eqn{M+1} 
when there are no ties.
The rows are winners and the columns are losers, e.g.,
the 2-1 element is now many times Competitor 2 has beaten
Competitor 1.
The matrices are best labelled with the competitors' names.
}
  \item{ties}{ Matrix of counts. This should be the same
dimension as \code{mat}. By default, there are no ties.
The matrix must be symmetric, and the diagonal should contain
\code{NA}s. 
}
  \item{string}{ Character. 
The matrices are labelled with the first value of the descriptor, e.g.,
\code{"NZ > Oz"} `means' NZ beats Australia in rugby.
Suggested alternatives include \code{" beats "} or \code{" wins against "}.
The second value is used to handle ties.

}
}
\details{
  In the \pkg{VGAM} package it is necessary for each matrix to
  be represented as a single row of data by \code{\link{brat}} and
  \code{\link{bratt}}.  Hence the non-diagonal elements of the \eqn{M+1}
  by \eqn{M+1} matrix are concatenated into \eqn{M(M+1)} values (no ties),
  while if there are ties, the non-diagonal elements of the \eqn{M}
  by \eqn{M} matrix are concatenated into \eqn{M(M-1)} values.

}
\value{
  A matrix with 1 row and either \eqn{M(M+1)} or \eqn{M(M-1)} columns. 
}
\references{
Agresti, A. (2002)
\emph{Categorical Data Analysis},
2nd ed. New York: Wiley.
}
\author{ T. W. Yee }
\note{
This is a data preprocessing function for
\code{\link{brat}} and \code{\link{bratt}}.

Yet to do: merge \code{InverseBrat} into \code{brat}.
}
\seealso{ 
\code{\link{brat}},
\code{\link{bratt}},
\code{InverseBrat}.
}
\examples{
journal = c("Biometrika", "Comm Statist", "JASA", "JRSS-B")
m = matrix(c( NA, 33, 320, 284,   730, NA, 813, 276,
             498, 68,  NA, 325,   221, 17, 142, NA), 4, 4)
dimnames(m) = list(winner = journal, loser = journal)
Brat(m)
vglm(Brat(m) ~ 1, brat, trace=TRUE)
}
\keyword{models}
\keyword{regression}


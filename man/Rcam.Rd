\name{Rcam}
\alias{Rcam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Mark the baseline of row and column on a matrix data 

}
\description{
  This function will grab the baseline of row and column for
  the rank-zero model of row-column association models or main
  effects models.

}
\usage{
  Rcam(mat, rbaseline = 1, cbaseline = 1) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{mat}{
  Matrix of \eqn{r} by \eqn{c}.
  The matrix is better labeled with row and column names.

}
\item{rbaseline}{
  Numeric (row number of the matrix y) or can be 
  character of row names of the matrix \code{y} that user
  wants it as row baseline or reference. 
  Similarly \code{cbaseline} for the column

}
\item{cbaseline}{
 Numeric of column number of the matrix \code{mat}, or can be 
 Character of column names of the matrix \code{mat} that user
 wants it as column baseline or reference.

}
}
\details{
  In some regression analysis, it is necessary to put the baseline
  (or reference) of the factor, and as well in rank-zero of
  row-column association model. It should be done by marking
  a baseline on a level of row and column.  This is a data
  preprocessing function for \code{\link{rcam}}.

}

\value{
  Matrix of the same dimension as the input, rearranged the first
  row or column.

}
\author{
Alfian F. Hadi.

}
\note{
This is a data preprocessing function for \code{\link{rcam}}.

}


\seealso{
  \code{\link{moffset}},
  \code{\link{rcam}},
  \code{\link{plotrcam0}}.

}
\examples{
(alcoff.e <- moffset(alcoff, roffset = "6"))
Rcam(alcoff.e, rbaseline = "11", cbaseline = "Monday")
}

\name{Rcim}
\alias{Rcim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Mark the Baseline of Row and Column on a Matrix data

}
\description{
  Rearrange the rows and columns of the input so
  that the first row and first column are baseline.
  This function is for rank-zero row-column interaction models
  (RCIMs; i.e., general main effects models).


}
\usage{
  Rcim(mat, rbaseline = 1, cbaseline = 1)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{mat}{
  Matrix, of dimension \eqn{r} by \eqn{c}.
  It is best that it is labelled with row and column names.


}
\item{rbaseline, cbaseline}{
  Numeric (row number of the matrix \code{mat}) or
  character (matching a row name of \code{mat}) that the user
  wants as the row baseline or reference level.
  Similarly \code{cbaseline} for the column.


}
}
\details{
  This is a data preprocessing function for \code{\link{rcim}}.
  For rank-zero row-column interaction models this function
  establishes the baseline (or reference) levels of the matrix
  response with respect to the row and columns---these become
  the new first row and column.


}

\value{
  Matrix of the same dimension as the input,
  with \code{rbaseline} and \code{cbaseline} specifying the
  first rows and columns.
  The default is no change in \code{mat}.


}
\author{
Alfian F. Hadi and T. W. Yee.


}
\note{
  This function is similar to \code{\link{moffset}}; see
  \code{\link{moffset}} for information about the differences.
  If numeric, the arguments
  \code{rbaseline}
  and
  \code{cbaseline}
  differ from arguments
  \code{roffset}
  and
  \code{coffset}
  in \code{\link{moffset}}
  by 1 (when elements of the matrix agree).


}


\seealso{
  \code{\link{moffset}},
  \code{\link{rcim}},
  \code{\link{plotrcim0}}.


}
\examples{
(alcoff.e <- moffset(alcoff, roffset = "6", postfix = "*"))
(aa <- Rcim(alcoff,    rbaseline = "11", cbaseline = "Sun"))
(bb <- moffset(alcoff,             "11",             "Sun", postfix = "*"))
aa - bb  # Note the difference!
}

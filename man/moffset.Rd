\name{moffset}
\alias{moffset}

\title{
Makes an effective day data set 

}
\description{
  Modify a matrix by shitting successive elements.

}
\usage{
moffset(mat, roffset = 1, coffset = 1)

}
\arguments{
\item{mat}{
  Data frame or matrix, which is considered to 24 hours of rows by 
  7 days a week of columns, or can be any matrix.

}
\item{roffset, coffset}{
  Numeric or character. The number or name of the row/column
  of matrix input. rows/column start as the first cell in
  matrix output. For example, for the \code{\link{alcoff}},
  put \code{roffset = "6"} means that we make an effective days
  dataset start from 6:00 am of the day, and from midnight to
  05.59 am is part the previous days. Similarly \code{coffset}
  for the column.

}
}

\details{
  In \pkg{VGAM} there are crash data (e.g., see
  \code{\link{crashi}}) and alcohol-related related
  datasets, and also function of \code{\link{rcam}} and
  \code{\link{plotrcam0}} available. In regression analysis of
  main effect model or rank zero of row-column association model,
  it becomes more realistic using the effective days dataset. Since
  in the effective days datasets for example, partying on Friday
  night includes Saturday morning.


}

\value{
  A matrix of the same dimensional as its input.

}


\author{
  Alfian F. Hadi,
  T. W. Yee.

}
\note{
  This function was specifically  24 x 7 dimensional matrix (24
  hours of days of row by 7 days a week), \code{\link{alcoff}}.
  Change into a matrix of "effective days" with flexible start of
  first hour of in the effective day. But now can more generaly
  applied into any dimensional matrix.


  This is a data preprocessing function for \code{\link{rcam}}
  and \code{\link{plotrcam0}}.  The differences between
  \code{\link{Rcam}} and \code{\link{moffset}} is that the
  \code{\link{Rcam}} only reordering the level of row and column
  and marking the row/column base line, there is no change in
  entire rows and column.  But here in \code{\link{moffset}}
  values in row/column before the rstar/cstar will be moved to
  previuous row/column. So if there is very similar matrix to
  \code{\link{alcoff}}, but on other dimension 7 by 24, for example,
  then the user can use this function very easy.


}
\seealso{
  \code{\link{Rcam}}, 
  \code{\link{rcam}}, 
  \code{\link{plotrcam0}}.

}
\examples{
alcoff.e <- moffset(alcoff, roffset = "6")
fit0.orig <- rcam(alcoff) #default baselines are first row & col
fit0 <- rcam(alcoff.e)    #default baselines are first roe & col 

#customise rowbaseline
fit0.based <- rcam(Rcam(alcoff.e, rbaseline="11", cbaseline="Monday"))
} 

\name{moffset}
\alias{moffset}

\title{
Matrix Offset

}
\description{
  Modify a matrix by shifting successive elements.

}
\usage{
moffset(mat, roffset = 0, coffset = 0, postfix = "",
        rprefix = "Row.", cprefix = "Col.")
}
\arguments{
\item{mat}{
  Data frame or matrix.
  This ought to have at least three rows and three columns.
  The elements are shifted in the order of \code{c(mat)},
  i.e., going down successive columns, as the columns go
  from left to right. Wrapping of values is done.


}
\item{roffset, coffset}{
  Numeric or character.
  If numeric, the amount of shift (offset) for each row and column.
  The default is no change to \code{mat}.
  If character, the offset is computed by matching with
  the row or column names.
  For example, for the \code{\link{alcoff}},
  put \code{roffset = "6"} means that we make an effective day's
  dataset start from 6:00 am, and this wraps around to
  include midnight to 05.59 am on the next day.


}
\item{postfix}{
  Character.
  Modified rows and columns are renamed by pasting this argument
  to the end of each name.
  The default is no change.


}
\item{rprefix, cprefix}{
  Same as \code{\link{rcim}}.


}
}

\details{
  This function allows a matrix to be rearranged so that
  element (\code{roffset} + 1, \code{coffset} + 1)
  becomes the (1, 1) element.
  The elements are assumed to be ordered in the same way
  as the elements of \code{c(mat)},


  This function is applicable to, e.g.,
  \code{\link{alcoff}},
  where it is useful to define the \emph{effective day} as starting
  at some other hour than midnight, e.g., 6.00am.
  This is because partying on Friday night continues on into
  Saturday morning, therefore it is more interpretable to use
  the effective day when considering a daily effect.


  This is a data preprocessing function for \code{\link{rcim}}
  and \code{\link{plotrcim0}}.  The differences between
  \code{\link{Rcim}} and \code{\link{moffset}} is that
  \code{\link{Rcim}} only reorders the level of the rows and columns
  so that the data is shifted but not moved.
  That is, a value in one row stays in that row,
  and ditto for column.
  But in \code{\link{moffset}}
  values in one column can be moved to a previous column.
  See the examples below.


}

\value{
  A matrix of the same dimensional as its input.


}


\author{
  T. W. Yee,
  Alfian F. Hadi.


}
\note{
% This function was originally for a 24 x 7 dimensional matrix
% (24 hours of the day by 7 days per week) such as \code{\link{alcoff}}.
% Of course, this function can be applied to any moderately
% large matrix.


  The input \code{mat} should have row names and column names.


}
\seealso{
  \code{\link{Rcim}},
  \code{\link{rcim}},
  \code{\link{plotrcim0}},
  \code{\link{alcoff}},
  \code{\link{crashi}}.


}
\examples{
moffset(alcoff, 3, 2, "*")  # Some day's data is moved to previous day.
Rcim(alcoff, 3 + 1, 2 + 1)  # Data does not move as much.
alcoff  # Original data
moffset(alcoff, 3, 2, "*") - Rcim(alcoff, 3+1, 2+1)  # Note the differences

# An 'effective day' data set:
alcoff.e <- moffset(alcoff, roffset = "6", postfix = "*")
fit.o <- rcim(alcoff)    # default baselines are first row and col
fit.e <- rcim(alcoff.e)  # default baselines are first row and col

\dontrun{ par(mfrow = c(2, 2), mar = c(9, 4, 2, 1))
plot(fit.o, rsub = "Not very interpretable", csub = "Not very interpretable")
plot(fit.e, rsub = "More interpretable", csub = "More interpretable")
}

# Some checking
all.equal(moffset(alcoff), alcoff)  # Should be no change
moffset(alcoff, 1, 1, "*")
moffset(alcoff, 2, 3, "*")
moffset(alcoff, 1, 0, "*")
moffset(alcoff, 0, 1, "*")
moffset(alcoff, "6", "Mon", "*")  # This one is good

# Customise row and column baselines
fit2 <- rcim(Rcim(alcoff.e, rbaseline = "11", cbaseline = "Mon*"))
}

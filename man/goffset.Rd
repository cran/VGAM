\name{goffset}
\alias{goffset}

\title{
GAITD Offset for the GTE Method

}
\description{
  Utility function to create a matrix of log-offset values,
  to help facilitate the Generally-Truncated-Expansion method

}
\usage{
goffset(mux, n,
        a.mix = NULL, i.mix = NULL, d.mix = NULL,
        a.mlm = NULL, i.mlm = NULL, d.mlm = NULL, par1or2 = 1)
}
\arguments{
\item{mux}{
  Multiplier. Usually a small positive integer.
  Must be positive.
  The value 1 means no change.



}
\item{n}{
  Number of rows.
  A positive integer, it should be the number of rows of
  the data frame containing the data.


}
\item{a.mix, i.mix, d.mix}{
  See, e.g., \code{\link{gaitdpoisson}}.



}
\item{a.mlm, i.mlm, d.mlm}{
  See, e.g., \code{\link{gaitdpoisson}}.



}
\item{par1or2}{
  Number of parameters of the parent distribution.
  Set \code{par1or2 = 2}  for \code{\link{gaitdnbinomial}},
  else the default value should be used.

}
}

\details{
  This function is intended to make the 
  Generally-Truncated-Expansion (GTE) method
  easier for the user.
  It only makes sense if the linear predictors(s) are
  log  of the mean of the parent distribution,
  which is the usual case for
  \code{\link{gaitdpoisson}} and
  \code{\link{gaitdnbinomial}}.
  However, for \code{\link{gaitdlog}}
  and  \code{\link{gaitdzeta}} one should be using
  \code{\link[VGAMextra]{logffMlink}} and
  \code{\link[VGAMextra]{zetaffMlink}}.



  Without this function, the user must do quite a lot
  of book-keeping to know which columns of the offset
  matrix is to be assigned \code{log(mux)}.
  This can be rather laborious.



  In the fictitional example below the response is underdispersed
  with respect to a Poisson distribution and
  doubling the response achieves approximate equidispersion.



}

\value{
  A matrix with \code{n} rows and the same number of
  columns that a GAITD regression would produce for
  its matrix of linear predictors.
  The matrix can be inputted into \code{\link{vglm}}
  by assigning the \code{offset} argument.


  

}


%\author{
%  T. W. Yee.



%}
\note{
  This function is still in a developmental stage.
  The order of the arguments might change, hence it's
  safest to invoke it with full specification.




}
\seealso{
  \code{\link{gaitdpoisson}},
  \code{\link{gaitdlog}},
  \code{\link{gaitdzeta}},
  \code{\link{gaitdnbinomial}},
  \code{\link{Trunc}},
  \code{\link[stats]{offset}}.



}
\examples{
i.mix <- c(5, 10, 15, 20); a.mlm <- 13; mymux <- 2
goffset(mymux, 10, i.mix = i.mix, a.mlm = a.mlm)
\dontrun{org1  <- with(gdata, range(y))  # Original range of the data
vglm(mymux * y ~ 1,
     offset = goffset(mymux, nrow(gdata), i.mix = i.mix, a.mlm = a.mlm),
     gaitdpoisson(a.mlm = mymux * a.mlm, i.mix = mymux * i.mix,
                  truncate = Trunc(org1, mymux)),
     data = gdata)
}}

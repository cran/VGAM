\name{UtilitiesVGAM}
\alias{UtilitiesVGAM}
\alias{param.names}
\alias{dimm}
\alias{interleave.VGAM}
\title{Utility Functions for the VGAM Package }
\description{
  A set of common utility functions used by
  \pkg{VGAM} family functions.

}
\usage{
param.names(string, S)
dimm(M, hbw = M)
interleave.VGAM(.M, M1, inverse = FALSE)
}
\arguments{
  \item{string}{
  Character.
  Name of the parameter.


  }
  \item{M, .M}{
  Numeric. The total number of linear/additive predictors, called
  \eqn{M}.
  By total, it is meant summed over the number of responses.
  Often, \eqn{M} is the total number of parameters to be estimated (but
  this is not the same as the number of regression coefficients, unless
  the RHS of the formula is an intercept-only).
  The use of \code{.M} is unfortunate, but it is a compromise solution
  to what is presented in Yee (2015).
  Ideally, \code{.M} should be just \code{M}.


  }
  \item{M1}{
  Numeric. The number of linear/additive predictors for one response, called
  \eqn{M_1}.
  This argument used to be called \code{M}, but is now renamed properly.


  }
  \item{inverse}{
  Logical. Useful for the inverse function of \code{interleave.VGAM()}.



  }
  \item{S}{
  Numeric. The number of responses.


  }
  \item{hbw}{
  Numeric. The half-bandwidth, which measures the number
  of bands emanating from the central diagonal band.


  }
}
\value{
  For \code{param.names()}, this function returns the parameter names
  for \eqn{S} responses,
  i.e., \code{string} is returned unchanged if \eqn{S=1},
  else \code{paste(string, 1:S, sep = "")}.


  For \code{dimm()}, this function returns the number of elements
  to be stored for each of the working weight matrices.
  They are represented as columns in the matrix \code{wz} in
  e.g., \code{vglm.fit()}.
  See  the \emph{matrix-band} format described in
  Section 18.3.5 of Yee (2015).



  For \code{interleave.VGAM()}, this function returns a reordering
  of the linear/additive predictors depending on the number of responses.
  The arguments presented in Table 18.5 may not be valid
  in your version of Yee (2015).


}
%\section{Warning }{
%  The \code{zero} argument is supplied for convenience but conflicts
%}

\details{
  See Yee (2015) for some details about some of these functions.



}

\references{



Yee, T. W. (2015)
Vector Generalized Linear and Additive Models:
With an Implementation in R.
New York, USA: \emph{Springer}.



}

\seealso{
  \code{\link{CommonVGAMffArguments}},
  \code{\link{VGAM-package}}.


}
\author{T. W. Yee.
  Victor Miranda added the \code{inverse} argument to \code{interleave.VGAM()}.


}

%\note{
%  See \code{\link{Links}} regarding a major change in
%
%}

\examples{
param.names("shape", 1)  # "shape"
param.names("shape", 3)  # c("shape1", "shape2", "shape3")

dimm(3, hbw = 1)  # Diagonal matrix; the 3 elements need storage.
dimm(3)  # A general 3 x 3 symmetrix matrix has 6 unique elements.
dimm(3, hbw = 2)  # Tridiagonal matrix; the 3-3 element is 0 and unneeded.

M1 <- 2; ncoly <- 3; M <- ncoly * M1
mynames1 <- param.names("location", ncoly)
mynames2 <- param.names("scale",    ncoly)
(parameters.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)])
# The  following is/was in Yee (2015) and has a poor/deceptive style:
(parameters.names <- c(mynames1, mynames2)[interleave.VGAM(M, M  = M1)])
parameters.names[interleave.VGAM(M, M1 = M1, inverse = TRUE)]
}
\keyword{distribution}
\keyword{regression}
\keyword{programming}
\keyword{models}


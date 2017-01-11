\name{coefvgam}
\alias{coefvgam}
\alias{coef,vgam-method}
\alias{coefficients,vgam-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract Model Coefficients of a vgam() Object}
\description{
  Extracts the estimated
  coefficients from vgam() objects.


}
\usage{
coefvgam(object, type = c("linear", "nonlinear"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A
    \code{\link{vgam}} object.


  }
  \item{type}{ Character.
    The default is the first choice.


  }
  \item{\ldots}{
    Optional arguments fed into
    \code{\link{coefvlm}}.


  }
}
\details{
  For VGAMs, because modified backfitting is performed,
  each fitted function is decomposed into a linear and nonlinear
  (smooth) part.
  The argument \code{type} is used to return which one is wanted.



}
\value{
  A vector if \code{type = "linear"}.
  A list if \code{type = "nonlinear"}, and each component of
  this list corresponds to an \code{\link{s}} term;
  the component contains an S4 object with slot names such as
  \code{"Bcoefficients"},
  \code{"knots"},
  \code{"xmin"},
  \code{"xmax"}.


}
%\references{
%
%
%}
\author{ Thomas W. Yee }

%\note{
%}

%\section{Warning }{

%}

\seealso{
   \code{\link{vgam}},
   \code{\link{coefvlm}},
   \code{\link[stats]{coef}}.


%  \code{\link{coef-method}},


}
\examples{
fit <- vgam(agaaus ~ s(altitude, df = 2), binomialff, data = hunua)
coef(fit)  # Same as coef(fit, type = "linear")
(ii <- coef(fit, type = "nonlinear"))
is.list(ii)
names(ii)
slotNames(ii[[1]])
}
\keyword{models}
\keyword{regression}


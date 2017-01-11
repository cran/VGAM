\name{Coef.rrvglm-class}
\docType{class}
\alias{Coef.rrvglm-class}
\title{Class ``Coef.rrvglm'' }
\description{  The most pertinent matrices and other quantities
pertaining to a RR-VGLM. }
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{Coef(object, ...)} where \code{object} is an object
of class \code{rrvglm} (see \code{\link{rrvglm-class}}).

In this document, \eqn{M} is the number of linear predictors
and \eqn{n} is the number of observations.
}
\section{Slots}{
  \describe{
    \item{\code{A}:}{Of class \code{"matrix"}, \bold{A}. }
    \item{\code{B1}:}{Of class \code{"matrix"}, \bold{B1}. }
    \item{\code{C}:}{Of class \code{"matrix"}, \bold{C}. }
    \item{\code{Rank}:}{The rank of the RR-VGLM. }
    \item{\code{colx1.index}:}{Index of the columns of the
          \code{"vlm"}-type model matrix corresponding to the variables
          in \bold{x1}. These correspond to \bold{B1}.
    }
    \item{\code{colx2.index}:}{
          Index of the columns of the
          \code{"vlm"}-type model matrix corresponding to the variables
          in \bold{x2}. These correspond to the reduced-rank regression.
    }
    \item{\code{Atilde}:}{Object of class \code{"matrix"}, the
          \bold{A} matrix with the corner rows removed. Thus each of the
          elements have been estimated. This matrix is returned only
          if corner constraints were used.
    }
}
}

%\section{Methods}{
%No methods defined with class "Coef.rrvglm" in the signature.
%}

\references{


Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


}
\author{ Thomas W. Yee }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{Coef.rrvglm}},
  \code{\link{rrvglm}},
  \code{\link{rrvglm-class}},
  \code{print.Coef.rrvglm}.


}

\examples{
# Rank-1 stereotype model of Anderson (1984)
pneumo <- transform(pneumo, let = log(exposure.time), x3 = runif(nrow(pneumo)))
fit <- rrvglm(cbind(normal, mild, severe) ~ let + x3, multinomial, data = pneumo)
coef(fit, matrix = TRUE)
Coef(fit)
# print(Coef(fit), digits = 3)
}
\keyword{classes}

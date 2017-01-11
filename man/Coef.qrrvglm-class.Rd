\name{Coef.qrrvglm-class}
\docType{class}
\alias{Coef.qrrvglm-class}
\title{Class ``Coef.qrrvglm'' }
\description{
   The most pertinent matrices and other quantities pertaining to a
   QRR-VGLM (CQO model).

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{Coef(object,
...)} where \code{object} is an object of class \code{"qrrvglm"}
(created by \code{\link{cqo}}).

In this document, \eqn{R} is the \emph{rank}, \eqn{M} is the number of
linear predictors and \eqn{n} is the number of observations.

}
\section{Slots}{
  \describe{
    \item{\code{A}:}{Of class \code{"matrix"}, \bold{A}, which are the
    linear `coefficients' of the matrix of latent variables.
    It is \eqn{M} by \eqn{R}. }
    \item{\code{B1}:}{Of class \code{"matrix"}, \bold{B1}.
    These correspond to terms of the argument \code{noRRR}. }
    \item{\code{C}:}{Of class \code{"matrix"}, \bold{C}, the
    canonical coefficients. It has \eqn{R} columns. }
    \item{\code{Constrained}:}{Logical. Whether the model is
    a constrained ordination model. }
    \item{\code{D}:}{Of class \code{"array"},
    \code{D[,,j]} is an order-\code{Rank} matrix, for
    \code{j} = 1,\dots,\eqn{M}.
    Ideally, these are negative-definite in order to make the response
    curves/surfaces bell-shaped.
    }
    \item{\code{Rank}:}{The rank (dimension, number of latent variables)
    of the RR-VGLM. Called \eqn{R}. }
    \item{\code{latvar}:}{\eqn{n} by \eqn{R} matrix
          of latent variable values. }
    \item{\code{latvar.order}:}{Of class \code{"matrix"}, the permutation
          returned when the function
          \code{\link{order}} is applied to each column of \code{latvar}.
          This enables each column of \code{latvar} to be easily sorted.
          }
    \item{\code{Maximum}:}{Of class \code{"numeric"}, the
          \eqn{M} maximum fitted values. That is, the fitted values
          at the optimums for \code{noRRR = ~ 1} models.
    If \code{noRRR} is not \code{~ 1} then these will be \code{NA}s. }
    \item{\code{NOS}:}{Number of species.}
    \item{\code{Optimum}:}{Of class \code{"matrix"}, the values
          of the latent variables where the optimums are.
          If the curves are not bell-shaped, then the value will
          be \code{NA} or \code{NaN}.}
    \item{\code{Optimum.order}:}{Of class \code{"matrix"}, the permutation
          returned when the function
          \code{\link{order}} is applied to each column of \code{Optimum}.
          This enables each row of \code{Optimum} to be easily sorted.
          }
%   \item{\code{Diagonal}:}{Vector of logicals: are the
%         \code{D[,,j]} diagonal? }
    \item{\code{bellshaped}:}{Vector of logicals: is each
          response curve/surface bell-shaped? }
    \item{\code{dispersion}:}{Dispersion parameter(s). }
    \item{\code{Dzero}:}{Vector of logicals, is each of the
          response curves linear in the latent variable(s)?
          It will be if and only if
          \code{D[,,j]} equals \bold{O}, for
          \code{j} = 1,\dots,\eqn{M} . }
    \item{\code{Tolerance}:}{Object of class \code{"array"},
          \code{Tolerance[,,j]} is an order-\code{Rank} matrix, for
          \code{j} = 1,\dots,\eqn{M}, being the matrix of
          tolerances (squared if on the diagonal).
          These are denoted by \bold{T} in Yee (2004).
    Ideally, these are positive-definite in order to make the response
    curves/surfaces bell-shaped.
The tolerance matrices satisfy
\eqn{T_s = -\frac12 D_s^{-1}}{T_s = -(0.5 D_s^(-1)}.

     }
  }
}

%\section{Methods}{
%No methods defined with class "Coef.qrrvglm" in the signature.
%}
\references{


Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


}
\author{ Thomas W. Yee }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{Coef.qrrvglm}},
    \code{\link{cqo}},
    \code{print.Coef.qrrvglm}.


% \code{qrrvglm-class},


}

\examples{
x2 <- rnorm(n <- 100)
x3 <- rnorm(n)
x4 <- rnorm(n)
latvar1 <- 0 + x3 - 2*x4
lambda1 <- exp(3 - 0.5 * ( latvar1-0)^2)
lambda2 <- exp(2 - 0.5 * ( latvar1-1)^2)
lambda3 <- exp(2 - 0.5 * ((latvar1+4)/2)^2)
y1 <- rpois(n, lambda1)
y2 <- rpois(n, lambda2)
y3 <- rpois(n, lambda3)
yy <- cbind(y1, y2, y3)
# vvv p1 <- cqo(yy ~ x2 + x3 + x4, fam = poissonff, trace = FALSE)
\dontrun{
lvplot(p1, y = TRUE, lcol = 1:3, pch = 1:3, pcol = 1:3)
}
# vvv print(Coef(p1), digits = 3)
}
\keyword{classes}

\name{Coef.qrrvglm}
\alias{Coef.qrrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns Important Matrices etc. of a QO Object }
\description{
  This methods function returns important matrices etc. of a
  QO object.


}
\usage{
Coef.qrrvglm(object, varI.latvar = FALSE, refResponse = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
% A CQO or UQO object.
  A CQO object.
  The former has class \code{"qrrvglm"}.


  }


  \item{varI.latvar}{
  Logical indicating whether to scale the site scores (latent variables)
  to have variance-covariance matrix equal to the rank-\eqn{R} identity
  matrix. All models have uncorrelated site scores (latent variables),
  and this option stretches or shrinks the ordination axes if \code{TRUE}.
  See below for further details.



  }
  \item{refResponse}{
    Integer or character.
    Specifies the \emph{reference response} or \emph{reference species}.
    By default, the reference
    species is found by searching sequentially starting from the first
    species until a positive-definite tolerance matrix is found. Then
    this tolerance matrix is transformed to the identity matrix. Then
    the sites scores (latent variables) are made uncorrelated.
    See below for further details.


%   If \code{eq.tolerances=FALSE}, then transformations occur so that
%   the reference species has a tolerance matrix equal to the rank-\eqn{R}
%   identity matrix.


  }
  \item{\dots}{ Currently unused. }
}
\details{

  If \code{I.tolerances=TRUE} or \code{eq.tolerances=TRUE} (and its
  estimated tolerance matrix is positive-definite) then all species'
  tolerances are unity by transformation or by definition, and the spread
  of the site scores can be compared to them. Vice versa, if one wishes
  to compare the tolerances with the sites score variability then setting
  \code{varI.latvar=TRUE} is more appropriate.


  For rank-2 QRR-VGLMs, one of the species can be chosen so that the
  angle of its major axis and minor axis is zero, i.e., parallel to
  the ordination axes.  This means the effect on the latent vars is
  independent on that species, and that its tolerance matrix is diagonal.
  The argument \code{refResponse} allows one to choose which is the reference
  species, which must have a positive-definite tolerance matrix, i.e.,
  is bell-shaped.  If \code{refResponse} is not specified, then the code will
  try to choose some reference species starting from the first species.
  Although the \code{refResponse} argument could possibly be offered as
  an option when fitting the model, it is currently available after
  fitting the model, e.g., in the functions \code{\link{Coef.qrrvglm}} and
  \code{\link{lvplot.qrrvglm}}.


}
\value{
  The \bold{A}, \bold{B1}, \bold{C},  \bold{T},  \bold{D} matrices/arrays
  are returned, along with other slots.
  The returned object has class \code{"Coef.qrrvglm"}
  (see \code{\link{Coef.qrrvglm-class}}).


% For UQO, \bold{C} is undefined.



}
\references{
Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.


}
\author{ Thomas W. Yee }
\note{
Consider an equal-tolerances Poisson/binomial CQO model with \code{noRRR = ~ 1}.
For \eqn{R=1} it has about \eqn{2S+p_2}{2*S+p2} parameters.
For \eqn{R=2} it has about \eqn{3S+2 p_2}{3*S+2*p_2} parameters.
Here, \eqn{S} is the number of species, and \eqn{p_2=p-1}{p2=p-1} is
the number of environmental variables making up the latent variable.
For an unequal-tolerances Poisson/binomial CQO model with
\code{noRRR = ~ 1}, it has about \eqn{3S -1 +p_2}{3*S-1+p2} parameters
for \eqn{R=1}, and about \eqn{6S -3 +2p_2}{6*S -3 +2*p2} parameters
for \eqn{R=2}.
Since the total number of data points is \eqn{nS}{n*S}, where
\eqn{n} is the number of sites, it pays to divide the number
of data points by the number of parameters to get some idea
about how much information the parameters contain.


}

% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{cqo}},
\code{\link{Coef.qrrvglm-class}},
\code{print.Coef.qrrvglm},
\code{\link{lvplot.qrrvglm}}.


}

\examples{
set.seed(123)
x2 <- rnorm(n <- 100)
x3 <- rnorm(n)
x4 <- rnorm(n)
latvar1 <- 0 + x3 - 2*x4
lambda1 <- exp(3 - 0.5 * ( latvar1-0)^2)
lambda2 <- exp(2 - 0.5 * ( latvar1-1)^2)
lambda3 <- exp(2 - 0.5 * ((latvar1+4)/2)^2)  # Unequal tolerances
y1 <- rpois(n, lambda1)
y2 <- rpois(n, lambda2)
y3 <- rpois(n, lambda3)
set.seed(111)
# vvv p1 <- cqo(cbind(y1, y2, y3) ~ x2 + x3 + x4, poissonff, trace = FALSE)
\dontrun{ lvplot(p1, y = TRUE, lcol = 1:3, pch = 1:3, pcol = 1:3)
}
# vvv Coef(p1)
# vvv print(Coef(p1), digits=3)
}
\keyword{models}
\keyword{regression}



\name{Hurea}
\alias{Hurea}
\alias{dhurea}
%\alias{phurea}
%\alias{qhurea}
%\alias{rhurea}
\title{The Husler-Reiss Angular Surface Distribution}
\description{
  Density
  for the Husler-Reiss angular surface distribution.

% distribution function,
% quantile function and random generation

}
\usage{
dhurea(x, shape, log = FALSE)
}
\arguments{
  \item{x}{
  Same as \code{\link[stats:Uniform]{Uniform}}.


  }
  \item{shape}{the positive (shape) parameter.
  It is often called \eqn{lambda}{\lambda} and
  it might not be a shape parameter at all.


  }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the
  density is returned.

  }
% \item{lower.tail, log.p}{
% Same meaning as in \code{\link[stats:Normal]{pnorm}}
% or \code{\link[stats:Normal]{qnorm}}.


% }

}
\value{
  \code{dhurea} gives the density.

% \code{phurea} gives the distribution function,
% \code{qhurea} gives the quantile function, and
% \code{rhurea} generates random deviates.


}
%\references{

%  Mhalla, L. and de Carvalho, M. and Chavez-Demoulin,
%  V. (2019).
%  Regression-type models for extremal dependence.
%  \emph{Scandinavian Journal of Statistics},
%  \bold{46}, 1141--1167.



%}
\author{ T. W. Yee }
\details{
  See \code{\link{hurea}}, the \pkg{VGAM}
  family function for
  estimating the (shape) parameter \eqn{s} by
  maximum likelihood
  estimation, for the formula of the
  probability density function.


}
\note{
  Difficulties are encountered as
  the shape parameter approaches 0 with
  respect to \code{\link[stats]{integrate}}
  because the density converges to a degenerate
  distribution with probability mass at 0 and 1.
  That is, when \eqn{s} is around 0.5 the
  density is ``u'' shaped and the area around the
  endpoints becomes concentrated at the
  two points.
  See the examples below.
  Approximately, the 
  density is ``u'' shaped for \eqn{s < 1}
  and unimodal for \eqn{s > 2}.

  

}
\section{Warning}{
  The cases \code{x == 0}, \code{x == 1},
  \code{shape == 0} and
  \code{shape == Inf}
  may not be handled correctly.



}
\seealso{
  \code{\link{hurea}}.



}
\examples{
integrate(dhurea, 0, 1, shape = 0.20)  # Incorrect
integrate(dhurea, 0, 1, shape = 0.35)  # struggling but okay
\dontrun{x <- seq(0, 1, length = 501)
par(mfrow = c(2, 2))
plot(x, dhurea(x, 0.7), col = "blue", type = "l")
plot(x, dhurea(x, 1.1), col = "blue", type = "l")
plot(x, dhurea(x, 1.4), col = "blue", type = "l")
plot(x, dhurea(x, 3.0), col = "blue", type = "l")
}
}
\keyword{distribution}



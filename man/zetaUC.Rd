\name{Zeta}
\alias{Zeta}
\alias{dzeta}
\alias{pzeta}
\alias{qzeta}
\alias{rzeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Zeta Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the zeta distribution.

}
\usage{
dzeta(x, shape, log = FALSE)
pzeta(q, shape, lower.tail = TRUE)
qzeta(p, shape)
rzeta(n, shape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{Same as \code{\link[stats]{Poisson}}. }
  \item{shape}{ The positive shape parameter \eqn{p}. }
  \item{lower.tail, log}{
  Same meaning as in \code{\link[stats]{Normal}}.


  }
}
\details{
  The density function of the zeta distribution is given by
  \deqn{y^{-s-1} / \zeta(s+1)}{%
        y^(-s-1) /  zeta(s+1)}
  where \eqn{s>0}, \eqn{y=1,2,\ldots}, and \eqn{\zeta}{zeta} is
  Riemann's zeta function.



}
\value{
  \code{dzeta} gives the density,
  \code{pzeta} gives the distribution function,
  \code{qzeta} gives the quantile function, and
  \code{rzeta} generates random deviates.


}
\references{


Johnson N. L., Kotz S., and Balakrishnan N. (1993)
\emph{Univariate Discrete Distributions},
2nd ed.
New York: Wiley.



% Lindsey, J. K. (2002zz)
% \emph{Applied Statistical Modelling}, 2nd ed.
% London: Chapman & Hall.zz


% Knight, K. (2002zz)
% Theory book.
% London: Chapman & Hall.zz

}
\author{ T. W. Yee }
\note{
    \code{qzeta()} runs slower and slower as \code{shape} approaches 0
    and \code{p} approaches 1.
    The \pkg{VGAM} family function \code{\link{zetaff}} estimates the
    shape parameter \eqn{s}.



}

%\section{Warning}{
%    These functions have not been fully tested.



%}

\seealso{
    \code{\link{zeta}},
    \code{\link{zetaff}},
    \code{\link{Oazeta}},
    \code{\link{Oizeta}},
    \code{\link{Otzeta}}.


}

\examples{
dzeta(1:20, shape = 2)
myshape <- 0.5
max(abs(pzeta(1:200, myshape) -
    cumsum(1/(1:200)^(1+myshape)) / zeta(myshape+1)))  # Should be 0

\dontrun{ plot(1:6, dzeta(1:6, 2), type = "h", las = 1,
               col = "orange", ylab = "Probability",
     main = "zeta probability function; orange: shape = 2; blue: shape = 1")
points(0.10 + 1:6, dzeta(1:6, 1), type = "h", col = "blue") }
}
\keyword{distribution}


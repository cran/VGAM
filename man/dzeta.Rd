\name{Zeta}
\alias{Zeta}
\alias{dzeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Zeta Distribution }
\description{
  Density for the zeta distribution.
}
% zz p is not a good argument name, esp. with qzeta(p, p)
\usage{
dzeta(x, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Numerical vector/matrix to evaluate the density. }
  \item{p}{ The parameter \eqn{p}. This must be greater than 1. }
}
\details{
  The density function of the zeta distribution is given by
  \deqn{y^{-p} / \zeta(p)}{%
        y^(-p) /  zeta(p)}
  where \eqn{p>1}, \eqn{y=1,2,\ldots}, and \eqn{\zeta}{zeta} is
  Riemann's zeta function.

}
\value{
  Returns the density evaluated at \code{x}.
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
    The \pkg{VGAM} family function \code{\link{zetaff}} estimates the
    parameter \eqn{p}.

}

\section{Warning}{
    This function has not been fully tested.
}

\seealso{
    \code{\link{zeta}},
    \code{\link{zetaff}}.
}

\examples{
dzeta(1:20, p=2)
\dontrun{
plot(1:6, dzeta(1:6, p=4), type="h", las=1, ylab="Probability",
     main="zeta probability function; black: p = 4; blue: p = 2")
points(0.10+1:6, dzeta(1:6, p=2), type="h", col="blue")
}
}
\keyword{distribution}


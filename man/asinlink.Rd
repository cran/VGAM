\name{asinlink}
\alias{asinlink}
%- Also NEED an '\alias' for EACH other
% topic documented here.
\title{ Arcsine Link Function}
\description{
  Computes the arcsine link,
  including its inverse and
  the first few derivatives.

}
\usage{
asinlink(theta, bvalue = NULL, inverse = FALSE,
   deriv = 0, short = TRUE, tag = FALSE, c10 = c(4, -pi))
}
%- maybe also 'usage' for other objs doced here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.
  }
  \item{bvalue}{
  See \code{\link{Links}}.


  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }
  \item{c10}{
    Similar to \code{\link{sqrtlink}}.
    The default is intended to match
    \code{\link{lcalogitlink}} for \code{\link{binomialff}}
    at binomial probabilities (\code{theta}) equal
    to 0.5.
    
  }

}
\details{
  Function \code{\link{alogitlink}}
  gives some motivation for this link.
  However, the problem with this link
  is that it is bounded by default
  between \code{(-pi, pi)}
  so that it can be unsuitable for regression.
  This link is a scaled and centred
  CDF of the arcsine distribution.
  The centring is chosen so that
  \code{asinlink(0.5)} is 0,
  and the scaling is chosen so that
  \code{asinlink(0.5, deriv = 1)} and
  \code{logitlink(0.5, deriv = 1)}
  are equal (the value 4 actually),
  hence this link will operate similar to the
  \code{\link{logitlink}}
  when close to 0.5.



}
\value{
  Similar to \code{\link{logitlink}}
  but using different formulas.



}
%\references{

% Yee, T. W. (2023).
% \emph{Constant information augmented link
% functions impervious
% to the Hauck--Donner effect in
% vector generalized linear models}.
% Under review.



%}
\author{ Thomas W. Yee }

\section{Warning }{
  It is possible that the scaling might change
  in the future.

}

%\note{
%  Numerical instability may occur when
%  \code{theta} is close to 1 or 0.  One way of
%  overcoming this is to use \code{bvalue}.

%}

\seealso{
  \code{\link{logitlink}},
  \code{\link{alogitlink}},
  \code{\link{Links}},
  \code{\link{probitlink}},
  \code{\link{clogloglink}},
  \code{\link{cauchitlink}},
  \code{\link{binomialff}},
  \code{\link{sloglink}},
  \code{\link{hdeff}}.



 }
\examples{
p <- seq(0.01, 0.99, length= 10)
asinlink(p)
max(abs(asinlink(asinlink(p), inv = TRUE) - p))  # 0?

\dontrun{
par(mfrow = c(2, 2), lwd = (mylwd <- 2))
y <- seq(-4, 4, length = 100)
p <- seq(0.01, 0.99, by = 0.01)

for (d in 0:1) {
  matplot(p, cbind(logitlink(p, deriv = d), probitlink(p, deriv = d)),
          type = "n", col = "blue", ylab = "transformation",
          log = ifelse(d == 1, "y", ""),
          las = 1, main = if (d == 0) "Some probability link functions"
          else "First derivative")
  lines(p,   logitlink(p, deriv = d), col = "green")
  lines(p,  probitlink(p, deriv = d), col = "blue")
  lines(p, clogloglink(p, deriv = d), col = "tan")
  lines(p,    asinlink(p, deriv = d), col = "red3")
  if (d == 0) {
    abline(v = 0.5, h = 0, lty = "dashed")
    legend(0, 4.5, c("logitlink", "probitlink", "clogloglink",
           "asinlink"), lwd = mylwd,
           col = c("green", "blue", "tan", "red3"))
  } else
    abline(v = 0.5, lwd = 0.5, col = "gray")
}

for (d in 0) {
  matplot(y, cbind( logitlink(y, deriv = d, inverse = TRUE),
                   probitlink(y, deriv = d, inverse = TRUE)),
          type  = "n", col = "blue", xlab = "transformation", ylab = "p",
          main = if (d == 0) "Some inverse probability link functions"
          else "First derivative", las=1)
  lines(y,   logitlink(y, deriv = d, inverse = TRUE), col = "green")
  lines(y,  probitlink(y, deriv = d, inverse = TRUE), col = "blue")
  lines(y, clogloglink(y, deriv = d, inverse = TRUE), col = "tan")
  lines(y,    asinlink(y, deriv = d, inverse = TRUE), col = "red3")
  if (d == 0) {
      abline(h = 0.5, v = 0, lwd = 0.5, col = "gray")
      legend(-4, 1, c("logitlink", "probitlink", "clogloglink",
             "asinlink"), lwd = mylwd,
             col = c("green", "blue", "tan", "red3"))
  }
}
par(lwd = 1)
}
}
\keyword{math}
\keyword{models}
\keyword{regression}

%plot(y, logitlink(y, inverse = TRUE), type = "l", col = "green",
%     xlab = "transformation", ylab = "p",
%     lwd=2, las=1, main = "Some inverse probability link functions")
%lines(y, probitlink(y, inverse = TRUE), col = "blue", lwd=2)
%lines(y, clogloglink(y, inverse = TRUE), col = "tan", lwd=2)
%abline(h=0.5, v = 0, lty = "dashed")





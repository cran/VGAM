\name{sqrtlink}
%\name{foldsqrtlink}   % orig. prior to 20231023
\alias{foldsqrtlink}
\alias{sqrtlink}
%\alias{sqrtpoislink}
% \alias{foldsqrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Square Root and Folded Square Root Link Functions }
\description{
  Computes the
  square root and
  folded square root
  transformations,
  including their inverse
  and their first two derivatives.

}
\usage{
foldsqrtlink(theta, min = 0, max = 1, mux = sqrt(2),
     inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
sqrtlink(theta, inverse = FALSE, deriv = 0, short = TRUE,
         tag = FALSE, c10 = c(2, -2))
}
%- maybe also 'usage' for other objects documented here.
%sqrtpoislink(theta, inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{min, max, mux}{
  These are called \eqn{L}, \eqn{U} and \eqn{K} below.

  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }
  \item{c10}{
  Numeric, 2-vector \code{c(c1, c0)} 
  for a linear transformation.
  The plain link is multiplied by \code{c1}
  and then \code{c0} is added so that \code{c1 = 1:0} is
  simply \code{\link[base:sqrt]{sqrt}}.
  The default is intended to match
  \code{\link{lcsloglink}} for \code{\link{poissonff}}
  at \code{lambda} (\code{theta}) equal to 1.



  }

}
\details{
  The folded square root link function can be applied to
  parameters that lie between \eqn{L} and \eqn{U}
  inclusive.
  Numerical values of \code{theta}
  out of range result in \code{NA} or \code{NaN}.



% While the ordinary square root link
% is \code{\link{sqrtlink}}
% (and is a special case of
% \code{\link{foldsqrtlink}})
% the link \code{\link{sqrtpoislink}}
% is a scaled version of \code{\link{sqrtlink}}
% specifically for Poisson regression.
% It is multiplied by two and has two subtracted
% from it.
% In particular, it is a
% constant information augmented (CIA) link function
% that allows a form of Poisson regression to be performed
% that combats the Hauck--Donner effect (HDE).
% That is, the resulting Wald statistic p-values are
% HDE-free.
% The scaled version matches the 0th and 1st
% derivatives of \code{\link{loglink}} at
% the null value of unity.
% Because \code{\link{sqrtpoislink}}
% has a restricted range,
% it can be unsuitable for regression,
% hence it is augmented with
% \code{\link{loglink}} in the form of a mixture,
% and this is implemented in \code{\link{sloglink}}.
  

  
  More general information
  can be found at \code{\link{alogitlink}}.


}
\value{
  For \code{foldsqrtlink} with \code{deriv = 0}:
  \eqn{K (\sqrt{\theta-L} - \sqrt{U-\theta})}{K *
          (sqrt(theta-L) - sqrt(U-theta))}
  or
  \code{mux * (sqrt(theta-min) - sqrt(max-theta))}
  when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then some more
  complicated function that returns a \code{NA} unless
  \code{theta} is between \code{-mux*sqrt(max-min)} and
  \code{mux*sqrt(max-min)}.


  For \code{sqrtlink} with \code{deriv = 0}
  and \code{c10 = 1:0}:
  \eqn{\sqrt{\theta}}{sqrt(theta)}
  when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then the square
  is returned.



  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a
  function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it
  returns the reciprocal.


}
%\references{
%
%}
\author{ Thomas W. Yee }

\note{
  For \code{foldsqrtlink},
  the default has, if \code{theta} is 0 or 1,
  the link function
  value is \code{-sqrt(2)} and \code{+sqrt(2)}
  respectively.
  These are finite values, therefore one cannot use
  this link function for
  general modelling of probabilities because
  of numerical problem,
  e.g., with \code{\link{binomialff}},
  \code{\link{cumulative}}. See
  the example below.


}

\seealso{
  \code{\link{Links}},
  \code{\link{poissonff}},
  \code{\link{sloglink}},
  \code{\link{hdeff}}.



% \url{https://www.cia.gov}.


  
}
\examples{
p <- seq(0.01, 0.99, by = 0.01)
foldsqrtlink(p)
max(abs(foldsqrtlink(foldsqrtlink(p), inverse = TRUE) - p))  # 0

p <- c(seq(-0.02, 0.02, by = 0.01), seq(0.97, 1.02, by = 0.01))
foldsqrtlink(p)  # Has NAs

\dontrun{
p <- seq(0.01, 0.99, by = 0.01)
par(mfrow = c(2, 2), lwd = (mylwd <- 2))
y <- seq(-4, 4, length = 100)
for (d in 0:1) {
  matplot(p, cbind(   logitlink(p, deriv = d),
                   foldsqrtlink(p, deriv = d)),
          col = "blue", ylab = "transformation",
          main = ifelse(d == 0, "Some probability links",
          "First derivative"), type = "n", las = 1)
  lines(p,    logitlink(p, deriv = d), col = "green")
  lines(p,   probitlink(p, deriv = d), col = "blue")
  lines(p,  clogloglink(p, deriv = d), col = "red")
  lines(p, foldsqrtlink(p, deriv = d), col = "tan")
  if (d == 0) {
    abline(v = 0.5, h = 0, lty = "dashed")
    legend(0, 4.5, c("logitlink", "probitlink",
                     "clogloglink", "foldsqrtlink"),
           lwd = 2, col = c("green", "blue",
                            "red", "tan"))
  } else
    abline(v = 0.5, lty = "dashed")
}

for (d in 0) {
  matplot(y,
          cbind(   logitlink(y, deriv = d, inverse = TRUE),
                foldsqrtlink(y, deriv = d, inverse = TRUE)),
          type = "n", col = "blue", xlab = "transformation",
          ylab = "p", lwd = 2, las = 1, main = if (d == 0)
          "Some inverse probability link functions" else
          "First derivative")
  lines(y,    logitlink(y, deriv=d, inverse=TRUE), col="green")
  lines(y,   probitlink(y, deriv=d, inverse=TRUE), col="blue")
  lines(y,  clogloglink(y, deriv=d, inverse=TRUE), col="red")
  lines(y, foldsqrtlink(y, deriv=d, inverse=TRUE), col="tan")
  if (d == 0) {
    abline(h = 0.5, v = 0, lty = "dashed")
    legend(-4, 1, c("logitlink", "probitlink",
                    "clogloglink", "foldsqrtlink"), lwd = 2, 
           col = c("green", "blue", "red", "tan"))
  }
}
par(lwd = 1)
}

# This is lucky to converge
fit.h <- vglm(agaaus ~ sm.bs(altitude),
              binomialff(foldsqrtlink(mux = 5)),
              hunua, trace = TRUE)
\dontrun{
plotvgam(fit.h, se = TRUE, lcol = "orange", scol = "orange",
         main = "Orange is Hunua, Blue is Waitakere") }
head(predict(fit.h, hunua, type = "response"))

\dontrun{
# The following fails.
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let,
       cumulative(foldsqrtlink(mux = 10), par = TRUE, rev = TRUE),
       data = pneumo, trace = TRUE, maxit = 200) }
}
\keyword{math}
\keyword{models}
\keyword{regression}


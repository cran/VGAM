\name{A1A2A3}
\alias{A1A2A3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The A1A2A3 Blood Group System }
\description{
  Estimates the three independent parameters of the
  the A1A2A3 blood group system.

}
\usage{
A1A2A3(link = "logit", inbreeding = FALSE, ip1 = NULL, ip2 = NULL, iF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to \code{p1}, \code{p2} and \code{f}.
  See \code{\link{Links}} for more choices.

  }
  \item{inbreeding}{
    Logical. Is there inbreeding?


%   Logical. Is the HWE assumption to be made?


  }
  \item{ip1, ip2, iF}{
  Optional initial value for \code{p1}, \code{p2} and \code{f}.

  }
}
\details{
  The parameters \code{p1} and \code{p2} are probabilities, so that
  \code{p3=1-p1-p2} is the third probability.
  The parameter \code{f} is the third independent parameter if
  \code{inbreeding = TRUE}.
  If \code{inbreeding = FALSE} then \eqn{f = 0} and Hardy-Weinberg
  Equilibrium (HWE) is assumed.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{

Lange, K. (2002)
\emph{Mathematical and Statistical Methods for Genetic Analysis},
2nd ed. New York: Springer-Verlag.


}
\author{ T. W. Yee }
\note{
  The input can be a 6-column matrix of counts,
  with columns corresponding to
  \code{A1A1},
  \code{A1A2},
  \code{A1A3},
  \code{A2A2},
  \code{A2A3},
  \code{A3A3} (in order).
  Alternatively, the input can be a 6-column matrix of
  proportions (so each row adds to 1) and the \code{weights}
  argument is used to specify the total number of counts for each row.


}

\seealso{
  \code{\link{AA.Aa.aa}},
  \code{\link{AB.Ab.aB.ab}},
  \code{\link{ABO}},
  \code{\link{MNSs}}.


% \code{\link{AB.Ab.aB.ab2}},



}
\examples{
ymat <- cbind(108, 196, 429, 143, 513, 559)
fit <- vglm(ymat ~ 1, A1A2A3(link = probit), trace = TRUE, crit = "coef")
fit <- vglm(ymat ~ 1, A1A2A3(link = logit, ip1 = 0.3, ip2 = 0.3, iF = 0.02),
            trace = TRUE, crit = "coef")
Coef(fit)  # Estimated p1 and p2
rbind(ymat, sum(ymat) * fitted(fit))
sqrt(diag(vcov(fit)))
}
\keyword{models}
\keyword{regression}

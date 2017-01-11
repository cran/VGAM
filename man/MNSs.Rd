\name{MNSs}
\alias{MNSs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The MNSs Blood Group System }
\description{
  Estimates the three independent parameters of the
  the MNSs blood group system.
}
\usage{
MNSs(link = "logit", imS = NULL, ims = NULL, inS = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to the three parameters.
  See \code{\link{Links}} for more choices.

  }
  \item{imS, ims, inS}{
  Optional initial value for \code{mS}, \code{ms}
  and \code{nS} respectively.
  A \code{NULL} means they are computed internally.

  }
}
\details{
  There are three independent
  parameters: \code{m_S}, \code{m_s}, \code{n_S}, say, so that
  \code{n_s = 1 - m_S - m_s - n_S}.
  We let the eta vector (transposed) be
  \code{(g(m_S), g(m_s), g(n_S))} where \code{g} is the
  link function.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{
  Elandt-Johnson, R. C. (1971)
  \emph{Probability Models and Statistical Methods in Genetics},
  New York: Wiley.


}
\author{ T. W. Yee }
\note{
  The input can be a 6-column matrix of counts, where the columns are
  MS, Ms, MNS, MNs, NS, Ns (in order).
  Alternatively, the input can be a 6-column matrix of
  proportions (so each row adds to 1) and the \code{weights}
  argument is used to specify the total number of counts for each row.


}
\seealso{
  \code{\link{AA.Aa.aa}},
  \code{\link{AB.Ab.aB.ab}},
  \code{\link{ABO}},
  \code{\link{A1A2A3}}.


% \code{\link{AB.Ab.aB.ab2}},


}
\examples{
# Order matters only:
y <- cbind(MS = 295, Ms = 107, MNS = 379, MNs = 322, NS = 102, Ns = 214)
fit <- vglm(y ~ 1, MNSs("logit", .25, .28, .08), trace = TRUE)
fit <- vglm(y ~ 1, MNSs(link = logit), trace = TRUE, crit = "coef")
Coef(fit)
rbind(y, sum(y)*fitted(fit))
sqrt(diag(vcov(fit)))
}
\keyword{models}
\keyword{regression}

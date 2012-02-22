\name{ABO}
\alias{ABO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The ABO Blood Group System }
\description{
  Estimates the two independent parameters of the 
  the ABO blood group system.

}
\usage{
ABO(link = "logit", earg=list(), ipA = NULL, ipO = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to \code{pA} and \code{pB}.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument applied to each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ipA, ipO}{
  Optional initial value for \code{pA} and \code{pO}.
  A \code{NULL} value means values are computed internally.
  }

}
\details{
  The parameters \code{pA} and \code{pB} are probabilities, so that
  \code{pO=1-pA-pB} is the third probability.
  The probabilities \code{pA} and \code{pB} correspond to A and B respectively,
  so that \code{pO} is the probability for O.
  It is easier to make use of initial values for \code{pO} than for \code{pB}.
  In documentation elsewhere I sometimes use
  \code{pA=p},
  \code{pB=q},
  \code{pO=r}.


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
  The input can be a 4-column matrix of counts, where the columns 
  are A, B, AB, O (in order).
  Alternatively, the input can be a 4-column matrix of 
  proportions (so each row adds to 1) and the \code{weights}
  argument is used to specify the total number of counts for each row.


}

\seealso{
  \code{\link{AA.Aa.aa}},
  \code{\link{AB.Ab.aB.ab}},
  \code{\link{AB.Ab.aB.ab2}},
  \code{\link{G1G2G3}},
  \code{\link{MNSs}}.


}
\examples{
ymatrix = cbind(A = 725, B = 258, AB = 72, O = 1073) # Order matters, not the name
fit = vglm(ymatrix ~ 1, ABO(link = identity), trace = TRUE, cri = "coef")
coef(fit, matrix = TRUE)
Coef(fit) # Estimated pA and pB
rbind(ymatrix, sum(ymatrix) * fitted(fit))
sqrt(diag(vcov(fit)))
}
\keyword{models}
\keyword{regression}

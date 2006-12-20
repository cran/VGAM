\name{ABO}
\alias{ABO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The ABO Blood Group System }
\description{
  Estimates the two independent parameters of the 
  the ABO blood group system.
}
\usage{
ABO(link = "logit", earg=list(), ir = NULL, ip = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to \code{p} and \code{q}.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument applied to each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ir, ip}{ Optional initial value for \code{r} and \code{p}.
  A \code{NULL} value means values are computed internally. }
}
\details{
  The parameters \code{p} and \code{q} are probabilities, so that
  \code{r=1-p-q} is the third probability.
  The probabilities \code{p} and \code{r} correspond to
  A and O respectively, so that \code{q} is the probability for B.
  It is easier to make use of initial values for \code{r} than for
  \code{q}.
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
y = cbind(A=725, B=258, AB=72, O=1073) # Order matters, not the name
fit = vglm(y ~ 1, ABO(link=logit), trace=TRUE, cri="coef")
fit = vglm(y ~ 1, ABO(link=identity), trace=TRUE, cri="coef")
coef(fit, matrix=TRUE)
Coef(fit) # Estimated p and q
rbind(y, sum(y)*fitted(fit))
sqrt(diag(vcov(fit)))
}
\keyword{models}
\keyword{regression}

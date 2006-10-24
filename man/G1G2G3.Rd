\name{G1G2G3}
\alias{G1G2G3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The G1G2G3 Blood Group System }
\description{
  Estimates the three independent parameters of the 
  the G1G2G3 blood group system.
}
\usage{
G1G2G3(link = "logit", ip1 = NULL, ip2 = NULL, iF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to \code{p1}, \code{p2} and \code{f}.
  See \code{\link{Links}} for more choices.

  }
  \item{ip1, ip2, iF}{
  Optional initial value for \code{p1}, \code{p2} and \code{f}.

  }
}
\details{
  The parameters \code{p1} and \code{p2} are probabilities, so that
  \code{p3=1-p1-p2} is the third probability.
  The parameter \code{f} is the third independent parameter.

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
\code{G_1/G_1},
\code{G_1/G_2},
\code{G_1/G_3},
\code{G_2/G_2},
\code{G_2/G_3},
\code{G_3/G_3} (in order). 
Alternatively, the input can be a 6-column matrix of 
proportions (so each row adds to 1) and the \code{weights}
argument is used to specify the total number of counts for each row.

}

\seealso{
\code{\link{AA.Aa.aa}},
\code{\link{AB.Ab.aB.ab}},
\code{\link{AB.Ab.aB.ab2}},
\code{\link{ABO}},
\code{\link{MNSs}}.
}
\examples{
y = cbind(108, 196, 429, 143, 513, 559)
fit = vglm(y ~ 1, G1G2G3(link=probit), trace=TRUE, crit="coef")
fit = vglm(y ~ 1, G1G2G3(link=logit, .3, .3, .02), trace=TRUE, crit="coef")
fit = vglm(y ~ 1, G1G2G3(link="identity"), trace=TRUE)
Coef(fit) # Estimated p1, p2 and f
rbind(y, sum(y)*fitted(fit))
sqrt(diag(vcov(fit)))
}
\keyword{models}
\keyword{regression}

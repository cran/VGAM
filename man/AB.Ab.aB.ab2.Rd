\name{AB.Ab.aB.ab2}
\alias{AB.Ab.aB.ab2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The AB-Ab-aB-ab2 Blood Group System }
\description{
  Estimates the parameter of the 
  the AB-Ab-aB-ab2 blood group system.

}
\usage{
AB.Ab.aB.ab2(link = "logit", earg=list(), init.p = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{ 
  Link function applied to \code{p}.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{init.p}{ Optional initial value for \code{p}. }
}
\details{
  This one parameter model involves a probability called \code{p}.

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
The input can be a 4-column matrix of counts.
Alternatively, the input can be a 4-column matrix of
proportions (so each row adds to 1) and the \code{weights}
argument is used to specify the total number of counts for each row.

}

\section{Warning}{
There may be a bug in the \code{deriv} and \code{weight} slot of the
family function.

}
\seealso{
\code{\link{AA.Aa.aa}},
\code{\link{AB.Ab.aB.ab}},
\code{\link{ABO}},
\code{\link{G1G2G3}},
\code{\link{MNSs}}.
}

\examples{
# See Elandt-Johnson, pp.430,427
# Estimated variance is approx 0.0021
y = cbind(68, 11, 13, 21)
fit = vglm(y ~ 1, AB.Ab.aB.ab2(link=logit), trace=TRUE, crit="coef")
fit = vglm(y ~ 1, AB.Ab.aB.ab2(link=cloglog), trace=TRUE, crit="coef")
Coef(fit)   # Estimated p
rbind(y, sum(y)*fitted(fit))
sqrt(diag(vcov(fit)))
}
\keyword{models}
\keyword{regression}

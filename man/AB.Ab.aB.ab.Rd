\name{AB.Ab.aB.ab}
\alias{AB.Ab.aB.ab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The AB-Ab-aB-ab Blood Group System }
\description{
   Estimates the parameter of the 
   AB-Ab-aB-ab blood group system.

}
\usage{
AB.Ab.aB.ab(link = "logit", earg=list(), init.p = NULL)
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
Lange, K. (2002)
\emph{Mathematical and Statistical Methods for Genetic Analysis},
2nd ed. New York: Springer-Verlag.

}
\author{ T. W. Yee }
\note{ 
The input can be a 4-column matrix of counts, where the columns 
are AB, Ab, aB and ab
(in order).
Alternatively, the input can be a 4-column matrix of 
proportions (so each row adds to 1) and the \code{weights}
argument is used to specify the total number of counts for each row.

}

\seealso{
\code{\link{AA.Aa.aa}},
\code{\link{AB.Ab.aB.ab2}},
\code{\link{ABO}},
\code{\link{G1G2G3}},
\code{\link{MNSs}}.

}

\examples{
y = cbind(AB=1997, Ab=906, aB=904, ab=32) # Data from Fisher (1925)
fit = vglm(y ~ 1, AB.Ab.aB.ab(link="identity", init.p=0.9), trace=TRUE)
fit = vglm(y ~ 1, AB.Ab.aB.ab, trace=TRUE)
rbind(y, sum(y)*fitted(fit))
Coef(fit) # Estimated p
p = sqrt(4*(fitted(fit)[,4]))
p*p
summary(fit)
}
\keyword{models}
\keyword{regression}

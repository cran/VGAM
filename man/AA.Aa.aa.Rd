\name{AA.Aa.aa}
\alias{AA.Aa.aa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The AA-Aa-aa Blood Group System }
\description{
   Estimates the parameter of the 
   AA-Aa-aa blood group system.
}
\usage{
AA.Aa.aa(link = "logit", earg=list(), init.pA = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to \code{pA}.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{init.pA}{ Optional initial value for \code{pA}. }
}
\details{
  This one parameter model involves a probability called \code{pA}.
  The probability of getting a count in the first column of the
  input (an AA) is \code{pA*pA}.
}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.
}
\references{ 
Weir, B. S. (1996)
\emph{Genetic Data Analysis II: Methods for Discrete Population
      Genetic Data},
Sunderland, MA: Sinauer Associates, Inc.
}
\author{ T. W. Yee }
\note{ 
  The input can be a 3-column matrix of counts, where the columns 
  are AA, Ab and aa
  (in order).
  Alternatively, the input can be a 3-column matrix of 
  proportions (so each row adds to 1) and the \code{weights}
  argument is used to specify the total number of counts for each row.


}
\seealso{
\code{\link{AB.Ab.aB.ab}},
\code{\link{AB.Ab.aB.ab2}},
\code{\link{ABO}},
\code{\link{G1G2G3}},
\code{\link{MNSs}}.

}
\examples{
y = cbind(53, 95, 38)
fit = vglm(y ~ 1, AA.Aa.aa(link="logit"), trace=TRUE)
fit = vglm(y ~ 1, AA.Aa.aa(link="probit"), trace=TRUE)
rbind(y, sum(y)*fitted(fit))
Coef(fit) # Estimated pA
summary(fit)
}
\keyword{models}
\keyword{regression}

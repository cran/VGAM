\name{AB.Ab.aB.ab}
\alias{AB.Ab.aB.ab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The AB-Ab-aB-ab Blood Group System }
\description{
   Estimates the parameter of the
   AB-Ab-aB-ab blood group system.

}
\usage{
AB.Ab.aB.ab(link = "logit", init.p = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to \code{p}.
  See \code{\link{Links}} for more choices.

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
  \code{\link{ABO}},
  \code{\link{A1A2A3}},
  \code{\link{MNSs}}.


% \code{\link{AB.Ab.aB.ab2}},



}

\examples{
ymat <- cbind(AB=1997, Ab=906, aB=904, ab=32)  # Data from Fisher (1925)
fit <- vglm(ymat ~ 1, AB.Ab.aB.ab(link = "identitylink"), trace = TRUE)
fit <- vglm(ymat ~ 1, AB.Ab.aB.ab, trace = TRUE)
rbind(ymat, sum(ymat)*fitted(fit))
Coef(fit)  # Estimated p
p <- sqrt(4*(fitted(fit)[, 4]))
p*p
summary(fit)
}
\keyword{models}
\keyword{regression}

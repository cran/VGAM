\name{AA.Aa.aa}
\alias{AA.Aa.aa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The AA-Aa-aa Blood Group System }
\description{
   Estimates the parameter of the
   AA-Aa-aa blood group system,
   with or without Hardy Weinberg equilibrium.

}
\usage{
AA.Aa.aa(linkp = "logit", linkf = "logit", inbreeding = FALSE,
         ipA = NULL, ifp = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{linkp, linkf}{
  Link functions applied to \code{pA} and \code{f}.
  See \code{\link{Links}} for more choices.

  }
  \item{ipA, ifp}{
  Optional initial values for \code{pA} and \code{f}.


  }
  \item{inbreeding}{
    Logical. Is there inbreeding?


%HWE assumption to be made?




  }
  \item{zero}{
    See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
  This one or two parameter model involves a probability called \code{pA}.
  The probability of getting a count in the first column of the
  input (an AA) is \code{pA*pA}.
  When \code{inbreeding = TRUE}, an additional parameter \code{f} is used.
  If \code{inbreeding = FALSE} then \eqn{f = 0} and Hardy-Weinberg
  Equilibrium (HWE) is assumed.
  The EIM is used if \code{inbreeding = FALSE}.



% With Hardy Weinberg equilibrium (HWE),
% Without the HWE assumption, an additional parameter \code{f} is used.



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

\section{Warning }{
  Setting \code{inbreeding = FALSE} makes estimation difficult
  with non-intercept-only models.
  Currently, this code seems to work with intercept-only models.


}


\seealso{
  \code{\link{AB.Ab.aB.ab}},
  \code{\link{ABO}},
  \code{\link{A1A2A3}},
  \code{\link{MNSs}}.


% \code{\link{AB.Ab.aB.ab2}},



}
\examples{
y <- cbind(53, 95, 38)
fit1 <- vglm(y ~ 1, AA.Aa.aa, trace = TRUE)
fit2 <- vglm(y ~ 1, AA.Aa.aa(inbreeding = TRUE), trace = TRUE)
rbind(y, sum(y) * fitted(fit1))
Coef(fit1)  # Estimated pA
Coef(fit2)  # Estimated pA and f
summary(fit1)
}
\keyword{models}
\keyword{regression}

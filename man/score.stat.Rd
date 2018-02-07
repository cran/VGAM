\name{score.stat}
\alias{score.stat}
\alias{score.stat.vlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Rao's Score Test
  Statistics Evaluated at the Null Values }
\description{
  Generic function that computes
  Rao's score test statistics evaluated at the null values
  (consequently they do not suffer from the Hauck-Donner effect).


}
\usage{
score.stat(object, ...)
score.stat.vlm(object, values0 = 0, subset = NULL, omit1s = TRUE,
          all.out = FALSE, iterate = TRUE, trace = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object, values0, subset}{
  Same as in \code{\link{wald.stat.vlm}}.


}
\item{omit1s, all.out, iterate}{
  Same as in \code{\link{wald.stat.vlm}}.


}
\item{trace}{
  Same as in \code{\link{wald.stat.vlm}}.


}
\item{\dots}{
  Ignored for now.


}
}
\details{
  The (Rao) \emph{score test}
  (also known as the \emph{Lagrange multiplier test} in econometrics)
  is a third general method for
  hypothesis testing under a likelihood-based framework
  (the others are the likelihood ratio test and
  Wald test; see \code{\link{lrt.stat}} and
  \code{\link{wald.stat}}).
  Asymptotically, the three tests are equivalent.
  The Wald test is not invariant to parameterization, and
  the usual Wald test statistics computed at the estimates
  make it vulnerable to the Hauck-Donner effect
  (HDE; see \code{\link{hdeff}}).
  This function is similar to \code{\link{wald.stat}} in that
  one coefficient is set to 0 (by default) and the other
  coefficients are iterated by IRLS to get their MLE subject to this
  constraint.


}
\value{
  By default the
  signed square root of the
  Rao score statistics are returned.
  If \code{all.out = TRUE} then a list is returned with the
  following components:
  \code{score.stat} the score statistic,
  \code{SE0} the standard error of that coefficient,
  \code{values0} the null values.
  Approximately, the default score statistics output are
  standard normal random variates if each null hypothesis is true.



}
%\references{
%
%}
\author{ Thomas W. Yee }

%\note{
%}

\section{Warning }{
  See \code{\link{wald.stat.vlm}}.


}


\seealso{
  \code{\link{wald.stat}},
  \code{\link{lrt.stat}},
  \code{\link{summaryvglm}},
  \code{\link[stats]{summary.glm}},
  \code{\link{vglm}},
  \code{\link{hdeff}}.



}
\examples{
set.seed(1)
pneumo <- transform(pneumo, let = log(exposure.time),
                            x3 = rnorm(nrow(pneumo)))
(pfit <- vglm(cbind(normal, mild, severe) ~ let + x3,
             propodds, data = pneumo))
score.stat(pfit)  # No HDE here; should be similar to the next line:
coef(summary(pfit))[, "z value"]  # Wald statistics computed at the MLE
summary(pfit, score0 = TRUE)
}
\keyword{models}
\keyword{regression}


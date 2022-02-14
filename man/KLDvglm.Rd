\name{KLD}
\alias{KLD}
\alias{KLDvglm}
%\alias{KLDvglm}
%\alias{KLDvgam}
%\alias{KLDrrvglm}  % 20211129
%\alias{KLDqrrvglm}
%\alias{KLDrrvgam}
%\alias{KLDc,vglm-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Kullback-Leibler Divergence }
\description{
  Calculates the Kullback-Leibler divergence
  for certain fitted model objects

}
\usage{
 KLD(object, \dots)
 KLDvglm(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Some \pkg{VGAM} object, for example, having
  class \code{\link{vglm-class}}.
  Currently \code{object} must be intercept-only.




  }
  \item{\dots}{
  Other possible arguments fed into
  \code{KLDvglm} in order to compute the KLD.


  }
}
\details{
  The \emph{Kullback-Leibler divergence} (KLD),
  or \emph{relative entropy},
  is a measure of how one probability distribution differs
  from a second reference probability distribution.
  Currently the \pkg{VGAM} package computes the KLD 
  for GAITD regression models
  (e.g., see \code{\link{gaitdpoisson}} and
  \code{\link{gaitdnbinomial}}) where the reference distribution
  is the (unscaled) parent or base distribution.
  For such, the formula for the KLD simplifies somewhat.
  Hence one can obtain a quantitative measure for the overall
  effect of altering, inflating, truncating and deflating certain
  (special) values.




}
\value{
  Returns a numeric nonnegative value with the corresponding KLD.
  A 0 value means no difference between an ordinary parent or base
  distribution.

  

}
\author{T. W. Yee. }
%\note{
%  Although it has had some testing and checking,
%  this function is still in the experimental phase.



% It is intended that improvements be made in the future
% such as to increase its efficiency.



%}

\references{

  Kullback, S. and Leibler, R. A. (1951).
  On information and sufficiency.
  \emph{Annals of Mathematical Statistics},
  \bold{22}, 79--86.



% \bold{22}(1), 79--86.



  M'Kendrick, A. G. (1925).
  Applications of mathematics to medical problems.
  \emph{Proc. Edinb. Math. Soc.},
  \bold{44}, 98--130.



}

\section{Warning }{
  Numerical problems might occur if any of the evaluated 
  probabilities
  of the unscaled parent distribution are very close to 0.


  
%  This code has been checked but is is still in the experimental phase.



}

\seealso{
  \code{\link{gaitdpoisson}},
  \code{\link{gaitdnbinomial}}.



}
\examples{
#  McKendrick (1925): Data from 223 Indian village households
cholera <- data.frame(ncases = 0:4,  # Number of cholera cases,
                      wfreq  = c(168, 32, 16, 6, 1))  # Frequencies
fit7 <- vglm(ncases ~ 1, gaitdpoisson(i.mlm = 0, ilambda.p = 1),
             weight = wfreq, data = cholera, trace = TRUE)
coef(fit7, matrix = TRUE)
KLD(fit7)
}
\keyword{models}
\keyword{regression}


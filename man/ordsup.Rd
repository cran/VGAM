\name{ordsup}
\alias{ordsup}
\alias{ordsup.vglm}
%\alias{score.stat}
%\alias{score.stat.vlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Ordinal Superiority
  Measures }
\description{
  Ordinal superiority measures for the linear model
  and cumulative link models: 
the probability that an observation from one
distribution falls above an independent observation
from the other distribution,
adjusted for explanatory variables in a model.



}
\usage{
ordsup(object, ...)
ordsup.vglm(object, all.vars = FALSE, confint = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{
  A \code{\link{vglm}} fit.
  Currently it must be one of:
  \code{\link{cumulative}},
  \code{\link{uninormal}}.
  The links for \code{\link{cumulative}} must be
  \code{\link{logitlink}} or \code{\link{probitlink}},
  and \code{parallel = TRUE} is also needed.
  For \code{\link{uninormal}} the mean must
  use \code{\link{identitylink}} and model the
  \code{sd} as intercept-only.


%  An object that is ideally an
%  \code{\link{vglm}} fit.


}
\item{all.vars}{
  Logical. The default is to use explanatory variables
  which are binary, but all variables are used (except the intercept)
  if set to \code{TRUE}.
  

}
\item{confint}{
  Logical.
  If \code{TRUE} then \code{\link{confintvglm}}
  is called to return confidence intervals for
  \eqn{\gamma}{gamma} and
  \eqn{\Delta}{Delta}.
  By default, Wald intervals are produced, but they
  can be replaced by profile intervals by setting
  \code{method = "profile"}.



  %Currently must be \code{FALSE}.

  

}
\item{\dots}{
  Parameters that can be fed into \code{\link{confintvglm}},
  e.g., \code{level = 0.95} and
  \code{method = c("wald", "profile")}.


}
}
\details{
  Details are given in Agresti and Kateri (2017) and this help
  file draws directly from this.
  This function returns two quantities for comparing two groups
  on an ordinal categorical response variable, while adjusting
  for other explanatory variables.
  They are called ``ordinal superiority'' measures, and
  the two groups can be compared without supplementary
  explanatory variables.
  Let \eqn{Y_1}{Y1} and \eqn{Y_2}{Y2} be independent random
  variables from groups A and B, say, for a quantitative ordinal
  categorical scale. Then
  \eqn{\Delta = P(Y_1 > Y_2) -
    P(Y_2 > Y_1)}{Delta = P(Y1 > Y2) - P(Y2 > Y1)} 
  summarizes their relative size.
  A second quantity is
  \eqn{\gamma = P(Y_1 > Y_2) -
    0.5 \times P(Y_2 = Y_1)}{gamma = P(Y1 > Y2) - 0.5 * P(Y2 = Y1)}.
  Then \eqn{\Delta=2 \times \gamma - 1}{Delta=2 * gamma -1}.
  whereas \eqn{\gamma=(\Delta + 1)/2}{gamma=(Delta + 1)/2}.
  The range of \eqn{\gamma}{gamma} is \eqn{[0, 1]}, while
  the range of \eqn{\Delta}{Delta} is \eqn{[-1, 1]}.
The examples below are based on that paper.
This function is currently implemented for a very limited
number of specific models.




}
\value{
  By default,
  a list with components
  \code{gamma} and
  \code{Delta},
  where each is a vector with elements corresponding to
  binary explanatory variables (i.e., 0 or 1),
  and if no explanatory variables are binary then a
  \code{NULL} is returned.
  If \code{confint = TRUE} then the list contains 4 more components:
  \code{lower.gamma},
  \code{upper.gamma},
  \code{Lower.Delta},
  \code{Upper.Delta}.



}
\references{

  Agresti, A. and Kateri, M. (2017).
  Ordinal probability effect measures for group
  comparisons in multinomial cumulative link models.
  \emph{Biometrics}, \bold{73}, 214--219.


}
\author{ Thomas W. Yee }

%\note{
%}


%\section{Warning }{
%  This function has not yet been thoroughly tested.


%}

\seealso{
  \code{\link{cumulative}},
  \code{\link{propodds}},
  \code{\link{uninormal}}.



}
\examples{
\dontrun{
Mental <- read.table("http://www.stat.ufl.edu/~aa/glm/data/Mental.dat",
                     header = TRUE)  # Make take a while to load in
Mental$impair <- ordered(Mental$impair)
pfit3 <- vglm(impair ~ ses + life, data = Mental,
         cumulative(link = "probitlink", reverse = FALSE, parallel = TRUE))
coef(pfit3, matrix = TRUE)
ordsup(pfit3)  # The 'ses' variable is binary

# Fit a crude LM
fit7 <- vglm(as.numeric(impair) ~ ses + life, uninormal, data = Mental)
coef(fit7, matrix = TRUE)  # 'sd' is estimated by MLE
ordsup(fit7)
ordsup(fit7, all.vars = TRUE)  # Some output may not be meaningful
ordsup(fit7, confint = TRUE, method = "profile")
}
}
\keyword{models}
\keyword{regression}


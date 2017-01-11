\name{brat}
\alias{brat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bradley Terry Model }
\description{
  Fits a Bradley Terry model (intercept-only model) by maximum
  likelihood estimation.

}
\usage{
brat(refgp = "last", refvalue = 1, ialpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refgp}{
  Integer whose value must be from the set
  \{1,\ldots,\eqn{M+1}\}, where there are \eqn{M+1}
  competitors. The default value indicates the last
  competitor is used---but don't input a character string,
  in general.


  }
  \item{refvalue}{
  Numeric. A positive value for the reference group.


  }
  \item{ialpha}{
  Initial values for the \eqn{\alpha}{alpha}s.
  These are recycled to the appropriate length.

  }

}
\details{
  The Bradley Terry model involves \eqn{M+1} competitors
  who either win or lose against each other (no draws/ties
  allowed in this implementation--see \code{\link{bratt}}
  if there are ties). The probability that Competitor
  \eqn{i} beats Competitor \eqn{j} is \eqn{\alpha_i /
  (\alpha_i+\alpha_j)}{alpha_i / (alpha_i + alpha_j)},
  where all the \eqn{\alpha}{alpha}s are positive.
  Loosely, the \eqn{\alpha}{alpha}s can be thought of as
  the competitors' `abilities'. For identifiability, one
  of the \eqn{\alpha_i}{alpha_i} is set to a known value
  \code{refvalue}, e.g., 1. By default, this function
  chooses the last competitor to have this reference value.
  The data can be represented in the form of a \eqn{M+1}
  by \eqn{M+1} matrix of counts, where winners are the
  rows and losers are the columns. However, this is not
  the way the data should be inputted (see below).


  Excluding the reference value/group, this function
  chooses \eqn{\log(\alpha_j)}{log(alpha_j)} as the
  \eqn{M} linear predictors. The log link ensures that
  the \eqn{\alpha}{alpha}s are positive.


  The Bradley Terry model can be fitted by logistic
  regression, but this approach is not taken here.
  The Bradley Terry model can be fitted with covariates,
  e.g., a home advantage variable, but unfortunately, this
  lies outside the VGLM theoretical framework and therefore
  cannot be handled with this code.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}.


}
\references{
Agresti, A. (2013)
\emph{Categorical Data Analysis},
3rd ed. Hoboken, NJ, USA: Wiley.


Stigler, S. (1994)
Citation patterns in the journals of statistics and probability.
\emph{Statistical Science},
\bold{9}, 94--108.


The \pkg{BradleyTerry2} package has more comprehensive capabilities
than this function.



}
\author{ T. W. Yee }
\note{
  The function \code{\link{Brat}} is useful for coercing
  a \eqn{M+1} by \eqn{M+1} matrix of counts into a one-row
  matrix suitable for \code{brat}. Diagonal elements are
  skipped, and the usual S order of \code{c(a.matrix)}
  of elements is used. There should be no missing values
  apart from the diagonal elements of the square matrix.
  The matrix should have winners as the rows, and losers
  as the columns. In general, the response should be a
  1-row matrix with \eqn{M(M+1)} columns.


  Only an intercept model is recommended with \code{brat}.
  It doesn't make sense really to include covariates because
  of the limited VGLM framework.


  Notationally, note that the \pkg{VGAM} family function
  \code{\link{brat}} has \eqn{M+1} contestants, while
  \code{bratt} has \eqn{M} contestants.


}
\section{Warning }{
  Presently, the residuals are wrong, and the prior weights
  are not handled correctly. Ideally, the total number of
  counts should be the prior weights, after the response has
  been converted to proportions. This would make it similar
  to family functions such as \code{\link{multinomial}}
  and \code{\link{binomialff}}.


}

\seealso{
  \code{\link{bratt}},
  \code{\link{Brat}},
  \code{\link{multinomial}},
  \code{\link{binomialff}}.


}
\examples{
# Citation statistics: being cited is a 'win'; citing is a 'loss'
journal <- c("Biometrika", "Comm.Statist", "JASA", "JRSS-B")
mat <- matrix(c( NA, 33, 320, 284,
                730, NA, 813, 276,
                498, 68,  NA, 325,
                221, 17, 142,  NA), 4, 4)
dimnames(mat) <- list(winner = journal, loser = journal)
fit <- vglm(Brat(mat) ~ 1, brat(refgp = 1), trace = TRUE)
fit <- vglm(Brat(mat) ~ 1, brat(refgp = 1), trace = TRUE, crit = "coef")
summary(fit)
c(0, coef(fit))  # Log-abilities (in order of "journal")
c(1, Coef(fit))  # Abilities (in order of "journal")
fitted(fit)     # Probabilities of winning in awkward form
(check <- InverseBrat(fitted(fit)))  # Probabilities of winning
check + t(check)  # Should be 1's in the off-diagonals
}
\keyword{models}
\keyword{regression}

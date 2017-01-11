\name{bratt}
\alias{bratt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bradley Terry Model With Ties }
\description{
  Fits a Bradley Terry model with ties (intercept-only model)
  by maximum likelihood estimation.

}
\usage{
bratt(refgp = "last", refvalue = 1, ialpha = 1, i0 = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refgp}{
  Integer whose value must be from the set \{1,\ldots,\eqn{M}\},
  where there are \eqn{M} competitors. The default value
  indicates the last competitor is used---but don't input
  a character string, in general.


  }
  \item{refvalue}{
  Numeric. A positive value for the reference group.


  }
  \item{ialpha}{
  Initial values for the \eqn{\alpha}{alpha}s.
  These are recycled to the appropriate length.


  }
  \item{i0}{
  Initial value for \eqn{\alpha_0}{alpha_0}.
  If convergence fails, try another positive value.


  }
}
\details{
  There are several models that extend the ordinary
  Bradley Terry model to handle ties. This family function
  implements one of these models.  It involves \eqn{M}
  competitors who either win or lose or tie against
  each other.  (If there are no draws/ties then use
  \code{\link{brat}}).  The probability that Competitor
  \eqn{i} beats Competitor \eqn{j} is \eqn{\alpha_i /
  (\alpha_i+\alpha_j+\alpha_0)}{alpha_i / (alpha_i +
  alpha_j + alpha_0)}, where all the \eqn{\alpha}{alpha}s
  are positive.  The probability that Competitor \eqn{i}
  ties with Competitor \eqn{j} is \eqn{\alpha_0 /
  (\alpha_i+\alpha_j+\alpha_0)}{alpha_0 / (alpha_i +
  alpha_j + alpha_0)}.  Loosely, the \eqn{\alpha}{alpha}s
  can be thought of as the competitors' `abilities',
  and \eqn{\alpha_0}{alpha_0} is an added parameter
  to model ties.  For identifiability, one of the
  \eqn{\alpha_i}{alpha_i} is set to a known value
  \code{refvalue}, e.g., 1.  By default, this function
  chooses the last competitor to have this reference value.
  The data can be represented in the form of a \eqn{M}
  by \eqn{M} matrix of counts, where winners are the rows
  and losers are the columns.  However, this is not the
  way the data should be inputted (see below).


  Excluding the reference value/group, this function
  chooses \eqn{\log(\alpha_j)}{log(alpha_j)} as the first
  \eqn{M-1} linear predictors.  The log link ensures that
  the \eqn{\alpha}{alpha}s are positive.  The last linear
  predictor is \eqn{\log(\alpha_0)}{log(alpha_0)}.


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

  Torsney, B. (2004)
  Fitting Bradley Terry models using a multiplicative algorithm.
  In: Antoch, J. (ed.)
  \emph{Proceedings in Computational Statistics COMPSTAT 2004},
  Physica-Verlag: Heidelberg. Pages 513--526.


}
\author{ T. W. Yee }
\note{
  The function \code{\link{Brat}} is useful for coercing
  a \eqn{M} by \eqn{M} matrix of counts into a one-row
  matrix suitable for \code{bratt}.  Diagonal elements
  are skipped, and the usual S order of \code{c(a.matrix)}
  of elements is used. There should be no missing values
  apart from the diagonal elements of the square matrix.
  The matrix should have winners as the rows, and losers
  as the columns.  In general, the response should be a
  matrix with \eqn{M(M-1)} columns.


  Also, a symmetric matrix of ties should be passed into
  \code{\link{Brat}}. The diagonal of this matrix should
  be all \code{NA}s.


  Only an intercept model is recommended with \code{bratt}.
  It doesn't make sense really to include covariates because
  of the limited VGLM framework.


  Notationally, note that the \pkg{VGAM} family function
  \code{\link{brat}} has \eqn{M+1} contestants, while
  \code{bratt} has \eqn{M} contestants.


}

\seealso{
  \code{\link{brat}},
  \code{\link{Brat}},
  \code{\link{binomialff}}.


}
\examples{
# citation statistics: being cited is a 'win'; citing is a 'loss'
journal <- c("Biometrika", "Comm.Statist", "JASA", "JRSS-B")
mat <- matrix(c( NA, 33, 320, 284,
                730, NA, 813, 276,
                498, 68,  NA, 325,
                221, 17, 142,  NA), 4, 4)
dimnames(mat) <- list(winner = journal, loser = journal)

# Add some ties. This is fictitional data.
ties <- 5 + 0 * mat
ties[2, 1] <- ties[1,2] <- 9

# Now fit the model
fit <- vglm(Brat(mat, ties) ~ 1, bratt(refgp = 1), trace = TRUE)
fit <- vglm(Brat(mat, ties) ~ 1, bratt(refgp = 1), trace = TRUE, crit = "coef")

summary(fit)
c(0, coef(fit))  # Log-abilities (in order of "journal"); last is log(alpha0)
c(1, Coef(fit))  #     Abilities (in order of "journal"); last is alpha0

fit@misc$alpha   # alpha_1,...,alpha_M
fit@misc$alpha0  # alpha_0

fitted(fit)  # Probabilities of winning and tying, in awkward form
predict(fit)
(check <- InverseBrat(fitted(fit)))    # Probabilities of winning
qprob <- attr(fitted(fit), "probtie")  # Probabilities of a tie
qprobmat <- InverseBrat(c(qprob), NCo = nrow(ties))  # Probabilities of a tie
check + t(check) + qprobmat  # Should be 1s in the off-diagonals
}
\keyword{models}
\keyword{regression}

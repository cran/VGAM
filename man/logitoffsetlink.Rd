\name{logitoffsetlink}
\alias{logitoffsetlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Logit-with-an-Offset Link Function }
\description{
  Computes the logitoffsetlink transformation, including its inverse and the
  first two derivatives.

}
\usage{
logitoffsetlink(theta, offset = 0, inverse = FALSE, deriv = 0,
      short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.


  }
  \item{offset}{
    The offset value(s), which must be non-negative.
    It is called \eqn{K} below.


  }

  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  This  link function allows for some asymmetry compared to the
  ordinary \code{\link{logit}} link.
  The formula is
  \deqn{\log(\theta/(1-\theta) - K)}{%
        log(theta/(1-theta) - K)}
  and the default value for the offset \eqn{K} is corresponds to the
  ordinary \code{\link{logit}} link.
  When \code{inverse = TRUE} will mean that the value will
  lie in the interval \eqn{(K / (1+K), 1)}.


}
\value{
  For \code{logitoffsetlink} with \code{deriv = 0}, the
  logitoffsetlink of \code{theta}, i.e.,
  \code{log(theta/(1-theta) - K)} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{(K + exp(theta))/(1 + exp(theta) + K)}.



  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.



  Here, all logarithms are natural logarithms, i.e., to base \emph{e}.


}
\references{
  Komori, O. and Eguchi, S. et al., 2016.
  An asymmetric logistic model for ecological data.
  \emph{Methods in Ecology and Evolution},
  \bold{7}.


}
\author{ Thomas W. Yee }

\note{
  This function is numerical less stability than
  \code{\link{logit}}.


}

\seealso{
    \code{\link{Links}},
    \code{\link{logit}}.


}
\examples{
p <- seq(0.05, 0.99, by = 0.01); myoff <- 0.05
logitoffsetlink(p, myoff)
max(abs(logitoffsetlink(logitoffsetlink(p, myoff),
                        myoff, inverse = TRUE) - p))  # Should be 0
}
\keyword{math}
\keyword{models}
\keyword{regression}







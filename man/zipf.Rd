\name{zipf}
\alias{zipf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zipf Distribution Family Function }
\description{
  Estimates the parameter of the Zipf distribution.

}
\usage{
zipf(N = NULL, lshape = "loge", ishape = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
  Number of elements, an integer satisfying \code{1 < N < Inf}.
  The default is to use the maximum value of the response.
  If given, \code{N} must be no less that the largest response value.
  If \code{N = Inf} and \eqn{s>1} then this is the zeta
  distribution (use \code{\link{zetaff}} instead).


  }
  \item{lshape}{
  Parameter link function applied to the (positive) shape parameter \eqn{s}.
  See \code{\link{Links}} for more choices.


  }
  \item{ishape}{
  Optional initial value for the parameter \eqn{s}.
  The default is to choose an initial value internally.
  If converge failure occurs use this argument to input a value.


  }
}
\details{
  The probability function for a response \eqn{Y} is
  \deqn{P(Y=y) = y^{-s} / \sum_{i=1}^N i^{-s},\ \ s>0,\ \ y=1,2,\ldots,N,}{%
        P(Y=y) = (y^(-s)) / sum((1:N)^(-s)), s>0, y=1,2,...,N,}
  where \eqn{s} is the exponent characterizing the distribution.
  The mean of \eqn{Y}, which are returned as the fitted values,
  is \eqn{\mu = H_{N,s-1} / H_{N,s}}{H(N,s-1) / H(N,s)}
  where \eqn{H_{n,m}= \sum_{i=1}^n i^{-m}}{H(n,m)=sum((1:n)^(-m))}
  is the \eqn{n}th generalized harmonic number.



  Zipf's law is an experimental law which is often applied
  to the study of the frequency of words in a corpus of
  natural language utterances. It states that the frequency
  of any word is inversely proportional to its rank in the
  frequency table. For example, \code{"the"} and \code{"of"}
  are first two most common words, and Zipf's law states
  that \code{"the"} is twice as common as \code{"of"}.
  Many other natural phenomena conform to Zipf's law.



}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as
  \code{\link{vglm}} and \code{\link{vgam}}.


}

\references{

  pp.526-- of Chapter 11 of
  Johnson N. L., Kemp, A. W. and Kotz S. (2005)
  \emph{Univariate Discrete Distributions},
  3rd edition,
  Hoboken, New Jersey, USA: Wiley.


}
\author{ T. W. Yee }
\note{
  Upon convergence, the \code{N} is stored as \code{@misc$N}.


}

\seealso{
  \code{\link{dzipf}},
  \code{\link{zetaff}},
  \code{\link{simulate.vlm}}.


}
\examples{
zdata <- data.frame(y = 1:5, ofreq = c(63, 14, 5, 1, 2))
zfit <- vglm(y ~ 1, zipf, data = zdata, trace = TRUE, weight = ofreq)
zfit <- vglm(y ~ 1, zipf(lshape = "identitylink", ishape = 3.4), data = zdata,
            trace = TRUE, weight = ofreq, crit = "coef")
zfit@misc$N
(shape.hat <- Coef(zfit))
with(zdata, weighted.mean(y, ofreq))
fitted(zfit, matrix = FALSE)
}
\keyword{models}
\keyword{regression}


%pp.465--471, Chapter 11 of
%Johnson N. L., Kotz S., and Kemp A. W. (1993)
%\emph{Univariate Discrete Distributions},
%2nd ed.
%New York: Wiley.

%http://www.math.uah.edu/stat/special/Zeta.html calls s 'shape'




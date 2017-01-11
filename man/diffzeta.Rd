\name{diffzeta}
\alias{diffzeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Differenced Zeta Distribution Family Function }
\description{
  Estimates the parameter of the differenced zeta distribution.

}
\usage{
diffzeta(start = 1, lshape = "loge", ishape = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, ishape}{
  Same as \code{\link{zetaff}}.


  }
  \item{start}{
    Smallest value of the support of the distribution.
    Must be a positive integer.


  }
}
\details{
  The PMF is
  \deqn{P(Y=y) = (a/y)^{s} - (a/(1+y))^{s},\ \ s>0,\ \ y=a,a+1,\ldots,}{%
        P(Y=y) = (a/y)^(s) - / (a/(1+y))^(s), s>0, y=a,a+1,...,}
  where \eqn{s} is the positive shape parameter, and \eqn{a} is \code{start}.
  According to   Moreno-Sanchez et al. (2016), this model
  fits quite well to about 40 percent of all the English books in the
  Project Gutenberg data base (about 30,000 texts).
  Multiple responses are handled.



}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.



}

\references{
Moreno-Sanchez, I. and Font-Clos, F. and Corral, A.
Large-Scale Analysis of Zipf's Law in English Texts,
2016.
PLoS ONE, \bold{11}(1), 1--19.



}
\author{ T. W. Yee }
%\note{
%  The \code{\link{zeta}} function may be used to compute values
%  of the zeta function.
%
%
%}

\seealso{
  \code{\link{Diffzeta}},
  \code{\link{zetaff}},
  \code{\link{zeta}},
  \code{\link{zipf}},
  \code{\link{zipf}}.


}
\examples{
odata <- data.frame(x2 = runif(nn <- 1000))  # Artificial data
odata <- transform(odata, shape = loge(-0.25 + x2, inverse = TRUE))
odata <- transform(odata, y1 = rdiffzeta(nn, shape))
with(odata, table(y1))
ofit <- vglm(y1 ~ x2, diffzeta, data = odata, trace = TRUE, crit = "coef")
coef(ofit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}
%

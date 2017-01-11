\name{rec.exp1}
\alias{rec.exp1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Upper Record Values from a 1-parameter Exponential Distribution }
\description{
  Maximum likelihood estimation of the rate parameter of a
  1-parameter exponential distribution when the observations are upper
  record values.

}
\usage{
rec.exp1(lrate = "loge", irate = NULL, imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lrate}{
  Link function applied to the rate parameter.
  See \code{\link{Links}} for more choices.


  }
  \item{irate}{
  Numeric. Optional initial values for the rate.
  The default value \code{NULL} means they are computed internally,
  with the help of \code{imethod}.


  }
  \item{imethod}{
  Integer, either 1 or 2 or 3. Initial method, three algorithms are
  implemented. Choose the another value if convergence fails, or use
  \code{irate}.


  }
}
\details{
  The response must be a vector or one-column matrix with strictly
  increasing values.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
    Arnold, B. C. and Balakrishnan, N. and Nagaraja, H. N. (1998)
    \emph{Records},
    New York: John Wiley & Sons.


}
\author{ T. W. Yee }
\note{
  By default, this family function has the intercept-only MLE as the
  initial value, therefore convergence may only take one iteration.
  Fisher scoring is used.


}

\seealso{
  \code{\link{exponential}}.


}
\examples{
rawy <- rexp(n <- 10000, rate = exp(1))
y <- unique(cummax(rawy))  # Keep only the records

length(y) / y[length(y)]   # MLE of rate

fit <- vglm(y ~ 1, rec.exp1, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{models}
\keyword{regression}

%# Keep only the records
%delete = c(FALSE, rep(TRUE, len = n-1))
%for (i in 2:length(rawy))
%    if (rawy[i] > max(rawy[1:(i-1)])) delete[i] = FALSE
%(y = rawy[!delete])




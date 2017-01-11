\name{studentt}
\alias{studentt}
\alias{studentt2}
\alias{studentt3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Student t Distribution }
\description{
  Estimating the parameters of a Student t distribution.


}
\usage{
studentt (ldf = "loglog", idf = NULL, tol1 = 0.1, imethod = 1)
studentt2(df = Inf, llocation = "identitylink", lscale = "loge",
          ilocation = NULL, iscale = NULL, imethod = 1, zero = "scale")
studentt3(llocation = "identitylink", lscale = "loge", ldf = "loglog",
          ilocation = NULL, iscale = NULL, idf = NULL,
          imethod = 1, zero = c("scale", "df"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale, ldf}{
  Parameter link functions for each parameter,
  e.g., for degrees of freedom \eqn{\nu}{nu}.
  See \code{\link{Links}} for more choices.
  The defaults ensures the parameters are in range.
  A \code{\link{loglog}} link keeps the degrees of freedom greater
  than unity; see below.


  }
  \item{ilocation, iscale, idf}{
  Optional initial values.
  If given, the values must be in range.
  The default is to compute an initial value internally.


  }

  \item{tol1}{
  A positive value, the tolerance for testing whether an
  initial value is 1.
  Best to leave this argument alone.


  }

  \item{df}{
  Numeric, user-specified degrees of freedom.
  It may be of length equal to the number of columns of a response
  matrix.


  }
  \item{imethod, zero}{
  See \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The Student t density function is
  \deqn{f(y;\nu) = \frac{\Gamma((\nu+1)/2)}{\sqrt{\nu \pi} \Gamma(\nu/2)}
               \left(1 + \frac{y^2}{\nu} \right)^{-(\nu+1)/2}}{%
        f(y;nu) = (gamma((nu+1)/2) / (sqrt(nu*pi) gamma(nu/2))) *
               (1 + y^2 / nu)^{-(nu+1)/2}}
  for all real \eqn{y}.
  Then \eqn{E(Y)=0} if \eqn{\nu>1}{nu>1} (returned as the fitted values),
  and \eqn{Var(Y)= \nu/(\nu-2)}{Var(Y)= nu/(nu-2)}
  for \eqn{\nu > 2}{nu > 2}.
  When \eqn{\nu=1}{nu=1} then the Student \eqn{t}-distribution
  corresponds to the standard Cauchy distribution,
  \code{\link{cauchy1}}.
  When \eqn{\nu=2}{nu=2} with a scale parameter of \code{sqrt(2)} then
  the Student \eqn{t}-distribution
  corresponds to the standard (Koenker) distribution,
  \code{\link{sc.studentt2}}.
  The degrees of freedom can be treated as a parameter to be estimated,
  and as a real and not an integer.
  The Student t distribution is used for a variety of reasons
  in statistics, including robust regression.


  Let \eqn{Y = (T - \mu) / \sigma}{Y = (T -  mu) /  sigma} where
  \eqn{\mu}{mu} and \eqn{\sigma}{sigma} are the location
  and scale parameters respectively.
  Then \code{studentt3} estimates the location, scale and
  degrees of freedom parameters.
  And \code{studentt2} estimates the location, scale parameters
  for a user-specified degrees of freedom, \code{df}.
  And \code{studentt} estimates the degrees of freedom parameter only.
  The fitted values are the location parameters.
  By default the linear/additive predictors are
  \eqn{(\mu, \log(\sigma), \log\log(\nu))^T}{(mu, log(sigma), log log(nu))^T}
  or subsets thereof.


  In general convergence can be slow, especially when there are
  covariates.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{


Student (1908)
The probable error of a mean.
\emph{Biometrika}, \bold{6}, 1--25.


Zhu, D. and Galbraith, J. W. (2010)
A generalized asymmetric Student-\emph{t} distribution with
application to financial econometrics.
\emph{Journal of Econometrics}, \bold{157}, 297--305.


}

\author{ T. W. Yee }
\note{
  \code{studentt3()} and \code{studentt2()} can handle multiple responses.


  Practical experience has shown reasonably good initial values are
  required. If convergence failure occurs try using arguments
  such as \code{idf}.
  Local solutions are also possible, especially when
  the degrees of freedom is close to unity or
  the scale parameter is close to zero.


  A standard normal distribution corresponds to a \emph{t} distribution
  with infinite degrees of freedom. Consequently, if the data is close
  to normal, there may be convergence problems; best to use
  \code{\link{uninormal}} instead.


}
\seealso{
  \code{\link{uninormal}},
  \code{\link{cauchy1}},
  \code{\link{logistic}},
  \code{\link{huber2}},
  \code{\link{sc.studentt2}},
  \code{\link[stats]{TDist}},
  \code{\link{simulate.vlm}}.


}
\examples{
tdata <- data.frame(x2 = runif(nn <- 1000))
tdata <- transform(tdata, y1 = rt(nn, df = exp(exp(0.5 - x2))),
                          y2 = rt(nn, df = exp(exp(0.5 - x2))))
fit1 <- vglm(y1 ~ x2, studentt, data = tdata, trace = TRUE)
coef(fit1, matrix = TRUE)

fit2 <- vglm(y1 ~ x2, studentt2(df = exp(exp(0.5))), data = tdata)
coef(fit2, matrix = TRUE)  # df inputted into studentt2() not quite right

fit3 <- vglm(cbind(y1, y2) ~ x2, studentt3, data = tdata, trace = TRUE)
coef(fit3, matrix = TRUE)
}
\keyword{models}
\keyword{regression}

%Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
%\emph{Statistical Distributions},
%Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.



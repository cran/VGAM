\name{freund61}
\alias{freund61}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Freund's (1961) Bivariate Extension of the Exponential Distribution }
\description{
  Estimate the four parameters of the Freund (1961) bivariate extension
  of the exponential distribution
  by maximum likelihood estimation.

}
\usage{
freund61(la = "loge",  lap = "loge",  lb = "loge", lbp = "loge",
         ia = NULL, iap = NULL, ib = NULL, ibp = NULL,
         independent = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{la, lap, lb, lbp}{
  Link functions applied to the (positive)
  parameters \eqn{\alpha}{alpha}, \eqn{\alpha'}{alpha'},
  \eqn{\beta}{beta} and \eqn{\beta'}{beta'}, respectively
  (the ``\code{p}'' stands for ``prime'').
  See \code{\link{Links}} for more choices.


  }
  \item{ia, iap, ib, ibp}{
  Initial value for the four parameters respectively.
  The default is to estimate them all internally.


  }
  \item{independent}{
  Logical. If \code{TRUE} then the parameters are constrained to satisfy
  \eqn{\alpha=\alpha'}{alpha=alpha'} and \eqn{\beta=\beta'}{beta=beta'},
  which implies that \eqn{y_1}{y1} and \eqn{y_2}{y2} are independent
  and each have an ordinary exponential distribution.


  }
  \item{zero}{
  A vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values can be from the set \{1,2,3,4\}.
  The default is none of them.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  This model represents one type of bivariate extension of the exponential
  distribution that is applicable to certain problems, in particular,
  to two-component systems which can function if one of the components
  has failed. For example, engine failures in two-engine planes, paired
  organs such as peoples' eyes, ears and kidneys.
  Suppose \eqn{y_1}{y1} and \eqn{y_2}{y2} are random variables
  representing the lifetimes of two components \eqn{A} and \eqn{B}
  in a two component system.
  The dependence between \eqn{y_1}{y1} and \eqn{y_2}{y2}
  is essentially such that the failure of the \eqn{B} component
  changes the parameter of the exponential life distribution
  of the \eqn{A}  component from \eqn{\alpha}{alpha} to
  \eqn{\alpha'}{alpha'}, while the failure of the \eqn{A}  component
  changes the parameter of the exponential life distribution
  of the \eqn{B}  component from \eqn{\beta}{beta} to
  \eqn{\beta'}{beta'}.


  The joint probability density function is given by
  \deqn{f(y_1,y_2) = \alpha \beta' \exp(-\beta' y_2 -
                      (\alpha+\beta-\beta')y_1) }{%
        f(y1,y2) = alpha * beta' * exp(-beta' * y2 -
                      (alpha+beta-beta') * y1) }
  for \eqn{0 < y_1 < y_2}{0 < y1 < y2}, and
  \deqn{f(y_1,y_2) = \beta \alpha' \exp(-\alpha' y_1 -
                      (\alpha+\beta-\alpha')y_2) }{%
        f(y1,y2) = beta * alpha' * exp(-alpha' * y1 -
                      (alpha+beta-alpha') * y2) }
  for \eqn{0 < y_2 < y_1}{0 < y2 < y1}.
  Here, all four parameters are positive, as well as the responses
  \eqn{y_1}{y1} and \eqn{y_2}{y2}.
  Under this model, the probability that component \eqn{A}
  is the first to fail is
  \eqn{\alpha/(\alpha+\beta)}{alpha/(alpha+beta)}.
  The time to the first failure is distributed as an
  exponential distribution with rate
  \eqn{\alpha+\beta}{alpha+beta}. Furthermore, the
  distribution of the time from first failure to failure
  of the other component is a mixture of
  Exponential(\eqn{\alpha'}{alpha'}) and
  Exponential(\eqn{\beta'}{beta'}) with proportions
  \eqn{\beta/(\alpha+\beta)}{beta/(alpha+beta)}
  and \eqn{\alpha/(\alpha+\beta)}{alpha/(alpha+beta)}
  respectively.


  The marginal distributions are, in general, not exponential.
  By default, the linear/additive predictors are
  \eqn{\eta_1=\log(\alpha)}{eta1=log(alpha)},
  \eqn{\eta_2=\log(\alpha')}{eta2=log(alpha')},
  \eqn{\eta_3=\log(\beta)}{eta3=log(beta)},
  \eqn{\eta_4=\log(\beta')}{eta4=log(beta')}.


  A special case is when \eqn{\alpha=\alpha'}{alpha=alpha'}
  and \eqn{\beta=\beta'}{beta'=beta'}, which means that
  \eqn{y_1}{y1} and \eqn{y_2}{y2} are independent, and
  both have an ordinary exponential distribution with means
  \eqn{1 / \alpha}{1/alpha} and \eqn{1 / \beta}{1/beta}
  respectively.


  Fisher scoring is used,
  and the initial values correspond to the MLEs of an intercept model.
  Consequently, convergence may take only one iteration.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}

\references{

Freund, J. E. (1961)
A bivariate extension of the exponential distribution.
\emph{Journal of the American Statistical Association},
\bold{56}, 971--977.


}
\author{ T. W. Yee }
\note{
  To estimate all four parameters, it is necessary to have some
  data where \eqn{y_1<y_2}{y1<y2} and \eqn{y_2<y_1}{y2<y1}.


  The response must be a two-column matrix, with columns
  \eqn{y_1}{y1} and \eqn{y_2}{y2}.
  Currently, the fitted value is a matrix with two columns; the
  first column has values
  \eqn{(\alpha'+\beta)/(\alpha' (\alpha+\beta))}{(alpha'+beta)/(alpha' * (alpha+beta))}
  for the mean of \eqn{y_1}{y1},
  while the second column has values
  \eqn{(\beta'+\alpha)/(\beta' (\alpha+\beta))}{(beta'+alpha)/(beta' * (alpha+beta))}
  for the mean of \eqn{y_2}{y2}.
  The variance of \eqn{y_1}{y1} is
  \deqn{ \frac{(\alpha')^2 + 2 \alpha \beta + \beta^2}{
               (\alpha')^2 (\alpha + \beta)^2}, }{%
              [(alpha')^2 + 2 * alpha * beta + beta^2]/
              [(alpha')^2 * (alpha + beta)^2], }
  the variance of \eqn{y_2}{y2} is
  \deqn{ \frac{(\beta')^2 + 2 \alpha \beta + \alpha^2 }{
               (\beta')^2 (\alpha + \beta)^2 }, }{%
              [(beta')^2 + 2 * alpha * beta + alpha^2]/
              [(beta')^2 * (alpha + beta)^2], }
  the covariance of \eqn{y_1}{y1} and \eqn{y_2}{y2} is
  \deqn{ \frac{\alpha' \beta' - \alpha \beta }{
               \alpha' \beta' (\alpha + \beta)^2}. }{%
              [alpha' * beta' - alpha * beta]/
              [alpha' * beta' * (alpha + beta)^2]. }


}
%\section{Warning}{
%}

\seealso{
  \code{\link{exponential}}.


}
\examples{
fdata <- data.frame(y1 = rexp(nn <- 1000, rate = exp(1)))
fdata <- transform(fdata, y2 = rexp(nn, rate = exp(2)))
fit1 <- vglm(cbind(y1, y2) ~ 1, fam = freund61, data = fdata, trace = TRUE)
coef(fit1, matrix = TRUE)
Coef(fit1)
vcov(fit1)
head(fitted(fit1))
summary(fit1)

# y1 and y2 are independent, so fit an independence model
fit2 <- vglm(cbind(y1, y2) ~ 1, freund61(indep = TRUE),
             data = fdata, trace = TRUE)
coef(fit2, matrix = TRUE)
constraints(fit2)
pchisq(2 * (logLik(fit1) - logLik(fit2)),  # p-value
       df = df.residual(fit2) - df.residual(fit1), lower.tail = FALSE)
lrtest(fit1, fit2)  # Better alternative
}
\keyword{models}
\keyword{regression}


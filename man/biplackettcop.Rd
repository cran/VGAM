\name{biplackettcop}
\alias{biplackettcop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plackett's Bivariate Copula Family Function }
\description{
  Estimate the association parameter of Plackett's bivariate distribution
  (copula)
  by maximum likelihood estimation.


}
\usage{
biplackettcop(link = "loge", ioratio = NULL, imethod = 1, nsimEIM = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to the (positive) odds ratio \eqn{\psi}{psi}.
  See \code{\link{Links}} for more choices
  and information.


  }
  \item{ioratio}{
  Numeric. Optional initial value for \eqn{\psi}{psi}.
  If a convergence failure occurs try assigning a value or a different value.


  }
  \item{imethod, nsimEIM}{
  See \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The defining equation is
  \deqn{\psi = H \times (1-y_1-y_2+H) / ((y_1-H) \times (y_2-H))}{%
         psi = H*(1-y1-y2+H) / ((y1-H)*(y2-H))}
  where
  \eqn{P(Y_1 \leq y_1, Y_2 \leq y_2) = H_{\psi}(y_1,y_2)}{P(Y1 <= y1, Y2 <= y2)=
       H(y1,y2)}
  is the cumulative distribution function.
  The density function is \eqn{h_{\psi}(y_1,y_2) =}{h(y1,y2) =}
  \deqn{\psi [1 + (\psi-1)(y_1 + y_2 - 2 y_1 y_2) ] / \left(
                     [1 + (\psi-1)(y_1 + y_2) ]^2 - 4 \psi
              (\psi-1) y_1 y_2 \right)^{3/2}}{%
              psi*[1 + (psi-1)*(y1 + y2 - 2*y1*y2) ] / (
                     [1 + (psi-1)*(y1 + y2)]^2 -
              4*psi*(psi-1)*y1*y2)^(3/2)}
  for \eqn{\psi > 0}{psi > 0}.
  Some writers call \eqn{\psi}{psi} the \emph{cross product ratio}
  but it is called the \emph{odds ratio} here.
  The support of the function is the unit square.
  The marginal distributions here are the standard uniform although
  it is commonly generalized to other distributions.


  If \eqn{\psi = 1}{psi=1} then
  \eqn{h_{\psi}(y_1,y_2) = y_1 y_2}{h(y1,y2) = y1*y2},
  i.e., independence.
  As the odds ratio tends to infinity one has \eqn{y_1=y_2}{y1=y2}.
  As the odds ratio tends to 0 one has \eqn{y_2=1-y_1}{y2=1-y1}.


  Fisher scoring is implemented using \code{\link{rbiplackcop}}.
  Convergence is often quite slow.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}

\references{


Plackett, R. L. (1965)
A class of bivariate distributions.
\emph{Journal of the American Statistical Association},
\bold{60}, 516--522.


}
\author{ T. W. Yee }
\note{
  The response must be a two-column matrix.  Currently, the fitted
  value is a 2-column matrix with 0.5 values
  because the marginal distributions correspond to a standard
  uniform distribution.


}

\seealso{
  \code{\link{rbiplackcop}},
  \code{\link{bifrankcop}}.


}
\examples{
\dontrun{
ymat <- rbiplackcop(n = 2000, oratio = exp(2))
plot(ymat, col = "blue")
fit <- vglm(ymat ~ 1, fam = biplackettcop, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
vcov(fit)
head(fitted(fit))
summary(fit)
}
}
\keyword{models}
\keyword{regression}


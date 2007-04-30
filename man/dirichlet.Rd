\name{dirichlet}
\alias{dirichlet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting a Dirichlet Distribution }
\description{
  Fits a Dirichlet distribution to a matrix of compositions.

}
\usage{
dirichlet(link = "loge", earg=list(), zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  In the following, the response is assumed to be a \eqn{M}-column
  matrix with positive values and whose rows each sum to unity. Such
  data can be thought of as compositional data.
  There are \eqn{M} linear/additive predictors \eqn{\eta_j}{eta_j}.

  \item{link}{
  Link function applied to each of the \eqn{M} (positive) shape
  parameters \eqn{\alpha_j}{alpha_j}.
  See \code{\link{Links}} for more choices.
  The default gives \eqn{\eta_j=\log(\alpha_j)}{eta_j=log(alpha_j)}.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The default is none of them.
  If used, choose values from the set \eqn{\{1,2,\ldots,M\}}.

  }
}
\details{
  The Dirichlet distribution is commonly used to model compositional
  data, including applications in genetics.
  Suppose \eqn{(Y_1,\ldots,Y_{M})^T}{(Y_1,\ldots,Y_M)^T} is
  the response. Then it has a Dirichlet distribution if
  \eqn{(Y_1,\ldots,Y_{M-1})^T}{(Y_1,\ldots,Y_{M-1})^T} has density
  \deqn{\frac{\Gamma(\alpha_{+})}
           {\prod_{j=1}^{M} \Gamma(\alpha_{j})}
           \prod_{j=1}^{M} y_j^{\alpha_{j} -1}}{%
           (Gamma(alpha_+) / prod_{j=1}^M gamma(alpha_j))
           prod_{j=1}^M y_j^(alpha_j -1)}
  where \eqn{\alpha_+=\alpha_1+\cdots+\alpha_M}{alpha_+= alpha_1 + \dots + alpha_M},
  \eqn{\alpha_j > 0}{alpha_j > 0},
  and the density is defined on the unit simplex
  \deqn{\Delta_{M} = \left\{
  (y_1,\ldots,y_{M})^T :
  y_1 > 0, \ldots, y_{M} > 0,
  \sum_{j=1}^{M} y_j = 1 \right\}. }{%
  Delta_M = {
  (y_1,\ldots,y_M)^T :
  y_1 > 0, \dots, y_M > 0,
  \sum_{j=1}^M y_j = 1 }.
  }
  One has \eqn{E(Y_j) = \alpha_j / \alpha_{+}}{E(Y_j) = alpha_j /
  alpha_{+}}, which are returned as the fitted values.
  For this distribution Fisher scoring corresponds to Newton-Raphson.

  The Dirichlet distribution can be motivated by considering the random variables
  \eqn{(G_1,\ldots,G_{M})^T}{(G_1,\ldots,G_M)^T} which are each independent
  and identically distributed as a gamma distribution with density 
  \eqn{f(g_j)=g_j^{\alpha_j - 1} e^{-g_j} / \Gamma(\alpha_j)}{f(g_j)=
              g_j^(alpha_j - 1) e^(-g_j) / gamma(alpha_j)}.
  Then the Dirichlet distribution arises when
  \eqn{Y_j=G_j / (G_1 + \cdots + G_M)}{Y_j = G_j / (G_1 + ... + G_M)}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.

  When fitted, the \code{fitted.values} slot of the object contains the
  \eqn{M}-column matrix of means. 

}
\references{
Lange, K. (2002)
\emph{Mathematical and Statistical Methods for Genetic Analysis},
2nd ed. New York: Springer-Verlag.

Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.

%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.

}
\author{ Thomas W. Yee }
\note{
  The response should be a matrix of positive values whose rows
  each sum to unity. Similar to this is count data, where probably a
  multinomial logit model (\code{\link{multinomial}}) may be appropriate.
  Another similar distribution to the Dirichlet is the
  Dirichlet-multinomial (see \code{\link{dirmultinomial}}).

}

\seealso{
  \code{\link{rdiric}},
  \code{\link{dirmultinomial}},
  \code{\link{multinomial}}.

}
\examples{
y = rdiric(n=1000, shape=c(3, 1, 4))
fit = vglm(y ~ 1, dirichlet, trace = TRUE, crit="c")
Coef(fit)
coef(fit, matrix=TRUE)
fitted(fit)[1:2,]
}
\keyword{models}
\keyword{regression}

\name{acat}
\alias{acat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Ordinal Regression with Adjacent Categories Probabilities }
\description{
  Fits an adjacent categories regression model to an ordered
  (preferably) factor response.

}
\usage{
acat(link = "loglink", parallel = FALSE, reverse = FALSE,
     zero = NULL, ynames = FALSE, Thresh = NULL, Trev = reverse,
     Tref = if (Trev) "M" else 1, whitespace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{link}{
  Link function applied to the ratios of the
  adjacent categories probabilities.
  See \code{\link{Links}} for more choices.

  }
  \item{parallel}{
  A logical, or formula specifying which terms have
  equal/unequal coefficients.

  }
  \item{reverse}{
  Logical.
  By default, the linear/additive predictors used are
  \eqn{\eta_j = \log(P[Y=j+1]/P[Y=j])}{eta_j = log(P[Y=j+1]/P[Y=j])}
  for \eqn{j=1,\ldots,M}.
  If \code{reverse} is \code{TRUE} then
  \eqn{\eta_j = \log(P[Y=j]/P[Y=j+1])}{eta_j=log(P[Y=j]/P[Y=j+1])}
  will be used.


  }
  \item{ynames}{
  See \code{\link{multinomial}} for information.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2,\ldots,\eqn{M}\}.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{Thresh, Trev, Tref}{
    See \code{\link{cumulative}} for information.
    These arguments apply to ordinal
    categorical regression models.
    

  }
  \item{whitespace}{
    See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
  In this help file the response \eqn{Y} is assumed to be a
  factor with ordered values \eqn{1,2,\ldots,M+1}, so that
  \eqn{M} is the number of linear/additive predictors
  \eqn{\eta_j}{eta_j}.
  By default, the log link is used because the ratio of
  two probabilities is positive.



  Internally, \code{\link[stats]{deriv3}} is called to
  perform symbolic differentiation and
  consequently this family function will struggle if
  \eqn{M} becomes too large.
  If this occurs, try combining levels so that
  \eqn{M} is effectively reduced.
  One idea is to aggregate levels with the fewest observations
  in them first.



}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{
Agresti, A. (2013).
\emph{Categorical Data Analysis},
3rd ed. Hoboken, NJ, USA: Wiley.
\cr
Tutz, G. (2012).
\emph{Regression for Categorical Data},
Cambridge: Cambridge University Press.
\cr
Yee, T. W. (2010).
The \pkg{VGAM} package for categorical data analysis.
\emph{Journal of Statistical Software},
\bold{32}, 1--34.
\doi{10.18637/jss.v032.i10}.




%\cr
%Simonoff, J. S. (2003).
%\emph{Analyzing Categorical Data},
%New York: Springer-Verlag.
%\cr



% \url{https://www.jstatsoft.org/article/view/v032i10/}.
% \url{https://www.jstatsoft.org/v32/i10/}.
%Documentation accompanying the \pkg{VGAM} package at
%\url{https://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.




}
\author{ Thomas W. Yee }
\note{
  The response should be either a matrix of counts
  (with row sums that are
  all positive), or an ordered factor. In both cases,
  the \code{y} slot returned
  by \code{vglm}/\code{vgam}/\code{rrvglm} is the
  matrix of counts.


  For a nominal (unordered) factor response,
  the multinomial logit model
  (\code{\link{multinomial}}) is more appropriate.


  
  Here is an example of the usage of the \code{parallel} argument.
  If there are covariates \code{x1}, \code{x2} and \code{x3}, then
  \code{parallel = TRUE ~ x1 + x2 -1} and
  \code{parallel = FALSE ~ x3} are equivalent.
  This would constrain the regression coefficients
  for \code{x1} and \code{x2} to be equal; those of the
  intercepts and \code{x3} would be different.


}
\section{Warning }{
  No check is made to verify that the response is ordinal if the
  response is a matrix;
  see \code{\link[base:factor]{ordered}}.


}

\seealso{
    \code{\link{cumulative}},
    \code{\link{cratio}},
    \code{\link{sratio}},
    \code{\link{multinomial}},
  \code{\link{CM.equid}},
  \code{\link{CommonVGAMffArguments}},
    \code{\link{margeff}},
    \code{\link{pneumo}},
    \code{\link{budworm}},
    \code{\link[stats]{deriv3}}.


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
(fit <- vglm(cbind(normal, mild, severe) ~ let, acat, pneumo))
coef(fit, matrix = TRUE)
constraints(fit)
model.matrix(fit)
}
\keyword{models}
\keyword{regression}
%pneumo$let <- log(pneumo$exposure.time)

%(fit2 <- vglm(cbind(normal, mild, severe) ~ let,
%              acat(Thresh = "equid"), pneumo))
%coef(fit2, matrix = TRUE)
%constraints(fit2)
%model.matrix(fit2)





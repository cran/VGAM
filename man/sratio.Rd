\name{sratio}
\alias{sratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Ordinal Regression with Stopping Ratios }
\description{
  Fits a stopping ratio logit/probit/cloglog/cauchit/...
  regression model to an ordered (preferably) factor response.
}
\usage{
sratio(link = "logitlink", parallel = FALSE, reverse = FALSE,
       zero = NULL, thresholds = c("unconstrained", "equidistant"),
       Treverse = reverse, Tref = if (Treverse) "M" else 1,
       whitespace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{link}{
    Link function applied to the \eqn{M}
    stopping ratio probabilities.
    See \code{\link{Links}} for more choices.

  }
  \item{parallel}{
    A logical, or formula specifying which terms have
    equal/unequal coefficients.

  }
  \item{reverse}{
  Logical.
  By default, the stopping ratios used are
  \eqn{\eta_j = logit(P[Y=j|Y \geq j])}{eta_j = logit(P[Y=j|Y>=j])}
  for \eqn{j=1,\dots,M}.
  If \code{reverse} is \code{TRUE}, then
\eqn{\eta_j = logit(P[Y=j+1|Y \leq j+1])}{eta_j = logit(P[Y=j+1|Y<=j+1])}
  will be used.

  }
  \item{zero}{
  Can be an integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2,\ldots,\eqn{M}\}.
  The default value means none are modelled as
  intercept-only terms.
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{thresholds, Treverse, Tref}{
    See \code{\link{cumulative}} for information.
    These arguments apply to ordinal
    categorical regression models.
    

  }
  \item{whitespace}{
  See \code{\link{CommonVGAMffArguments}} for information.

  }
}
\details{
  In this help file the response \eqn{Y} is assumed to be a factor
  with ordered values \eqn{1,2,\dots,M+1}, so that
  \eqn{M} is the number of linear/additive predictors
  \eqn{\eta_j}{eta_j}.

  

  There are a number of definitions for the \emph{continuation ratio}
  in the literature. To make life easier, in the \pkg{VGAM} package,
  we use \emph{continuation} ratios (see \code{\link{cratio}})
  and \emph{stopping} ratios.
  Continuation ratios deal with quantities such as
  \code{logitlink(P[Y>j|Y>=j])}.


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



Boersch-Supan, P. H. (2021).
Modeling insect phenology using ordinal
regression and continuation ratio models.
\emph{ReScience C},
\bold{7.1}, 1--14.
\doi{10.18637/jss.v032.i10}.

%\bold{7.1}(#5), 1--14.



%Simonoff, J. S. (2003).
%\emph{Analyzing Categorical Data},
%New York, USA: Springer-Verlag.


McCullagh, P. and Nelder, J. A. (1989).
\emph{Generalized Linear Models},
2nd ed. London: Chapman & Hall.


Tutz, G. (2012).
\emph{Regression for Categorical Data},
Cambridge: Cambridge University Press.



Yee, T. W. (2010).
The \pkg{VGAM} package for categorical data analysis.
\emph{Journal of Statistical Software},
\bold{32}, 1--34.
\doi{10.18637/jss.v032.i10}.



% \url{https://www.jstatsoft.org/article/view/v032i10/}.
% \url{https://www.jstatsoft.org/v32/i10/}.
%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.



}
\author{ Thomas W. Yee }
\note{
  The response should be either a matrix of counts
  (with row sums that
  are all positive), or a factor. In both cases,
  the \code{y} slot
  returned by \code{vglm}/\code{vgam}/\code{rrvglm}
  is the matrix
  of counts.


  For a nominal (unordered) factor response, the multinomial
  logit model (\code{\link{multinomial}}) is more appropriate.


  Here is an example of the usage of the \code{parallel} argument.
  If there are covariates \code{x1}, \code{x2} and \code{x3}, then
  \code{parallel = TRUE ~ x1 + x2 -1} and
  \code{parallel = FALSE ~ x3} are equivalent. This would constrain
  the regression coefficients for \code{x1} and \code{x2} to be
  equal; those of the intercepts and \code{x3} would be different.


}
\section{Warning }{
  No check is made to verify that the response is ordinal if the
  response is a matrix;
  see \code{\link[base:factor]{ordered}}.




% 20220822:

  Boersch-Supan (2021) considers a sparse data set
  (called \code{\link{budworm}})
  and the numerical problems encountered when
  fitting models such as
  \code{\link{cratio}},
  \code{\link{sratio}},
  \code{\link{cumulative}}.
  Although improvements to links such as
  \code{\link{clogloglink}} have been made,
  currently these family functions have not been
  properly adapted to handle sparse data as well as they could.



}

\seealso{
  \code{\link{cratio}},
  \code{\link{acat}},
  \code{\link{cumulative}},
  \code{\link{multinomial}},
  \code{\link{margeff}},
  \code{\link{pneumo}},
  \code{\link{budworm}},
  \code{\link{logitlink}},
  \code{\link{probitlink}},
  \code{\link{clogloglink}},
  \code{\link{cauchitlink}}.


}

\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
(fit <- vglm(cbind(normal, mild, severe) ~ let,
             sratio(parallel = TRUE), data = pneumo))
coef(fit, matrix = TRUE)
constraints(fit)
predict(fit)
predict(fit, untransform = TRUE)
}
\keyword{models}
\keyword{regression}


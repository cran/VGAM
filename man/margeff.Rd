\name{margeff}
\alias{margeff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Marginal effects for several categorical response models }
\description{
  Marginal effects for the multinomial logit model and
  cumulative logit/probit/... models and
  continuation ratio models and
  stopping ratio models and
  adjacent categories models:
  the derivative of the fitted probabilities with respect to
  each explanatory variable.


}
\usage{
margeff(object, subset = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A \code{\link{vglm}} object,
  with one of the following family functions:
  \code{\link{multinomial}},
  \code{\link{cumulative}},
  \code{\link{cratio}},
  \code{\link{sratio}}
  or
  \code{\link{acat}}.



  }
  \item{subset}{
  Numerical or logical vector, denoting the required observation(s).
  Recycling is used if possible.
  The default means all observations.


  }

  \item{\dots}{
  further arguments passed into the other methods functions.
% e.g., \code{subset}.


  }
}
\details{
  Computes the derivative of the fitted probabilities
  of the categorical response model
  with respect to each explanatory variable.
  Formerly one big function, this function now uses S4
  dispatch to break up the computations.



% 20151215
  The function \code{margeff()} is \emph{not} generic. However, it
  calls the function \code{margeffS4VGAM()} which \emph{is}.
  This is based on the class of the \code{VGAMff} argument, and
  it uses the S4 function \code{\link[methods]{setMethod}} to
  correctly dispatch to the required methods function.
  The inheritance is given by the \code{vfamily} slot of the
  \pkg{VGAM} family function.



}
\value{
  A \eqn{p} by \eqn{M+1} by \eqn{n} array, where \eqn{p} is the number
  of explanatory variables and the (hopefully) nominal response has
  \eqn{M+1} levels, and there are \eqn{n} observations.


  In general, if
  \code{is.numeric(subset)}
  and
  \code{length(subset) == 1} then a
  \eqn{p} by \eqn{M+1} matrix is returned.


}
% \references{ ~put references to the literature/web site here ~ }
\author{ T. W. Yee,
with some help and motivation from Stasha Rmandic.



}

\section{Warning }{
  Care is needed in interpretation, e.g., the change is not
  universally accurate for a unit change in each explanatory
  variable because eventually the `new' probabilities may become
  negative or greater than unity. Also, the `new' probabilities
  will not sum to one.


  This function is not applicable for models with
  data-dependent terms such as \code{\link{bs}} and
  \code{\link{poly}}.
  Also the function should not be applied to models with any terms that
  have generated more than one column of the LM model matrix,
  such as \code{\link{bs}} and \code{\link{poly}}.
  For such try using numerical methods such as finite-differences.
  The \code{formula} in \code{object} should comprise of simple terms
  of the form \code{ ~ x2 + x3 + x4}, etc.



  Some numerical problems may occur if the fitted values are
  close to 0 or 1 for the
  \code{\link{cratio}} and
  \code{\link{sratio}} models.
  Models with offsets may result in an incorrect answer.



}

\note{
  For \code{\link{multinomial}}
  this function should handle any value of \code{refLevel} and also
  any constraint matrices.
  However, it does not currently handle
  the \code{xij} or \code{form2} arguments,
  nor \code{\link{vgam}} objects.



% 20151211; this is now false, so can delete this:
% For \code{\link{multinomial}},
% if \code{subset} is numeric then the function uses a \code{for} loop over
% the observations (slow).
% The default computations use vectorization; this uses more memory than a
% \code{for} loop but is faster.



  Some other limitations are imposed, e.g.,
  for \code{\link{acat}} models
  only a \code{\link{loge}} link is allowed.



}
\seealso{
  \code{\link{multinomial}},
  \code{\link{cumulative}},
  \code{\link{propodds}},
  \code{\link{acat}},
  \code{\link{cratio}},
  \code{\link{sratio}},
  \code{\link{vglm}}.


}

\examples{
# Not a good example for multinomial() because the response is ordinal!!
ii <- 3; hh <- 1/100
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let, multinomial, data = pneumo)
fit <- vglm(cbind(normal, mild, severe) ~ let,
            cumulative(reverse = TRUE,  parallel = TRUE),
            data = pneumo)
fitted(fit)[ii, ]

mynewdata <- with(pneumo, data.frame(let = let[ii] + hh))
(newp <- predict(fit, newdata = mynewdata, type = "response"))

# Compare the difference. Should be the same as hh --> 0.
round(digits = 3, (newp-fitted(fit)[ii, ])/hh)  # Finite-difference approxn
round(digits = 3, margeff(fit, subset = ii)["let",])

# Other examples
round(digits = 3, margeff(fit))
round(digits = 3, margeff(fit, subset = 2)["let",])
round(digits = 3, margeff(fit, subset = c(FALSE, TRUE))["let",,])  # recycling
round(digits = 3, margeff(fit, subset = c(2, 4, 6, 8))["let",,])
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}

% set \code{i=1:n}.
% hh * margeff(fit, i=ii)["let",]
%          cumulative(reverse=TRUE,  parallel=TRUE),
%          cumulative(reverse=FALSE, parallel=TRUE),
%          cumulative(reverse=TRUE,  parallel=FALSE),
%          cumulative(reverse=FALSE, parallel=FALSE),




\name{posbernUC}
\alias{posbernUC}
\alias{dposbern}
\alias{rposbern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Bernoulli Sequence Model }
\description{
  Density, and random generation
  for multiple Bernoulli responses where
  each row in the response matrix has at least one success.

}
\usage{
rposbern(n, nTimePts = 5, pvars = length(xcoeff), xcoeff = c(-2, 1, 2),
         Xmatrix = NULL, cap.effect = 1, is.popn = FALSE,
         link = "logit", earg.link = FALSE)
dposbern(x, prob, prob0 = prob, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  response vector or matrix.
  Should only have 0 and 1 values, at least two columns, and each row
  should have at least one 1.


  }

  \item{nTimePts}{Number of sampling occasions.
  Called \eqn{\tau} in \code{\link{posbernoulli.b}}
  and \code{\link{posbernoulli.t}}.


  }
  \item{n}{number of observations.
  Usually a single positive integer, else the length of the vector is used.
  See argument \code{is.popn}.


  }

  \item{is.popn}{
  Logical.
  If \code{TRUE} then argument \code{n} is the population size
  and what is returned may have substantially less rows than \code{n}.
  That is, if an animal has at least one one in its sequence then
  it is returned, else that animal is not returned because it
  never was captured.
% Put in other words, only animals captured at least once are
% returned in the sample.


  }
  \item{Xmatrix}{
  Optional \bold{X} matrix.
  If given, the \bold{X} matrix is not generated internally.


  }
  \item{cap.effect}{
  Numeric, the capture effect.
  Added to the linear predictor if captured previously.
  A positive or negative value corresponds to
  a trap-happy and trap-shy effect respectively.


  }

% \item{double.ch}{
% Logical.
% If \code{TRUE} then the values of \code{ch0}, \code{ch1}, \ldots are
% 2 or 0, else 1 or 0.
% Setting this argument \code{TRUE} means that a model can be fitted
% with half the capture history in both denominator and numerator
% (this is a compromise of the Huggins (1991) model where the full
% capture history only appears in the numerator).
% }

  \item{pvars}{ Number of other numeric covariates that make up
  the  linear predictor.
  Labelled \code{x1}, \code{x2}, \ldots,
  where the first is an intercept, and the others are
  independent standard \code{\link[stats:Uniform]{runif}} random variates.
  The first \code{pvars} elements of \code{xcoeff} are used.


  }

  \item{xcoeff}{
  The regression coefficients of the linear predictor.
  These correspond to \code{x1}, \code{x2}, \ldots,
  and the first is for the intercept.
  The length of \code{xcoeff} must be at least \code{pvars}.


  }

  \item{link, earg.link}{
  The former is used to generate the probabilities for capture
  at each occasion.
  Other details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{prob, prob0}{
    Matrix of probabilities for the numerator and denominators
    respectively.
    The default does \emph{not} correspond to the
    \eqn{M_b} model since the \eqn{M_b} model has a denominator
    which involves the capture history.


  }
  \item{log}{
  Logical. Return the logarithm of the answer?


  }

}
\details{
  The form of the conditional likelihood is described in
  \code{\link{posbernoulli.b}} and/or
  \code{\link{posbernoulli.t}} and/or
  \code{\link{posbernoulli.tb}}.
  The denominator is equally shared among the elements of
  the matrix \code{x}.


}
\value{
  \code{rposbern} returns a data frame with some attributes.
  The function generates random deviates
  (\eqn{\tau} columns labelled \code{y1}, \code{y2}, \ldots)
  for the response.
  Some indicator columns are also included
  (those starting with \code{ch} are for previous capture history).
  The default setting corresponds to a \eqn{M_{bh}} model that
  has a single trap-happy effect.
  Covariates \code{x1}, \code{x2}, \ldots have the same
  affect on capture/recapture at every sampling occasion
  (see the argument \code{parallel.t} in, e.g.,
  \code{\link{posbernoulli.tb}}).



% and these are useful for the \code{xij} argument.


  The function \code{dposbern} gives the density,


}
%\references{ }
\author{ Thomas W. Yee. }
\note{
  The \code{r}-type function is experimental only and does not follow the
  usual conventions of \code{r}-type R functions.
  It may change a lot in the future.
  The \code{d}-type function is more conventional and is less
  likely to change.


}

\seealso{
  \code{\link{posbernoulli.tb}},
  \code{\link{posbernoulli.b}},
  \code{\link{posbernoulli.t}}.


% \code{\link{huggins91}},


}
\examples{
rposbern(n = 10)
attributes(pdata <- rposbern(n = 100))
M.bh <- vglm(cbind(y1, y2, y3, y4, y5) ~ x2 + x3, posbernoulli.b(I2 = FALSE),
             data = pdata, trace = TRUE)
constraints(M.bh)
summary(M.bh)
}
\keyword{distribution}
\keyword{datagen}


%double.ch = FALSE,
% and those starting with \code{z} are zero.




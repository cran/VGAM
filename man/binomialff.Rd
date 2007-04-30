\name{binomialff}
%\alias{binomial}
\alias{binomialff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Binomial Family Function }
\description{
  Family function for fitting generalized linear models to binomial
  responses, where the dispersion parameter may be known or unknown.

}
\usage{
binomialff(link = "logit", earg = list(),
           dispersion = 1, mv = FALSE, onedpar = !mv,
           parallel = FALSE, zero = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  The notation \eqn{M} is used to denote the number of linear/additive
  predictors.

  \item{link}{ 
  Link function. See \code{\link{Links}} for more choices.

  }
  \item{earg}{ 
  Extra argument optionally used by the link function.
  See \code{\link{Links}} for more information.

  }
  \item{dispersion}{ 
  Dispersion parameter. By default, maximum likelihood is used to
  estimate the model because it is known.  However, the user can specify
  \code{dispersion = 0} to have it estimated, or else specify a known
  positive value (or values if \code{mv} is \code{TRUE}).

  }
  \item{mv}{ 
  Multivariate response? If \code{TRUE}, then the response is interpreted
  as \eqn{M} binary responses, where \eqn{M} is the number of columns
  of the response matrix. In this case, the response matrix should have
  zero/one values only.

  If \code{FALSE} and the response is a (2-column) matrix, then the number
  of successes is given in the first column, and the second column is
  the number of failures.

  }
  \item{onedpar}{ 
  One dispersion parameter? If \code{mv}, then a separate dispersion
  parameter will be computed for each response (column), by default.
  Setting \code{onedpar=TRUE} will pool them so that there is only one
  dispersion parameter to be estimated.

  }
  \item{parallel}{ 
  A logical or formula. Used only if \code{mv} is \code{TRUE}.  This
  argument allows for the parallelism assumption whereby the regression
  coefficients for a variable is constrained to be equal over the \eqn{M}
  linear/additive predictors.

  }
  \item{zero}{ 
  An integer-valued vector specifying which linear/additive predictors
  are modelled as intercepts only.  The values must be from the set
  \{1,2,\ldots,\eqn{M}\}, where \eqn{M} is the number of columns of the
  matrix response.

  }
}
\details{
  If the dispersion parameter is unknown, then the resulting estimate
  is not fully a maximum likelihood estimate (see pp.124--8 of McCullagh
  and Nelder, 1989).

  A dispersion parameter that is less/greater than unity corresponds to
  under-/over-dispersion relative to the binomial model.  Over-dispersion
  is more common in practice.

  Setting \code{mv=TRUE} is necessary when fitting a Quadratic RR-VGLM
  (see \code{\link{cqo}}) because the response is a matrix of \eqn{M}
  columns (e.g., one column per species). Then there will be \eqn{M}
  dispersion parameters (one per column of the response matrix).

  When used with \code{\link{cqo}} and \code{\link{cao}}, it may be
  preferable to use the \code{\link{cloglog}} link.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as
  \code{\link{vglm}},
  \code{\link{vgam}},
  \code{\link{rrvglm}},
  \code{\link{cqo}},
  and \code{\link{cao}}.

}
\references{
 McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.
}
\author{ Thomas W. Yee }

\note{
  If \code{mv} is \code{FALSE} (default), then the response can be of one
  of three formats: a factor (first level taken as success), a vector of
  proportions of success, or a 2-column matrix (first column = successes)
  of counts.  The argument \code{weights} in the modelling function can
  also be specified. In particular, for a general vector of proportions,
  you will need to specify \code{weights} because the number of trials
  is needed.

  If \code{mv} is \code{TRUE}, then the matrix response can only be of
  one format: a matrix of 1's and 0's (1=success).

  The call \code{binomialff(dispersion=0, ...)} is equivalent to
  \code{quasibinomialff(...)}.  The latter was written so that R users
  of \code{quasibinomial()} would only need to add a  ``\code{ff}''
  to the end of the family function name.

  Regardless of whether the dispersion parameter is to be estimated or
  not, its value can be seen from the output from the \code{summary()}
  of the object.

% With the introduction of name spaces for the \pkg{VGAM} package,
% \code{"ff"} can be dropped for this family function.

}
\seealso{
    \code{\link{quasibinomialff}},
    \code{\link{Links}},
    \code{\link{rrvglm}},
    \code{\link{cqo}},
    \code{\link{cao}},
    \code{\link{zibinomial}},
    \code{\link[stats:Binomial]{binomial}}.
}
\section{Warning }{
    With a multivariate response, assigning a known dispersion parameter
    for \emph{each} response is not handled well yet.  Currently, only
    a single known dispersion parameter is handled well.

}

\examples{
quasibinomialff()
quasibinomialff(link="probit")

data(hunua)
fit = vgam(agaaus ~ poly(altitude, 2), binomialff(link=cloglog), hunua)
\dontrun{
attach(hunua)
plot(altitude, agaaus, col="blue", ylab="P(agaaus=1)",
     main="Presence/absence of Agathis australis", las=1)
o = order(altitude)
lines(altitude[o], fitted(fit)[o], col="red", lwd=2)
detach(hunua)
}
}
\keyword{models}
\keyword{regression}


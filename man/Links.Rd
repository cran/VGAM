\name{Links}
\alias{Links}
\alias{TypicalVGAMlinkFunction}
\title{Link functions for VGLM/VGAM/etc. families}
\description{
  The \pkg{VGAM} package provides a number of (parameter) link functions
  which are described in general here. Collectively, they offer the user
  considerable flexibility for modelling data.

}
\usage{
TypicalVGAMlinkFunction(theta, earg=list(), inverse=FALSE,
                        deriv=0, short=TRUE, tag=FALSE)
}
\arguments{
  \item{theta}{
  Numeric or character.
  Actually this can be \eqn{\theta}{theta} (default) or \eqn{\eta}{eta},
  depending on the other arguments.
  If \code{theta} is character then \code{inverse} and 
  \code{deriv} are ignored.

  }
  \item{earg}{
  List.
  Extra argument allowing for additional information, specific to the
  link function.  For example, for \code{\link{logoff}}, this will
  contain the offset value.  The argument \code{earg} is
  always a list with \emph{named} components. See each specific link
  function to find the component names for the list.

  Almost all \pkg{VGAM} family functions with a single link
  function have an argument (often called \code{earg}) which will
  allow parameters to be inputted for that link function.
  For \pkg{VGAM} family functions with more than one link
  function there usually will be an \code{earg}-type argument for
  each link. For example, if there are two links called
  \code{lshape} and \code{lscale} then 
  the \code{earg}-type arguments for these might be called
  \code{eshape} and \code{escale}, say.

  }
  \item{inverse}{
  Logical. If \code{TRUE} the inverse link value
  \eqn{\theta}{theta} is returned, hence the argument
  \code{theta} is really \eqn{\eta}{eta}.

  }
  \item{deriv}{
  Integer. Either 0, 1, or 2 specifying the order of the derivative.

  }
  \item{short, tag}{
  Logical.
  Used for labelling the \code{blurb} slot of a
  \code{\link{vglmff-class}} object.
  Used only if \code{theta} is character, and gives the formula
  for the link in character form.
  If \code{tag=TRUE} then the result contains a little more information.

  }
}
\value{

  Returns one of the link function value or its first or second derivative,
  the inverse link or its first or second derivative,
  or a character description of the link.

  Here are the general details.
  If \code{inverse=FALSE} and \code{deriv=0} (default) then the
  ordinary link
  function \eqn{\eta=g(\theta)}{eta=g(theta)} is returned.
  If \code{inverse=FALSE} and \code{deriv=1} then it is
  \eqn{d\theta / d\eta}{d theta / d eta} 
  \emph{as a function of} \eqn{\theta}{theta}.
  If \code{inverse=FALSE} and \code{deriv=2} then it is
  \eqn{d^2\theta / d\eta^2}{d^2 theta / d eta^2}
  \emph{as a function of} \eqn{\theta}{theta}.

  If \code{inverse=TRUE} and \code{deriv=0} then the inverse
  link function is returned, hence \code{theta} is really
  \eqn{\eta}{eta}.
  If \code{inverse=TRUE} and \code{deriv} is positive then the 
  \emph{reciprocal} of the same link function with
  \code{(theta=theta, earg=earg, inverse=TRUE, deriv=deriv)}
  is returned.

}
\details{
  Almost all \pkg{VGAM} link functions have something similar to
  the argument list as given above.
  That is, there is a matching \code{earg} for each \code{link} argument.
  In this help file
  we have \eqn{\eta=g(\theta)}{eta=g(theta)}
  where \eqn{g} is the link function, \eqn{\theta}{theta} is the parameter
  and \eqn{\eta}{eta} is the linear/additive predictor.





  The following is a brief enumeration of all \pkg{VGAM} link functions.

  For parameters lying between 0 and 1 (e.g., probabilities):
  \code{\link{logit}},
  \code{\link{probit}},
  \code{\link{cloglog}},
  \code{\link{cauchit}},
  \code{\link{loglog}},
  \code{\link{fsqrt}},
  \code{\link{logc}},
  \code{\link{golf}},
  \code{\link{polf}},
  \code{\link{nbolf}}.

  For positive parameters (i.e., greater than 0):
  \code{\link{loge}},
  \code{\link{nloge}},
  \code{\link{powl}}.

  For parameters greater than 1:
  \code{\link{loglog}}.

  For parameters between \eqn{-1} and \eqn{1}:
  \code{\link{fisherz}},
  \code{\link{rhobit}}.

  For parameters between \eqn{A} and \eqn{B}:
  \code{\link{elogit}},
  \code{\link{logoff}} (\eqn{B=\infty}{B=Inf}).

  For unrestricted parameters (i.e., any value):
  \code{\link{identity}},
  \code{\link{nidentity}},
  \code{\link{reciprocal}},
  \code{\link{nreciprocal}}.

% Other links:

}
\references{
    McCullagh, P. and Nelder, J. A. (1989)
    \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.
}

\seealso{
  \code{\link{TypicalVGAMfamilyFunction}},
  \code{\link{vglm}},
  \code{\link{vgam}},
  \code{\link{rrvglm}}.
  \code{\link{cqo}},
  \code{\link{cao}},
  \code{\link{uqo}}.

}
\author{T. W. Yee}
\note{
  \pkg{VGAM} link functions are generally not compatible with other functions
  outside the package.
  In particular, they won't work with
  \code{\link[stats]{glm}} or any other package for fitting GAMs.

  From October 2006 onwards,
  all \pkg{VGAM} family functions will only
  contain one default value for each link argument rather than giving a
  vector of choices.
  For example, rather than 
  \code{binomialff(link=c("logit", "probit", "cloglog",
  "cauchit", "identity"), ...)}
  it is now
  \code{binomialff(link="logit", ...)}
  No checking will be done to see if the user's choice is reasonable.
  This means that the user can write his/her own \pkg{VGAM} link function
  and use it within any \pkg{VGAM} family function.
  Altogether this provides greater flexibility.
  The downside is that the user must specify the \emph{full} name of the
  link function, by either assigning the link argument the full name as
  a character string, or just the name itself.  See the examples below.

}
\examples{
logit("a")
logit("a", short=FALSE)
logit("a", short=FALSE, tag=TRUE)

logoff(1:5, earg=list(offset=1))  # Same as log(1:5 + 1)
powl(1:5, earg=list(power=2))     # Same as (1:5)^2

fit1 = vgam(agaaus ~ altitude, binomialff(link=cloglog), hunua)    # ok
fit2 = vgam(agaaus ~ altitude, binomialff(link="cloglog"), hunua)  # ok

\dontrun{
# This no longer works since "clog" is not a valid VGAM link function:
fit3 = vgam(agaaus ~ altitude, binomialff(link="clog"), hunua) # not ok


# No matter what the link, the estimated var-cov matrix is the same
y = rbeta(n=1000, shape1=exp(0), shape2=exp(1))
fit1 = vglm(y ~ 1, beta.ab(lshape1="identity", lshape2="identity"),
            trace = TRUE, crit="c")
fit2 = vglm(y ~ 1, beta.ab(lshape1=logoff, eshape1=list(offset=1.1),
                           lshape2=logoff, eshape2=list(offset=1.1)),
            trace = TRUE, crit="c")
vcov(fit1, untran=TRUE)
vcov(fit1, untran=TRUE)-vcov(fit2, untran=TRUE)  # Should be all 0s
fit1@misc$earg   # No 'special' parameters
fit2@misc$earg   # Some 'special' parameters are here


par(mfrow=c(2,2))
p = seq(0.01, 0.99, len=200)
x = seq(-4, 4, len=200)
plot(p, logit(p), type="l", col="blue")
plot(x, logit(x, inverse=TRUE), type="l", col="blue")
plot(p, logit(p, deriv=1), type="l", col="blue") # reciprocal!
plot(p, logit(p, deriv=2), type="l", col="blue") # reciprocal!
}
}
\keyword{models}


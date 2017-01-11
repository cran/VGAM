\name{Links}
\alias{Links}
\alias{TypicalVGAMlink}
\title{Link functions for VGLM/VGAM/etc. families}
\description{
  The \pkg{VGAM} package provides a number of (parameter) link functions
  which are described in general here. Collectively, they offer the user
  considerable choice and flexibility for modelling data.

}
\usage{
TypicalVGAMlink(theta, someParameter = 0, bvalue = NULL, inverse = FALSE,
                deriv = 0, short = TRUE, tag = FALSE)
}
\arguments{
  \item{theta}{
  Numeric or character.
  This is usually \eqn{\theta}{theta} (default) but can sometimes
  be \eqn{\eta}{eta},
  depending on the other arguments.
  If \code{theta} is character then \code{inverse} and
  \code{deriv} are ignored.
  The name \code{theta} should always be the name of the first argument.

  }

  \item{someParameter}{
  Some parameter, e.g., an offset.


  }
  \item{bvalue}{
  Boundary value, positive if given.
  If \code{0 < theta} then
  values of \code{theta} which are less than or equal to 0 can be
  replaced by \code{bvalue}
  before computing the link function value.
  Values of \code{theta} which are greater than or equal to 1 can be
  replaced by 1 minus \code{bvalue}
  before computing the link function value.
  The value \code{bvalue = .Machine$double.eps} is sometimes a reasonable
  value, or something slightly higher.


  }


% \item{earg}{
% List.
% Extra argument allowing for additional information, specific to the
% link function.  For example, for \code{\link{logoff}}, this will
% contain the offset value.  The argument \code{earg} is
% always a list with \emph{named} components. See each specific link
% function to find the component names for the list.
%
%
% Almost all \pkg{VGAM} family functions with a single link
% function have an argument (often called \code{earg}) which will
% allow parameters to be inputted for that link function.
% For \pkg{VGAM} family functions with more than one link
% function there usually will be an \code{earg}-type argument for
% each link. For example, if there are two links called
% \code{lshape} and \code{lscale} then
% the \code{earg}-type arguments for these might be called
% \code{eshape} and \code{escale}, say.
%
% }




  \item{inverse}{
    Logical. If \code{TRUE} and \code{deriv = 0}  then
    the inverse link value
  \eqn{\theta}{theta} is returned, hence the argument
  \code{theta} is really \eqn{\eta}{eta}.
  In all other cases, the argument \code{theta} is
  really \eqn{\theta}{theta}.


  }
  \item{deriv}{
  Integer. Either 0, 1, or 2, specifying the order of the derivative.


  }
  \item{short, tag}{
  Logical.
  These are used for labelling the \code{blurb} slot of a
  \code{\link{vglmff-class}} object.
  These arguments are used only if \code{theta} is character,
  and gives the formula for the link in character form.
  If \code{tag = TRUE} then the result is preceeded by a little
  more information.


  }
}
\value{

  Returns one of:
  the link function value or its first or second derivative,
  the inverse link or its first or second derivative,
  or a character description of the link.


  Here are the general details.
  If \code{inverse = FALSE} and \code{deriv = 0} (default) then the
  ordinary link
  function \eqn{\eta = g(\theta)}{eta = g(theta)} is returned.


  If \code{inverse = TRUE} and \code{deriv = 0} then the inverse
  link function value is returned, hence \code{theta} is really
  \eqn{\eta}{eta} (the only occasion this happens).


  If \code{inverse = FALSE} and \code{deriv = 1} then it is
  \eqn{d\eta / d\theta}{d eta / d theta}
  \emph{as a function of} \eqn{\theta}{theta}.
  If \code{inverse = FALSE} and \code{deriv = 2} then it is
  \eqn{d^2\eta / d\theta^2}{d^2 eta / d theta^2}
  \emph{as a function of} \eqn{\theta}{theta}.


  If \code{inverse = TRUE} and \code{deriv = 1} then it is
  \eqn{d\theta / d\eta}{d theta / d eta}
  \emph{as a function of} \eqn{\theta}{theta}.
  If \code{inverse = TRUE} and \code{deriv = 2} then it is
  \eqn{d^2\theta / d\eta^2}{d^2 theta / d eta^2}
  \emph{as a function of} \eqn{\theta}{theta}.



  It is only when \code{deriv = 1} that
  \code{linkfun(theta, deriv = 1, inverse = TRUE)}
  and
  \code{linkfun(theta, deriv = 1, inverse = FALSE)}
    are \emph{reciprocals} of each other.
    In particular,
  \code{linkfun(theta, deriv = 2, inverse = TRUE)}
  and
  \code{linkfun(theta, deriv = 2, inverse = FALSE)}
  are \emph{not} reciprocals of each other in general.



% Prior to 20150711; this was what it was:

% If \code{inverse = FALSE} and \code{deriv = 1} then it is
% \eqn{d\theta / d\eta}{d theta / d eta}
% \emph{as a function of} \eqn{\theta}{theta}.
% If \code{inverse = FALSE} and \code{deriv = 2} then it is
% \eqn{d^2\theta / d\eta^2}{d^2 theta / d eta^2}
% \emph{as a function of} \eqn{\theta}{theta}.


% If \code{inverse = TRUE} and \code{deriv = 0} then the inverse
% link function is returned, hence \code{theta} is really
% \eqn{\eta}{eta}.
% If \code{inverse = TRUE} and \code{deriv} is positive then the
% \emph{reciprocal} of the same link function with
% \code{(theta = theta, someParameter, inverse = TRUE, deriv = deriv)}
% is returned.


}
\details{
  Almost all \pkg{VGAM} link functions have something similar to
  the argument list as given above.
  In this help file we have \eqn{\eta = g(\theta)}{eta = g(theta)}
  where \eqn{g} is the link function, \eqn{\theta}{theta} is the
  parameter and \eqn{\eta}{eta} is the linear/additive predictor.
  The link \eqn{g} must be strictly monotonic and
  twice-differentiable in its range.



% The arguments \code{short} and \code{tag} are used only if
% \code{theta} is character.




% That is, there is a matching \code{earg} for each \code{link} argument.



  The following is a brief enumeration of all \pkg{VGAM} link functions.



  For parameters lying between 0 and 1 (e.g., probabilities):
  \code{\link{logit}},
  \code{\link{probit}},
  \code{\link{cloglog}},
  \code{\link{cauchit}},
  \code{\link{foldsqrt}},
  \code{\link{logc}},
  \code{\link{golf}},
  \code{\link{polf}},
  \code{\link{nbolf}}.


  For positive parameters (i.e., greater than 0):
  \code{\link{loge}},
  \code{\link{negloge}},
  \code{\link{powerlink}}.


  For parameters greater than 1:
  \code{\link{loglog}}.


  For parameters between \eqn{-1} and \eqn{1}:
  \code{\link{fisherz}},
  \code{\link{rhobit}}.


  For parameters between \eqn{A} and \eqn{B}:
  \code{\link{extlogit}},
  \code{\link{logoff}} (\eqn{B = \infty}{B = Inf}).


  For unrestricted parameters (i.e., any value):
  \code{\link{identity}},
  \code{\link{negidentity}},
  \code{\link{reciprocal}},
  \code{\link{negreciprocal}}.


% Other links:

}
\references{
    McCullagh, P. and Nelder, J. A. (1989)
    \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}


\section{Warning }{
  The output of link functions changed at \pkg{VGAM} \code{0.9-9}
  (date was around 2015-07).
  Formerly, \code{linkfun(theta, deriv = 1)} is now
  \code{linkfun(theta, deriv = 1, inverse = TRUE)}, or equivalently,
  \code{1 / linkfun(theta, deriv = 1, inverse = TRUE)}.
  Also, formerly, \code{linkfun(theta, deriv = 2)} was
  \code{1 / linkfun(theta, deriv = 2, inverse = TRUE)}.
  This was a bug.
  Altogether, these are big changes and the user should beware!



  One day in the future, \emph{all} \pkg{VGAM} link functions
  may be renamed to end in the characters \code{"link"}.



}

\seealso{
  \code{\link{TypicalVGAMfamilyFunction}},
  \code{\link{linkfun}},
  \code{\link{vglm}},
  \code{\link{vgam}},
  \code{\link{rrvglm}}.
  \code{\link{cqo}},
  \code{\link{cao}}.


% \code{\link{uqo}}.


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
  \code{binomialff(link = c("logit", "probit", "cloglog",
  "cauchit", "identitylink"), ...)}
  it is now
  \code{binomialff(link = "logit", ...)}.
  No checking will be done to see if the user's choice is reasonable.
  This means that the user can write his/her own \pkg{VGAM} link function
  and use it within any \pkg{VGAM} family function.
  Altogether this provides greater flexibility.
  The downside is that the user must specify the \emph{full} name of the
  link function, by either assigning the link argument the full name as
  a character string, or just the name itself.  See the examples below.


  From August 2012 onwards, a major change in link functions
  occurred.
  Argument \code{esigma} (and the like such as \code{earg})
  used to be in \pkg{VGAM} prior to version 0.9-0 (released
  during the 2nd half of 2012).
  The major change is that arguments such as \code{offset} that used to
  be passed in via those arguments can done directly through
  the link function. For example,
  \code{gev(lshape = "logoff", eshape = list(offset = 0.5))}
  is replaced by
  \code{gev(lshape = logoff(offset = 0.5))}.
  The \code{@misc} slot no longer
  has \code{link} and \code{earg} components, but two other
  components replace these. Functions such as \code{dtheta.deta()},
  \code{d2theta.deta2()}, \code{eta2theta()}, \code{theta2eta()}
  are modified.






}
\examples{
logit("a")
logit("a", short = FALSE)
logit("a", short = FALSE, tag = TRUE)

logoff(1:5, offset = 1)  # Same as log(1:5 + 1)
powerlink(1:5, power = 2)  # Same as (1:5)^2

\dontrun{ # This is old and no longer works:
logoff(1:5, earg = list(offset = 1))
powerlink(1:5, earg = list(power = 2))
}

fit1 <- vgam(agaaus ~ altitude, binomialff(link = "cloglog"), hunua)  # best
fit2 <- vgam(agaaus ~ altitude, binomialff(link =  cloglog ), hunua)  # okay

\dontrun{
# This no longer works since "clog" is not a valid VGAM link function:
fit3 <- vgam(agaaus ~ altitude, binomialff(link = "clog"), hunua)  # not okay


# No matter what the link, the estimated var-cov matrix is the same
y <- rbeta(n = 1000, shape1 = exp(0), shape2 = exp(1))
fit1 <- vglm(y ~ 1, betaR(lshape1 = "identitylink", lshape2 = "identitylink"),
             trace = TRUE, crit = "coef")
fit2 <- vglm(y ~ 1, betaR(lshape1 = logoff(offset = 1.1),
                          lshape2 = logoff(offset = 1.1)), trace = TRUE)
vcov(fit1, untransform = TRUE)
vcov(fit1, untransform = TRUE) - vcov(fit2, untransform = TRUE)  # Should be all 0s
\dontrun{ # This is old:
fit1@misc$earg  # Some 'special' parameters
fit2@misc$earg  # Some 'special' parameters are here
}


par(mfrow = c(2, 2))
p <- seq(0.05, 0.95, len = 200)  # A rather restricted range
x <- seq(-4, 4, len = 200)
plot(p, logit(p), type = "l", col = "blue")
plot(x, logit(x, inverse = TRUE), type = "l", col = "blue")
plot(p, logit(p, deriv = 1), type = "l", col = "blue")  # 1 / (p*(1-p))
plot(p, logit(p, deriv = 2), type = "l", col = "blue")  # (2*p-1)/(p*(1-p))^2
}
}
\keyword{models}


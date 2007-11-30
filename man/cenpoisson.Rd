\name{cenpoisson}
%\alias{cenpoisson}
\alias{cenpoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Censored Poisson Family Function }
\description{
  Family function for a censored
  Poisson response.

}
\usage{
cenpoisson(link = "loge", earg = list(), imu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to the mean.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  Extra argument optionally used by the link function.
  See \code{\link{Links}} for more information.

  }
  \item{imu}{
    Optional initial value.
    See \code{\link{CommonVGAMffArguments}} for more information.

  }
}
\details{
  Often a table of Poisson counts has an entry \emph{J+} meaning
  \eqn{\ge J}.
  This family function is similar to \code{\link{poissonff}} but handles
  such censored data. The input requires \code{\link{Surv}}.
  Only a univariate response is allowed.
  The Newton-Raphson algorithm is used.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as
  \code{\link{vglm}} and
  \code{\link{vgam}}.

}
\references{
  See \pkg{survival} for background.

}

\author{ Thomas W. Yee }

\note{
  The function \code{\link{poissonff}} should be used
  when there are no censored observations.
  Also, \code{NA}s are not permitted with \code{\link{Surv}},
  nor is \code{type="counting"}.

}
\section{Warning }{
  As the response is discrete,
  care is required with \code{\link{Surv}}, especially with
  \code{"interval"} censored data because of the
  \code{(start, end]} format.
  See the examples below.
  The examples have
  \code{y < L} as left censored and
  \code{y >= U} (formatted as \code{U+}) as right censored observations,
  therefore
  \code{L <= y <  U} is for uncensored and/or interval censored observations.
  Consequently the input must be tweaked to conform to the
  \code{(start, end]} format.

}

\seealso{
    \code{\link{Surv}},
    \code{\link{poissonff}},
    \code{\link{Links}}.

}
\examples{
# Example 1: right censored data
set.seed(123)
y = rpois(n <- 100, exp(3))
U = 20
cy = pmin(U, y)
rcensored = y >= U
table(cy)
table(rcensored)
status = ifelse(rcensored, 0, 1)
table(i <- print(Surv(cy, status)))  # Check; U+ means >= U
fit = vglm(Surv(cy, status) ~ 1, cenpoisson, trace=TRUE)
coef(fit, mat=TRUE)
table(print(fit@y))  # Another check; U+ means >= U


# Example 2: left censored data
L = 15
cy = pmax(L, y)
lcensored = y <  L   # Note y < L, not cy == L or y <= L
table(cy)
table(lcensored)
status = ifelse(lcensored, 0, 1)
table(i <- print(Surv(cy, status, type="left")))  # Check
fit = vglm(Surv(cy, status, type="left") ~ 1, cenpoisson, trace=TRUE)
coef(fit, mat=TRUE)


# Example 3: interval censored data
Lvec = rep(L, len=n)
Uvec = rep(U, len=n)
icensored = Lvec <= y & y < Uvec  # Neither lcensored or rcensored
table(icensored)
status = rep(3, n)                    # 3 means interval censored
status = ifelse(rcensored, 0, status) # 0 means right censored
status = ifelse(lcensored, 2, status) # 2 means left  censored
# Have to adjust Lvec and Uvec because of the (start, end] format:
Lvec[icensored] = Lvec[icensored] - 1
Uvec[icensored] = Uvec[icensored] - 1
Lvec[lcensored] = Lvec[lcensored]  # Remains unchanged
Lvec[rcensored] = Uvec[rcensored]  # Remains unchanged
table(i <- print(Surv(Lvec, Uvec, status, type="interval")))  # Check

fit = vglm(Surv(Lvec, Uvec, status, type="interval") ~ 1,
           cenpoisson, trace=TRUE)
coef(fit, mat=TRUE)
table(print(fit@y))  # Another check


# Example 4: Add in some uncensored observations
index = (1:n)[icensored]
index = index[1:min(4,length(index))]
status[index] = 1 # actual or uncensored value
Lvec[index] = y[index]
table(i <- print(Surv(Lvec, Uvec, status, type="interval")))  # Check

fit = vglm(Surv(Lvec, Uvec, status, type="interval") ~ 1,
           cenpoisson, trace=TRUE, crit="c")
coef(fit, mat=TRUE)
table(print(fit@y))  # Another check
}
\keyword{models}
\keyword{regression}


\name{cens.poisson}
%\alias{cens.poisson}
\alias{cens.poisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Censored Poisson Family Function }
\description{
  Family function for a censored
  Poisson response.

}
\usage{
cens.poisson(link = "loge", imu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to the mean;
  see \code{\link{Links}} for more choices.

  }
  \item{imu}{
    Optional initial value;
    see \code{\link{CommonVGAMffArguments}} for more information.

  }
}
\details{
  Often a table of Poisson counts has an entry \emph{J+} meaning
  \eqn{\ge J}{>= J}.
  This family function is similar to \code{\link{poissonff}} but handles
  such censored data. The input requires \code{\link{SurvS4}}.
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
  Also, \code{NA}s are not permitted with \code{\link{SurvS4}},
  nor is \code{type = "counting"}.


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
    \code{\link{SurvS4}},
    \code{\link{poissonff}},
    \code{\link{Links}}.


}
\examples{
# Example 1: right censored data
set.seed(123); U <- 20
cdata <- data.frame(y = rpois(N <- 100, exp(3)))
cdata <- transform(cdata, cy = pmin(U, y),
                          rcensored = (y >= U))
cdata <- transform(cdata, status = ifelse(rcensored, 0, 1))
with(cdata, table(cy))
with(cdata, table(rcensored))
with(cdata, table(ii <- print(SurvS4(cy, status))))  # Check; U+ means >= U
fit <- vglm(SurvS4(cy, status) ~ 1, cens.poisson, data = cdata, trace = TRUE)
coef(fit, matrix = TRUE)
table(print(depvar(fit)))  # Another check; U+ means >= U


# Example 2: left censored data
L <- 15
cdata <- transform(cdata,
                   cY = pmax(L, y),
                   lcensored = y <  L)  # Note y < L, not cY == L or y <= L
cdata <- transform(cdata, status = ifelse(lcensored, 0, 1))
with(cdata, table(cY))
with(cdata, table(lcensored))
with(cdata, table(ii <- print(SurvS4(cY, status, type = "left"))))  # Check
fit <- vglm(SurvS4(cY, status, type = "left") ~ 1, cens.poisson,
            data = cdata, trace = TRUE)
coef(fit, matrix = TRUE)


# Example 3: interval censored data
cdata <- transform(cdata, Lvec = rep(L, len = N),
                          Uvec = rep(U, len = N))
cdata <-
  transform(cdata,
            icensored = Lvec <= y & y < Uvec)  # Not lcensored or rcensored
with(cdata, table(icensored))
cdata <- transform(cdata, status = rep(3, N))  # 3 means interval censored
cdata <-
  transform(cdata,
            status = ifelse(rcensored, 0, status))  # 0 means right censored
cdata <-
  transform(cdata,
            status = ifelse(lcensored, 2, status))  # 2 means left  censored
# Have to adjust Lvec and Uvec because of the (start, end] format:
cdata$Lvec[with(cdata, icensored)] <- cdata$Lvec[with(cdata, icensored)] - 1
cdata$Uvec[with(cdata, icensored)] <- cdata$Uvec[with(cdata, icensored)] - 1
# Unchanged:
cdata$Lvec[with(cdata, lcensored)] <- cdata$Lvec[with(cdata, lcensored)]
cdata$Lvec[with(cdata, rcensored)] <- cdata$Uvec[with(cdata, rcensored)]
with(cdata,
     table(ii <- print(SurvS4(Lvec, Uvec, status, type = "interval"))))  # Check

fit <- vglm(SurvS4(Lvec, Uvec, status, type = "interval") ~ 1,
            cens.poisson, data = cdata, trace = TRUE)
coef(fit, matrix = TRUE)
table(print(depvar(fit)))  # Another check


# Example 4: Add in some uncensored observations
index <- (1:N)[with(cdata, icensored)]
index <- head(index, 4)
cdata$status[index] <- 1  # actual or uncensored value
cdata$Lvec[index] <- cdata$y[index]
with(cdata, table(ii <- print(SurvS4(Lvec, Uvec, status,
                                     type = "interval"))))  # Check

fit <- vglm(SurvS4(Lvec, Uvec, status, type = "interval") ~ 1,
            cens.poisson, data = cdata, trace = TRUE, crit = "c")
coef(fit, matrix = TRUE)
table(print(depvar(fit)))  # Another check
}
\keyword{models}
\keyword{regression}


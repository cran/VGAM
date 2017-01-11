\name{lms.bcg}
\alias{lms.bcg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LMS Quantile Regression with a Box-Cox transformation to a Gamma Distribution }
\description{
  LMS quantile regression with the Box-Cox transformation
  to the gamma distribution.
}
\usage{
lms.bcg(percentiles = c(25, 50, 75), zero = c("lambda", "sigma"),
        llambda = "identitylink", lmu = "identitylink", lsigma = "loge",
        idf.mu = 4, idf.sigma = 2, ilambda = 1, isigma = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{percentiles}{
  A numerical vector containing values between 0 and 100,
  which are the quantiles. They will be returned as `fitted values'.

  }

  \item{zero}{
  See \code{\link{lms.bcn}}.

  }

  \item{llambda, lmu, lsigma}{
  See \code{\link{lms.bcn}}.

  }

  \item{idf.mu, idf.sigma}{
  See \code{\link{lms.bcn}}.

  }
  \item{ilambda, isigma}{
  See \code{\link{lms.bcn}}.

  }

}
\details{
 Given a value of the covariate, this function applies a Box-Cox
 transformation to the response to best obtain a gamma distribution.
 The parameters
 chosen to do this are estimated by maximum likelihood or penalized
 maximum likelihood.
 Similar details can be found at \code{\link{lms.bcn}}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{

Lopatatzidis A. and Green, P. J. (unpublished manuscript)
Semiparametric quantile regression using the gamma distribution.


Yee, T. W. (2004)
Quantile regression via vector generalized additive models.
\emph{Statistics in Medicine}, \bold{23}, 2295--2315.


%Documentation accompanying the \pkg{VGAM} package at
%\url{https://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.


}
\author{ Thomas W. Yee }
\note{
 Similar notes can be found at \code{\link{lms.bcn}}.


}
\section{Warning }{
  This \pkg{VGAM} family function comes with the same
  warnings as \code{\link{lms.bcn}}.
  Also, the expected value of the second derivative with
  respect to lambda may be incorrect (my calculations do
  not agree with the Lopatatzidis and Green manuscript.)


}
\seealso{
  \code{\link{lms.bcn}},
  \code{\link{lms.yjn}},
  \code{\link{qtplot.lmscreg}},
  \code{\link{deplot.lmscreg}},
  \code{\link{cdf.lmscreg}},
  \code{\link{bmi.nz}},
  \code{\link{amlexponential}}.


}

\examples{
# This converges, but deplot(fit) and qtplot(fit) do not work
fit0 <- vglm(BMI ~ sm.bs(age, df = 4), lms.bcg, data = bmi.nz, trace = TRUE)
coef(fit0, matrix = TRUE)
\dontrun{
par(mfrow = c(1, 1))
plotvgam(fit0, se = TRUE)  # Plot mu function (only)
}

# Use a trick: fit0 is used for initial values for fit1.
fit1 <- vgam(BMI ~ s(age, df = c(4, 2)), etastart = predict(fit0),
             lms.bcg(zero = 1), bmi.nz, trace = TRUE)

# Difficult to get a model that converges.
# Here, we prematurely stop iterations because it fails near the solution.
fit2 <- vgam(BMI ~ s(age, df = c(4, 2)), maxit = 4,
             lms.bcg(zero = 1, ilam = 3), bmi.nz, trace = TRUE)
summary(fit1)
head(predict(fit1))
head(fitted(fit1))
head(bmi.nz)
# Person 1 is near the lower quartile of BMI amongst people his age
head(cdf(fit1))

\dontrun{
# Quantile plot
par(bty = "l", mar=c(5, 4, 4, 3) + 0.1, xpd = TRUE)
qtplot(fit1, percentiles=c(5, 50, 90, 99), main = "Quantiles",
       xlim = c(15, 90), las = 1, ylab = "BMI", lwd = 2, lcol = 4)

# Density plot
ygrid <- seq(15, 43, len = 100)  # BMI ranges
par(mfrow = c(1, 1), lwd = 2)
(aa <- deplot(fit1, x0 = 20, y = ygrid, xlab = "BMI", col = "black",
  main = "Density functions at Age = 20 (black), 42 (red) and 55 (blue)"))
aa <- deplot(fit1, x0=42, y=ygrid, add=TRUE, llty=2, col="red")
aa <- deplot(fit1, x0=55, y=ygrid, add=TRUE, llty=4, col="blue", Attach=TRUE)
aa@post$deplot  # Contains density function values
}
}
\keyword{models}
\keyword{regression}


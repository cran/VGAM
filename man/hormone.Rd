\name{hormone}
\alias{hormone}
\docType{data}
\title{
  Hormone Data

}
\description{
  A data set described in Carroll and Ruppert (1988) 
  concerning hormone assay.

%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(hormone)}
\format{
  A data frame with 85 observations on the following 2 variables.
  
  \describe{
    \item{\code{X}}{a numeric vector, suitable as the x-axis in
    a scatter plot.

    }
    \item{\code{Y}}{a numeric vector, suitable as the y-axis in
    a scatter plot.

    }
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~

The data is described in
Carroll and Ruppert (1988).


}
%\source{

% Originally,

%}
\references{

  Carroll, R. J. and Ruppert, D. (1988) 
  \emph{Transformation and Weighting in Regression}.
  New York, USA: Chapman & Hall.


  Yee, T. W. (2012)
  Two-parameter reduced-rank vector generalized linear models.
  \emph{In preparation}.

}

\seealso{
  \code{\link{normal1}},
  \code{\link{rrvglm}}.


}



\examples{
data(hormone)
summary(hormone)

modelI <-rrvglm(Y ~ 1 + X, data = hormone, trace = TRUE,
                normal1(zero = NULL, lsd = "identity", imethod = 2))

# Alternative way to fit modelI
modelI.other <- vglm(Y ~ 1 + X, data = hormone, trace = TRUE,
                     normal1(zero = NULL, lsd = "identity"))

# Inferior to modelI
modelII <- vglm(Y ~ 1 + X, data = hormone, trace = TRUE,
                family = normal1(zero = NULL))

logLik(modelI)
logLik(modelII) # Less than logLik(modelI)


# Reproduce Equations (1)--(3) on p.65 of Carroll and Ruppert (1988)

# Equation (1)
hormone <- transform(hormone, rX = 1 / X)
clist <- list("(Intercept)" = diag(2), X = diag(2), rX = rbind(0, 1))
fit1 <- vglm(Y ~ 1 + X + rX, family = normal1(zero = NULL),
             constraints = clist, data = hormone, trace = TRUE)
coef(fit1, matrix = TRUE)
summary(fit1) # Actually, the intercepts do not seem significant
\dontrun{ plot(Y ~ X, hormone, col = "blue")
lines(fitted(fit1) ~ X, hormone, col = "orange") }

# Equation (2)
fit2 <- rrvglm(Y ~ 1 + X, normal1(zero = NULL), hormone, trace = TRUE)
coef(fit2, matrix = TRUE)
\dontrun{ plot(Y ~ X, hormone, col = "blue")
lines(fitted(fit2) ~ X, hormone, col = "red")
# Add +- 2 SEs
lines(fitted(fit2) + 2 * exp(predict(fit2)[, "log(sd)"]) ~ X,
      hormone, col = "orange")
lines(fitted(fit2) - 2 * exp(predict(fit2)[, "log(sd)"]) ~ X,
      hormone, col = "orange") }

# Equation (3)
# Does not fit well because the loge link for the mean is not good.
fit3 <- rrvglm(Y ~ 1 + X, maxit = 300, data = hormone, trace = TRUE,
               normal1(lmean = "loge", zero = NULL))
coef(fit3, matrix = TRUE)
\dontrun{ plot(Y ~ X, hormone, col = "blue") # Does not look okay.
lines(exp(predict(fit3)[, 1]) ~ X, hormone, col = "red")
# Add +- 2 SEs
lines(fitted(fit3) + 2 * exp(predict(fit3)[, "log(sd)"]) ~ X,
      hormone, col = "orange")
lines(fitted(fit3) - 2 * exp(predict(fit3)[, "log(sd)"]) ~ X,
      hormone, col = "orange") }
}
\keyword{datasets}

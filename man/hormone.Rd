\name{hormone}
\alias{hormone}
\docType{data}
\title{
  Hormone Assay Data

}
\description{
  A hormone assay data set from Carroll and Ruppert (1988).


%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(hormone)}
\format{
  A data frame with 85 observations on the following 2 variables.

  \describe{
    \item{\code{X}}{a numeric vector, suitable as the x-axis in
    a scatter plot.
    The reference method.

    }
    \item{\code{Y}}{a numeric vector, suitable as the y-axis in
    a scatter plot.
    The test method.

    }
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~

The data is given in Table 2.4 of
Carroll and Ruppert (1988), and was downloaded
from \url{http://www.stat.tamu.edu/~carroll}.
The book describes the data as follows.
The data are the results of two assay methods for hormone
data; the scale of the data as presented is not
particularly meaningful, and the original source
of the data refused permission to divulge further
information. As in a similar example of
Leurgans (1980), the old or reference method is
being used to predict the new or test method.
The overall goal is to see whether we can reproduce
the test-method measurements with the reference-method
measurements.
Thus calibration might be of interest for the data.



}
%\source{

% Originally,

%}


\references{

  Carroll, R. J. and Ruppert, D. (1988)
  \emph{Transformation and Weighting in Regression}.
  New York, USA: Chapman & Hall.


  Leurgans, S. (1980)
  Evaluating laboratory measurement techniques.
  \emph{Biostatistics Casebook}.
  Eds.: Miller, R. G. Jr., and Efron, B. and
  Brown, B. W. Jr., and Moses, L.
  New York, USA: Wiley.


  Yee, T. W. (2014)
  Reduced-rank vector generalized linear models with two linear predictors.
  \emph{Computational Statistics and Data Analysis},
  \bold{71}, 889--902.


}

\seealso{
  \code{\link[VGAM]{uninormal}},
  \code{\link[VGAM]{rrvglm}}.


}



\examples{
\dontrun{
data(hormone)
summary(hormone)

modelI <-rrvglm(Y ~ 1 + X, data = hormone, trace = TRUE,
                uninormal(zero = NULL, lsd = "identitylink", imethod = 2))

# Alternative way to fit modelI
modelI.other <- vglm(Y ~ 1 + X, data = hormone, trace = TRUE,
                     uninormal(zero = NULL, lsd = "identitylink"))

# Inferior to modelI
modelII <- vglm(Y ~ 1 + X, data = hormone, trace = TRUE,
                family = uninormal(zero = NULL))

logLik(modelI)
logLik(modelII)  # Less than logLik(modelI)


# Reproduce the top 3 equations on p.65 of Carroll and Ruppert (1988).
# They are called Equations (1)--(3) here.

# Equation (1)
hormone <- transform(hormone, rX = 1 / X)
clist <- list("(Intercept)" = diag(2), X = diag(2), rX = rbind(0, 1))
fit1 <- vglm(Y ~ 1 + X + rX, family = uninormal(zero = NULL),
             constraints = clist, data = hormone, trace = TRUE)
coef(fit1, matrix = TRUE)
summary(fit1)  # Actually, the intercepts do not seem significant
plot(Y ~ X, hormone, col = "blue")
lines(fitted(fit1) ~ X, hormone, col = "orange")

# Equation (2)
fit2 <- rrvglm(Y ~ 1 + X, uninormal(zero = NULL), hormone, trace = TRUE)
coef(fit2, matrix = TRUE)
plot(Y ~ X, hormone, col = "blue")
lines(fitted(fit2) ~ X, hormone, col = "red")
# Add +- 2 SEs
lines(fitted(fit2) + 2 * exp(predict(fit2)[, "loge(sd)"]) ~ X,
      hormone, col = "orange")
lines(fitted(fit2) - 2 * exp(predict(fit2)[, "loge(sd)"]) ~ X,
      hormone, col = "orange")

# Equation (3)
# Does not fit well because the loge link for the mean is not good.
fit3 <- rrvglm(Y ~ 1 + X, maxit = 300, data = hormone, trace = TRUE,
               uninormal(lmean = "loge", zero = NULL))
coef(fit3, matrix = TRUE)
plot(Y ~ X, hormone, col = "blue")  # Does not look okay.
lines(exp(predict(fit3)[, 1]) ~ X, hormone, col = "red")
# Add +- 2 SEs
lines(fitted(fit3) + 2 * exp(predict(fit3)[, "loge(sd)"]) ~ X,
      hormone, col = "orange")
lines(fitted(fit3) - 2 * exp(predict(fit3)[, "loge(sd)"]) ~ X,
      hormone, col = "orange")
}
}
\keyword{datasets}


% from \url{http://www.stat.tamu.edu/~carroll/data/hormone_data.txt}.



\name{huggins91}
\alias{huggins91}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Huggins (1991) Capture-recapture Model Family Function
        (approximation only) }
\description{
  Fits a Huggins (1991) capture-recapture model to a matrix of 0s
  and 1s: animals sampled on several occasions and individual
  animals caught at least once.

}
\usage{
huggins91(link = "logit", earg = list(), parallel = TRUE,
          iprob = NULL, eim.not.oim = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, earg, parallel, iprob}{
  See \code{\link{CommonVGAMffArguments}} for information.
  The \code{parallel} argument should generally be left alone since
  parallelism is assumed by Huggins (1991).

  }
  \item{eim.not.oim}{
  Logical. If \code{TRUE} use the EIM, else the OIM.

  }
}
\details{
  This model operates on a response matrix of 0s and 1s.  Each of
  at least two columns is an occasion where animals are potentially
  captured (e.g., a field trip), and each row is an individual
  animal.  Capture is a 1, else a 0. Each row has at least one
  capture.  It is well-known that animals are affected by capture,
  e.g., trap-shy or trap-happy.  This \pkg{VGAM} family function
  attempts to allow the capture history to be modelled. This
  involves the use of the \code{xij} argument.  Ignoring capture
  history effects would mean \code{\link{posbinomial}} could
  be used by aggregating over the sampling occasions.


  Huggins (1991) suggests a model involving maximizing a
  conditional likelihood.
  The form of this is a numerator divided by a denominator,
  where the true model has part of the linear/additive predictor
  modelling capture history applying to the numerator only,
  so that part is set to zero in the denominator.
  The numerator of the conditional likelihood corresponds
  to a sequence of Bernoulli trials,
  with at least one success,
  for each animal.


  Unfortunately the Huggins model is too difficult to fit in this
  package, and one can only use the \emph{same} linear/additive
  predictor in the numerator as the denominator. Hence this
  \pkg{VGAM} family function does \emph{not} implement the model
  properly.


  The number of linear/additive predictors is twice the number
  of sampling occasions, i.e., \eqn{M = 2T}, say.
  The first two correspond to the first sampling occasion,
  the next two correspond to the second sampling occasion, etc.
  Even-numbered linear/additive predictors should correspond to what
  would happen if no capture had occurred (they belong to
  the denominator.)
  Odd-numbered linear/additive predictors correspond to
  what actually happened (they belong to the numerator.)


  The fitted value for column \eqn{t} is the \eqn{t}th
  numerator probability divided by the denominator.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

Huggins, R. M. (1991)
Some practical aspects of a conditional likelihood
approach to capture experiments.
\emph{Biometrics},
\bold{47}, 725--732.


}
\author{ Thomas W. Yee }

\note{
  The \code{weights} argument of \code{\link{vglm}} need not be
  assigned, and the default is just a matrix of ones.


  This \pkg{VGAM} family function is currently more complicated
  than it needs to be, e.g., it is possible to
  simplify \eqn{M = T}, say.


}

\section{Warning }{
  This \pkg{VGAM} family function is experimental and does not
  work properly because the linear/additive predictor in the
  numerator and denominator must be the same.
  The parameter estimates of the Huggins (1991) model ought to
  be similar (probably in between, in some sense) to two models:
  Model 1 is where the capture history variable is included,
  Model 2 is where the capture history variable is not included.
  See the example below.
  A third model, called Model 3, allows for 'half' the capture
  history to be put in both numerator and denominator. This
  might be thought of as a compromise between Models 1 and 2,
  and may be useful as a crude approximation.


  Under- or over-flow may occur if the data is ill-conditioned.


}
\seealso{ 
    \code{\link{vglm.control}} for \code{xij},
    \code{\link{dhuggins91}},
    \code{\link{rhuggins91}}.
    \code{\link{posbinomial}}.

}

\examples{
set.seed(123); nTimePts = 5
hdata = rhuggins91(n = 1000, nTimePts = nTimePts, pvars = 2)

# The truth: xcoeffs are c(-2, 1, 2) and capeffect = -1

# Model 1 is where capture history information is used
model1  = vglm(cbind(y1, y2, y3, y4, y5) ~ x2 + Chistory,
               huggins91, data = hdata, trace  = TRUE,
               xij = list(Chistory ~ ch0 + zch0 +
                                     ch1 + zch1 + ch2 + zch2 +
                                     ch3 + zch3 + ch4 + zch4 - 1),
               form2 = ~ 1 + x2 + Chistory +
                          ch0 +  ch1 +  ch2 +  ch3 +  ch4 +
                         zch0 + zch1 + zch2 + zch3 + zch4)

coef(model1, matrix = TRUE)  # Biased!!
summary(model1)
head(fitted(model1))
head(model.matrix(model1, type = "vlm"), 21)
head(hdata)

# Model 2 is where no capture history information is used
model2  = vglm(cbind(y1, y2, y3, y4, y5) ~ x2,
               huggins91, data = hdata, trace  = TRUE)
coef(model2, matrix = TRUE)  # Biased!!
summary(model2)

# Model 3 is where half the capture history is used in both
# the numerator and denominator
set.seed(123); nTimePts = 5
hdata2 = rhuggins91(n = 1000, nTimePts = nTimePts, pvars = 2,
                    double.ch = TRUE)
head(hdata2)  # 2s have replaced the 1s in hdata
model3  = vglm(cbind(y1, y2, y3, y4, y5) ~ x2 + Chistory,
               huggins91, data = hdata2, trace  = TRUE,
               xij = list(Chistory ~ ch0 + zch0 +
                                     ch1 + zch1 + ch2 + zch2 +
                                     ch3 + zch3 + ch4 + zch4 - 1),
               form2 = ~ 1 + x2 + Chistory +
                          ch0 +  ch1 +  ch2 +  ch3 +  ch4 +
                         zch0 + zch1 + zch2 + zch3 + zch4)
coef(model3, matrix = TRUE)  # Biased!!
}
\keyword{models}
\keyword{regression}


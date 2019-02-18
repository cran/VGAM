\name{genipoisson}
\alias{genipoisson}
%\alias{genipoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Inflated Poisson Distribution Family Function}
\description{
  Fits a generally-inflated Poisson distribution based on a mixture
  model involving a multinomial distribution
  and a Poisson distribution.



  
}
\usage{
genipoisson(inflate = 0, zero = NULL, llambda = "loglink",
           type.fitted = c("mean", "lambda", "pstr.i", "onempstr.i"),
           imethod = 1, mux.inflate = 0.5, ipstr0 = NULL,
           ilambda = NULL, ishrinkage = 0.95, probs.y = 0.35)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inflate}{
  Vector of inflated values, i.e., nonnegative integers.
  Must have unique values only.
  May be a \code{NULL}, which stands for empty set
  (same as \code{\link{poissonff}}).
  The default is the 0-inflated Poisson.
  See below for recommendations of its use.


% Must be sorted and have unique values only.
  }
  \item{llambda}{
  See \code{\link{Links}}
  for more choices and information.

  }

  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}} for information.
  The choice \code{"pstr.i"} is the sum of the structural probabilities
  and \code{"onempstr.i"} is its complement.
  

  }

 
  \item{mux.inflate}{
    Numeric, between 0 and 1.
    For the initial values of the multinomial probabilities,
    the sample proportions are multiplied by
    this multiplication factor to estimate the amount of inflation.
    The default is to assume that half of all the inflated
    values are due to inflation, and the other half from the
    Poisson distribution.
    A value close to 1 is recommended for severe amounts of inflation.
    


  }
  \item{imethod, ipstr0, ilambda}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}} for information.
  Setting \code{zero = "pstr"} will model the multinomial
  probabilities as simple as possible (intercept-only), hence
  should be more numerically stable than the default, and
  this is recommended for many analyses especially when there
  are many explanatory variables.


  }
}
\details{
  The generally-inflated Poisson distribution is a mixture
  of an ordinary Poisson
  distribution and the probability of certain values
  (\code{inflate} argument) being inflated by a multinomial distribution.
  The latter is modelled by a multinomial logit model.
  Thus each inflated value has two sources: either from the Poisson
  distribution or the multinomial distribution.
  The 0-inflated Poisson distribution is a special case
  and is the default.


  This family function can be easily misused.
  The data should exhibit \emph{strong} evidence of inflation in
  certain values before the \code{inflate} argument is assigned these
  values.
  Deflation or very little inflation will result in numerical problems.
  The \code{inflate} argument should not be assigned a vector
  having too many values.
  Probably \code{\link{genapoisson}} is recommended in general
  instead as it handles deflation and inflation.

  

  This function implements Fisher scoring and currently
  does not handle multiple responses.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  
  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted},
  returns the mean (default).
  The choice \code{type.fitted = "pstr.i"} returns the sum of all
  the inflated probabilities.



% returns the mean \eqn{\mu}{mu} (default).





}
%\references{
%}



%20111123; this has been fixed up with proper FS using EIM.
%\section{Warning }{
%  Inference obtained from \code{summary.vglm}
%  and \code{summary.vgam} may or may not be correct.
%  In particular, the p-values, standard errors and degrees of
%  freedom may need adjustment. Use simulation on artificial
%  data to check that these are reasonable.
%
%
%}


\author{ T. W. Yee and Chenchen Ma}
\note{
  This family function does not yet have the robustness of
  \code{\link{multinomial}} when computing the working weight matrices.
  Several improvements are needed, e.g.,
  better labelling and initial values and handling multiple responses.
  And \code{\link{summaryvglm}} does not work yet on these fitted
  models.



% yettodo: see lines just above.


  

% This family function effectively combines \code{\link{poissonff}}
% and \code{\link{multinomial}} into one family function.
% This family function can handle multiple responses,
% e.g., more than one species.


}

\seealso{
  \code{\link{Genipois}},
  \code{\link{genapoisson}},
  \code{\link{gentpoisson}},
  \code{\link{zipoisson}},
  \code{\link{zipoissonff}},
  \code{\link{multinomial}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.


% \code{\link{rgenipois}},



}

\examples{
ivec <- c(0, 3)  # Inflate these values
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, x3 = runif(nn))
gdata <- transform(gdata, lambda1   = exp(1 + 2 * x2 + 0.5 * x3),
                          lambda2   = exp(1 - 1 * x2 + 0.5 * x3),
                          lambda3   = exp(1))
gdata <- transform(gdata,
  y1 = rgenipois(nn, lambda1, pstr.i = c(0.2, 0.3), byrow = TRUE,
                 inflate = ivec),
  y2 = rgenipois(nn, lambda2, pstr.i = c(0.2, 0.3), byrow = TRUE,
                 inflate = ivec),
  y3 = rgenipois(nn, lambda3, pstr.i = c(0.2, 0.3), byrow = TRUE,
                 inflate = ivec))
genipoisson(inflate = ivec)
with(gdata, table(y1))
with(gdata, table(y2))
with(gdata, table(y3))

fit1 <- vglm(y1 ~ x2 + x3, genipoisson(inflate = ivec),
             crit = "coef", trace = TRUE, data = gdata)
head(fitted(fit1))
head(predict(fit1))
coef(fit1, matrix = TRUE)
# summary(fit1)  # Currently does not work!!

# Another example ------------------------------
fit3 <- vglm(y3 ~ 1, genipoisson(inflate = ivec), gdata, trace = TRUE)
coef(fit3, matrix = TRUE)
}
\keyword{models}
\keyword{regression}


%genipoisson(lpobs0 = "logitlink", llambda = "loglink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%genipoissonff(llambda = "loglink", lonempobs0 = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")




\name{genapoisson}
\alias{genapoisson}
%\alias{genapoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered Poisson Distribution Family Function}
\description{
  Fits a generally-altered Poisson distribution based on a conditional
  model involving a multinomial distribution
  and a generally-truncated Poisson distribution.



  
}
\usage{
genapoisson(alter = 0, zero = NULL, llambda = "loglink", 
            type.fitted = c("mean", "lambda", "pobs.a", "onempobs.a"),
            imethod = 1, ilambda = NULL, ishrinkage = 0.95, probs.y = 0.35)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{alter}{
  Vector of altered values, i.e., nonnegative integers.
  Must have unique values only.
  May be a \code{NULL}, which stands for empty set
  (same as \code{\link{poissonff}}).
  The default is the 0-altered Poisson.


% Must be sorted and have unique values only.
  }
  \item{llambda}{
  See \code{\link{Links}}
  for more choices and information.

  }

  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}} for information.
  The choice \code{"pobs.a"} is the probability of an altered value,
  and \code{"onempobs.a"} is its complement.
  

  }

 
  \item{imethod, ilambda}{
  See \code{\link{CommonVGAMffArguments}} for information.


% ipobs0, 



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}} for information.
  Setting \code{zero = "pobs"} will model the multinomial
  probabilities as simple as possible (intercept-only), hence
  should be more numerically stable than the default, and
  this is recommended for many analyses especially when there
  are many explanatory variables..


  }
}
\details{
  The generally-altered Poisson distribution is the ordinary Poisson
  distribution with the probability of certain values
  (\code{alter} argument) being modelled using a multinomial
  logit model.
  The other values are modelled using a
  \emph{generally-truncated Poisson} distribution.
  The 0-altered Poisson distribution is a special case
  and is the default
  (it is called a \emph{hurdle} model by some people).



  This function implements Fisher scoring and currently
  does not handle multiple responses.
  Compared to \code{\link{genipoisson}}
  this family function handles deflation and inflation,
  therefore handles a wider range of data.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  
  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted},
  returns the mean \eqn{\mu}{mu} (default).
  The choice \code{type.fitted = "pobs.a"} returns the sum of all
  the altered probabilities.



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


  

  This family function effectively combines \code{\link{gentpoisson}}
  and \code{\link{multinomial}} into one family function.
% This family function can handle multiple responses,
% e.g., more than one species.


}

\seealso{
  \code{\link{Genapois}},
  \code{\link{gentpoisson}},
  \code{\link{genipoisson}},
  \code{\link{zapoisson}},
  \code{\link{multinomial}},
  \code{\link{rpospois}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.





}

\examples{
avec <- c(0, 3)  # Alter these values
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, x3 = runif(nn))
gdata <- transform(gdata, lambda1   = exp(1 + 2 * x2 + 0.5 * x3),
                          lambda2   = exp(1 - 1 * x2 + 0.5 * x3),
                          lambda3   = exp(1))
gdata <- transform(gdata,
  y1 = rgenapois(nn, lambda = lambda1, pobs.a = c(0.2, 0.3), byrow = TRUE,
                 alter = avec),
  y2 = rgenapois(nn, lambda = lambda2, pobs.a = c(0.2, 0.3), byrow = TRUE,
                 alter = avec),
  y3 = rgenapois(nn, lambda = lambda3, pobs.a = c(0.2, 0.3), byrow = TRUE,
                 alter = avec))
genapoisson(alter = avec)
with(gdata, table(y1))
with(gdata, table(y2))
with(gdata, table(y3))

fit1 <- vglm(y1 ~ x2 + x3, genapoisson(alter = avec),
             crit = "coef", trace = TRUE, data = gdata)
head(fitted(fit1))
head(predict(fit1))
coef(fit1, matrix = TRUE)
# summary(fit1)  # Currently does not work!!

# Another example ------------------------------
fit3 <- vglm(y3 ~ 1, genapoisson(alter = avec), gdata, trace = TRUE)
coef(fit3, matrix = TRUE)
}
\keyword{models}
\keyword{regression}


%genapoisson(lpobs0 = "logitlink", llambda = "loglink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%genapoissonff(llambda = "loglink", lonempobs0 = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")




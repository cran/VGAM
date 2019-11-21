\name{gabinomial.mlm}
\alias{gabinomial.mlm}
%\alias{gabinomial.mlmff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered Binomial Regression Family Function
        (multinomial logit model based; GA-binom-MLM)

}
\description{
  Fits a generally-altered binomial regression
 (using a multinomial logit model for the altered values). 



  
}
\usage{
gabinomial.mlm(alter = 0, zero = NULL, lprob = "logitlink", 
    type.fitted = c("mean", "prob", "pobs.a", "Pobs.a"),
    imethod = 1, iprob = NULL)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
% ishrinkage = 0.95,, probs.y = 0.35
\arguments{
  \item{alter}{
  Vector of altered values, i.e., nonnegative integers.
  Must have unique values only.
  May be a \code{NULL}, which stands for empty set
  (same as \code{\link{binomialff}}).
  The default is the 0-altered binomial


% Must be sorted and have unique values only.
  }
  \item{lprob}{
  See \code{\link{Links}}
  for more choices and information.

  }

  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}} for information.
  The choice \code{"pobs.a"} is the probability of an altered value.
  See \code{\link{gatpoisson.mlm}} for more details.

% and \code{"onempobs.a"} is its complement.


  }

 
  \item{imethod, iprob}{
  See \code{\link{CommonVGAMffArguments}} for information.


% ipobs0, 



  }
% \item{probs.y}{
% See \code{\link{CommonVGAMffArguments}} for information.
% }



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
  The generally-altered binomial distribution is the ordinary binomial
  distribution with the probability of certain values
  (\code{alter} argument) being modelled using a multinomial
  logit model.
  The other values are modelled using a
  \emph{generally-truncated binomial} distribution.
  The 0-altered binomial distribution is a special case
  and is the default
  (it is called a \emph{hurdle} model by some people).



  This function implements Fisher scoring and currently
  does not handle multiple responses.
  Compared to \code{gibinomial} (not complete yet)
  this family function handles deflation and inflation,
  therefore handles a wider range of data.


% Compared to \code{\link{gibinomial}}



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
  See \code{\link{gatpoisson.mlm}} for more details.



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



  The default value of \code{zero} is not the best for most
  data sets; so assign it a value (see above).



% And \code{\link{summaryvglm}} does not work yet on these fitted
% models.



% yettodo: see lines just above.


  

  This family function effectively combines \code{\link{gtbinomial}}
  and \code{\link{multinomial}} into one family function.
% This family function can handle multiple responses,
% e.g., more than one species.


}

\seealso{
  \code{\link{Gaitbinom.mlm}},
  \code{\link{gtbinomial}},
  \code{\link{zabinomial}},
  \code{\link{multinomial}},
  \code{\link{rposbinom}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.



% \code{\link{gibinomial}},


}

\examples{
avec <- c(0, 5)  # Alter these values
Size <- 10
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, x3 = runif(nn))
gdata <- transform(gdata,
           prob1 = logitlink(1 + 2 * x2 + .5 * x3, inverse = TRUE),
           prob2 = logitlink(1 - 1 * x2 + .5 * x3, inverse = TRUE),
           prob3 = logitlink(1, inverse = TRUE))
gdata <- transform(gdata,
  y1 = rgaitbinom.mlm(nn, prob = prob1, Size, pobs.a = c(0.2, 0.3),
                      byrow = TRUE, alter = avec),
  y2 = rgaitbinom.mlm(nn, prob = prob2, Size, pobs.a = c(0.2, 0.3),
                      byrow = TRUE, alter = avec),
  y3 = rgaitbinom.mlm(nn, prob = prob3, Size, pobs.a = c(0.2, 0.3),
                      byrow = TRUE, alter = avec))
gabinomial.mlm(alter = avec)
with(gdata, table(y1))
with(gdata, table(y2))
with(gdata, table(y3))

fit1 <- vglm(cbind(y1, Size-y1) ~ x2 + x3, gabinomial.mlm(alter = avec),
             crit = "coef", trace = TRUE, data = gdata)
head(fitted(fit1))
head(predict(fit1))
coef(fit1, matrix = TRUE)
summary(fit1)

# Another example ------------------------------
fit3 <- vglm(cbind(y3, Size-y3) ~ 1, gabinomial.mlm(alter = avec),
             gdata, trace = TRUE)
coef(fit3, matrix = TRUE)
}
\keyword{models}
\keyword{regression}


%gabinomial.mlm(lpobs0 = "logitlink", lprob = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%gabinomial.mlmff(lprob = "logitlink", lonempobs0 = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")




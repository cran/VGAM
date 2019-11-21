\name{gatnbinomial.mlm}
\alias{gatnbinomial.mlm}
%\alias{ganbinomialff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered and -Truncated Negative Binomial
  Regression Family Function (GAT--NB--MLM Variant)

}
\description{
  Fits a generally-altered and -truncated negative binomial
  regression
  (using a multinomial logit model for the altered values).



% The truncation may include values in the upper tail.


%  based on a conditional
%  model involving a multinomial distribution
%  and a generally-truncated negative binomial distribution.


}
\usage{
gatnbinomial.mlm(alter = NULL, truncate = NULL, zero = "size",
   lmunb = "loglink", lsize = "loglink", type.fitted = c("mean",
   "pobs.a", "Pobs.a", "prob.a", "prob.t"),
   imethod = 1, imunb = NULL, isize = exp(1), ishrinkage = 0.95,
   probs.y = 0.35, cutoff.prob = 0.999, eps.trig = 1e-7,
   max.chunk.MB = 30)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{alter, truncate}{
    Same as \code{\link{gatpoisson.mlm}}.



% Must be sorted and have unique values only.
  }
  \item{lmunb, lsize}{
  See \code{\link{Links}}
  for more choices and information.
  Similar to \code{\link{negbinomial}}.

  }

  \item{type.fitted, zero}{
    See \code{\link{CommonVGAMffArguments}} and
    \code{\link{gatpoisson.mlm}}
    for information.


  }

 
  \item{imethod, imunb, isize}{
  See \code{\link{CommonVGAMffArguments}} for information.


% ipobs0, 



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{cutoff.prob, eps.trig}{
  See \code{\link{negbinomial}} for information.



  }
  \item{max.chunk.MB}{
  See \code{\link{negbinomial}} for information.



  }
}
\details{
  The generally-truncated (GT) negative binomial distribution is an
  ordinary negative binomial
  distribution with the probability of certain values
  (given by the \code{truncate} argument) being zero.
  Thus the other probabilities are scaled up.
  The (0-truncated) positive-negative binomial distribution is a
  special case
  (\code{\link{posnegbinomial}}).


  
  The generally-altered (GA) negative binomial distribution is an
  ordinary negative binomial
  distribution with the probability of certain values
  (given by the \code{alter} argument) being modelled using a multinomial
  logit model
  (see \code{\link{multinomial}}).
  The 0-altered negative binomial distribution is a special case
  (\code{\link{zanegbinomial}})
  and it is called a \emph{hurdle} model by some people.

  
  
% The other values are modelled using a
% \emph{generally-truncated negative binomial} distribution.


This function can fit both the GA and GT models simultaneously,
called the GAT-NB-MLM.
It might be applied to heaped data.
That is, each special value can be altered or truncated but
not both.
The default settings make this family function
equivalent to \code{\link{negbinomial}}.



  This function implements Fisher scoring and currently
  does not handle multiple responses for GT or the GA model.
  Compared to what could be \code{ginbinomial.mlm}
  this family function handles deflation and inflation,
  therefore handles a wider range of data.
  For further details please see
  \code{\link{Gaitnbinom.mlm}}.
  

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
  See \code{\link{gatpoisson.mlm}} for more information.


}
%\references{
%}



\section{Warning }{
  Due to its flexibility, it is easy to misuse this function;
  the \code{truncate} vector should ideally be not very long and
  have values that can be justified by the application on hand.
  Likewise, the \code{alter} vector should be short too, and
  each value should have good justification for being included.
  Adding unnecessary values to these two arguments willy-nilly
  is a recipe for disaster.



  Regarding truncation,
  under- or over-flow may occur if the data is ill-conditioned.
  The response is checked to see that no values equal
  any values of the \code{truncate} vector.



  Compared to \code{\link{gatpoisson.mlm}}
  this family function is even more difficult to fit because
  it is more flexible and involves approximating an infinite
  series when computing the expected information matrix.

    
% See \code{\link{gatpoisson.mlm}} for more information;




}


\author{ T. W. Yee }
\note{
  See \code{\link{gatpoisson.mlm}} for more information
  about robustness.


% yettodo: see lines just above.


  

This family function effectively combines what could have been
\code{ganbinomial.mlm()} and
\code{gtnbinomial.mlm()} together.
The former would have effectively included
\code{\link{multinomial}} inside it.
% This family function can handle multiple responses,
% e.g., more than one species.
It is possible in the near future that
  \code{\link{zanegbinomial}},
  \code{\link{zanegbinomialff}},
  \code{\link{Zanegbin}},
  \code{\link{posnegbinomial}},
  \code{\link{Posnegbin}}
be moved into \pkg{VGAMdata}.
Likewise for
  \code{\link{zapoisson}},
  \code{\link{zapoissonff}},
  \code{\link{Zapois}}.



An argument \code{max.support} is not implemented because
the upper RHS tail of the NB distribution does not admit
practical formulas that can be computed.
However, instead  of using something like
\code{max.support = A}
one could try something like
\code{truncate = (A+1):(A+large)}
for some positive integer \code{large} that is large enough,
but within reason.
For example, instead of
\code{max.support = 20}
one could try
\code{truncate = 21:30}
or
\code{truncate = 21:40}.


}

\seealso{
  \code{\link{Gaitnbinom.mlm}},
  \code{\link{gatpoisson.mlm}},
  \code{\link{gatnbinomial.mix}},
  \code{\link{zanegbinomial}},
  \code{\link{multinomial}},
  \code{\link{Posnegbin}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.


% \code{\link{gtnbinomial}},



}

\examples{
\dontrun{
avec <- c(10, 20, 30)    # Alter these values
tvec <- 0  # Truncate this value
pobs.a <- logitlink(-(2:4), inverse = TRUE)  # Between 0.02 and 0.12
size1 <- exp(1)
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, lambda1 = exp(2 + 0.5 * x2))
gdata <- transform(gdata,
  y1 = rgaitnbinom.mlm(nn, size1, mu = lambda1, pobs.a = pobs.a,
                       truncate = tvec, byrow = TRUE, alter = avec))
gatnbinomial.mlm(alter = avec)
(ty1 <- with(gdata, table(y1)))
propn1 <- c(ty1) / sum(ty1)
plot(as.numeric(names(ty1)), propn1, las = 1, xlab = "y",
     yaxs = "i", ylim = c(0, max(propn1) * 1.1), main = "Heaped data",
     ylab = "Proportion", lwd = 3, type = "h", col = "blue")
fit1 <- vglm(y1 ~ x2, trace = TRUE, data = gdata, crit = "coef",
             gatnbinomial.mlm(alter = avec, truncate = tvec,
                              zero = c("size", "pobs")))
head(fitted(fit1))
head(predict(fit1))
coef(fit1, matrix = TRUE)
summary(fit1) }
}
\keyword{models}
\keyword{regression}


%ganbinomial(lpobs0 = "logitlink", llambda = "loglink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%ganbinomialff(llambda = "loglink", lonempobs0 = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")


%                cutoff.prob = 0.999, eps.trig = 1e-7,



\name{gatnbinomial.mix}
\alias{gatnbinomial.mix}
%\alias{gapoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered and -Truncated Negative Binomial
  Regression Family Function (GAT--NB--NB Mixture Variant)

}
\description{
  Fits a generally-altered and -truncated negative binomial regression
  (mixture of 2 NBs on differing support).
  The truncation may include values in the upper tail.


}
\usage{
gatnbinomial.mix(alter = NULL, truncate = NULL,
    zero = c("pobs.a", "size"), parallel = FALSE,
    lpobs.a = "logitlink", lmunb.p = "loglink", lsize.p = "loglink",
    lmunb.a = "loglink", lsize.a = "loglink",
    type.fitted = c("mean", "pobs.a", "Pobs.a", "prob.a", "prob.t"),
    imethod = 1, imunb.p = NULL, isize.p = NULL,
    imunb.a = imunb.p, isize.a = isize.p, ishrinkage = 0.95,
    probs.y = 0.35, cutoff.prob = 0.999, eps.trig = 1e-7,
    max.chunk.MB = 30)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{alter, truncate}{
    Very similar to \code{\link{gatpoisson.mix}}.
\emph{Note:} \code{alter}
\emph{must be assigned a vector of length 3 or more for this function to work}.
  In general, the regression tends to be more stable when
  \code{length(alter)} is not too small and the values of
  \code{alter} are spread out over the entire support.

  

% Must be sorted and have unique values only.
  }
  \item{lpobs.a, lmunb.p, lmunb.a}{
    Link functions;
    the \code{.p} and \code{.a} refer to the parent
    and altered distributions respectively.
    See \code{\link{Links}} for more choices and information.

  }
  \item{lsize.p, lsize.a}{
    Same as above.

  }

  \item{parallel, type.fitted}{
    Very similar to \code{\link{gatpoisson.mix}};
    the former concerns: constrain the mean and size parameters
    to be equal?
    Setting \code{parallel = TRUE} is probably a good idea
    for many data sets, especially when \code{length(alter)} is low
    and the values of \code{alter} are not spread out.


  }

 
  \item{imethod, imunb.p, isize.p}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }


  \item{imunb.a, isize.a}{
  See \code{\link{CommonVGAMffArguments}} for information.

  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{zero}{
    See \code{\link{CommonVGAMffArguments}}
    and \code{\link{gatpoisson.mix}}
    for information.


  }
  \item{cutoff.prob, eps.trig, max.chunk.MB}{
    See \code{\link{negbinomial}}
    for information.


  }
}
\details{
  This distribution,
  also known as the GAT--NB--NB mixture,
  is an alternative to
  \code{\link{gatpoisson.mix}} and allows for overdispersion
  relative to the Poisson distribution.
  Because the \code{size} parameters are estimated this
  family function is more difficult to fit and is numerically
  fraught, in comparison.


  The distribution is also a more structured model compared
  to \code{\link{gatnbinomial.mlm}} because the
  (outer) distribution of the altered values is also an NB too.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  
  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted},
  are similar to \code{\link{gatnbinomial.mlm}}.




}
\references{
Yee, T. W. and Ma, C. C. (2019)
Generally-altered, -inflated and -truncated count regression,
with application to heaped and seeped data.
\emph{In preparation}.


%, \bold{3}, 15--41.


}



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


\section{Warning }{
  The same caution needed
  for \code{\link{gatnbinomial.mlm}}
  and \code{\link{gatpoisson.mix}}
  applies here,
  however this function
  is a bit more parametric (structured) in comparison to
  \code{.mlm} functions,
  especially when \code{parallel = TRUE}.

    
}








\author{ T. W. Yee}
\note{
  Convergence is rather slow because of an infinite series that is
  approximated with finite terms.
  This family function is quite expensive.
  See also \code{\link{gatpoisson.mix}}.
  This family function offers potential for regressing
  heaped count data.

  

% This is not true, as 'alter' needs a 2-vector at least:
%  This family function effectively
%  renders the following functions as obsolete
%  (or rather, they are just special cases):
%  \code{\link{pospoisson}},
%  \code{\link{zapoisson}}.



}

\seealso{
  \code{\link{Gaitnbinom.mix}},
  \code{\link{gatnbinomial.mlm}},
  \code{\link{gatpoisson.mix}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.



% \code{\link{rposnegbin}},
% \code{\link{multinomial}},
% \code{\link{zapoisson}},
% \code{\link{gatnbinomial.mlm}},
% \code{\link{gipoisson}},



}

\examples{
avec <- seq(5, 25, by = 5)  # Alter these values
tvec <- c(4, 7, 22)   # Truncate these values
pobs.a <- logitlink(-0.5, inverse = TRUE)  # About 0.5
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, munb.p = exp(2 + 0 * x2), size.p = exp(1))
gdata <- transform(gdata,
  y1 = rgaitnbinom.mix(nn, size.p = size.p, munb.p = munb.p,
         pobs.a = pobs.a, truncate = tvec, alter = avec))
gatnbinomial.mix(alter = avec)
(ty1 <- with(gdata, table(y1)))
\dontrun{
plot(as.numeric(names(ty1)), c(ty1) / sum(ty1), xlab = "y",
     ylab = "Proportion", las = 1, type = "h", col = "blue")
fit1 <- vglm(y1 ~ 1, trace = TRUE, data = gdata,
             gatnbinomial.mix(alter = avec, truncate = tvec,
                              parallel = TRUE))
head(fitted(fit1, type.fitted = "Pobs.a"))
head(predict(fit1))
coef(fit1, matrix = TRUE)
summary(fit1) }
}
\keyword{models}
\keyword{regression}


%gapoisson(lpobs0 = "logitlink", llambda = "loglink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%gapoissonff(llambda = "loglink", lonempobs0 = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")




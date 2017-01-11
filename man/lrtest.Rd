\name{lrtest}
\alias{lrtest}
\alias{lrtest_vglm}
%\alias{update_formula}
%\alias{update_default}
\title{Likelihood Ratio Test of Nested Models}

\description{
 \code{lrtest} is a generic function for carrying out likelihood ratio tests.
 The default method can be employed for comparing nested VGLMs
 (see details below).

}

\usage{
 lrtest(object, \dots)

 lrtest_vglm(object, \dots, no.warning = FALSE, name = NULL)

}
%\method{lrtest}{default}(object, \dots, name = NULL)

%\method{lrtest}{formula}(object, \dots, data = list())



\arguments{
  \item{object}{
  a \code{\link{vglm}} object.
  See below for details.


  }
  \item{\dots}{
  further object specifications passed to methods.
  See below for details.


  }
  \item{no.warning}{
  logical; if \code{TRUE} then no warning is issued.
  For example, setting \code{TRUE} might be a good idea when testing
  for linearity of a variable for a \code{"pvgam"} object.



  }
  \item{name}{
  a function for extracting a suitable name/description from
  a fitted model object.
  By default the name is queried by calling \code{\link{formula}}.


  }



% \item{data}{
% a data frame containing the variables in the model.
%
% }


}

\details{
  \code{lrtest} is intended to be a generic function for
  comparisons of models via asymptotic likelihood ratio
  tests. The default method consecutively compares the
  fitted model object \code{object} with the models passed
  in \code{\dots}. Instead of passing the fitted model
  objects in \code{\dots}, several other specifications
  are possible. The updating mechanism is the same as for
  \code{\link[lmtest]{waldtest}}: the models in \code{\dots}
  can be specified as integers, characters (both for terms
  that should be eliminated from the previous model),
  update formulas or fitted model objects. Except for
  the last case, the existence of an \code{\link[stats]{update}}
  method is assumed. See \code{\link[lmtest]{waldtest}} for details.


  Subsequently, an asymptotic likelihood ratio test for each
  two consecutive models is carried out: Twice the difference
  in log-likelihoods (as derived by the \code{\link[stats]{logLik}}
  methods) is compared with a Chi-squared distribution.


% The \code{"formula"} method fits a \code{\link{lm}}
% first and then calls the default method.


}

\note{
  The code was adapted directly from \pkg{lmtest} (written by
  T. Hothorn, A. Zeileis, G. Millo, D. Mitchell)
  and made to work for VGLMs and S4.
  This help file also was adapted from \pkg{lmtest}.


  \emph{Approximate} LRTs might be applied to VGAMs, as
  produced by \code{\link{vgam}}, but it is probably better in inference
  to use \code{\link{vglm}} with regression splines
  (\code{\link[splines]{bs}} and
   \code{\link[splines]{ns}}).
  This methods function should not be applied to other models such as
  those produced
  by \code{\link{rrvglm}},
  by \code{\link{cqo}},
  by \code{\link{cao}}.


}

\section{Warning }{
  Several \pkg{VGAM} family functions implement distributions
  which do not satisfying the usual regularity conditions needed for
  the LRT to work. No checking or warning is given for these.

}


\value{
 An object of class \code{"VGAManova"} which contains a slot
 with the
 log-likelihood, degrees of freedom, the difference in
 degrees of freedom, likelihood ratio Chi-squared statistic
 and corresponding p value.
 These are printed by \code{stats:::print.anova()};
 see \code{\link[stats]{anova}}.


}

\seealso{
  \pkg{lmtest},
  \code{\link{vglm}}.


% \code{\link{waldtest}}
% \code{update_default},
% \code{update_formula}.



}

\examples{
set.seed(1)
pneumo <- transform(pneumo, let = log(exposure.time), x3 = runif(nrow(pneumo)))
fit1 <- vglm(cbind(normal, mild, severe) ~ let     , propodds, data = pneumo)
fit2 <- vglm(cbind(normal, mild, severe) ~ let + x3, propodds, data = pneumo)
fit3 <- vglm(cbind(normal, mild, severe) ~ let     , cumulative, data = pneumo)
# Various equivalent specifications of the LR test for testing x3
(ans1 <- lrtest(fit2, fit1))
ans2 <- lrtest(fit2, 2)
ans3 <- lrtest(fit2, "x3")
ans4 <- lrtest(fit2, . ~ . - x3)
c(all.equal(ans1, ans2), all.equal(ans1, ans3), all.equal(ans1, ans4))

# Doing it manually
(testStatistic <- 2 * (logLik(fit2) - logLik(fit1)))
(mypval <- pchisq(testStatistic, df = length(coef(fit2)) - length(coef(fit1)),
                  lower.tail = FALSE))

(ans4 <- lrtest(fit3, fit1))  # Test proportional odds (parallelism) assumption
}

\keyword{htest}


%(testStatistic <- 2 * (logLik(fit3) - logLik(fit1)))
%(mypval <- pchisq(testStatistic, df = length(coef(fit3)) - length(coef(fit1)),
%                  lower.tail = FALSE))


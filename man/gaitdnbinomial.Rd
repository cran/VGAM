\name{gaitdnbinomial}
\alias{gaitdnbinomial}
% 20211124; gaitdpoisson.Rd was the template
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated, --Truncated and Deflated
        Negative Binomial Regression


}
\description{
  Fits a generally--altered, --inflated --truncated and deflated
  negative binomial regression by MLE.
  The GAITD combo model having 7 types of special values is implemented.
  This allows mixtures of negative binomial distributions
  on nested and/or partitioned
  support as well as a multinomial logit model for
  (nonparametric) altered, inflated and deflated values.



% Truncation may include the upper tail.



}
%    eq.ap = FALSE, eq.ip = FALSE, eq.dp = FALSE,
\usage{
 gaitdnbinomial(a.mix = NULL, i.mix = NULL, d.mix = NULL,
           a.mlm = NULL, i.mlm = NULL, d.mlm = NULL,
           truncate = NULL, zero = c("size", "pobs", "pstr", "pdip"),
           eq.ap = TRUE, eq.ip = TRUE, eq.dp = TRUE,
           parallel.a = FALSE, parallel.i = FALSE, parallel.d = FALSE,
           lmunb.p = "loglink",
           lmunb.a = lmunb.p, lmunb.i = lmunb.p, lmunb.d = lmunb.p,
           lsize.p = "loglink",
           lsize.a = lsize.p, lsize.i = lsize.p, lsize.d = lsize.p,
           type.fitted = c("mean", "munbs", "sizes", "pobs.mlm",
           "pstr.mlm", "pdip.mlm", "pobs.mix", "pstr.mix", "pdip.mix",
           "Pobs.mix", "Pstr.mix", "Pdip.mix", "nonspecial", "Numer",
           "Denom.p", "sum.mlm.i", "sum.mix.i",
           "sum.mlm.d", "sum.mix.d", "ptrunc.p", "cdf.max.s"),
           gpstr.mix = ppoints(7) / 3,
           gpstr.mlm = ppoints(7) / (3 + length(i.mlm)),
           imethod = 1, mux.init = c(0.75, 0.5, 0.75, 0.5),
           imunb.p = NULL, imunb.a = imunb.p,
           imunb.i = imunb.p, imunb.d = imunb.p,
           isize.p = NULL,  isize.a = isize.p,
           isize.i = isize.p, isize.d = isize.p,
           ipobs.mix = NULL, ipstr.mix = NULL,
           ipdip.mix = NULL, ipobs.mlm = NULL,
           ipstr.mlm = NULL, ipdip.mlm = NULL,
           byrow.aid = FALSE, ishrinkage = 0.95, probs.y = 0.35,
           nsimEIM = 500, cutoff.prob = 0.999, eps.trig = 1e-7,
           nbd.max.support = 4000, max.chunk.MB = 30)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
%    gpdip.mix = ppoints(7) / 3,
%    gpdip.mlm = ppoints(7) / (3 + length(d.mlm)),
\arguments{
% \item{truncate, max.support}
  \item{truncate}{
    See \code{\link{gaitdpoisson}}.
    


  }
\item{a.mix, i.mix, d.mix}{
    See \code{\link{gaitdpoisson}}.


  }

  \item{a.mlm, i.mlm, d.mlm}{
    See \code{\link{gaitdpoisson}}.

    

}

  \item{lmunb.p, lmunb.a, lmunb.i, lmunb.d}{
    Link functions pertaining to the mean parameters.
    See \code{\link{gaitdpoisson}} where \code{llambda.p} etc. are
    the equivalent.


    
  }
  \item{lsize.p, lsize.a, lsize.i, lsize.d}{
    Link functions pertaining to the \code{size} parameters.
    See \code{\link[stats]{NegBinomial}}.


    
  }




  \item{eq.ap, eq.ip, eq.dp}{
    See \code{\link{gaitdpoisson}}.
    These apply to both \code{munb} and \code{size} parameters
    simultaneously.
    See \code{\link[stats]{NegBinomial}} also.



  }
  \item{parallel.a, parallel.i, parallel.d}{
    See \code{\link{gaitdpoisson}}.



  }

  \item{type.fitted}{
    See \code{\link{gaitdpoisson}}.


  }


  \item{gpstr.mix, gpstr.mlm}{
    See \code{\link{gaitdpoisson}}.
  
    
}
%\item{gpdip.mix, gpdip.mlm}{
%  Similar to \code{gpstr.mix} and \code{gpstr.mlm}.
%}



  \item{imethod, ipobs.mix, ipstr.mix, ipdip.mix}{
    See \code{\link{gaitdpoisson}} and
    \code{\link{CommonVGAMffArguments}}.



  }

  \item{ipobs.mlm, ipstr.mlm, ipdip.mlm}{
    See \code{\link{gaitdpoisson}}.




  }
  \item{mux.init}{
    Numeric, of length 4.
    General downward multiplier for initial values for
    the sample proportions (MLEs actually).
    See \code{\link{gaitdpoisson}}.
    The fourth value corresponds to \code{size}.
    

    
}
  
  \item{imunb.p, imunb.a, imunb.i, imunb.d}{
    See \code{\link{gaitdpoisson}};
    \code{imunb.p} is similar to  \code{ilambda.p}, etc.



  }
  \item{isize.p, isize.a, isize.i, isize.d}{
    See \code{\link{gaitdpoisson}};
    \code{isize.p} is similar to  \code{ilambda.p}, etc.



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{byrow.aid}{
  Details are at \code{\link{Gaitdpois}}.


  }
  \item{zero}{
    See \code{\link{gaitdpoisson}} and
    \code{\link{CommonVGAMffArguments}}.



  }
  \item{nsimEIM, cutoff.prob, eps.trig}{
    See \code{\link{negbinomial}}.



  }
  \item{nbd.max.support, max.chunk.MB}{
    See \code{\link{negbinomial}}.



  }
}
\details{
  The GAITD--NB combo model is the pinnacle of GAITD regression
  for counts because it potentially handles
  underdispersion, 
  equidispersion and
  overdispersion relative to the Poisson,
  as well as
  alteration,
  inflation,
  deflation and
  truncation at arbitrary support points.
  In contrast, \code{\link{gaitdpoisson}} cannot handle
  overdispersion so well.
  The GAITD--NB is so flexible that it can accommodate up to
  seven modes.



  The full
  GAITD--NB--NB--MLM--NB-MLM--NB-MLM combo model
  may be fitted with this family function.
  There are seven types of special values and all arguments for these
  may be used in a single model.
  Here, the MLM represents the nonparametric while the NB
  refers to the negative binomial mixtures.
  The defaults for this function correspond to an
  ordinary negative binomial
  regression so that \code{\link{negbinomial}} is called instead.




  While much of the documentation here draws upon
  \code{\link{gaitdpoisson}}, there are additional
  details here because the NBD is a \emph{two} parameter
  distribution that handles \emph{overdispersion} relative
  to the Possion.
  Consequently, this family function is exceeding flexible
  and there are many more pitfalls to avoid.
 




  
  The order of the linear/additive predictors is best explained by
  an example.
  Suppose a combo model has
  \code{length(a.mix) > 3} and
  \code{length(i.mix) > 3},
  \code{length(d.mix) > 3},
  \code{a.mlm = 3:5},
  \code{i.mlm = 6:9} and
  \code{d.mlm = 10:12}, say.
  Then \code{loglink(munb.p)} and \code{loglink(size.p)}
  are the first two.
  The third is \code{multilogitlink(pobs.mix)} followed
  by \code{loglink(munb.a)}
  and \code{loglink(size.a)}
  because \code{a.mix} is long enough.
  The sixth is \code{multilogitlink(pstr.mix)} followed
  by \code{loglink(munb.i)}
  and \code{loglink(size.i)}
  because \code{i.mix} is long enough.
  The ninth is \code{multilogitlink(pdip.mix)} followed
  by \code{loglink(munb.d)}
  and \code{loglink(size.d)}
  because \code{d.mix} is long enough.
  Next are the probabilities for the \code{a.mlm} values.
  Then are the probabilities for the \code{i.mlm} values.
  Lastly are the probabilities for the \code{d.mlm} values.
  All the probabilities are estimated by one big MLM and effectively
  the \code{"(Others)"} column of left over probabilities is
  associated with the nonspecial values.
  These might be called the
  \emph{nonspecial baseline probabilities} (NBP)
  or reserve probabilities.
  The dimension of the vector of linear/additive predictors here
  is \eqn{M=21}.


  
% 11 + length(c(3:12))







Apart from the order of the linear/additive predictors,
the following are (or should be) equivalent:
\code{gaitdnbinomial()} and \code{negbinomial()},
\code{gaitdnbinomial(a.mix = 0)} and \code{zanegbinomial(zero = "pobs0")},
\code{gaitdnbinomial(i.mix = 0)} and \code{zinegbinomial(zero = "pstr0")},
\code{gaitdnbinomial(truncate = 0)} and \code{posnegbinomial()}.
Likewise, if
\code{a.mix} and \code{i.mix} are assigned a scalar then
it effectively moves that scalar to \code{a.mlm} and \code{i.mlm}
because there is no
parameters such as \code{munb.i} being estimated.
Thus
\code{gaitdnbinomial(a.mix = 0)} and \code{gaitdnbinomial(a.mlm = 0)}
are the effectively same, and ditto for
\code{gaitdnbinomial(i.mix = 0)} and \code{gaitdnbinomial(i.mlm = 0)}.







}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.




  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted},
  returns the mean \eqn{\mu}{mu} by default.
  See the information above on \code{type.fitted}.



  
}
\references{
Yee, T. W. and Ma, C. (2022).
Generally--altered, --inflated, --truncated and --deflated
regression, with application to heaped and seeped data.
\emph{In preparation}.



% count data.
%, \bold{3}, 15--41.
% with application to heaped and seeped counts.


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
  See \code{\link{gaitdpoisson}}.
  Also, having \code{eq.ap = TRUE}, \code{eq.ip = TRUE}
  and \code{eq.dp = TRUE} is often needed to obtain
  initial values that are good enough because they borrow
  strength across the different operators.
  It is usually easy to relax these assumptions later.



  This family function is under constant development and
  future changes will occur.




}





\author{ T. W. Yee}
\note{
  If \code{length(a.mix)} is 1 then effectively this becomes a
  value of \code{a.mlm}.
  If \code{length(a.mix)} is 2 then an error message will be issued
  (overfitting really).
  If \code{length(a.mix)} is 3 then this is almost overfitting too.
  Hence \code{length(a.mix)} should be 4 or more.
  Ditto for \code{length(i.mix)} and \code{length(d.mix)}.


  
  See \code{\link{gaitdpoisson}} for notes about numerical
  problems that can easily arise. With the NBD there is
  even more potential trouble that can occur.
  In particular, good initial values are more necessary so
  it pays to experiment with arguments such as
  \code{imunb.p} and \code{isize.p}, as well as
  fitting an intercept-only model first before adding
  covariates and using \code{etastart}.


 
  Currently \code{max.support} is missing because only
  \code{Inf} is handled. This might change later.


  

  
}

\seealso{
  \code{\link{Gaitdnbinom}},
  \code{\link{multinomial}},
  \code{\link{rootogram4}},
  \code{\link{specials}},
  \code{\link{plotdgaitd}},
  \code{\link{spikeplot}},
  \code{\link{meangaitd}},
  \code{\link{KLD}},
  \code{\link{gaitdpoisson}},
  \code{\link{gaitdlog}},
  \code{\link{gaitdzeta}},
  \code{\link{multilogitlink}},
  \code{\link{multinomial}},
  \code{\link{goffset}},
  \code{\link{Trunc}},
  \code{\link{negbinomial}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.



% \code{\link{zapoisson}},
% \code{\link{zipoisson}},
% \code{\link{pospoisson}},
% \code{\link{CommonVGAMffArguments}},
% \code{\link{simulate.vlm}}.


% \code{\link{Trunc}},
% \code{\link{gaitlog.mix}},



% \code{\link{gatnbinomial.mix}},
% \code{\link{gatnbinomial.mlm}},
% \code{\link{gatpoisson.mix}},
% \code{\link{multinomial}},
% \code{\link{zapoisson}},
% \code{\link{gipoisson}},



}

\examples{
i.mix <- c(5, 10, 12, 16)  # Inflate these values parametrically
i.mlm <- c(14, 15)  # Inflate these values
a.mix <- c(1, 6, 13, 20)  # Alter these values
tvec <- c(3, 11)   # Truncate these values
pstr.mlm <- 0.1  # So parallel.i = TRUE
pobs.mix <- pstr.mix <- 0.1; set.seed(1)
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, munb.p = exp(2 + 0.0 * x2), size.p = exp(1))
gdata <- transform(gdata,
  y1 = rgaitdnbinom(nn, size.p, munb.p, a.mix = a.mix, i.mix = i.mix,
                    pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                    i.mlm = i.mlm, pstr.mlm = pstr.mlm,
                    truncate = tvec))
gaitdnbinomial(a.mix = a.mix, i.mix = i.mix, i.mlm = i.mlm)
with(gdata, table(y1))
fit1 <- vglm(y1 ~ 1, crit = "coef", trace = TRUE, data = gdata,
             gaitdnbinomial(a.mix = a.mix, i.mix = i.mix, i.mlm = i.mlm,
                            parallel.i = TRUE, eq.ap = TRUE,
                            eq.ip = TRUE, truncate = tvec))
head(fitted(fit1, type.fitted = "Pstr.mix"))
head(predict(fit1))
t(coef(fit1, matrix = TRUE))  # Easier to see with t()
summary(fit1)
\dontrun{ spikeplot(with(gdata, y1), lwd = 2)
plotdgaitd(fit1, new.plot = FALSE, offset.x = 0.2, all.lwd = 2)  }
}
\keyword{models}
\keyword{regression}


%gapoisson(lpobs0 = "logitlink", lmunb = "loglink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%gapoissonff(lmunb = "loglink", lonempobs0 = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")




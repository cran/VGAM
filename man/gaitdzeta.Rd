\name{gaitdzeta}
\alias{gaitdzeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally Altered, Inflated, Truncated
  and Deflated
  Zeta Regression

}
\description{
  Fits a generally altered, inflated, truncated
  and deflated zeta regression by MLE.
  The GAITD combo model having 7 types of special
  values is implemented.  This allows mixtures of
  zetas on nested and/or partitioned support as
  well as a multinomial logit model for altered,
  inflated and deflated values.


% The truncation may include values in the upper tail.


}
\usage{
gaitdzeta(a.mix = NULL, i.mix = NULL, d.mix = NULL,
         a.mlm = NULL, i.mlm = NULL, d.mlm = NULL,
         truncate = NULL, max.support = Inf,
         zero = c("pobs", "pstr", "pdip"), eq.ap = TRUE, eq.ip = TRUE,
         eq.dp = TRUE, parallel.a = FALSE,
         parallel.i = FALSE, parallel.d = FALSE,
         lshape.p = "loglink", lshape.a = lshape.p,
         lshape.i = lshape.p, lshape.d = lshape.p,
         type.fitted = c("mean", "shapes", "pobs.mlm", "pstr.mlm",
         "pdip.mlm", "pobs.mix", "pstr.mix", "pdip.mix", "Pobs.mix",
         "Pstr.mix", "Pdip.mix", "nonspecial",
         "Numer", "Denom.p", "sum.mlm.i", "sum.mix.i", "sum.mlm.d",
         "sum.mix.d", "ptrunc.p", "cdf.max.s"),
         gshape.p = -expm1(-ppoints(7)), gpstr.mix = ppoints(7) / 3,
         gpstr.mlm = ppoints(7) / (3 + length(i.mlm)),
         imethod = 1, mux.init = c(0.75, 0.5, 0.75),
         ishape.p = NULL, ishape.a = ishape.p,
         ishape.i = ishape.p, ishape.d = ishape.p,
         ipobs.mix = NULL, ipstr.mix = NULL, ipdip.mix = NULL,
         ipobs.mlm = NULL, ipstr.mlm = NULL, ipdip.mlm = NULL,
         byrow.aid = FALSE, ishrinkage = 0.95, probs.y = 0.35)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{truncate, max.support}{
    See \code{\link{gaitdpoisson}}.
    Only \code{max.support = Inf} is allowed
    because some equations are intractable.

    
  }
  \item{a.mix, i.mix, d.mix}{
    See \code{\link{gaitdpoisson}}.

  }
  \item{a.mlm, i.mlm, d.mlm}{
    See \code{\link{gaitdpoisson}}.

  }

  \item{lshape.p, lshape.a, lshape.i, lshape.d}{
    Link functions.
    See \code{\link{gaitdpoisson}}
    and \code{\link{Links}} for more choices
    and information.  Actually, it is usually
    a good idea to set these arguments equal to
    \code{\link[VGAMextra]{zetaffMlink}} because
    the log-mean is the first linear/additive
    predictor so it is like a Poisson regression.

    
  }

  \item{eq.ap, eq.ip, eq.dp}{
    Single logical each.
    See \code{\link{gaitdpoisson}}



  }


  \item{parallel.a, parallel.i, parallel.d}{
    Single logical each.
    See \code{\link{gaitdpoisson}}.



  }

  \item{type.fitted, mux.init}{
    See \code{\link{gaitdpoisson}}.


  }


  \item{imethod, ipobs.mix, ipstr.mix, ipdip.mix}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{gaitdpoisson}} for information.



% ipobs0,



  }
  \item{ipobs.mlm, ipstr.mlm, ipdip.mlm, byrow.aid}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{gaitdpoisson}} for information.
  



  }
  \item{gpstr.mix, gpstr.mlm}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{gaitdpoisson}} for information.


  }
% \item{gpdip.mix, gpdip.mlm}{
% See \code{\link{CommonVGAMffArguments}}
% and \code{\link{gaitdpoisson}} for information.


% }
  \item{gshape.p, ishape.p}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{gaitdpoisson}} for information.
  The former is used only if the latter is
  not given.  Practical experience has shown
  that good initial values are needed, so
  if convergence is not obtained then try a
  finer grid.


% ipobs0,



  }
  \item{ishape.a, ishape.i, ishape.d}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{gaitdpoisson}} for information.


% ipobs0,



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}}
  for information.



  }
  \item{zero}{
    See \code{\link{gaitdpoisson}}
    and \code{\link{CommonVGAMffArguments}}
    for information.



  }
}
\details{
  Many details to this family function can be
  found in \code{\link{gaitdpoisson}} because it
  is also a 1-parameter discrete distribution.
  This function currently does not handle
  multiple responses.  Further details are at
  \code{\link{Gaitdzeta}}.




  As alluded to above, when there are covariates
  it is much more interpretable to model
  the mean rather than the shape parameter.
  Hence \code{\link[VGAMextra]{zetaffMlink}}
  is recommended.  (This might become the default
  in the future.)  So installing \pkg{VGAMextra}
  is a good idea.





Apart from the order of the linear/additive predictors,
the following are (or should be) equivalent:
\code{gaitdzeta()} and \code{zetaff()},
\code{gaitdzeta(a.mix = 1)} and \code{oazeta(zero = "pobs1")},
\code{gaitdzeta(i.mix = 1)} and \code{oizeta(zero = "pstr1")},
\code{gaitdzeta(truncate = 1)} and \code{otzeta()}.
The functions
\code{\link[VGAMdata]{oazeta}},
\code{\link[VGAMdata]{oizeta}} and
\code{\link[VGAMdata]{otzeta}}
have been placed in \pkg{VGAMdata}.





}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such
  as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.




% The \code{fitted.values} slot of the fitted object,
% which should be extracted by the generic function \code{fitted},
% are similar to \code{\link{gaitdzeta.mlm}}.




}
% \references{
%Yee, T. W. and Ma, C. (2023)
%Generally altered, inflated, truncated and deflated regression.
%\emph{In preparation}.

%% , \bold{3}, 15--41.


% }




\section{Warning }{
    See \code{\link{gaitdpoisson}}.



}



\author{ T. W. Yee}
\note{
  See \code{\link{gaitdpoisson}}.

  

% See \code{\link{gaitdzeta.mlm}} for other general details.


}

\seealso{
  \code{\link{Gaitdzeta}},
  \code{\link{zetaff}},
  \code{\link[VGAMextra]{zetaffMlink}},
  \code{\link{Gaitdpois}},
  \code{\link{gaitdpoisson}},
  \code{\link{gaitdlog}},
  \code{\link{spikeplot}},
  \code{\link{goffset}},
  \code{\link{Trunc}},
  \code{\link[VGAMdata]{oazeta}},
  \code{\link[VGAMdata]{oizeta}},
  \code{\link[VGAMdata]{otzeta}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{rootogram4}},
  \code{\link{simulate.vlm}}.


% \code{\link{gaitdzeta.mlm}},





}

\examples{
\dontrun{
avec <- c(5, 10)  # Alter these values parametrically
ivec <- c(3, 15)  # Inflate these values
tvec <- c(6, 7)   # Truncate these values
set.seed(1); pobs.a <- pstr.i <- 0.1
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, shape.p = logitlink(2, inverse = TRUE))
gdata <- transform(gdata,
  y1 = rgaitdzeta(nn, shape.p, a.mix = avec, pobs.mix = pobs.a,
                  i.mix = ivec, pstr.mix = pstr.i, truncate = tvec))
gaitdzeta(a.mix = avec, i.mix = ivec)
with(gdata, table(y1))
spikeplot(with(gdata, y1), las = 1)
fit7 <- vglm(y1 ~ 1, trace = TRUE, data = gdata, crit = "coef",
             gaitdzeta(i.mix = ivec, truncate = tvec,
                       a.mix = avec, eq.ap = TRUE, eq.ip = TRUE))
head(fitted(fit7, type.fitted = "Pstr.mix"))
head(predict(fit7))
t(coef(fit7, matrix = TRUE))  # Easier to see with t()
summary(fit7)
spikeplot(with(gdata, y1), lwd = 2, ylim = c(0, 0.6), xlim = c(0, 20))
plotdgaitd(fit7, new.plot = FALSE, offset.x = 0.2, all.lwd = 2)
}}
\keyword{models}
\keyword{regression}



%                     truncate = tvec, max.support = max.support
%                               eq.ip = TRUE, max.support = max.support



\name{gaitlog.mix}
\alias{gaitlog.mix}
%\alias{galogff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated and --Truncated Logarithmic
  Regression Family Function
  (GAIT--Log--Log--Log Mixture Variant)

}
\description{
  Fits a generally--altered, --inflated and --truncated logarithmic
  regression
  (mixtures of logarithmic distributions on
  nested and/or partitioned supports).


% The truncation may include values in the upper tail.


}
\usage{
gaitlog.mix(alter = NULL, inflate = NULL, truncate = NULL,
    max.support = Inf, zero = c("pobs.a", "pstr.i"),
    eq.ap = FALSE, eq.ip = FALSE, lshape.p = "logitlink",
    lpobs.a = "logitlink", lshape.a = "logitlink",
    lpstr.i = "logitlink", lshape.i = "logitlink",
    type.fitted = c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
    "prob.a", "prob.i", "prob.t", "lhs.prob"),
    imethod = 1, gshape.p = -expm1(-7 * ppoints(12)), ishape.p = NULL,
    ishape.a = NULL, ishape.i = NULL,
    ipobs.a = NULL, ipstr.i = NULL, ishrinkage = 0.95, probs.y = 0.35)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{alter, inflate, truncate}{
    See \code{\link{gaitpoisson.mix}}.

  }
  \item{lshape.p, lshape.a, lshape.i}{
    Link functions.
    See \code{\link{gaitpoisson.mix}}
    and \code{\link{Links}} for more choices and information.

  }
  \item{lpobs.a, lpstr.i}{
    Link functions.
    See \code{\link{gaitpoisson.mix}}
    and \code{\link{Links}} for more choices and information.

  }

  \item{eq.ap, eq.ip}{
    Single logical each.
    See \code{\link{gaitpoisson.mix}}



  }


  \item{type.fitted, max.support}{
    See \code{\link{gaitpoisson.mix}}.


  }


  \item{imethod, ipobs.a, ipstr.i}{
  See \code{\link{CommonVGAMffArguments}} for information.


% ipobs0,



  }
  \item{gshape.p, ishape.p}{
  See \code{\link{CommonVGAMffArguments}} for information.
  The former is used only if the latter is not given.
  Practical experience has shown that good initial values are needed,
  so if convergence is not obtained then try a finer grid.


% ipobs0,



  }
  \item{ishape.a, ishape.i}{
  See \code{\link{CommonVGAMffArguments}} for information.


% ipobs0,



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{zero}{
    See \code{\link{gaitpoisson.mix}}
    and \code{\link{CommonVGAMffArguments}} for information.



  }
}
\details{
  Although the full GAIT--Log--Log--Log model may be fitted,
  the two submodels that may be fitted can be abbreviated
  GAT--Log--Log or
  GIT--Log--Log,
  which is where the inner distribution for
  ordinary values is the logarithmic distribution, and
  the outer distribution for the altered or inflated values
  is another logarithmic distribution with a different shape parameter
  by default.
  Thus for the GAT model
  the distribution being fitted is a (spliced) mixture
  of two logarithmic with differing (partitioned) support.
  Likewise, for the GIT model
  the distribution being fitted is a mixture
  of two logarithmic with nested support.
  Two shape parameters at a time may be constrained to be equal using
  \code{eq.ap} or \code{eq.ip}.



  For the GIT model,
  by default, a logistic regression models the (structural)
  probability \code{pstr.i} that the response is inflated.



  This function currently does not handle multiple responses.
  Further details are at \code{\link{Gaitlog}}.
% An alternative variant of this distribution,
% more unstructured in nature, is based
% on the multinomial logit model---see \code{\link{gaitlog.mlm}}.




  For the GIT model,
  the ordering of the linear/additive predictors corresponds to
  \code{length(inflate)} equalling 0, 1, and more than 1;
  the dimension grows accordingly.
  The same idea holds for the GAT model.



  Apart from the order of the linear/additive predictors,
  the following are (or should be) equivalent:
\code{gaitlog.mix()} and \code{logff()},
\code{gaitlog.mix(alter = 1)} and \code{oalog(zero = "pobs1")},
\code{gaitlog.mix(inflate = 1)} and \code{oilog(zero = "pstr1")},
\code{gaitlog.mix(truncate = 1)} and \code{otlog()}.




}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.



% The \code{fitted.values} slot of the fitted object,
% which should be extracted by the generic function \code{fitted},
% are similar to \code{\link{gaitlog.mlm}}.




}
% \references{
% Yee, T. W. and Ma, C. (2020).
% Generally--altered, --inflated and --truncated regression,
% with application to heaped and seeped count data.
% \emph{In preparation}.


%% , \bold{3}, 15--41.


% }




\section{Warning }{
    See \code{\link{gaitpoisson.mix}}



}



\author{ T. W. Yee}
\note{
  See \code{\link{gaitpoisson.mix}}


% See \code{\link{gaitlog.mlm}} for other general details.


}

\seealso{
  \code{\link{Gaitlog}},
  \code{\link{gaitpoisson.mix}},
  \code{\link{gaitpoisson.mlm}},
  \code{\link{oalog}},
  \code{\link{oilog}},
  \code{\link{otlog}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{rootogram4}},
  \code{\link{simulate.vlm}}.


% \code{\link{gaitlog.mlm}},





}

\examples{
avec <- c(5, 10)  # Alter these values
ivec <- c(3, 15)  # Inflate these values
tvec <- c(6, 7)   # Truncate these values
max.support <- 20; set.seed(1)
pobs.a <- logitlink(-1, inverse = TRUE)  # About 0.27
pstr.i <- logitlink(-1, inverse = TRUE)  # About 0.27
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, shape.p = logitlink(2+0.5*x2, inverse = TRUE))
gdata <- transform(gdata,
  y1 = rgaitlog(nn, shape.p, alter.mix = avec, pobs.mix.a = pobs.a,
                inflate.mix = ivec, pstr.mix.i = pstr.i, truncate = tvec,
                max.support = max.support))
gaitlog.mix(alter = avec, inflate = ivec, max.support = max.support)
with(gdata, table(y1))
gaitlxfit <- vglm(y1 ~ x2, crit = "coef", trace = TRUE, data = gdata,
                  gaitlog.mix(inflate = ivec, truncate = tvec,
                              max.support = max.support,
                              alter = avec, eq.ap = TRUE, eq.ip = TRUE))
head(fitted(gaitlxfit, type.fitted = "Pstr.i"))
head(predict(gaitlxfit))
coef(gaitlxfit, matrix = TRUE)
summary(gaitlxfit)
}
\keyword{models}
\keyword{regression}



%                     truncate = tvec, max.support = max.support
%                               eq.ip = TRUE, max.support = max.support



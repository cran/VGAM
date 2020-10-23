\name{gaitzeta}
\alias{gaitzeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated and --Truncated Zeta
  Regression

}
\description{
  Fits a generally--altered, --inflated and --truncated zeta
  regression by MLE.
  The GAIT combo model having 5 types of special values is implemented.
  This allows mixtures of zetas on nested and/or partitioned
  support as well as a multinomial logit model for altered and
  inflated values.


% The truncation may include values in the upper tail.


}
\usage{
gaitzeta(alt.mix = NULL, inf.mix = NULL, alt.mlm = NULL,
     inf.mlm = NULL, truncate = NULL, max.support = Inf,
     zero = c("pobs", "pstr"), eq.ap = FALSE, eq.ip = FALSE,
     parallel.ap = FALSE, parallel.ip = FALSE, lshape.p = "loglink",
     lshape.a = "loglink", lshape.i = "loglink",
     type.fitted = c("mean", "shapes", "pobs.mlm", "pstr.mlm",
     "pobs.mix", "pstr.mix", "Pobs.mix", "Pstr.mix", "nonspecial",
     "Numer", "Denom.p", "sum.mlm.i", "sum.mix.i", "ptrunc.p",
     "cdf.max.s"), gshape.p = 1 + exp(-seq(7)),
     gpstr.mix = ppoints(9) / 2,
     gpstr.mlm = ppoints(9) / (2 + length(inf.mlm)), imethod = 1,
     imux = 0.5, ishape.p = NULL, ishape.a = ishape.p,
     ishape.i = ishape.p, ipobs.mix = NULL, ipstr.mix = NULL,
     ipobs.mlm = NULL, ipstr.mlm = NULL, byrow.ai = FALSE,
     ishrinkage = 0.95, probs.y = 0.35)



}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{truncate, max.support}{
    See \code{\link{gaitpoisson}}.
    Only \code{max.support = Inf} is allowed
    because some equations are intractable.

    
  }
  \item{alt.mix, inf.mix}{
    See \code{\link{gaitpoisson}}.

  }
  \item{alt.mlm, inf.mlm}{
    See \code{\link{gaitpoisson}}.

  }

  \item{lshape.p, lshape.a, lshape.i}{
    Link functions.
    See \code{\link{gaitpoisson}}
    and \code{\link{Links}} for more choices and information.

    
  }

  \item{eq.ap, eq.ip}{
    Single logical each.
    See \code{\link{gaitpoisson}}



  }


  \item{parallel.ap, parallel.ip}{
    Single logical each.
    See \code{\link{gaitpoisson}}.



  }

  \item{type.fitted, imux}{
    See \code{\link{gaitpoisson}}.


  }


  \item{imethod, ipobs.mix, ipstr.mix}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{gaitpoisson}} for information.



% ipobs0,



  }
  \item{ipobs.mlm, ipstr.mlm, byrow.ai}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{gaitpoisson}} for information.
  



  }
  \item{gpstr.mix, gpstr.mlm}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{gaitpoisson}} for information.


  }
  \item{gshape.p, ishape.p}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{gaitpoisson}} for information.
  The former is used only if the latter is not given.
  Practical experience has shown that good initial values are needed,
  so if convergence is not obtained then try a finer grid.


% ipobs0,



  }
  \item{ishape.a, ishape.i}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{gaitpoisson}} for information.


% ipobs0,



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{zero}{
    See \code{\link{gaitpoisson}}
    and \code{\link{CommonVGAMffArguments}} for information.



  }
}
\details{
  Many details to this family function can be found
  in \code{\link{gaitpoisson}} because it is also a 1-parameter
  discrete distribution.
  This function currently does not handle multiple responses.
  Further details are at \code{\link{Gaitzeta}}.






  Apart from the order of the linear/additive predictors,
  the following are (or should be) equivalent:
\code{gaitzeta()} and \code{zetaff()},
\code{gaitzeta(alt.mix = 1)} and \code{oazeta(zero = "pobs1")},
\code{gaitzeta(inf.mix = 1)} and \code{oizeta(zero = "pstr1")},
\code{gaitzeta(truncate = 1)} and \code{otzeta()}.
The functions
\code{\link[VGAMdata]{oazeta}},
\code{\link[VGAMdata]{oizeta}} and
\code{\link[VGAMdata]{otzeta}}
have been placed in \pkg{VGAMdata}.





}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.



% The \code{fitted.values} slot of the fitted object,
% which should be extracted by the generic function \code{fitted},
% are similar to \code{\link{gaitzeta.mlm}}.




}
% \references{
% Yee, T. W. and Ma, C. (2020).
% Generally--altered, --inflated and --truncated regression,
% with application to heaped and seeped counts.  % count data.
% \emph{In preparation}.


%% , \bold{3}, 15--41.


% }




\section{Warning }{
    See \code{\link{gaitpoisson}}.



}



\author{ T. W. Yee}
\note{
  See \code{\link{gaitpoisson}}.

  

% See \code{\link{gaitzeta.mlm}} for other general details.


}

\seealso{
  \code{\link{Gaitzeta}},
  \code{\link{zetaff}},
  \code{\link{Gaitpois}},
  \code{\link{gaitpoisson}},
  \code{\link{gaitlog}},
  \code{\link[VGAMdata]{oazeta}},
  \code{\link[VGAMdata]{oizeta}},
  \code{\link[VGAMdata]{otzeta}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{rootogram4}},
  \code{\link{simulate.vlm}}.


% \code{\link{gaitzeta.mlm}},





}

\examples{
\dontrun{
avec <- c(5, 10)  # Alter these values parametrically
ivec <- c(3, 15)  # Inflate these values
tvec <- c(6, 7)   # Truncate these values
set.seed(1); pobs.mix <- pstr.mix <- 0.1
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, shape.p = loglink(0.1+0.1*x2, inverse = TRUE))
gdata <- transform(gdata,
  y1 = rgaitzeta(nn, shape.p, alt.mix = avec, pobs.mix = pobs.mix,
                 inf.mix = ivec, pstr.mix = pstr.mix, truncate = tvec))
gaitzeta(alt.mix = avec, inf.mix = ivec, max.support = max.support)
with(gdata, table(y1))
gaitzxfit <- vglm(y1 ~ x2, crit = "coef", trace = TRUE, data = gdata,
                  gaitzeta(inf.mix = ivec, truncate = tvec,
                           alt.mix = avec, eq.ap = TRUE, eq.ip = TRUE))
head(fitted(gaitzxfit, type.fitted = "Pstr.mix"))
head(predict(gaitzxfit))
t(coef(gaitzxfit, matrix = TRUE))  # Easier to see with t()
summary(gaitzxfit, HDEtest = FALSE)  # summary(gaitzxfit) is better
}}
\keyword{models}
\keyword{regression}



%                     truncate = tvec, max.support = max.support
%                               eq.ip = TRUE, max.support = max.support



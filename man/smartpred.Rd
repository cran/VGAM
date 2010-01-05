\name{smartpred}
\alias{smartpred}
\title{ Smart Prediction }
\description{
  Data-dependent parameters in formula terms
  can cause problems in when predicting.
  The \pkg{smartpred} package for \R and S-PLUS saves
  data-dependent parameters on the object so that the bug is fixed.
  The \code{\link[stats]{lm}} and \code{\link[stats]{glm}} functions have
  been fixed properly. Note that the \pkg{VGAM} package by T. W. Yee
  automatically comes with smart prediction.

}
%\usage{
%lm()
%glm()
%ns()
%bs()
%poly()
%scale()
%vglm()
%rrvglm()
%vgam()
%cao()
%cqo()
%uqo()
%}
\value{
  Returns the usual object, but with one list/slot component called
  \code{smart.prediction} containing any data-dependent parameters.
}
\section{Side Effects}{
  The variables
  \code{.max.smart},
  \code{.smart.prediction} and 
  \code{.smart.prediction.counter}
  are created while the model is being fitted.
  In \R they are created in a new environment called \code{smartpredenv}.
  In S-PLUS they are created in frame 1.
  These variables are deleted after the model has been fitted.
  However, in \R,
  if there is an error in the model fitting function or the fitting
  model is killed (e.g., by typing control-C) then these variables will
  be left in \code{smartpredenv}.  At the beginning of model fitting,
  these variables are deleted if present in \code{smartpredenv}.

  During prediction, the variables
  \code{.smart.prediction} and 
  \code{.smart.prediction.counter}
  are reconstructed and read by the smart functions when the model
  frame is re-evaluated. 
  After prediction, these variables are deleted. 

  If the modelling function is used with argument \code{smart=FALSE}
  (e.g., \code{vglm(..., smart=FALSE)}) then smart prediction will not
  be used, and the results should match with the original \R or S-PLUS
  functions.

}
\details{
  \R version 1.6.0 introduced a partial fix for the prediction
  problem because it does not work all the time,
  e.g., for terms such as
  \code{I(poly(x, 3))}, 
  \code{poly(c(scale(x)), 3)},
  \code{bs(scale(x), 3)}, 
  \code{scale(scale(x))}.
  See the examples below.
  Smart prediction, however, will always work.

% albeit, not so elegantly.

  The basic idea is that the functions in the formula are now smart, and the
  modelling functions make use of these smart functions.  Smart prediction
  works in two ways: using \code{\link{smart.expression}}, or using a
  combination of \code{\link{put.smart}} and \code{\link{get.smart}}.

}

\author{T. W. Yee and T. J. Hastie}
\note{
  In S-PLUS you will need to load in the \pkg{smartpred} library with
  the argument \code{first=T}, e.g.,
  \code{library(smartpred, lib="./mys8libs", first=T)}.
  Here, \code{mys8libs} is the name of a directory of installed packages.
  To install the smartpred package in Linux/Unix, type something like
  \code{Splus8 INSTALL -l ./mys8libs ./smartpred_0.8-2.tar.gz}.

}

%\note{
%  In \R and
%  prior to the \pkg{VGAM} package using name spaces, the location of the
%  variables was the workspace.  The present use of \code{smartpredenv}
%  is superior, and is somewhat similar to the S-PLUS implementation in
%  that the user is more oblivious to its existence.
%
%}

\seealso{
  \code{\link{get.smart.prediction}},
  \code{\link{get.smart}},
  \code{\link{put.smart}},
  \code{\link{smart.expression}},
  \code{\link{smart.mode.is}},
  \code{\link{setup.smart}},
  \code{\link{wrapup.smart}}.
  Commonly used data-dependent functions include
  \code{\link[base]{scale}}, 
  \code{\link[stats]{poly}}, 
  \code{\link[splines]{bs}}, 
  \code{\link[splines]{ns}}.
  In \R, 
  the functions \code{\link[splines]{bs}}
  and \code{\link[splines]{ns}} are in the
  \pkg{splines} package, and this library is automatically
  loaded in because it contains compiled code that 
  \code{\link[splines]{bs}} and \code{\link[splines]{ns}} call.

  The website \url{http://www.stat.auckland.ac.nz/~yee}
  contains more information such as how to write a
  smart function, and other technical details.

  The functions \code{\link[VGAM]{vglm}},
  \code{\link[VGAM]{vgam}},
  \code{\link[VGAM]{rrvglm}}
  and
  \code{\link[VGAM]{cqo}}
  in T. W. Yee's \pkg{VGAM}
  package are examples of modelling functions that employ smart prediction.

}
\section{WARNING }{
  In S-PLUS,
  if the \code{"bigdata"} library is loaded then it is
  \code{detach()}'ed. This is done because
  \code{scale} cannot be made smart if \code{"bigdata"} is loaded
  (it is loaded by default in the Windows version of
  Splus 8.0, but not in Linux/Unix).
  The function \code{\link[base]{search}} tells what is
  currently attached.

  In \R and S-PLUS the functions
  \code{\link[splines]{predict.bs}} and
  \code{predict.ns}
  are not smart.
  That is because they operate on objects that contain attributes only
  and do not have list components or slots.
  In \R the function
  \code{\link[stats:poly]{predict.poly}} is not smart.

}

\examples{
# Create some data first
n = 20
set.seed(86) # For reproducibility of the random numbers
x = sort(runif(n)) 
y = sort(runif(n))
\dontrun{if(is.R()) library(splines)   # To get ns() in R
}

# This will work for R 1.6.0 and later, but fail for S-PLUS
fit = lm(y ~ ns(x, df=5))
\dontrun{plot(x, y)
lines(x, fitted(fit))
newx = seq(0, 1, len=n)
points(newx, predict(fit, data.frame(x=newx)), type="b", col=2, err=-1)
}

# The following fails for R 1.6.x and later but works with smart prediction
\dontrun{fit = lm(y ~ ns(scale(x), df=5))
fit$smart.prediction
plot(x, y)
lines(x, fitted(fit))
newx = seq(0, 1, len=n)
points(newx, predict(fit, data.frame(x=newx)), type="b", col=2, err=-1)
}

# The following requires the VGAM package to be loaded 
\dontrun{library(VGAM)
fit = vlm(y ~ ns(scale(x), df=5))
fit@smart.prediction
plot(x, y)
lines(x, fitted(fit))
newx = seq(0, 1, len=n)
points(newx, predict(fit, data.frame(x=newx)), type="b", col=2, err=-1)
}
}
%\keyword{smart}
\keyword{models}
\keyword{regression}
\keyword{programming}



%lm(..., smart=TRUE)
%glm(..., smart=TRUE)
%ns()
%bs()
%poly()
%scale()
%vglm(..., smart=TRUE)
%rrvglm(..., smart=TRUE)
%vgam(..., smart=TRUE)
%cao(..., smart=TRUE)
%cqo(..., smart=TRUE)
%uqo(..., smart=TRUE)

%library(smartpred, lib="./mys8libs", first=T)





\name{smartpred}
\alias{smartpred}
\alias{sm.bs}
\alias{sm.ns}
\alias{sm.scale}
\alias{sm.scale.default}
\alias{sm.poly}
\title{ Smart Prediction }
\description{
  Data-dependent parameters in formula terms
  can cause problems in when predicting.
  The \pkg{smartpred} package
  saves
  data-dependent parameters on the object so that the bug is fixed.
  The \code{\link[stats]{lm}} and \code{\link[stats]{glm}} functions have
  been fixed properly. Note that the \pkg{VGAM} package by T. W. Yee
  automatically comes with smart prediction.




}
\usage{
sm.bs(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE, 
      Boundary.knots = range(x))
sm.ns(x, df = NULL, knots = NULL, intercept = FALSE,
      Boundary.knots = range(x))
sm.poly(x, ..., degree = 1, coefs = NULL, raw = FALSE) 
sm.scale(x, center = TRUE, scale = TRUE) 
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



\arguments{
  \item{x}{
  The \code{x} argument is actually common to them all.


  }

  \item{df, knots, intercept, Boundary.knots}{
  See \code{\link[splines]{bs}} and/or
      \code{\link[splines]{ns}}.


  }
  \item{degree, \dots, coefs, raw}{
  See \code{\link[stats]{poly}}.


  }
  \item{center, scale}{
  See \code{\link[base]{scale}}.


  }

}




\value{
  The usual value returned by
  \code{\link[splines]{bs}}, 
  \code{\link[splines]{ns}},
  \code{\link[stats]{poly}} and
  \code{\link[base]{scale}},
  When used with functions such as \code{\link[VGAM]{vglm}}
  the data-dependent parameters are  saved on one slot component called
  \code{smart.prediction}.


}
\section{Side Effects}{
  The variables
  \code{.max.smart},
  \code{.smart.prediction} and 
  \code{.smart.prediction.counter}
  are created while the model is being fitted.
  They are created in a new environment called \code{smartpredenv}.
  These variables are deleted after the model has been fitted.
  However,
  if there is an error in the model fitting function or the fitting
  model is killed (e.g., by typing control-C) then these variables will
  be left in \code{smartpredenv}.  At the beginning of model fitting,
  these variables are deleted if present in \code{smartpredenv}.


% In S-PLUS they are created in frame 1.



  During prediction, the variables
  \code{.smart.prediction} and 
  \code{.smart.prediction.counter}
  are reconstructed and read by the smart functions when the model
  frame is re-evaluated. 
  After prediction, these variables are deleted. 


  If the modelling function is used with argument \code{smart = FALSE}
  (e.g., \code{vglm(..., smart = FALSE)}) then smart prediction will not
  be used, and the results should match with the original \R functions.


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


%\note{
% In S-PLUS you will need to load in the \pkg{smartpred} library with
% the argument \code{first = T}, e.g.,
% \code{library(smartpred, lib = "./mys8libs", first = T)}.
% Here, \code{mys8libs} is the name of a directory of installed packages.
% To install the smartpred package in Linux/Unix, type something like
% \code{Splus8 INSTALL -l ./mys8libs ./smartpred_0.8-2.tar.gz}.

%}




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
  For \code{\link[VGAM]{vgam}} in \pkg{VGAM},
  \code{\link[VGAM]{sm.ps}} is important.
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


% The website \url{http://www.stat.auckland.ac.nz/~yee}
% contains more information such as how to write a
% smart function, and other technical details.


  The functions \code{\link[VGAM]{vglm}},
  \code{\link[VGAM]{vgam}},
  \code{\link[VGAM]{rrvglm}}
  and
  \code{\link[VGAM]{cqo}}
  in T. W. Yee's \pkg{VGAM}
  package are examples of modelling functions that employ smart prediction.



}
\section{WARNING }{
% In S-PLUS,
% if the \code{"bigdata"} library is loaded then it is
% \code{detach()}'ed. This is done because
% \code{scale} cannot be made smart if \code{"bigdata"} is loaded
% (it is loaded by default in the Windows version of
% Splus 8.0, but not in Linux/Unix).
% The function \code{\link[base]{search}} tells what is
% currently attached.


% In \R and S-PLUS


  The functions
  \code{\link[splines]{bs}},
  \code{\link[splines]{ns}},
  \code{\link[stats]{poly}} and
  \code{\link[base]{scale}}
  are now left alone (from 2014-05 onwards) and no longer smart.
  They work via safe prediction.
  The smart versions of these functions have been renamed and
  they begin with \code{"sm."}.




  The functions
  \code{\link[splines]{predict.bs}} and
  \code{predict.ns}
  are not smart.
  That is because they operate on objects that contain attributes only
  and do not have list components or slots.
  The function
  \code{\link[stats:poly]{predict.poly}} is not smart.






}

\examples{
# Create some data first
n <- 20
set.seed(86)  # For reproducibility of the random numbers
ldata <- data.frame(x2 = sort(runif(n)), y = sort(runif(n)))
library("splines")  # To get ns() in R

# This will work for R 1.6.0 and later
fit <- lm(y ~ ns(x2, df = 5), data = ldata)
\dontrun{
plot(y ~ x2, data = ldata)
lines(fitted(fit) ~ x2, data = ldata)
new.ldata <- data.frame(x2 = seq(0, 1, len = n))
points(predict(fit, new.ldata) ~ x2, new.ldata, type = "b", col = 2, err = -1)
}

# The following fails for R 1.6.x and later. It can be
# made to work with smart prediction provided
# ns is changed to sm.ns and scale is changed to sm.scale:
fit1 <- lm(y ~ ns(scale(x2), df = 5), data = ldata)
\dontrun{
plot(y ~ x2, data = ldata, main = "Safe prediction fails")
lines(fitted(fit1) ~ x2, data = ldata)
points(predict(fit1, new.ldata) ~ x2, new.ldata, type = "b", col = 2, err = -1)
}

# Fit the above using smart prediction
\dontrun{
library("VGAM")  # The following requires the VGAM package to be loaded 
fit2 <- vglm(y ~ sm.ns(sm.scale(x2), df = 5), uninormal, data = ldata)
fit2@smart.prediction
plot(y ~ x2, data = ldata, main = "Smart prediction")
lines(fitted(fit2) ~ x2, data = ldata)
points(predict(fit2, new.ldata, type = "response") ~ x2, data = new.ldata,
       type = "b", col = 2, err = -1)
}
}
%\keyword{smart}
\keyword{models}
\keyword{regression}
\keyword{programming}



%lm(..., smart = TRUE)
%glm(..., smart = TRUE)
%ns()
%bs()
%poly()
%scale()
%vglm(..., smart = TRUE)
%rrvglm(..., smart = TRUE)
%vgam(..., smart = TRUE)
%cao(..., smart = TRUE)
%cqo(..., smart = TRUE)
%uqo(..., smart = TRUE)

%library(smartpred, lib = "./mys8libs", first = T)





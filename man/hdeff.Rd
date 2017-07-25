\name{hdeff}
\alias{hdeff}
\alias{hdeff.vglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hauck-Donner effects: detection test for Wald tests }
\description{
  A detection test for
  Hauck-Donner effects of each regression coefficient
  in a VGLM regression model.


}
\usage{
hdeff(object, ...)
hdeff.vglm(object, derivative = NULL, se.arg = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A \code{\link{vglm}} object.
  Currently only a limited number of family functions have
  the HDE detection test:
  \code{\link{binomialff}},
  \code{\link{cumulative}},
  \code{\link{erlang}},
  \code{\link{poissonff}},
  \code{\link{topple}},
  \code{\link{uninormal}},
  \code{\link{zipoissonff}},
  and
  \code{\link{zipoisson}}.
  More will be implemented in the short future!



% \code{\link{multinomial}},



  }

  \item{derivative}{
  Numeric. Either 1 or 2.
  Currently only a few models having one linear predictor are handled
  when \code{derivative = 2}, e.g.,
  \code{\link{binomialff}},
  \code{\link{poissonff}}.
  


  }
  \item{se.arg}{
  Logical. If \code{TRUE} then the derivatives of the standard errors
  are returned as well, otherwise the derivatives are of the
  Wald statistics.


  }
  \item{\dots}{
  further arguments passed into the other methods functions.
% e.g., \code{subset}.


  }
}
\details{
Hauck and Donner (1977)
first observed an
aberration of the Wald test statistic not monotonically increasing as a
function of increasing distance between the parameter estimate and the
null value (called the Hauck-Donner effect, or HDE, here).
This "disturbing" and "undesirable" underappreciated effect has since
been observed in other regression models by various authors.  This
function computes the first, and possibly second, derivative of the Wald
statistic for each regression coefficient.  A negative value of the
first derivative is indicative of the HDE being present.



By default this function returns a labelled logical vector;
a \code{TRUE} means the HDE is affirmative for that coefficient.
  Hence ideally all values are \code{FALSE}.
  Any \code{TRUE} values suggests that the MLE is
  near the boundary of the parameter space,
  and that the p-value for that regression coefficient
  is biased upwards.



  % and  that a likelihood ratio test is recommended.



}
\value{
  By default, a vector of logicals.
  Setting \code{deriv = 1} returns a vector of first
  derivatives  of the Wald statistics.
  Setting \code{deriv = 2} returns a 2-column matrix of first
  and second derivatives of the Wald statistics.
  Setting \code{se.arg = TRUE} returns an additional 1 or 2 columns.



  For those \pkg{VGAM} family functions whose HDE test has not yet
  been implmented a \code{NULL} is returned.


}
\references{


Hauck, J. W. W. and A. Donner (1977)
Wald's test as applied to hypotheses in logit analysis. 
\emph{Journal of the American Statistical Association},
\bold{72}, 851--853.
Corrigenda: JASA, \bold{75}, 482.



% \textit{JASA 72(360): 851--3}] 75 (370), 482



Yee, T. W. (2017)
\emph{Detecting the Hauck-Donner effect in Wald tests}
(in preparation).



}
\author{ T. W. Yee.  }

\section{Warning }{
  Some 2nd derivatives are \code{NA}, meaning that they
  have not been programmed in yet.



}

\note{
  The function \code{\link{summaryvglm}} conducts the HDE
  detection test if possible and prints out a
  modified Wald table if some HDEs are detected.



  This function is currently a little beyond the
  experimental stage and may change quite a bit
  in the short future.



}
\seealso{
  \code{\link{summaryvglm}},
  \code{\link{vglm}},
  \code{\link{lrp.vglm}}.



% \code{\link{multinomial}},
% \code{\link{cumulative}},



}

\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let, data = pneumo,
            cumulative(reverse = TRUE,  parallel = TRUE))
hdeff(fit)
hdeff(fit, deriv = 1)
hdeff(fit, deriv = 2)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}





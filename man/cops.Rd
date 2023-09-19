\name{cops}
\alias{cops}
\alias{copsvglm}
\alias{cops,vglm-method}
%- Also NEED an '\alias' for EACH other topic
%- documented here.
\title{ Centre of the Parameter Space }
\description{
  Returns a vector similar to coefs() comprising
  the centre of the parameter space (COPS) values,
  given a fitted VGLM regression.


}
\usage{
cops(object, ...)
copsvglm(object, beta.range = c(-5, 6),
         tol = .Machine$double.eps^0.25,
         dointercepts = TRUE, trace. = FALSE,
         slowtrain = FALSE, ...)
}
%- maybe also 'usage' for other objects
%- documented here.
\arguments{
  \item{object}{
  A \code{\link{vglm}} object.
  However, this function will not work for all
  such objects.



  }

  \item{beta.range}{
    Numeric.
    Interval for the numerical search.
    After a little scaling, it is effectively
    fed into \code{interval} in
    \code{\link[stats]{optimize}}.
    Convergence failure may occur if this
    argument is too wide so it is a good
    idea to vary this argument.


  }
  \item{tol}{
    Numeric.
    Fed into \code{tol} in
    \code{\link[stats]{optimize}}.


  }
  \item{dointercepts}{
    Logical.
    Compute the COPS for the intercepts?
    This should be set to \code{FALSE}
    for models such as
    \code{\link{propodds}} and
    \code{\link{cumulative}}.
  

  }
  \item{trace.}{
    Logical.
    Print a running log?
    This may or may not work properly.

    

  }
  \item{slowtrain}{
    Logical.
    If \code{TRUE} then all columns of a
    matrix is computed.
    If \code{FALSE} then only one column of a
    matrix is computed, and this is the only
    column needed.
    

  }

  \item{\dots}{
    currently unused but may be used in the
    future for further arguments passed into
    the other methods functions.
% e.g., \code{subset}.


  }
}
\details{
  For many models, some COPS values will be
  \code{Inf} or \code{-Inf}
  so that manual checking is needed,
  for example, \code{\link{poissonff}}.
  Each value returned may be effectively
  that of \code{beta.range}
  or \code{NA}.
  The answers returned by this function only
  make sense if the COPSs are in the
  interior of the parameter space.

  

}
\value{
  A named vector, similar to \code{\link{coefvlm}}.
  If \code{trace.} then a list is returned,
  having a componennt comprising a
  matrix of function evaluations used by
  \code{\link[stats]{optimize}}.


}
\references{



Yee, T. W. (2021).
Some new results concerning the Hauck-Donner effect.
\emph{Manuscript in preparation}.



%Yee, T. W. (2022).
%Some new results concerning the Wald tests and
%the parameter space.
%\emph{In preparation}.




}
\author{ Thomas W. Yee.  }

%\section{Warning }{



%}

\note{


  This function is experimental and
  can be made to run more efficiently
  in the future.


}
\seealso{
  \code{\link{hdeff}}.
  \code{\link{coefvlm}},
  \code{\link[stats]{coef}}.



}

\examples{
\dontrun{data("xs.nz", package = "VGAMdata")
data1 <-
    na.omit(xs.nz[, c("age", "cancer", "sex")])
fit1 <- vglm(cancer ~ age + sex, binomialff,
             data = data1)
cops(fit1)
}}

% Add >=1 standard keywords, see file 'KEYWORDS'
% in the R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{htest}





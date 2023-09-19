\name{Influence}
\alias{Influence}
\alias{Influence.vglm}
\title{
  Influence Function (S4 generic)
  of a Fitted Model

}
\description{
  Returns a matrix containing the influence function
  of a fitted model, e.g., a "vglm" object.


}
\usage{
Influence(object, \dots)
Influence.vglm(object, weighted = TRUE, \dots)
}
\arguments{
  \item{object}{
    an object, especially that of class \code{"vglm"}---see
    \code{\link{vglm-class}}.
    Currently other classes such as \code{"vgam"} are not yet
    implemented.



  }
  \item{weighted}{
    Logical. Include the prior weights?
    Currently only \code{TRUE} is accepted.
    This might change in the future and/or
    the default value might change.


  }
  \item{\dots}{
    any additional arguments such as to
    allow or disallow the prior weights.


  }
}
\value{
    An \code{n} by \code{p.vlm} matrix.


%   \code{\link[countreg]{rootogram}};


}
\details{
  Influence functions are useful in fields such
  as sample survey theory,
  e.g.,
  \pkg{survey},
  \pkg{svyVGAM}.
  For each \eqn{i=1,\ldots,n}{i=1,...,n},
  the formula is approximately \eqn{-I U}
  where \eqn{I} is the weighted Fisher
  information matrix and U is
  the \eqn{i}th score vector.
  


  
}
% \note{
%
%


%}
\section{Warning}{
  This function is currently experimental and
  defaults may change.
  Use with caution!
  The functions here should not be confused with
  \code{\link[stats]{lm.influence}}.



}
\seealso{
  \code{\link{vglm}},
  \code{\link{vglm-class}},
  \pkg{survey}.



}
%\references{
%
%   
%}



%\author{
%  T. W. Yee.
%
%}




\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let, acat, data = pneumo)
coef(fit)  # 8-vector
Influence(fit)  # 8 x 4
all(abs(colSums(Influence(fit))) < 1e-6)  # TRUE
}
\keyword{models}





%\donttest{}
%\dontshow{utils::example("lm", echo = FALSE)}
%           (   pkg = c("countreg", "vcd"), \dots)

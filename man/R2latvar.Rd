\name{R2latvar}
\alias{R2latvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ R-squared for Latent Variable Models }
\description{
  R-squared goodness of fit for latent variable models,
  such as cumulative link models.
  Some software such as Stata call the quantity
  the McKelvey--Zavoina R-squared, which was proposed
  in their 1975 paper for cumulative probit models.

  

}
\usage{
R2latvar(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{
  A \code{\link{cumulative}} or
  \code{\link{binomialff}} fit using
  \code{\link{vglm}}.
  Only a few selected link functions are currently permitted:
  \code{\link{logitlink}},
  \code{\link{probitlink}},
  \code{\link{clogloglink}}.
  For models with more than one linear predictor,
  a parallelism assumption is needed also, i.e.,
  the constraint matrices must be a 1-column matrix of 1s
  (except for the intercept).
  The model is assumed to have an intercept term.


}
}
\details{
  Models such as the proportional odds model have
  a latent variable interpretation
  (see, e.g., Section 6.2.6 of Agresti (2018),
  Section 14.4.1.1 of Yee (2015),
  Section 5.2.2 of McCullagh and Nelder (1989)).
  It is possible to summarize the predictive power of
  the model by computing \eqn{R^2} on the transformed
  scale, e.g., on a standard normal distribution for
  a \code{\link{probitlink}} link.
  For more details see Section 6.3.7 of Agresti (2018).



}
\value{
  The \eqn{R^2} value.
  Approximately, that amount is the variability in the
  latent variable of the model explained by all the explanatory
  variables.
  Then taking the positive square-root gives an approximate
  multiple correlation \eqn{R}.



}
\references{

% Agresti, A. (2007)
% \emph{An Introduction to Categorical Data Analysis, 2nd ed.},
% New York: John Wiley & Sons.
% Page 38.



  Agresti, A. (2018)
  \emph{An Introduction to Categorical Data Analysis, 3rd ed.},
  New York: John Wiley & Sons.


  
  McKelvey, R. D. and W. Zavoina (1975)
  A statistical model for the analysis of
  ordinal level dependent variables.
  \emph{The Journal of Mathematical Sociology}, \bold{4},
  103--120.


}

\author{ Thomas W. Yee }

%\note{
%}

%\section{Warning }{
%  This


%}

\seealso{
  \code{\link{vglm}},
  \code{\link{cumulative}},
  \code{\link{propodds}},
  \code{\link{logitlink}},
  \code{\link{probitlink}},
  \code{\link{clogloglink}},
  \code{\link[stats]{summary.lm}}.


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
(fit <- vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo))
R2latvar(fit)
}
\keyword{models}
\keyword{regression}


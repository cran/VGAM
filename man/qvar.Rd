\name{qvar}
\alias{qvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quasi-variances Extraction Function

%%  ~~function to do ... ~~
}
\description{
  Takes a \code{\link{rcim}} fit of the appropriate format and
  returns either the quasi-variances or quasi-standard errors.


%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qvar(object, se = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A \code{\link{rcim}} object that has family function
  \code{\link{uninormal}} with the
  \code{\link{explink}} link.
  See below for an example.


}
\item{se}{
  Logical. If \code{FALSE} then the quasi-variances are returned,
  else the square root of them, called quasi-standard errors.


}

\item{\ldots}{
  Currently unused.


}
}
\details{

  This simple function is ad hoc and simply is equivalent to
  computing the quasi-variances
  by \code{diag(predict(fit1)[, c(TRUE, FALSE)]) / 2}.
  This function is for convenience only.
  Serious users of quasi-variances ought to understand why and how this
  function works.


}
\value{
  A vector of quasi-variances  or quasi-standard errors.


}
%\references{
%
%}

\author{

  T. W. Yee.


}

%\note{
%  This is an adaptation of \code{qvcalc()} in \pkg{qvcalc}.
%
%
%}

%\section{Warning }{
%  N
%
%
%}


\seealso{
  \code{\link{rcim}},
  \code{\link{uninormal}},
  \code{\link{explink}},
  \code{\link{Qvar}},
  \code{\link[MASS]{ships}}.


%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("ships", package = "MASS")
Shipmodel <- vglm(incidents ~ type + year + period,
                  quasipoissonff, offset = log(service),
                  data = ships, subset = (service > 0))

# Easiest form of input
fit1 <- rcim(Qvar(Shipmodel, "type"), uninormal("explink"), maxit = 99)
qvar(fit1)             # Quasi-variances
qvar(fit1, se = TRUE)  # Quasi-standard errors

# Manually compute them:
(quasiVar <- exp(diag(fitted(fit1))) / 2)                 # Version 1
(quasiVar <- diag(predict(fit1)[, c(TRUE, FALSE)]) / 2)   # Version 2
(quasiSE  <- sqrt(quasiVar))

\dontrun{ qvplot(fit1, col = "green", lwd = 3, scol = "blue", slwd = 2, las = 1) }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
% \code{\link[qvcalc:qvcalc]{qvcalc}} in \pkg{qvcalc}

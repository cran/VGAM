\name{plotdgaitd.vglm}
\alias{plotdgaitd}
\alias{plotdgaitd.vglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting the GAITD Combo Density from a GAITD Regression Object

}
\description{
  Given a GAITD regression object,
  plots the probability mass function.


  %%  ~~ A concise (1-5 lines) description of what the function does. ~~


}
\usage{
plotdgaitd(object, ...)
plotdgaitd.vglm(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A fitted GAITD combo regression, e.g.,
    \code{\link{gaitdpoisson}}.

  }
  \item{\dots}{
    Graphical arguments passed into \code{\link{dgaitdplot}}.

    
}
}
\details{
  This is meant to be a more convenient function for plotting
  the PMF of the GAITD combo model from a fitted regression model.
  The fit should be intercept-only and the distribution
  should have 1 or 2 parameters.
  Currently it should work for a \code{\link{gaitdpoisson}} fit.
  As much information as needed
  such as the special values
  is extracted from the object
  and fed into \code{\link{dgaitdplot}}.





}
\value{
  Same as \code{\link{dgaitdplot}}.

}


%\references{
%% ~put references to the literature/web site here ~
%}
\author{
T. W. Yee.
}

\note{
This function is subject to change.


  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{dgaitdplot}},
  \code{\link{spikeplot}},
  \code{\link{gaitdpoisson}}.


% \code{\link{Gaitgenpois1}}.



}
\examples{
\dontrun{
example(gaitdpoisson)
gaitpfit2 <-
  vglm(y1 ~ 1, crit = "coef", trace = TRUE, data = gdata,
       gaitdpoisson(a.mix = a.mix, i.mix = i.mix,
                    i.mlm = i.mlm, eq.ap = TRUE, eq.ip = TRUE,
                    truncate = tvec, max.support = max.support))
plotdgaitd(gaitpfit2)
}}
%\keyword{graphs}
%\keyword{models}
\keyword{regression}
\keyword{hplot}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

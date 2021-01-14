\name{meangait}
\alias{meangait}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean of the GAIT Combo Density
}
\description{
  Returns the mean of a 1- or 2-parameter GAIT combo probability mass function.


  %%  ~~ A concise (1-5 lines) description of what the function does. ~~


}
\usage{
meangait(theta.p, fam = c("pois", "log", "zeta"),
    alt.mix = NULL, inf.mix = NULL,
    alt.mlm = NULL, inf.mlm = NULL, truncate = NULL,
    max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
    pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
    theta.a = theta.p, theta.i = theta.p, deflation = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta.p}{
  Same as \code{\link{dgaitplot}};
  usually of length 1 but may be of length 2.




}
  \item{fam}{
  Same as \code{\link{dgaitplot}}.
  The default is the first one.
  All other choices are listed in that vector.


    
}
  \item{alt.mix, inf.mix, alt.mlm, inf.mlm}{
  Same as \code{\link{dgaitplot}}.


}
  \item{truncate, max.support}{
  Same as \code{\link{dgaitplot}}.

}
  \item{pobs.mix, pobs.mlm, pstr.mix, pstr.mlm, byrow.ai}{
  Same as \code{\link{dgaitplot}}.

}
\item{theta.a, theta.i, deflation}{
  Same as \code{\link{dgaitplot}}.

  
}
  \item{\dots}{
   Currently unused.



  }
}
\details{
  This function returns the mean of the PMF of
  the GAIT combo model.
  Many of its arguments are the same as \code{\link{dgaitplot}}.
  More functionality may be added in the future, such as
  returning the variance.





}
\value{
  The mean.

}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
T. W. Yee.
}

\note{
This utility function may change a lot in the future.


  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{dgaitplot}},
  \code{\link{Gaitpois}},
  \code{\link{gaitpoisson}}.



% \code{\link{Gaitgenpois1}}.



}
\examples{
i.mix <- seq(0, 15, by = 5)
lambda.p <- 10
meangait(lambda.p, alt.mix = i.mix + 1, inf.mix = i.mix,
         max.support = 17, pobs.mix = 0.1, pstr.mix = 0.1)
}
\keyword{models}
\keyword{regression}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.







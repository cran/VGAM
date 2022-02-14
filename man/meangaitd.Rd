\name{meangaitd}
\alias{meangaitd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean of the GAITD Combo Density

}
\description{
  Returns the mean of a 1- or 2-parameter GAITD combo
  probability mass function.


}
\usage{
meangaitd(theta.p, fam = c("pois", "log", "zeta"),
    a.mix = NULL, i.mix = NULL, d.mix = NULL,
    a.mlm = NULL, i.mlm = NULL, d.mlm = NULL,
    truncate = NULL, max.support = Inf,
    pobs.mix = 0, pobs.mlm = 0,
    pstr.mix = 0, pstr.mlm = 0,
    pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
    theta.a = theta.p, theta.i = theta.p, theta.d = theta.p, ...)
}
%   theta.d = theta.p, deflation = FALSE, ...
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta.p}{
  Same as \code{\link{dgaitdplot}};
  usually of length 1 but may be of length 2.




}
  \item{fam}{
  Same as \code{\link{dgaitdplot}}.
  The default is the first one.
  All other choices are listed in that vector.


    
}
  \item{a.mix, i.mix, a.mlm, i.mlm}{
  Same as \code{\link{dgaitdplot}}.


}
  \item{d.mix, d.mlm}{
  Same as \code{\link{dgaitdplot}}.


}
  \item{truncate, max.support}{
  Same as \code{\link{dgaitdplot}}.

}
  \item{pobs.mix, pobs.mlm, byrow.aid}{
  Same as \code{\link{dgaitdplot}}.

}
  \item{pstr.mix, pstr.mlm, pdip.mix, pdip.mlm}{
  Same as \code{\link{dgaitdplot}}.

}
\item{theta.a, theta.i, theta.d}{
  Same as \code{\link{dgaitdplot}}.

  
}
  \item{\dots}{
   Currently unused.



  }
}
\details{
  This function returns the mean of the PMF of
  the GAITD combo model.
  Many of its arguments are the same as \code{\link{dgaitdplot}}.
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
  \code{\link{dgaitdplot}},
  \code{\link{Gaitdpois}},
  \code{\link{gaitdpoisson}}.



% \code{\link{Gaitgenpois1}}.



}
\examples{
i.mix <- seq(0, 15, by = 5)
lambda.p <- 10
meangaitd(lambda.p, a.mix = i.mix + 1, i.mix = i.mix,
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







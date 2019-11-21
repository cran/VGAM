\name{rootogram4}
\alias{rootogram4}
\alias{rootogram4vglm}
\title{
  Rootograms (S4 generic)
  for Assessing Goodness of Fit of Probability Models

}
\description{
  A graphical technique for comparing the observed and
  fitted counts from a probability model, on a
  square root scale.


}
\usage{
rootogram4(object, \dots)
rootogram4vglm(object, newdata = NULL, breaks = NULL, max = NULL,
               xlab = NULL, main = NULL, width = NULL, \dots)
}
\arguments{
  \item{object}{
    an object of class \code{"vglm"}.
    zz This includes \code{"vgam"} because
     \code{"vlm"} handles both VGLM and VGAM objects.


%    It is strongly recommended that this be the full model
%    because a backward direction is taken first.


  }
  \item{newdata}{ Data upon which to base the calculations.
    The default is the one used to fit the model.


  }
  \item{breaks}{numeric. Breaks for the histogram intervals.}
  \item{max}{maximum count displayed.}
  \item{xlab, main}{graphical parameters.}
  \item{width}{numeric. Widths of the histogram bars.}
%  \item{pkg}{which package to call \code{rootogram()}.
%    The first is the default. }

  \item{\dots}{
    any additional arguments to
    \code{rootogram.default} and
    \code{plot.rootogram} in \pkg{countreg}.


%   \code{\link[countreg]{rootogram.default}}
%    and \code{\link[countreg]{plot.rootogram}}.


    
%   and its associated functions.


  }
}
\value{
    See
    \code{rootogram} in \pkg{countreg};
    an object of class \code{"rootogram0"}
    inheriting from \code{"data.frame"} with
    about 8 variables.


%   \code{\link[countreg]{rootogram}};


}
\details{
  Rootograms are a useful graphical technique for
  comparing the observed counts with the expected
  counts given a probability model.

  
% on \code{\link[countreg]{rootogram}}


  This S4 implementation is based very heavily
  on \code{rootogram} coming from
  \pkg{countreg}. This package is primarily written by
  A. Zeileis and
  C. Kleiber.
  That package is currently on R-Forge but not CRAN, and
  it is based on S3.
  Since \pkg{VGAM} is written using S4, it was necessary
  to define an S4 generic function called
  \code{rootogram4()} which dispatches appropriately for
  S4 objects.


% The second package is \pkg{vcd} is on CRAN and is written
% by David Meyer [aut, cre],
% Achim Zeileis ORCID iD [aut],
% Kurt Hornik [aut],
% Florian Gerber [ctb],
% Michael Friendly [ctb].


  Currently, only a selected number of \pkg{VGAM} family functions
  are implemented. Over time, hopefully more and more will be
  completed.



}
\note{
  The function names used coming from \pkg{countreg} have
  been renamed slightly to avoid conflict.

  
%  Ditto for \pkg{vcd}.



}
\section{Warning}{
  This function is rudimentary and based totally
  on the implementation in \pkg{countreg}.



%   and \pkg{vcd}.



}
\seealso{
  \code{\link{vglm}},
  \code{\link{vgam}},
  \code{\link[stats]{glm}},
  \code{\link{zipoisson}},
  \code{\link{zapoisson}},
  \code{rootogram} in \pkg{countreg}.


%  \code{rootogram} in \pkg{vcd}.



% \code{\link[countreg]{rootogram}}.



}
\references{
  Friendly, M. and Meyer, D. (2016)
  \emph{Discrete Data Analysis with R: Visualization and
    Modeling Techniques for Categorical and Count Data},
Boca Raton, FL, USA: Chapman & Hall/CRC Press.



Kleiber, C. and Zeileis, A. (2016)
  \dQuote{Visualizing Count Data Regressions Using Rootograms.} 
  \emph{The American Statistician},
  \bold{70}(3), 296--303.
  \doi{10.1080/00031305.2016.1173590}.

  
  Tukey, J. W. (1977)
  \emph{Exploratory Data Analysis},
   Reading, MA, USA: Addison-Wesley.

   
}




%\references{
%  Hastie, T. J. and Pregibon, D. (1992)
%  \emph{Generalized linear models.}
%  Chapter 6 of \emph{Statistical Models in S}
%  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

%  Venables, W. N. and Ripley, B. D. (2002)
%  \emph{Modern Applied Statistics with S.}
%  New York: Springer (4th ed).
%}



\author{
  Package \pkg{countreg} is primarily written by
  A. Zeileis and
  C. Kleiber.
  Function \code{rootogram4()} is based very heavily
  on \pkg{countreg}.


% and \pkg{vcd}.


  
%  Package \pkg{vcd} is written
%  by David Meyer [aut, cre],
%  Achim Zeileis ORCID iD [aut],
%  Kurt Hornik [aut],
%  Florian Gerber [ctb],
%  Michael Friendly [ctb].

}




\examples{
\dontrun{
data("hspider", package = "VGAM")  # Count responses
hs.p   <- vglm(Pardlugu ~ CoveHerb,   poissonff, data = hspider)
hs.nb  <- vglm(Pardlugu ~ CoveHerb, negbinomial, data = hspider)
hs.zip <- vglm(Pardlugu ~ CoveHerb,   zipoisson, data = hspider)
hs.zap <- vglm(Pardlugu ~ CoveHerb,   zapoisson, data = hspider)

opar <- par(mfrow = c(2, 2))  # Plot the rootograms
rootogram4(hs.p,   max = 15, main = "poissonff")
rootogram4(hs.nb,  max = 15, main = "negbinomial")
rootogram4(hs.zip, max = 15, main = "zipoisson")
rootogram4(hs.zap, max = 15, main = "zapoisson")
par(opar)
}
}
\keyword{models}





%\donttest{}
%\dontshow{utils::example("lm", echo = FALSE)}
%           (   pkg = c("countreg", "vcd"), \dots)

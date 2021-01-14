\name{dgaitplot}
\alias{dgaitplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting the GAIT Combo Density
}
\description{
  Plots a 1- or 2-parameter GAIT combo probability mass function.


  %%  ~~ A concise (1-5 lines) description of what the function does. ~~


}
\usage{
dgaitplot(theta.p, fam = "pois", alt.mix = NULL, inf.mix = NULL,
    alt.mlm = NULL, inf.mlm = NULL, truncate = NULL,
    max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
    pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
    theta.a = theta.p, theta.i = theta.p, deflation = FALSE,
    plot.it = TRUE, new.plot = TRUE, offset.x = ifelse(new.plot, 0, 0.25),
    type.plot = "h", xlim = c(0, min(100, max.support + 2)),
    ylim = NULL, xlab = "", ylab = "Probability", main = "",
    cex.main = 1.2, posn.main = NULL,
    lty.p = "solid", lty.a.mix = "longdash", lty.a.mlm = "longdash",
    lty.i.mix = "dashed", lty.i.mlm = "dashed", col.p = "pink2",
    col.a.mix = "#007FFF", col.a.mlm = "blue",
    col.i.mix = "#3F00FF", col.i.mlm = "purple",
    col.t = "tan", cex.p = 1, lwd.p = NULL, lwd.a = NULL, lwd.i = NULL,
    iontop = TRUE, las = 0, lend = "round",
    axes.x = TRUE, axes.y = TRUE, Plot.trunc = TRUE, cex.t = 1,
    pch.t = 1, baseparams.argnames = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta.p}{
Numeric, usually scalar but may have length 2.
This matches with, e.g., \code{lambda.p} for \code{\link{Gaitpois}}.
A length 2 example is
\code{c(mean.p, dispind.p)} for \code{Gaitgenpois1}, in
which case \code{fam = "genpois1"}.


}
  \item{fam}{
    Character, \code{paste0("dgait", fam)} should be a
    \code{d}-type function returning the PMF.
    The default is for the GAIT Poisson combo.


    
}
  \item{alt.mix, inf.mix, alt.mlm, inf.mlm}{
  See \code{\link{Gaitpois}} and \code{\link{gaitpoisson}}.


}
  \item{truncate, max.support}{
  See \code{\link{Gaitpois}} and \code{\link{gaitpoisson}}.

}
  \item{pobs.mix, pobs.mlm, pstr.mix, pstr.mlm, byrow.ai}{
  See \code{\link{Gaitpois}} and \code{\link{gaitpoisson}}.

}
\item{theta.a, theta.i}{
  Similar to \code{theta.p}, and they should have the same length
  too.

  
}
  \item{deflation}{
  Logical. Allow deflation?
}
  \item{plot.it}{
  Logical. Plot the PMF?
}
  \item{new.plot, offset.x}{
  If \code{new.plot} then \code{\link[graphics]{plot}} is called.
  If multiple plots are desired then use \code{offset.x} to
  shift the lines.

}
  \item{xlim, ylim, xlab, ylab}{
    See \code{\link[graphics]{par}} and \code{\link[graphics]{plot.default}}.
    
}
  \item{main, cex.main, posn.main}{
    Character, size and position of \code{main} for the title.
    See \code{\link[graphics]{title}},
    \code{\link[graphics]{par}} and \code{\link[graphics]{plot.default}}.
    The position is used if it is a 2-vector.



}
  \item{lty.p, lty.a.mix, lty.a.mlm, lty.i.mix, lty.i.mlm}{
    Line type for parent, altered and inflated.
    See \code{\link[graphics]{par}} and \code{\link[graphics]{plot.default}}.


    
}
  \item{col.p, col.a.mix, col.a.mlm, col.i.mix, col.i.mlm}{
    Line colour for parent (nonspecial), altered and inflated
    values.
    See \code{\link[graphics]{par}} and \code{\link[graphics]{plot.default}}.
    Roughly, by default,
    the parent is pink, the altered are azure (blue)
    and the inflated are indigo (purplish/violet).

    

    
}
  \item{col.t}{
    Point colour for truncated values, the default is \code{"tan"}.

    
}
  \item{type.plot, cex.p}{
    The former matches 'type' argument in \code{\link[graphics]{plot.default}}.
    The latter is the size of the point if \code{type.plot = "p"}
    or \code{type.plot = "b"}, etc.


    
}
  \item{lwd.p, lwd.a, lwd.i}{
    Line width for parent, altered and inflated.
    See \code{\link[graphics]{par}} and \code{\link[graphics]{plot.default}}.
    By default \code{par()\$lwd} is used for all of them.





}
  \item{las, lend}{
    See \code{\link[graphics]{par}}.



}
  \item{iontop}{
    Logical. Draw the inflated bars on top?
    The default is to draw the spikes on top, but if
    \code{FALSE} then the spikes are drawn from the bottom---this
    makes it easier to see their distribution.

    

}
  \item{axes.x, axes.y}{
    Logical. Plot axes?
    See \code{\link[graphics]{par}} and \code{\link[graphics]{plot.default}}.

    
}
  \item{Plot.trunc, cex.t, pch.t}{
    Logical. Plot the truncated values?
    If so, then specify the size and plotting character.
    See \code{\link[graphics]{par}} and \code{\link[graphics]{plot.default}}.

    
}
  \item{baseparams.argnames}{
   Character string specifying the argument name for the generic
   parameter \code{theta}, e.g.,
   \code{"lambda"} for \code{\link[VGAM]{gaitpoisson}}, 
   By appending \code{.p}, there is an argument called
   \code{lambda.p} in \code{\link[VGAM]{dgaitpois}}.
   Another example is for \code{\link[VGAM]{gaitlog}}:
   \code{"shape"} appended with \code{.p} means that
   \code{\link[VGAM]{dgaitlog}} should have an argument
   called \code{shape.p}.
   This argument is optional and increases the reliability of the
   \code{\link[base]{do.call}} call internally.



}
  \item{\dots}{
   Currently unused but there is provision for passing graphical
   arguments in in the future;
   see \code{\link[graphics]{par}}.



  }
}
\details{
  This is meant to be a crude function to plot the PMF of
  the GAIT combo model.
  Some flexibility is offered via many graphical arguments, but
  there are still many improvements that could be done.





}
\value{
  A list is returned invisibly. The components are:

  \item{x}{The integer values between the values of \code{xlim}. }

  \item{pmf.z}{The value of the PMF, by calling the
  \code{d}-type function with all the arguments fed in. }


  \item{mid.pmf}{
  For inflated values, the value where the spikes begin,
  which are at the same level as the scaled parent
  distribution.
  This is a convenient way to obtain them as it is quite
  cumbersome to compute manually.
  For non-inflated values, they are the value of \code{pmf.z}.
}
%% ...
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
  \code{\link{plotdgait}},
  \code{\link{meangait}},
  \code{\link{Gaitpois}},
  \code{\link{gaitpoisson}}.



% \code{\link{Gaitgenpois1}}.



}
\examples{
\dontrun{
i.mix <- seq(0, 25, by = 5)
mean.p <- 10 
dispind.p <- 8^2 / mean.p  # Var(Y) = dispind.p * mean.p
dgaitplot(c(mean.p, dispind.p), fam = "genpois1",
  alt.mix = i.mix + 1, inf.mix = i.mix, max.support = 33,
  pobs.mix = 0.1, pstr.mix = 0.1, lwd.p = 2, lwd.a = 2, lwd.i = 2) }
}
\keyword{graphs}
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


%   col.t = "tan", lwd.p = 1, lwd.a = 1, lwd.i = 1,


% 20201026; lots of
%   See \code{\link[graphics]{par}} and \code{\link[base]{plot}}.






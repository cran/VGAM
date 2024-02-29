\name{dgaitdplot}
\alias{dgaitdplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting the GAITD Combo Density
}
\description{
  Plots a 1- or 2-parameter GAITD combo
  probability mass function.

}
\usage{
dgaitdplot(theta.p, fam = "pois", a.mix = NULL, i.mix = NULL,
    d.mix = NULL, a.mlm = NULL, i.mlm = NULL,
    d.mlm = NULL, truncate = NULL, max.support = Inf,
    pobs.mix = 0, pobs.mlm = 0,
    pstr.mix = 0, pstr.mlm = 0,
    pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
    theta.a = theta.p, theta.i = theta.p, theta.d = theta.p,
    deflation = FALSE, plot.it = TRUE, new.plot = TRUE,
    offset.x = ifelse(new.plot, 0, 0.25), type.plot = "h",
    xlim = c(0, min(100, max.support + 2)),
    ylim = NULL, xlab = "", ylab = "Probability", main = "",
    cex.main = 1.2, posn.main = NULL,
    all.col = NULL, all.lty = NULL, all.lwd = NULL, 
    lty.p = "solid", lty.a.mix = "longdash", lty.a.mlm = "longdash",
    lty.i.mix = "dashed", lty.i.mlm = "dashed",
    lty.d.mix = "solid", lty.d.mlm = "solid", lty.d.dip = "dashed",
    col.p = "pink2",
    col.a.mix = artichoke.col, col.a.mlm = asparagus.col,
    col.i.mix = indigo.col, col.i.mlm = iris.col,
    col.d.mix = deer.col, col.d.mlm = dirt.col, col.d.dip = desire.col,
    col.t = turquoise.col, cex.p = 1, lwd.p = NULL, lwd.a = NULL,
    lwd.i = NULL, lwd.d = NULL, iontop = TRUE, dontop = TRUE,
    las = 0, lend = "round", axes.x = TRUE, axes.y = TRUE,
    Plot.trunc = TRUE, cex.t = 1, pch.t = 1,
    baseparams.argnames = NULL, nparams = 1, flip.args = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta.p}{
Numeric, usually scalar but may have length 2.
This matches with, e.g., \code{lambda.p} for
\code{\link{Gaitdpois}}.
A length 2 example is \code{c(size.p, munb.p)} for
\code{\link{Gaitdnbinom}},
in which case \code{fam = "nbinom"}.
Another length 2 example is
\code{c(mean.p, dispind.p)} for \code{Gaitgenpois1},
in which case \code{fam = "genpois1"}.



}
  \item{fam}{
    Character, \code{paste0("dgait", fam)}
    should be a \code{d}-type function returning
    the PMF.  The default is for the GAITD
    Poisson combo.



%   Setting \code{fam = "nbinom"}
%   will handle the NBD.


%   and \code{nparams} and \code{flip.args}


    
}
  \item{a.mix, i.mix, a.mlm, i.mlm}{
  See \code{\link{Gaitdpois}} and \code{\link{gaitdpoisson}}.


}
  \item{d.mix, d.mlm}{
  See \code{\link{Gaitdpois}} and \code{\link{gaitdpoisson}}.


}
  \item{truncate, max.support}{
  See \code{\link{Gaitdpois}} and \code{\link{gaitdpoisson}}.

}
  \item{pobs.mix, pobs.mlm, byrow.aid}{
  See \code{\link{Gaitdpois}} and \code{\link{gaitdpoisson}}.

}
  \item{pstr.mix, pstr.mlm, pdip.mix, pdip.mlm}{
  See \code{\link{Gaitdpois}} and \code{\link{gaitdpoisson}}.

}
\item{theta.a, theta.i, theta.d}{
  Similar to \code{theta.p}, and they should
  have the same length too.

  
}
  \item{deflation}{
    Logical.
    Plot the deflation (dip) probabilities?


}
  \item{plot.it}{
  Logical. Plot the PMF?
}
  \item{new.plot, offset.x}{
  If \code{new.plot} then
  \code{\link[graphics]{plot}} is called.
  If multiple plots are desired then use
  \code{offset.x} to shift the lines.

}
\item{xlim, ylim, xlab, ylab}{
  See \code{\link[graphics]{par}} and
  \code{\link[graphics]{plot.default}}.
  Argument \code{xlim} should be integer-valued.
 

}
  \item{main, cex.main, posn.main}{
  Character, size and position of \code{main} for the title.
  See \code{\link[graphics]{title}},
  \code{\link[graphics]{par}} and
  \code{\link[graphics]{plot.default}}.
  The position is used if it is a 2-vector.



}
\item{all.col, all.lty, all.lwd}{
  These arguments allow all the colours,
  line types and line widths arguments to be
  assigned to these values, i.e., so that they
  are the same for all values of the support.
  For example, if \code{all.lwd = 2} then this
  sets \code{lwd.p}, \code{lwd.a}, \code{lwd.i}
  and \code{lwd.d} all equal to 2.


  }
  \item{lty.p, lty.a.mix, lty.a.mlm, lty.i.mix, lty.i.mlm}{
  Line type for parent, altered and inflated.
  See \code{\link[graphics]{par}} and
  \code{\link[graphics]{plot.default}}.


    
}
  \item{col.p, col.a.mix, col.a.mlm, col.i.mix, col.i.mlm}{
    Line colour for parent (nonspecial), altered, inflated,
    truncated and deflated values.
    See \code{\link[graphics]{par}} and
    \code{\link[graphics]{plot.default}}.
    Roughly, by default and currently,
    the parent is pink-like,
    the altered are greenish,
    the inflated are purplish/violet,
    the truncated are light blue,
    and the deflated are brownish with the dip
    probabilities being reddish.
    The proper colour names are similar to
    being acrostic.  For each operator, the
    colours of \code{"mix"} vs \code{"mlm"} are
    similar but different---this is intentional.
    Warning: the default colours might change,
    depending on style!

    
}
  \item{lty.d.mix, lty.d.mlm, lty.d.dip}{
    Similar to above.
    Used when \code{deflation = TRUE}.


}
  \item{col.d.mix, col.d.mlm, col.d.dip}{
    Similar to above.
    Used when \code{deflation = TRUE}.
    The website \url{https://www.spycolor.com} was used to
    choose some of the default colours;
    the first two are also called \code{"dirt"}
    and \code{"deer"} respectively, which are
    both brownish.



}
  \item{col.t}{
    Point colour for truncated values, the default is
    \code{"tan"}.

    
}
  \item{type.plot, cex.p}{
    The former matches 'type' argument in
    \code{\link[graphics]{plot.default}}.
    The latter is the size of the point if
    \code{type.plot = "p"} or \code{type.plot = "b"},
    etc.


    
}
  \item{lwd.p, lwd.a, lwd.i, lwd.d}{
   Line width for parent, altered and inflated.
   See \code{\link[graphics]{par}} and
   \code{\link[graphics]{plot.default}}.
   By default \code{par()\$lwd} is used for all of them.


}
  \item{las, lend}{
    See \code{\link[graphics]{par}}.



}
  \item{iontop, dontop}{
    Logicals.
    Draw the inflated and deflated bars on top?
    The default is to draw the spikes on top, but if
    \code{FALSE} then the spikes are drawn
    from the bottom---this makes it easier
    to see their distribution.  Likewise, if
    \code{deflation = TRUE} then \code{dontop}
    is used to position the deflation (dip)
    probabilities.

    

}
  \item{axes.x, axes.y}{
    Logical. Plot axes?
    See \code{\link[graphics]{par}} and
    \code{\link[graphics]{plot.default}}.

    
}
  \item{Plot.trunc, cex.t, pch.t}{
    Logical. Plot the truncated values?
    If so, then specify the size and plotting character.
    See \code{\link[graphics]{par}} and
    \code{\link[graphics]{plot.default}}.

    
}
  \item{baseparams.argnames}{
   Character string specifying the argument name for the generic
   parameter \code{theta}, e.g.,
   \code{"lambda"} for \code{\link[VGAM]{gaitdpoisson}}, 
   By appending \code{.p}, there is an argument called
   \code{lambda.p} in \code{\link[VGAM]{dgaitdpois}}.
   Another example is for \code{\link[VGAM]{gaitdlog}}:
   \code{"shape"} appended with \code{.p} means that
   \code{\link[VGAM]{dgaitdlog}} should have an argument
   called \code{shape.p}.
   This argument is optional and increases the reliability of the
   \code{\link[base]{do.call}} call internally.



  }
  \item{nparams, flip.args}{
    Not for use by the user. It is used
    internally to handle the NBD.

  }
  \item{\dots}{
   Currently unused but there is provision for
   passing graphical arguments in in the future;
   see \code{\link[graphics]{par}}.



  }
}
\details{
  This is meant to be a crude function to
  plot the PMF of the GAITD combo model.
  Some flexibility is offered via many
  graphical arguments, but there are still many
  improvements that could be done.





}
\value{
  A list is returned invisibly. The components
  are:

  \item{x}{The integer values between the values
  of \code{xlim}.

  }

  \item{pmf.z}{The value of the PMF, by
  calling the \code{d}-type function with all
  the arguments fed in.

  }


  \item{sc.parent}{
  The same level as the scaled parent
  distribution.  Thus for inflated values,
  the value where the spikes begin.  And for
  deflated values, the value at the top of
  the dips.  This is a convenient way to obtain
  them as it is quite cumbersome to compute
  them manually.  For any nonspecial value,
  such as non-inflated and non-deflated values,
  they are equal to \code{pmf.z}.



% 20211008 was called mid.pmf

}
  \item{unsc.parent}{
  Unscaled parent distribution.
  If there is no alteration, inflation,
  deflation and truncation then this is
  the basic PMF stipulated by the parent
  distribution only.  Usually this is FYI only.



% 20211027; renamed to sc.parent and unsc.parent.

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
  Because this function is called by a \pkg{shiny} app,
  if any parameter values lie outside the 
  parameter space then \code{\link[base:stop]{stop}}
  will be called.
  For example, too much deflation results in
  \code{NaN} values returned by
  \code{\link{dgaitdnbinom}}.


% 20240216; parameter space info added.

  
  
}

%% Make other sections like Warning with \section{Warning }{....}

\seealso{
  \code{\link{plotdgaitd}},
  \code{\link{spikeplot}},
  \code{\link{meangaitd}},
  \code{\link{Gaitdpois}},
  \code{\link{gaitdpoisson}},
  \code{\link{Gaitdnbinom}},
  \code{\link{multilogitlink}}.



% \code{\link{Gaitgenpois1}}.



}
\examples{
\dontrun{  i.mix <- seq(0, 25, by = 5)
mean.p <- 10; size.p <- 8
dgaitdplot(c(size.p, mean.p), fam = "nbinom", xlim = c(0, 25),
     a.mix = i.mix + 1, i.mix = i.mix, pobs.mix = 0.1,
     pstr.mix = 0.1, lwd.i = 2,lwd.p = 2, lwd.a = 2)
}}
%\keyword{graphs}
%\keyword{models}
%\keyword{regression}
\keyword{dplot}
\keyword{hplot}
\keyword{distribution}



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
% See \code{\link[graphics]{par}} and \code{\link[base]{plot}}.






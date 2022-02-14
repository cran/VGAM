\name{spikeplot}
\alias{spikeplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spike Plot

}
\description{
  Produces a spike plot of a numeric vector.


}
\usage{
spikeplot(x, freq = FALSE, as.table = FALSE, col = par("col"),
    lty = par("lty"), lwd = par("lwd"), lend = par("lend"),
    type = "h", xlab = deparse1(substitute(x)), ylab = NULL,
    capped = FALSE, cex = sqrt(lwd) / 2, pch = 19, pcol = col, scol = NULL,
    slty = NULL, slwd = NULL, new.plot = TRUE, offset.x = 0, ymux = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Numeric, passed into \code{\link[base]{table}}.

  }
  \item{freq}{
    Logical. If \code{TRUE} then the y-axis measures
    the frequencies, else the sample proportions.
    Intended to be as \code{\link[graphics]{hist}}.

  }
  \item{as.table}{
    Logical.
    If \code{TRUE} then the call to \code{\link[base]{plot}} is
    closer to \code{plot(table(x), ...)}, meaning the labelling
    differs from \code{as.table = FALSE}.
    The default is to convert \code{table(x)} into a numeric
    vector which is then passed into \code{\link[base]{plot}}
    so that the labelling is more uniform along the x-axis.


  
  }
  \item{col, type, lty, lwd}{
  See \code{\link[graphics]{par}}.


  }
  \item{lend, xlab, ylab}{
  See \code{\link[graphics]{par}}.


  }
  \item{capped, cex, pch, pcol}{
    First argument is logical.
    If \code{TRUE} then the others argument are used to place
    points at the top using \code{\link[graphics]{points}}
    with \code{pcol} being its colour.  
    See \code{\link[graphics]{par}}.


  }
  \item{scol, slty, slwd}{
    Similar to \code{col}, \code{lty} and \code{lwd} but
    apply to some selected values.
    The input may be a named list such as
    \code{scol = list("green" = c(2, 4, 6), "blue" = 5)},
    \code{slty = list("dashed" = c(2, 4, 6), "dotted" = 5)},
    \code{slwd = list("2" = c(2, 4, 6), "3" = 5)},
    else a named vector such as
  \code{scol = c("green" = 2, "green" = 4, "green" = 6, "blue" = 5)},
  \code{slty = c("dashed" = 2, "dashed" = 4, "dashed" = 6, "dotted" = 5)},
  \code{slwd = c("2" = 2, "2" = 4, "2" = 6, "3" = 5)}.
  The three arguments are ignored if \code{as.table = TRUE}.


  
%  Not true since the code handles both cases:
%  Note: \code{names(slty)} should be character rather than numeric.
%  \code{slty = list("2" = c(2, 4, 6), "3" = 5)} is now handled.




}
  \item{new.plot, offset.x}{
   Logical and numeric.
   Add to an existing plot? If so, set \code{new.plot = FALSE}
   and it is useful for
   the spikes to be shifted by some amount \code{offset.x}. 


}
  \item{ymux}{
  Numeric, y-multiplier. The response is multiplied by \code{ymux}.
  This can be useful when plotting subsets side-by-side so that
  the constituent proportions add up to the overall proportion.


}
  \item{\dots}{
    Additional graphical arguments passed into an ordinary
    \code{\link[base]{plot}}, for example,
    \code{xlim}, \code{las}, \code{main}.


}
}
\details{
  \emph{Heaping} is a very commonly occurring phenomenon in
  retrospective self-reported survey data.
  Also known as \emph{digit preference} data,
  it is often characterized by an excess of multiples of 10 or 5
  upon rounding.
  For this type of data
  this simple function is meant to be convenient for
  plotting the frequencies or sample proportions of
  a vector \code{x} representing a discrete random variable.
  This type of plot
  is known as a \emph{spike plot} in STATA circles.
  If \code{table(x)} works then this function should hopefully
  work.
  The default for \code{type} means that any heaping and
  \emph{seeping} should easily be seen. 
  If such features exist then \emph{GAITD regression} is
  potentially useful---see \code{\link{gaitdpoisson}} etc.
  Currently missing values are ignored totally because
  \code{table(x)} is used without further arguments;
  this might change in the future.





}
\value{
  Returns invisibly \code{table(x)}.

}


%\references{
%% ~put references to the literature/web site here ~
%}

\author{
T. W. Yee.
}

%\note{
%This function is subject to change.


  
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~


%\section{Warning }{This function needs R 4.0.0 or later, as it uses
%  \code{\link[base]{deparse1}}.
%
%
% This is no longer true as deparse1 is now included inside spikeplot().
% yettodo: remove it when DESCRIPTION has R 4.0.0 or above is needed.
%
%
%
%
%}


\seealso{
  \code{\link[base]{table}},
  \code{\link[base]{plot}},
  \code{\link[graphics]{par}},
  \code{\link[base]{deparse1}},
  \code{\link{dgaitdplot}},
  \code{\link{plotdgaitd}},
  \code{\link{gaitdpoisson}}.


% \code{\link{Gaitgenpois1}}.



}
\examples{
\dontrun{
spikeplot(with(marital.nz, age), col = "pink2", lwd = 2)
}}
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
%  # data(marital.nz, package = "VGAMdata")


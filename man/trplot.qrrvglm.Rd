\name{trplot.qrrvglm}
\alias{trplot.qrrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Trajectory plot for QRR-VGLMs }
\description{
Produces a trajectory plot for
\emph{quadratic reduced-rank vector generalized linear models}
(QRR-VGLMs).
It is only applicable for rank-1 models with argument
\code{noRRR = ~ 1}.

}
\usage{
trplot.qrrvglm(object, which.species = NULL, add = FALSE, show.plot = TRUE,
               label.sites = FALSE, sitenames = rownames(object@y),
               axes.equal = TRUE, cex = par()$cex,
               col = 1:(nos * (nos - 1)/2), log = "",
               lty = rep_len(par()$lty, nos * (nos - 1)/2),
               lwd = rep_len(par()$lwd, nos * (nos - 1)/2),
               tcol = rep_len(par()$col, nos * (nos - 1)/2),
               xlab = NULL, ylab = NULL,
               main = "", type = "b", check.ok = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class \code{"qrrvglm"}, i.e., a CQO object. }
  \item{which.species}{ Integer or character vector specifying the
  species to be plotted. If integer, these are the columns of the
  response matrix. If character, these must match exactly with the
  species' names.
  The default is to use all species.


  }
  \item{add}{ Logical. Add to an existing plot? If \code{FALSE} (default),
  a new plot is made. }
  \item{show.plot}{ Logical. Plot it? }
  \item{label.sites}{ Logical. If \code{TRUE}, the points on the
  curves/trajectories are labelled with the \code{sitenames}. }
  \item{sitenames}{ Character vector. The names of the sites. }
  \item{axes.equal}{ Logical. If \code{TRUE}, the x- and y-axes
  will be on the same scale.


  }
  \item{cex}{ Character expansion of the labelling of the site names.
  Used only if \code{label.sites} is \code{TRUE}.
  See the \code{cex} argument in \code{\link[graphics]{par}}.


  }
  \item{col}{Color of the lines.
  See the \code{col} argument in \code{\link[graphics]{par}}.
  Here, \code{nos} is the number of species.


  }
  \item{log}{ Character, specifying which (if any) of the x- and
  y-axes are to be on a logarithmic scale.
  See the \code{log} argument in \code{\link[graphics]{par}}.


  }
  \item{lty}{  Line type.
  See the \code{lty} argument of \code{\link[graphics]{par}}.


  }
  \item{lwd}{  Line width.
  See the \code{lwd} argument of \code{\link[graphics]{par}}.


  }
  \item{tcol}{Color of the text for the site names.
  See the \code{col} argument in \code{\link[graphics]{par}}.
  Used only if \code{label.sites} is \code{TRUE}.


  }
  \item{xlab}{Character caption for the x-axis.
  By default, a suitable caption is found.
  See the \code{xlab} argument in \code{\link[graphics]{plot}}
  or \code{\link[graphics]{title}}.


  }
  \item{ylab}{Character caption for the y-axis.
  By default, a suitable caption is found.
  See the \code{xlab} argument in \code{\link[graphics]{plot}}
  or \code{\link[graphics]{title}}.


  }
  \item{main}{ Character, giving the title of the plot.
  See the \code{main} argument in \code{\link[graphics]{plot}}
  or \code{\link[graphics]{title}}.


  }
  \item{type}{ Character, giving the type of plot. A common
  option is to use \code{type="l"} for lines only.
  See the \code{type} argument of \code{\link[graphics]{plot}}.


  }
  \item{check.ok}{ Logical. Whether a check is performed to see
  that \code{noRRR = ~ 1} was used.
  It doesn't make sense to have a trace plot unless this is so.


  }
  \item{\dots}{ Arguments passed into the \code{plot} function
  when setting up the entire plot. Useful arguments here include
  \code{xlim} and \code{ylim}.


  }
}
\details{
 A trajectory plot plots the fitted values of a `second' species
 against a `first' species. The argument \code{which.species} must
 therefore contain at least two species. By default, all of the
 species that were fitted in \code{object} are plotted.
 With more than a few species
 the resulting plot will be very congested, and so it is recommended
 that only a few species be selected for plotting.


In the above, \eqn{M} is the number of species selected for plotting,
so there will be \eqn{M(M-1)/2}{M*(M-1)/2} curves/trajectories in total.


A trajectory plot will be fitted only if \code{noRRR = ~ 1} because
otherwise the trajectory will not be a smooth function of the latent
variables.


}
\value{
  A list with the following components.
  \item{species.names}{
   A matrix of characters giving the `first' and `second' species. The
   number of different combinations of species is given by the number
   of rows.  This is useful for creating a legend.

   }
  \item{sitenames}{A character vector of site names, sorted by
   the latent variable (from low to high).


  }
}
\references{


Yee, T. W. (2012)
On constrained and unconstrained
quadratic ordination.
\emph{Manuscript in preparation}.


}

\author{ Thomas W. Yee }

\note{
   Plotting the axes on a log scale is often a good idea.
   The use of \code{xlim} and \code{ylim} to control the axis limits
   is also a good idea, so as to limit the extent of the curves at low
   abundances or probabilities.
   Setting \code{label.sites = TRUE} is a good idea only if the number of
   sites is small, otherwise there is too much clutter.


}

\seealso{
  \code{\link{cqo}},
  \code{\link[graphics]{par}},
  \code{\link[graphics]{title}}.


}

\examples{\dontrun{ set.seed(111)  # This leads to the global solution
# hspider[,1:6] <- scale(hspider[,1:6])  # Standardize the environmental variables
p1 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi, Auloalbi,
                Pardlugu, Pardmont, Pardnigr, Pardpull, Trocterr, Zoraspin) ~
          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
          poissonff, data = hspider, trace = FALSE)

trplot(p1, which.species = 1:3, log = "xy", type = "b", lty = 1,
       main = "Trajectory plot of three hunting spiders species",
       col = c("blue","red","green"), lwd = 2, label = TRUE) -> ii
legend(0.00005, 0.3, lwd = 2, lty = 1, col = c("blue", "red", "green"),
       with(ii, paste(species.names[,1], species.names[,2], sep = " and ")))
abline(a = 0, b = 1, lty = "dashed", col = "grey")  # Useful reference line
}
}
\keyword{models}
\keyword{regression}
\keyword{graphs}

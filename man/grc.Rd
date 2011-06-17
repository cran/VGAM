\name{grc}
\alias{grc}
\alias{rcam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Row-Column Association Models including Goodman's RC Association Model }
\description{
  Fits a Goodman's RC association model to a matrix of counts,
  and more generally, a sub-class of row-column association models.

}
\usage{
grc(y, Rank = 1, Index.corner = 2:(1 + Rank),
    szero = 1, summary.arg = FALSE, h.step = 1e-04, ...)
rcam(y, family = poissonff, Rank = 0, Musual = NULL,
     Index.corner = if (!Rank) NULL else 1 + Musual * (1:Rank),
     rprefix = "Row.", cprefix = "Col.",
     szero = if (!Rank) NULL else {
       if (Musual == 1) 1 else setdiff(1:(Musual * ncol(y)),
                    c(1 + (1:ncol(y)) * Musual, Index.corner))
     },
     summary.arg = FALSE, h.step = 0.0001,
     rbaseline = 1, cbaseline = 1, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{
  For \code{grc} a matrix of counts.
  For \code{rcam} a general matrix response depending on \code{family}.
  Output from \code{table()} is acceptable; it is converted into a matrix.
  Note that \code{y} should be at least 3 by 3 in dimension.

  }
  \item{family}{
  A \pkg{VGAM} family function.
  The first linear/additive predictor is fitted using main effects plus
  an optional rank-\code{Rank} interaction term.
  Not all family functions are suitable or make sense.
  All other linear/additive predictors are fitted using an intercept-only,
  so it has a common value over all rows and columns.
  For example,
  \code{\link{zipoissonff}} may be suitable for counts but not
  \code{\link{zipoisson}} because of the ordering of the
  linear/additive predictors.
  If the \pkg{VGAM} family function does not have an \code{infos}
  slot then \code{Musual} needs to be inputted (the number of
  linear predictors for an ordinary (usually univariate) response,
  aka \eqn{M}).
  The \pkg{VGAM} family function also needs to be able to
  handle multiple responses; and not all of them can do this.


  }
  \item{Rank}{
  An integer from the set
  \{0,\ldots,\code{min(nrow(y), ncol(y))}\}.
  This is the dimension of the fit in terms of the interaction.
  For \code{grc()} this argument must be positive.
  A value of 0 means no interactions (i.e., main effects only);
  each row and column is represented by an indicator variable.

  }
  \item{Index.corner}{
  A vector of \code{Rank} integers.
  These are used to store the \code{Rank} by \code{Rank}
  identity matrix in the
  \code{A} matrix; corner constraints are used.


  }
  \item{rprefix, cprefix}{ 
  Character, for rows and columns resp.
  For labelling the indicator variables.

  }
  \item{szero}{ 
  An integer from the set \{1,\ldots,\code{min(nrow(y), ncol(y))}\},
  specifying the row that is used as the structural zero.

  }
  \item{summary.arg}{
  Logical. If \code{TRUE}, a summary is returned.
  If \code{TRUE}, \code{y} may be the output (fitted
  object) of \code{grc()}.

  }
  \item{h.step}{
  A small positive value that is passed into
  \code{summary.rrvglm()}. Only used when \code{summary.arg = TRUE}. }
  \item{\dots}{ Arguments that are passed into \code{rrvglm.control()}.

  }
  \item{Musual}{
  The number of linear predictors of the \pkg{VGAM} \code{family} function
  for an ordinary (univariate) response.
  Then the number of linear predictors of the \code{rcam()} fit is
  usually the number of columns of \code{y} multiplied by \code{Musual}.
  The default is to evaluate the \code{infos} slot of the
  \pkg{VGAM} \code{family} function to try to evaluate it;
  see \code{\link{vglmff-class}}.
  If this information is not yet supplied by the family function then
  the value needs to be inputted manually using this argument.


  }
  \item{rbaseline, cbaseline}{
  Baseline reference levels for the rows and columns.
  Currently stored on the object but not used.

  }
}
\details{
  Goodman's RC association model fits a reduced-rank approximation
  to a table of counts.  The log of each cell mean is decomposed as an
  intercept plus a row effect plus a column effect plus a reduced-rank
  component. The latter can be collectively written \code{A \%*\% t(C)},
  the product of two `thin' matrices.
  Indeed, \code{A} and \code{C} have \code{Rank} columns.
By default, the first column and row of the interaction matrix
\code{A \%*\% t(C)} is chosen 
to be structural zeros, because \code{szero = 1}.
This means the first row of \code{A} are all zeros.   


This function uses \code{options()$contrasts} to set up the row and 
column indicator variables.
In particular, Equation (4.5) of Yee and Hastie (2003) is used.
These are called \code{Row.} and \code{Col.} (by default) followed
by the row or column number.


The function \code{rcam()} is more general than \code{grc()}.
Its default is a no-interaction model of \code{grc()}, i.e.,
rank-0 and a Poisson distribution. This means that each
row and column has a dummy variable associated with it.
The first row and column is baseline.
The power of \code{rcam()} is that many \pkg{VGAM} family functions
can be assigned to its \code{family} argument.
For example, 
\code{\link{normal1}} fits something in between a 2-way
ANOVA with and without interactions,
\code{\link{alaplace2}} with \code{Rank = 0} is something like
\code{\link[stats]{medpolish}}.
Others include
\code{\link{zipoissonff}},
\code{\link{negbinomial}}.
Hopefully one day \emph{all} \pkg{VGAM} family functions will
work when assigned to the \code{family} argument, although the
result may not have meaning.


}
\value{
  An object of class \code{"grc"}, which currently is the same as
  an \code{"rrvglm"} object.
  Currently,
  a rank-0 \code{rcam()} object is of class \code{\link{rcam0-class}},
  else of class \code{"rcam"} (this may change in the future).

% Currently,
% a rank-0 \code{rcam()} object is of class \code{\link{vglm-class}},
% but it may become of class \code{"rcam"} one day.


}
\references{
Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


Goodman, L. A. (1981)
Association models and canonical correlation in the analysis
of cross-classifications having ordered categories.
\emph{Journal of the American Statistical Association},
\bold{76}, 320--334.


Documentation accompanying the \pkg{VGAM} package at 
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information about the setting up of the
indicator variables.

 
}
\author{
Thomas W. Yee, with
assistance from Alfian F. Hadi.


}
\note{
  These functions set up the indicator variables etc. before calling
  \code{\link{rrvglm}}
  or
  \code{\link{vglm}}.
  The \code{...} is passed into \code{\link{rrvglm.control}} or
  \code{\link{vglm.control}},
  This means, e.g., \code{Rank = 1} is default for \code{grc()}.


  The data should be labelled with \code{\link[base]{rownames}} and
  \code{\link[base]{colnames}}.
  Setting \code{trace = TRUE} is recommended for monitoring
  convergence.
  Using \code{criterion = "coefficients"} can result in slow convergence.


  If \code{summary = TRUE}, then \code{y} can be a
  \code{"grc"} object, in which case a summary can be
  returned. That is, \code{grc(y, summary = TRUE)} is
  equivalent to \code{summary(grc(y))}.


}

\section{Warning}{
  The function \code{rcam()} is experimental at this stage and
  may have bugs.
  Quite a lot of expertise is needed when fitting and in its
  interpretion thereof. For example, the constraint
  matrices applies the reduced-rank regression to the first linear
  predictor and the other linear predictors are intercept-only and
  have a common value throughout the entire data set.
  This means that \code{family =} \code{\link{zipoissonff}} is
  appropriate but not
  \code{family =} \code{\link{zipoisson}}.
  To understand what is going on, do examine the constraint
  matrices of the fitted object, and reconcile this with Equations
  (4.3) to (4.5) of Yee and Hastie (2003).


  The functions temporarily create a permanent data frame
  called \code{.grc.df} or \code{.rcam.df}, which used
  to be needed by \code{summary.rrvglm()}.  Then these
  data frames are deleted before exiting the function.
  If an error occurs, then the data frames may be present
  in the workspace.



}

\seealso{
  \code{\link{rrvglm}},
  \code{\link{rrvglm.control}},
  \code{\link{rrvglm-class}},
  \code{summary.grc},
  \code{\link{moffset}},
  \code{\link{Rcam}},
  \code{\link{plotrcam0}},
  \code{\link{alcoff}},
  \code{\link{crashi}},
  \code{\link{auuc}},
  \code{\link{olympic}},
  \code{\link{poissonff}}.

}

\examples{
# Some undergraduate student enrolments at the University of Auckland in 1990
grc1 <- grc(auuc)
fitted(grc1)
summary(grc1)

grc2 <- grc(auuc, Rank = 2, Index.corner = c(2, 5))
fitted(grc2)
summary(grc2)


# 2008 Summer Olympic Games in Beijing
top10 <- head(olympic, n = 10)
oly1 <- with(top10, grc(cbind(gold, silver, bronze)))
round(fitted(oly1))
round(resid(oly1, type = "response"), dig = 1)  # Response residuals
summary(oly1)
Coef(oly1)


# Roughly median polish
rcam0 <- rcam(auuc, fam = alaplace2(tau = 0.5, intparloc = TRUE), trace = TRUE)
round(fitted(rcam0), dig = 0)
round(100 * (fitted(rcam0) - auuc) / auuc, dig = 0)   # Discrepancy
rcam0@y
round(coef(rcam0, matrix = TRUE), dig = 2)
print(Coef(rcam0, matrix = TRUE), dig = 3)
# constraints(rcam0)
names(constraints(rcam0))

# Compare with medpolish():
(med.a <- medpolish(auuc))
fv <- med.a$overall + outer(med.a$row, med.a$col, "+")
round(100 * (fitted(rcam0) - fv) / fv)   # Hopefully should be all 0s
}
\keyword{models}
\keyword{regression}
% plot(oly1)
% oly2 <- with(top10, grc(cbind(gold,silver,bronze), Rank = 2)) # Saturated model
% round(fitted(oly2))
% round(fitted(oly2)) - with(top10, cbind(gold,silver,bronze))
% summary(oly2) # Saturated model



% zz 20100927 unsure
% Then \code{.grc.df} is deleted before exiting the function.





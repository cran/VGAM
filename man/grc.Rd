\name{grc}
\alias{grc}
\alias{rcim}
\alias{uqo}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Row-Column Interaction Models including Goodman's RC Association
  Model and Unconstrained Quadratic Ordination }
\description{
  Fits a Goodman's RC association model (GRC) to a matrix of counts,
  and more generally, row-column interaction models (RCIMs).
  RCIMs allow for unconstrained quadratic ordination (UQO).



}
\usage{
grc(y, Rank = 1, Index.corner = 2:(1 + Rank),
    str0 = 1, summary.arg = FALSE, h.step = 1e-04, ...)
rcim(y, family = poissonff, Rank = 0, M1 = NULL,
     weights = NULL, which.linpred = 1,
     Index.corner = ifelse(is.null(str0), 0, max(str0)) + 1:Rank,
     rprefix = "Row.", cprefix = "Col.", iprefix = "X2.",
     offset = 0, str0 = if (Rank) 1 else NULL,
     summary.arg = FALSE, h.step = 0.0001,
     rbaseline = 1, cbaseline = 1,
     has.intercept = TRUE, M = NULL,
     rindex = 2:nrow(y), cindex = 2:ncol(y), iindex = 2:nrow(y), ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{
  For \code{grc()}: a matrix of counts.
  For \code{rcim()}: a general matrix response depending on \code{family}.
  Output from \code{table()} is acceptable; it is converted into a matrix.
  Note that \code{y} should be at least 3 by 3 in dimension.


  }
  \item{family}{
  A \pkg{VGAM} family function.
  By default, the first linear/additive predictor
  is fitted
  using main effects plus an optional rank-\code{Rank}
  interaction term.
  Not all family functions are suitable or make sense.
  All other linear/additive predictors are fitted using an intercept-only,
  so it has a common value over all rows and columns.
  For example,
  \code{\link{zipoissonff}} may be suitable for counts but not
  \code{\link{zipoisson}} because of the ordering of the
  linear/additive predictors.
  If the \pkg{VGAM} family function does not have an \code{infos}
  slot then \code{M1} needs to be inputted (the number of
  linear predictors for an ordinary (usually univariate) response,
  aka \eqn{M}).
  The \pkg{VGAM} family function also needs to be able to
  handle multiple responses (currently not all of them can do this).


  }
  \item{Rank}{
  An integer from the set
  \{0,\ldots,\code{min(nrow(y), ncol(y))}\}.
  This is the dimension of the fit in terms of the interaction.
  For \code{grc()} this argument must be positive.
  A value of 0 means no interactions (i.e., main effects only);
  each row and column is represented by an indicator variable.


  }
  \item{weights}{
  Prior weights. Fed into
  \code{\link{rrvglm}}
  or
  \code{\link{vglm}}.


  }
  \item{which.linpred}{
  Single integer.
  Specifies which linear predictor is modelled as the sum of an
  intercept, row effect, column effect plus an optional interaction
  term. It should be one value from the set \code{1:M1}.


  }
  \item{Index.corner}{
  A vector of \code{Rank} integers.
  These are used to store the \code{Rank} by \code{Rank}
  identity matrix in the
  \code{A} matrix; corner constraints are used.


  }
  \item{rprefix, cprefix, iprefix}{
  Character, for rows and columns and interactions respectively.
  For labelling the indicator variables.


  }
  \item{offset}{
  Numeric. Either a matrix of the right dimension, else
  a single numeric expanded into such a matrix.


  }
  \item{str0}{
  Ignored if \code{Rank = 0}, else
  an integer from the set \{1,\ldots,\code{min(nrow(y), ncol(y))}\},
  specifying the row that is used as the structural zero.
  Passed into \code{\link{rrvglm.control}} if \code{Rank > 0}.
  Set \code{str0 = NULL} for none.


  }
  \item{summary.arg}{
  Logical. If \code{TRUE} then a summary is returned.
  If \code{TRUE} then \code{y} may be the output (fitted
  object) of \code{grc()}.


  }
  \item{h.step}{
  A small positive value that is passed into
  \code{summary.rrvglm()}. Only used when \code{summary.arg = TRUE}.


  }
  \item{\dots}{ Arguments that are passed
  into \code{rrvglm.control()}.


  }
  \item{M1}{
  The number of linear predictors of the \pkg{VGAM} \code{family} function
  for an ordinary (univariate) response.
  Then the number of linear predictors of the \code{rcim()} fit is
  usually the number of columns of \code{y} multiplied by \code{M1}.
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
  \item{has.intercept}{
  Logical. Include an intercept?


  }
  \item{M, cindex}{
  \eqn{M}  is the usual \pkg{VGAM} \eqn{M}, viz. the number of linear/additive
  predictors in total.
  Also, \code{cindex} means column index, and these point to the columns
  of \code{y} which are part of the vector of linear/additive predictor
  \emph{main effects}.


  For \code{family = multinomial} it is necessary to input these arguments
  as \code{M = ncol(y)-1} and
  \code{cindex = 2:(ncol(y)-1)}.


% except for the possibly the first one (due to identifiability constraints).


  }
  \item{rindex, iindex}{
  \code{rindex} means row index, and these are similar to \code{cindex}.
  \code{iindex} means interaction index, and these are similar to \code{cindex}.


  }
}
\details{
  Goodman's RC association model fits a reduced-rank approximation
  to a table of counts.
  A Poisson model is assumed.
  The log of each cell mean is decomposed as an
  intercept plus a row effect plus a column effect plus a reduced-rank
  component. The latter can be collectively written \code{A \%*\% t(C)},
  the product of two `thin' matrices.
  Indeed, \code{A} and \code{C} have \code{Rank} columns.
By default, the first column and row of the interaction matrix
\code{A \%*\% t(C)} is chosen
to be structural zeros, because \code{str0 = 1}.
This means the first row of \code{A} are all zeros.


This function uses \code{options()$contrasts} to set up the row and
column indicator variables.
In particular, Equation (4.5) of Yee and Hastie (2003) is used.
These are called \code{Row.} and \code{Col.} (by default) followed
by the row or column number.


The function \code{rcim()} is more general than \code{grc()}.
Its default is a no-interaction model of \code{grc()}, i.e.,
rank-0 and a Poisson distribution. This means that each
row and column has a dummy variable associated with it.
The first row and first column are baseline.
The power of \code{rcim()} is that many \pkg{VGAM} family functions
can be assigned to its \code{family} argument.
For example,
\code{\link{uninormal}} fits something in between a 2-way
ANOVA with and without interactions,
\code{\link{alaplace2}} with \code{Rank = 0} is something like
\code{\link[stats]{medpolish}}.
Others include
\code{\link{zipoissonff}} and
\code{\link{negbinomial}}.
Hopefully one day \emph{all} \pkg{VGAM} family functions will
work when assigned to the \code{family} argument, although the
result may not have meaning.



\emph{Unconstrained quadratic ordination} (UQO) can be performed
using \code{rcim()} and \code{grc()}.
This has been called \emph{unconstrained Gaussian ordination}
in the literature, however the word \emph{Gaussian} has two
meanings which is confusing; it is better to use
\emph{quadratic} because the bell-shape response surface is meant.
UQO is similar to CQO (\code{\link{cqo}}) except there are
no environmental/explanatory variables.
Here, a GLM is fitted to each column (species)
that is a quadratic function of hypothetical latent variables
or gradients.
Thus each row of the response has an associated site score,
and each column of the response has an associated optimum
and tolerance matrix.
  UQO can be performed here under the assumption that all species
  have the same tolerance matrices.
  See Yee and Hadi (2014) for details.
  It is not recommended that presence/absence data be inputted
  because the information content is so low for each site-species
  cell.
  The example below uses Poisson counts.



}
\value{
  An object of class \code{"grc"}, which currently is the same as
  an \code{"rrvglm"} object.
  Currently,
  a rank-0 \code{rcim()} object is of class \code{\link{rcim0-class}},
  else of class \code{"rcim"} (this may change in the future).


% Currently,
% a rank-0 \code{rcim()} object is of class \code{\link{vglm-class}},
% but it may become of class \code{"rcim"} one day.


}
\references{

Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


Yee, T. W. and Hadi, A. F. (2014)
Row-column interaction models, with an R implementation.
\emph{Computational Statistics},
\bold{29}, 1427--1445.


Goodman, L. A. (1981)
Association models and canonical correlation in the analysis
of cross-classifications having ordered categories.
\emph{Journal of the American Statistical Association},
\bold{76}, 320--334.


%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information about the setting up of the
%indicator variables.


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
  Setting \code{trace = TRUE} is recommended to monitor
  convergence.
  Using \code{criterion = "coefficients"} can result in slow convergence.


  If \code{summary = TRUE} then \code{y} can be a
  \code{"grc"} object, in which case a summary can be returned.
  That is, \code{grc(y, summary = TRUE)} is
  equivalent to \code{summary(grc(y))}.
  It is not possible to plot a
   \code{grc(y, summary = TRUE)} or
  \code{rcim(y, summary = TRUE)} object.


}

\section{Warning}{
  The function \code{rcim()} is experimental at this stage and
  may have bugs.
  Quite a lot of expertise is needed when fitting and in its
  interpretion thereof. For example, the constraint
  matrices applies the reduced-rank regression to the first
  (see \code{which.linpred})
  linear predictor and the other linear predictors are intercept-only and
  have a common value throughout the entire data set.
  This means that, by default, \code{family =} \code{\link{zipoissonff}} is
  appropriate but not
  \code{family =} \code{\link{zipoisson}}.
  Else set \code{family =} \code{\link{zipoisson}} and \code{which.linpred = 2}.
  To understand what is going on, do examine the constraint
  matrices of the fitted object, and reconcile this with
  Equations (4.3) to (4.5) of Yee and Hastie (2003).


  The functions temporarily create a permanent data frame
  called \code{.grc.df} or \code{.rcim.df}, which used
  to be needed by \code{summary.rrvglm()}. Then these
  data frames are deleted before exiting the function.
  If an error occurs then the data frames may be present
  in the workspace.



}

\seealso{
  \code{\link{rrvglm}},
  \code{\link{rrvglm.control}},
  \code{\link{rrvglm-class}},
  \code{summary.grc},
  \code{\link{moffset}},
  \code{\link{Rcim}},
  \code{\link{Select}},
  \code{\link{Qvar}},
  \code{\link{plotrcim0}},
  \code{\link{cqo}},
  \code{\link{multinomial}},
  \code{\link{alcoff}},
  \code{\link{crashi}},
  \code{\link{auuc}},
  \code{\link[VGAM:olym08]{olym08}},
  \code{\link[VGAM:olym12]{olym12}},
  \code{\link{poissonff}},
  \code{\link[stats]{medpolish}}.


}

\examples{
# Example 1: Undergraduate enrolments at Auckland University in 1990
fitted(grc1 <- grc(auuc))
summary(grc1)

grc2 <- grc(auuc, Rank = 2, Index.corner = c(2, 5))
fitted(grc2)
summary(grc2)

model3 <- rcim(auuc, Rank = 1, fam = multinomial,
               M = ncol(auuc)-1, cindex = 2:(ncol(auuc)-1), trace = TRUE)
fitted(model3)
summary(model3)

# Median polish but not 100 percent reliable. Maybe call alaplace2()...
\dontrun{
rcim0 <- rcim(auuc, fam = alaplace1(tau = 0.5), trace = FALSE, maxit = 500)
round(fitted(rcim0), digits = 0)
round(100 * (fitted(rcim0) - auuc) / auuc, digits = 0)  # Discrepancy
depvar(rcim0)
round(coef(rcim0, matrix = TRUE), digits = 2)
Coef(rcim0, matrix = TRUE)
# constraints(rcim0)
names(constraints(rcim0))

# Compare with medpolish():
(med.a <- medpolish(auuc))
fv <- med.a$overall + outer(med.a$row, med.a$col, "+")
round(100 * (fitted(rcim0) - fv) / fv)  # Hopefully should be all 0s
}


# Example 2: 2012 Summer Olympic Games in London
\dontrun{ top10 <- head(olym12, 10)
grc1.oly12 <- with(top10, grc(cbind(gold, silver, bronze)))
round(fitted(grc1.oly12))
round(resid(grc1.oly12, type = "response"), digits = 1)  # Response residuals
summary(grc1.oly12)
Coef(grc1.oly12)
}


# Example 3: Unconstrained quadratic ordination (UQO); see Yee and Hadi (2014)
\dontrun{
n <- 100; p <- 5; S <- 10
pdata <- rcqo(n, p, S, es.opt = FALSE, eq.max = FALSE,
              eq.tol = TRUE, sd.latvar = 0.75)  # Poisson counts
true.nu <- attr(pdata, "latvar")  # The 'truth'; site scores
attr(pdata, "tolerances")  # The 'truth'; tolerances

Y <- Select(pdata, "y", sort = FALSE)  # Y matrix (n x S); the "y" vars
uqo.rcim1 <- rcim(Y, Rank = 1,
                  str0 = NULL,  # Delta covers entire n x M matrix
                  iindex = 1:nrow(Y),  # RRR covers the entire Y
                  has.intercept = FALSE)  # Suppress the intercept

# Plot 1
par(mfrow = c(2, 2))
plot(attr(pdata, "optimums"), Coef(uqo.rcim1)@A,
     col = "blue", type = "p", main = "(a) UQO optimums",
     xlab = "True optimums", ylab = "Estimated (UQO) optimums")
mylm <- lm(Coef(uqo.rcim1)@A ~ attr(pdata, "optimums"))
abline(coef = coef(mylm), col = "orange", lty = "dashed")

# Plot 2
fill.val <- NULL  # Choose this for the new parameterization
plot(attr(pdata, "latvar"), c(fill.val, concoef(uqo.rcim1)),
     las = 1, col = "blue", type = "p", main = "(b) UQO site scores",
     xlab = "True site scores", ylab = "Estimated (UQO) site scores" )
mylm <- lm(c(fill.val, concoef(uqo.rcim1)) ~ attr(pdata, "latvar"))
abline(coef = coef(mylm), col = "orange", lty = "dashed")

# Plots 3 and 4
myform <- attr(pdata, "formula")
p1ut <- cqo(myform, family = poissonff,
            eq.tol = FALSE, trace = FALSE, data = pdata)
c1ut <- cqo(Select(pdata, "y", sort = FALSE) ~ scale(latvar(uqo.rcim1)),
            family = poissonff, eq.tol = FALSE, trace = FALSE, data = pdata)
lvplot(p1ut, lcol = 1:S, y = TRUE, pcol = 1:S, pch = 1:S, pcex = 0.5,
       main = "(c) CQO fitted to the original data",
       xlab = "Estimated (CQO) site scores")
lvplot(c1ut, lcol = 1:S, y = TRUE, pcol = 1:S, pch = 1:S, pcex = 0.5,
       main = "(d) CQO fitted to the scaled UQO site scores",
       xlab = "Estimated (UQO) site scores")
}
}
\keyword{models}
\keyword{regression}
% plot(grc.oly1)
% oly2 <- with(top10, grc(cbind(gold,silver,bronze), Rank = 2))  # Saturated model
% round(fitted(oly2))
% round(fitted(oly2)) - with(top10, cbind(gold,silver,bronze))
% summary(oly2)  # Saturated model



% zz 20100927 unsure
% Then \code{.grc.df} is deleted before exiting the function.


% print(Coef(rcim0, matrix = TRUE), digits = 3)


% Prior to 201310:
%    str0 = if (!Rank) NULL else {
%      if (M1 == 1) 1 else setdiff(1:(M1 * ncol(y)),
%                   c(1 + (1:ncol(y)) * M1, Index.corner))
%    },
%    str0 = if (Rank > 0) 1 else NULL,
%    Index.corner = if (!Rank) NULL else 1 + M1 * (1:Rank),





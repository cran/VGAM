\name{qrrvglm.control}
\alias{qrrvglm.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Control Function for QRR-VGLMs (CQO) }
\description{
  Algorithmic constants and parameters for a constrained quadratic
  ordination (CQO), by fitting a \emph{quadratic reduced-rank vector
  generalized linear model} (QRR-VGLM), are set using this function.
  It is the control function for \code{\link{cqo}}.

}
\usage{
qrrvglm.control(Rank = 1,
                Bestof = if (length(Cinit)) 1 else 10,
                checkwz = TRUE,
                Cinit = NULL,
                Crow1positive = TRUE,
                epsilon = 1.0e-06,
                EqualTolerances = NULL,
                eq.tolerances = TRUE,
                Etamat.colmax = 10,
                FastAlgorithm = TRUE,
                GradientFunction = TRUE,
                Hstep = 0.001,
                isd.latvar = rep_len(c(2, 1, rep_len(0.5, Rank)), Rank),
                iKvector = 0.1,
                iShape = 0.1,
                ITolerances = NULL,
                I.tolerances = FALSE,
                maxitl = 40,
                imethod = 1,
                Maxit.optim = 250,
                MUXfactor = rep_len(7, Rank),
                noRRR = ~ 1, Norrr = NA,
                optim.maxit = 20,
                Parscale = if (I.tolerances) 0.001 else 1.0,
                sd.Cinit = 0.02,
                SmallNo = 5.0e-13,
                trace = TRUE,
                Use.Init.Poisson.QO = TRUE,
                wzepsilon = .Machine$double.eps^0.75, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
    In the following, \eqn{R} is the \code{Rank}, \eqn{M} is the number
    of linear predictors, and \eqn{S} is the number of responses
    (species).
    Thus \eqn{M=S} for binomial and Poisson responses, and
    \eqn{M=2S} for the negative binomial and 2-parameter gamma distributions.

  \item{Rank}{
    The numerical rank \eqn{R} of the model, i.e., the
    number of ordination axes. Must be an element from the set
    \{1,2,\ldots,min(\eqn{M},\eqn{p_2}{p2})\} where the vector of explanatory
    variables \eqn{x} is partitioned into (\eqn{x_1},\eqn{x_2}), which is
    of dimension \eqn{p_1+p_2}{p1+p2}. The variables making up \eqn{x_1}
    are given by the terms in the \code{noRRR} argument, and the rest
    of the terms comprise \eqn{x_2}.

    }
  \item{Bestof}{
    Integer. The best of \code{Bestof} models fitted is returned.
    This argument helps guard against local solutions by (hopefully)
    finding the global solution from many fits. The argument has value
    1 if an initial value for \eqn{C} is inputted using \code{Cinit}.

    }
  \item{checkwz}{ logical indicating whether the diagonal elements of
    the working weight matrices should be checked whether they are
    sufficiently positive, i.e., greater than \code{wzepsilon}. If not,
    any values less than \code{wzepsilon} are replaced with this value.

  }
  \item{Cinit}{
    Optional initial \eqn{C} matrix, which must be a \eqn{p_2}{p2} by \eqn{R}
    matrix. The default is to apply \code{.Init.Poisson.QO()} to obtain
    initial values.

  }
  \item{Crow1positive}{
      Logical vector of length \code{Rank} (recycled if necessary): are
      the elements of the first row of \eqn{C} positive? For example,
      if \code{Rank} is 4, then specifying \code{Crow1positive = c(FALSE,
      TRUE)} will force \eqn{C[1,1]} and \eqn{C[1,3]} to be negative,
      and \eqn{C[1,2]} and \eqn{C[1,4]} to be positive. This argument
      allows for a reflection in the ordination axes because the
      coefficients of the latent variables are unique up to a sign.

  }
  \item{epsilon}{
    Positive numeric. Used to test for convergence for GLMs fitted in C.
    Larger values mean a loosening of the convergence criterion.
    If an error code of 3 is reported, try increasing this value.

  }
  \item{eq.tolerances}{
    Logical indicating whether each (quadratic) predictor will
    have equal tolerances. Having \code{eq.tolerances = TRUE}
    can help avoid numerical problems, especially with binary data.
    Note that the estimated (common) tolerance matrix may or may
    not be positive-definite. If it is  then it can be scaled to
    the \eqn{R} by \eqn{R} identity matrix, i.e., made equivalent
    to \code{I.tolerances = TRUE}. Setting \code{I.tolerances = TRUE}
    will \emph{force} a common \eqn{R} by \eqn{R} identity matrix as
    the tolerance matrix to the data even if it is not appropriate.
    In general, setting \code{I.tolerances = TRUE} is
    preferred over \code{eq.tolerances = TRUE} because,
    if it works, it is much faster and uses less memory.
    However, \code{I.tolerances = TRUE} requires the
    environmental variables to be scaled appropriately.
    See \bold{Details} for more details.

  }
  \item{EqualTolerances}{
    Defunct argument.
    Use \code{eq.tolerances} instead.

  }











% \item{Eta.range}{ Numerical vector of length 2 or \code{NULL}.
%   Gives the lower and upper bounds on the values that can be taken
%   by the quadratic predictor (i.e., on the eta-scale).
%   Since \code{FastAlgorithm = TRUE}, this argument should be ignored.
% }
  \item{Etamat.colmax}{
    Positive integer, no smaller than \code{Rank}. Controls the amount
    of memory used by \code{.Init.Poisson.QO()}. It is the maximum
    number of columns allowed for the pseudo-response and its weights.
    In general, the larger the value, the better the initial value.
    Used only if \code{Use.Init.Poisson.QO = TRUE}.

  }

  \item{FastAlgorithm}{
   Logical. Whether a new fast algorithm is to be used. The fast
   algorithm results in a large speed increases compared to Yee
   (2004).
   Some details of the fast algorithm are found in Appendix A of Yee (2006).
   Setting \code{FastAlgorithm = FALSE} will give an error.

  }
  \item{GradientFunction}{
   Logical. Whether \code{\link[stats]{optim}}'s argument \code{gr}
   is used or not, i.e., to compute gradient values. Used only if
   \code{FastAlgorithm} is \code{TRUE}. The default value is usually
   faster on most problems.

  }
  \item{Hstep}{
   Positive value. Used as the step size in the finite difference
   approximation to the derivatives by \code{\link[stats]{optim}}.

%  Used only if \code{FastAlgorithm} is \code{TRUE}.
  }
  \item{isd.latvar}{
   Initial standard deviations for the latent variables (site scores).
   Numeric, positive and of length \eqn{R} (recycled if necessary).
   This argument is used only if \code{I.tolerances = TRUE}. Used by
   \code{.Init.Poisson.QO()} to obtain initial values for the constrained
   coefficients \eqn{C} adjusted to a reasonable value. It adjusts the
   spread of the site scores relative to a common species tolerance of 1
   for each ordination axis. A value between 0.5 and 10 is recommended;
   a value such as 10 means that the range of the environmental space is
   very large relative to the niche width of the species. The successive
   values should decrease because the first ordination axis should have
   the most spread of site scores, followed by the second ordination
   axis, etc.

  }
  \item{iKvector, iShape}{
   Numeric, recycled to length \eqn{S} if necessary.
   Initial values used for estimating the positive \eqn{k} and
   \eqn{\lambda}{lambda} parameters of the negative binomial and
   2-parameter gamma distributions respectively. For further information
   see \code{\link{negbinomial}} and \code{\link{gamma2}}.
   These arguments override the \code{ik} and \code{ishape}
   arguments in \code{\link{negbinomial}} and \code{\link{gamma2}}.

  }

  \item{I.tolerances}{
   Logical. If \code{TRUE} then the (common) tolerance matrix is the
   \eqn{R} by \eqn{R} identity matrix by definition. Note that having
   \code{I.tolerances = TRUE} implies \code{eq.tolerances = TRUE}, but
   not vice versa. Internally, the quadratic terms will be treated as
   offsets (in GLM jargon) and so the models can potentially be fitted
   very efficiently. \emph{However, it is a very good idea to center
   and scale all numerical variables in the \eqn{x_2} vector}.
   See \bold{Details} for more details.
   The success of \code{I.tolerances = TRUE} often
   depends on suitable values for \code{isd.latvar} and/or
   \code{MUXfactor}.


  }
  \item{ITolerances}{
    Defunct argument.
    Use \code{I.tolerances} instead.


 }




 \item{maxitl}{
    Maximum number of times the optimizer is called or restarted.
    Most users should ignore this argument.

  }
  \item{imethod}{
    Method of initialization. A positive integer 1 or 2 or 3 etc.
    depending on the \pkg{VGAM} family function.
    Currently it is used for \code{\link{negbinomial}} and
    \code{\link{gamma2}} only, and used within the C.

  }
  \item{Maxit.optim}{
    Positive integer. Number of iterations given to the function
    \code{\link[stats]{optim}} at each of the \code{optim.maxit}
    iterations.

  }

  \item{MUXfactor}{
   Multiplication factor for detecting large offset values. Numeric,
   positive and of length \eqn{R} (recycled if necessary). This argument
   is used only if \code{I.tolerances = TRUE}. Offsets are \eqn{-0.5}
   multiplied by the sum of the squares of all \eqn{R} latent variable
   values. If the latent variable values are too large then this will
   result in numerical problems. By too large, it is meant that the
   standard deviation of the latent variable values are greater than
   \code{MUXfactor[r] * isd.latvar[r]} for \code{r=1:Rank} (this is why
   centering and scaling all the numerical predictor variables in
   \eqn{x_2} is recommended). A value about 3 or 4 is recommended.
   If failure to converge occurs, try a slightly lower value.

}
  \item{optim.maxit}{
    Positive integer. Number of times \code{\link[stats]{optim}}
    is invoked. At iteration \code{i}, the \code{i}th value of
    \code{Maxit.optim} is fed into \code{\link[stats]{optim}}.

  }

  \item{noRRR}{
    Formula giving terms that are \emph{not} to be included in the
    reduced-rank regression (or formation of the latent variables),
    i.e., those belong to \eqn{x_1}.
    Those variables which do not make up the latent variable (reduced-rank
    regression) correspond to the \eqn{B_1}{B_1} matrix.
    The default is to omit the intercept term from the latent variables.


  }
  \item{Norrr}{
  Defunct. Please use \code{noRRR}.
  Use of \code{Norrr} will become an error soon.


  }


  \item{Parscale}{
    Numerical and positive-valued vector of length \eqn{C}
   (recycled if necessary).
   Passed into \code{optim(..., control = list(parscale = Parscale))};
   the elements of \eqn{C} become \eqn{C} / \code{Parscale}.
   Setting \code{I.tolerances = TRUE} results in line searches that
   are very large, therefore \eqn{C} has to be scaled accordingly
   to avoid large step sizes.
   See \bold{Details} for more information.
   It's probably best to leave this argument alone.
  }
  \item{sd.Cinit}{
      Standard deviation of the initial values for the elements
      of \eqn{C}.
      These are normally distributed with mean zero.
      This argument is used only if \code{Use.Init.Poisson.QO = FALSE}
      and \eqn{C} is not inputted using \code{Cinit}.

  }
  \item{trace}{
    Logical indicating if output should be produced for
    each iteration. The default is \code{TRUE} because the
    calculations are numerically intensive, meaning it may take
    a long time, so that the user might think the computer has
    locked up if \code{trace = FALSE}.

}

%  \item{Kinit}{ Initial values for the index parameters \code{k} in the
%   negative binomial distribution (one per species).
%   In general, a smaller number is preferred over a larger number.
%   The vector is recycled to the number of responses (species).
%   The argument is currently not used.
%   }

%  \item{Dzero}{ Integer vector specifying which squared terms
%      are to be zeroed. These linear predictors will correspond to
%      a RR-VGLM.
%      The values must be elements from the set \{1,2,\ldots,\eqn{M}\}.
%      Used only if \code{Quadratic = TRUE} and \code{FastAlgorithm = FALSE}.
%  }
\item{SmallNo}{
  Positive numeric between \code{.Machine$double.eps} and \code{0.0001}.
  Used to avoid under- or over-flow in the IRLS algorithm.
  Used only if \code{FastAlgorithm} is \code{TRUE}.
  }
  \item{Use.Init.Poisson.QO}{
    Logical. If \code{TRUE} then the function \code{.Init.Poisson.QO()} is
    used to obtain initial values for the canonical coefficients \eqn{C}.
    If \code{FALSE} then random numbers are used instead.
  }
  \item{wzepsilon}{
  Small positive number used to test whether the diagonals of the working
  weight matrices are sufficiently positive.

  }
\item{\dots}{ Ignored at present. }
}
\details{
  Recall that the central formula for CQO is
  \deqn{\eta = B_1^T x_1 + A \nu +
               \sum_{m=1}^M (\nu^T D_m \nu) e_m}{%
         eta = B_1^T x_1 + A nu +
         sum_{m=1}^M (nu^T D_m nu) e_m}
  where \eqn{x_1}{x_1} is a vector (usually just a 1 for an intercept),
  \eqn{x_2}{x_2} is a vector of environmental variables, \eqn{\nu=C^T
  x_2}{nu=C^T x_2} is a \eqn{R}-vector of latent variables, \eqn{e_m} is
  a vector of 0s but with a 1 in the \eqn{m}th position.
    QRR-VGLMs are an extension of RR-VGLMs and allow for maximum
    likelihood solutions to constrained quadratic ordination (CQO) models.


%    For the fitting of QRR-VGLMs, the default is that the \eqn{C} matrix
%    (containing the \emph{canonical} or \emph{constrained coefficients}
%    corresponding to \eqn{x_2})
%    is constrained by forcing the latent variables to have sample
%    variance-covariance matrix equalling \code{diag(Rank)}, i.e.,
%    unit variance and uncorrelated. The tolerance matrices are, in
%    general, diagonal under such a constraint.


   Having \code{I.tolerances = TRUE} means all the tolerance matrices
   are the order-\eqn{R} identity matrix, i.e., it \emph{forces}
   bell-shaped curves/surfaces on all species. This results in a
   more difficult optimization problem (especially for 2-parameter
   models such as the negative binomial and gamma) because of overflow
   errors and it appears there are more local solutions. To help avoid
   the overflow errors, scaling \eqn{C} by the factor \code{Parscale}
   can help enormously. Even better, scaling \eqn{C} by specifying
   \code{isd.latvar} is more understandable to humans. If failure to
   converge occurs, try adjusting \code{Parscale}, or better, setting
   \code{eq.tolerances = TRUE} (and hope that the estimated tolerance
   matrix is positive-definite). To fit an equal-tolerances model, it
   is firstly best to try setting \code{I.tolerances = TRUE} and varying
   \code{isd.latvar} and/or \code{MUXfactor} if it fails to converge.
   If it still fails to converge after many attempts, try setting
   \code{eq.tolerances = TRUE}, however this will usually be a lot slower
   because it requires a lot more memory.


   With a \eqn{R > 1} model, the latent variables are always uncorrelated,
   i.e., the variance-covariance matrix of the site scores is a diagonal
   matrix.


   If setting \code{eq.tolerances = TRUE} is used and the common
   estimated tolerance matrix is positive-definite then that model is
   effectively the same as the \code{I.tolerances = TRUE} model (the two are
   transformations of each other). In general, \code{I.tolerances = TRUE}
   is numerically more unstable and presents a more difficult problem
   to optimize; the arguments \code{isd.latvar} and/or \code{MUXfactor} often
   must be assigned some good value(s) (possibly found by trial and error)
   in order for convergence to occur. Setting \code{I.tolerances = TRUE}
   \emph{forces} a bell-shaped curve or surface onto all the species data,
   therefore this option should be used with deliberation. If unsuitable,
   the resulting fit may be very misleading. Usually it is a good idea
   for the user to set \code{eq.tolerances = FALSE} to see which species
   appear to have a bell-shaped curve or surface. Improvements to the
   fit can often be achieved using transformations, e.g., nitrogen
   concentration to log nitrogen concentration.


   Fitting a CAO model (see \code{\link{cao}}) first is a good idea for
   pre-examining the data and checking whether it is appropriate to fit
   a CQO model.


%Suppose \code{FastAlgorithm = FALSE}. In theory (if
%\code{Eta.range = NULL}), for QRR-VGLMs, the predictors have the values of
%a quadratic form. However, when \code{Eta.range} is assigned a numerical
%vector of length 2 (giving the endpoints of an interval), then those
%values lying outside the interval are assigned the closest boundary
%value.  The \code{Eta.range} argument is provided to help avoid
%numerical problems associated with the inner minimization problem.  A
%consequence of this is that the fitted values are bounded, e.g., between
%\code{1/(1+exp(-Eta.range[1]))} and \code{1/(1+exp(-Eta.range[2]))} for
%binary data (logit link), and greater than \code{exp(Eta.range[1])} for
%Poisson data (log link).  It is suggested that, for binary responses,
%\code{c(-16, 16)} be used, and for Poisson responses, \code{c(-16, Inf)}
%be used.  The value \code{NULL} corresponds to \code{c(-Inf, Inf)}.

}
\value{
  A list with components matching the input names.


}
\references{
Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology},
\bold{87}, 203--213.


}
\author{ Thomas W. Yee }
\note{
  When \code{I.tolerances = TRUE} it is a good idea to apply
  \code{\link[base]{scale}} to all the numerical variables that make up
  the latent variable, i.e., those of \eqn{x_2}{x_2}. This is to make
  them have mean 0, and hence avoid large offset values which cause
  numerical problems.


  This function has many arguments that are common with
  \code{\link{rrvglm.control}} and \code{\link{vglm.control}}.


  It is usually a good idea to try fitting a model with
  \code{I.tolerances = TRUE} first, and if convergence is unsuccessful,
  then try \code{eq.tolerances = TRUE} and \code{I.tolerances = FALSE}.
  Ordination diagrams with
  \code{eq.tolerances = TRUE} have a natural interpretation, but
  with \code{eq.tolerances = FALSE} they are more complicated and
  requires, e.g., contours to be overlaid on the ordination diagram
  (see \code{\link{lvplot.qrrvglm}}).


% and/or use the \code{Eta.range} argument.


  In the example below, an equal-tolerances CQO model is fitted to the
  hunting spiders data. Because \code{I.tolerances = TRUE}, it is a good idea
  to center all the \eqn{x_2} variables first. Upon fitting the model,
  the actual standard deviation of the site scores are computed. Ideally,
  the \code{isd.latvar} argument should have had this value for the best
  chances of getting good initial values.  For comparison, the model is
  refitted with that value and it should run more faster and reliably.
}

\section{Warning }{

  The default value of \code{Bestof} is a bare minimum for many datasets,
  therefore it will be necessary to increase its value to increase the
  chances of obtaining the global solution.


%Suppose \code{FastAlgorithm = FALSE}.
%The fitted values of QRR-VGLMs can be restricted to lie between two values
%in order to help make the computation numerically stable.  For some data
%sets, it may be necessary to use \code{Eta.range} to obtain convergence;
%however, the fitted values etc. will no longer be accurate, especially at
%small and/or large values.  Convergence is slower when \code{Eta.range}
%is used to restrict values.

}

\seealso{
  \code{\link{cqo}},
  \code{\link{rcqo}},
  \code{\link{Coef.qrrvglm}},
  \code{\link{Coef.qrrvglm-class}},
  \code{\link[stats]{optim}},
  \code{\link{binomialff}},
  \code{\link{poissonff}},
  \code{\link{negbinomial}},
  \code{\link{gamma2}},
  \code{\link{gaussianff}}.


% \code{\link{rrvglm}},
% \code{\link{rrvglm.control}},
% \code{\link{rrvglm.optim.control}},


}

\examples{
\dontrun{ # Poisson CQO with equal tolerances
set.seed(111)  # This leads to the global solution
hspider[,1:6] <- scale(hspider[,1:6])  # Good idea when I.tolerances = TRUE
p1 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi, Auloalbi,
               Pardlugu, Pardmont, Pardnigr, Pardpull, Trocterr, Zoraspin) ~
          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
          quasipoissonff, data = hspider, eq.tolerances = TRUE)
sort(deviance(p1, history = TRUE))  # A history of all the iterations

(isd.latvar <- apply(latvar(p1), 2, sd))  # Should be approx isd.latvar

# Refit the model with better initial values
set.seed(111)  # This leads to the global solution
p1 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi, Auloalbi,
               Pardlugu, Pardmont, Pardnigr, Pardpull, Trocterr, Zoraspin) ~
          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
          I.tolerances = TRUE, quasipoissonff, data = hspider,
          isd.latvar = isd.latvar)  # Note the use of isd.latvar here
sort(deviance(p1, history = TRUE))  # A history of all the iterations
}
}
\keyword{models}
\keyword{regression}



%\dontrun{
%# 20120221; withdrawn for a while coz it creates a lot of error messages.
%# Negative binomial CQO; smallest deviance is about 275.389
%set.seed(1234)  # This leads to a reasonable (but not the global) solution?
%nb1 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi, Auloalbi,
%                Pardlugu, Pardmont, Pardnigr, Pardpull, Trocterr, Zoraspin) ~
%          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
%          I.tol = FALSE, eq.tol = TRUE,  # A good idea for negbinomial
%          fam = negbinomial, data = hspider)
%sort(deviance(nb1, history = TRUE))  # A history of all the iterations
%summary(nb1)
%}
%\dontrun{ lvplot(nb1, lcol = 1:12, y = TRUE, pcol = 1:12) }




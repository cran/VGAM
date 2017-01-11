\name{cao.control}
\alias{cao.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Control Function for RR-VGAMs (CAO) }
\description{
  Algorithmic constants and parameters for a constrained additive
  ordination (CAO), by fitting a \emph{reduced-rank vector generalized
  additive model} (RR-VGAM), are set using this function.
  This is the control function for \code{\link{cao}}.

}
\usage{
cao.control(Rank = 1, all.knots = FALSE, criterion = "deviance", Cinit = NULL,
            Crow1positive = TRUE, epsilon = 1.0e-05, Etamat.colmax = 10,
            GradientFunction = FALSE, iKvector = 0.1, iShape = 0.1,
            noRRR = ~ 1, Norrr = NA,
            SmallNo = 5.0e-13, Use.Init.Poisson.QO = TRUE,
            Bestof = if (length(Cinit)) 1 else 10, maxitl = 10,
            imethod = 1, bf.epsilon = 1.0e-7, bf.maxit = 10,
            Maxit.optim = 250, optim.maxit = 20, sd.sitescores = 1.0,
            sd.Cinit = 0.02, suppress.warnings = TRUE,
            trace = TRUE, df1.nl = 2.5, df2.nl = 2.5,
            spar1 = 0, spar2 = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{Rank}{
    The numerical rank \eqn{R} of the model, i.e., the number of latent
    variables.  Currently only \code{Rank = 1} is implemented.

  }
  \item{all.knots}{
    Logical indicating if all distinct points of the smoothing
    variables are to be used as knots.  Assigning the value
    \code{FALSE} means fewer knots are chosen when the number
    of distinct points is large, meaning less computational
    expense. See \code{\link{vgam.control}} for details.

  }
  \item{criterion}{
    Convergence criterion. Currently, only one is supported:
    the deviance is minimized.

  }
  \item{Cinit}{
    Optional initial \bold{C} matrix which may speed up
    convergence.

  }
  \item{Crow1positive}{
    Logical vector of length \code{Rank} (recycled if
    necessary): are the elements of the first row of \bold{C}
    positive?  For example, if \code{Rank} is 4, then specifying
    \code{Crow1positive = c(FALSE, TRUE)} will force \bold{C[1,1]}
    and \bold{C[1,3]} to be negative, and \bold{C[1,2]} and
    \bold{C[1,4]} to be positive.

  }
    \item{epsilon}{
    Positive numeric. Used to test for convergence for GLMs
    fitted in FORTRAN.  Larger values mean a loosening of the
    convergence criterion.

%   Used only if \code{FastAlgorithm} is \code{TRUE}.

    }
  \item{Etamat.colmax}{
    Positive integer, no smaller than \code{Rank}.  Controls the
    amount of memory used by \code{.Init.Poisson.QO()}.  It is the
    maximum number of columns allowed for the pseudo-response and
    its weights.  In general, the larger the value, the better the
    initial value.  Used only if \code{Use.Init.Poisson.QO = TRUE}.

  }

% \item{FastAlgorithm}{
%   Logical.
%   Whether compiled code is used.
%   For \code{\link{cao}} this must be \code{TRUE}.
%
% }

\item{GradientFunction}{
  Logical. Whether \code{\link[stats]{optim}}'s argument \code{gr}
  is used or not, i.e., to compute gradient values.  Used only if
  \code{FastAlgorithm} is \code{TRUE}.  Currently, this argument
  must be set to \code{FALSE}.



}
  \item{iKvector, iShape}{
  See \code{\link{qrrvglm.control}}.

}

%  \item{Hstep}{ Positive value. Used as the step size in the
%   finite difference approximation to the derivatives by
%   \code{\link[stats]{optim}}.
%   Used only if \code{GradientFunction} is \code{TRUE}.
%
%  }

% \item{Kinit}{
%  Initial values for the index parameters \code{k} in the negative
%  binomial distribution (one per species).  In general, a smaller number
%  is preferred over a larger number.  The vector is recycled to the
%  number of responses (species).

%}

  \item{noRRR}{
  Formula giving terms that are \emph{not} to be included in the
  reduced-rank regression (or formation of the latent variables).
  The default is to omit the intercept term from the latent
  variables.  Currently, only \code{noRRR = ~ 1} is implemented.


  }
  \item{Norrr}{
  Defunct. Please use \code{noRRR}.
  Use of \code{Norrr} will become an error soon.


  }



% \item{Parscale}{
%  Numerical and positive-valued vector of length \bold{C}
%  (recycled if necessary).  Passed into \code{optim(...,
%  control = list(parscale = Parscale))}; the elements of \bold{C} become
%  \bold{C} / \code{Parscale}.  Setting \code{I.tolerances = TRUE} results
%  in line searches that are very large, therefore \bold{C} has to be
%  scaled accordingly to avoid large step sizes.

% }

  \item{SmallNo}{
   Positive numeric between \code{.Machine$double.eps} and
   \code{0.0001}.  Used to avoid under- or over-flow in the
   IRLS algorithm.


%  Used only if \code{FastAlgorithm} is \code{TRUE}.


  }
  \item{Use.Init.Poisson.QO }{
  Logical. If \code{TRUE} then the function
  \code{.Init.Poisson.QO} is used to obtain initial values
  for the canonical coefficients \bold{C}.  If \code{FALSE}
  then random numbers are used instead.


  }

  \item{Bestof}{
    Integer. The best of \code{Bestof} models fitted is returned. This
    argument helps guard against local solutions by (hopefully) finding
    the global solution from many fits. The argument works only when
    the function generates its own initial value for \bold{C}, i.e.,
    when \bold{C} are \emph{not} passed in as initial values.
    The default is only a convenient minimal number and users are urged
    to increase this value.


  }
  \item{maxitl}{
    Positive integer. Maximum number of
    Newton-Raphson/Fisher-scoring/local-scoring iterations allowed.


  }
  \item{imethod}{
  See \code{\link{qrrvglm.control}}.

  }
  \item{bf.epsilon}{
    Positive numeric. Tolerance used by the modified vector backfitting
    algorithm for testing convergence.

  }
  \item{bf.maxit}{
    Positive integer.
    Number of backfitting iterations allowed in the compiled code.
  }
  \item{Maxit.optim}{
    Positive integer.
    Number of iterations given to the function \code{\link[stats]{optim}}
    at each of the \code{optim.maxit} iterations.

  }
  \item{optim.maxit}{
    Positive integer.
    Number of times \code{\link[stats]{optim}} is invoked.

%   At iteration \code{i}, the \code{i}th value of \code{Maxit.optim}
%   is fed into \code{\link[stats]{optim}}.

    }
% \item{se.fit}{
%   Logical indicating whether approximate
%   pointwise standard errors are to be saved on the object.
%   Currently this argument must have the value \code{FALSE}.

% }
  \item{sd.sitescores}{
    Numeric. Standard deviation of the
    initial values of the site scores, which are generated from
    a normal distribution.
    Used when \code{Use.Init.Poisson.QO} is \code{FALSE}.


    }
  \item{sd.Cinit}{
    Standard deviation of the initial values for the elements
    of \bold{C}.
    These are normally distributed with mean zero.
    This argument is used only if \code{Use.Init.Poisson.QO = FALSE}.

  }
  \item{suppress.warnings}{
    Logical. Suppress warnings?


  }
  \item{trace}{
    Logical indicating if output should be produced for each
    iteration. Having the value \code{TRUE} is a good idea for large
    data sets.

  }
  \item{df1.nl, df2.nl}{
    Numeric and non-negative, recycled to length \emph{S}.
    Nonlinear degrees
    of freedom for smooths of the first and second latent variables.
    A value of 0 means the smooth is linear.  Roughly, a value between
    1.0 and 2.0 often has the approximate flexibility of a quadratic.
    The user should not assign too large a value to this argument, e.g.,
    the value 4.0 is probably too high.  The argument \code{df1.nl} is
    ignored if \code{spar1} is assigned a positive value or values. Ditto
    for \code{df2.nl}.

  }
  \item{spar1, spar2}{
    Numeric and non-negative, recycled to length \emph{S}.
    Smoothing parameters of the
    smooths of the first and second latent variables. The larger the value, the
    more smooth (less wiggly) the fitted curves. These arguments are an
    alternative to specifying \code{df1.nl} and \code{df2.nl}.  A value
    0 (the default) for \code{spar1} means that \code{df1.nl} is used.
    Ditto for \code{spar2}.
    The values are on a scaled version of the latent variables.
    See Green and Silverman (1994) for more information.

  }
  \item{\dots}{ Ignored at present. }
}
\details{
  Many of these arguments are identical to \code{\link{qrrvglm.control}}.
  Here, \eqn{R} is the \code{Rank}, \eqn{M} is the number
  of additive predictors, and \eqn{S} is the number of responses
  (species).
  Thus \eqn{M=S} for binomial and Poisson responses, and
  \eqn{M=2S} for the negative binomial and 2-parameter gamma distributions.


  Allowing the smooths too much flexibility means the CAO optimization
  problem becomes more difficult to solve. This is because the number
  of local solutions increases as the nonlinearity of the smooths
  increases. In situations of high nonlinearity, many initial values
  should be used, so that \code{Bestof} should be assigned a larger
  value. In general, there should be a reasonable value of \code{df1.nl}
  somewhere between 0 and about 3 for most data sets.

}
\value{
  A list with the components corresponding to its arguments, after some
  basic error checking.

}
\references{


Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.


Green, P. J. and Silverman, B. W. (1994)
\emph{Nonparametric Regression and Generalized Linear Models:
      A Roughness Penalty Approach},
London: Chapman & Hall.


}
\author{T. W. Yee}
\note{
  The argument \code{df1.nl} can be inputted in the format \code{c(spp1 = 2,
  spp2 = 3, 2.5)}, say, meaning the default value is 2.5, but two species
  have alternative values.


  If \code{spar1 = 0} and \code{df1.nl = 0} then this represents fitting
  linear functions (CLO). Currently, this is handled in the awkward
  manner of setting \code{df1.nl} to be a small positive value, so that
  the smooth is almost linear but not quite.
  A proper fix to this special case should done in the short future.


}
\seealso{
  \code{\link{cao}}.


}

\examples{\dontrun{
hspider[,1:6] <- scale(hspider[,1:6])  # Standardized environmental vars
set.seed(123)
ap1 <- cao(cbind(Pardlugu, Pardmont, Pardnigr, Pardpull, Zoraspin) ~
           WaterCon + BareSand + FallTwig +
           CoveMoss + CoveHerb + ReflLux,
           family = poissonff, data = hspider,
           df1.nl = c(Zoraspin = 2.3, 2.1),
           Bestof = 10, Crow1positive = FALSE)
sort(deviance(ap1, history = TRUE))  # A history of all the iterations

Coef(ap1)

par(mfrow = c(2, 3))  # All or most of the curves are unimodal; some are
plot(ap1, lcol = "blue")  # quite symmetric. Hence a CQO model should be ok

par(mfrow = c(1, 1), las = 1)
index <- 1:ncol(depvar(ap1))  # lvplot is jagged because only 28 sites
lvplot(ap1, lcol = index, pcol = index, y = TRUE)

trplot(ap1, label = TRUE, col = index)
abline(a = 0, b = 1, lty = 2)

persp(ap1, label = TRUE, col = 1:4)
}
}
\keyword{models}
\keyword{regression}





%cao.control(Rank = 1, all.knots = FALSE,
%            criterion = "deviance", Cinit = NULL,
%            Crow1positive = TRUE, epsilon = 1e-05,
%            Etamat.colmax = 10,
%%           FastAlgorithm = TRUE,
%%                           is.loaded(symbol.For("cqo2f")),
%%           GradientFunction = FALSE,
%            iKvector = 0.1,
%            iShape = 0.1,
%            noRRR = ~1,
%%           Parscale = 1,
%            SmallNo = 5e-13,
%            Use.Init.Poisson.QO = TRUE,
%            Bestof = if (length(Cinit)) 1 else 10, maxitl = 40,
%            bf.epsilon = 1.0e-7, bf.maxit = 40,
%            Maxit.optim = 250, optim.maxit = 20,
%%           se.fit = FALSE,
%            sd.sitescores = 1,
%            sd.Cinit = 0.02, trace = TRUE,
%%            df1.nl = 2.5, spar1 = 0, ...)

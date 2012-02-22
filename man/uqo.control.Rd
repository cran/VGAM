\name{uqo.control}
\alias{uqo.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Control Function for UQO models }
\description{
  Algorithmic constants and parameters for an
  unconstrained quadratic ordination (UQO) model, by fitting a
  \emph{quadratic unconstrained vector generalized additive model}
  (QU-VGLM), are set using this function.
  It is the control function of \code{\link{uqo}}.
  
}
\usage{
uqo.control(Rank=1, Bestof = if (length(lvstart) &&
            !jitter.sitescores) 1 else 10, CA1 = FALSE, Crow1positive
            = TRUE, epsilon = 1.0e-07, EqualTolerances = ITolerances,
            Etamat.colmax = 10, GradientFunction=TRUE, Hstep = 0.001,
            isdlv = rep(c(2, 1, rep(0.5, len=Rank)), len=Rank),
            ITolerances = FALSE, lvstart = NULL, jitter.sitescores
            = FALSE, maxitl = 40, Maxit.optim = 250, MUXfactor =
            rep(3, length=Rank), optim.maxit = 20, nRmax = 250,
            SD.sitescores = 1.0, SmallNo = 5.0e-13, trace = TRUE,
            Use.Init.Poisson.QO=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Rank}{ The numerical rank \eqn{R} of the model,
    i.e., the number of latent variables or ordination axes.
    Currently only \eqn{R=1} is recommended.
  }
  \item{Bestof}{ Integer. The best of \code{Bestof} models fitted is
    returned. This argument helps guard against local solutions by
    (hopefully) finding the global solution from many fits.
    The argument has value 1 if an initial value for the site scores is
    inputted using \code{lvstart}.

  }
  \item{CA1}{ 
    Logical. If \code{TRUE} the site scores from a correspondence analysis
    (CA) are computed and used on the first axis as initial values.
    Both \code{CA1} and \code{Use.Init.Poisson.QO} cannot both be
    \code{TRUE}.

  }
  \item{Crow1positive}{ 
    Logical vector of length \code{Rank} (recycled if necessary):
    are the elements of the first row of the latent variable matrix
    \eqn{\nu}{nu} positive?
      For example, if \code{Rank} is 2, then specifying
      \code{Crow1positive=c(FALSE, TRUE)} will force the first
      site score's first element to be negative, and the first site
      score's second element to be positive.  Note that there is no
      \eqn{C} matrix with UQO, but the argument's name comes from
      \code{\link{qrrvglm.control}} and is left unchanged for convenience.

  }
    \item{epsilon}{
      Positive numeric. Used to test for convergence for GLMs fitted
      in FORTRAN.  Larger values mean a loosening of the convergence
      criterion.

    }
    \item{EqualTolerances}{
      Logical indicating whether each (quadratic) predictor will have
      equal tolerances. Setting \code{EqualTolerances=TRUE} can
      help avoid numerical problems, especially with binary data.
      Note that the estimated (common) tolerance matrix may or may not be
      positive-definite. If it is, then it can be scaled to the \eqn{R}
      x \eqn{R} identity matrix.  Setting \code{ITolerances=TRUE} will
      fit a common \eqn{R} x \eqn{R} identity matrix as the tolerance
      matrix to the data, but this is model-driven rather than being
      data-driven because it \emph{forces} bell-shaped curves/surfaces
      onto the data.  If the estimated (common) tolerance matrix happens
      to be positive-definite, then this model is essentially equivalent
      to the model with \code{ITolerances=TRUE}.
      See \bold{Details} in \code{\link{cqo}} and \code{\link{qrrvglm.control}}
      for more details.

    }
  \item{Etamat.colmax}{
    Positive integer, no smaller than \code{Rank}.  Controls the amount
    of memory used by \code{.Init.Poisson.QO()}.  It is the maximum
    number of columns allowed for the pseudo-response and its weights.
    In general, the larger the value, the better the initial value.
    Used only if \code{Use.Init.Poisson.QO=TRUE}.

  }

  \item{GradientFunction}{ 
   Logical. Whether \code{\link[stats]{optim}}'s argument \code{gr} is
   used or not, i.e., to compute gradient values.  The default value is
   usually faster on most problems.

  }
  \item{Hstep}{ 
   Positive value. Used as the step size in the finite difference
   approximation to the derivatives by \code{\link[stats]{optim}}.

  }
  \item{isdlv}{
   Initial standard deviations for the latent variables (site scores).
   Numeric, positive and of length \eqn{R} (recycled if necessary).
   This argument is used only if \code{ITolerances=TRUE}.  Used by
   \code{.Init.Poisson.QO()} to obtain initial values for the constrained
   coefficients \eqn{C} adjusted to a reasonable value. It adjusts the
   spread of the site scores relative to a common species tolerance of 1
   for each ordination axis.  A value between 0.5 and 10 is recommended;
   a value such as 10 means that the range of the environmental space is
   very large relative to the niche width of the species.  The successive
   values should decrease because the first ordination axis should have
   the most spread of site scores, followed by the second ordination
   axis, etc.

 }
  \item{ITolerances}{
   Logical. If \code{TRUE} then the (common) tolerance matrix is
   the \eqn{R} x \eqn{R} identity matrix by definition.  Note that
   \code{ITolerances=TRUE} implies \code{EqualTolerances=TRUE}, but
   not vice versa.  Internally, the quadratic terms will be treated
   as offsets (in GLM jargon) and so the models can potentially be
   fitted very efficiently. 
   See \bold{Details} in \code{\link{cqo}} and \code{\link{qrrvglm.control}}
   for more details.
   more details.  The success of \code{ITolerances=TRUE} often depends
   on suitable values for \code{isdlv} and/or \code{MUXfactor}.

 }
 \item{lvstart}{ 
   Optional matrix of initial values of the site scores. If given, the
   matrix must be \eqn{n} by \eqn{R}, where \eqn{n} is the number
   of sites and \eqn{R} is the rank.  This argument overrides the
   arguments \code{Use.Init.Poisson.QO} and \code{CA1}.
   Good possibilities for \code{lvstart} are the site scores from a
   constrained ordination, e.g., from \code{\link{cqo}}.

  }
  \item{jitter.sitescores}{ Logical.
   If \code{TRUE} the initial values for the site scores are jittered
   to add a random element to the starting values.

  }

  \item{maxitl}{ 
    Positive integer.  Number of iterations allowed for the IRLS algorithm
    implemented in the compiled code.

    }
  \item{Maxit.optim}{ 
    Positive integer.  Number of iterations given to the function
    \code{\link[stats]{optim}} at each of the \code{optim.maxit}
    iterations.

  }
  \item{MUXfactor}{
   Multiplication factor for detecting large offset values.  Numeric,
   positive and of length \eqn{R} (recycled if necessary).  This argument
   is used only if \code{ITolerances=TRUE}. Offsets are \eqn{-0.5}
   multiplied by the sum of the squares of all \eqn{R} latent variable
   values. If the latent variable values are too large then this will
   result in numerical problems. By too large, it is meant that the
   standard deviation of the latent variable values are greater than
   \code{MUXfactor[r] * isdlv[r]} for \code{r=1:Rank} (this is why
   centering and scaling all the numerical predictor variables in
   \eqn{x_2} is recommended).  A value about 3 or 4 is recommended.
   If failure to converge occurs, try a slightly lower value.

}
  \item{optim.maxit}{ 
    Positive integer.  Number of times \code{\link[stats]{optim}}
    is invoked.

%   At iteration \code{i}, the \code{i}th value of \code{Maxit.optim}
%   is fed into \code{\link[stats]{optim}}.

  }
  \item{nRmax}{ 
    Positive integer.  If the number of parameters making up the latent
    variable values (\eqn{n} multiplied by \eqn{R}) is greater than this
    argument then a conjugate-gradients algorithm is used, otherwise a
    quasi-Newton algorithm is used by \code{\link[stats]{optim}}. The
    conjugate-gradients method is more suitable when the number of
    parameters is large because it requires less memory.

    }
  \item{SD.sitescores}{ Numeric. Standard deviation of the
    initial values of the site scores, which are generated from
    a normal distribution.
    }
% \item{Dzero}{ Integer vector specifying which squared terms
%     are to be zeroed. These linear predictors will correspond to
%     a RR-VGLM.
%     The values must be elements from the set \{1,2,\ldots,\eqn{M}\}.
%     Used only if \code{Quadratic=TRUE} and \code{FastAlgorithm=FALSE}.
% }
  \item{SmallNo}{ Positive numeric between \code{.Machine$double.eps} and
      \code{0.0001}.
      Used to avoid under- or over-flow in the IRLS algorithm.
  }
  \item{trace}{ Logical indicating if output should be produced for
    each iteration.
  }
%  \item{Kinit}{ Initial values for the index parameters \code{k} in the
%   negative binomial distribution (one per species).
%   In general, a smaller number is preferred over a larger number.
%   The vector is recycled to the number of responses (species).
% }
  \item{Use.Init.Poisson.QO}{
    Logical. If \code{TRUE} then the function \code{.Init.Poisson.QO()} is
    used to obtain initial values for the site scores.  If \code{FALSE}
    then random numbers are used instead.  Both \code{CA1} and
    \code{Use.Init.Poisson.QO} cannot both be \code{TRUE}.

  }
  \item{\dots}{ Ignored at present. }
}
\details{
   The algorithm currently used by \code{\link{uqo}} is unsophisticated
   and fails often. Improvements will hopefully be made soon.

   See \code{\link{cqo}} and \code{\link{qrrvglm.control}} for more details
   that are equally pertinent to UQO.

% zz site scores are centered. Possibly uncorrelated too?

   To reduce the number of parameters being estimated, setting
   \code{ITolerances = TRUE} or \code{EqualTolerances = TRUE} is advised.

}
\value{
  A list with the components corresponding to its arguments, after
  some basic error checking.
}
\references{


Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.


%Yee, T. W. (2012)
%On constrained and unconstrained quadratic ordination.
%\emph{Manuscript in preparation}.


}
\author{T. W. Yee}
\note{
  This is a difficult optimization problem, and the current
  algorithm needs to be improved.


}
\seealso{
  \code{\link{uqo}}.


}

\section{Warning }{
  This function is currently very sensitive to initial values. Setting
  \code{Bestof} some reasonably large integer is recommended.


}

\examples{
uqo.control()
}
\keyword{models}
\keyword{regression}


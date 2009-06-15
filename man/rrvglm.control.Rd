\name{rrvglm.control}
\alias{rrvglm.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Control function for rrvglm }
\description{
  Algorithmic constants and parameters for running \code{rrvglm} are set
  using this function.
}
\usage{
rrvglm.control(Rank = 1, Algorithm = c("alternating", "derivative"),
    Corner = TRUE, Uncor = FALSE, Wmat = NULL, Svd.arg = FALSE, 
    Index.corner = if (length(Structural.zero)) 
    head((1:1000)[-Structural.zero], Rank) else 1:Rank,
    Alpha = 0.5, Bestof = 1, Cinit = NULL,
    Etamat.colmax = 10,
    SD.Cinit = 0.02, Structural.zero = NULL,
    Norrr = ~1, trace = FALSE, Use.Init.Poisson.QO = FALSE, 
    checkwz = TRUE, wzepsilon = .Machine$double.eps^0.75, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{Rank}{
    The numerical rank \eqn{R} of the model.
    Must be an element from the set \{1,2,\ldots,min(\eqn{M},\emph{p2})\}. 
    Here, the vector of explanatory variables \bold{x} is partitioned into
    (\bold{x1},\bold{x2}), which is of dimension \emph{p1}+\emph{p2}.
    The variables making up \bold{x1} are given by the terms in
    \code{Norrr} argument, and the rest of the terms comprise \bold{x2}.

  }
  \item{Algorithm}{
    Character string indicating what algorithm is
    to be used. The default is the first one.

  }
  \item{Corner}{
    Logical indicating whether corner constraints are
    to be used. This is one method for ensuring a unique solution.
    If \code{TRUE}, \code{Index.corner} specifies the \eqn{R} rows
    of the constraint matrices that are use as the corner constraints,
    i.e., they hold an order-\eqn{R} identity matrix.

  }
\item{Uncor}{
  Logical indicating whether uncorrelated latent variables are to be used.
  This is normalization forces the variance-covariance
  matrix of the latent variables to be \code{diag(Rank)}, i.e., unit
  variance and uncorrelated. This constraint does not lead to a unique
  solution because it can be rotated.

}
  \item{Wmat}{ Yet to be done. }
  \item{Svd.arg}{
    Logical indicating whether a singular value decomposition
    of the outer product is to computed.  This is another normalization
    which ensures uniqueness.  See the argument \code{Alpha} below.

}
  \item{Index.corner}{
    Specifies the \eqn{R} rows of the constraint matrices that are
    used for the corner constraints, i.e., they hold an order-\eqn{R}
    identity matrix.

  }
  \item{Alpha}{
    The exponent in the singular value decomposition that is used in
    the first part: if the SVD is
    \eqn{U D V^T}{ U \%*\% D \%*\% t(V) } then the
    first and second parts are 
    \eqn{U D^{\alpha}}{ U \%*\% D^Alpha} 
    and
    \eqn{D^{1-\alpha}  V^T}{D^(1-Alpha) \%*\% t(V)} respectively.
    A value of 0.5 is `symmetrical'. 
    This argument is used only when \code{Svd.arg=TRUE}.

  }
  \item{Bestof}{
    Integer. The best of \code{Bestof} models fitted is
    returned. This argument helps guard against local solutions by
    (hopefully) finding the global solution from many fits. The
    argument works only when the function generates its own initial value
    for \bold{C},
    i.e., when \bold{C} is \emph{not} passed in as initial values.

  }
  \item{Cinit}{
    Initial \bold{C} matrix which may speed up convergence.
    It must be of the correct dimension.

  }
  \item{Etamat.colmax}{
    Positive integer, no smaller than \code{Rank}.  Controls the amount
    of memory used by \code{.Init.Poisson.QO()}.  It is the maximum
    number of columns allowed for the pseudo-response and its weights.
    In general, the larger the value, the better the initial value.
    Used only if \code{Use.Init.Poisson.QO=TRUE}.

  }

% \item{Quadratic}{
%   Logical indicating whether a \emph{Quadratic}
%   RR-VGLM is to be fitted. If \code{TRUE}, an object of class
%   \code{"qrrvglm"} will be returned, otherwise \code{"rrvglm"}.
% }
  \item{Structural.zero}{
      Integer vector specifying which rows
      of the constraint matrices are to be all zeros.
  }
  \item{SD.Cinit}{
      Standard deviation of the initial values for the elements
      of \bold{C}.
      These are normally distributed with mean zero.  
      This argument is used only if \code{Use.Init.Poisson.QO = FALSE}.
  }
% \item{ppar}{ Ignore this. }
  \item{Norrr}{
    Formula giving terms that are not to be included in the reduced-rank
    regression. These variables constitute the \eqn{\bold{B}_1}{\bold{B}1}
    matrix in the Yee and Hastie paper.  Those variables which
    are subject to the reduced-rank regression correspond to the
    \eqn{\bold{B}_2}{\bold{B}2} matrix.

  }
  \item{trace}{
    Logical indicating if output should be produced for
    each iteration.

%   Useful when \code{Quadratic=TRUE} because QRR-VGLMs are
%   computationally expensive and it's good to see that the program
%   is working!

  }
  \item{Use.Init.Poisson.QO}{
    Logical indicating whether the \code{.Init.Poisson.QO()} should
    be used to obtain initial values for the \bold{C}.  The function
    uses a new method that can work well if the data are Poisson counts
    coming from an equal-tolerances QRR-VGLM (CQO).  This option is less
    realistic for RR-VGLMs compared to QRR-VGLMs.

  }
  \item{checkwz}{
    logical indicating whether the diagonal elements of
    the working weight matrices should be checked whether they are
    sufficiently positive, i.e., greater than \code{wzepsilon}. If not,
    any values less than \code{wzepsilon} are replaced with this value.

  }
  \item{wzepsilon}{
  Small positive number used to test whether the diagonals of the working
  weight matrices are sufficiently positive.

  }
  \item{\dots}{
    Variables in \dots are passed into
    \code{\link{vglm.control}}. If the derivative algorithm is used, then
    \dots are also passed into \code{\link{rrvglm.optim.control}}.
  }
  In the above, \eqn{R} is the \code{Rank} and 
  \eqn{M} is the number of linear predictors.
  
}
\details{
%   QRR-VGLMs are an extension of RR-VGLMs and are useful for constrained
%   ordination. QRR-VGLMs fitted with \pkg{VGAM} allow a maximum
%   likelihood solution to constrained quadratic ordination (CQO;
%   formerly called canonical Gaussian ordination) models.

%   For QRR-VGLMs, if \code{EqualTolerances=TRUE} and
%   \code{ITolerances=FALSE} then the default is that the \bold{C}
%   matrix is constrained by forcing the latent variables to have sample
%   variance-covariance matrix equalling \code{diag(Rank)}, i.e., unit
%   variance and uncorrelated.

    \pkg{VGAM} supports three normalizations to ensure a unique
    solution. Of these, only corner constraints will work with
    \code{summary} of RR-VGLM objects.

}
\value{
  A list with components matching the input names.
  Some error checking is done, but not much.
}
\references{
  Yee, T. W. and Hastie, T. J. (2003)
  Reduced-rank vector generalized linear models.
  \emph{Statistical Modelling},
  \bold{3}, 15--41.

}
\author{ Thomas W. Yee }
\note{

% The function call \code{cqo(...)} is equivalent to
% \code{rrvglm(..., Quadratic=TRUE)}, and hence uses this function.
% For QRR-VGLMs, the function \code{\link{qrrvglm.control}} is called too.

  The arguments in this function begin with an upper case letter to help
  avoid interference with those of \code{\link{vglm.control}}.

  In the example below a rank-1 \emph{stereotype} model (Anderson, 1984)
  is fitted.

}

%- \section{Warning }{ }

\seealso{
  \code{\link{rrvglm}},
  \code{\link{rrvglm.optim.control}},
  \code{\link{rrvglm-class}},
  \code{\link{vglm}},
  \code{\link{vglm.control}},
  \code{\link{cqo}}.
}

\examples{
set.seed(111)
pneumo = transform(pneumo, let=log(exposure.time),
            x1 = runif(nrow(pneumo)))  # x1 is some unrelated covariate
fit = rrvglm(cbind(normal, mild, severe) ~ let + x1,
             multinomial, pneumo, Rank=1, Index.corner=2)
constraints(fit)
vcov(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}


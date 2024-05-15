\name{rrvglm.control}
\alias{rrvglm.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Control Function for rrvglm() }
\description{
  Algorithmic constants and parameters for
  running \code{rrvglm} are set using this
  function.
  Doubly constrained RR-VGLMs (DRR-VGLMs) are
  also catered for.
  
}
\usage{
rrvglm.control(Rank = 1, Corner = TRUE,
    Index.corner = head(setdiff(seq(length(str0) +
    Rank), str0), Rank), noRRR = ~ 1, str0 = NULL,
    Crow1positive = NULL, trace = FALSE, Bestof = 1,
    H.A.thy = list(), H.C = list(),
    Ainit = NULL, Cinit = NULL, sd.Cinit = 0.02,
    Algorithm = "alternating", Etamat.colmax = 10,
    noWarning = FALSE, Use.Init.Poisson.QO = FALSE,
    checkwz = TRUE, Check.rank = TRUE, Check.cm.rank = TRUE,
    wzepsilon = .Machine$double.eps^0.75, ...)
}
%- maybe also `usage' for other objects documented here.
%
%   noRRR = ~ 1, str0 = NULL, Crow1positive = TRUE,
% 20240327; Defunct arguments:
%   scaleA = FALSE,
%   Alpha = 0.5,
%   Uncorrelated.latvar = FALSE, Wmat = NULL,
%   Svd.arg = FALSE, Norrr = NA,
\arguments{
  \item{Rank}{
    The numerical rank \eqn{R} of the
    model.  Must be an element from the set
    \{1,2,\ldots,min(\eqn{M},\emph{p2})\}.
    Here, the vector of explanatory
    variables \bold{x} is partitioned into
    (\bold{x1}, \bold{x2}), which is of dimension
    \emph{p1}+\emph{p2}.  The variables making
    up \bold{x1} are given by the terms in
    \code{noRRR} argument, and the rest of the
    terms comprise \bold{x2}.


  }
  \item{Corner}{
    Logical indicating whether corner
    constraints are to be used.
    Strongly recommended as the \emph{only}
    method for fitting RR-VGLMs and
    DRR-VGLMs.
    This is one
    method for ensuring a unique solution
    and the availability of standard errors.
    If \code{TRUE} then \code{Index.corner}
    specifies the position of the corner
    constraints.
    Also used by DRR-VGLMs, and this is
    known as RCCs (see below).


  }
  \item{Index.corner}{
    Specifies the \eqn{R} rows of
    the corner constraints, e.g., they hold an
    order-\eqn{R} identity matrix for RR-VGLMs.
    The default are the first \eqn{R} positive
    integers that exclude \code{str0}.


    
% 20231223:

For DRR-VGLMs one needs
to have (restricted) corner constraints.
Then argument \code{str0} will be ignored.
If there are structural zeros then they
should be built into the constraint matrices
(see \code{H.A.thy}).


}
  \item{noRRR}{
    Formula giving terms that are \emph{not}
    to be included in the reduced-rank
    regression.  That is, \code{noRRR}
    specifes which explanatory variables
    are in the \bold{x1} vector of
    \code{\link{rrvglm}}, and the rest go into
    \bold{x2}.
    I will write \bold{x1} as \eqn{x_1}{x1}
    and \bold{x2} as \eqn{x_2}{x2} hereafter.
    The \eqn{x_1}{x1} variables
    constitute the \eqn{\bold{B}_1}{\bold{B}1}
    matrix in Yee and Hastie (2003).  Those
    \eqn{x_2}{x2} variables which are subject
    to the reduced-rank regression correspond
    to the \eqn{\bold{B}_2}{\bold{B}2}
    matrix.  Set \code{noRRR = NULL} for the
    reduced-rank regression to be applied to
    every explanatory variable including the
    intercept.


  }
  \item{str0}{
  Integer vector specifying which rows of the
  estimated constraint matrices (\bold{A})
  are to be all zeros.  These are called
  \emph{structural zeros}.  Must not have
  any common value with \code{Index.corner},
  and be a subset of the vector \code{1:M}.
  The default, \code{str0 = NULL}, means no
  structural zero rows at all.
  This argument is ignored by DRR-VGLMs.



  }
  \item{Crow1positive}{
    Currently this argument has no effect.
    In the future, it may be a
    logical vector of length \code{Rank}
    (recycled if necessary): are
    the elements of the first
    row of \eqn{C} positive?
    Because of corner constraints (\code{diag(Rank)})
    this argument currently must be \code{NULL}.
    Having a logical vector would mean some of the diagonals
    might be \eqn{1} and others \eqn{-1}. 


%Possibly this might be relaxed in the future or
%else the argument withdrawn.




% 20240327; making this argument a lame duck.
%   For example,
%   if \code{Rank} is 4, then
%   specifying \code{Crow1positive = c(FALSE,
%   TRUE)} will force \eqn{C[1,1]}
%   and \eqn{C[1,3]}
%   to be negative,
%   and \eqn{C[1,2]} and \eqn{C[1,4]}
%   to be positive.
%   This argument therefore
%   allows the user to determine the direction
%   of the latent variables since they are
%   unique up to a sign.
%   This argument certainly works for RR-VGLMs
%   but may not be applicable to DRR-VGLMs
%   because the constraint matrices may 
%   control their sign.


  }


  \item{trace}{
    Logical indicating if output should be
    produced for each iteration.


%Useful when \code{Quadratic=TRUE} because
% QRR-VGLMs are computationally expensive and
% it's good to see that the program is working!

  }




  \item{Bestof}{
    Integer. The best of \code{Bestof}
    models fitted is returned. This argument
    helps guard against local solutions by
    (hopefully) finding the global solution
    from many fits. The argument works only
    when the function generates its own initial
    value for \bold{C}, i.e., when \bold{C}
    is \emph{not} passed in as initial values.


  }
  \item{H.A.thy, H.C}{
  Lists.
  DRR-VGLMs are \emph{Doubly constrained}
  RR-VGLMs where \bold{A} has
  \code{Rank}
  constraint matrices
  (one for each column)
  in a list called \code{H.A.thy},
  and \bold{C} has a
  constraint matrix for each row,
  i.e., for each explanatory
  variable making up it.
  The class
  \code{"drrvglm"} may arise for such models.
  So \code{H.C}
  should be a named list of \eqn{p_2}
  constraint matrices,
  each one for a different row of \bold{C},
  i.e., \eqn{p_2} is the number of variables
  making up the latent variable.
  Note that if \code{H.C} has
  \code{\link[base:names]{names}} then
  matching is done with that,
  and the components of \code{H.C}
  are reordered if they are not sorted
  according to the terms in \code{formula}.
  If they are not named, then
  their \emph{order} is used,
  for example, \code{H.C[[1]]}
  and \code{H.C[[2]]} are taken as
  the constraint matrices for the first two
  variables of the latent variable(s).


% keyword: setp2 and Pos456


% Argument \code{drrvglm} is a logical and
% needs to be set \code{TRUE} for \code{H.A.thy}
% and \code{H.C} to work.  Then \code{H.C}
% is a named list of constraint matrices,

  }
  \item{Algorithm}{
    Character string indicating what algorithm is
    to be used. The default is the first one.
    The choice \code{"derivative"} has been
    withdrawn in \pkg{VGAM} 1.1-10.
    

  }
  \item{Ainit, Cinit}{
    Initial \bold{A} and \bold{C} matrices which
    may speed up convergence.  They must be of
    the correct dimension.


  }
  \item{sd.Cinit}{
  Standard deviation of the initial values
  for the elements of \bold{C}.
  These are normally distributed with
  mean zero.  This argument is used only if
  \code{Use.Init.Poisson.QO = FALSE}.


  }
  \item{Etamat.colmax}{
    Positive integer, no smaller than
    \code{Rank}.  Controls the amount of
    memory used by \code{.Init.Poisson.QO()}.
    It is the maximum number of columns allowed
    for the pseudo-response and its weights.
    In general, the larger the value, the
    better the initial value.  Used only if
    \code{Use.Init.Poisson.QO = TRUE}.


  }
%\item{Uncorrelated.latvar}{
%  Logical indicating whether uncorrelated
%  latent variables are to be used.  This is
%  normalization forces the variance-covariance
%  matrix of the latent variables to be
%  \code{diag(Rank)}, i.e., unit variance and
%  uncorrelated. This constraint does not lead to
%  a unique solution because it can be rotated.
%  \bold{Update during 2023/2024:}
%  setting this argument to be \code{TRUE}
%  is not recommended and might not work anymore.
%}


% \item{Wmat}{ Yet to be done. }
% \item{Svd.arg}{
%   Logical indicating whether a singular value
%   decomposition of the outer product is to
%   computed.  This is another normalization
%   which ensures uniqueness.  See the argument
%   \code{Alpha} below.
% \bold{Update during 2023/2024:}
% setting this argument to be \code{TRUE}
% is not recommended and might not work anymore.
%}


% \item{Alpha}{
%   The exponent in the singular value
%   decomposition that is used in the first
%   part: if the SVD is
%   \eqn{U D V^T}{ U \%*\% D \%*\% t(V) } then the
%   first and second parts are
%   \eqn{U D^{\alpha}}{ U \%*\% D^Alpha}
%   and
%   \eqn{D^{1-\alpha}  V^T}{D^(1-Alpha) \%*\% t(V)}
%   respectively.
%   A value of 0.5 is `symmetrical'.
%   This argument is used only when
%   \code{Svd.arg=TRUE}.
% \bold{Update during 2023/2024:}
% using this argument
% is not recommended and might not work anymore.
% }

% \item{Quadratic}{
%   Logical indicating whether a \emph{Quadratic}
%   RR-VGLM is to be fitted. If \code{TRUE}, an object of class
%   \code{"qrrvglm"} will be returned, otherwise \code{"rrvglm"}.
% }
% \item{Norrr}{
% Defunct. Please use \code{noRRR}.
% Its use results in an error.
% The argument may be removed soon.
% }
% \item{ppar}{ Ignore this. }
  \item{Use.Init.Poisson.QO}{
    Logical indicating whether the
    \code{.Init.Poisson.QO()} should be used
    to obtain initial values for the \bold{C}.
    The function uses a new method that can
    work well if the data are Poisson counts
    coming from an equal-tolerances QRR-VGLM
    (CQO).  This option is less realistic for
    RR-VGLMs compared to QRR-VGLMs.

  }
  \item{checkwz}{
    logical indicating whether the diagonal
    elements of the working weight matrices
    should be checked whether they are
    sufficiently positive, i.e., greater than
    \code{wzepsilon}. If not, any values less
    than \code{wzepsilon} are replaced with
    this value.


  }
  \item{noWarning, Check.rank, Check.cm.rank}{
    Same as \code{\link{vglm.control}}.
    Ignored for \pkg{VGAM} 0.9-7 and higher.



  }
  \item{wzepsilon}{
  Small positive number used to test whether
  the diagonals of the working weight matrices
  are sufficiently positive.

  }
% \item{scaleA}{
% Logical.
% This is not recommended because it has been
% replaced by
% \emph{restricted corner constraints}
% (RCCs)---see
% \code{\link{rrvglm}}.



  
% Another uniqueness constraint to obtain a
% unique \bold{A} and \bold{C}.
% If \code{H.A.thy} and/or \code{H.C} are
% inputted then sometimes one wants to
% preserve the structure in \bold{A},
% e.g., \code{\link{CM.qnorm}}.
% Here, \code{A <- scale(A, center = FALSE)}
% so that only the columns are multiplicatively
% scaled.
% Note that the estimates of the elements of
% \bold{A} and \bold{C} are unique, up to
% their sign.
% Also note that ideally the attributes
% \code{attr(,"scaled:scale")} should
% be unity upon convergence so that
% if they differ substantially from that
% then this suggests some misbehaviour
% in convergence.
% }



  \item{\dots}{
    Variables in \dots are passed into
    \code{\link{vglm.control}}.
    If the derivative algorithm is used
    then \dots are also passed into
    \code{\link{rrvglm.optim.control}};
    and if the alternating algorithm is
    used then \dots are also passed into
    \code{\link{valt0.control}}.


  }
  In the above, \eqn{R} is the \code{Rank} and
  \eqn{M} is the number of linear predictors.

}
\details{
% QRR-VGLMs are an extension of RR-VGLMs and are useful for constrained
% ordination. QRR-VGLMs fitted with \pkg{VGAM} allow a maximum
% likelihood solution to constrained quadratic ordination (CQO;
% formerly called canonical Gaussian ordination) models.


% For QRR-VGLMs, if \code{eq.tolerances=TRUE} and
% \code{I.tolerances=FALSE} then the default is that the \bold{C} matrix
% is constrained by forcing the latent variables to have sample
% variance-covariance matrix equalling \code{diag(Rank)}, i.e., unit
% variance and uncorrelated.



  \pkg{VGAM} supported three normalizations
  to ensure a unique solution.
  But currently,
  only corner constraints will work with
  \code{summary} of RR-VGLM
  and DRR-VGLM objects.
  \bold{Update during late-2023/early-2024:}
  with ongoing work implementing
  the \code{"drrvglm"} class, there may
  be disruption and changes to other
  normalizations. However, corner
  constraints should be fully supported
  and have the greatest priority.



}
\value{
  A list with components matching the input
  names.  Some error checking is done, but
  not much.


}
\references{
  Yee, T. W. and Hastie, T. J. (2003).
  Reduced-rank vector generalized linear models.
  \emph{Statistical Modelling},
  \bold{3}, 15--41.

}
\author{ Thomas W. Yee }
\note{

% The function call \code{cqo(...)} is equivalent to
% \code{rrvglm(..., Quadratic=TRUE)}, and hence uses this function.
% For QRR-VGLMs, the function \code{\link{qrrvglm.control}} is called too.



  In \pkg{VGAM} 1.1-11 and higher,
  the following arguments are no longer supported:
  \code{Wmat}, \code{Norrr}, \code{Svd.arg},
  \code{Uncorrelated.latvar}, \code{scaleA}.
  Users should use corner constraints only.



  The arguments in this function begin with an
  upper case letter to help avoid interference
  with those of \code{\link{vglm.control}}.


  In the example below a rank-1
  \emph{stereotype} model (Anderson, 1984)
  is fitted,
  however, the intercepts are completely
  unconstrained rather than sorted.


}

%- \section{Warning }{ }

\seealso{
  \code{\link{rrvglm}},
  \code{\link{rrvglm-class}},
  \code{\link{summary.drrvglm}},
  \code{\link{rrvglm.optim.control}},
  \code{\link{vglm}},
  \code{\link{vglm.control}},
  \code{\link{TypicalVGAMfamilyFunction}},
  \code{\link{CM.qnorm}},
  \code{\link{cqo}}.


}

\examples{
\dontrun{
set.seed(111)
pneumo <- transform(pneumo, let = log(exposure.time),
                            x3 = runif(nrow(pneumo)))  # Unrelated
fit <- rrvglm(cbind(normal, mild, severe) ~ let + x3,
              multinomial, pneumo, Rank = 1, Index.corner = 2)
constraints(fit)
vcov(fit)
summary(fit)
}
}
\keyword{optimize}
\keyword{models}
\concept{Reduced-Rank Vector Generalized Linear Model}
%\keyword{regression}

% Was:
% if (length(str0))  head((1:1000)[-str0], Rank) else 1:Rank,
% vglm.control(Rank = 1, Algorithm = c("alternating", "derivative"),







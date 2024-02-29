\name{CM.equid}
\alias{CM.equid}
\alias{CM.free}
\alias{CM.ones}
\alias{CM.symm0}
\alias{CM.symm1}
\alias{CM.qnorm}
\alias{CM.qlogis}
%\alias{CM.symm1}
%\alias{CM.symm1}
%
\title{ Constraint Matrices for Symmetry, Order,
    Parallelism, etc.
%%  ~~function to do ... ~~
}
\description{
    Given \emph{M} linear/additive predictors,
    construct the constraint matrices to allow
    symmetry, (linear and normal) ordering, etc.
    in terms such as the intercept.

}
\usage{
CM.equid(M, Trev = FALSE, Tref = 1)
 CM.free(M, Trev = FALSE, Tref = 1)
 CM.ones(M, Trev = FALSE, Tref = 1)
CM.symm0(M, Trev = FALSE, Tref = 1)
CM.symm1(M, Trev = FALSE, Tref = 1)
CM.qnorm(M, Trev = FALSE, Tref = 1)
}
%- maybe also 'usage' 4 other objs docted here.
\arguments{
  \item{M}{ Number of linear/additive predictors,
        usually \eqn{>1}.

}
\item{Tref}{ Reference level for the threshold,
  this should be a single value from \code{1:M}.
  This argument is ignored by some of the above
  functions.


}
\item{Trev}{ Logical.
  Apply reverse direction
  for the thresholds direction?
  This argument is ignored by some of the above
  functions.

}




}
\details{
  A constraint matrix is \eqn{M \times R} where
  \eqn{R} is its rank and usually the elements are
  0, 1 or \eqn{-1}.
  There is a constraint matrix for each column
  of the LM matrix used to fit the
  \code{\link{vglm}}.
  They are used to apportion the regression
  coefficients to the linear predictors, e.g.,
  parallelism, exchangeability, etc.
  The functions described here are intended
  to construct
  constraint matrices easily for
  symmetry constraints and
  linear ordering etc.
  They are potentially useful for categorical data
  analysis (e.g., \code{\link{cumulative}},
  \code{\link{multinomial}}), especially for the
  intercept term.
  When applied to \code{\link{cumulative}},
  they are sometimes called
  \emph{structured thresholds},
  e.g., \pkg{ordinal}.



  One example is the stereotype model proposed
  by Anderson (1984)
  (see \code{\link{multinomial}} and
  \code{\link{rrvglm}}) where the elements of
  the \bold{A} matrix are ordered.
  This is not fully possible in \pkg{VGAM} but
  some special cases can be fitted, e.g.,
  use \code{\link{CM.equid}} to create
  a linear ordering.
  And \code{\link{CM.symm1}} might result in
  fully ordered estimates too, etc.



  \code{\link{CM.free}} creates
  \emph{free} or unconstrained estimates.
  It is almost always the case for VGLMs,
  and is simply \code{diag(M)}.



  \code{\link{CM.ones}} creates
  \emph{equal} estimates,
  which is also known as the \emph{parallelism}
  assumption in models such as
  \code{\link{cumulative}}.
  It gets its name because the constraint matrix
  is simply \code{matrix(1, M, 1)}.




  \code{\link{CM.equid}} creates
  \emph{equid}istant estimates. This is a
  linear scaling, and the direction and
  origin are controlled by \code{Treverse}
  and \code{Tref} respectively.



  \code{\link{CM.qnorm}} and
  \code{\link{CM.qlogis}} are based on
  \code{\link[stats]{qnorm}} and
  \code{\link[stats]{qlogis}}.
  For example, \code{CM.qnorm(M)} is essentially
  \code{cbind(qnorm(seq(M) / (M + 1)))}.
  This might be useful with a model with
  \code{\link{probitlink}} applied to multiple
  intercepts.
  


  Further details can be found at
  \code{\link{cumulative}} and
  \code{\link{CommonVGAMffArguments}},



%% If nec, more details than the descrip above
}
\value{
    A constraint matrix.


}
%\references{
%% ~put ref to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% \section{Warning }{....}

\seealso{
  \code{\link{CommonVGAMffArguments}},
  \code{\link{cumulative}},
  \code{\link{acat}},
  \code{\link{cratio}},
  \code{\link{sratio}},
  \code{\link{multinomial}}.

}
\examples{
CM.equid(4)
CM.equid(4, Trev = TRUE, Tref = 3)
CM.symm1(5)
CM.symm0(5)
CM.qnorm(5)
}
%\keyword{regression}
\keyword{models}
\keyword{utilities}
% Add std keywords, see file 'KEYWORDS' in the
% (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-std keywords, use \concept instead:
\concept{Constraint matrices}
% \concept{ ~cpt2 }
% Use only one concept per line.

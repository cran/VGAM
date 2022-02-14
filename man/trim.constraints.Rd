\name{trim.constraints}
\alias{trim.constraints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Trimmed Constraint Matrices }
\description{
  Deletes statistically nonsignficant regression coefficients via
  their constraint matrices, for future refitting.



}
\usage{
trim.constraints(object, sig.level = 0.05, max.num = Inf,
                 intercepts = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Some \pkg{VGAM} object, especially having
  class \code{\link{vglmff-class}}.
  It has not yet been tested on non-\code{"vglm"} objects.




  }
  \item{sig.level}{
    Significance levels, with values in \eqn{[0, 1]}.
    Columns of constraint matices whose p-values are larger than
    this argument are deleted.
    With terms that generate more than one column of the
    \code{"lm"} model matrix, all p-values must be greater
    than this argument for deletion.
    This argument is recycled to the total number of
    regression coefficients of \code{object}.


  }
  \item{max.num}{
    Numeric, positive and integer-valued.
    Maximum number of regression coefficients allowable for deletion.
    This allows one to limit the number of deleted coefficients.
    For example,
    if \code{max.num = 1} then only the largest p-value is used
    for the deletion,
    provided it is larger than \code{sig.level}.
    The default is to delete all those coefficients whose
    p-values are greater than \code{sig.level}.
    With a finite value,
    this argument will probably not work properly
    when there are terms that
    generate more than one column of the LM model matrix.
    Having a value greater than unity might be unsuitable
    in the presence of multicollinearity because all
    correlated variables might be eliminated at once.



    
%   ---this corresponds
%   to the \emph{backward elimination} technique if \eqn{M=1}.


    
  }
  \item{intercepts}{
  Logical. Trim the intercept term?
  If \code{FALSE} then the constraint matrix for the
  \code{"(Intercept)"} term is left unchanged.



  }
  \item{\dots}{
  Unused but for provision in the future.


  }

}
\details{
  This utility function is intended to simplify an existing
  \code{\link{vglm}} object having
  variables (terms) that affect unnecessary parameters.
  Suppose the explanatory variables in the formula
  includes a simple numeric covariate called \code{x2}.
  This variable will affect every linear predictor if
  \code{zero = NULL} in the \pkg{VGAM} family function.
  This situation may correspond to the constraint matrices having
  unnecessary columns because their regression coefficients are
  statistically nonsignificant.
  This function attempts to delete those columns and
  return a possibly simplified list of constraint matrices
  that can make refitting a simpler model easy to do.
  P-values obtained from \code{\link{summaryvglm}}
  (with \code{HDEtest = FALSE} for increased speed)
  are compared to \code{sig.level} to test for
  statistical significance.


  
  For terms that generate more than one column of the
  \code{"lm"} model matrix,
  such as \code{\link[splines]{bs}} and \code{\link[stats]{poly}},
  the column is deleted if all regression coefficients
  are statistically nonsignificant.
  Incidentally, users should instead use
  \code{\link{sm.bs}},
  \code{\link{sm.ns}},
  \code{\link{sm.poly}},
  etc.,
  for smart and safe prediction.


  
  One can think of this function as facilitating
  \emph{backward elimination} for variable selection,
  especially if \code{max.num = 1} and \eqn{M=1},
  however usually more than one regression coefficient is deleted
  here by default.



% It would be possible to adjust \code{sig.level} manually
% so that only one column is deleted at a time, however this would be
% labourious.
 
  

}
\value{
  A list of possibly simpler constraint matrices
  that can be fed back into the model using the
  \code{constraints} argument
  (usually \code{zero = NULL} is needed to avoid a warning).
  Consequently, they are required to be of the \code{"term"}-type.
  After the model is refitted, applying
  \code{\link{summaryvglm}} should result in
  regression coefficients that are `all' statistically
  significant.



}
\author{T. W. Yee }
\section{Warning }{
  This function has not been tested thoroughly.
  One extreme is that a term is totally deleted because
  none of its regression coefficients are needed,
  and that situation has not yet been finalized.
  Ideally, \code{object} only contains terms where at least
  one regression coefficient has a p-value less than
  \code{sig.level}.
  For ordered factors and other situations, deleting
  certain columns may not make sense and destroy interpretability.



  As stated above, \code{max.num} may not work properly
  when there are terms that
  generate more than one column of the LM model matrix.
  However, this limitation may change in the future.




}

\note{
  This function is experimental and may be replaced by
  some other function in the future.
  This function does not use S4 object oriented programming
  but may be converted to such in the future.


}

%\references{
%}


\seealso{
  \code{\link{constraints}},
  \code{\link{vglm}},
  \code{\link{summaryvglm}},
  \code{\link{model.matrixvlm}},
  \code{\link{drop1.vglm}},
  \code{\link{step4vglm}},
  \code{\link{sm.bs}},
  \code{\link{sm.ns}},
  \code{\link{sm.poly}}.



}
\examples{
\dontrun{ data("xs.nz", package = "VGAMdata")
fit1 <-
  vglm(cbind(worry, worrier) ~ bs(age) + sex + ethnicity + cat + dog,
       binom2.or(zero = NULL), data = xs.nz, trace = TRUE)
summary(fit1, HDEtest = FALSE)  # 'cat' is not significant at all
dim(constraints(fit1, matrix = TRUE))
(tclist1 <- trim.constraints(fit1))  # No 'cat'
fit2 <-  # Delete 'cat' manually from the formula:
  vglm(cbind(worry, worrier) ~ bs(age) + sex + ethnicity +       dog,
       binom2.or(zero = NULL), data = xs.nz,
       constraints = tclist1, trace = TRUE)
summary(fit2, HDEtest = FALSE)  # A simplified model
dim(constraints(fit2, matrix = TRUE))  # Fewer regression coefficients
}}
\keyword{models}
\keyword{regression}


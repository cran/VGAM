\name{cumulative}
\alias{cumulative}
%\alias{scumulative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Ordinal Regression with Cumulative Probabilities }
\description{
  Fits a cumulative link
  regression model to a (preferably ordered) factor response.

}
\usage{
cumulative(link = "logitlink", parallel = FALSE,
    reverse = FALSE, multiple.responses = FALSE,
    ynames = FALSE, Thresh = NULL, Trev = reverse,
    Tref = if (Trev) "M" else 1, whitespace = FALSE)
}
%apply.parint = FALSE,
%scumulative(link = "logitlink",
%            lscale = "loglink", escale = list(),
%            parallel = FALSE, sparallel = TRUE, reverse = FALSE,
%            iscale = 1)
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{link}{
  Link function applied to the \eqn{J} cumulative probabilities.
  See \code{\link{Links}} for more choices,
  e.g., for the cumulative
\code{\link{probitlink}}/\code{\link{clogloglink}}/\ldots
  models.


  }
% \item{lscale}{
% Link function applied to the \eqn{J} scaling parameters.
% See \code{\link{Links}} for more choices.
%
% }


  \item{parallel}{
  A logical or formula specifying which terms have
  equal/unequal coefficients.
  See below for more information about the parallelism
  assumption.
  The default results in what some people call the
  \emph{generalized ordered logit model} to be fitted.
  If \code{parallel = TRUE} then it does not apply to the
  intercept.


  
  The \emph{partial proportional odds model} can be
  fitted by assigning this argument something like
  \code{parallel = TRUE ~ -1 + x3 + x5} so that
  there is one regression coefficient
  for \code{x3} and \code{x5}.
  Equivalently, setting
  \code{parallel = FALSE ~ 1 + x2 + x4} means \eqn{M}
  regression coefficients for the intercept
  and \code{x2} and \code{x4}.
  It is important that the intercept is never parallel.
  See \code{\link{CommonVGAMffArguments}} for more
  information.




  }
% \item{sparallel}{
% For the scaling parameters.
% A logical, or formula specifying which terms have
% equal/unequal coefficients.
% This argument is not applied to the intercept.
% The \code{scumulative()} function requires covariates; for
% intercept models use \code{cumulative()}.

% }
  \item{reverse}{
  Logical.
  By default, the cumulative probabilities used are
  \eqn{P(Y\leq 1)}{P(Y<=1)}, \eqn{P(Y\leq 2)}{P(Y<=2)},
  \dots, \eqn{P(Y\leq J)}{P(Y<=J)}.
  If \code{reverse} is \code{TRUE} then
  \eqn{P(Y\geq 2)}{P(Y>=2)}, \eqn{P(Y\geq 3)}{P(Y>=3)}, \dots,
  \eqn{P(Y\geq J+1)}{P(Y>=J+1)} are used.

  

% This should be set to \code{TRUE} for \code{link=}
% \code{\link{gordlink}},
% \code{\link{pordlink}},
% \code{\link{nbordlink}}.
% For these links the cutpoints must be an increasing sequence;
% if \code{reverse = FALSE} for then the cutpoints must be an
% decreasing sequence.


  
  }
  \item{ynames}{
  See \code{\link{multinomial}} for information.

  }
  \item{multiple.responses}{
  Logical.
  Multiple responses?
  If \code{TRUE} then the input should be
  a matrix with values \eqn{1,2,\dots,L},
  where \eqn{L=J+1} is the
  number of levels.
  Each column of the matrix is a response,
  i.e., multiple responses.
  A suitable matrix can be obtained from
  \code{Cut}.



  }
%  \item{apply.parint}{
%  Logical.
%  Whether the \code{parallel} argument should be applied to
%  the intercept term.  This should be set to \code{TRUE}
%  for \code{link=}
%  \code{\link{gordlink}},
%  \code{\link{pordlink}},
%  \code{\link{nbordlink}}.
%  See \code{\link{CommonVGAMffArguments}} for more information.
%
%
%  }
% \item{iscale}{
% Numeric. Initial values for the scale parameters.


% }


  \item{Thresh}{
  Character.
  The choices concern constraint matrices applied
  to the intercepts.
  They can be constrained to be
  equally-spaced (\emph{equid}istant)
  etc.
  See \code{\link{CommonVGAMffArguments}} and
  \code{\link{constraints}} for general
  information.
  Basically, the choice is pasted to the end of
  \code{"CM."} and that function is called.
  This means users can easily write their own
  \code{CM.}-type function.



%   The first choice is the default and causes the intercepts
%   to be estimated in an unconstrained (flexible) manner.
%   Actually, for this model, they will be sorted either in
%   ascending (default) or descending order, depending
%   on \code{reverse}.



% Orig.: "equidistant", "symmetric1",
%   "symmetric0", "qnorm"




If equally-spaced then the direction and the
reference level are controlled
by \code{Trev} and \code{Tref}, and
the constraint matrix
will be \eqn{M} by 2, with the second
column corresponding to the distance between
the thresholds.



If \code{"symm1"} then the fitted intercepts are
\emph{symmetric} about the median (\eqn{M} odd)
intercept.
If \eqn{M} is even then the median is the mean of
the two most inner and adjacent intercepts.
For this, \code{\link{CM.symm1}} is used to
construct the appropriate constraint matrix.



If \code{"symm0"} then the median intercept is 0
by definition and the symmetry occurs about the
origin. Thus the intercepts comprise pairs
that differ by sign only.
The appropriate constraint matrix is as with
\code{"symm1"} but with the first column deleted.
The choices \code{"symm1"} and \code{"symm0"}
are effectively equivalent to
\code{"symmetric"} and \code{"symmetric2"}
respectively in \pkg{ordinal}.




For \code{"qnorm"} then
\code{\link{CM.qnorm}} uses the
\code{qnorm((1:M)/(M+1))} quantiles
of the standard normal.



  }
  \item{Trev, Tref}{
  Support arguments for \code{Thresh}
  for equally-spaced intercepts.
  The logical argument \code{Trev} is
  applied first to give the direction
  (i.e., ascending or descending)
  before
  row \code{Tref} (ultimately numeric)
  of the first (intercept) constraint matrix
  is set to the reference level.
  See \code{\link{constraints}} for information.
    


  }

  \item{whitespace}{
  See \code{\link{CommonVGAMffArguments}} for
  information.


  }

}
\details{
  In this help file the response \eqn{Y} is assumed
  to be a factor with ordered values \eqn{1,2,\dots,J+1}.
  Hence \eqn{M} is the number of linear/additive
  predictors \eqn{\eta_j}{eta_j};
  for \code{cumulative()} one has \eqn{M=J}.
% and for \code{scumulative()} \eqn{M=2J}.



  This \pkg{VGAM} family function fits the class of
  \emph{cumulative link models} to (hopefully)
  an ordinal response.
  By default, the \emph{non-parallel} cumulative logit model
  is fitted, i.e.,
  \deqn{\eta_j = logit(P[Y \leq j])}{%
         eta_j = logit(P[Y<=j])}
  where \eqn{j=1,2,\dots,M} and
  the \eqn{\eta_j}{eta_j} are not constrained to be parallel.
  This is also known as the \emph{non-proportional odds model}.
  If the logit link is replaced by a complementary log-log link
  (\code{\link{clogloglink}}) then
  this is known as the \emph{proportional-hazards model}.



  In almost all the literature, the constraint matrices
  associated with this family of models are known.
  For example, setting
  \code{parallel = TRUE} will make all constraint matrices
  (except for the intercept) equal to a vector of \eqn{M} 1's.
  If the constraint matrices are equal, unknown and to
  be estimated,
  then this can be achieved by fitting the model as a
  reduced-rank vector generalized
  linear model (RR-VGLM; see \code{\link{rrvglm}}).
  Currently, reduced-rank vector generalized additive models
  (RR-VGAMs) have not been implemented here.



% The scaled version of \code{cumulative()},
% called \code{scumulative()},
% has \eqn{J} positive scaling factors.
% They are described in pages 154 and 177 of
% McCullagh and Nelder (1989);
% see their equation (5.4) in particular,
% which they call the \emph{generalized rational model}.



}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

Agresti, A. (2013).
\emph{Categorical Data Analysis},
3rd ed. Hoboken, NJ, USA: Wiley.


Agresti, A. (2010).
\emph{Analysis of Ordinal Categorical Data},
2nd ed. Hoboken, NJ, USA: Wiley.




%Boersch-Supan, P. H. (2021).
%Modeling insect phenology using ordinal
%regression and continuation ratio models.
%\emph{ReScience C},
%\bold{7.1}, 1--14.
%\doi{10.18637/jss.v032.i10}.

%\bold{7.1}(#5), 1--14.



%Dobson, A. J. and Barnett, A. (2008).
%\emph{An Introduction to Generalized Linear Models},
%3rd ed. Boca Raton, FL, USA: Chapman & Hall/CRC Press.


McCullagh, P. and Nelder, J. A. (1989).
\emph{Generalized Linear Models}, 2nd ed.
London: Chapman & Hall.


%Simonoff, J. S. (2003).
%\emph{Analyzing Categorical Data},
%New York: Springer-Verlag.


Tutz, G. (2012).
\emph{Regression for Categorical Data},
Cambridge: Cambridge University Press.



Tutz, G. and Berger, M. (2022).
Sparser ordinal regression models based
on parametric and additive location-shift
approaches.
\emph{International Statistical Review},
\bold{90}, 306--327.
\doi{10.1111/insr.12484}.



Yee, T. W. (2010).
The \pkg{VGAM} package for categorical data analysis.
\emph{Journal of Statistical Software},
\bold{32}, 1--34.
\doi{10.18637/jss.v032.i10}.



% \url{https://www.jstatsoft.org/article/view/v032i10/}.
% \url{https://www.jstatsoft.org/v32/i10/}.


Yee, T. W. and Wild, C. J. (1996).
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society,
  Series B, Methodological},
\bold{58}, 481--493.


%Further information and examples on categorical data analysis
%by the \pkg{VGAM} package can be found at
%\url{www.stat.auckland.ac.nz/~yee/VGAM/doc/categorical.pdf}.


}
\author{ Thomas W. Yee }
\note{
  The response should be either a matrix of counts
  (with row sums that
  are all positive), or a factor. In both cases,
  the \code{y} slot
  returned by
  \code{\link{vglm}}/\code{\link{vgam}}/\code{\link{rrvglm}}
  is the matrix
  of counts.
  The formula must contain an intercept term.
  Other \pkg{VGAM} family functions for an ordinal response
  include
  \code{\link{acat}},
  \code{\link{cratio}},
  \code{\link{sratio}}.
  For a nominal (unordered) factor response, the multinomial
  logit model (\code{\link{multinomial}}) is more appropriate.

  

  With the logit link, setting \code{parallel =
  TRUE} will fit a proportional odds model. Note
  that the \code{TRUE} here does not apply to
  the intercept term.  In practice, the validity
  of the \emph{proportional odds assumption}
  needs to be checked, e.g., by a likelihood
  ratio test (LRT).  If acceptable on the data,
  then numerical problems are less likely
  to occur during the fitting, and there are
  less parameters. Numerical problems occur
  when the linear/additive predictors cross,
  which results in probabilities outside of
  \eqn{(0,1)}; setting \code{parallel = TRUE}
  will help avoid this problem.

  

  Here is an example of the usage of the \code{parallel}
  argument.
  If there are covariates \code{x2}, \code{x3} and
  \code{x4}, then
  \code{parallel = TRUE ~ x2 + x3 -1} and
  \code{parallel = FALSE ~ x4} are equivalent.
  This would constrain the regression coefficients
  for \code{x2} and \code{x3} to be equal;
  those of the intercepts and \code{x4} would be different.

  

  If the data is inputted in \emph{long} format
  (not \emph{wide} format, as in \code{\link{pneumo}} below)
  and the self-starting initial values are not
  good enough then try using
  \code{mustart},
  \code{coefstart} and/or
  \code{etatstart}.
  See the example below.


  
  To fit the proportional odds model one can use the
  \pkg{VGAM} family function \code{\link{propodds}}.
  Note that \code{propodds(reverse)} is equivalent to
  \code{cumulative(parallel = TRUE, reverse = reverse)}
  (which is equivalent to
  \code{cumulative(parallel =}
  \code{TRUE, reverse = reverse, link = "logitlink")}).
  It is for convenience only. A call to
  \code{cumulative()} is preferred since it reminds the user
  that a parallelism assumption is made, as well as
  being a lot more flexible.


  
% In the future, this family function may be renamed to
% ``\code{cups}''
% (for \bold{cu}mulative \bold{p}robabilitie\bold{s})
% or ``\code{cute}''
% (for \bold{cu}mulative probabili\bold{t}i\bold{e}s).
% Please let me know if you strongly agree or disagree
% about this.



Category specific effects may be modelled using
the \code{xij}-facility; see
  \code{\link{vglm.control}} and \code{\link{fill1}}.




With most \code{Thresh}old choices,
the first few fitted regression coefficients
need care in their interpretation.  For example,
some values could be the distance away from
the median intercept.  Typing something
like \code{constraints(fit)[[1]]} gives the
constraint matrix of the intercept term.




}
\section{Warning }{
  No check is made to verify that the response is ordinal
  if the response is a matrix;
  see \code{\link[base:factor]{ordered}}.



  Boersch-Supan (2021) looks at sparse data and
  the numerical problems that result;
  see \code{\link{sratio}}.




}

\seealso{
  \code{\link{propodds}},
  \code{\link{constraints}},
  \code{\link{CM.ones}},
  \code{\link{CM.equid}},
  \code{\link{R2latvar}},
  \code{\link{ordsup}},
  \code{\link{prplot}},
  \code{\link{margeff}},
  \code{\link{acat}},
  \code{\link{cratio}},
  \code{\link{sratio}},
  \code{\link{multinomial}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{pneumo}},
  \code{\link{budworm}},
  \code{\link{Links}},
  \code{\link{hdeff.vglm}},
  \code{\link{logitlink}},
  \code{\link{probitlink}},
  \code{\link{clogloglink}},
  \code{\link{cauchitlink}},
  \code{\link{logistic1}}.


% \code{\link{pordlink}},
% \code{\link{gordlink}},
% \code{\link{nbordlink}},





}
\examples{
# Proportional odds model (p.179) of McCullagh and Nelder (1989)
pneumo <- transform(pneumo, let = log(exposure.time))
(fit <- vglm(cbind(normal, mild, severe) ~ let,
             cumulative(parallel = TRUE, reverse = TRUE), pneumo))
depvar(fit)  # Sample proportions (good technique)
fit@y        # Sample proportions (bad technique)
weights(fit, type = "prior")  # Number of observations
coef(fit, matrix = TRUE)
constraints(fit)  # Constraint matrices
apply(fitted(fit), 1, which.max)  # Classification
apply(predict(fit, newdata = pneumo, type = "response"),
      1, which.max)  # Classification
R2latvar(fit)

# Check that the model is linear in let ----------------------
fit2 <- vgam(cbind(normal, mild, severe) ~ s(let, df = 2),
             cumulative(reverse = TRUE), data = pneumo)
\dontrun{
 plot(fit2, se = TRUE, overlay = TRUE, lcol = 1:2, scol = 1:2) }

# Check the proportional odds assumption with a LRT ----------
(fit3 <- vglm(cbind(normal, mild, severe) ~ let,
              cumulative(parallel = FALSE, reverse = TRUE), pneumo))
pchisq(2 * (logLik(fit3) - logLik(fit)), df =
       length(coef(fit3)) - length(coef(fit)), lower.tail = FALSE)
lrtest(fit3, fit)  # More elegant

# A factor() version of fit ----------------------------------
# This is in long format (cf. wide format above)
Nobs <- round(depvar(fit) * c(weights(fit, type = "prior")))
sumNobs <- colSums(Nobs)  # apply(Nobs, 2, sum)

pneumo.long <-
  data.frame(symptoms = ordered(rep(rep(colnames(Nobs), nrow(Nobs)),
                                        times = c(t(Nobs))),
                                levels = colnames(Nobs)),
             let = rep(rep(with(pneumo, let), each = ncol(Nobs)),
                       times = c(t(Nobs))))
with(pneumo.long, table(let, symptoms))  # Should be same as pneumo


(fit.long1 <- vglm(symptoms ~ let, data = pneumo.long, trace = TRUE,
                   cumulative(parallel = TRUE, reverse = TRUE)))
coef(fit.long1, matrix = TRUE)  # cf. coef(fit, matrix = TRUE)
# Could try using mustart if fit.long1 failed to converge.
mymustart <- matrix(sumNobs / sum(sumNobs),
                    nrow(pneumo.long), ncol(Nobs), byrow = TRUE)
fit.long2 <- vglm(symptoms ~ let, mustart = mymustart,
                  cumulative(parallel = TRUE, reverse = TRUE),
                  data = pneumo.long, trace = TRUE)
coef(fit.long2, matrix = TRUE)  # cf. coef(fit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}

% pneumo$let <- log(pneumo$exposure.time)

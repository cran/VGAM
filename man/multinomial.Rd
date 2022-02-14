\name{multinomial}
\alias{multinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multinomial Logit Model }
\description{
  Fits a multinomial logit model (MLM) to a (preferably unordered)
  factor response.
% The deflated--altered MLM (DAML model, or DAMLM) may also be
% fitted.

  
}
\usage{
multinomial(zero = NULL, parallel = FALSE, nointercept = NULL,
            refLevel = "(Last)", imethod = 1, imu = NULL,
            byrow.arg = FALSE, whitespace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
%           refLevel = "(Last)", d.mlm = NULL, imethod = 1,
\arguments{

  \item{zero}{
  Can be an integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Any values must be from the set \{1,2,\ldots,\eqn{M}\}.
  The default value means none are modelled as intercept-only terms.
  See \code{\link{CommonVGAMffArguments}} for more information.



  }
  \item{parallel}{
  A logical, or formula specifying which terms have
  equal/unequal coefficients.


  }
  \item{nointercept, whitespace}{
  See \code{\link{CommonVGAMffArguments}} for details.


  }
  \item{imu, byrow.arg}{
  See \code{\link{CommonVGAMffArguments}} for details.


  }
  \item{refLevel}{
  Either a (1) single positive integer or (2) a value of the factor
  or (3) a character string.
  If inputted as an integer then it specifies which
  column of the response matrix is the reference or baseline level.
  The default is the \emph{last} one (the \eqn{(M+1)}th one).
  If used, this argument will be usually assigned the value \code{1}.
  If inputted as a value of a factor then beware of missing values
  of certain levels of the factor
  (\code{drop.unused.levels = TRUE} or
   \code{drop.unused.levels = FALSE}).
  See the example below.
  If inputted as a character string then this should be
  equal to (A) one of the levels of the factor response, else (B) one
  of the column names of the matrix response of counts;
  e.g.,
  \code{vglm(cbind(normal, mild, severe) ~ let,}
  \code{multinomial(refLevel = "severe"), data = pneumo)}
  if it was (incorrectly because the response is ordinal)
  applied to the \code{\link{pneumo}} data set.
  Another example is
  \code{vglm(ethnicity ~ age,}
  \code{multinomial(refLevel = "European"), data = xs.nz)}
  if it was applied to the \code{\link[VGAMdata]{xs.nz}} data set.



  }
  \item{imethod}{
    Choosing 2 will use the mean sample proportions of each
    column of the response matrix, which corresponds to
    the MLEs for intercept-only models.
    See \code{\link{CommonVGAMffArguments}} for more details.



%    Initial values for the DAMLM often need to be better than
%    the ordinary MLM, hence choosing various values of
%    this argument is encouraged.



  }
% \item{d.mlm}{
%   This argument helps implement GAITD regression by fitting
%   the DAMLM.
%   The argument is the same as \code{\link{gaitdpoisson}}
%   and is the set of deflated values (levels)
%   from \{1,2,\ldots,\eqn{M+1}\}.
%   The ordinary MLM may be called the AMLM because the
%   probabilities are modelled as \emph{altered} values.

    


%   The second argument may be assigned higher values if
%   necessary in order to ensure that the probabilities
%   are positive---the vector is recycled to \code{length(d.mlm)}.
%   If the baseline group is chosen to be the one with the
%   largest  probability or probabilities
%   then \code{eta.d.max = 0} should suffice,
%   however, it might be necessary to increase its value
%   slightly, e.g., by 1 or 2 maybe.



% }
}
\details{
  In this help file the response \eqn{Y} is assumed to be
  a factor with unordered values \eqn{1,2,\dots,M+1}, so
  that \eqn{M} is the number of linear/additive predictors
  \eqn{\eta_j}{eta_j}.


  
  The default model can be written
  \deqn{\eta_j = \log(P[Y=j]/ P[Y=M+1])}{%
         eta_j =  log(P[Y=j]/ P[Y=M+1])}
  where \eqn{\eta_j}{eta_j} is the \eqn{j}th
  linear/additive predictor.
  Here, \eqn{j=1,\ldots,M}, and \eqn{\eta_{M+1}}{eta_{M+1}}
  is 0 by definition. That is, the last level of the factor,
  or last column of the response matrix, is taken as the
  reference level or baseline---this is for identifiability
  of the parameters. The reference or baseline level can
  be changed with the \code{refLevel} argument.


  
  In almost all the literature, the constraint matrices
  associated with this family of models are known.
  For example, setting \code{parallel = TRUE} will make all
  constraint matrices
  (including the intercept)
  equal to a vector of \eqn{M} 1's;
  to suppress the intercepts from being parallel then
  set \code{parallel = FALSE ~ 1}.
  If the constraint matrices are
  unknown and to be estimated, then this can be achieved by
  fitting the model as a reduced-rank vector generalized
  linear model (RR-VGLM; see \code{\link{rrvglm}}).
  In particular, a multinomial logit model with unknown
  constraint matrices is known as a \emph{stereotype}
  model (Anderson, 1984), and can be fitted with
  \code{\link{rrvglm}}.



% Pre 20170816; Stuart Coles picked up an error:
% For example, setting \code{parallel = TRUE} will make all
% constraint matrices
% (except for the intercept)
% equal to a vector of \eqn{M} 1's.



The above details correspond to the ordinary MLM where
all the levels are \emph{altered} (in the terminology
of GAITD regression).


%The \eqn{j}th linear/additive predictor is different
%from above if level \eqn{j} is deflated: it is negated
%so that increasing the linear/additive predictor results
%in a decreasing probability.
%Thus the DAMLM simply has a sign change for those levels
%that are deflated.





%DAMLMs require skill to fit because handling
%data that is extreme (high or low fitted values)
%means that \code{eta.d.max} may need to be increased.
%It pays to fit an intercept-only model first and to
%monitor convergence carefully.
%One disadvantage of DAMLMs is that the regression
%coefficients depend on \code{eta.d.max} so there is
%no unique solution.




}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.


}
\references{



Yee, T. W. (2010).
The \pkg{VGAM} package for categorical data analysis.
\emph{Journal of Statistical Software},
\bold{32}, 1--34.
\doi{10.18637/jss.v032.i10}.



% \url{https://www.jstatsoft.org/article/view/v032i10/}.
% \url{https://www.jstatsoft.org/v32/i10/}.



Yee, T. W. and Hastie, T. J. (2003).
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.



Yee, T. W. and Ma, C. (2021).
Generally--altered, --inflated, --truncated and --altered
count regression.
\emph{In preparation}.



McCullagh, P. and Nelder, J. A. (1989).
\emph{Generalized Linear Models}, 2nd ed.
London: Chapman & Hall.


Agresti, A. (2013).
\emph{Categorical Data Analysis},
3rd ed. Hoboken, NJ, USA: Wiley.


Hastie, T. J., Tibshirani, R. J. and Friedman, J. H. (2009).
\emph{The Elements of Statistical Learning: Data Mining,
      Inference and Prediction},
2nd ed.
New York, USA: Springer-Verlag.


Simonoff, J. S. (2003).
\emph{Analyzing Categorical Data},
New York, USA: Springer-Verlag.


Anderson, J. A. (1984).
Regression and ordered categorical variables.
\emph{Journal of the Royal Statistical Society, Series B,
  Methodological},
\bold{46}, 1--30.


Tutz, G. (2012).
\emph{Regression for Categorical Data},
Cambridge University Press.



%Further information and examples on categorical data analysis
%by the \pkg{VGAM} package can be found at
%\url{http://www.stat.auckland.ac.nz/~yee/VGAM/doc/categorical.pdf}.


}
\author{ Thomas W. Yee }

\note{
  The response should be either a matrix of counts
  (with row sums that are all positive), or a
  factor. In both cases, the \code{y} slot returned by
  \code{\link{vglm}}/\code{\link{vgam}}/\code{\link{rrvglm}}
  is the matrix of sample proportions.


  The multinomial logit model is more appropriate for a nominal
  (unordered) factor response than for an ordinal (ordered) factor
  response.
  Models more suited for the latter include those based on
  cumulative probabilities, e.g., \code{\link{cumulative}}.


  \code{multinomial} is prone to numerical difficulties if
  the groups are separable and/or the fitted probabilities
  are close to 0 or 1. The fitted values returned
  are estimates of the probabilities \eqn{P[Y=j]} for
  \eqn{j=1,\ldots,M+1}. See \pkg{safeBinaryRegression}
  for the logistic regression case.


  Here is an example of the usage of the \code{parallel}
  argument. If there are covariates \code{x2}, \code{x3}
  and \code{x4}, then \code{parallel = TRUE ~ x2 + x3 -
  1} and \code{parallel = FALSE ~ x4} are equivalent. This
  would constrain the regression coefficients for \code{x2}
  and \code{x3} to be equal; those of the intercepts and
  \code{x4} would be different.


  In Example 4 below, a conditional logit model is
  fitted to an artificial data set that explores how
  cost and travel time affect people's decision about
  how to travel to work. Walking is the baseline group.
  The variable \code{Cost.car} is the difference between
  the cost of travel to work by car and walking, etc. The
  variable \code{Time.car} is the difference between
  the travel duration/time to work by car and walking,
  etc. For other details about the \code{xij} argument see
  \code{\link{vglm.control}} and \code{\link{fill}}.


  The \code{\link[nnet]{multinom}} function in the
  \pkg{nnet} package uses the first level of the factor as
  baseline, whereas the last level of the factor is used
  here. Consequently the estimated regression coefficients
  differ.


}

% In the future, this family function may be renamed to
% ``\code{mum}'' (for \bold{mu}ltinomial logit \bold{m}odel).
% Please let me know if you strongly agree or disagree about this.

\section{Warning }{
  No check is made to verify that the response is nominal.


% Using \code{d.mlm} may not work with other options,
% e.g., the Hauck-Donner effect detection using derivatives.


  See \code{\link{CommonVGAMffArguments}} for more warnings.


}
\seealso{
    \code{\link{multilogitlink}},
    \code{\link{margeff}},
    \code{\link{cumulative}},
    \code{\link{acat}},
    \code{\link{cratio}},
    \code{\link{sratio}},
    \code{\link{dirichlet}},
    \code{\link{dirmultinomial}},
    \code{\link{rrvglm}},
    \code{\link{fill1}},
    \code{\link[stats:Multinom]{Multinomial}},
    \code{\link{gaitdpoisson}},
    \code{\link{Gaitdpois}},
    \code{\link[datasets]{iris}}.



%   \code{\link{gatnbinomial.mlm}},
% The author's homepage has further documentation about
% categorical data analysis using \pkg{VGAM}.


}
%   \code{\link[base:Multinom]{rmultinom}}
%   \code{\link{pneumo}},


\examples{
# Example 1: fit a multinomial logit model to Edgar Anderson's iris data
data(iris)
\dontrun{ fit <- vglm(Species ~ ., multinomial, iris)
coef(fit, matrix = TRUE) }

# Example 2a: a simple example
ycounts <- t(rmultinom(10, size = 20, prob = c(0.1, 0.2, 0.8)))
fit <- vglm(ycounts ~ 1, multinomial)
head(fitted(fit))   # Proportions
fit@prior.weights   # NOT recommended for extraction of prior weights
weights(fit, type = "prior", matrix = FALSE)  # The better method
depvar(fit)         # Sample proportions; same as fit@y
constraints(fit)    # Constraint matrices

# Example 2b: Different reference level used as the baseline
fit2 <- vglm(ycounts ~ 1, multinomial(refLevel = 2))
coef(fit2, matrix = TRUE)
coef(fit , matrix = TRUE)  # Easy to reconcile this output with fit2

# Example 3: The response is a factor.
nn <- 10
dframe3 <- data.frame(yfac = gl(3, nn, labels = c("Ctrl", "Trt1", "Trt2")),
                      x2   = runif(3 * nn))
myrefLevel <- with(dframe3, yfac[12])
fit3a <- vglm(yfac ~ x2, multinomial(refLevel = myrefLevel), dframe3)
fit3b <- vglm(yfac ~ x2, multinomial(refLevel = 2), dframe3)
coef(fit3a, matrix = TRUE)  # "Trt1" is the reference level
coef(fit3b, matrix = TRUE)  # "Trt1" is the reference level
margeff(fit3b)

# Example 4: Fit a rank-1 stereotype model
fit4 <- rrvglm(Country ~ Width + Height + HP, multinomial, data = car.all)
coef(fit4)  # Contains the C matrix
constraints(fit4)$HP       # The A matrix
coef(fit4, matrix = TRUE)  # The B matrix
Coef(fit4)@C               # The C matrix
concoef(fit4)              # Better to get the C matrix this way
Coef(fit4)@A               # The A matrix
svd(coef(fit4, matrix = TRUE)[-1, ])$d  # This has rank 1; = C %*% t(A)
# Classification (but watch out for NAs in some of the variables):
apply(fitted(fit4), 1, which.max)  # Classification
colnames(fitted(fit4))[apply(fitted(fit4), 1, which.max)]  # Classification
apply(predict(fit4, car.all, type = "response"), 1, which.max)  # Ditto


# Example 5: The use of the xij argument (aka conditional logit model)
set.seed(111)
nn <- 100  # Number of people who travel to work
M <- 3  # There are M+1 models of transport to go to work
ycounts <- matrix(0, nn, M+1)
ycounts[cbind(1:nn, sample(x = M+1, size = nn, replace = TRUE))] = 1
dimnames(ycounts) <- list(NULL, c("bus","train","car","walk"))
gotowork <- data.frame(cost.bus  = runif(nn), time.bus  = runif(nn),
                       cost.train= runif(nn), time.train= runif(nn),
                       cost.car  = runif(nn), time.car  = runif(nn),
                       cost.walk = runif(nn), time.walk = runif(nn))
gotowork <- round(gotowork, digits = 2)  # For convenience
gotowork <- transform(gotowork,
                      Cost.bus   = cost.bus   - cost.walk,
                      Cost.car   = cost.car   - cost.walk,
                      Cost.train = cost.train - cost.walk,
                      Cost       = cost.train - cost.walk,  # for labelling
                      Time.bus   = time.bus   - time.walk,
                      Time.car   = time.car   - time.walk,
                      Time.train = time.train - time.walk,
                      Time       = time.train - time.walk)  # for labelling
fit <- vglm(ycounts ~ Cost + Time,
            multinomial(parall = TRUE ~ Cost + Time - 1),
            xij = list(Cost ~ Cost.bus + Cost.train + Cost.car,
                       Time ~ Time.bus + Time.train + Time.car),
            form2 =  ~ Cost + Cost.bus + Cost.train + Cost.car +
                       Time + Time.bus + Time.train + Time.car,
            data = gotowork, trace = TRUE)
head(model.matrix(fit, type = "lm"))   # LM model matrix
head(model.matrix(fit, type = "vlm"))  # Big VLM model matrix
coef(fit)
coef(fit, matrix = TRUE)
constraints(fit)
summary(fit)
max(abs(predict(fit) - predict(fit, new = gotowork)))  # Should be 0
}
\keyword{models}
\keyword{regression}



% 20100915; this no longer works:
% # Example 2c: Different input to Example 2a but same result
% w <- apply(ycounts, 1, sum)  # Prior weights
% yprop <- ycounts / w    # Sample proportions
% fitprop <- vglm(yprop ~ 1, multinomial, weights=w)
% head(fitted(fitprop))   # Proportions
% weights(fitprop, type="prior", matrix=FALSE)
% fitprop@y # Same as the input




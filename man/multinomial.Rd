\name{multinomial}
\alias{multinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multinomial Logit Model }
\description{
  Fits a multinomial logit model to an unordered factor response.
}
\usage{
multinomial(zero = NULL, parallel = FALSE, nointercept = NULL,
            refLevel = "last")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  In the following, the response \eqn{Y} is assumed to be a factor with
  unordered values \eqn{1,2,\dots,M+1}, so that \eqn{M} is the number
  of linear/additive predictors \eqn{\eta_j}{eta_j}.

  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2,\ldots,\eqn{M}\}.
  The default value means none are modelled as intercept-only terms.

  }
  \item{parallel}{
  A logical, or formula specifying which terms have
  equal/unequal coefficients.

  }
  \item{nointercept}{
  An integer-valued vector specifying which
  linear/additive predictors have no intercepts.
  The values must be from the set \{1,2,\ldots,\eqn{M}\}.

  }
  \item{refLevel}{
  Either a single positive integer or a value of a factor.
  If an integer then it specifies which
  column of the response matrix is the reference or baseline level.
  The default is the last one (the \eqn{(M+1)}th one).
  If used, this argument will be often assigned the value \code{1}.
  If inputted as a value of a factor then beware of missing values
  of certain levels of the factor
  (\code{drop.unused.levels=TRUE} or
   \code{drop.unused.levels=FALSE}).
  See the example below.

  }
}
\details{
  The default model can be written
  \deqn{\eta_j = \log(P[Y=j]/ P[Y=M+1])}{%
         eta_j =  log(P[Y=j]/ P[Y=M+1])}
  where \eqn{\eta_j}{eta_j} is the \eqn{j}th linear/additive predictor.
  Here, \eqn{j=1,\ldots,M} and \eqn{\eta_{M+1}}{eta_{M+1}} is 0 by
  definition.  That is, the last level of the factor, or last column of
  the response matrix, is taken as the reference level or baseline---this
  is for identifiability of the parameters.
  The reference or baseline level can be chosen by using the
  \code{refLevel} argument.

  In almost all the literature, the constraint matrices associated
  with this family of models are known. For example, setting
  \code{parallel=TRUE} will make all constraint matrices (except for
  the intercept) equal to a vector of \eqn{M} 1's.  If the constraint
  matrices are unknown and to be estimated, then this can be achieved
  by fitting the model as a reduced-rank vector generalized linear model
  (RR-VGLM; see \code{\link{rrvglm}}).  In particular, a multinomial logit
  model with unknown constraint matrices is known as a stereotype model
  (Anderson, 1984), and can be fitted with \code{\link{rrvglm}}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.

}
\references{
Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling}, 
\bold{3}, 15--41.

McCullagh, P. and Nelder, J. A. (1989)
\emph{Generalized Linear Models}, 2nd ed.
London: Chapman & Hall.

Agresti, A. (2002)
\emph{Categorical Data Analysis},
2nd ed. New York: Wiley.

Simonoff, J. S. (2003)
\emph{Analyzing Categorical Data},
New York: Springer-Verlag.

Anderson, J. A. (1984)
Regression and ordered categorical variables. 
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{46}, 1--30.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.

}
\author{ Thomas W. Yee }

\note{
  The response should be either a matrix of counts (with row sums that are
  all positive), or a factor. In both cases, the \code{y} slot returned
  by \code{\link{vglm}}/\code{\link{vgam}}/\code{\link{rrvglm}} is the
  matrix of counts.

  The multinomial logit model is more appropriate for a nominal
  (unordered) factor response. For an ordinal (ordered) factor
  response, models such as those based on cumulative probabilities
  (see \code{\link{cumulative}}) are more suited.

  \code{multinomial} is prone to numerical difficulties if the groups
  are separable and/or the fitted probabilities are close to 0 or 1.
  The fitted values returned are estimates of the probabilities
  \eqn{P[Y=j]} for \eqn{j=1,\ldots,M+1}.

  Here is an example of the usage of the \code{parallel} argument.
  If there are covariates \code{x1}, \code{x2} and \code{x3}, then
  \code{parallel = TRUE ~ x1 + x2 -1} and
  \code{parallel = FALSE ~ x3} are equivalent. This would constrain
  the regression coefficients for \code{x1} and \code{x2} to be
  equal; those of the intercepts and \code{x3} would be different.

  In Example 4 below, a conditional logit model is fitted to a artificial
  data set that explores how cost and travel time affect people's
  decision about how to travel to work.  Walking is the baseline group.
  The variable \code{Cost.car} is the difference between the cost of
  travel to work by car and walking, etc.  The variable \code{Durn.car}
  is the difference between the travel duration/time to work by car and
  walking, etc.  For other details about the \code{xij} argument see
  \code{\link{vglm.control}} and \code{\link{fill}}.

  The \code{\link[nnet]{multinom}} function in the \pkg{nnet} package
  uses the first level of the factor as baseline, whereas the last
  level of the factor is used here. Consequently the estimated 
  regression coefficients differ.

}

% In the future, this family function may be renamed to
% ``\code{mum}'' (for \bold{mu}ltinomial logit \bold{m}odel).
% Please let me know if you strongly agree or disagree about this.

\section{Warning }{
  The arguments \code{zero} and \code{nointercept} can be inputted
  with values that fail. For example, \code{multinomial(zero=2,
  nointercept=1:3)} means the second linear/additive predictor is
  identically zero, which will cause a failure.

  Be careful about the use of other potentially contradictory constraints,
  e.g., \code{multinomial(zero=2, parallel = TRUE ~ x3)}.  If in doubt,
  apply \code{constraints()} to the fitted object to check.

  No check is made to verify that the response is nominal.

}
\seealso{
    \code{\link{acat}},
    \code{\link{cumulative}},
    \code{\link{cratio}},
    \code{\link{sratio}},
    \code{\link{dirichlet}},
    \code{\link{dirmultinomial}},
    \code{\link{rrvglm}},
    \code{\link[stats:Multinom]{Multinomial}},
    \code{\link[base]{iris}}.

}
%   \code{\link[base:Multinom]{rmultinom}}

\examples{
# Example 1: fit a multinomial logit model to Edgar Anderson's iris data
data(iris)
\dontrun{
fit = vglm(Species ~ ., multinomial, iris)
coef(fit, matrix=TRUE) 
}


# Example 2a: a simple example 
y = t(rmultinom(10, size = 20, prob=c(0.1,0.2,0.8))) # Counts
fit = vglm(y ~ 1, multinomial)
fitted(fit)[1:4,]   # Proportions
fit@prior.weights # Not recommended for extraction of prior weights
weights(fit, type="prior", matrix=FALSE) # The better method
fit@y   # Sample proportions
constraints(fit)   # Constraint matrices

# Example 2b: Different reference level used as the baseline 
fit2 = vglm(y ~ 1, multinomial(refLevel=2))
coef(fit2, matrix=TRUE)
coef(fit , matrix=TRUE) # Easy to reconcile this output with fit2

# Example 2c: Different input to Example 2a but same result
w = apply(y, 1, sum) # Prior weights
yprop = y / w    # Sample proportions
fitprop = vglm(yprop ~ 1, multinomial, weights=w)
fitted(fitprop)[1:4,]   # Proportions
weights(fitprop, type="prior", matrix=FALSE)
fitprop@y # Same as the input


# Example 3: The response is a factor.
nn = 10
yfactor = gl(3, nn, labels = c("Control", "Treatment1", "Treatment2"))
x = runif(3 * nn)
fit3a = vglm(yfactor ~ x, multinomial(refLevel=yfactor[12]))
fit3b = vglm(yfactor ~ x, multinomial(refLevel=2))
coef(fit3a, matrix=TRUE)  # "Treatment1" is the reference level
coef(fit3b, matrix=TRUE)  # "Treatment1" is the reference level


# Example 4: Fit a rank-1 stereotype model 
data(car.all)
fit4 = rrvglm(Country ~ Width + Height + HP, multinomial, car.all, Rank=1)
coef(fit4)   # Contains the C matrix
constraints(fit4)$HP     # The A matrix 
coef(fit4, matrix=TRUE)  # The B matrix
Coef(fit4)@C             # The C matrix 
ccoef(fit4)              # Better to get the C matrix this way
Coef(fit4)@A             # The A matrix 
svd(coef(fit4, matrix=TRUE)[-1,])$d    # This has rank 1; = C %*% t(A) 


# Example 5: The use of the xij argument (conditional logit model)
set.seed(111)
n = 100  # Number of people who travel to work
M = 3  # There are M+1 models of transport
ymat = matrix(0, n, M+1)
ymat[cbind(1:n, sample(x=M+1, size=n, replace=TRUE))] = 1
dimnames(ymat) = list(NULL, c("bus","train","car","walk"))
transport = data.frame(cost.bus=runif(n), cost.train=runif(n),
                       cost.car=runif(n), cost.walk=runif(n),
                       durn.bus=runif(n), durn.train=runif(n),
                       durn.car=runif(n), durn.walk=runif(n))
transport = round(transport, dig=2) # For convenience
transport = transform(transport,
                      Cost.bus   = cost.bus   - cost.walk,
                      Cost.car   = cost.car   - cost.walk,
                      Cost.train = cost.train - cost.walk,
                      Durn.bus   = durn.bus   - durn.walk,
                      Durn.car   = durn.car   - durn.walk,
                      Durn.train = durn.train - durn.walk)
fit = vglm(ymat ~ Cost.bus + Cost.train + Cost.car + 
                  Durn.bus + Durn.train + Durn.car,
           fam = multinomial,
           xij = list(Cost ~ Cost.bus + Cost.train + Cost.car,
                      Durn ~ Durn.bus + Durn.train + Durn.car),
           data=transport)
model.matrix(fit, type="lm")[1:7,]   # LM model matrix
model.matrix(fit, type="vlm")[1:7,]  # Big VLM model matrix
coef(fit)
coef(fit, matrix=TRUE)
coef(fit, matrix=TRUE, compress=FALSE)
summary(fit)
}
\keyword{models}
\keyword{regression}


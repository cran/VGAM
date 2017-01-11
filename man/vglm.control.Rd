\name{vglm.control}
\alias{vglm.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Control Function for vglm() }
\description{
  Algorithmic constants and parameters for running \code{vglm} are set
  using this function.

}
\usage{
vglm.control(checkwz = TRUE, Check.rank = TRUE, Check.cm.rank = TRUE,
             criterion = names(.min.criterion.VGAM),
             epsilon = 1e-07, half.stepsizing = TRUE,
             maxit = 30, noWarning = FALSE,
             stepsize = 1, save.weights = FALSE,
             trace = FALSE, wzepsilon = .Machine$double.eps^0.75,
             xij = NULL, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{checkwz}{
  logical indicating whether the diagonal elements
  of the working weight matrices should be checked
  whether they are sufficiently positive, i.e., greater
  than \code{wzepsilon}. If not, any values less than
  \code{wzepsilon} are replaced with this value.


  }
  \item{Check.rank}{
  logical indicating whether the rank of the VLM matrix
  should be checked. If this is not of full column rank then
  the results are not to be trusted.
  The default is to give an error message if the VLM
  matrix is not of full column rank.


  }
  \item{Check.cm.rank}{
  logical indicating whether the rank of each constraint matrix
  should be checked. If this is not of full column rank then
  an error will occur. Under no circumstances should any
  constraint matrix have a rank less than the number of columns.



  }
  \item{criterion}{
  character variable describing what criterion is to be
  used to test for convergence. The possibilities are
  listed in \code{.min.criterion.VGAM}, but most family
  functions only implement a few of these.


  }
  \item{epsilon}{
  positive convergence tolerance epsilon. Roughly speaking,
  the Newton-Raphson/Fisher-scoring iterations are assumed
  to have converged when two successive \code{criterion}
  values are within \code{epsilon} of each other.


  }
  \item{half.stepsizing}{
  logical indicating if half-stepsizing is allowed. For
  example, in maximizing a log-likelihood, if the next
  iteration has a log-likelihood that is less than
  the current value of the log-likelihood, then a half
  step will be taken.  If the log-likelihood is still
  less than at the current position, a quarter-step
  will be taken etc. Eventually a step will be taken
  so that an improvement is made to the convergence
  criterion.  \code{half.stepsizing} is ignored if
  \code{criterion == "coefficients"}.


  }
  \item{maxit}{
  maximum number of (usually Fisher-scoring) iterations
  allowed.  Sometimes Newton-Raphson is used.


  }
  \item{noWarning}{
  logical indicating whether to suppress a warning if
  convergence is not obtained within \code{maxit} iterations.
  This is ignored if \code{maxit = 1} is set.


  }
  \item{stepsize}{
  usual step size to be taken between each
  Newton-Raphson/Fisher-scoring iteration. It should be a
  value between 0 and 1, where a value of unity corresponds
  to an ordinary step.  A value of 0.5 means half-steps are
  taken.  Setting a value near zero will cause convergence
  to be generally slow but may help increase the chances
  of successful convergence for some family functions.


  }
  \item{save.weights}{
  logical indicating whether the \code{weights} slot of a
  \code{"vglm"} object will be saved on the object. If not,
  it will be reconstructed when needed, e.g., \code{summary}.
  Some family functions have \code{save.weights = TRUE} and
  others have \code{save.weights = FALSE} in their control
  functions.


  }
  \item{trace}{
  logical indicating if output should be produced for each
  iteration.  Setting \code{trace = TRUE} is recommended in
  general because \pkg{VGAM} fits a very broad variety of
  models and distributions, and for some of them, convergence
  is intrinsically more difficult. Monitoring convergence
  can help check that the solution is reasonable or that
  a problem has occurred.  It may suggest better initial
  values are needed, the making of invalid assumptions,
  or that the model is inappropriate for the data, etc.


  }
  \item{wzepsilon}{
  small positive number used to test whether the diagonals
  of the working weight matrices are sufficiently positive.


  }
  \item{xij}{
  A formula or a list of formulas.
  Each formula has a RHS giving \eqn{M} terms making up a
  covariate-dependent term (whose name is the response).
  That is, it creates a variable
  that takes on different values for each linear/additive predictor,
  e.g., the ocular pressure of each eye.
  The \eqn{M} terms must be unique;
  use \code{\link{fill1}}, \code{fill2}, \code{fill3}, etc. if necessary.
  Each formula should have a response which is taken as the
  name of that variable, and the \eqn{M} terms are enumerated
  in sequential order.  Each of the \eqn{M} terms multiply
  each successive row of the constraint matrix.
  When \code{xij} is used, the use of \code{form2} is also required
  to give \emph{every} term used by the model.



  The function \code{\link{Select}} can be used to
  select variables beginning with the same character string.



  }
% \item{jix}{
% A formula or a list of formulas specifying
% which explanatory variables are to be plotted for each \code{xij} term.
% For example, in the code below,
% the term \code{BS(dumm)} could be plotted against either
% \code{dum1} or \code{dum2}, therefore
% either \code{jix=dum1} or \code{jix=dum2} are ok.
% This argument is made use of by \code{plotvgam()}.
% Each formula has a RHS giving \eqn{r_k} unique terms,
% one for each column of the constraint matrix.
% Each formula should have a response that matches the
% \code{formula} argument.
% The argument \code{jix} is a reversal of \code{xij} to emphasize
% the same framework for handling terms involving covariates that have
% different values for each linear/additive predictor.
%
% }

  \item{\dots}{
  other parameters that may be picked up from control
  functions that are specific to the \pkg{VGAM} family function.


  }
}
\details{
  Most of the control parameters are used within
  \code{vglm.fit} and you will have to look at that to
  understand the full details.


  Setting \code{save.weights = FALSE} is useful for some
  models because the \code{weights} slot of the object
  is the largest and so less memory is used to store the
  object. However, for some \pkg{VGAM} family function,
  it is necessary to set \code{save.weights = TRUE} because
  the \code{weights} slot cannot be reconstructed later.


  }
\value{
  A list with components matching the input names. A little error
  checking is done, but not much.
  The list is assigned to the \code{control} slot of
  \code{vglm} objects.


}
\references{

  Yee, T. W. and Hastie, T. J. (2003)
  Reduced-rank vector generalized linear models.
  \emph{Statistical Modelling},
  \bold{3}, 15--41.


}
\author{ Thomas W. Yee}
\note{
  Reiterating from above,
  setting \code{trace = TRUE} is recommended in general.


  In Example 2 below there are two covariates that have
  linear/additive predictor specific values.  These are
  handled using the \code{xij} argument.


}

\section{Warning}{
  For some applications the default convergence criterion should
  be tightened.
  Setting something like \code{criterion = "coef", epsilon = 1e-09}
  is one way to achieve this, and also add
  \code{trace = TRUE} to monitor the convergence.
  Setting  \code{maxit} to some higher number is usually not
  needed, and needing to do so suggests something is wrong, e.g.,
  an ill-conditioned model, over-fitting or under-fitting.


}



\seealso{
  \code{\link{vglm}},
  \code{\link{fill1}}.
  The author's homepage has further documentation about
  the \code{xij} argument;
  see also \code{\link{Select}}.


}

\examples{
# Example 1.
pneumo <- transform(pneumo, let = log(exposure.time))
vglm(cbind(normal, mild, severe) ~ let, multinomial, data = pneumo,
     crit = "coef", step = 0.5, trace = TRUE, epsil = 1e-8, maxit = 40)


# Example 2. The use of the xij argument (simple case).
ymat <- rdiric(n <- 1000, shape = rep(exp(2), len = 4))
mydat <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n),
                    x4 = runif(n),
                    z1 = runif(n), z2 = runif(n), z3 = runif(n),
                    z4 = runif(n))
mydat <- transform(mydat, X = x1, Z = z1)
mydat <- round(mydat, digits = 2)
fit2 <- vglm(ymat ~ X + Z,
             dirichlet(parallel = TRUE), data = mydat, trace = TRUE,
             xij = list(Z ~ z1 + z2 + z3 + z4,
                        X ~ x1 + x2 + x3 + x4),
             form2 = ~  Z + z1 + z2 + z3 + z4 +
                        X + x1 + x2 + x3 + x4)
head(model.matrix(fit2, type =  "lm"))  # LM model matrix
head(model.matrix(fit2, type = "vlm"))  # Big VLM model matrix
coef(fit2)
coef(fit2, matrix = TRUE)
max(abs(predict(fit2)-predict(fit2, new = mydat)))  # Predicts correctly
summary(fit2)
\dontrun{
# plotvgam(fit2, se = TRUE, xlab = "x1", which.term = 1)  # Bug!
# plotvgam(fit2, se = TRUE, xlab = "z1", which.term = 2)  # Bug!
plotvgam(fit2, xlab = "x1")  # Correct
plotvgam(fit2, xlab = "z1")  # Correct
}


# Example 3. The use of the xij argument (complex case).
set.seed(123)
coalminers <- transform(coalminers,
                        Age = (age - 42) / 5,
                        dum1 = round(runif(nrow(coalminers)), digits = 2),
                        dum2 = round(runif(nrow(coalminers)), digits = 2),
                        dum3 = round(runif(nrow(coalminers)), digits = 2),
                        dumm = round(runif(nrow(coalminers)), digits = 2))
BS <- function(x, ..., df = 3)
  sm.bs(c(x,...), df = df)[1:length(x),,drop = FALSE]
NS <- function(x, ..., df = 3)
  sm.ns(c(x,...), df = df)[1:length(x),,drop = FALSE]

# Equivalently...
BS <- function(x, ..., df = 3)
  head(sm.bs(c(x,...), df = df), length(x), drop = FALSE)
NS <- function(x, ..., df = 3)
  head(sm.ns(c(x,...), df = df), length(x), drop = FALSE)

fit3 <- vglm(cbind(nBnW,nBW,BnW,BW) ~ Age + NS(dum1, dum2),
             fam = binom2.or(exchangeable = TRUE, zero = 3),
             xij = list(NS(dum1, dum2) ~ NS(dum1, dum2) +
                                         NS(dum2, dum1) +
                                         fill(NS( dum1))),
             form2 = ~  NS(dum1, dum2) + NS(dum2, dum1) + fill(NS(dum1)) +
                        dum1 + dum2 + dum3 + Age + age + dumm,
             data = coalminers, trace = TRUE)
head(model.matrix(fit3, type = "lm"))   # LM model matrix
head(model.matrix(fit3, type = "vlm"))  # Big VLM model matrix
coef(fit3)
coef(fit3, matrix = TRUE)
\dontrun{
plotvgam(fit3, se = TRUE, lcol = "red", scol = "blue", xlab = "dum1")
}
}
\keyword{models}
\keyword{regression}





% zz 20090506 put elsewhere:
%
%
%# Example 4. The use of the xij argument (complex case).
%# Here is one method to handle the xij argument with a term that
%# produces more than one column in the model matrix.
%# The constraint matrix for 'op' has one column.
%POLY3 <- function(x, ...) {
%    # A cubic; ensures that the basis functions are the same.
%    poly(c(x,...), 3)[1:length(x),]
%}
%
%\dontrun{
%fit4 <- vglm(cbind(leye,reye) ~ POLY3(op), trace=TRUE,
%            fam = binom2.or(exchangeable=TRUE, zero=3),  data=eyesdata,
%            xij = list(POLY3(op) ~ POLY3(lop,rop) + POLY3(rop,lop) +
%                                   fill(POLY3(lop,rop))),
%            form2 =  ~ POLY3(op) + POLY3(lop,rop) + POLY3(rop,lop) +
%                                   fill(POLY3(lop,rop)))
%coef(fit4)
%coef(fit4, matrix=TRUE)
%head(predict(fit4))
%}


\name{fill}
\alias{fill}
\alias{fill1}
%- \alias{fill2}
%- \alias{fill3}
%- \alias{fill4}
%- \alias{fill5}
%- \alias{fill6}
%- \alias{fill7}
%- \alias{fill8}
%- \alias{fill9}
%- \alias{fill10}
%- \alias{fill11}
%- \alias{fill12}
%- \alias{fill13}
%- \alias{fill14}
%- \alias{fill15}
%- \alias{fill16}
%- \alias{fill17}
%- \alias{fill18}
%- \alias{fill19}
%- \alias{fill20}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Creates a Matrix of Appropriate Dimension }
\description{
  A support function for the argument \code{xij}, it generates a matrix
  of an appropriate dimension.


}
\usage{
fill(x, values = 0, ncolx = ncol(x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A vector or matrix which is used to determine the dimension of the
    answer, in particular, the number of rows.  After converting \code{x}
    to a matrix if necessary, the answer is a matrix of \code{values}
    values, of dimension \code{nrow(x)} by \code{ncolx}.


  }
  \item{values}{
    Numeric.
    The answer contains these values,
    which are recycled \emph{columnwise} if necessary, i.e.,
    as \code{matrix(values, ..., byrow=TRUE)}.


  }
  \item{ncolx}{
    The number of columns of the returned matrix.
    The default is the number of columns of \code{x}.


  }
}
\details{
  The \code{xij} argument for \code{\link{vglm}} allows the user to input
  variables specific to each linear/additive predictor.
  For example, consider
  the bivariate logit model where the first/second linear/additive
  predictor is the logistic regression of the first/second binary response
  respectively. The third linear/additive predictor is \code{log(OR) =
  eta3}, where \code{OR} is the odds ratio.  If one has ocular pressure
  as a covariate in this model then \code{xij} is required to handle the
  ocular pressure for each eye, since these will be different in general.
  [This contrasts with a variable such as \code{age}, the age of the
  person, which has a common value for both eyes.]  In order to input
  these data into \code{\link{vglm}} one often finds that functions
  \code{fill}, \code{fill1}, etc. are useful.


  All terms in the \code{xij}
  and \code{formula} arguments in \code{\link{vglm}}
  must appear in the \code{form2} argument too.


}
\value{
  \code{matrix(values, nrow=nrow(x), ncol=ncolx)}, i.e., a matrix
  consisting of values \code{values}, with the number of rows matching
  \code{x}, and the default number of columns is the number of columns
  of \code{x}.


}
%\references{
%  More information can be found at
%  \url{http://www.stat.auckland.ac.nz/~yee}.
%
%
%}

% \section{Warning }{
% Care is needed in such cases.
% See the examples below.
%
%}

\author{ T. W. Yee }
\note{
  The effect of the \code{xij} argument is after other arguments such as
  \code{exchangeable} and \code{zero}.
  Hence \code{xij} does not affect constraint matrices.


  Additionally, there are currently 3 other identical \code{fill}
  functions, called \code{fill1}, \code{fill2} and \code{fill3};
  if you need more then assign \code{fill4 = fill5 = fill1} etc.
  The reason for this is that if more than one \code{fill} function is
  needed then they must be unique.
  For example, if \eqn{M=4} then
  \code{xij = op ~ lop + rop + fill(mop) + fill(mop)} would reduce to
  \code{xij = op ~ lop + rop + fill(mop)}, whereas
  \code{xij = op ~ lop + rop + fill1(mop) + fill2(mop)} would retain
  all \eqn{M} terms, which is needed.


% The constraint matrices, as returned by \code{constraints}, do not
% have a different meaning when \code{xij} is used.


  In Examples 1 to 3 below, the \code{xij} argument illustrates covariates
  that are specific to a linear predictor. Here, \code{lop}/\code{rop} are
  the ocular pressures of the left/right eye in an artificial dataset,
  and \code{mop} is their mean. Variables \code{leye} and \code{reye}
  might be the presence/absence of a particular disease on the LHS/RHS
  eye respectively.

%
% Examples 1 and 2 are deliberately misspecified.
% The output from, e.g., \code{coef(fit, matrix=TRUE)}, looks wrong but
% is correct because the coefficients are multiplied by the zeros
% produced from \code{fill}.


  In Example 3,
  the \code{xij} argument illustrates fitting the (exchangeable) model
  where there
  is a common smooth function of the ocular pressure. One should use
  regression splines since \code{\link{s}} in \code{\link{vgam}} does not
  handle the \code{xij} argument.  However, regression splines such as
  \code{\link[splines]{bs}} and \code{\link[splines]{ns}} need to have
  the same basis functions here for both functions, and Example 3 illustrates
  a trick involving a function \code{BS} to obtain this, e.g., same knots.
  Although regression splines create more than a single column per term
  in the model matrix, \code{fill(BS(lop,rop))} creates the required
  (same) number of columns.


}

\seealso{
  \code{\link{vglm.control}},
  \code{\link{vglm}},
  \code{\link{multinomial}},
  \code{\link{Select}}.


}
\examples{
fill(runif(5))
fill(runif(5), ncol = 3)
fill(runif(5), val = 1, ncol = 3)

# Generate eyes data for the examples below. Eyes are independent (OR=1).
nn <- 1000  # Number of people
eyesdata <- data.frame(lop = round(runif(nn), 2),
                       rop = round(runif(nn), 2),
                       age = round(rnorm(nn, 40, 10)))
eyesdata <- transform(eyesdata,
    mop = (lop + rop) / 2,        # Mean ocular pressure
    op  = (lop + rop) / 2,        # Value unimportant unless plotting
#   op  =  lop,                   # Choose this if plotting
    eta1 = 0 - 2*lop + 0.04*age,  # Linear predictor for left eye
    eta2 = 0 - 2*rop + 0.04*age)  # Linear predictor for right eye
eyesdata <- transform(eyesdata,
    leye = rbinom(nn, size = 1, prob = logit(eta1, inverse = TRUE)),
    reye = rbinom(nn, size = 1, prob = logit(eta2, inverse = TRUE)))

# Example 1
# All effects are linear
fit1 <- vglm(cbind(leye,reye) ~ op + age,
             family = binom2.or(exchangeable = TRUE, zero = 3),
             data = eyesdata, trace = TRUE,
             xij = list(op ~ lop + rop + fill(lop)),
             form2 =  ~ op + lop + rop + fill(lop) + age)
head(model.matrix(fit1, type = "lm"))   # LM model matrix
head(model.matrix(fit1, type = "vlm"))  # Big VLM model matrix
coef(fit1)
coef(fit1, matrix = TRUE)  # Unchanged with 'xij'
constraints(fit1)
max(abs(predict(fit1)-predict(fit1, new = eyesdata)))  # Predicts correctly
summary(fit1)
\dontrun{
plotvgam(fit1, se = TRUE)  # Wrong, e.g., because it plots against op, not lop.
# So set op = lop in the above for a correct plot.
}



# Example 2
# Model OR as a linear function of mop
fit2 <- vglm(cbind(leye,reye) ~ op + age, data = eyesdata, trace = TRUE,
            binom2.or(exchangeable = TRUE, zero = NULL),
            xij   = list(op ~ lop + rop + mop),
            form2 =    ~ op + lop + rop + mop + age)
head(model.matrix(fit2, type = "lm"))   # LM model matrix
head(model.matrix(fit2, type = "vlm"))  # Big VLM model matrix
coef(fit2)
coef(fit2, matrix = TRUE)  # Unchanged with 'xij'
max(abs(predict(fit2) - predict(fit2, new = eyesdata)))  # Predicts correctly
summary(fit2)
\dontrun{
plotvgam(fit2, se = TRUE)  # Wrong because it plots against op, not lop.
}


# Example 3. This model uses regression splines on ocular pressure.
# It uses a trick to ensure common basis functions.
BS <- function(x, ...)
  sm.bs(c(x,...), df = 3)[1:length(x), , drop = FALSE]  # trick

fit3 <- vglm(cbind(leye,reye) ~ BS(lop,rop) + age,
             family = binom2.or(exchangeable = TRUE, zero = 3),
             data = eyesdata, trace = TRUE,
             xij = list(BS(lop,rop) ~ BS(lop,rop) +
                                      BS(rop,lop) +
                                      fill(BS(lop,rop))),
             form2 = ~  BS(lop,rop) + BS(rop,lop) + fill(BS(lop,rop)) +
                        lop + rop + age)
head(model.matrix(fit3, type =  "lm"))  # LM model matrix
head(model.matrix(fit3, type = "vlm"))  # Big VLM model matrix
coef(fit3)
coef(fit3, matrix = TRUE)
summary(fit3)
fit3@smart.prediction
max(abs(predict(fit3) - predict(fit3, new = eyesdata)))  # Predicts correctly
predict(fit3, new = head(eyesdata))  # Note the 'scalar' OR, i.e., zero=3
max(abs(head(predict(fit3)) -
             predict(fit3, new = head(eyesdata))))  # Should be 0
\dontrun{
plotvgam(fit3, se = TRUE, xlab = "lop")  # Correct
}


# Example 4. Capture-recapture model with ephemeral and enduring
# memory effects. Similar to Yang and Chao (2005), Biometrics.
deermice <- transform(deermice, Lag1 = y1)
M.tbh.lag1 <-
  vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight + Lag1,
       posbernoulli.tb(parallel.t = FALSE ~ 0,
                       parallel.b = FALSE ~ 0,
                       drop.b = FALSE ~ 1),
       xij = list(Lag1 ~ fill(y1) + fill(y2) + fill(y3) + fill(y4) +
                         fill(y5) + fill(y6) +
                         y1 + y2 + y3 + y4 + y5),
       form2 = ~ sex + weight + Lag1 +
                 fill(y1) + fill(y2) + fill(y3) + fill(y4) +
                 fill(y5) + fill(y6) +
                 y1 + y2 + y3 + y4 + y5 + y6,
       data = deermice, trace = TRUE)
coef(M.tbh.lag1)
}
\keyword{models}
\keyword{regression}


%This function is unrelated to the \code{zero} argument found in many
%\pkg{VGAM} family functions.  [zz implies I should call it
%\code{fill1(x, value=0, ncolx=ncol(x))} and create .Rd file for
%\code{zero} argument.]

%eyesdata$leye <- ifelse(runif(n) < exp(eta1)/(1+exp(eta1)), 1, 0)
%eyesdata$reye <- ifelse(runif(n) < exp(eta2)/(1+exp(eta2)), 1, 0)

%   \deqn{logit P(Y_k=1) = f_k(x_{ijk}) }{%
%         logit P(Y_k=1) = f_k(x_{ijk}) }
% for \code{k=1,2}.
%                 fill1(lop, ncol=ncol(BS(lop,rop,mop))), data=eyesdata)

% Models using the \code{xij} argument may or may not predict correctly,
% and inference obtained using \code{summary} may be incorrect.


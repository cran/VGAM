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
    Numeric. The answer contains these values which are recycled if
    necessary.

  }
  \item{ncolx}{ 
    The number of columns of the returned matrix.
    The default is the number of columns of \code{x}.

  }
}
\details{
  The \code{xij} argument for \code{\link{vglm}} allows the user to input
  variables specific to each linear predictor.  For example, consider
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

  All terms in the \code{xij} argument must appear in the main
  \code{formula} argument in \code{\link{vglm}}.

}
\value{
  \code{matrix(values, nrow=nrow(x), ncol=ncolx)}, i.e., a matrix
  consisting of values \code{values}, with the number of rows matching
  \code{x}, and the default number of columns is the number of columns
  of \code{x}.

}
\references{ 
  More information can be found at
  \url{http://www.stat.auckland.ac.nz/~yee}.

}
\section{Warning }{
  The use of the \code{xij} argument overrides other arguments such as
  \code{exchangeable} and \code{zero}.  Care is needed in such cases.
  See the examples below.

}

\author{ T. W. Yee }
\note{ 
  Additionally, there are currently 3 other identical \code{fill}
  functions, called \code{fill1}, \code{fill2} and \code{fill3}; if you
  need more then assign \code{fill4 = fill5 = fill1} etc.
  The reason for this is that if more than one \code{fill} function is
  needed then they must be unique.
  For example, if \eqn{M=4} then 
  \code{xij = op ~ lop + rop + fill(mop) + fill(mop)} would reduce to 
  \code{xij = op ~ lop + rop + fill(mop)}, whereas
  \code{xij = op ~ lop + rop + fill1(mop) + fill2(mop)} would retain
  \eqn{M} terms, which is needed.

  The constraint matrices, as returned by \code{constraints}, have a
  different meaning when \code{xij} is used.

  In Examples 1 to 3 below, the \code{xij} argument illustrates covariates
  that are specific to a linear predictor. Here, \code{lop}/\code{rop} are
  the ocular pressures of the left/right eye in an artificial dataset,
  and \code{mop} is their mean.  Variables \code{leye} and \code{reye}
  might be the presence/absence of a particular disease on the LHS/RHS
  eye respectively.  Examples 1 and 2 are deliberately misspecified.
  The output from, e.g., \code{coef(fit, matrix=TRUE)}, looks wrong but
  is correct because the coefficients are multiplied by the zeros 
  produced from \code{fill}.

  In Example 4,
  the \code{xij} argument illustrates fitting the model where there
  is a common smooth function of the ocular pressure.  One should use
  regression splines since \code{\link{s}} in \code{\link{vgam}} does not
  handle the \code{xij} argument.  However, regression splines such as
  \code{\link[splines]{bs}} and \code{\link[splines]{ns}} need to have
  the same knots here for both functions, and Example 4 illustrates
  a trick involving a function \code{BS} to obtain this. Although
  regression splines create more than a single column per term in the
  model matrix, \code{fill(BS(lop,rop,mop))} creates the required (same)
  number of columns.

}

\seealso{ 
  \code{\link{vglm}},
  \code{\link{vglm.control}}.
}
\examples{
fill(runif(5))
fill(runif(5), ncol=3)
fill(runif(5), val=1, ncol=3)

# Generate eyes data for the examples below. Eyes are independent (OR=1).
set.seed(123)
n = 2000  # Number of people
eyes = data.frame(lop = round(runif(n), 2),
                  rop = round(runif(n), 2),
                  age = round(rnorm(n, 40, 10)))
eyes = transform(eyes,
    mop = (lop + rop) / 2, # mean ocular pressure
    eta1 = 0 - 2*lop + 0.04*age, # Linear predictor for left eye
    eta2 = 0 - 2*rop + 0.04*age) # Linear predictor for right eye
eyes = transform(eyes,
    leye = rbinom(n, size=1, prob=exp(eta1)/(1+exp(eta1))),
    reye = rbinom(n, size=1, prob=exp(eta2)/(1+exp(eta2))))

# Example 1
# Non-exchangeable errors (misspecified model)
fit1 = vglm(cbind(leye,reye) ~ lop + rop + fill(lop) + age,
            family = binom2.or(exchangeable=FALSE, zero=NULL),
            xij = op ~ lop + rop + fill(lop), data=eyes)
model.matrix(fit1, type="lm")[1:7,]   # LM model matrix
model.matrix(fit1, type="vlm")[1:7,]  # Big VLM model matrix
coef(fit1)
coef(fit1, matrix=TRUE)  # Looks wrong but is correct
coef(fit1, matrix=TRUE, compress=FALSE)  # Looks wrong but is correct
constraints(fit1)
max(abs(predict(fit1)-predict(fit1, new=eyes))) # Predicts correctly
summary(fit1)


# Example 2
# Nonexchangeable errors (misspecified model), OR is a function of mop
fit2 = vglm(cbind(leye,reye) ~ lop + rop + mop + age,
            family = binom2.or(exchangeable=FALSE, zero=NULL),
            xij = op ~ lop + rop + mop, data=eyes)
model.matrix(fit2, type="lm")[1:7,]   # LM model matrix
model.matrix(fit2, type="vlm")[1:7,]  # Big VLM model matrix
coef(fit2)
coef(fit2, matrix=TRUE)  # correct
coef(fit2, matrix=TRUE, compress=FALSE)  # correct
max(abs(predict(fit2)-predict(fit2, new=eyes))) # Predicts correctly
summary(fit2)


# Example 3. This model is correctly specified.
# Exchangeable errors
fit3 = vglm(cbind(leye,reye) ~ lop + rop + fill(lop) + age,
            family = binom2.or(exchangeable=TRUE, zero=3),
            xij = op ~ lop + rop + fill(lop), data=eyes)
model.matrix(fit3, type="lm")[1:7,]   # LM model matrix
model.matrix(fit3, type="vlm")[1:7,]  # Big VLM model matrix
coef(fit3)
coef(fit3, matrix=TRUE) # Looks wrong but is correct
coef(fit3, matrix=TRUE, compress=FALSE) # Looks wrong but is correct
predict(fit3, new=eyes[1:4,])  # Note the 'scalar' OR, i.e., zero=3
max(abs(predict(fit3)-predict(fit3, new=eyes))) # Predicts correctly
summary(fit3)


# Example 4. This model uses regression splines on ocular pressure.
# It assumes exchangeable errors.
BS = function(x, ...) bs(c(x,...), df=3)[1:length(x),]
fit4 = vglm(cbind(leye,reye) ~ BS(lop,rop,mop) + BS(rop,lop,mop) +
            fill(BS(lop,rop,mop)) + age,
            family = binom2.or(exchangeable=TRUE, zero=3),
            xij = BS(op) ~ BS(lop,rop,mop) + BS(rop,lop,mop) +
                  fill(BS(lop,rop,mop)), data=eyes)
model.matrix(fit4, type="lm")[1:7,]   # LM model matrix
model.matrix(fit4, type="vlm")[1:7,]  # Big VLM model matrix
coef(fit4)
coef(fit4, matrix=TRUE) # Looks wrong but is correct
coef(fit4, matrix=TRUE, compress=FALSE) # Looks wrong but is correct
predict(fit4, new=eyes[1:4,])  # Note the 'scalar' OR, i.e., zero=3
max(abs(predict(fit4)-predict(fit4, new=eyes))) # Predicts correctly
summary(fit4)
}
\keyword{models}
\keyword{regression}


%This function is unrelated to the \code{zero} argument found in many
%\pkg{VGAM} family functions.  [zz implies I should call it
%\code{fill1(x, value=0, ncolx=ncol(x))} and create .Rd file for
%\code{zero} argument.]

%eyes$leye = ifelse(runif(n) < exp(eta1)/(1+exp(eta1)), 1, 0)
%eyes$reye = ifelse(runif(n) < exp(eta2)/(1+exp(eta2)), 1, 0)

%   \deqn{logit P(Y_k=1) = f_k(x_{ijk}) }{%
%         logit P(Y_k=1) = f_k(x_{ijk}) }
% for \code{k=1,2}.
%                 fill1(lop, ncol=ncol(BS(lop,rop,mop))), data=eyes)

% Models using the \code{xij} argument may or may not predict correctly,
% and inference obtained using \code{summary} may be incorrect.


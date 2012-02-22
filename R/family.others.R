# These functions are Copyright (C) 1998-2012 T. W. Yee   All rights reserved.


# family.others.R
# This file contains functions written by other people.


# Last modified:
# 20110317: (a) James Lauder files.
# This file is in one part: see (a), and later (b), (c).







# ----------------------------------------------------------------------
# (a) James Lauder code put here.
# ----------------------------------------------------------------------
# Edited from james.familyfuncs.2.R on 20110317
# ----------------------------------------------------------------------


# 13/12/10; [drpq]exppois() and exppois().
# reference: Karlis CSDA 53 (2009) pg 894 and 
# Kus CSDA 51 (2007) pg 4497
# everything functioning except for hypergeometric function
# (see R package "hypergeo")


# ref: Kus, section 4.1, pg 4500
# updated on 22/15/2010
dexppois <- function(x, lambda, betave = 1, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  N <- max(length(x), length(lambda), length(betave))
  x <- rep(x, len = N); lambda = rep(lambda, len = N);
  betave <- rep(betave, len = N)

  logdensity <- rep(log(0), len = N)
  xok <- (0 < x)
 
  logdensity[xok] <- log(lambda[xok]) + log(betave[xok]) -
                     log1p(-exp(-lambda[xok])) - lambda[xok] - 
                     betave[xok] * x[xok] + lambda[xok] * 
                     exp(-betave[xok] * x[xok])
   
  logdensity[lambda <= 0] <- NaN
  logdensity[betave <= 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}


# ref: calculated from F(x) from Kus, pg 4499
# updated and working on 22/15/2010
qexppois<- function(p, lambda, betave = 1) {
  ans <- -log(log(p * -(expm1(lambda)) +
         exp(lambda)) / lambda) / betave
  ans[(lambda <= 0) | (betave <= 0)] = NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans
}



# ref: Kus, eqn 2, pg 4499
# Updated on 22/12/2010
pexppois<- function(q, lambda, betave = 1) {
  ans <-(exp(lambda * exp(-betave * q)) - exp(lambda)) / -expm1(lambda)  
  ans[q <= 0] <- 0
  ans[(lambda <= 0) | (betave <= 0)] <- NaN
  ans
}



# ref: calculated from F(x) from Kus, pg 4499
# updated and working on 22/15/2010
rexppois <- function(n, lambda, betave = 1) {
  ans <- -log(log(runif(n) * -(expm1(lambda)) +
         exp(lambda)) / lambda) / betave
  ans[(lambda <= 0) | (betave <= 0)] <- NaN
  ans
}






###################
# the family function
# reference: Karlis CSDA 53 (2009) pg 894 and 
# Kus CSDA 51 (2007) pg 4497
#
# Notes:
# 1. Requires the \pkg{hypergeo} package
# (to use their \code{\link[hypergeo]{genhypergeo}} function).

 exppoisson = function (llambda = "loge", lbetave = "loge",
                        elambda = list(), ebetave = list(),
                        ilambda = 1.1,   ibetave = 2.0,
                        zero = NULL) {

  if (mode(llambda) != "character" && mode(llambda) != "name")
    llambda = as.character(substitute(llambda))
  if (mode(lbetave) != "character" && mode(lbetave) != "name")
    lbetave = as.character(substitute(lbetave))

  if (length(zero) && !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")
  if (length(ilambda) && !is.Numeric(ilambda, positive = TRUE))
    stop("bad input for argument 'ilambda'")
  if (length(ibetave) && !is.Numeric(ibetave, positive = TRUE))
    stop("bad input for argument 'ibetave'")

  ilambda[abs(ilambda - 1) < 0.01] = 1.1
  if (!is.list(ebetave))
    ebetave = list()
  if (!is.list(elambda))
    elambda = list()

#print("hi4, 20110319")
  new("vglmff",
  blurb = c("Exponential Poisson distribution \n \n",
            "Links:    ",
            namesof("lambda", llambda, earg = elambda), ", ",
            namesof("betave", lbetave, earg = ebetave), "\n",
            "Mean:     lambda/(expm1(lambda) * betave)) * ",
                      "genhypergeo(c(1,1),c(2,2),lambda)"),

# genhypergeo() from package: hypergeo
# ref = mean from Kus pg 4499

  constraints = eval(substitute(expression({
    constraints = cm.zero.vgam(constraints, x, .zero , M)
    }), list( .zero = zero))),

  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")

    predictors.names = c(
      namesof("lambda", .llambda, earg = .elambda, short = TRUE),
      namesof("betave", .lbetave, earg = .ebetave, short = TRUE))

    if (!length(etastart)) {
#MLE for lambda from Kus eqn(6) pg 4500
      betave.init = if (length( .ibetave ))
                      rep( .ibetave , len = n) else
                      stop("Need to input a value into argument 'ibetave'")
                   ## (lambda.init/(expm1(lambda.init) * (y + 1/8))) *
                   ##  genhypergeo(c(1,1),c(2,2),lambda.init)
      lambda.init = if (length( .ilambda ))
                      rep( .ilambda , len = n) else
                      (1/betave.init - mean(y)) / ((y * 
                      exp(-betave.init * y))/n)

# supply inital values for now to get function working

      betave.init = rep(weighted.mean(betave.init, w = w), len = n)
#print("head(lambda.init)")
#print( head(lambda.init) )
#print("head(betave.init)")
#print( head(betave.init) )
      
      etastart = cbind(theta2eta(lambda.init, .llambda ,earg = .elambda ),
                       theta2eta(betave.init, .lbetave ,earg = .ebetave ))

#print("head(etastart, 3)")
#print( head(etastart, 3) )
    }
   }), list( .llambda = llambda, .lbetave = lbetave, 
             .ilambda = ilambda, .ibetave = ibetave, 
             .elambda = elambda, .ebetave = ebetave))), 

  linkinv = eval(substitute(function(eta, extra = NULL) {
    lambda = eta2theta(eta[, 1], .llambda , earg = .elambda )
    betave = eta2theta(eta[, 2], .lbetave , earg = .ebetave )

# warning("returning dud means")
#   mu
#   runif(nrow(eta))

#  20110319; not sure about the following:
#  20110319; and not in .Rd file.
   -lambda * genhypergeo(c(1, 1), c(2, 2), lambda) / (expm1(-lambda) *
    betave)
  }, list( .llambda = llambda, .lbetave = lbetave, 
           .elambda = elambda, .ebetave = ebetave))), 

  last = eval(substitute(expression({
    misc$link =    c(lambda = .llambda , betave = .lbetave )
    misc$earg = list(lambda = .elambda , betave = .ebetave )
    misc$expected = TRUE
    
  }), list( .llambda = llambda, .lbetave = lbetave,
            .elambda = elambda, .ebetave = ebetave))), 

  loglikelihood = eval(substitute(function(mu, y, w, 
                  residuals = FALSE, eta, extra = NULL) {
    lambda = eta2theta(eta[, 1], .llambda , earg = .elambda )
    betave = eta2theta(eta[, 2], .lbetave , earg = .ebetave )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(w * dexppois(x = y, lambda = lambda, betave = betave,
                       log = TRUE))
    }
  }, list( .lbetave = lbetave , .llambda = llambda , 
           .elambda = elambda , .ebetave = ebetave ))), 

  vfamily = c("exppoisson"),

# Updated 22/12/2010
  deriv = eval(substitute(expression({
    lambda = eta2theta(eta[, 1], .llambda , earg = .elambda )
    betave = eta2theta(eta[, 2], .lbetave , earg = .ebetave )
    dl.dbetave = 1/betave - y - y * lambda * exp(-betave * y)
    dl.dlambda = 1/lambda - 1/expm1(lambda) - 1 + exp(-betave * y)
    dbetave.deta = dtheta.deta(betave, .lbetave , earg = .ebetave )
    dlambda.deta = dtheta.deta(lambda, .llambda , earg = .elambda )
    c(w) * cbind(dl.dlambda * dlambda.deta, dl.dbetave * dbetave.deta)
  }), list( .llambda = llambda, .lbetave = lbetave,
            .elambda = elambda, .ebetave = ebetave ))), 

  weight = eval(substitute(expression({
    
# Updated 22/12/2010
    temp1 = -expm1(-lambda)
    
# Ref: Kus, pg 4502, J11
    ed2l.dlambda2 = (1 + exp(2 * lambda) - lambda^2 * exp(lambda) - 2 *
                    exp(lambda)) / (lambda * temp1)^2


# Ref: Kus, pg 4502, J22
    ed2l.dbetave2 = 1 / betave^2 - (lambda^2 * exp(-lambda) / (4 * 
                    betave^2 * temp1)) * 
                    genhypergeo(c(2,2,2),c(3,3,3),lambda) 

# Ref: Kus, pg 4502,J12 
    ed2l.dbetavelambda = (lambda * exp(-lambda) / (4 * betave * temp1)) *
                         genhypergeo(c(2,2),c(3,3),lambda)   

    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M)] = dlambda.deta^2 * ed2l.dlambda2
    wz[, iam(2, 2, M)] = dbetave.deta^2 * ed2l.dbetave2
    wz[, iam(1, 2, M)] = dbetave.deta * dlambda.deta * ed2l.dbetavelambda
    c(w) * wz
  }), list( .zero = zero ))))
}






#=======================================================================
# 14/12/10 [drpq]genray() and genrayleigh().
# References: Kundu and Raqab, CSDA 49 (2005) pg 187 and
# Raqab and Kundu "Burr Type X Distribution Revisited"

# Updated by Thomas 10/01/2011
# Notes:
# 1. scale = 1 / \lambda here, = \delta, say.
# 2. My calculations showed EIM_{12} did not agree with Kundu and
#    Raqab, (2005). So am using nsimEIM.



# Ref: Kundu pg 188
#  Updated 22/12/10
dgenray <- function(x, shape, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  N <- max(length(x), length(shape), length(scale))
  x <- rep(x, len = N)
  shape <- rep(shape, len = N)
  scale <- rep(scale, len = N)

  logdensity <- rep(log(0), len = N)
  if (any(xok <- (x > 0))) {
    temp1 <- x[xok] / scale[xok]
    logdensity[xok] <- log(2) + log(shape[xok]) + log(x[xok]) -
                       2 * log(scale[xok]) - temp1^2  +
                       (shape[xok] - 1) * log1p(-exp(-temp1^2))
  }
  logdensity[(shape <= 0) | (scale <= 0)] <- NaN
  if (log.arg) {
    logdensity
  } else {
     exp(logdensity)
  }
}


#  Ref: Kundu pg 188
#  Updated 22/12/10
pgenray <- function(q, shape, scale = 1) {
  ans <- (-expm1(-(q/scale)^2))^shape
  ans[q <= 0] <- 0
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}


# Ref: Kundu pg 193
# Updated 22/12/10  
qgenray <- function(p, shape, scale = 1) {
  ans <- scale * sqrt(-log1p(-(p^(1/shape))))
  ans[(shape <= 0) | (scale <= 0)] = NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans
}




# Ref: Kundu pg 193
#  Updated 22/12/10
rgenray <- function(n, shape, scale = 1) {
  ans <- qgenray(runif(n), shape = shape, scale = scale)
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}



###################
# The family function
# References: Kundu CSDA 49 (2005) pg 187 & Raqab and 
# Kundu "Burr Type X Distribution Revisited"
# updated 05/01/2011
# updated by Thomas 10/01/2011

genrayleigh.control <- function(save.weight = TRUE, ...)
{
# Because of nsimEIM in @weight
    list(save.weight = save.weight)
}

 genrayleigh = function (lshape = "loge", lscale = "loge",
                         eshape = list(), escale = list(),
                         ishape = NULL,   iscale = NULL,
                         tol12 = 1.0e-05, 
                         nsimEIM = 300, zero = 1) {

  if (mode(lshape) != "character" && mode(lshape) != "name")
    lshape = as.character(substitute(lshape))
  if (mode(lscale) != "character" && mode(lscale) != "name")
    lscale = as.character(substitute(lscale))

  if (length(ishape) && !is.Numeric(ishape, positive = TRUE))
    stop("bad input for argument 'ishape'")
  if (length(iscale) && !is.Numeric(iscale, positive = TRUE)) 
    stop("bad input for argument 'iscale'")

  if (length(zero) && !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")
  if (!is.Numeric(nsimEIM, allowable.length = 1, integer.valued = TRUE) || nsimEIM <= 50)
      stop("'nsimEIM' should be an integer greater than 50")

  if (!is.list(escale))
    escale = list()
  if (!is.list(eshape))
    eshape = list()


  new("vglmff",
  blurb = c("Generalized Rayleigh distribution\n",
            "Links:    ",
            namesof("shape", lshape, earg = eshape), ", ",
            namesof("scale", lscale, earg = escale), "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),

  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1) 
      stop("response must be a vector or a one-column matrix")
       
    predictors.names = c(
      namesof("shape", .lshape , earg = .eshape , short = TRUE),
      namesof("scale", .lscale , earg = .escale , short = TRUE))

# Following getmaxmin method implemented on 06/01/11
    if (!length(etastart)) {
      genrayleigh.Loglikfun = function(scale, y, x, w, extraargs) {
        temp1 <- y / scale
# Equation (7) from Kundu and Raqab, p.190, which derives from their (2).
# It gives the MLE of shape, given Scale.
        shape = -1 / weighted.mean(log1p(-exp(-temp1^2)), w = w)

        ans <- sum(w * (log(2) + log(shape) + log(y) - 2 * log(scale) -
                   temp1^2  + (shape - 1) * log1p(-exp(-temp1^2))))
#print("c(scale, ans)")
#print( c(scale, ans) )
        ans
      }
# Note: problems occur if scale values too close to zero:
      scale.grid = seq(0.2 * stats::sd(y), 5 * stats::sd(y), len = 29)
      scale.init = if (length( .iscale )) .iscale else
                     getMaxMin(scale.grid, objfun = genrayleigh.Loglikfun,
                               y = y, x = x, w = w)
#print("head(scale.init)")
#print( head(scale.init) )
      scale.init = rep(scale.init, length = length(y))
 
      shape.init = if (length( .ishape )) .ishape else
                   -1 / weighted.mean(log1p(-exp(-(y/scale.init)^2)),
                    w = w)
#print("head(shape.init)")
#print( head(shape.init) )
      shape.init = rep(shape.init, length = length(y))

      etastart = cbind(theta2eta(shape.init, .lshape, earg = .eshape),
                       theta2eta(scale.init, .lscale, earg = .escale))
#print(",,,,,,,,,,,,,,,,,,,")
        }
    }), list( .lscale = lscale, .lshape = lshape,
              .iscale = iscale, .ishape = ishape,
              .escale = escale, .eshape = eshape))), 

  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape = eta2theta(eta[, 1], .lshape , earg = .eshape )
    Scale = eta2theta(eta[, 2], .lscale , earg = .escale )
# zz yet to do: Find expression for mean
# Much easier to return the median rather than the mean:
    qgenray(p = 0.5, shape = shape, scale = Scale)
  }, list( .lshape = lshape, .lscale = lscale, 
           .eshape = eshape, .escale = escale ))),

  last = eval(substitute(expression({
    misc$link =    c(shape = .lshape , scale = .lscale )
    misc$earg = list(shape = .eshape , scale = .escale )
    misc$expected = TRUE
    misc$nsimEIM = .nsimEIM
  }), list( .lshape = lshape, .lscale = lscale,
            .eshape = eshape, .escale = escale,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(function(mu, y, w, 
                  residuals = FALSE, eta, extra = NULL) {

    shape = eta2theta(eta[, 1], .lshape , earg = .eshape )
    Scale = eta2theta(eta[, 2], .lscale , earg = .escale )
#print("head(shape, 3)")
#print( head(shape, 3) )
#print("head(Scale, 3)")
#print( head(Scale, 3) )

    if (residuals) stop("loglikelihood residuals",
                        "not implemented yet") else {
      sum(w * dgenray(x = y, shape = shape, scale = Scale, log = TRUE))
    }
  }, list( .lshape = lshape , .lscale = lscale , 
           .eshape = eshape , .escale = escale ))), 
      
  vfamily = c("genrayleigh"),

  deriv = eval(substitute(expression({
    shape = eta2theta(eta[, 1], .lshape , earg = .eshape )
    Scale = eta2theta(eta[, 2], .lscale , earg = .escale )
#print("head(shape, 3)")
#print( head(shape, 3) )
#print("head(Scale, 3)")
#print( head(Scale, 3) )
    dshape.deta = dtheta.deta(shape, .lshape , earg = .eshape )
    dscale.deta = dtheta.deta(Scale, .lscale , earg = .escale )
    dthetas.detas = cbind(dshape.deta, dscale.deta)

# Note: singularities wrt derivatives at shape==0 and zz:
    temp1 <- y / Scale
    temp2 <- exp(-temp1^2)
    temp3 <- temp1^2 / Scale
    AAA   <- 2 * temp1^2 / Scale  # 2 * y^2 / Scale^3
    BBB   <- -expm1(-temp1^2)     # denominator
    dl.dshape = 1/shape + log1p(-temp2)
    dl.dscale = -2 / Scale + AAA * (1 - (shape - 1) * temp2 / BBB)

# Special fixup:
    dl.dshape[!is.finite(dl.dshape)] = max(dl.dshape[is.finite(dl.dshape)])

    answer <- c(w) * cbind(dl.dshape, dl.dscale) * dthetas.detas
#print("summary(answer)")
#print( summary(answer) )
#print("head(answer, 3)")
#print( head(answer, 3) )
    answer
  }), list( .lshape = lshape , .lscale = lscale,
            .eshape = eshape,  .escale = escale ))),

  weight = eval(substitute(expression({
# 20110108; I disagree with EIM_{12} of pg 190 of Kundu and Raqab.
# So am using simulated Fisher scoring.

# Notes:
# 1. Inf occurs (albeit infequently) for dl.dshape when ysim is close to 0
#    Special fixup to handle this.

    run.varcov = 0
    ind1 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    for(ii in 1:( .nsimEIM )) {
        ysim = rgenray(n = n, shape = shape, scale = Scale)

        temp1 <- ysim / Scale
        temp2 <- exp(-temp1^2)  # May be 1 if ysim is very close to 0.
        temp3 <- temp1^2 / Scale
        AAA   <- 2 * temp1^2 / Scale  # 2 * y^2 / Scale^3
        BBB   <- -expm1(-temp1^2)     # denominator
        dl.dshape = 1/shape + log1p(-temp2)
        dl.dscale = -2 / Scale + AAA * (1 - (shape - 1) * temp2 / BBB)

# Special fixup:
        dl.dshape[!is.finite(dl.dshape)] = max(
        dl.dshape[is.finite(dl.dshape)])

        temp3 = cbind(dl.dshape, dl.dscale)
        run.varcov = run.varcov + temp3[, ind1$row.index] *
                                  temp3[, ind1$col.index]
    }
    run.varcov = run.varcov / .nsimEIM

    wz = if (intercept.only)
        matrix(colMeans(run.varcov, na.rm = FALSE),
               n, ncol(run.varcov), byrow = TRUE) else run.varcov
    wz = wz * dthetas.detas[, ind1$row] * dthetas.detas[, ind1$col]
#print("summary(run.varcov)")
#print( summary(run.varcov) )
#print("summary(wz)")
#print( summary(wz) )
#print("head(wz,3)")
#print( head(wz,3) )
    c(w) * wz
  }), list( .lshape = lshape , .lscale = lscale,
            .eshape = eshape,  .escale = escale,
            .tol12 = tol12, .nsimEIM = nsimEIM ))))
}







#=======================================================================
# 20/01/10; [drpq]expgeom() and expgeometric().
# Reference: Adamidis and Loukas, SPL 39 (1998) pg 35--42 

# Notes:
# Scale is the reciprocal of scale in Adamidis.
# Updated and working 03/02/2011


# Ref: Adamidis pg.36
dexpgeom <- function(x, scale = 1, shape, log = FALSE) {
# 20110201; looks okay.
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  N <- max(length(x), length(scale), length(shape))
  x <- rep(x, len = N)
  scale <- rep(scale, len = N)
  shape <- rep(shape, len = N)

  logdensity <- rep(log(0), len = N)
  if (any(xok <- (x > 0))) {
    temp1 <- -x[xok] / scale[xok]
    logdensity[xok] <- -log(scale[xok]) + log1p(-shape[xok]) + 
                       temp1 - 2 * log1p(-shape[xok] * exp(temp1))
  }

  logdensity[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  if (log.arg) {
    logdensity
  } else {
     exp(logdensity)
  }
}


# Ref: Adamidis p.37, (3.1)
pexpgeom <- function(q, scale = 1, shape) {
# 20110201; looks okay.
  temp1 <- -q / scale
  ans <- -expm1(temp1) / (1 - shape * exp(temp1))
  ans[q <= 0] <- 0
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}

 
qexpgeom <- function(p, scale = 1, shape) {
# 20110201; looks okay.
  ans <- (-scale) * log((p - 1) / (p * shape - 1))
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans
}


rexpgeom <- function(n, scale = 1, shape) {
  ans <- qexpgeom(runif(n), shape = shape, scale = scale)
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}



#=================================================================
# Exponential geometric family function.
# Reference: Adamidis & Loukas, SPL 39 (1998) pg 35--42 
# All derivatives etc copied directly from article
# Updated and working 03/02/2011

# Notes:
# Scale is the reciprocal of scale in Adamidis.

expgeometric.control <- function(save.weight = TRUE, ...)
{
# Because of nsimEIM in @weight
    list(save.weight = save.weight)
}


 expgeometric = function (lscale = "loge", lshape = "logit",
                          escale = list(), eshape = list(),
                          iscale = NULL,   ishape = NULL, 
                          tol12 = 1.0e-05, zero = 1,
                          nsimEIM = 400) {

# 20110102; modified by TWYee. Works.
# Yet to do: get proper Fisher scoring going.

  if (mode(lshape) != "character" && mode(lshape) != "name")
    lshape = as.character(substitute(lshape))
  if (mode(lscale) != "character" && mode(lscale) != "name")
    lscale = as.character(substitute(lscale))

  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE) || any(ishape >= 1))
      stop("bad input for argument 'ishape'")

  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")

  if (length(zero) && !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")

  if (!is.list(escale))
    escale = list()
  if (!is.list(eshape))
    eshape = list()

  if (!is.Numeric(nsimEIM, allowable.length = 1, integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
      stop("'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Exponential geometric distribution\n\n",
            "Links:    ",
            namesof("Scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     ", "(shape - 1) * log(1 - ",
            "shape) / (shape / Scale)"), 
# mean = Adamidis eqn. (3.2)
                           
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
 

  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")

    predictors.names = c(
      namesof("Scale", .lscale , earg = .escale , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {

      scale.init = if (is.Numeric( .iscale , positive = TRUE)) {
                     rep( .iscale , len = n)
                   } else {
# The scale parameter should be 
# the standard deviation of y.
                      stats::sd(y)  # The papers scale parameter beta
                   }
#print("head(scale.init)")
#print( head(scale.init) )

      shape.init = if (is.Numeric( .ishape , positive = TRUE)) {
                     rep( .ishape , len = n)
                   } else {
# Use the formula for the median:
                      rep(2 - exp(median(y)/scale.init), len = n)
                   }
# But avoid extremes:
      shape.init[shape.init >= 0.95] = 0.95
      shape.init[shape.init <= 0.05] = 0.05

#print("head(shape.init)")
#print( head(shape.init) )
      
      etastart = cbind(theta2eta(scale.init, .lscale , earg = .escale ),
                       theta2eta(shape.init, .lshape , earg = .eshape ))

#print("head(etastart, 3)")
#print( head(etastart, 3) )
    }
   }), list( .lscale = lscale, .lshape = lshape, 
             .iscale = iscale, .ishape = ishape, 
             .escale = escale, .eshape = eshape))), 

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale = eta2theta(eta[, 1], .lscale , earg = .escale )
    shape = eta2theta(eta[, 2], .lshape , earg = .eshape )
    
# Return the mean as fitted value; Adamidis Equation (3.2)
    (shape - 1) * log1p(-shape) / (shape / Scale)

  }, list( .lscale = lscale, .lshape = lshape, 
           .escale = escale, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link =    c(Scale = .lscale , shape = .lshape )
    misc$earg = list(Scale = .escale , shape = .eshape )
    misc$expected = TRUE
    misc$nsimEIM = .nsimEIM
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(function(mu, y, w, 
                  residuals = FALSE, eta, extra = NULL) {

    Scale = eta2theta(eta[, 1], .lscale , earg = .escale )
    shape = eta2theta(eta[, 2], .lshape , earg = .eshape )
#print("head(shape, 3)")
#print( head(shape, 3) )
#print("head(Scale, 3)")
#print( head(Scale, 3) )
    
    if (residuals) stop("loglikelihood residuals",
                        "not implemented yet") else {
      sum(w * dexpgeom(x = y, scale = Scale, shape = shape, log = TRUE))     
    }
  }, list( .lscale = lscale , .lshape = lshape , 
           .escale = escale , .eshape = eshape ))), 
      
  vfamily = c("expgeometric"),

  deriv = eval(substitute(expression({
    Scale = eta2theta(eta[, 1], .lscale , earg = .escale )
    shape = eta2theta(eta[, 2], .lshape , earg = .eshape )

# JGL calculated:
     temp2 <- exp(-y / Scale)
     temp3 <- shape * temp2
     temp4 <- y / Scale^2
     dl.dscale =  -1 / Scale + temp4 + 2 * temp4 * temp3 / (1 - temp3)
     dl.dshape = -1 / (1 - shape)    + 2 * temp2 / (1 - temp3)

    dscale.deta = dtheta.deta(Scale, .lscale , earg = .escale )            
    dshape.deta = dtheta.deta(shape, .lshape , earg = .eshape )
    dthetas.detas = cbind(dscale.deta, dshape.deta)

    answer <- c(w) * cbind(dl.dscale, dl.dshape) * dthetas.detas
#print("summary(answer)")
#print( summary(answer) )
#print("head(answer, 3)")
#print( head(answer, 3) )
    answer
  }), list( .lscale = lscale , .lshape = lshape,
            .escale = escale,  .eshape = eshape ))),

#######################
  weight = eval(substitute(expression({
  
#EIM copied exactly as Adamidis article page 40
# Yet to do: get this proper Fisher scoring going.

# gls package function "dilog()" used for polylog function..check up
# on this.
# if (FALSE) {

#   ed2l.dscale2 = (3 * shape - 2 * (shape - (1 - shape) *
#                   (gsl::dilog(shape,2)$val))) / (3 * Scale^2 * shape)

#   ed2l.dshape2 = (1 - shape)^(-2) / 3

#   ed2l.dscaleshape = (4 * shape^2 - shape + (1 - shape)^2 *
#                     log1p(-shape)) / (3 * Scale * shape^2 * (1 - shape))

#   wz <- matrix(0, n, dimm(M))
#   wz[, iam(1, 1, M)] = dscale.deta^2 * ed2l.dscale2
#   wz[, iam(2, 2, M)] = dshape.deta^2 * ed2l.dshape2
#   wz[, iam(1, 2, M)] = dscale.deta * dshape.deta * ed2l.dscaleshape
#   c(w) * wz
# }


# 5/10/07: Use simulation to estimate the EIM
# Use an updating formula for the mean and variance
# Ref.: Hastie and Tibshirani, 1990, GAM book, p.35.
# Here, the variance has 'n' in denominator, not 'n-1'.

        run.varcov = 0
        ind1 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)

        if (length( .nsimEIM )) {
# Simulated FS used only if nsimEIM was specified.
            for(ii in 1:( .nsimEIM )) {
                ysim = rexpgeom(n, scale=Scale, shape=shape)

# Now compute some quantities
                temp2 <- exp(-ysim / Scale)
                temp3 <- shape * temp2
                temp4 <- ysim / Scale^2
                dl.dscale =  -1 / Scale + temp4 + 
                             2 * temp4 * temp3 / (1 - temp3)
                dl.dshape = -1 / (1 - shape) + 
                             2 * temp2 / (1 - temp3)

                temp6 = cbind(dl.dscale, dl.dshape)
#print("temp6[1:3,]")
#print( temp6[1:3,] )
                run.varcov = run.varcov +
                           temp6[,ind1$row.index] * temp6[,ind1$col.index]
            }

            run.varcov = run.varcov / .nsimEIM

# Can do even better if it is an intercept-only model
            wz = if (intercept.only)
                matrix(colMeans(run.varcov),
                       n, ncol(run.varcov), byrow = TRUE) else run.varcov

#print("wz[1:3,]")
#print( wz[1:3,] )
            wz = wz * dthetas.detas[, ind1$row] *
                      dthetas.detas[, ind1$col]
#print("using simulation")
        }
#print("wz[1:3,]")
#print( wz[1:3,] )

    c(w) * wz      
  }), list( .nsimEIM = nsimEIM ))))
}






#=======================================================================
# 16/02/10; [drpq]explog() and explogarithmic().
# Reference: Tahmasabi and Rezaei, CSDA 52 (2008) pg 3889--3901

# Notes:
# Scale is the reciprocal of scale in Tahmasabi.


# Ref: Tahmasabi pg.3890
dexplog <- function(x, scale = 1, shape, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  N <- max(length(x), length(scale), length(shape))
  x <- rep(x, len = N)
  scale <- rep(scale, len = N)
  shape <- rep(shape, len = N)

  logdensity <- rep(log(0), len = N)
  if (any(xok <- (x > 0))) {
    temp1 <- -x[xok] / scale[xok]
    logdensity[xok] <- -log(-log(shape[xok])) - log(scale[xok]) + 
                       log1p(-shape[xok]) + temp1 - 
                       log1p(-(1-shape[xok]) * exp(temp1))
  }

  logdensity[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  if (log.arg) {
    logdensity
  } else {
     exp(logdensity)
  }
}


# Ref: Tahmasabi pg. 3890
pexplog <- function(q, scale = 1, shape) {
  ans <- 1 - log1p(-(1-shape) * exp(-q / scale)) / log(shape)
  ans[q <= 0] <- 0
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}



#ref:  Tahmasabi pg. 3892
# 20110319; this was wrong. Corrected by TWY.
qexplog <- function(p, scale = 1, shape) {

# orig is wrong:
# ans <- scale * log((1 - shape) / (1 - shape^p))

# 20110319, twy picked up an error:
  ans <- -scale * (log1p(-shape^(1.0 - p)) - log1p(-shape))

  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans
}



#ref:  Tahmasabi pg. 3892
rexplog <- function(n, scale = 1, shape) {
  ans <- qexplog(runif(n), scale = scale, shape = shape)
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}






#=================================================================
# Exponential logarithmic.
# Reference: Tahmasbi and Rezaei, CSDA 52 (2008) pg 3889--3901

# Notes:
# Scale is the reciprocal of scale in Tahmasabi.

#updated and working 27/02/11
explogarithmic.control <- function(save.weight = TRUE, ...)
{
# Because of nsimEIM in @weight
    list(save.weight = save.weight)
}

 explogarithmic = function (lscale = "loge", lshape = "logit",
                            escale = list(), eshape = list(),
                            iscale = NULL,   ishape = NULL,
                            tol12 = 1.0e-05, zero = 1,
                            nsimEIM = 400) {

  if (mode(lshape) != "character" && mode(lshape) != "name")
    lshape = as.character(substitute(lshape))
  if (mode(lscale) != "character" && mode(lscale) != "name")
    lscale = as.character(substitute(lscale))

  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE) || any(ishape >= 1))
      stop("bad input for argument 'ishape'")

  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")

  if (length(zero) && !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")

  if (!is.list(escale))
    escale = list()
  if (!is.list(eshape))
    eshape = list()

  if (!is.Numeric(nsimEIM, allowable.length = 1, integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
      stop("'nsimEIM' should be an integer greater than 50")

  new("vglmff",
  blurb = c("Exponential logarithmic distribution\n\n",
            "Links:    ",
            namesof("Scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     ", "(-polylog(2, 1 - p) * Scale) / log(shape)"),
# mean = Tahmabasi pg. 3891

  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),

  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")

    predictors.names = c(
      namesof("Scale", .lscale , earg = .escale , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {

      scale.init = if (is.Numeric( .iscale , positive = TRUE)) {
                     rep( .iscale , len = n)
                   } else {
# The scale parameter should be
# the standard deviation of y.
                     stats::sd(y)  
                   }

      shape.init = if (is.Numeric( .ishape , positive = TRUE)) {
                     rep( .ishape , len = n)
                   } else {
# Use the formula for the median (Tahmasabi pg. 3891):
                      rep((exp(median(y)/scale.init) - 1)^2, len = n)
                   }
# But avoid extremes:
      shape.init[shape.init >= 0.95] = 0.95
      shape.init[shape.init <= 0.05] = 0.05

#print("head(scale.init)")
#print( head(scale.init) )
#print("head(shape.init)")
#print( head(shape.init) )

      etastart = cbind(theta2eta(scale.init, .lscale , earg = .escale ),
                       theta2eta(shape.init, .lshape , earg = .eshape ))

#print("head(etastart, 3)")
#print( head(etastart, 3) )
    }
   }), list( .lscale = lscale, .lshape = lshape,
             .iscale = iscale, .ishape = ishape,
             .escale = escale, .eshape = eshape))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale = eta2theta(eta[, 1], .lscale , earg = .escale )
    shape = eta2theta(eta[, 2], .lshape , earg = .eshape )

#  warning("returning dud means")
#    runif(nrow(eta))        
# zz yet to do: Find polylog function

# mean should be the fitted value; Tahmasabi pg. 3891
#    (-polylog(2, 1 - p) * Scale) / log(shape)
# mean contains polylog function therefore return median for now:

    qexplog(p = 0.5, shape = shape, scale = Scale)  

  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link =    c(Scale = .lscale , shape = .lshape )
    misc$earg = list(Scale = .escale , shape = .eshape )
    misc$expected = TRUE
    misc$nsimEIM = .nsimEIM
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(function(mu, y, w,
                  residuals = FALSE, eta, extra = NULL) {

    Scale = eta2theta(eta[, 1], .lscale , earg = .escale )
    shape = eta2theta(eta[, 2], .lshape , earg = .eshape )
    
#print("head(shape, 3)")
#print( head(shape, 3) )
#print("head(Scale, 3)")
#print( head(Scale, 3) )

    if (residuals) stop("loglikelihood residuals",
                        "not implemented yet") else {
      sum(w * dexplog(x = y, scale = Scale, shape = shape, log = TRUE))
    }
  }, list( .lscale = lscale , .lshape = lshape ,
           .escale = escale , .eshape = eshape ))),

  vfamily = c("explogarithmic"),

  deriv = eval(substitute(expression({
    Scale = eta2theta(eta[, 1], .lscale , earg = .escale )
    shape = eta2theta(eta[, 2], .lshape , earg = .eshape )

# JGL calculated:
     temp2 <- exp(-y / Scale)
     temp3 <- y / Scale^2
     temp4 <- 1 - shape
     dl.dscale = (-1 / Scale) + temp3 + (temp4 * temp3 *
                 temp2) / (1 - temp4 * temp2)
     dl.dshape = -1 / (shape * log(shape)) - 1 / temp4 -
                 temp2 / (1 - temp4 * temp2)

    dscale.deta = dtheta.deta(Scale, .lscale , earg = .escale )
    dshape.deta = dtheta.deta(shape, .lshape , earg = .eshape )
    dthetas.detas = cbind(dscale.deta, dshape.deta)

    answer <- c(w) * cbind(dl.dscale, dl.dshape) * dthetas.detas
#print("summary(answer)")
#print( summary(answer) )
#print("head(answer, 3)")
#print( head(answer, 3) )
    answer
  }), list( .lscale = lscale , .lshape = lshape,
            .escale = escale,  .eshape = eshape ))),

#######################
  weight = eval(substitute(expression({


# 5/10/07: Use simulation to estimate the EIM

        run.varcov = 0
        ind1 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)

        if (length( .nsimEIM )) {
# Simulated FS used only if nsimEIM was specified.
            for(ii in 1:( .nsimEIM )) {
                ysim = rexplog(n, scale=Scale, shape=shape)

# Now compute some quantities
                temp2 <- exp(-ysim / Scale)
                temp3 <- ysim / Scale^2
                temp4 <- 1 - shape
                dl.dscale = (-1 / Scale) + temp3 + (temp4 * temp3 *
                             temp2) / (1 - temp4 * temp2)
                dl.dshape = -1 / (shape * log(shape)) - 1 / temp4 -
                            temp2 / (1 - temp4 * temp2)

                temp6 = cbind(dl.dscale, dl.dshape)
#print("temp6[1:3,]")
#print( temp6[1:3,] )
                run.varcov = run.varcov +
                           temp6[,ind1$row.index] * temp6[,ind1$col.index]
            }

            run.varcov = run.varcov / .nsimEIM

# Can do even better if it is an intercept-only model
            wz = if (intercept.only)
                matrix(colMeans(run.varcov),
                       n, ncol(run.varcov), byrow = TRUE) else run.varcov

#print("wz[1:3,]")
#print( wz[1:3,] )
            wz = wz * dthetas.detas[, ind1$row] *
                      dthetas.detas[, ind1$col]
#print("using simulation")
        }
#print("wz[1:3,]")
#print( wz[1:3,] )

    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))
}







#=======================================================================
# 09/02/10; [drpq]weibull3()
# Reference "The Weibull distribution - A Handbook" by Horst Rinne


# 20110319; withdrawing [dpqrt]weibull3() due to regularity conditions not
# being met.


  
#Ref: pg. 30
#working 10/02/2010
dweibull3 <- function(x, location = 0, scale = 1, shape, log = FALSE) {

    log.arg = log
    rm(log)
    dweibull(x = x - location, shape = shape, scale = scale, log = log.arg)
}

# Ref: pg 43
# working 10/02/2010
pweibull3 <- function(q, location = 0, scale = 1, shape) {
  pweibull(q = q - location, scale = scale, shape = shape)
}


# Ref: pg 68
# updated and working 18/02/2010
qweibull3 <- function(p, location = 0, scale = 1, shape) {
  location + qweibull(p = p, shape = shape, scale = scale)
}


# Ref: pg 68
# working 11/02/2010
rweibull3 <- function(n, location = 0, scale = 1, shape) {
  location + rweibull(n = n, shape = shape, scale = scale)
}


#=====================================
# 3-parameter Weibull function
# 07/02/2011

# This code is based on the 2-parameter Weibull function weibull()
# Does not accomodate censoring yet.
# Reference "The Weibull distribution - A Handbook" by Horst Rinne

if (FALSE)
 weibull3    = function(llocation = "identity", lscale = "loge",
                       lshape = "loge", elocation = list(),
                       escale = list(), eshape = list(),
                       ilocation = NULL, iscale = NULL, ishape = NULL,
                       imethod = 1, zero = c(2, 3))
{

  llocat = llocation
  elocat = elocation
  ilocat = ilocation

  if (mode(llocat) != "character" && mode(llocat) != "name")
    llocat = as.character(substitute(llocat))
  if (mode(lscale) != "character" && mode(lscale) != "name")
    lscale = as.character(substitute(lscale))
  if (mode(lshape) != "character" && mode(lshape) != "name")
    lshape = as.character(substitute(lshape))
   
  if (length(zero) && !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
      stop("argument 'imethod' must be 1, 2 or 3")

  if (!is.list(elocat)) elocat = list()
  if (!is.list(eshape)) eshape = list()
  if (!is.list(escale)) escale = list()

  new("vglmff",
  blurb = c("3-parameter Weibull distribution\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale",    lscale, earg = escale), ", ",
            namesof("shape",    lshape, earg = eshape), "\n",
            "Mean:     location + scale * gamma(1 + 1/shape)\n",
            "Variance: scale^2 * (gamma(1 + 2/shape) - ",
                      "gamma(1 + 1/shape)^2)"),
#Ref: Rinne (Mean - pg 77 eqn. 2.64b; Var - pg 89 eqn. 2.88a)
  constraints = eval(substitute(expression({
    constraints = cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
  
  initialize = eval(substitute(expression({
    y = cbind(y)
    if (ncol(y) > 1)
      stop("the response must be a vector or a 1-column matrix")

    if (is.SurvS4(y))
      stop("only uncensored observations are allowed; don't use SurvS4()")

    predictors.names =
      c(namesof("location", .llocat, earg = .elocat, tag = FALSE),
        namesof("scale",    .lscale, earg = .escale, tag = FALSE),
        namesof("shape",    .lshape, earg = .eshape, tag = FALSE))


    if (!length(etastart)) {
#Assigning shape.init, scale.init, locat.init
                
        if ( .imethod == 1) {
# method of moments - Rinne page 464
# working - 22/02/2011
                                  
          if(length( .ishape )) {
             shape.init = rep( .ishape , len = n )
          } else {  
# approximating equation for shape
# eqn (12.10b)
             alpha3 = ((1/n) *(sum((y - mean(y))^3)))/((1/n) * (sum((y - 
                       mean(y))^2)))^(3/2)
# eqn (12.10d)
             temp2 = (alpha3 + 1.14)          
             shape.init = rep(-0.729268 - 0.338679 * alpha3 + 4.96077 * 
                          temp2^(-1.0422) + 0.683609 * 
                          (log(temp2))^2, len = n)
#valid for (0.52 <= shape.init <= 100)
          }                                  
                                 
#eqn (12.9b)
          scale.init = if(length( .iscale )) {
                         rep( .iscale , len = n )
                       } else {
                         rep(stats::sd(y) / sqrt(gamma(1 + 2/shape.init) - 
                             gamma(1 + 1/shape.init)^2) , len = n)
                       } 
                                                           
#eqn (12.8b)
          locat.init = if(length( .ilocat )) { 
                          rep( .ilocat , len = n )
                       } else {
                          rep(mean(y) - scale.init * gamma(1 + 1/shape.init),
                              len = n)
                       }
#location = just below min value if smaller than MOM locat.init                       
          locat.init = pmin(min(y) - 0.05 * diff(range(y)), locat.init)             
        }
        if ( .imethod == 2 || .imethod == 3) {  
        #least squares method for scale and shape
        #with two separate methods for locat
                                       
          #code from weibull (2-parameter) for least squares method                             
          if (!length( .ishape ) || !length( .iscale )) {
            anyc = FALSE  # extra$leftcensored | extra$rightcensored
            i11 = if ( .imethod == 2 || .imethod == 3) anyc else 
                  FALSE 
            # can be all data
            qvec = c(.25, .5, .75)  # Arbitrary; could be made an argument
            init.shape = if (length( .ishape )) .ishape else 1
            ###init.shape???  should be shape.init?
            xvec = log(-log1p(-qvec))
            fit0 = lsfit(x = xvec, y=log(quantile(y[!i11], qvec)))
          }
                                       
            shape.init = rep(if(length( .ishape )) .ishape else 
                         1/fit0$coef["X"], len = n)
            scale.init = rep(if(length( .iscale )) .iscale else 
                         exp(fit0$coef["Intercept"]), len = n)
            locat.init = rep(if(length( .ilocat )) .ilocat else                                   
                           if ( .imethod == 2) {
                             ifelse(min(y)>0, 0.75, 1.25) * min(y)
                           } else {
                             min(y) - 0.05 * diff(range(y))
                           } 
                           , len = n)
         }
#print("min(y)")         
#print( min(y) )                            
#print("head(locat.init)")
#print( head(locat.init) )
#print("head(scale.init)")
#print( head(scale.init) )
#print("head(shape.init)")
#print( head(shape.init) )
                                           
        etastart =
        cbind(theta2eta(locat.init, .llocat, earg = .elocat ),
              theta2eta(scale.init, .lscale, earg = .escale ),
              theta2eta(shape.init, .lshape, earg = .eshape ))

 print("head(etastart, 3)")
 print( head(etastart, 3) )

    }
    }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
              .elocat = elocat, .escale = escale, .eshape = eshape,
              .ilocat = ilocat, .iscale = iscale, .ishape = ishape,
              .imethod = imethod) )),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    locat = eta2theta(eta[, 1], .llocat, earg = .elocat )
    scale = eta2theta(eta[, 2], .lscale, earg = .escale )
    shape = eta2theta(eta[, 3], .lshape, earg = .eshape )

# fitted value = mean (pg.77 eqn. 2.64b)
    locat + scale * gamma(1 + 1/shape)
    
  }, list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
           .elocat = elocat, .escale = escale, .eshape = eshape ) )),
  last = eval(substitute(expression({

# From 2-parameter Weibull code:
    if (regnotok <- any(shape <= 2))
      warning("MLE regularity conditions are violated",
              "(shape <= 2) at the final iteration")

# Putting the MLE warning here good because it could possibly be violated
# only in the early iterations.
# Putting the MLE warning here is bad  because vcov() gets no warning.

    misc$link =    c(location = .llocat, scale = .lscale, shape = .lshape)
    misc$earg = list(location = .elocat, scale = .escale, shape = .eshape)
    misc$RegCondOK = !regnotok   # Save this for later
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape
           ) )),
 
  loglikelihood = eval(substitute(
          function(mu, y, w, residuals = FALSE,eta, extra = NULL) {
    locat = eta2theta(eta[, 1], .llocat, earg = .elocat )      
    scale = eta2theta(eta[, 2], .lscale, earg = .escale )
    shape = eta2theta(eta[, 3], .lshape, earg = .eshape )


# 20110319; Some of this code comes from gev().
    if (any(bad <- (y <= locat))) {
        cat("There are", sum(bad), "range violations in @loglikelihood\n")
        flush.console()
    }
    old.answer =
            sum(bad) * (-1.0e10) + ifelse(any(!bad),
            sum(w[!bad] * dweibull3(x = y[!bad], location = locat[!bad],
                                    scale = scale[!bad],
                                    shape = shape[!bad], log = TRUE)), 0)

#   ell2 = dweibull3(x = y, location = locat, scale = scale, 
#                    shape = shape, log = TRUE)

#pg 405 eqn. 11.4b

#    temp3 = y - locat
#    ell1 = log(shape) - shape * log(scale) + (shape-1) * log(temp3) - 
#            (temp3/scale)^shape

#print("max(abs(ell1 - ell2))")
#print( max(abs(ell1 - ell2)) )

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
#     sum(w * ell2)
      old.answer
    }
  }, list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
           .elocat = elocat, .escale = escale, .eshape = eshape ) )),
  vfamily = c("weibull3"),
  deriv = eval(substitute(expression({
 print("in @deriv")
 print("head(eta, 3)")
 print( head(eta, 3) )
    locat = eta2theta(eta[, 1], .llocat, earg = .elocat )
    scale = eta2theta(eta[, 2], .lscale, earg = .escale )
    shape = eta2theta(eta[, 3], .lshape, earg = .eshape )

    dlocat.deta = dtheta.deta(locat, .llocat, earg = .elocat )
    dshape.deta = dtheta.deta(shape, .lshape, earg = .eshape )
    dscale.deta = dtheta.deta(scale, .lscale, earg = .escale )

# equations from pg 405
    temp4 = shape / scale
    zedd = (y - locat) / scale
 print("min(zedd)")
 print( min(zedd) )
    
    if (min(zedd) <= 0)
       warning("Boundary problem. Taking evasive action.")
        
    dl.dlocat =  (1 - shape) / (y - locat) + temp4 * zedd^(shape - 1)
    dl.dscale = temp4 * (-1 + zedd^shape)
    dl.dshape = 1 / shape + log(abs(zedd)) - log(abs(zedd)) * zedd^shape

    c(w) * cbind( dl.dlocat * dlocat.deta,
                  dl.dscale * dscale.deta,
                  dl.dshape * dshape.deta)
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape ) )),
  weight = eval(substitute(expression({
 print("in @weight of weibull3()")
#   EulerM = 0.57721566490153286 
    EulerM = -digamma(1.0)

 print("head(locat)")
 print( head(locat) )
 print("head(scale)")
 print( head(scale) )
 print("head(shape)")
 print( head(shape) )


    wz = matrix(as.numeric(NA), n, dimm(M)) 

# equations involving location parameter from Horst Rinne pg 410
    temp6 = 1 - 1 / shape
    ed2l.dlocat2 = gamma(1 - 2/shape) * ((shape - 1) / scale)^2
#   ed2l.dshape2 = ((1 - EulerM)^2 + (pi^2)/6)/shape^2 # Kleiber&Kotz (2003)
#   ed2l.dshape2: modified from the 2-parameter weibull code:
    ed2l.dscale2 = (shape/scale)^2
    ed2l.dshape2 = (6 * (EulerM - 1)^2 + pi^2)/(6 * shape^2)
    ed2l.dlocatscale = -gamma(2 - 1/shape) * (shape/scale)^2
    ed2l.dlocatshape = -(1/scale) * temp6 * gamma(temp6) *
                       (1 + digamma(temp6))
    ed2l.dshapescale = (EulerM - 1) / scale
    
    wz[, iam(1,1,M)] = ed2l.dlocat2 * dlocat.deta^2
    wz[, iam(2,2,M)] = ed2l.dscale2 * dscale.deta^2
    wz[, iam(3,3,M)] = ed2l.dshape2 * dshape.deta^2
    wz[, iam(1,2,M)] = ed2l.dlocatscale * dlocat.deta * dscale.deta
    wz[, iam(1,3,M)] = ed2l.dlocatshape * dlocat.deta * dshape.deta
    wz[, iam(2,3,M)] = ed2l.dshapescale * dshape.deta * dscale.deta
    
# Putting the MLE warning here is bad because could possibly be violated
# only in the early iterations.
# Putting MLE warning here is good because vcov() gets another warning.

 print("head(wz)")
 print( head(wz) )

    wz = c(w) * wz
    wz
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape ))))
}


# End of James Lauder code here

#=========================================================================




# ----------------------------------------------------------------------
# (b) Arash code
# 20110615
# TPN.R
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

   ### Two-piece normal (TPN) family 

################      dtpn       ##################################

dtpn <- function(x, location = 0, scale = 1, skewpar = 0.5,
                 log.arg = FALSE) {

# Reference: Arash handnotes

  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
           na.rm = TRUE))
    stop("some parameters out of bound")

# Recycle the vectors to equal lengths
  LLL = max(length(x), length(location), length(scale),
            length(skewpar))
  if (length(x) != LLL) x = rep(x, length = LLL)
  if (length(location) != LLL) location = rep(location, length = LLL)
  if (length(scale) != LLL) scale = rep(scale, length = LLL)
  if (length(skewpar) != LLL) skewpar = rep(skewpar, length = LLL)
    
  zedd <- (x - location) / scale

  log.s1 <-  -zedd^2 / (8 * skewpar^2)
  log.s2 <-  -zedd^2 / (8 * (1 - skewpar)^2)
            
  logdensity <- log.s1
  logdensity[zedd > 0] <- log.s2[zedd > 0]
  
  logdensity <- logdensity -log(scale) - log(sqrt(2 * pi))

  if (log.arg) logdensity else exp(logdensity)
}

################      ptpn         ################################
ptpn <- function(q, location = 0, scale = 1, skewpar = 0.5) {

  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
          na.rm = TRUE))
    stop("some parameters out of bound")

# Reference: Arash handnotes

 zedd <- (q - location) / scale

  s1 <- 2 * skewpar * pnorm(zedd, sd = 2 * skewpar) #/ scale
  s2 <- skewpar + (1 - skewpar) * pgamma(zedd^2 / (8 * (1-skewpar)^2), 0.5)
 
ans <- rep(0.0, length(zedd))
ans[zedd <= 0] <- s1[zedd <= 0]
ans[zedd > 0] <- s2[zedd > 0]

ans
}



##################### qtpn  ############################################
pos <- function(x) ifelse(x > 0, x, 0.0)
 

qtpn <- function(p, location = 0, scale = 1, skewpar = 0.5){

  pp = p
  if (any(pp      <= 0 |
          pp      >= 1 |
          skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
             na.rm = TRUE))
    stop("some parameters out of bound")
    # Recycle the vectors to equal lengths
  LLL = max(length(pp), length(location), length(scale),
            length(skewpar))
  if (length(pp) != LLL) pp = rep(pp, length = LLL)
  if (length(location) != LLL) location = rep(location, length = LLL)
  if (length(scale) != LLL) scale = rep(scale, length = LLL)
  if (length(skewpar) != LLL) skewpar = rep(skewpar, length = LLL)
       
  qtpn <- rep(as.numeric(NA), length(LLL))
  qtpn <- qnorm(pp / (2 * skewpar), sd = 2 * skewpar)
  qtpn[pp > skewpar] <- sqrt(8 * ( 1 - skewpar)^2 * 
                        qgamma(pos( pp - skewpar) / ( 
                        1 - skewpar),.5))[pp > skewpar]
        
   qtpn * scale + location
  
}




###########     rast     ##########################################

rtpn <- function(n, location = 0, scale = 1, skewpar = 0.5) {


  qtpn(p = runif(n), location = location, scale = scale, skewpar = skewpar)
}


### Two-piece normal family function via VGAM

tpnff <- function(llocation = "identity", lscale = "loge",
                 elocation = list(), escale = list(), 
                  pp = 0.5, method.init = 1,  zero = 2)
{
# Arash : At the moment, I am working on two important(In Quant. Reg.)
# parameters of the TPN distribution, I am not worry about the skew 
#  parameter p.     
# Note :  pp = Skewparameter
  if (!is.Numeric(method.init, allowable.length = 1, integer.valued = TRUE, positive = TRUE) ||
      method.init > 4)
       stop("'imethod' must be 1 or 2 or 3 or 4")

  if (!is.Numeric(pp, allowable.length = 1, positive = TRUE))
      stop("bad input for argument 'pp'")

  if (mode(llocation)  !=  "character" && mode(llocation) != "name")
       llocation = as.character(substitute(llocation))
  if (mode(lscale)  !=  "character" && mode(lscale) != "name")
       lscale = as.character(substitute(lscale))
  if (length(zero) && !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
       stop("bad input for argument 'zero'")
  if (!is.list(elocation)) elocation = list()
  if (!is.list(escale)) escale = list()

  new("vglmff",
    blurb = c("Two-piece normal distribution \n\n",
              "Links: ",
              namesof("location",  llocation,  earg = elocation), ", ",
              namesof("scale",     lscale,     earg = escale), "\n\n",
              "Mean: "),
    constraints = eval(substitute(expression({
            constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
      predictors.names <-
         c(namesof("location", .llocat, earg = .elocat, tag = FALSE),
           namesof("scale",    .lscale, earg = .escale, tag = FALSE))
      if (ncol(y <- cbind(y)) != 1)
           stop("response must be a vector or a one-column matrix")
      if (!length(etastart)) {
          junk = lm.wfit(x = x, y = y, w = w)
          scale.y.est <- sqrt( sum(w * junk$resid^2) / junk$df.residual )
          location.init <- if ( .llocat == "loge") pmax(1/1024, y) else {
            if ( .method.init == 3) {
              rep(weighted.mean(y, w), len = n)
            } else if ( .method.init == 2) {
              rep(median(rep(y, w)), len = n)
            } else if ( .method.init == 1) {
              junk$fitted
            } else {
              y
            }
          }
          etastart <- cbind(
               theta2eta(location.init,  .llocat, earg = .elocat),
               theta2eta(scale.y.est,    .lscale, earg = .escale))
      }
    }), list( .llocat = llocation, .lscale = lscale,
              .elocat = elocation, .escale = escale,
              .method.init=method.init ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
      eta2theta(eta[,1], .llocat, earg = .elocat)
    }, list( .llocat = llocation,
             .elocat = elocation, .escale = escale ))),
    last = eval(substitute(expression({
      misc$link     <-    c("location" = .llocat, "scale" = .lscale)
      misc$earg     <- list("location" = .elocat, "scale" = .escale)
      misc$expected <- TRUE
      misc$pp       <- .pp
      misc$method.init <- .method.init
    }), list( .llocat = llocation, .lscale = lscale,
              .elocat = elocation, .escale = escale,
              .pp     = pp,        .method.init = method.init ))),
   loglikelihood = eval(substitute(
     function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
     location <- eta2theta(eta[,1], .llocat, earg = .elocat)
     myscale  <- eta2theta(eta[,2], .lscale, earg = .escale)
     ppay      <- .pp
     if (residuals) stop("loglikelihood residuals not ",
                         "implemented yet") else {
       sum(w * dtpn(y, skewpar = ppay, location = location,  scale = myscale,
                      log.arg = TRUE))
     }
   }, list( .llocat = llocation, .lscale = lscale,
            .elocat = elocation, .escale = escale,
            .pp      = pp ))),
    vfamily = c("tpnff"),
    deriv = eval(substitute(expression({
      mylocat <- eta2theta(eta[,1], .llocat,  earg = .elocat)
      myscale <- eta2theta(eta[,2], .lscale,  earg = .escale)
      mypp    <- .pp

      zedd <- (y - mylocat) / myscale
 #     cond1 <-    (zedd <= 0)
       cond2 <-    (zedd > 0)

      dl.dlocat        <-  zedd / (4 * mypp^2)  # cond1
      dl.dlocat[cond2] <- (zedd / (4 * (1 - mypp)^2))[cond2]
      dl.dlocat        <- dl.dlocat / myscale

      dl.dscale        <-  zedd^2 / (4 * mypp^2)
      dl.dscale[cond2] <- (zedd^2 / (4 * (1 - mypp)^2))[cond2]
      dl.dscale        <- (-1 + dl.dscale) / myscale

      #dl.dpp        <-  zedd^2 /  (4 * mypp^3)
      #dl.dpp[cond2] <- -zedd^2 /  (4 * (1 - mypp)^3)[cond2]
      


      dlocat.deta <- dtheta.deta(mylocat, .llocat, earg = .elocat)
      dscale.deta <- dtheta.deta(myscale, .lscale, earg = .escale)
      ans <-
      w * cbind(dl.dlocat * dlocat.deta,
                dl.dscale * dscale.deta)
      ans
    }), list( .llocat = llocation, .lscale = lscale,
              .elocat = elocation, .escale = escale,
              .pp      = pp ))),
    weight = eval(substitute(expression({
      wz   <- matrix(as.numeric(NA), n, M) # diag matrix; y is one-col too
      temp10 <- mypp * (1 - mypp)
  ed2l.dlocat2        <- 1 / ((4 * temp10) * myscale^2)
  ed2l.dscale2        <- 2 /  myscale^2
# ed2l.dskewpar       <- 1 / temp10
# ed2l.dlocatdskewpar <- (-2 * sqrt(2)) / (temp10 * sqrt(pi) * myscale)
     

    wz[, iam(1,1,M)] <- ed2l.dlocat2 * dlocat.deta^2
    wz[, iam(2,2,M)] <- ed2l.dscale2 * dscale.deta^2
  # wz[, iam(3,3,M)] <- ed2l.dskewpar2 * dskewpa.deta^2
  # wz[, iam(1,3,M)] <-  ed2l.dlocatdskewpar * dskewpar.deta * dlocat.deta
        ans
      w * wz
    })
    
    )))
}



  ########################################################################
# Two-piece normal family function via VGAM (All 3 parameters will estimate)

tpnff3 <- function(llocation = "identity", elocation = list(),
                    lscale   = "loge",     escale    = list(),
                    lskewpar = "identity", eskewpar  = list(),
                    method.init = 1,  zero = 2)
{
  if (!is.Numeric(method.init, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      method.init > 4)
    stop("'imethod' must be 1 or 2 or 3 or 4")

 # if (!is.Numeric(pp, allowable.length = 1, positive = TRUE))
  #    stop("bad input for argument 'pp'")

  if (mode(llocation)  !=  "character" && mode(llocation) != "name")
     llocation = as.character(substitute(llocation))
  if (mode(lscale)  !=  "character" && mode(lscale) != "name")
     lscale = as.character(substitute(lscale))
  if (mode(lskewpar)  !=  "character" && mode(lskewpar) != "name")
     lscale = as.character(substitute(lscale))
  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
     stop("bad input for argument 'zero'")

  if (!is.list(elocation)) elocation = list()
  if (!is.list(escale))    escale    = list()
  if (!is.list(eskewpar))  eskewpar = list()

  new("vglmff",
    blurb = c("Two-piece normal distribution \n\n",
              "Links: ",
              namesof("location", llocation,  earg = elocation), ", ",
              namesof("scale",    lscale,     earg = escale),  ", ",
              namesof("skewpar",  lscale,     earg = eskewpar),  "\n\n",
              "Mean: "),
    constraints = eval(substitute(expression({
            constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
      predictors.names <-
         c(namesof("location", .llocat, earg = .elocat, tag = FALSE),
           namesof("scale",    .lscale, earg = .escale, tag = FALSE),
           namesof("skewpar",  .lskewpar, earg = .eskewpar, tag = FALSE))
      if (ncol(y <- cbind(y)) != 1)
           stop("response must be a vector or a one-column matrix")
      if (!length(etastart)) {
          junk = lm.wfit(x = x, y = y, w = w)
          scale.y.est <- sqrt( sum(w * junk$resid^2) / junk$df.residual )
          location.init <- if ( .llocat == "loge") pmax(1/1024, y) else {
            if ( .method.init == 3) {
              rep(weighted.mean(y, w), len = n)
            } else if ( .method.init == 2) {
              rep(median(rep(y, w)), len = n)
            } else if ( .method.init == 1) {
              junk$fitted
            } else {
              y
            }
          }
          skew.l.in <- sum((y < location.init)) / length(y)
          etastart <- cbind(
               theta2eta(location.init, .llocat,   earg = .elocat),
               theta2eta(scale.y.est,   .lscale,   earg = .escale),
               theta2eta(skew.l.in,     .lskewpar, earg = .escale))
      }
    }), list( .llocat = llocation, .lscale = lscale, .lskewpar = lskewpar,
              .elocat = elocation, .escale = escale, .eskewpar = eskewpar,
              
              .method.init=method.init ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
      eta2theta(eta[,1], .llocat, earg = .elocat)
    }, list( .llocat = llocation,
             .elocat = elocation, .escale = escale ))),
    last = eval(substitute(expression({
      misc$link     <- c("location" = .llocat, "scale" = .lscale, 
                                    "skewpar" = .lskewpar)
      misc$earg     <- list( "location" = .elocat, "scale" = .escale,
                             "skewpar"  = .eskewpar)
      misc$expected <- TRUE
           misc$method.init <- .method.init
    }), list( .llocat = llocation, .lscale = lscale, .lskewpar = lskewpar,
              .elocat = elocation, .escale = escale, .eskewpar = lskewpar,
                      .method.init = method.init ))),
   loglikelihood = eval(substitute(
     function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
     location <- eta2theta(eta[,1], .llocat, earg = .elocat)
     myscale  <- eta2theta(eta[,2], .lscale, earg = .escale)
     myskew  <- eta2theta(eta[,3], .lskewpar, earg = .eskewpar)
 
      if (residuals) stop("loglikelihood residuals not ",
                         "implemented yet") else {
       sum(w * dtpn(y, location = location,  scale = myscale,
                      skewpar = myskew, log.arg = TRUE))
     }
   }, list( .llocat = llocation, .lscale = lscale, .lskewpar = lskewpar,
            .elocat = elocation, .escale = escale, .eskewpar = eskewpar
             ))),
    vfamily = c("tpnff3"),
     deriv = eval(substitute(expression({
      mylocat <- eta2theta(eta[,1], .llocat,  earg = .elocat)
      myscale <- eta2theta(eta[,2], .lscale,  earg = .escale)
     myskew   <- eta2theta(eta[,3], .lskewpar,  earg = .eskewpar)
    

      zedd <- (y - mylocat) / myscale
 #     cond1 <-    (zedd <= 0)
       cond2 <-    (zedd > 0)

      dl.dlocat        <-  zedd / (4 * myskew^2)  # cond1
      dl.dlocat[cond2] <- (zedd / (4 * (1 - myskew)^2))[cond2]
      dl.dlocat        <- dl.dlocat / myscale

      dl.dscale        <-  zedd^2 / (4 * myskew^2)
      dl.dscale[cond2] <- (zedd^2 / (4 * (1 - myskew)^2))[cond2]
      dl.dscale        <- (-1 + dl.dscale) / myscale

      dl.dskewpar      <-     zedd^2 /  (4 * myskew^3)
      dl.dskewpar[cond2] <- (-zedd^2 /  (4 * (1 - myskew)^3))[cond2]
      


      dlocat.deta <- dtheta.deta(mylocat, .llocat, earg = .elocat)
      dscale.deta <- dtheta.deta(myscale, .lscale, earg = .escale)
      dskewpar.deta <- dtheta.deta(myskew, .lskewpar, earg = .eskewpar)
      ans <-
      w * cbind(dl.dlocat * dlocat.deta,
                dl.dscale * dscale.deta,
                dl.dskewpar * dskewpar.deta
                )
      ans
    }), list( .llocat = llocation, .lscale = lscale, .lskewpar = lskewpar,
              .elocat = elocation, .escale = escale, .eskewpar = eskewpar
              ))),
    weight = eval(substitute(expression({
      wz   <- matrix(as.numeric(NA), n, dimm(M)) # diag matrix; y is one-col too
     
      temp10 <- myskew * (1 - myskew)
  ed2l.dlocat2        <- 1 / ((4 * temp10) * myscale^2)
  ed2l.dscale2        <- 2 /  myscale^2
  ed2l.dskewpar2      <- 3 / temp10
  ed2l.dlocatdskewpar <- (-2 * sqrt(2)) / (temp10 * sqrt(pi) * myscale)
     
    print("hello")
   wz[, iam(1,1,M)] <- ed2l.dlocat2 * dlocat.deta^2
   wz[, iam(2,2,M)] <- ed2l.dscale2 * dscale.deta^2
   wz[, iam(3,3,M)] <- ed2l.dskewpar2 * dskewpar.deta^2
   wz[, iam(1,3,M)] <- ed2l.dlocatdskewpar * dskewpar.deta * dlocat.deta
   
       ans
      w * wz
    })
    
    )))
}


# ----------------------------------------------------------------------
# (c) Not yet assigned
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------





# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.











dkumar <- function(x, shape1, shape2, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  N <- max(length(x), length(shape1), length(shape2))
  x <- rep(x, len = N); shape1 <- rep(shape1, len = N);
  shape2 <- rep(shape2, len = N)

  logdensity <- rep(log(0), len = N)
  xok <- (0 <= x & x <= 1)
  logdensity[xok] <- log(shape1[xok]) + log(shape2[xok]) +
                     (shape1[xok] - 1) * log(x[xok]) +
                     (shape2[xok] - 1) * log1p(-x[xok]^shape1[xok])

  logdensity[shape1 <= 0] <- NaN
  logdensity[shape2 <= 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}


rkumar <- function(n, shape1, shape2) {
  ans <- (1 - (runif(n))^(1/shape2))^(1/shape1)
  ans[(shape1 <= 0) | (shape2 <= 0)] <- NaN
  ans
}


qkumar <- function(p, shape1, shape2) {


  ans <- (1.0 - (1.0 - p)^(1/shape2))^(1/shape1)
  ans[(shape1 <= 0) | (shape2 <= 0)] = NaN
  ans[p <  0] <- NaN
  ans[p >  1] <- NaN
  ans
}


pkumar = function(q, shape1, shape2) {

  ans <- 1.0 - (1.0 - q^shape1)^shape2
  ans[q <= 0] <- 0
  ans[q >= 1] <- 1
  ans[(shape1 <= 0) | (shape2 <= 0)] <- NaN
  ans
}





 kumar <- function(lshape1 = "loge", lshape2 = "loge",
                   eshape1 = list(), eshape2 = list(),
                   ishape1 = NULL,   ishape2 = NULL,
                   grid.shape1 = c(0.4, 6.0),
                   tol12 = 1.0e-4, zero = NULL)
{
  if (mode(lshape1) != "character" && mode(lshape1) != "name")
    lshape1 <- as.character(substitute(lshape1))
  if (mode(lshape2) != "character" && mode(lshape2) != "name")
    lshape2 <- as.character(substitute(lshape2))
  if (length(ishape1) &&
     (!is.Numeric(ishape1, allowable.length = 1, positive = TRUE)))
      stop("bad input for argument 'ishape1'")
  if (length(ishape2) && !is.Numeric(ishape2))
    stop("bad input for argument 'ishape2'")

  if (!is.list(eshape1)) eshape1 = list()
  if (!is.list(eshape2)) eshape2 = list()

  if (!is.Numeric(tol12, allowable.length = 1, positive = TRUE))
    stop("bad input for argument 'tol12'")
  if (!is.Numeric(grid.shape1, allowable.length = 2, positive = TRUE))
    stop("bad input for argument 'grid.shape1'")

  new("vglmff",
  blurb = c("Kumaraswamy distribution\n\n",
            "Links:    ",
            namesof("shape1", lshape1, earg = eshape1, tag = FALSE), ", ",
            namesof("shape2", lshape2, earg = eshape2, tag = FALSE), "\n",
            "Mean:     ",
            "shape2 * beta(1+1/shape1, shape2)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    if (ncol(y <- cbind(y)) != 1)
      stop("the response must be a vector or one-column matrix")
    if (any((y <= 0) | (y >=1)))
      stop("the response must be in (0,1)")

    predictors.names <- c(
        namesof("shape1", .lshape1 , earg = .eshape1 , tag = FALSE),
        namesof("shape2", .lshape2 , earg = .eshape2 , tag = FALSE))
    if (!length(etastart)) {



      kumar.Loglikfun <- function(shape1, y, x, w, extraargs) {




           medy <- weighted.mean(y, w)


          shape2 <- log(0.5) / log1p(-(medy^shape1))
          sum(w * (log(shape1) + log(shape2) + (shape1-1)*log(y) +
                  (shape2-1)*log1p(-y^shape1)))
      }

      shape1.grid <- seq( .grid.shape1[1], .grid.shape1[2], len = 19)
      shape1.init <- if (length( .ishape1 )) .ishape1 else
      getMaxMin(shape1.grid, objfun = kumar.Loglikfun, y = y,  x = x, w = w)
      shape1.init <- rep(shape1.init, length = length(y))

       medy <- weighted.mean(y, w)



      shape2.init <- if (length( .ishape2 )) .ishape2 else
        log(0.5) / log1p(-(medy^shape1.init))
      shape2.init <- rep(shape2.init, length = length(y))
      etastart <- cbind(
            theta2eta(shape1.init, .lshape1 , earg = .eshape1 ),
            theta2eta(shape2.init, .lshape2 , earg = .eshape2 ))
    }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .ishape1 = ishape1, .ishape2 = ishape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .grid.shape1 = grid.shape1 ))),
  linkinv = eval(substitute(function(eta, extra = NULL){
    shape1 <- eta2theta(eta[,1], link = .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[,2], link = .lshape2 , earg = .eshape2 )
    shape2 * (base::beta(1 + 1/shape1, shape2))
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  last = eval(substitute(expression({
    misc$link <-    c("shape1" = .lshape1, "shape2" = .lshape2)
    misc$earg <- list("shape1" = .eshape1, "shape2" = .eshape2)
    misc$expected = TRUE
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  loglikelihood = eval(substitute(
          function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    shape1 <- eta2theta(eta[,1], link = .lshape1, earg = .eshape1)
    shape2 <- eta2theta(eta[,2], link = .lshape2, earg = .eshape2)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      sum(w * dkumar(x=y, shape1 = shape1, shape2 = shape2, log = TRUE))
    }
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  vfamily = c("kumar"),
  deriv = eval(substitute(expression({
    shape1 <- eta2theta(eta[,1], link = .lshape1, earg = .eshape1)
    shape2 <- eta2theta(eta[,2], link = .lshape2, earg = .eshape2)
    dshape1.deta <- dtheta.deta(shape1, link = .lshape1, earg = .eshape1)
    dshape2.deta <- dtheta.deta(shape2, link = .lshape2, earg = .eshape2)

    dl.dshape1 <- 1 / shape1 + log(y) - (shape2 - 1) * log(y) *
                  (y^shape1) / (1 - y^shape1)
    dl.dshape2 <- 1 / shape2 + log1p(-y^shape1)

    c(w) * cbind(dl.dshape1 * dshape1.deta,
                 dl.dshape2 * dshape2.deta)
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  weight = eval(substitute(expression({
    ed2l.dshape11 <- (1 + (shape2 / (shape2 - 2)) *
           ((digamma(shape2) -  digamma(2))^2 -
           (trigamma(shape2) - trigamma(2)))) / shape1^2
    ed2l.dshape22 <- 1.0 / shape2^2
    ed2l.dshape12 <-
       -((digamma(1 + shape2) - digamma(2)) / (shape2 - 1.0)) / shape1

    index1 <- (abs(shape2 - 1.0) < .tol12)
    if (any(index1))
      ed2l.dshape12[index1] <- -trigamma(2) / shape1[index1]

    index2 <- (abs(shape2 - 2.0) < .tol12)
    if (any(index2))
      ed2l.dshape11[index2] <-
          (1.0 - 2.0 * psigamma(2.0, deriv = 2)) / shape1[index2]^2

    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M = M)] <- ed2l.dshape11 * dshape1.deta^2
    wz[, iam(2, 2, M = M)] <- ed2l.dshape22 * dshape2.deta^2
    wz[, iam(1, 2, M = M)] <- ed2l.dshape12 * dshape1.deta * dshape2.deta

    c(w) * wz
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .tol12 = tol12 ))))
}




drice <- function(x, vee, sigma, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(vee), length(sigma))
  x <- rep(x, len = N); vee <- rep(vee, len = N);
  sigma <- rep(sigma, len = N)

  logdensity <- rep(log(0), len = N)
  xok <- (x > 0)
  x.abs <- abs(x[xok] * vee[xok] / sigma[xok]^2)
  logdensity[xok] <- log(x[xok]) - 2 * log(sigma[xok]) +
                     (-(x[xok]^2+vee[xok]^2)/(2*sigma[xok]^2)) +
                     log(besselI(x.abs, nu=0, expon.scaled = TRUE)) + x.abs
  logdensity[sigma <= 0] <- NaN
  logdensity[vee < 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}


rrice <- function(n, vee, sigma) {
  if (!is.Numeric(n, integer.valued = TRUE, allowable.length = 1))
    stop("bad input for argument 'n'")
  theta <- 1 # any number
  X <- rnorm(n, mean = vee * cos(theta), sd = sigma)
  Y <- rnorm(n, mean = vee * sin(theta), sd = sigma)
  sqrt(X^2 + Y^2)
}



riceff.control <- function(save.weight = TRUE, ...) {
    list(save.weight = save.weight)
}


 riceff = function(lvee = "loge", lsigma = "loge",
                   evee = list(), esigma = list(),
                   ivee = NULL, isigma = NULL,
                   nsimEIM = 100, zero = NULL)
{
  if (mode(lvee) != "character" && mode(lvee) != "name")
    lvee = as.character(substitute(lvee))
  if (mode(lsigma) != "character" && mode(lsigma) != "name")
    lsigma = as.character(substitute(lsigma))
  if (length(ivee) && !is.Numeric(ivee, positive = TRUE))
    stop("bad input for argument 'ivee'")
  if (length(isigma) && !is.Numeric(isigma, positive = TRUE))
    stop("bad input for argument 'isigma'")
  if (!is.list(evee)) evee = list()
  if (!is.list(esigma)) esigma = list()
  if (!is.Numeric(nsimEIM, allowable.length = 1, integer.valued = TRUE) || nsimEIM <= 50)
    stop("'nsimEIM' should be an integer greater than 50")

  new("vglmff",
  blurb = c("Rice distribution\n\n",
            "Links:    ",
            namesof("vee",   lvee,   earg = evee,   tag = FALSE), ", ", 
            namesof("sigma", lsigma, earg = esigma, tag = FALSE), "\n",
            "Mean:     ",
            "sigma*sqrt(pi/2)*exp(z/2)*((1-z)*",
            "besselI(-z/2,nu=0)-z*besselI(-z/2,nu=1)) where z=-vee^2/(2*sigma^2)"),
  constraints = eval(substitute(expression({
    constraints = cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    if (ncol(y <- cbind(y)) != 1)
      stop("the response must be a vector or one-column matrix")
    if (any((y <= 0)))
      stop("the response must be in (0,Inf)")
    predictors.names = c(
                     namesof("vee", .lvee, earg = .evee, tag = FALSE),
                     namesof("sigma", .lsigma, earg = .esigma, tag = FALSE))
    if (!length(etastart)) {
        riceff.Loglikfun = function(vee, y, x, w, extraargs) {
            sigma.init = sd(rep(y, w))
            sum(w * (log(y) - 2*log(sigma.init) +
                     log(besselI(y*vee/sigma.init^2, nu=0)) -
                     (y^2 + vee^2)/(2*sigma.init^2)))
        }
        vee.grid = seq(quantile(rep(y, w), probs = seq(0, 1, 0.2))["20%"],
                       quantile(rep(y, w), probs = seq(0, 1, 0.2))["80%"], len=11)
        vee.init = if (length( .ivee )) .ivee else
            getMaxMin(vee.grid, objfun = riceff.Loglikfun, y = y,  x = x, w = w)
        vee.init = rep(vee.init, length = length(y))
        sigma.init = if (length( .isigma )) .isigma else
            sqrt(max((weighted.mean(y^2, w) - vee.init^2)/2, 0.001))
        sigma.init = rep(sigma.init, length = length(y))
        etastart = cbind(theta2eta(vee.init,   .lvee,   earg = .evee),
                         theta2eta(sigma.init, .lsigma, earg = .esigma))
    }
  }), list( .lvee = lvee, .lsigma = lsigma,
            .ivee = ivee, .isigma = isigma,
            .evee = evee, .esigma = esigma ))),
  linkinv = eval(substitute(function(eta, extra = NULL){
    vee   = eta2theta(eta[, 1], link = .lvee,   earg = .evee)
    sigma = eta2theta(eta[, 2], link = .lsigma, earg = .esigma)
    temp9 = -vee^2 / (2*sigma^2)


      sigma * sqrt(pi/2) * ((1-temp9) * besselI(-temp9/2,nu=0,expon = TRUE) -
                               temp9 * besselI(-temp9/2,nu=1,expon = TRUE))
  }, list( .lvee = lvee, .lsigma = lsigma,
           .evee = evee, .esigma = esigma ))),
  last = eval(substitute(expression({
    misc$link <-    c("vee" = .lvee, "sigma" = .lsigma)
    misc$earg <- list("vee" = .evee, "sigma" = .esigma)
    misc$expected = TRUE
    misc$nsimEIM = .nsimEIM
  }), list( .lvee = lvee, .lsigma = lsigma,
            .evee = evee, .esigma = esigma, .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
          function(mu,y, w,residuals = FALSE,eta,extra = NULL) {
      vee = eta2theta(eta[,1], link = .lvee, earg = .evee)
      sigma = eta2theta(eta[,2], link = .lsigma, earg = .esigma)
      if (residuals)
        stop("loglikelihood residuals not implemented yet") else {
          sum(w * drice(x=y, vee = vee, sigma = sigma, log = TRUE))
      }
  }, list( .lvee = lvee, .lsigma = lsigma,
           .evee = evee, .esigma = esigma ))),
  vfamily = c("riceff"),
  deriv = eval(substitute(expression({
    vee = eta2theta(eta[,1], link = .lvee, earg = .evee)
    sigma = eta2theta(eta[,2], link = .lsigma, earg = .esigma)
    dvee.deta = dtheta.deta(vee, link = .lvee, earg = .evee)
    dsigma.deta = dtheta.deta(sigma, link = .lsigma, earg = .esigma)
    temp8 = y * vee / sigma^2
    dl.dvee = -vee/sigma^2 + (y/sigma^2) *
              besselI(temp8, nu=1) / besselI(temp8, nu=0)
    dl.dsigma = -2/sigma + (y^2 + vee^2)/(sigma^3) - (2 * temp8 / sigma) *
                besselI(temp8, nu=1) / besselI(temp8, nu=0)
    c(w) * cbind(dl.dvee * dvee.deta,
                 dl.dsigma * dsigma.deta)
  }), list( .lvee = lvee, .lsigma = lsigma,
            .evee = evee, .esigma = esigma, .nsimEIM = nsimEIM ))),
  weight = eval(substitute(expression({
    run.var = run.cov = 0
    for(ii in 1:( .nsimEIM )) {
      ysim = rrice(n, vee = vee, sigma = sigma)
      temp8 = ysim * vee / sigma^2
      dl.dvee = -vee/sigma^2 + (ysim/sigma^2) *
                besselI(temp8, nu=1) / besselI(temp8, nu=0)
      dl.dsigma = -2/sigma + (ysim^2 + vee^2)/(sigma^3) -
                  (2 * temp8 / sigma) *
                  besselI(temp8, nu=1) / besselI(temp8, nu=0)

      rm(ysim)
      temp3 = cbind(dl.dvee, dl.dsigma)
      run.var = ((ii-1) * run.var + temp3^2) / ii
      run.cov = ((ii-1) * run.cov + temp3[,1] * temp3[,2]) / ii
    }
    wz = if (intercept.only)
        matrix(colMeans(cbind(run.var, run.cov)),
               n, dimm(M), byrow = TRUE) else cbind(run.var, run.cov)

    dtheta.detas = cbind(dvee.deta, dsigma.deta)
    index0 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    wz = wz * dtheta.detas[,index0$row] * dtheta.detas[,index0$col]
    c(w) * wz
  }), list( .lvee = lvee, .lsigma = lsigma,
            .evee = evee, .esigma = esigma, .nsimEIM = nsimEIM ))))
}




dskellam = function(x, mu1, mu2, log = FALSE) {
    log.arg = log; rm(log)
    if ( !is.logical( log.arg ) || length( log.arg )!=1 )
        stop("bad input for 'log.arg'")

    L = max(length(x), length(mu1), length(mu2))
    x = rep(x, len = L);
    mu1 = rep(mu1, len = L);
    mu2 = rep(mu2, len = L);
    ok2 <- is.finite(mu1) && is.finite(mu2) & (mu1 >= 0) & (mu2 >= 0)
    ok3 <- (mu1 == 0) & (mu2 >  0)
    ok4 <- (mu1 >  0) & (mu2 == 0)
    ok5 <- (mu1 == 0) & (mu2 == 0)
    if (log.arg) {
        ans = -mu1 - mu2 + 2 * sqrt(mu1*mu2) +
              0.5 * x * log(mu1) - 0.5 * x * log(mu2) +
              log(besselI(2 * sqrt(mu1*mu2), nu = x, expon.scaled = TRUE))
        ans[ok3] = dpois(x = -x[ok3], lambda = mu2[ok3], log = TRUE)
        ans[ok4] = dpois(x = -x[ok4], lambda = mu1[ok4], log = TRUE)
        ans[ok5] = dpois(x =  x[ok5], lambda = 0.0,      log = TRUE)
        ans[x != round(x)] = log(0.0)
    } else {
        ans = (mu1/mu2)^(x/2) * exp(-mu1-mu2 + 2 * sqrt(mu1*mu2)) *
              besselI(2 * sqrt(mu1*mu2), nu = x, expon.scaled = TRUE)
        ans[ok3] = dpois(x = -x[ok3], lambda = mu2[ok3])
        ans[ok4] = dpois(x = -x[ok4], lambda = mu1[ok4])
        ans[ok5] = dpois(x =  x[ok5], lambda = 0.0)
        ans[x != round(x)] = 0.0
    }
    ans[!ok2] = NaN
    ans
}






rskellam = function(n, mu1, mu2) {
    rpois(n, mu1) - rpois(n, mu2)
}



skellam.control <- function(save.weight = TRUE, ...) {
  list(save.weight = save.weight)
}


 skellam = function(lmu1 = "loge", lmu2 = "loge",
                    emu1 = list(), emu2= list(),
                    imu1 = NULL, imu2 = NULL,
                    nsimEIM = 100, parallel = FALSE, zero = NULL)
{
    if (mode(lmu1) != "character" && mode(lmu1) != "name")
        lmu1 = as.character(substitute(lmu1))
    if (mode(lmu2) != "character" && mode(lmu2) != "name")
        lmu2 = as.character(substitute(lmu2))
    if (length(imu1) && !is.Numeric(imu1, positive = TRUE))
        stop("bad input for argument 'imu1'")
    if (length(imu2) && !is.Numeric(imu2, positive = TRUE))
        stop("bad input for argument 'imu2'")
    if (!is.list(emu1)) emu1 = list()
    if (!is.list(emu2)) emu2 = list()
    if (!is.Numeric(nsimEIM, allowable.length = 1, integer.valued = TRUE) ||
        nsimEIM <= 50)
        stop("'nsimEIM' should be an integer greater than 50")

    new("vglmff",
    blurb = c("Skellam distribution\n\n",
           "Links:    ",
           namesof("mu1", lmu1, earg = emu1, tag = FALSE), ", ", 
           namesof("mu2", lmu2, earg = emu2, tag = FALSE), "\n",
           "Mean:     mu1-mu2", "\n",
           "Variance: mu1+mu2"),
    constraints = eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints,
                              intercept.apply = TRUE)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(y <- cbind(y)) != 1)
            stop("the response must be a vector or one-column matrix")
        if (any((y != round(y))))
            stop("the response should be integer-valued")
        predictors.names = c(
                       namesof("mu1", .lmu1, earg = .emu1, tag = FALSE),
                       namesof("mu2", .lmu2, earg = .emu2, tag = FALSE))
        if (!length(etastart)) {
            junk = lm.wfit(x = x, y = y, w = w)
            var.y.est = sum(w * junk$resid^2) / junk$df.residual
            mean.init = weighted.mean(y, w)
            mu1.init = max((var.y.est + mean.init)/2, 0.01)
            mu2.init = max((var.y.est - mean.init)/2, 0.01)
            mu1.init = rep(if(length( .imu1)) .imu1 else mu1.init, length=n)
            mu2.init = rep(if(length( .imu2)) .imu2 else mu2.init, length=n)
            etastart = cbind(theta2eta(mu1.init, .lmu1, earg = .emu1),
                             theta2eta(mu2.init, .lmu2, earg = .emu2))
        }
    }), list( .lmu1 = lmu1, .lmu2 = lmu2,
              .imu1=imu1, .imu2=imu2,
              .emu1 = emu1, .emu2 = emu2 ))),
    linkinv = eval(substitute(function(eta, extra = NULL){
        mu1 = eta2theta(eta[,1], link = .lmu1, earg = .emu1)
        mu2 = eta2theta(eta[,2], link = .lmu2, earg = .emu2)
        mu1 - mu2
    }, list( .lmu1 = lmu1, .lmu2 = lmu2,
             .emu1 = emu1, .emu2 = emu2 ))),
    last = eval(substitute(expression({
        misc$link <-    c("mu1" = .lmu1, "mu2" = .lmu2)
        misc$earg <- list("mu1" = .emu1, "mu2" = .emu2)
        misc$expected = TRUE
        misc$nsimEIM = .nsimEIM
    }), list( .lmu1 = lmu1, .lmu2 = lmu2,
              .emu1 = emu1, .emu2 = emu2, .nsimEIM = nsimEIM ))),
    loglikelihood = eval(substitute(
            function(mu,y, w,residuals = FALSE,eta,extra = NULL) {
        mu1 = eta2theta(eta[,1], link = .lmu1, earg = .emu1)
        mu2 = eta2theta(eta[,2], link = .lmu2, earg = .emu2)
        if (residuals)
          stop("loglikelihood residuals not implemented yet") else {




            if ( is.logical( .parallel ) && length( .parallel )== 1 &&
                .parallel )
                sum(w * log(besselI(2*mu1, nu=y, expon = TRUE))) else
                sum(w * (-mu1 - mu2 +
                        0.5 * y * log(mu1) -
                        0.5 * y * log(mu2) +
                        2 * sqrt(mu1*mu2) +  # Use this when expon = TRUE
                        log(besselI(2 * sqrt(mu1*mu2), nu=y, expon = TRUE))))
            }
    }, list( .lmu1 = lmu1, .lmu2 = lmu2,
             .parallel=parallel,
             .emu1 = emu1, .emu2 = emu2 ))),
    vfamily = c("skellam"),
    deriv = eval(substitute(expression({
        mu1 = eta2theta(eta[,1], link = .lmu1, earg = .emu1)
        mu2 = eta2theta(eta[,2], link = .lmu2, earg = .emu2)
        dmu1.deta = dtheta.deta(mu1, link = .lmu1, earg = .emu1)
        dmu2.deta = dtheta.deta(mu2, link = .lmu2, earg = .emu2)
        temp8 = 2 * sqrt(mu1*mu2)
        temp9 = besselI(temp8, nu=y, expon = TRUE)
        temp7 = (besselI(temp8, nu=y-1, expon = TRUE) +
                 besselI(temp8, nu=y+1, expon = TRUE)) / 2
        temp6 = temp7 / temp9
        dl.dmu1 = -1 + 0.5 * y / mu1 + sqrt(mu2/mu1) * temp6
        dl.dmu2 = -1 - 0.5 * y / mu2 + sqrt(mu1/mu2) * temp6
        c(w) * cbind(dl.dmu1 * dmu1.deta,
                     dl.dmu2 * dmu2.deta)
    }), list( .lmu1 = lmu1, .lmu2 = lmu2,
              .emu1 = emu1, .emu2 = emu2, .nsimEIM = nsimEIM ))),
    weight = eval(substitute(expression({
        run.var = run.cov = 0
        for(ii in 1:( .nsimEIM )) {
            ysim = rskellam(n, mu1=mu1, mu2=mu2)
            temp9 = besselI(temp8, nu=ysim, expon = TRUE)
            temp7 = (besselI(temp8, nu=ysim-1, expon = TRUE) +
                     besselI(temp8, nu=ysim+1, expon = TRUE)) / 2
            temp6 = temp7 / temp9
            dl.dmu1 = -1 + 0.5 * ysim/mu1 + sqrt(mu2/mu1) * temp6
            dl.dmu2 = -1 - 0.5 * ysim/mu2 + sqrt(mu1/mu2) * temp6
            rm(ysim)
            temp3 = cbind(dl.dmu1, dl.dmu2)
            run.var = ((ii-1) * run.var + temp3^2) / ii
            run.cov = ((ii-1) * run.cov + temp3[,1] * temp3[,2]) / ii
        }
        wz = if (intercept.only)
            matrix(colMeans(cbind(run.var, run.cov)),
                   n, dimm(M), byrow = TRUE) else cbind(run.var, run.cov)

        dtheta.detas = cbind(dmu1.deta, dmu2.deta)
        index0 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)
        wz = wz * dtheta.detas[,index0$row] * dtheta.detas[,index0$col]
        c(w) * wz
    }), list( .lmu1 = lmu1, .lmu2 = lmu2,
              .emu1 = emu1, .emu2 = emu2, .nsimEIM = nsimEIM ))))
}




dyules = function(x, rho, log = FALSE) {
  log.arg = log
  rm(log)
  if ( !is.logical( log.arg ) || length( log.arg )!=1 )
    stop("bad input for 'log.arg'")
  if ( log.arg ) {
    ans = log(rho) + lbeta(abs(x), rho+1)
    ans[(x != round(x)) | (x < 1)] = log(0)
  } else {
    ans = rho * beta(x, rho+1)
    ans[(x != round(x)) | (x < 1)] = 0
  }
  ans[!is.finite(rho) | (rho <= 0) | (rho <= 0)] = NA
  ans
}


ryules = function(n, rho) {
  if (!is.Numeric(n, integer.valued = TRUE, allowable.length = 1))
    stop("bad input for argument 'n'")
  rgeom(n, prob = exp(-rexp(n, rate=rho))) + 1
}


pyules = function(q, rho) {
  tq = trunc(q)
  ans = 1 - tq * beta(abs(tq), rho+1)
  ans[q<1] = 0
  ans[(rho <= 0) | (rho <= 0)] = NA
  ans
}




yulesimon.control <- function(save.weight = TRUE, ...) {
  list(save.weight = save.weight)
}


 yulesimon = function(link = "loge", earg = list(), irho = NULL, nsimEIM = 200)
{
    if (length(irho) && !is.Numeric(irho, positive = TRUE))
        stop("argument 'irho' must be > 0")
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()
    if (!is.Numeric(nsimEIM, allowable.length = 1, integer.valued = TRUE) ||
        nsimEIM <= 50)
        stop("'nsimEIM' should be an integer greater than 50")

    new("vglmff",
    blurb = c("Yule-Simon distribution f(y) = rho*beta(y,rho+1), rho>0, y=1,2,..\n\n",
            "Link:    ",
            namesof("p", link, earg =earg), "\n\n",
            "Mean:     rho/(rho-1), provided rho>1\n",
            "Variance: rho^2 / ((rho-1)^2 * (rho-2)), provided rho>2"),
    initialize = eval(substitute(expression({
        y = as.numeric(y)
        if (any(y < 1))
            stop("all y values must be in 1,2,3,...")
        if (any(y != round(y )))
            warning("y should be integer-valued")
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("rho", .link, earg =.earg, tag = FALSE) 

        if (!length(etastart)) {
            wmeany = weighted.mean(y, w) + 1/8
            rho.init = wmeany / (wmeany - 1)
            rho.init = rep( if (length( .irho )) .irho else rho.init, len = n)
            etastart = theta2eta(rho.init, .link, earg =.earg)
        }
    }), list( .link=link, .earg =earg, .irho=irho ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        ans = rho = eta2theta(eta, .link, earg =.earg)
        ans[rho>1] = rho / (rho - 1)
        ans[rho<=1] = NA
        ans
    }, list( .link=link, .earg =earg ))),
    last = eval(substitute(expression({
        misc$link <-    c(rho = .link)
        misc$earg <- list(rho = .earg)
        misc$expected = TRUE
        misc$nsimEIM = .nsimEIM
    }), list( .link=link, .earg =earg, .nsimEIM = nsimEIM ))),
    loglikelihood = eval(substitute(
        function(mu,y, w,residuals= FALSE,eta, extra = NULL) {
        rho = eta2theta(eta, .link, earg =.earg)
        if (residuals)
          stop("loglikelihood residuals not implemented yet") else {
            sum(w * dyules(x=y, rho=rho, log = TRUE))
        }
    }, list( .link=link, .earg =earg ))),
    vfamily = c("yulesimon"),
    deriv = eval(substitute(expression({
        rho = eta2theta(eta, .link, earg =.earg)
        dl.drho = 1/rho + digamma(1+rho) - digamma(1+rho+y)
        drho.deta = dtheta.deta(rho, .link, earg =.earg)
        w * dl.drho * drho.deta
    }), list( .link=link, .earg =earg ))),
    weight = eval(substitute(expression({
        run.var = 0
        for(ii in 1:( .nsimEIM )) {
            ysim = ryules(n, rho=rho)
            dl.drho = 1/rho + digamma(1+rho) - digamma(1+rho+ysim)
            rm(ysim)
            temp3 = dl.drho
            run.var = ((ii-1) * run.var + temp3^2) / ii
        }
        wz = if (intercept.only)
            matrix(colMeans(cbind(run.var)),
                   n, dimm(M), byrow = TRUE) else cbind(run.var)

        wz = wz * drho.deta^2


        c(w) * wz
    }), list( .nsimEIM = nsimEIM ))))
}





dslash <- function(x, mu = 0, sigma = 1, log = FALSE,
                   smallno =.Machine$double.eps*1000){
  log.arg = log
  rm(log)
  if (!is.Numeric(sigma) || any(sigma <= 0))
    stop("'sigma' must be positive")
  L = max(length(x), length(mu), length(sigma))
  x = rep(x, len = L); mu = rep(mu, len = L); sigma = rep(sigma, len = L)
  zedd = (x-mu)/sigma
  if (log.arg)
    ifelse(abs(zedd)<smallno, -log(2*sigma*sqrt(2*pi)),
    log1p(-exp(-zedd^2/2)) - log(sqrt(2*pi)*sigma*zedd^2)) else
    ifelse(abs(zedd)<smallno, 1/(2*sigma*sqrt(2*pi)),
    -expm1(-zedd^2/2)/(sqrt(2*pi)*sigma*zedd^2))
}

pslash <- function(q, mu = 0, sigma = 1){
    if (!is.Numeric(sigma) || any(sigma <= 0))
      stop("'sigma' must be positive")
    L = max(length(q), length(mu), length(sigma))
    q = rep(q, len = L); mu = rep(mu, len = L); sigma = rep(sigma, len = L)
    ans = q * NA
    for (ii in 1:L) {
        temp = integrate(dslash, lower = -Inf, upper = q[ii])
        if (temp$message != "OK") {
            warning("integrate() failed")
        } else
            ans[ii] = temp$value
    }
    ans
}

rslash <- function (n, mu = 0, sigma = 1){
    if (!is.Numeric(n, positive = TRUE, integer.valued = TRUE, allowable.length = 1))
      stop("bad input for argument 'n'")
    if (any(sigma <= 0))
      stop("argument 'sigma' must be positive")
    rnorm(n = n, mean=mu, sd=sigma) / runif(n = n)
}

slash.control <- function(save.weight = TRUE, ...)
{
    list(save.weight = save.weight)
}

 slash = function(lmu = "identity", lsigma = "loge",
                  emu = list(), esigma = list(),
                  imu = NULL, isigma = NULL,
                  iprobs = c(0.1, 0.9),
                  nsimEIM = 250, zero = NULL,
                  smallno = .Machine$double.eps*1000)
{
    if (mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if (mode(lsigma) != "character" && mode(lsigma) != "name")
        lsigma = as.character(substitute(lsigma))
    if (length(isigma) && !is.Numeric(isigma, positive = TRUE))
        stop("'isigma' must be > 0")
    if (length(zero) && !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(emu)) emu = list()
    if (!is.list(esigma)) esigma = list()
    if (!is.Numeric(nsimEIM, allowable.length = 1, integer.valued = TRUE) || nsimEIM <= 50)
        stop("'nsimEIM' should be an integer greater than 50")
    if (!is.Numeric(iprobs, positive = TRUE) || max(iprobs) >= 1 ||
       length(iprobs)!=2)
        stop("bad input for argument 'iprobs'")
    if (!is.Numeric(smallno, positive = TRUE) || smallno > 0.1)
        stop("bad input for argument 'smallno'")

    new("vglmff",
    blurb = c("Slash distribution\n\n",
           "Links:    ",
           namesof("mu",    lmu,    earg = emu,    tag = FALSE), ", ",
           namesof("sigma", lsigma, earg = esigma, tag = FALSE), "\n",
           paste(
           "1-exp(-(((y-mu)/sigma)^2)/2))/(sqrt(2*pi)*sigma*((y-mu)/sigma)^2)",
           "\ty!=mu",
           "\n1/(2*sigma*sqrt(2*pi))",
           "\t\t\t\t\t\t\ty=mu\n")),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(y <- cbind(y)) != 1)
            stop("the response must be a vector or one-column matrix")
        predictors.names = c(
            namesof("mu",    .lmu,    earg = .emu,    tag = FALSE),
            namesof("sigma", .lsigma, earg = .esigma, tag = FALSE))
        if (!length(etastart)) {

            slash.Loglikfun = function(mu, y, x, w, extraargs) {
                sigma = if (is.Numeric(.isigma)) .isigma else
                  max(0.01,
                     ((quantile(rep(y, w), prob = 0.75)/2)-mu)/qnorm(0.75))
                zedd = (y-mu)/sigma
                sum(w * ifelse(abs(zedd)<.smallno,
                               -log(2*sigma*sqrt(2*pi)),
                               log1p(-exp(-zedd^2/2)) -
                               log(sqrt(2*pi) * sigma * zedd^2)))
            }
            iprobs = .iprobs
            mu.grid = quantile(rep(y, w), probs=iprobs)
            mu.grid = seq(mu.grid[1], mu.grid[2], length=100)
            mu.init = if (length( .imu )) .imu else
                      getMaxMin(mu.grid, objfun = slash.Loglikfun,
                                y = y,  x = x, w = w)
            sigma.init = if (is.Numeric(.isigma)) .isigma else
              max(0.01, ((quantile(rep(y, w), prob = 0.75)/2)-mu.init)/qnorm(0.75))
            mu.init = rep(mu.init, length = length(y))
            etastart = matrix(0, n, 2)
            etastart[,1] = theta2eta(mu.init, .lmu, earg =.emu)
            etastart[,2] = theta2eta(sigma.init, .lsigma, earg =.esigma)
        }
    }), list( .lmu = lmu, .lsigma = lsigma,
              .imu = imu, .isigma = isigma,
              .emu = emu, .esigma = esigma,
              .iprobs=iprobs, .smallno = smallno))),
    linkinv = eval(substitute(function(eta, extra = NULL){
        NA * eta2theta(eta[,1], link = .lmu, earg = .emu)
    }, list( .lmu = lmu, .emu = emu ))),
    last = eval(substitute(expression({
        misc$link <-    c("mu" = .lmu, "sigma" = .lsigma)
        misc$earg <- list("mu" = .emu, "sigma" = .esigma)
        misc$expected = TRUE
        misc$nsimEIM = .nsimEIM
    }), list( .lmu = lmu, .lsigma = lsigma,
              .emu = emu, .esigma = esigma, .nsimEIM = nsimEIM ))),
    loglikelihood = eval(substitute(
            function(mu,y, w,residuals = FALSE,eta,extra = NULL) {
        mu = eta2theta(eta[,1], link = .lmu, earg = .emu)
        sigma = eta2theta(eta[,2], link = .lsigma, earg = .esigma)
        zedd = (y - mu) / sigma
        if (residuals)
          stop("loglikelihood residuals not implemented yet") else {
            sum(w * dslash(x=y, mu = mu, sigma = sigma, log = TRUE,
                           smallno = .smallno))
        }
    }, list( .lmu = lmu, .lsigma = lsigma,
             .emu = emu, .esigma = esigma, .smallno = smallno ))),
    vfamily = c("slash"),
    deriv = eval(substitute(expression({
        mu = eta2theta(eta[,1], link = .lmu, earg = .emu)
        sigma = eta2theta(eta[,2], link = .lsigma, earg = .esigma)
        dmu.deta = dtheta.deta(mu, link = .lmu, earg = .emu)
        dsigma.deta = dtheta.deta(sigma, link = .lsigma, earg = .esigma)
        zedd = (y-mu)/sigma
        d3 = deriv3(~ w * log(1-exp(-(((y-mu)/sigma)^2)/2))-
                    log(sqrt(2*pi)*sigma*((y-mu)/sigma)^2),
                    c("mu", "sigma"))
        eval.d3 = eval(d3)
        dl.dthetas =  attr(eval.d3, "gradient")
        dl.dmu = dl.dthetas[,1]
        dl.dsigma = dl.dthetas[,2]
        ind0 = (abs(zedd) < .smallno)
        dl.dmu[ind0] = 0
        dl.dsigma[ind0] = -1/sigma[ind0]
        ans =  c(w) * cbind(dl.dmu    * dmu.deta,
                            dl.dsigma * dsigma.deta)
        ans
    }), list( .lmu = lmu, .lsigma = lsigma,
              .emu = emu, .esigma = esigma, .smallno = smallno ))),
    weight=eval(substitute(expression({
        run.varcov = 0
        ind1 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)
        sd3 = deriv3(~ w * log(1-exp(-(((ysim-mu)/sigma)^2)/2))-
                     log(sqrt(2*pi)*sigma*((ysim-mu)/sigma)^2),
                     c("mu", "sigma"))
        for(ii in 1:( .nsimEIM )) {
            ysim = rslash(n, mu = mu, sigma = sigma)
            seval.d3 = eval(sd3)

            dl.dthetas =  attr(seval.d3, "gradient")
            dl.dmu = dl.dthetas[,1]
            dl.dsigma = dl.dthetas[,2]



            temp3 = cbind(dl.dmu, dl.dsigma)
            run.varcov = ((ii-1) * run.varcov +
                       temp3[,ind1$row.index]*temp3[,ind1$col.index]) / ii
        }
        wz = if (intercept.only)
            matrix(colMeans(run.varcov, na.rm = FALSE),
                   n, ncol(run.varcov), byrow = TRUE) else run.varcov
        dthetas.detas = cbind(dmu.deta, dsigma.deta)
        wz = wz * dthetas.detas[,ind1$row] * dthetas.detas[,ind1$col]
        c(w) * wz
    }), list( .lmu = lmu, .lsigma = lsigma,
              .emu = emu, .esigma = esigma,
              .nsimEIM = nsimEIM, .smallno = smallno ))))
}




dnefghs = function(x, tau, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    N = max(length(x), length(tau))
    x = rep(x, len = N); tau = rep(tau, len = N);

    logdensity = log(sin(pi*tau)) + (1-tau)*x - log(pi) - log1p(exp(x))
    logdensity[tau < 0] = NaN
    logdensity[tau > 1] = NaN
    if (log.arg) logdensity else exp(logdensity)
}



 nefghs <- function(link = "logit", earg = list(), itau = NULL,
                    imethod = 1)
{
    if (length(itau) && !is.Numeric(itau, positive = TRUE) || any(itau >= 1))
        stop("argument 'itau' must be in (0,1)")
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()
    if (!is.Numeric(imethod, allowable.length = 1, integer.valued = TRUE, positive = TRUE) ||
       imethod > 2)
        stop("argument 'imethod' must be 1 or 2")

    new("vglmff",
    blurb = c("Natural exponential family generalized hyperbolic ",
            "secant distribution\n",
            "f(y) = sin(pi*tau)*exp((1-tau)*y)/(pi*(1+exp(y))\n\n",
            "Link:    ",
            namesof("tau", link, earg =earg), "\n\n",
            "Mean:     pi / tan(pi * tau)\n"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("tau", .link, earg =.earg, tag = FALSE) 

        if (!length(etastart)) {
            wmeany = if ( .imethod == 1) weighted.mean(y, w) else
                     median(rep(y, w))
            if (abs(wmeany) < 0.01) wmeany = 0.01
            tau.init = atan(pi / wmeany) / pi + 0.5
            tau.init[tau.init < 0.03] = 0.03
            tau.init[tau.init > 0.97] = 0.97
            tau.init = rep( if (length( .itau )) .itau else tau.init, len = n)
            etastart = theta2eta(tau.init, .link, earg =.earg)
        }
    }), list( .link = link, .earg = earg, .itau = itau,
              .imethod = imethod ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        tau = eta2theta(eta, .link, earg =.earg)
        pi / tan(pi * tau)
    }, list( .link=link, .earg =earg ))),
    last = eval(substitute(expression({
        misc$link <-    c(tau = .link)
        misc$earg <- list(tau = .earg)
        misc$expected = TRUE
        misc$imethod= .imethod
    }), list( .link=link, .earg =earg, .imethod = imethod ))),
    loglikelihood = eval(substitute(
        function(mu,y, w,residuals= FALSE,eta, extra = NULL) {
        tau = eta2theta(eta, .link, earg =.earg)
        if (residuals)
          stop("loglikelihood residuals not implemented yet") else {
            sum(w * dnefghs(x=y, tau=tau, log = TRUE))
        }
    }, list( .link=link, .earg =earg ))),
    vfamily = c("nefghs"),
    deriv = eval(substitute(expression({
        tau = eta2theta(eta, .link, earg =.earg)
        dl.dtau = pi / tan(pi * tau) - y
        dtau.deta = dtheta.deta(tau, .link, earg =.earg)
        w * dl.dtau * dtau.deta
    }), list( .link=link, .earg =earg ))),
    weight = eval(substitute(expression({
        d2l.dtau2 = (pi / sin(pi * tau))^2
        wz = d2l.dtau2 * dtau.deta^2
        c(w) * wz
    }), list( .link = link ))))
}




dlogF = function(x, shape1, shape2, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    logdensity = -shape2*x - lbeta(shape1, shape2) -
                (shape1 + shape2) * log1p(exp(-x))
    if (log.arg) logdensity else exp(logdensity)
}




 logF = function(lshape1 = "loge", lshape2 = "loge",
                 eshape1 = list(), eshape2 = list(),
                 ishape1 = NULL, ishape2 = 1,
                 imethod = 1)
{
  if (length(ishape1) && !is.Numeric(ishape1, positive = TRUE))
    stop("argument 'ishape1' must be positive")
  if ( # length(ishape2) &&
   !is.Numeric(ishape2, positive = TRUE))
    stop("argument 'ishape2' must be positive")
  if (mode(lshape1) != "character" && mode(lshape1) != "name")
    lshape1 = as.character(substitute(lshape1))
  if (mode(lshape2) != "character" && mode(lshape2) != "name")
    lshape2 = as.character(substitute(lshape2))
  if (!is.list(eshape1)) eshape1 = list()
  if (!is.list(eshape2)) eshape2 = list()
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
      stop("argument 'imethod' must be 1 or 2")

  new("vglmff",
  blurb = c("log F distribution\n",
          "f(y) = exp(-shape2*y)/(beta(shape1,shape2)*",
          "(1+exp(-y))^(shape1+shape2))\n\n",
          "Link:    ",
          namesof("shape1", lshape1, earg =eshape1),
          ", ",
          namesof("shape2", lshape2, earg =eshape2),
          "\n\n",
          "Mean:     digamma(shape1) - digamma(shape2)"),
  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")
    predictors.names = c(
      namesof("shape1", .lshape1, earg = .eshape1, tag = FALSE),
      namesof("shape2", .lshape2, earg = .eshape2, tag = FALSE))

    if (!length(etastart)) {
      wmeany = if ( .imethod == 1) weighted.mean(y, w) else
               median(rep(y, w))


      shape1.init = shape2.init = rep( .ishape2, len = n)
      shape1.init = if (length( .ishape1))
                            rep( .ishape1, len = n) else {
                index1 = (y > wmeany)
                shape1.init[index1] = shape2.init[index1] + 1/1
                shape1.init[!index1] = shape2.init[!index1] - 1/1
                shape1.init = pmax(shape1.init, 1/8)
                shape1.init
              }
      etastart =
          cbind(theta2eta(shape1.init, .lshape1, earg = .eshape1),
                theta2eta(shape2.init, .lshape2, earg = .eshape2))
    }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .ishape1 = ishape1, .ishape2 = ishape2,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape1 = eta2theta(eta[,1], .lshape1, earg = .eshape1)
    shape2 = eta2theta(eta[,2], .lshape2, earg = .eshape2)
    digamma(shape1) - digamma(shape2)
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1 = .lshape1, shape2 = .lshape2)
    misc$earg <- list(shape1 = .eshape1, shape2 = .eshape2)
    misc$expected = TRUE
    misc$imethod= .imethod
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu,y, w,residuals= FALSE,eta, extra = NULL) {
    shape1 = eta2theta(eta[,1], .lshape1, earg = .eshape1)
    shape2 = eta2theta(eta[,2], .lshape2, earg = .eshape2)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
        sum(w * dlogF(x=y, shape1 = shape1, shape2 = shape2, log = TRUE))
    }
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  vfamily = c("logF"),
  deriv = eval(substitute(expression({
    shape1 = eta2theta(eta[,1], .lshape1, earg = .eshape1)
    shape2 = eta2theta(eta[,2], .lshape2, earg = .eshape2)
    tmp888 = digamma(shape1 + shape2) - log1p(exp(-y))
    dl.dshape1 = tmp888 - digamma(shape1)
    dl.dshape2 = tmp888 - digamma(shape2) - y
    dshape1.deta = dtheta.deta(shape1, .lshape1, earg = .eshape1)
    dshape2.deta = dtheta.deta(shape2, .lshape2, earg = .eshape2)
    c(w) * cbind(dl.dshape1 * dshape1.deta,
                 dl.dshape2 * dshape2.deta)
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  weight = eval(substitute(expression({
    tmp888 = trigamma(shape1 + shape2)
    d2l.dshape12 = trigamma(shape1) - tmp888
    d2l.dshape22 = trigamma(shape2) - tmp888
    d2l.dshape1shape2 = -tmp888
    wz = matrix(0, n, dimm(M))
    wz[,iam(1,1,M = M)] = d2l.dshape12 * dshape1.deta^2
    wz[,iam(2,2,M = M)] = d2l.dshape22 * dshape2.deta^2
    wz[,iam(1,2,M = M)] = d2l.dshape1shape2 * dshape1.deta *
                                                dshape2.deta
    c(w) * wz
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))))
}





dbenf <- function(x, ndigits = 1, log = FALSE) {
  if (!is.Numeric(ndigits, allowable.length = 1, positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
      stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  log.arg <- log; rm(log)
  ans <- x * NA
  indexTF <- is.finite(x) & (x >= lowerlimit)

  ans[indexTF] <- log10(1 + 1/x[indexTF])
  ans[!is.na(x) & !is.nan(x) &
     ((x < lowerlimit) |
      (x > upperlimit) |
      (x != round(x)))] <- 0.0
  if (log.arg) log(ans) else ans
}


rbenf <- function(n, ndigits = 1) {
  if (!is.Numeric(ndigits, allowable.length = 1, positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
      stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE, allowable.length = 1, positive = TRUE)) 
               stop("bad input for argument 'n'") else n
  myrunif <- runif(use.n)

  ans <- rep(lowerlimit, length = use.n)
  for(ii in (lowerlimit+1):upperlimit) {
      indexTF <- (pbenf(ii-1, ndigits = ndigits) < myrunif) &
                 (myrunif <= pbenf(ii, ndigits = ndigits))
      ans[indexTF] <- ii
  }
  ans
}


pbenf <- function(q, ndigits = 1, log.p = FALSE) {
  if (!is.Numeric(ndigits, allowable.length = 1, positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
      stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)

  ans <- q * NA
  floorq <- floor(q)
  indexTF <- is.finite(q) & (floorq >= lowerlimit)
  ans[indexTF] <- log10(1 + floorq[indexTF]) - ifelse(ndigits == 1, 0, 1)
  ans[!is.na(q) & !is.nan(q) & (q >= upperlimit)] <- 1
  ans[!is.na(q) & !is.nan(q) & (q <  lowerlimit)] <- 0
  if (log.p) log(ans) else ans
}




qbenf <- function(p, ndigits = 1) {
  if (!is.Numeric(ndigits, allowable.length = 1, positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
      stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  bad <- !is.na(p) & !is.nan(p) & ((p < 0) | (p > 1))
  if (any(bad))
    stop("bad input for 'p'")

  ans <- rep(lowerlimit, length = length(p))
  for(ii in (lowerlimit+1):upperlimit) {
    indexTF <- is.finite(p) &
              (pbenf(ii-1, ndigits = ndigits) < p) &
              (p <= pbenf(ii, ndigits = ndigits))
    ans[indexTF] <- ii
  }

  ans[ is.na(p) |  is.nan(p)] <- NA
  ans[!is.na(p) & !is.nan(p) & (p == 0)] <- lowerlimit
  ans[!is.na(p) & !is.nan(p) & (p == 1)] <- upperlimit
  ans
}








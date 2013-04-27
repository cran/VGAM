# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.
















dexppois <- function(x, lambda, betave = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
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


qexppois<- function(p, lambda, betave = 1) {
  ans <- -log(log(p * -(expm1(lambda)) +
         exp(lambda)) / lambda) / betave
  ans[(lambda <= 0) | (betave <= 0)] = NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans
}



pexppois<- function(q, lambda, betave = 1) {
  ans <-(exp(lambda * exp(-betave * q)) -
         exp(lambda)) / -expm1(lambda)  
  ans[q <= 0] <- 0
  ans[(lambda <= 0) | (betave <= 0)] <- NaN
  ans
}



rexppois <- function(n, lambda, betave = 1) {
  ans <- -log(log(runif(n) * -(expm1(lambda)) +
         exp(lambda)) / lambda) / betave
  ans[(lambda <= 0) | (betave <= 0)] <- NaN
  ans
}







 exppoisson <- function(llambda = "loge", lbetave = "loge",
                        ilambda = 1.1,   ibetave = 2.0,
                        zero = NULL) {

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  lbetave <- as.list(substitute(lbetave))
  ebetave <- link2list(lbetave)
  lbetave <- attr(ebetave, "function.name")





  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")

  if (length(ilambda) &&
      !is.Numeric(ilambda, positive = TRUE))
    stop("bad input for argument 'ilambda'")
  if (length(ibetave) &&
      !is.Numeric(ibetave, positive = TRUE))
    stop("bad input for argument 'ibetave'")

  ilambda[abs(ilambda - 1) < 0.01] = 1.1


  new("vglmff",
  blurb = c("Exponential Poisson distribution \n \n",
            "Links:    ",
            namesof("lambda", llambda, earg = elambda), ", ",
            namesof("betave", lbetave, earg = ebetave), "\n",
            "Mean:     lambda/(expm1(lambda) * betave)) * ",
                      "genhypergeo(c(1, 1),c(2, 2),lambda)"),


  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
    }), list( .zero = zero))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <- c(
      namesof("lambda", .llambda, earg = .elambda, short = TRUE),
      namesof("betave", .lbetave, earg = .ebetave, short = TRUE))

    if (!length(etastart)) {
      betave.init <- if (length( .ibetave ))
              rep( .ibetave , len = n) else
              stop("Need to input a value into argument 'ibetave'")
      lambda.init <- if (length( .ilambda ))
                      rep( .ilambda , len = n) else
                      (1/betave.init - mean(y)) / ((y * 
                      exp(-betave.init * y))/n)


      betave.init <- rep(weighted.mean(betave.init, w = w), len = n)
      
      etastart <-
        cbind(theta2eta(lambda.init, .llambda ,earg = .elambda ),
              theta2eta(betave.init, .lbetave ,earg = .ebetave ))

    }
  }), list( .llambda = llambda, .lbetave = lbetave, 
            .ilambda = ilambda, .ibetave = ibetave, 
            .elambda = elambda, .ebetave = ebetave))), 

  linkinv = eval(substitute(function(eta, extra = NULL) {
    lambda <- eta2theta(eta[, 1], .llambda , earg = .elambda )
    betave <- eta2theta(eta[, 2], .lbetave , earg = .ebetave )


    -lambda * genhypergeo(c(1, 1), c(2, 2), lambda) / (expm1(-lambda) *
    betave)
  }, list( .llambda = llambda, .lbetave = lbetave, 
           .elambda = elambda, .ebetave = ebetave))), 

  last = eval(substitute(expression({
    misc$link <-    c(lambda = .llambda , betave = .lbetave )

    misc$earg <- list(lambda = .elambda , betave = .ebetave )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .llambda = llambda, .lbetave = lbetave,
            .elambda = elambda, .ebetave = ebetave))), 

  loglikelihood = eval(substitute(function(mu, y, w, 
                  residuals = FALSE, eta, extra = NULL) {
    lambda <- eta2theta(eta[, 1], .llambda , earg = .elambda )
    betave <- eta2theta(eta[, 2], .lbetave , earg = .ebetave )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dexppois(x = y, lambda = lambda, betave = betave,
                          log = TRUE))
    }
  }, list( .lbetave = lbetave , .llambda = llambda , 
           .elambda = elambda , .ebetave = ebetave ))), 

  vfamily = c("exppoisson"),

  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta[, 1], .llambda , earg = .elambda )
    betave <- eta2theta(eta[, 2], .lbetave , earg = .ebetave )
    dl.dbetave <- 1/betave - y - y * lambda * exp(-betave * y)
    dl.dlambda <- 1/lambda - 1/expm1(lambda) - 1 + exp(-betave * y)
    dbetave.deta <- dtheta.deta(betave, .lbetave , earg = .ebetave )
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    c(w) * cbind(dl.dlambda * dlambda.deta,
                 dl.dbetave * dbetave.deta)
  }), list( .llambda = llambda, .lbetave = lbetave,
            .elambda = elambda, .ebetave = ebetave ))), 

  weight = eval(substitute(expression({
    
    temp1 <- -expm1(-lambda)
    
    ned2l.dlambda2 <- (1 + exp(2 * lambda) - lambda^2 * exp(lambda) - 2 *
                    exp(lambda)) / (lambda * temp1)^2


    ned2l.dbetave2 <- 1 / betave^2 - (lambda^2 * exp(-lambda) / (4 * 
                    betave^2 * temp1)) * 
                    genhypergeo(c(2, 2, 2),c(3, 3, 3),lambda) 

    ned2l.dbetavelambda <- (lambda * exp(-lambda) / (4 * betave * temp1)) *
                         genhypergeo(c(2, 2),c(3, 3),lambda)   

    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M)] <- dlambda.deta^2 * ned2l.dlambda2
    wz[, iam(2, 2, M)] <- dbetave.deta^2 * ned2l.dbetave2
    wz[, iam(1, 2, M)] <- dbetave.deta * dlambda.deta * ned2l.dbetavelambda
    c(w) * wz
  }), list( .zero = zero ))))
}










dgenray <- function(x, shape, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
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


pgenray <- function(q, shape, scale = 1) {
  ans <- (-expm1(-(q/scale)^2))^shape
  ans[q <= 0] <- 0
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}


qgenray <- function(p, shape, scale = 1) {
  ans <- scale * sqrt(-log1p(-(p^(1/shape))))
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans
}




rgenray <- function(n, shape, scale = 1) {
  ans <- qgenray(runif(n), shape = shape, scale = scale)
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}




genrayleigh.control <- function(save.weight = TRUE, ...) {
    list(save.weight = save.weight)
}


 genrayleigh <- function(lshape = "loge", lscale = "loge",
                         ishape = NULL,   iscale = NULL,
                         tol12 = 1.0e-05, 
                         nsimEIM = 300, zero = 1) {

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  if (length(ishape) &&
      !is.Numeric(ishape, positive = TRUE))
    stop("bad input for argument 'ishape'")
  if (length(iscale) &&
      !is.Numeric(iscale, positive = TRUE)) 
    stop("bad input for argument 'iscale'")

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")
  if (!is.Numeric(nsimEIM, allowable.length = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
      stop("argument 'nsimEIM' should be an integer greater than 50")



  new("vglmff",
  blurb = c("Generalized Rayleigh distribution\n",
            "Links:    ",
            namesof("shape", lshape, earg = eshape), ", ",
            namesof("scale", lscale, earg = escale), "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y





    predictors.names <- c(
      namesof("shape", .lshape , earg = .eshape , short = TRUE),
      namesof("scale", .lscale , earg = .escale , short = TRUE))

    if (!length(etastart)) {
      genrayleigh.Loglikfun <- function(scale, y, x, w, extraargs) {
        temp1 <- y / scale
        shape <- -1 / weighted.mean(log1p(-exp(-temp1^2)), w = w)

        ans <- sum(c(w) * (log(2) + log(shape) + log(y) -
                           2 * log(scale) - temp1^2  +
                           (shape - 1) * log1p(-exp(-temp1^2))))
        ans
      }
      scale.grid <- seq(0.2 * stats::sd(c(y)),
                        5.0 * stats::sd(c(y)), len = 29)
      scale.init <- if (length( .iscale )) .iscale else
                    getMaxMin(scale.grid, objfun = genrayleigh.Loglikfun,
                               y = y, x = x, w = w)
      scale.init <- rep(scale.init, length = length(y))
 
      shape.init <- if (length( .ishape )) .ishape else
                    -1 / weighted.mean(log1p(-exp(-(y/scale.init)^2)),
                     w = w)
      shape.init <- rep(shape.init, length = length(y))

      etastart <- cbind(theta2eta(shape.init, .lshape, earg = .eshape),
                        theta2eta(scale.init, .lscale, earg = .escale))
        }
    }), list( .lscale = lscale, .lshape = lshape,
              .iscale = iscale, .ishape = ishape,
              .escale = escale, .eshape = eshape))), 

  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta[, 1], .lshape , earg = .eshape )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    qgenray(p = 0.5, shape = shape, scale = Scale)
  }, list( .lshape = lshape, .lscale = lscale, 
           .eshape = eshape, .escale = escale ))),

  last = eval(substitute(expression({
    misc$link <-    c(shape = .lshape , scale = .lscale )

    misc$earg <- list(shape = .eshape , scale = .escale )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lshape = lshape, .lscale = lscale,
            .eshape = eshape, .escale = escale,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(function(mu, y, w, 
                  residuals = FALSE, eta, extra = NULL) {

    shape <- eta2theta(eta[, 1], .lshape , earg = .eshape )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )

    if (residuals) stop("loglikelihood residuals",
                        "not implemented yet") else {
      sum(c(w) * dgenray(x = y, shape = shape,
                         scale = Scale, log = TRUE))
    }
  }, list( .lshape = lshape , .lscale = lscale , 
           .eshape = eshape , .escale = escale ))), 
      
  vfamily = c("genrayleigh"),

  deriv = eval(substitute(expression({
    shape <- eta2theta(eta[, 1], .lshape , earg = .eshape )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
    dthetas.detas <- cbind(dshape.deta, dscale.deta)

    temp1 <- y / Scale
    temp2 <- exp(-temp1^2)
    temp3 <- temp1^2 / Scale
    AAA   <- 2 * temp1^2 / Scale  # 2 * y^2 / Scale^3
    BBB   <- -expm1(-temp1^2)     # denominator
    dl.dshape <- 1/shape + log1p(-temp2)
    dl.dscale <- -2 / Scale + AAA * (1 - (shape - 1) * temp2 / BBB)

    dl.dshape[!is.finite(dl.dshape)] =
      max(dl.dshape[is.finite(dl.dshape)])

    answer <- c(w) * cbind(dl.dshape, dl.dscale) * dthetas.detas
    answer
  }), list( .lshape = lshape , .lscale = lscale,
            .eshape = eshape,  .escale = escale ))),

  weight = eval(substitute(expression({


    run.varcov <- 0
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    for(ii in 1:( .nsimEIM )) {
        ysim <- rgenray(n = n, shape = shape, scale = Scale)

        temp1 <- ysim / Scale
        temp2 <- exp(-temp1^2)  # May be 1 if ysim is very close to 0.
        temp3 <- temp1^2 / Scale
        AAA   <- 2 * temp1^2 / Scale  # 2 * y^2 / Scale^3
        BBB   <- -expm1(-temp1^2)     # denominator
        dl.dshape <- 1/shape + log1p(-temp2)
        dl.dscale <- -2 / Scale + AAA * (1 - (shape - 1) * temp2 / BBB)

        dl.dshape[!is.finite(dl.dshape)] <- max(
        dl.dshape[is.finite(dl.dshape)])

        temp3 <- cbind(dl.dshape, dl.dscale)
        run.varcov <- run.varcov + temp3[, ind1$row.index] *
                                   temp3[, ind1$col.index]
    }
    run.varcov <- run.varcov / .nsimEIM

    wz <- if (intercept.only)
        matrix(colMeans(run.varcov, na.rm = FALSE),
               n, ncol(run.varcov), byrow = TRUE) else run.varcov
    wz <- wz * dthetas.detas[, ind1$row] * dthetas.detas[, ind1$col]
    c(w) * wz
  }), list( .lshape = lshape , .lscale = lscale,
            .eshape = eshape,  .escale = escale,
            .tol12 = tol12, .nsimEIM = nsimEIM ))))
}










dexpgeom <- function(x, scale = 1, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
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


pexpgeom <- function(q, scale = 1, shape) {
  temp1 <- -q / scale
  ans <- -expm1(temp1) / (1 - shape * exp(temp1))
  ans[q <= 0] <- 0
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}

 
qexpgeom <- function(p, scale = 1, shape) {
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





expgeometric.control <- function(save.weight = TRUE, ...) {
  list(save.weight = save.weight)
}


 expgeometric <- function(lscale = "loge", lshape = "logit",
                          iscale = NULL,   ishape = NULL, 
                          tol12 = 1.0e-05, zero = 1,
                          nsimEIM = 400) {


  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE) || any(ishape >= 1))
      stop("bad input for argument 'ishape'")

  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  if (!is.Numeric(nsimEIM, allowable.length = 1,
                  integer.valued = TRUE))
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
                           
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
 

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y







    predictors.names <- c(
      namesof("Scale", .lscale , earg = .escale , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {

      scale.init <- if (is.Numeric( .iscale , positive = TRUE)) {
                      rep( .iscale , len = n)
                    } else {
                      stats::sd(c(y)) # The papers scale parameter beta
                    }

      shape.init <- if (is.Numeric( .ishape , positive = TRUE)) {
                      rep( .ishape , len = n)
                    } else {
                      rep(2 - exp(median(y)/scale.init), len = n)
                    }
      shape.init[shape.init >= 0.95] <- 0.95
      shape.init[shape.init <= 0.05] <- 0.05

      
      etastart <-
        cbind(theta2eta(scale.init, .lscale , earg = .escale ),
              theta2eta(shape.init, .lshape , earg = .eshape ))

    }
   }), list( .lscale = lscale, .lshape = lshape, 
             .iscale = iscale, .ishape = ishape, 
             .escale = escale, .eshape = eshape))), 

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    
    (shape - 1) * log1p(-shape) / (shape / Scale)

  }, list( .lscale = lscale, .lshape = lshape, 
           .escale = escale, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link <-    c(Scale = .lscale , shape = .lshape )

    misc$earg <- list(Scale = .escale , shape = .eshape )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(function(mu, y, w, 
                  residuals = FALSE, eta, extra = NULL) {

    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    
    if (residuals) stop("loglikelihood residuals",
                        "not implemented yet") else {
      sum(c(w) * dexpgeom(x = y, scale = Scale, shape = shape,
                          log = TRUE))
    }
  }, list( .lscale = lscale , .lshape = lshape , 
           .escale = escale , .eshape = eshape ))), 
      
  vfamily = c("expgeometric"),

  deriv = eval(substitute(expression({
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

     temp2 <- exp(-y / Scale)
     temp3 <- shape * temp2
     temp4 <- y / Scale^2
     dl.dscale <-  -1 / Scale + temp4 + 2 * temp4 * temp3 / (1 - temp3)
     dl.dshape <- -1 / (1 - shape)    + 2 * temp2 / (1 - temp3)

    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )            
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    dthetas.detas <- cbind(dscale.deta, dshape.deta)

    answer <- c(w) * cbind(dl.dscale, dl.dshape) * dthetas.detas
    answer
  }), list( .lscale = lscale , .lshape = lshape,
            .escale = escale,  .eshape = eshape ))),

  weight = eval(substitute(expression({
  








        run.varcov <- 0
        ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)

        if (length( .nsimEIM )) {
            for(ii in 1:( .nsimEIM )) {
                ysim <- rexpgeom(n, scale=Scale, shape=shape)

                temp2 <- exp(-ysim / Scale)
                temp3 <- shape * temp2
                temp4 <- ysim / Scale^2
                dl.dscale <-  -1 / Scale + temp4 + 
                             2 * temp4 * temp3 / (1 - temp3)
                dl.dshape <- -1 / (1 - shape) + 
                             2 * temp2 / (1 - temp3)

                temp6 <- cbind(dl.dscale, dl.dshape)
                run.varcov <- run.varcov +
                    temp6[,ind1$row.index] * temp6[,ind1$col.index]
            }

            run.varcov <- run.varcov / .nsimEIM

            wz <- if (intercept.only)
                matrix(colMeans(run.varcov),
                       n, ncol(run.varcov), byrow = TRUE) else run.varcov

            wz <- wz * dthetas.detas[, ind1$row] *
                      dthetas.detas[, ind1$col]
        }

    c(w) * wz      
  }), list( .nsimEIM = nsimEIM ))))
}









dexplog <- function(x, scale = 1, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
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


pexplog <- function(q, scale = 1, shape) {
  ans <- 1 - log1p(-(1-shape) * exp(-q / scale)) / log(shape)
  ans[q <= 0] <- 0
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}



qexplog <- function(p, scale = 1, shape) {


  ans <- -scale * (log1p(-shape^(1.0 - p)) - log1p(-shape))

  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans
}



rexplog <- function(n, scale = 1, shape) {
  ans <- qexplog(runif(n), scale = scale, shape = shape)
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}








explogarithmic.control <- function(save.weight = TRUE, ...) {
    list(save.weight = save.weight)
}

 explogarithmic <- function(lscale = "loge", lshape = "logit",
                            iscale = NULL,   ishape = NULL,
                            tol12 = 1.0e-05, zero = 1,
                            nsimEIM = 400) {

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE) ||
        any(ishape >= 1))
      stop("bad input for argument 'ishape'")

  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE,
                  positive = TRUE))
    stop("bad input for argument 'zero'")


  if (!is.Numeric(nsimEIM, allowable.length = 1,
                  integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
      stop("argument 'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Exponential logarithmic distribution\n\n",
            "Links:    ",
            namesof("Scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     ", "(-polylog(2, 1 - p) * Scale) / log(shape)"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <- c(
      namesof("Scale", .lscale , earg = .escale , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {

      scale.init <- if (is.Numeric( .iscale , positive = TRUE)) {
                     rep( .iscale , len = n)
                   } else {
                     stats::sd(c(y))  
                   }

      shape.init <- if (is.Numeric( .ishape , positive = TRUE)) {
                     rep( .ishape , len = n)
                   } else {
                      rep((exp(median(y)/scale.init) - 1)^2, len = n)
                   }
      shape.init[shape.init >= 0.95] <- 0.95
      shape.init[shape.init <= 0.05] <- 0.05


      etastart <-
        cbind(theta2eta(scale.init, .lscale , earg = .escale ),
              theta2eta(shape.init, .lshape , earg = .eshape ))

    }
   }), list( .lscale = lscale, .lshape = lshape,
             .iscale = iscale, .ishape = ishape,
             .escale = escale, .eshape = eshape))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )



    qexplog(p = 0.5, shape = shape, scale = Scale)  

  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link <-    c(Scale = .lscale , shape = .lshape )

    misc$earg <- list(Scale = .escale , shape = .eshape )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(function(mu, y, w,
                  residuals = FALSE, eta, extra = NULL) {

    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    

    if (residuals) stop("loglikelihood residuals",
                        "not implemented yet") else {
      sum(c(w) * dexplog(x = y, scale = Scale,
                         shape = shape, log = TRUE))
    }
  }, list( .lscale = lscale , .lshape = lshape ,
           .escale = escale , .eshape = eshape ))),

  vfamily = c("explogarithmic"),

  deriv = eval(substitute(expression({
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

     temp2 <- exp(-y / Scale)
     temp3 <- y / Scale^2
     temp4 <- 1 - shape
     dl.dscale <- (-1 / Scale) + temp3 + (temp4 * temp3 *
                  temp2) / (1 - temp4 * temp2)
     dl.dshape <- -1 / (shape * log(shape)) - 1 / temp4 -
                  temp2 / (1 - temp4 * temp2)

    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    dthetas.detas <- cbind(dscale.deta, dshape.deta)

    answer <- c(w) * cbind(dl.dscale, dl.dshape) * dthetas.detas
    answer
  }), list( .lscale = lscale , .lshape = lshape,
            .escale = escale,  .eshape = eshape ))),

  weight = eval(substitute(expression({



        run.varcov <- 0
        ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)

        if (length( .nsimEIM )) {
            for(ii in 1:( .nsimEIM )) {
                ysim <- rexplog(n, scale=Scale, shape=shape)

                temp2 <- exp(-ysim / Scale)
                temp3 <- ysim / Scale^2
                temp4 <- 1 - shape
                dl.dscale <- (-1 / Scale) + temp3 + (temp4 * temp3 *
                             temp2) / (1 - temp4 * temp2)
                dl.dshape <- -1 / (shape * log(shape)) - 1 / temp4 -
                             temp2 / (1 - temp4 * temp2)

                temp6 <- cbind(dl.dscale, dl.dshape)
                run.varcov <- run.varcov +
                           temp6[,ind1$row.index] *
                           temp6[,ind1$col.index]
            }

            run.varcov <- run.varcov / .nsimEIM

            wz <- if (intercept.only)
                matrix(colMeans(run.varcov),
                       n, ncol(run.varcov), byrow = TRUE) else run.varcov

            wz <- wz * dthetas.detas[, ind1$row] *
                      dthetas.detas[, ind1$col]
        }

    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))
}











  
dweibull3 <- function(x, location = 0, scale = 1, shape,
                      log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  dweibull(x = x - location, shape = shape,
           scale = scale, log = log.arg)
}

pweibull3 <- function(q, location = 0, scale = 1, shape) {
  pweibull(q = q - location, scale = scale, shape = shape)
}


qweibull3 <- function(p, location = 0, scale = 1, shape) {
  location + qweibull(p = p, shape = shape, scale = scale)
}


rweibull3 <- function(n, location = 0, scale = 1, shape) {
  location + rweibull(n = n, shape = shape, scale = scale)
}









   ### Two-piece normal (TPN) family 


dtpn <- function(x, location = 0, scale = 1, skewpar = 0.5,
                 log.arg = FALSE) {


  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
           na.rm = TRUE))
    stop("some parameters out of bound")

  LLL <- max(length(x), length(location), length(scale),
            length(skewpar))
  if (length(x) != LLL) x <- rep(x, length = LLL)
  if (length(location) != LLL) location <- rep(location, length = LLL)
  if (length(scale) != LLL) scale <- rep(scale, length = LLL)
  if (length(skewpar) != LLL) skewpar <- rep(skewpar, length = LLL)
    
  zedd <- (x - location) / scale

  log.s1 <-  -zedd^2 / (8 * skewpar^2)
  log.s2 <-  -zedd^2 / (8 * (1 - skewpar)^2)
            
  logdensity <- log.s1
  logdensity[zedd > 0] <- log.s2[zedd > 0]
  
  logdensity <- logdensity -log(scale) - log(sqrt(2 * pi))

  if (log.arg) logdensity else exp(logdensity)
}

ptpn <- function(q, location = 0, scale = 1, skewpar = 0.5) {

  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
          na.rm = TRUE))
    stop("some parameters out of bound")


 zedd <- (q - location) / scale

  s1 <- 2 * skewpar * pnorm(zedd, sd = 2 * skewpar) #/ scale
  s2 <- skewpar + (1 - skewpar) *
        pgamma(zedd^2 / (8 * (1-skewpar)^2), 0.5)
 
ans <- rep(0.0, length(zedd))
ans[zedd <= 0] <- s1[zedd <= 0]
ans[zedd > 0] <- s2[zedd > 0]

ans
}



pos <- function(x) ifelse(x > 0, x, 0.0)
 

qtpn <- function(p, location = 0, scale = 1, skewpar = 0.5) {

  pp = p
  if (any(pp      <= 0 |
          pp      >= 1 |
          skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
             na.rm = TRUE))
    stop("some parameters out of bound")
    # Recycle the vectors to equal lengths
  LLL <- max(length(pp), length(location), length(scale),
            length(skewpar))
  if (length(pp) != LLL) pp <- rep(pp, length = LLL)
  if (length(location) != LLL) location <- rep(location, length = LLL)
  if (length(scale) != LLL) scale <- rep(scale, length = LLL)
  if (length(skewpar) != LLL) skewpar <- rep(skewpar, length = LLL)
       
  qtpn <- rep(as.numeric(NA), length(LLL))
  qtpn <- qnorm(pp / (2 * skewpar), sd = 2 * skewpar)
  qtpn[pp > skewpar] <- sqrt(8 * ( 1 - skewpar)^2 * 
                        qgamma(pos( pp - skewpar) / ( 
                        1 - skewpar),.5))[pp > skewpar]
        
   qtpn * scale + location
  
}





rtpn <- function(n, location = 0, scale = 1, skewpar = 0.5) {


  qtpn(p = runif(n), location = location,
       scale = scale, skewpar = skewpar)
}





tpnff <- function(llocation = "identity", lscale = "loge",
                  pp = 0.5, method.init = 1,  zero = 2)
{
  if (!is.Numeric(method.init, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      method.init > 4)
     stop("argument 'imethod' must be 1 or 2 or 3 or 4")

  if (!is.Numeric(pp, allowable.length = 1, positive = TRUE))
    stop("bad input for argument 'pp'")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
     stop("bad input for argument 'zero'")



  new("vglmff",
  blurb = c("Two-piece normal distribution \n\n",
            "Links: ",
            namesof("location",  llocat,  earg = elocat), ", ",
            namesof("scale",     lscale,  earg = escale), "\n\n",
            "Mean: "),
  constraints = eval(substitute(expression({
          constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <-
       c(namesof("location", .llocat, earg = .elocat, tag = FALSE),
         namesof("scale",    .lscale, earg = .escale, tag = FALSE))




    if (!length(etastart)) {
        junk <- lm.wfit(x = x, y = c(y), w = c(w))
        scale.y.est <-
          sqrt( sum(c(w) * junk$resid^2) / junk$df.residual )
        location.init <- if ( .llocat == "loge")
          pmax(1/1024, y) else {

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
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .method.init=method.init ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat, earg = .elocat)
  }, list( .llocat = llocat,
           .elocat = elocat, .escale = escale ))),
  last = eval(substitute(expression({
    misc$link     <-    c("location" = .llocat, "scale" = .lscale)

    misc$earg     <- list("location" = .elocat, "scale" = .escale)

    misc$expected <- TRUE
    misc$pp       <- .pp
    misc$method.init <- .method.init
    misc$multipleResponses <- FALSE
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .pp     = pp,        .method.init = method.init ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    location <- eta2theta(eta[, 1], .llocat, earg = .elocat)
    myscale  <- eta2theta(eta[, 2], .lscale, earg = .escale)
    ppay     <- .pp
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dtpn(y, skewpar = ppay, location = location,
                      scale = myscale, log.arg = TRUE))
    }
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale,
           .pp      = pp ))),
  vfamily = c("tpnff"),
  deriv = eval(substitute(expression({
    mylocat <- eta2theta(eta[, 1], .llocat,  earg = .elocat)
    myscale <- eta2theta(eta[, 2], .lscale,  earg = .escale)
    mypp    <- .pp

    zedd <- (y - mylocat) / myscale
 #   cond1 <-    (zedd <= 0)
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

    ans <- c(w) * cbind(dl.dlocat * dlocat.deta,
                        dl.dscale * dscale.deta)
    ans
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .pp      = pp ))),
  weight = eval(substitute(expression({
    wz   <- matrix(as.numeric(NA), n, M) # diag matrix; y is one-col too
    temp10 <- mypp * (1 - mypp)
    ned2l.dlocat2        <- 1 / ((4 * temp10) * myscale^2)
    ned2l.dscale2        <- 2 /  myscale^2
     

    wz[, iam(1, 1,M)] <- ned2l.dlocat2 * dlocat.deta^2
    wz[, iam(2, 2,M)] <- ned2l.dscale2 * dscale.deta^2
  # wz[, iam(3, 3,M)] <- ned2l.dskewpar2 * dskewpa.deta^2
  # wz[, iam(1, 3,M)] <- ned2l.dlocatdskewpar * dskewpar.deta * dlocat.deta
      ans
    c(w) * wz
  }))))
}



  ########################################################################


tpnff3 <- function(llocation = "identity",
                    lscale   = "loge",
                    lskewpar = "identity",
                    method.init = 1,  zero = 2)
{
  if (!is.Numeric(method.init, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      method.init > 4)
    stop("argument 'imethod' must be 1 or 2 or 3 or 4")



  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lskewp <- as.list(substitute(lskewpar))
  eskewp <- link2list(lskewp)
  lskewp <- attr(eskewp, "function.name")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")



  new("vglmff",
  blurb = c("Two-piece normal distribution \n\n",
            "Links: ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale",    lscale, earg = escale),  ", ",
            namesof("skewpar",  lscale, earg = eskewp),  "\n\n",
            "Mean: "),
  constraints = eval(substitute(expression({
          constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
       c(namesof("location", .llocat, earg = .elocat, tag = FALSE),
         namesof("scale",    .lscale, earg = .escale, tag = FALSE),
         namesof("skewpar",  .lskewp, earg = .eskewp, tag = FALSE))

    if (!length(etastart)) {
      junk = lm.wfit(x = x, y = c(y), w = c(w))
      scale.y.est <- sqrt(sum(c(w) * junk$resid^2) / junk$df.residual)
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
           theta2eta(skew.l.in,     .lskewp, earg = .escale))
    }
  }), list( .llocat = llocat, .lscale = lscale, .lskewp = lskewp,
            .elocat = elocat, .escale = escale, .eskewp = eskewp,
            
            .method.init=method.init ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat, earg = .elocat)
  }, list( .llocat = llocat,
           .elocat = elocat, .escale = escale ))),
  last = eval(substitute(expression({
    misc$link     <-     c("location" = .llocat,
                           "scale"    = .lscale, 
                           "skewpar"  = .lskewp)

    misc$earg     <-  list("location" = .elocat,
                           "scale"    = .escale,
                           "skewpar"  = .eskewp)

    misc$expected <- TRUE
         misc$method.init <- .method.init
  }), list( .llocat = llocat, .lscale = lscale, .lskewp = lskewp,
            .elocat = elocat, .escale = escale, .eskewp = eskewp,
                    .method.init = method.init ))),
 loglikelihood = eval(substitute(
   function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
   locat <- eta2theta(eta[, 1], .llocat, earg = .elocat)
   myscale  <- eta2theta(eta[, 2], .lscale, earg = .escale)
   myskew   <- eta2theta(eta[, 3], .lskewp, earg = .eskewp)

    if (residuals) stop("loglikelihood residuals not ",
                       "implemented yet") else {
     sum(c(w) * dtpn(y, location = locat,  scale = myscale,
                     skewpar = myskew, log.arg = TRUE))
   }
 }, list( .llocat = llocat, .lscale = lscale, .lskewp = lskewp,
          .elocat = elocat, .escale = escale, .eskewp = eskewp
           ))),
  vfamily = c("tpnff3"),
  deriv = eval(substitute(expression({
    mylocat <- eta2theta(eta[, 1], .llocat,   earg = .elocat)
    myscale <- eta2theta(eta[, 2], .lscale,   earg = .escale)
    myskew  <- eta2theta(eta[, 3], .lskewp, earg = .eskewp)
  

    zedd <- (y - mylocat) / myscale
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
    dskewpar.deta <- dtheta.deta(myskew, .lskewp, earg = .eskewp)
    ans <-
    c(w) * cbind(dl.dlocat * dlocat.deta,
              dl.dscale * dscale.deta,
              dl.dskewpar * dskewpar.deta
              )
    ans
  }), list( .llocat = llocat, .lscale = lscale, .lskewp = lskewp,
            .elocat = elocat, .escale = escale, .eskewp = eskewp
            ))),
  weight = eval(substitute(expression({
    wz <- matrix(as.numeric(NA), n, dimm(M)) # diag matrix; y is one-col too
   
    temp10 <- myskew * (1 - myskew)

    ned2l.dlocat2        <- 1 / ((4 * temp10) * myscale^2)
    ned2l.dscale2        <- 2 /  myscale^2
    ned2l.dskewpar2      <- 3 / temp10
    ned2l.dlocatdskewpar <- (-2 * sqrt(2)) / (temp10 * sqrt(pi) *
                             myscale)
     
    wz[, iam(1, 1,M)] <- ned2l.dlocat2 * dlocat.deta^2
    wz[, iam(2, 2,M)] <- ned2l.dscale2 * dscale.deta^2
    wz[, iam(3, 3,M)] <- ned2l.dskewpar2 * dskewpar.deta^2
    wz[, iam(1, 3,M)] <- ned2l.dlocatdskewpar * dskewpar.deta *
                         dlocat.deta
  
    ans
    c(w) * wz
  }))))
}







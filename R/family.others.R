# These functions are
# Copyright (C) 1998-2021 T.W. Yee, University of Auckland.
# All rights reserved.


















dexppois <- function(x, rate = 1, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(shape), length(rate))
  if (length(x)      != N) x      <- rep_len(x,     N)
  if (length(shape)  != N) shape  <- rep_len(shape, N)
  if (length(rate)   != N) rate   <- rep_len(rate,  N)
  logdensity <- rep_len(log(0), N)

  xok <- (0 < x)

  logdensity[xok] <- log(shape[xok]) + log(rate[xok]) -
                     log1p(-exp(-shape[xok])) - shape[xok] -
                     rate[xok] * x[xok] + shape[xok] *
                     exp(-rate[xok] * x[xok])

  logdensity[shape <= 0] <- NaN
  logdensity[rate <= 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}






qexppois<- function(p, rate = 1, shape,
                    lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- -log(log(exp(ln.p) *
                      (-expm1(shape)) + exp(shape)) / shape) / rate
      ans[ln.p > 0] <- NaN
    } else {
      ans <- -log(log(p * (-expm1(shape)) + exp(shape)) / shape) / rate
      ans[p < 0] <- NaN
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- -log(log(expm1(ln.p) *
                      expm1(shape) + exp(shape)) / shape) / rate
      ans[ln.p > 0] <- NaN
    } else {
      ans <- -log(log(p * expm1(shape) + 1) / shape) / rate
      ans[p < 0] <- NaN
      ans[p > 1] <- NaN
    }
  }
  ans[(shape <= 0) | (rate <= 0)] <- NaN
  ans
}














pexppois<- function(q, rate = 1, shape,
                    lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ans <- log((exp(shape * exp(-rate * q)) -
                    exp(shape)) / -expm1(shape))
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- (exp(shape * exp(-rate * q)) - exp(shape)) / (-expm1(shape))
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- log(expm1(shape * exp(-rate * q)) / expm1(shape))
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      ans <- expm1(shape * exp(-rate * q)) / expm1(shape)
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  }
  ans[(shape <= 0) | (rate <= 0)] <- NaN
  ans
}




rexppois <- function(n, rate = 1, shape) {
  ans <- -log(log(runif(n) * (-expm1(shape)) +
         exp(shape)) / shape) / rate
  ans[(shape <= 0) | (rate <= 0)] <- NaN
  ans
}









 exppoisson <- function(lrate = "loglink", lshape = "loglink",
                        irate = 2.0, ishape = 1.1,
                        zero = NULL) {

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lratee <- as.list(substitute(lrate))
  eratee <- link2list(lratee)
  lratee <- attr(eratee, "function.name")


  iratee <- irate




  if (length(ishape) &&
      !is.Numeric(ishape, positive = TRUE))
    stop("bad input for argument 'ishape'")
  if (length(iratee) &&
      !is.Numeric(iratee, positive = TRUE))
    stop("bad input for argument 'irate'")

  ishape[abs(ishape - 1) < 0.01] = 1.1


  new("vglmff",
  blurb = c("Exponential Poisson distribution \n \n",
            "Links:    ",
            namesof("rate",  lratee, earg = eratee), ", ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     shape/(expm1(shape) * rate)) * ",
                      "genhypergeo(c(1, 1), c(2, 2), shape)"),


  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("rate", "shape"),
         lrate  = .lratee ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .zero = zero, .lratee = lratee, .lshape = lshape ))),


  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <- c(
      namesof("rate",  .lratee , earg = .eratee , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {
      ratee.init <- if (length( .iratee ))
              rep_len( .iratee , n) else
              stop("Need to input a value into argument 'iratee'")
      shape.init <- if (length( .ishape ))
                      rep_len( .ishape , n) else
                      (1/ratee.init - mean(y)) / ((y *
                      exp(-ratee.init * y))/n)


      ratee.init <- rep_len(weighted.mean(ratee.init, w = w), n)

      etastart <-
        cbind(theta2eta(ratee.init, .lratee , earg = .eratee ),
              theta2eta(shape.init, .lshape , earg = .eshape ))


    }
  }), list( .lshape = lshape, .lratee = lratee,
            .ishape = ishape, .iratee = iratee,
            .eshape = eshape, .eratee = eratee))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    ratee <- eta2theta(eta[, 1], .lratee , earg = .eratee )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )





    qexppois(p = 0.5, rate = ratee, shape = shape)
  }, list( .lshape = lshape, .lratee = lratee,
           .eshape = eshape, .eratee = eratee))),

  last = eval(substitute(expression({
    misc$link <-    c( rate = .lratee , shape = .lshape )
    misc$earg <- list( rate = .eratee , shape = .eshape )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lshape = lshape, .lratee = lratee,
            .eshape = eshape, .eratee = eratee))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    ratee <- eta2theta(eta[, 1], .lratee , earg = .eratee )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dexppois(x = y, shape = shape, rate = ratee,
                                 log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lratee = lratee , .lshape = lshape ,
           .eshape = eshape , .eratee = eratee ))),

  vfamily = c("exppoisson"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    ratee <- eta2theta(eta[, 1], .lratee , earg = .eratee )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    okay1 <- all(is.finite(ratee)) && all(0 < ratee) &&
             all(is.finite(shape)) && all(0 < shape)
    okay1
  }, list( .lratee = lratee , .lshape = lshape ,
           .eshape = eshape , .eratee = eratee ))),

  deriv = eval(substitute(expression({
    ratee <- eta2theta(eta[, 1], .lratee , earg = .eratee )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

    dl.dratee <- 1/ratee - y - y * shape * exp(-ratee * y)
    dl.dshape <- 1/shape - 1/expm1(shape) - 1 + exp(-ratee * y)
    dratee.deta <- dtheta.deta(ratee, .lratee , earg = .eratee )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * cbind(dl.dratee * dratee.deta,
                 dl.dshape * dshape.deta)
  }), list( .lshape = lshape, .lratee = lratee,
            .eshape = eshape, .eratee = eratee ))),

  weight = eval(substitute(expression({

    temp1 <- -expm1(-shape)

    ned2l.dshape2 <- (1 + exp(2 * shape) - shape^2 * exp(shape) - 2 *
                      exp(shape)) / (shape * temp1)^2


    ned2l.dratee2 <- 1 / ratee^2 - (shape^2 * exp(-shape) / (4 *
                    ratee^2 * temp1)) *
                    genhypergeo(c(2, 2, 2), c(3, 3, 3), shape)

    ned2l.drateeshape <- (shape * exp(-shape) / (4 * ratee * temp1)) *
                           genhypergeo(c(2, 2), c(3, 3), shape)

    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M)] <- dratee.deta^2 * ned2l.dratee2
    wz[, iam(1, 2, M)] <- dratee.deta * dshape.deta * ned2l.drateeshape
    wz[, iam(2, 2, M)] <- dshape.deta^2 * ned2l.dshape2
    c(w) * wz
  }), list( .zero = zero ))))
}










dgenray <- function(x, scale = 1, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  N <- max(length(x), length(shape), length(scale))
  if (length(x)      != N) x      <- rep_len(x,     N)
  if (length(shape)  != N) shape  <- rep_len(shape, N)
  if (length(scale)  != N) scale  <- rep_len(scale, N)
  logdensity <- rep_len(log(0), N)

  if (any(xok <- (x > 0))) {
    temp1 <- x[xok] / scale[xok]
    logdensity[xok] <- log(2) + log(shape[xok]) + log(x[xok]) -
                       2 * log(scale[xok]) - temp1^2  +
                       (shape[xok] - 1) * log1p(-exp(-temp1^2))
  }
  logdensity[(shape <= 0) | (scale <= 0)] <- NaN
  logdensity[is.infinite(x)] <- log(0)  # 20141209 KaiH
  if (log.arg) {
    logdensity
  } else {
     exp(logdensity)
  }
}






pgenray <- function(q, scale = 1, shape,
                    lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ans <- log((-expm1(-(q/scale)^2))^shape)
      ans[q <= 0 ] <- -Inf
    } else {
      ans <- (-expm1(-(q/scale)^2))^shape
      ans[q <= 0] <- 0
    }
  } else {
    if (log.p) {
      ans <- log(-expm1(shape*log(-expm1(-(q/scale)^2))))
      ans[q <= 0] <- 0
    } else {
      ans <- -expm1(shape*log(-expm1(-(q/scale)^2)))
      ans[q <= 0] <- 1
    }
  }
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}









qgenray <- function(p, scale = 1, shape,
                    lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- scale * sqrt(-log1p(-(exp(ln.p)^(1/shape))))
      ans[ln.p > 0] <- NaN
    } else {
      ans <- scale * sqrt(-log1p(-(p^(1/shape))))
      ans[p < 0] <- NaN
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- scale * sqrt(-log1p(-((-expm1(ln.p))^(1/shape))))
      ans[ln.p > 0] <- NaN
    } else {
      ans <- scale * sqrt(-log1p(-exp((1/shape)*log1p(-p))))
      ans[p < 0] <- NaN
      ans[p > 1] <- NaN
    }
  }
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}






rgenray <- function(n, scale = 1, shape) {
  ans <- qgenray(runif(n), shape = shape, scale = scale)
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}




genrayleigh.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}


 genrayleigh <-
  function(lscale = "loglink", lshape = "loglink",
           iscale = NULL,   ishape = NULL,
           tol12 = 1.0e-05,
           nsimEIM = 300, zero = 2) {

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

  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
      stop("argument 'nsimEIM' should be an integer greater than 50")



  new("vglmff",
  blurb = c("Generalized Rayleigh distribution\n",
            "Links:    ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("scale", "shape"),
         nsimEIM = .nsimEIM ,
         lscale = .lscale ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .zero = zero, .lscale = lscale, .lshape = lshape,
           .nsimEIM = nsimEIM ))),


  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y





    predictors.names <- c(
      namesof("scale", .lscale , earg = .escale , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

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
                    grid.search(scale.grid, objfun = genrayleigh.Loglikfun,
                                y = y, x = x, w = w)
      scale.init <- rep_len(scale.init, length(y))

      shape.init <- if (length( .ishape )) .ishape else
                    -1 / weighted.mean(log1p(-exp(-(y/scale.init)^2)),
                     w = w)
      shape.init <- rep_len(shape.init, length(y))

      etastart <- cbind(theta2eta(scale.init, .lscale , earg = .escale ),
                        theta2eta(shape.init, .lshape , earg = .eshape ))

        }
    }), list( .lscale = lscale, .lshape = lshape,
              .iscale = iscale, .ishape = ishape,
              .escale = escale, .eshape = eshape))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    qgenray(p = 0.5, shape = shape, scale = Scale)
  }, list( .lshape = lshape, .lscale = lscale,
           .eshape = eshape, .escale = escale ))),

  last = eval(substitute(expression({
    misc$link <-    c(scale = .lscale, shape = .lshape )

    misc$earg <- list(scale = .escale, shape = .eshape )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lshape = lshape, .lscale = lscale,
            .eshape = eshape, .escale = escale,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {

    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dgenray(x = y, shape = shape,
                                scale = Scale, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape , .lscale = lscale ,
           .eshape = eshape , .escale = escale ))),

  vfamily = c("genrayleigh"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    okay1 <- all(is.finite(Scale)) && all(0 < Scale) &&
             all(is.finite(shape)) && all(0 < shape)
    okay1
  }, list( .lshape = lshape , .lscale = lscale ,
           .eshape = eshape , .escale = escale ))),

  deriv = eval(substitute(expression({
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    dthetas.detas <- cbind(dscale.deta, dshape.deta)

    temp1 <- y / Scale
    temp2 <- exp(-temp1^2)
    temp3 <- temp1^2 / Scale
    AAA   <- 2 * temp1^2 / Scale  # 2 * y^2 / Scale^3
    BBB   <- -expm1(-temp1^2)     # denominator
    dl.dshape <- 1/shape + log1p(-temp2)
    dl.dscale <- -2 / Scale + AAA * (1 - (shape - 1) * temp2 / BBB)

    dl.dshape[!is.finite(dl.dshape)] =
      max(dl.dshape[is.finite(dl.dshape)])

    answer <- c(w) * cbind(dl.dscale, dl.dshape) * dthetas.detas
    answer
  }), list( .lshape = lshape , .lscale = lscale,
            .eshape = eshape,  .escale = escale ))),

  weight = eval(substitute(expression({


    run.varcov <- 0
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    for (ii in 1:( .nsimEIM )) {
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

      temp3 <- cbind(dl.dscale, dl.dshape)
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
  if (length(x)      != N) x      <- rep_len(x,     N)
  if (length(scale)  != N) scale  <- rep_len(scale, N)
  if (length(shape)  != N) shape  <- rep_len(shape, N)
  logdensity <- rep_len(log(0), N)

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






expgeometric.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 expgeometric <- function(lscale = "loglink", lshape = "logitlink",
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



  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
      stop("'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Exponential geometric distribution\n\n",
            "Links:    ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     ", "(shape - 1) * log(1 - ",
            "shape) / (shape / scale)"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("scale", "shape"),
         nsimEIM = .nsimEIM ,
         lscale = .lscale ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .zero = zero, .lscale = lscale, .lshape = lshape,
           .nsimEIM = nsimEIM ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <- c(
      namesof("scale", .lscale , earg = .escale , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {

      scale.init <- if (is.Numeric( .iscale , positive = TRUE)) {
                      rep_len( .iscale , n)
                    } else {
                      stats::sd(c(y))  # The papers scale parameter beta
                    }

      shape.init <- if (is.Numeric( .ishape , positive = TRUE)) {
                      rep_len( .ishape , n)
                    } else {
                      rep_len(2 - exp(median(y)/scale.init), n)
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
    misc$link <-    c(scale = .lscale , shape = .lshape )

    misc$earg <- list(scale = .escale , shape = .eshape )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {

    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dexpgeom(x = y, scale = Scale, shape = shape,
                                 log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape))),

  vfamily = c("expgeometric"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    okay1 <- all(is.finite(Scale)) && all(0 < Scale) &&
             all(is.finite(shape)) && all(0 < shape & shape < 1)
    okay1
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape))),

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
        for (ii in 1:( .nsimEIM )) {
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
  if (length(x)      != N) x      <- rep_len(x,     N)
  if (length(scale)  != N) scale  <- rep_len(scale, N)
  if (length(shape)  != N) shape  <- rep_len(shape, N)

  logdensity <- rep_len(log(0), N)
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









explogff.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}


 explogff <- function(lscale = "loglink", lshape = "logitlink",
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



  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
      stop("argument 'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Exponential logarithmic distribution\n\n",
            "Links:    ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     ", "(-polylog(2, 1 - p) * scale) / log(shape)"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("scale", "shape"),
         nsimEIM = .nsimEIM ,
         lscale = .lscale ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .zero = zero, .lscale = lscale, .lshape = lshape,
           .nsimEIM = nsimEIM ))),


  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <- c(
      namesof("scale", .lscale , earg = .escale , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {

      scale.init <- if (is.Numeric( .iscale , positive = TRUE)) {
                     rep_len( .iscale , n)
                   } else {
                     stats::sd(c(y))
                   }

      shape.init <- if (is.Numeric( .ishape , positive = TRUE)) {
                     rep_len( .ishape , n)
                   } else {
                      rep_len((exp(median(y)/scale.init) - 1)^2, n)
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
    scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )



    qexplog(p = 0.5, shape = shape, scale = scale)

  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link <-    c(scale = .lscale , shape = .lshape )

    misc$earg <- list(scale = .escale , shape = .eshape )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dexplog(x = y, scale = Scale,
                                shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape))),

  vfamily = c("explogff"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    okay1 <- all(is.finite(Scale)) && all(0 < Scale) &&
             all(is.finite(shape)) && all(0 < shape & shape < 1)
    okay1
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape))),

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
          for (ii in 1:( .nsimEIM )) {
            ysim <- rexplog(n, scale = Scale, shape = shape)

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
                       n, ncol(run.varcov), byrow = TRUE) else
                run.varcov

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

  N <- max(length(x), length(location), length(scale),
           length(skewpar))
  if (length(x)        != N) x         <- rep_len(x,        N)
  if (length(scale)    != N) scale     <- rep_len(scale,    N)
  if (length(location) != N) location  <- rep_len(location, N)
  if (length(skewpar)  != N) skewpar   <- rep_len(skewpar,  N)

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

  s1 <- 2 * skewpar * pnorm(zedd, sd = 2 * skewpar)  #/ scale
  s2 <- skewpar + (1 - skewpar) *
        pgamma(zedd^2 / (8 * (1-skewpar)^2), 0.5)

ans <- rep_len(0.0, length(zedd))
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
  if (length(pp)       != LLL) pp       <- rep_len(pp,       LLL)
  if (length(location) != LLL) location <- rep_len(location, LLL)
  if (length(scale)    != LLL) scale    <- rep_len(scale,    LLL)
  if (length(skewpar)  != LLL) skewpar  <- rep_len(skewpar,  LLL)

  qtpn <- rep_len(NA_real_, length(LLL))
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





tpnff <- function(llocation = "identitylink", lscale = "loglink",
                  pp = 0.5, method.init = 1,  zero = 2) {
  if (!is.Numeric(method.init, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      method.init > 4)
     stop("argument 'imethod' must be 1 or 2 or 3 or 4")

  if (!is.Numeric(pp, length.arg = 1, positive = TRUE))
    stop("bad input for argument 'pp'")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  new("vglmff",
  blurb = c("Two-piece normal distribution \n\n",
            "Links: ",
            namesof("location",  llocat,  earg = elocat), ", ",
            namesof("scale",     lscale,  earg = escale), "\n\n",
            "Mean: "),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("location", "scale"),
         llocation = .llocat ,
         lscale    = .lscale ,
         zero = .zero )
  }, list( .zero = zero,
           .llocat = llocat,
           .lscale = lscale ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <-
       c(namesof("location", .llocat , earg = .elocat , tag = FALSE),
         namesof("scale",    .lscale , earg = .escale , tag = FALSE))




    if (!length(etastart)) {
        junk <- lm.wfit(x = x, y = c(y), w = c(w))
        scale.y.est <-
          sqrt( sum(c(w) * junk$resid^2) / junk$df.residual )
        location.init <- if ( .llocat == "loglink")
          pmax(1/1024, y) else {

        if ( .method.init == 3) {
          rep_len(weighted.mean(y, w), n)
        } else if ( .method.init == 2) {
          rep_len(median(rep(y, w)), n)
        } else if ( .method.init == 1) {
          junk$fitted
        } else {
          y
        }
      }
      etastart <- cbind(
           theta2eta(location.init,  .llocat , earg = .elocat ),
           theta2eta(scale.y.est,    .lscale , earg = .escale ))
    }
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .method.init = method.init ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat , earg = .elocat )
  }, list( .llocat = llocat,
           .elocat = elocat, .escale = escale ))),
  last = eval(substitute(expression({
    misc$link     <-    c("location" = .llocat , "scale" = .lscale )

    misc$earg     <- list("location" = .elocat , "scale" = .escale )

    misc$expected <- TRUE
    misc$pp       <- .pp
    misc$method.init <- .method.init
    misc$multipleResponses <- FALSE
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .pp     = pp,        .method.init = method.init ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    location <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    myscale  <- eta2theta(eta[, 2], .lscale , earg = .escale )
    ppay     <- .pp
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dtpn(y, skewpar = ppay, location = location,
                             scale = myscale, log.arg = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale,
           .pp      = pp ))),
  vfamily = c("tpnff"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    mylocat <- eta2theta(eta[, 1], .llocat ,  earg = .elocat )
    myscale <- eta2theta(eta[, 2], .lscale ,  earg = .escale )
    okay1 <- all(is.finite(mylocat)) &&
             all(is.finite(myscale)) && all(0 < myscale)
    okay1
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale,
           .pp      = pp ))),
  deriv = eval(substitute(expression({
    mylocat <- eta2theta(eta[, 1], .llocat ,  earg = .elocat )
    myscale <- eta2theta(eta[, 2], .lscale ,  earg = .escale )
    mypp    <- .pp

    zedd <- (y - mylocat) / myscale
    cond2 <-    (zedd > 0)

    dl.dlocat        <-  zedd / (4 * mypp^2)  # cond1
    dl.dlocat[cond2] <- (zedd / (4 * (1 - mypp)^2))[cond2]
    dl.dlocat        <- dl.dlocat / myscale

    dl.dscale        <-  zedd^2 / (4 * mypp^2)
    dl.dscale[cond2] <- (zedd^2 / (4 * (1 - mypp)^2))[cond2]
    dl.dscale        <- (-1 + dl.dscale) / myscale


    dlocat.deta <- dtheta.deta(mylocat, .llocat, earg = .elocat)
    dscale.deta <- dtheta.deta(myscale, .lscale, earg = .escale)

    ans <- c(w) * cbind(dl.dlocat * dlocat.deta,
                        dl.dscale * dscale.deta)
    ans
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .pp      = pp ))),
  weight = eval(substitute(expression({
    wz   <- matrix(0, n, M)  # diag matrix; y is one-col too
    temp10 <- mypp * (1 - mypp)
    ned2l.dlocat2        <- 1 / ((4 * temp10) * myscale^2)
    ned2l.dscale2        <- 2 /  myscale^2


    wz[, iam(1, 1, M)] <- ned2l.dlocat2 * dlocat.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dscale2 * dscale.deta^2
  # wz[, iam(3, 3, M)] <- ned2l.dskewpar2 * dskewpa.deta^2
  # wz[, iam(1, 3, M)] <- ned2l.dlocatdskewpar * dskewpar.deta * dlocat.deta
    c(w) * wz
  }))))
}



  ########################################################################


tpnff3 <- function(llocation = "identitylink",
                    lscale   = "loglink",
                    lskewpar = "identitylink",
                    method.init = 1,  zero = 2)
{
  if (!is.Numeric(method.init, length.arg = 1,
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






  new("vglmff",
  blurb = c("Two-piece normal distribution \n\n",
            "Links: ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale",    lscale, earg = escale),  ", ",
            namesof("skewpar",  lskewp, earg = eskewp),  "\n\n",
            "Mean: "),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("location", "scale", "skewpar"),
         llocation = .llocat ,
         lscale    = .lscale ,
         lskewpar  = .lskewp ,
         zero = .zero )
  }, list( .zero = zero,
           .llocat = llocat,
           .lscale = lscale,
           .lskewp = lskewp ))),


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
      location.init <- if ( .llocat == "loglink") pmax(1/1024, y) else {
        if ( .method.init == 3) {
          rep_len(weighted.mean(y, w), n)
        } else if ( .method.init == 2) {
          rep_len(median(rep(y, w)), n)
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
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

   locat    <- eta2theta(eta[, 1], .llocat , earg = .elocat )
   myscale  <- eta2theta(eta[, 2], .lscale , earg = .escale )
   myskew   <- eta2theta(eta[, 3], .lskewp , earg = .eskewp )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dtpn(y, location = locat,  scale = myscale,
                             skewpar = myskew, log.arg = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
   }
  }, list( .llocat = llocat, .lscale = lscale, .lskewp = lskewp,
           .elocat = elocat, .escale = escale, .eskewp = eskewp
           ))),
  vfamily = c("tpnff3"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    mylocat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    myscale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    myskew  <- eta2theta(eta[, 3], .lskewp , earg = .eskewp )
    okay1 <- all(is.finite(mylocat)) &&
             all(is.finite(myscale)) && all(0 < myscale) &&
             all(is.finite(myskew ))
    okay1
  }, list( .llocat = llocat, .lscale = lscale, .lskewp = lskewp,
           .elocat = elocat, .escale = escale, .eskewp = eskewp
         ))),
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
    wz <- matrix(NA_real_, n, dimm(M))  # diag matrix; y is one-col too

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





dzoabeta <- function(x, shape1, shape2, pobs0 = 0,
                     pobs1 = 0, log = FALSE, tol = .Machine$double.eps) {
  log.arg <- log
  rm(log)
  LLL <- max(length(x), length(shape1),
             length(shape2), length(pobs0), length(pobs1))
  if (LLL != length(x))      x      <- rep_len(x,      LLL)
  if (LLL != length(shape1)) shape1 <- rep_len(shape1, LLL)
  if (LLL != length(shape2)) shape2 <- rep_len(shape2, LLL)
  if (LLL != length(pobs0))  pobs0  <- rep_len(pobs0,  LLL)
  if (LLL != length(pobs1))  pobs1  <- rep_len(pobs1,  LLL)
  ans <- rep_len(NA_real_, LLL)

  k1 <- (pobs0 < -tol | pobs1 < -tol |
    (pobs0 + pobs1) > (1 + tol))
  k4 <- is.na(pobs0) | is.na(pobs1)
  ans[!k4 & !k1] <- dbeta(x[!k4 & !k1],
                          shape1[!k4 & !k1],
                          shape2[!k4 & !k1], log = TRUE) +
                    log1p(-(pobs0[!k4 & !k1] + pobs1[!k4 & !k1]))
  k2 <- x == 0 & pobs0 > 0 & !is.na(x)
  k3 <- x == 1 & pobs1 > 0 & !is.na(x)
  ans[k2 & !k4 & !k1] <- log(pobs0[k2 & !k4 & !k1])
  ans[k3 & !k4 & !k1] <- log(pobs1[k3 & !k4 & !k1])
  if (!log.arg) ans <- exp(ans)
  if (any(k1 & !k4)) {
    ans[k1 & !k4] <- NaN
    warning("NaNs produced")
  }
  ans
}



rzoabeta <- function(n, shape1, shape2, pobs0 = 0, pobs1 = 0,
                     tol = .Machine$double.eps) {
  use.n <- if ((length.n <- length(n)) > 1) {
    length.n
  } else {
    if (!is.Numeric(n, integer.valued = TRUE, length.arg = 1,
                    positive = TRUE)) {
      stop("bad input for argument 'n'")
    } else {
      n
    }
  }
  shape1 <- rep_len(shape1, use.n)
  shape2 <- rep_len(shape2, use.n)
  pobs0  <- rep_len(pobs0,  use.n)
  pobs1  <- rep_len(pobs1,  use.n)
  random.number <- runif(use.n)
  ans <- rep_len(NA_real_, use.n)
  k5 <- (pobs0 < -tol | pobs1 < -tol |
           (pobs0 + pobs1) > (1 + tol))
  k4 <- is.na(pobs0) | is.na(pobs1)
  ans[!k4] <- qzoabeta(random.number[!k4], shape1 = shape1,
                       shape2 = shape2, pobs0 = pobs0,
                       pobs1 = pobs1)
  if (any(k5 & !k4)) {
    ans[k5 & !k4] <- NaN
    warning("NaNs produced")
  }
  ans
}



pzoabeta <- function(q, shape1, shape2, pobs0 = 0, pobs1 = 0,
                     lower.tail = TRUE, log.p = FALSE,
                     tol = .Machine$double.eps) {
  LLL <- max(length(q), length(shape1),
             length(shape2), length(pobs0), length(pobs1))
  if (LLL != length(q))      q      <- rep_len(q,      LLL)
  if (LLL != length(shape1)) shape1 <- rep_len(shape1, LLL)
  if (LLL != length(shape2)) shape2 <- rep_len(shape2, LLL)
  if (LLL != length(pobs0))  pobs0  <- rep_len(pobs0,  LLL)
  if (LLL != length(pobs1))  pobs1  <- rep_len(pobs1,  LLL)
  k3 <- (pobs0 < -tol | pobs1 < -tol |
        (pobs0 + pobs1) > (1 + tol))
  k4 <- is.na(pobs0) | is.na(pobs1)
  ans <- rep_len(NA_real_, LLL)
  ans[!k3 & !k4] <- pbeta(q[!k3 & !k4],
                          shape1[!k3 & !k4],
                          shape2[!k3 & !k4], log.p = TRUE) +
    log1p(-(pobs0[!k3 & !k4] + pobs1[!k3 & !k4]))
  ans <- exp(ans)
  k1 <- q >= 0 & !is.na(q)
  k2 <- q >= 1 & !is.na(q)
  ans[k1 & !k3 & !k4] <- ans[k1 & !k3 & !k4] +
    pobs0[k1 & !k3 & !k4]
  ans[k2 & !k3 & !k4] <- ans[k2 & !k3 & !k4] +
    pobs1[k2 & !k3 & !k4]
  if (!lower.tail & log.p) {
    ans <- log1p(-ans)
  } else {
    if (!lower.tail)
      ans <- 1 - ans
    if (log.p)
      ans <- log(ans)
  }
  if (any(k3 & !k4)) {
    ans[k3 & !k4] <- NaN
    warning("NaNs produced")
  }
  ans
}



qzoabeta <- function(p, shape1, shape2, pobs0 = 0, pobs1 = 0,
                     lower.tail = TRUE, log.p = FALSE,
                     tol = .Machine$double.eps) {
  LLL <- max(length(p), length(shape1),
             length(shape2), length(pobs0), length(pobs1))
  if (LLL != length(p))      p      <- rep_len(p,      LLL)
  if (LLL != length(shape1)) shape1 <- rep_len(shape1, LLL)
  if (LLL != length(shape2)) shape2 <- rep_len(shape2, LLL)
  if (LLL != length(pobs0))  pobs0  <- rep_len(pobs0,  LLL)
  if (LLL != length(pobs1))  pobs1  <- rep_len(pobs1,  LLL)
  k0 <- (pobs0 < -tol | pobs1 < -tol |
        (pobs0 + pobs1) > (1 + tol))
  k4 <- is.na(pobs0) | is.na(pobs1)
  ans <- rep_len(NA_real_, LLL)
  if (!lower.tail & log.p) {
    p <- -expm1(p)
  } else{
    if (!lower.tail)
      p <- 1 - p
    if (log.p) {
      p <- exp(p)
    }
  }
  k1 <- p >= 0 & p <= pobs0 & !is.na(p)
  k2 <- p > pobs0 & p < (1 - pobs1) & !is.na(p)
  k3 <- p >= (1 - pobs1) & p <= 1 & !is.na(p)
  ans[k1 & !k0 & !k4] <- 0
  ans[k2 & !k0 & !k4] <-
    qbeta((p[k2 & !k0 & !k4] -
           pobs0[k2 & !k0 & !k4]) / (1 - pobs0[k2 & !k0 & !k4] -
           pobs1[k2 & !k0 & !k4]),
           shape1 = shape1[k2 & !k0 & !k4],
           shape2 = shape2[k2 & !k0 & !k4])
  ans[k3 & !k0 & !k4] <- 1
  if (any(k0 & !k4)) {
    ans[k3 & !k4] <- NaN
    warning("NaNs produced")
  }
  ans
}





log1mexp <- function(x) {
  if (any(x < 0 & !is.na(x)))
    stop("Inputs need to be non-negative!")
  ifelse(x <= log(2), log(-expm1(-x)), log1p(-exp(-x)))
}


log1pexp <- function(x){

  ifelse(x <= -37, exp(x),
         ifelse(x <= 18, log1p(exp(x)),
                ifelse(x <= 33, x + exp(-x), x)))
}






dzoibetabinom.ab <- function(x, size, shape1, shape2, pstr0 = 0,
                             pstrsize = 0, log = FALSE) {
  log.arg <- log
  rm(log)
  LLL <- max(length(x), length(size), length(shape1),
             length(shape2), length(pstr0), length(pstrsize))
  if (LLL != length(x))        x        <- rep_len(x,        LLL)
  if (LLL != length(size))     size     <- rep_len(size,     LLL)
  if (LLL != length(shape1))   shape1   <- rep_len(shape1,   LLL)
  if (LLL != length(shape2))   shape2   <- rep_len(shape2,   LLL)
  if (LLL != length(pstr0))    pstr0    <- rep_len(pstr0,    LLL)
  if (LLL != length(pstrsize)) pstrsize <- rep_len(pstrsize, LLL)
  ans <- rep_len(NA_real_, LLL)
  k1 <- pstr0 < 0 | pstrsize < 0 |
           (pstr0 + pstrsize) > 1
  k <- is.na(size) | is.na(shape1) | is.na(shape2) |
    is.na(pstr0) | is.na(pstrsize) | is.na(x)
  if (sum(!k & !k1) > 0) {
    ans[!k & !k1] <-
      dbetabinom.ab(x[!k & !k1], size[!k & !k1], shape1[!k & !k1],
                    shape2[!k & !k1], log = TRUE) +
      log1p(-(pstr0[!k & !k1]+pstrsize[!k & !k1]))
    if (!log.arg) ans <- exp(ans)
  }
  k2 <- x == 0 & pstr0 > 0
  k3 <- x == size & pstrsize > 0
  if (sum(k2 & !k & !k1) > 0)
    ans[k2 & !k & !k1] <- pstr0[k2 & !k & !k1] +
      ans[k2 & !k & !k1]
  if (sum(k3 & !k & !k1) > 0)
    ans[k3 & !k & !k1] <- pstrsize[k3 & !k & !k1] +
      ans[k3 & !k & !k1]
  if (any(k1 & !k)) {
    ans[k1 & !k] <- NaN
    warning("NaNs produced")
  }
  ans
}



dzoibetabinom <- function(x, size, prob, rho = 0, pstr0 = 0,
                          pstrsize = 0, log = FALSE) {
  dzoibetabinom.ab(x, size, shape1 = prob * (1 - rho) / rho,
                   shape2 = (1 - prob) * (1 - rho) / rho,
                   pstr0 = pstr0, pstrsize = pstrsize, log = log)
}



rzoibetabinom.ab <- function(n, size, shape1, shape2,
                             pstr0 = 0, pstrsize = 0) {
  use.n <- if ((length.n <- length(n)) > 1) {
    length.n
  } else {
    if (!is.Numeric(n, integer.valued = TRUE, length.arg = 1,
                    positive = TRUE)) {
      stop("bad input for argument 'n'")
    } else {
      n
    }
  }
  size     <- rep_len(size,     use.n)
  shape1   <- rep_len(shape1,   use.n)
  shape2   <- rep_len(shape2,   use.n)
  pstr0    <- rep_len(pstr0,    use.n)
  pstrsize <- rep_len(pstrsize, use.n)
  ans      <- rep_len(NA_real_, use.n)
  k <- is.na(size) | is.na(shape1) | is.na(shape2) |
       is.na(pstr0) | is.na(pstrsize)
  k1 <- pstr0 < 0 | pstrsize < 0 |
    (pstr0 + pstrsize) > 1
  random.number <- runif(use.n)
  k2 <- random.number[!k] < pstr0[!k]
  k3 <- pstr0[!k] <= random.number[!k] &
    random.number[!k] <= (1 - pstrsize[!k])
  k4 <- (1 - pstrsize[!k]) < random.number[!k]
  if (sum(k2 & !k1 & !k) > 0)
    ans[k2 & !k1 & !k] <- 0
  if (sum(k3 & !k1 & !k) > 0)
    ans[k3 & !k1 & !k] <- rbetabinom.ab(sum(k3 & !k1 & !k),
                                        size =  size[k3 & !k1 & !k],
                                        shape1 = shape1[k3 & !k1 & !k],
                                        shape2 = shape2[k3 & !k1 & !k])
  if (sum(k4 & !k1 & !k) > 0)
    ans[k4 & !k1 & !k] <- size[k4 & !k1 & !k]
  ans
}



rzoibetabinom <- function(n, size, prob, rho = 0, pstr0 = 0,
                          pstrsize = 0) {
  rzoibetabinom.ab(n, size, shape1 = prob * (1 - rho) / rho,
                   shape2 = (1 - prob) * (1 - rho) / rho,
                   pstr0 = pstr0,
                   pstrsize = pstrsize)
}



pzoibetabinom.ab <- function(q, size, shape1, shape2, pstr0 = 0,
                             pstrsize = 0, lower.tail = TRUE,
                             log.p = FALSE) {
  LLL <- max(length(q), length(size), length(shape1),
             length(shape2), length(pstr0), length(pstrsize))
  if (LLL != length(q))        q        <- rep_len(q,        LLL)
  if (LLL != length(size))     size     <- rep_len(size,     LLL)
  if (LLL != length(shape1))   shape1   <- rep_len(shape1,   LLL)
  if (LLL != length(shape2))   shape2   <- rep_len(shape2,   LLL)
  if (LLL != length(pstr0))    pstr0    <- rep_len(pstr0,    LLL)
  if (LLL != length(pstrsize)) pstrsize <- rep_len(pstrsize, LLL)
  ans <- rep_len(NA_real_, LLL)
  k <- is.na(size) | is.na(shape1) | is.na(shape2) |
    is.na(pstr0) | is.na(pstrsize) | is.na(q)
  k1 <- pstr0 < 0 | pstrsize < 0 |
    (pstr0 + pstrsize) > 1
  if (sum(!k1 & !k) > 0)
    ans[!k & !k1] <-
      pbetabinom.ab(q[!k & !k1], size[!k & !k1],
                    shape1[!k & !k1], shape2[!k & !k1], log.p = TRUE) +
      log1p(-(pstr0[!k & !k1] + pstrsize[!k & !k1]))
  ans <- exp(ans)
  k2 <- q >= 0
  k3 <- q >= size
  if (sum(k2 & !k1 & !k) > 0)
    ans[k2 & !k & !k1] <- ans[k2 & !k & !k1] +
      pstr0[k2 & !k & !k1]
  if (sum(k3 & !k1 & !k) > 0)
    ans[k3 & !k & !k1] <- ans[k3 & !k & !k1] +
      pstrsize[k3 & !k & !k1]
  if (!lower.tail & log.p) {
    ans <- log1p(-ans)
  } else {
    if (!lower.tail)
      ans <- 1 - ans
    if (log.p)
      ans <- log(ans)
  }
  if (any(!k & k1)) {
    ans[!k & k1] <- NaN
    warning("NaNs produced")
  }
  ans
}


pzoibetabinom <- function(q, size, prob, rho,
                          pstr0 = 0, pstrsize = 0,
                        lower.tail = TRUE, log.p = FALSE) {
  pzoibetabinom.ab(q, size, shape1 = prob * (1 - rho) / rho,
                 shape2 = (1 - prob) * (1 - rho) / rho,
                 pstr0 = pstr0, pstrsize = pstrsize,
                 lower.tail = lower.tail, log.p = log.p)
}














  AR1EIM<- function(x = NULL,
                    var.arg   = NULL,
                    p.drift   = NULL,
                    WNsd      = NULL,
                    ARcoeff1  = NULL,
                    eps.porat = 1e-2) {

    if (!is.matrix(x))
      stop("Argument 'x' must be a matrix.")

    yy   <- x
    M    <- 3
    nn   <- nrow(x)
    nn0  <- numeric(0)
    NOS  <- ncol(x)

    if (!is.matrix(WNsd))
      WNsd <- matrix(WNsd, nrow = nn, ncol = NOS, byrow = TRUE)

    if (!is.matrix(ARcoeff1))
      ARcoeff1 <- matrix(ARcoeff1, nrow = nn, ncol = NOS, byrow = TRUE)

    if (!is.Numeric(eps.porat, length.arg = 1) || eps.porat < 0 ||
        eps.porat > 1e-2)
      stop("Bad input for argument 'eps.porat'.")

    sdTSR   <- colMeans(WNsd)
    sdTSv   <- colMeans(WNsd)
    drift.v <- rep(p.drift, NOS)[1:NOS]
    Aux11   <- (NOS > 1)
    the1v   <- colMeans(ARcoeff1)
    JFin    <- array(0.0, dim = c(nn, NOS, M + (M - 1) + (M - 2) ))

    for (spp in 1:NOS) {

      x <- yy[, spp]
      the1    <- the1v[spp]
      drift.p <- drift.v[spp]
      sdTS    <- sdTSv[spp]

      r <- numeric(nn)
      r <- AR1.gammas(x = x, lags = nn - 1)
      r[nn] <- r[1]

      s0 <- numeric(nn)
      s1 <- numeric(nn)
      s1 <- if (var.arg) (the1^(0:(nn - 1))) / (1 - the1^2) else
             2 * (the1^(0:(nn - 1))) * sdTS / (1 - the1^2)

      s2    <- numeric(nn)
      help1 <- c(0:(nn - 1))
      s2    <- help1 * (the1^(help1 - 1)) * (sdTS^2) / (1 - the1^2) +
                   2 * (sdTS^2) * (the1^(help1 + 1)) / (1 - the1^2)^2
      sMat <- cbind(s0, s1, s2)

      J  <- array(NA_real_,
                  dim = c(length(the1) + 2, length(the1) + 2, nn))
      Jp <- array(NA_real_,
                  dim = c(length(the1) + 2, length(the1) + 2, nn))

      alpha    <- numeric(nn)
      alpha[1] <- 1
      delta    <- r[1]
      eta      <- matrix(NA_real_, nrow = nn, ncol = M)
      eta[1, ] <- cbind(s0[1], s1[1], s2[1])

      psi      <- matrix(0, nrow = nn, ncol = length(the1) + 2)
      psi[1, ] <- cbind(s0[1], s1[1], s2[1]) / r[1]

      u0  <- rep(1/(1 - sign(the1v[1]) * min(0.975, abs(the1v[1]))), nn )
      u1  <- rep(drift.p/(1 - the1)^2, nn)
      uMat <- cbind(u0, rep(0, nn), u1)

      aux1 <- matrix(sMat[1, ],
                     nrow = 2 + length(the1),
                     ncol = 2 + length(the1), byrow = TRUE)
      diag(aux1) <- sMat[1, ]
      J[, , 1]   <- Jp[, , 1] <- aux1 * t(aux1) / (2 * r[1]^2)
      J[1, 1, 1] <- Jp[1, 1, 1] <- 1 / sdTS^2
      JFin[1, spp, 1:M] <- Jp[, , 1][row(Jp[, , 1]) == col(Jp[, , 1])]
      Neps.porat <- 1.819*eps.porat*(1e-10)

      dk   <- matrix(NA_real_, nrow = 1, ncol = length(the1) + 2)
      eR   <- matrix(NA_real_, nrow = 1, ncol = length(the1) + 2)
      cAux2 <- d55 <- numeric(nn); d55[1] <- 0.1

      for (jay in 1:(nn - 1)) {

        cAux <- as.numeric(alpha[1:jay] %*%
                    r[2:(jay + 1)][length(r[2:(jay + 1)]):1])/delta

        dk <- alpha[1:jay] %*%
         sMat[2:(jay + 1), , drop = FALSE][length(sMat[2:(jay + 1)]):1, ]

        delta        <- delta * (1 - cAux^2)
        d55[jay + 1] <- cAux^2

        if ((d55[jay + 1] < eps.porat*1e-2) || (jay > 1e1)) {
          nn0 <- jay
          break
        }

        eta[jay + 1, ] <- dk
        tAux <- numeric(jay + 1)
        tAux <- alpha[1:(jay + 1)] -
                      cAux * alpha[1:(jay + 1)][(jay + 1):1]
        alpha[1:(jay + 1)] <- tAux[1:(jay + 1)]

        eR <- alpha[1:(jay + 1)][(jay + 1):1] %*%
          eta[1:(jay + 1), , drop = FALSE]

        tAux <- eta[1:(jay + 1), ] -
          cAux * eta[1:(jay + 1), ][(jay + 1):1, ]

        eta[1:(jay + 1), ] <- tAux

        AuxE <- matrix(eR, nrow = jay + 1, ncol = M, byrow = TRUE)
        Aux3 <- matrix(alpha[1:(jay + 1)][(jay + 1):1],
                       nrow = jay + 1, ncol = M, byrow = FALSE)
        Aux4 <- matrix(alpha[1:(jay + 1)],
                       nrow = jay + 1, ncol = M, byrow = FALSE)
                  tAux <- psi[1:(jay + 1), ] -
                      cAux * psi[1:(jay + 1), ][(jay + 1):1, ] +
                         AuxE * (Aux3 - cAux * Aux4) / delta

        if (any(dim(psi[1:(jay + 1), ])) != any(dim(tAux)) )
          stop("Invalids 'psi' and 'tAux'.")

        psi[1:(jay + 1), ] <- tAux
        fk <- alpha[1:(jay + 1)] %*% eta[1:(jay + 1), ]
        gk <- alpha[1:(jay + 1)][(jay + 1):1] %*% uMat[1:(jay + 1), ]

        Auxf <- matrix(fk, nrow = M, ncol = M, byrow = FALSE)
        Auxg <- matrix(gk, nrow = M, ncol = M, byrow = FALSE)
        J[, , jay + 1] <-
          J[, , jay] + t(eta[1:(jay + 1), ]) %*% psi[1:(jay + 1), ] /
                      delta - 0.5 * Auxf * t(Auxf) / delta^2 +
                          Auxg * t(Auxg) / delta

        Jp[, , jay + 1] <- J[, , jay + 1] - J[, , jay]
        JFin[jay + 1, spp , 1:M ] <-
          Jp[, , jay + 1][col(Jp[, , jay + 1]) == row(Jp[, , jay + 1])]

        helpC <- numeric(0)
        for (kk in 1:(M - 1))  {
          TF1 <- ( col(Jp[, , jay + 1]) >= row(Jp[, , jay + 1]) )
          TF2 <- (abs(col(Jp[, , jay + 1]) - row(Jp[, , jay + 1])) == kk )
          helpC <- c(helpC, Jp[, , jay + 1][TF1 & TF2])
        }
        rm(TF1, TF2)

        JFin[jay + 1, spp , -(1:M) ] <- helpC
      }

      if (length(nn0))
        for (kk in nn0:(nn - 1)) {
          J[, , kk + 1] <- J[, , nn0] + (kk - nn0 + 1) * Jp[, , nn0]
          Jp[, , kk + 1] <- J[, , kk + 1] - J[, , kk]

          JFin[kk + 1, spp , 1:M ] <-
            Jp[, , kk + 1][col(Jp[, , kk + 1]) == row(Jp[, , kk + 1])]

          helpC <- numeric(0)
          for (ll in 1:(M - 1))  {
           TF1 <- ( col(Jp[, , kk + 1]) >= row(Jp[, , kk + 1]) )
           TF2 <- (abs(col(Jp[, , kk + 1]) - row(Jp[, , kk + 1])) == ll)
            helpC <- c(helpC, Jp[, , kk + 1][TF1 & TF2])
          }
          rm(TF1, TF2)
          JFin[kk + 1, spp , -(1:M) ] <- helpC
        }
      JFin[which(JFin <= Neps.porat)] <-
        abs( JFin[which(JFin <= Neps.porat)])
    }

    JFin

  } # End






AR1.gammas <- function(x, y = NULL, lags = 1) {
  xx  <- matrix(x, ncol = 1)
  nx  <- nrow(xx)

  if (lags < 0 || !(is.Numeric(lags, integer.valued = TRUE)))
    stop("'lags' must be a positive integer.")

  if (length(y)) {
    yy <- matrix(y, ncol = 1)
    ny <- nrow(yy)
    if (nx != ny)
      stop("Number of rows differs.") else
        n <- nx
  } else {
    yy <- xx
    n  <- nrow(xx)
  }

  myD <- numeric(lags + 1)
  myD[1] <- if (length(y)) cov(xx, yy) else cov(xx, xx)  # i.e. var(xx)
  if (lags > 0)
    for (ii in 1:lags)
      myD[ii + 1]  <- cov(xx[-(1:ii), 1], yy[1:(n - ii) , 1])

  myD
}





dzipfmb <- function(x, shape, start = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(start), length(shape))
  if (length(x)     != LLL) x     <- rep_len(x,     LLL)
  if (length(start) != LLL) start <- rep_len(start, LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  
  bad0 <- !is.finite(shape) | !is.finite(start) |
          shape <= 0 | 1 <= shape |
          start != round(start) | start < 1
  bad <- bad0 | !is.finite(x) | x < start | x != round(x)

  logpdf <- x + shape + start
  
  if (any(!bad)) {
    logpdf[!bad] <- lbeta(    x[!bad] - shape[!bad], shape[!bad] + 1) - 
                    lbeta(start[!bad] - shape[!bad], shape[!bad])
  }
  
  logpdf[!bad0 & is.infinite(x)] <- log(0)
  logpdf[!bad0 & x < start     ] <- log(0)
  logpdf[!bad0 & x != round(x) ] <- log(0)
  logpdf[ bad0] <- NaN

  if (log.arg) logpdf else exp(logpdf)
} # dzipfmb()



pzipfmb <-
  function(q, shape, start = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log'")
  
  LLL <- max(length(shape), length(q), length(start))
  if (length(q)     != LLL) q     <- rep_len(q,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(start) != LLL) start <- rep_len(start, LLL)
  
  q.use <- pmax(start, floor(q + 1))
  
  bad0 <- !is.finite(shape) | !is.finite(start) |
          shape <= 0 | 1 <= shape |
          start != round(start) | start < 1
  bad <- bad0 | !is.finite(q.use)

  ans <- q + shape + start

  log.S.short <- lgamma(start[!bad]) + lgamma(q.use[!bad] - shape[!bad]) -
                 lgamma(q.use[!bad]) - lgamma(start[!bad] - shape[!bad])
  
  ans[!bad] <-
  if (lower.tail) {
    if (log.p) { # lower.tail = T, log.p = T
      log1p(-exp(log.S.short))
    } else {  # lower.tail = T, log.p = F
      -expm1(log.S.short)
    }
  } else { 
    if (log.p) {  # lower.tail = F, log.p = T
      log.S.short
    } else {      # lower.tail = F, log.p = F
      exp(log.S.short)
    }
  }

  ans[!bad0 & is.infinite(q.use) & start < q.use] <-
    if (lower.tail) {if (log.p) log(1) else 1} else
    {if (log.p) log(0) else 0}

  ans[bad0] <- NaN
  ans
}  # pzipfmb()






qzipfmb <- function(p, shape, start = 1) {

  LLL <- max(length(p), length(shape), length(start))
  if (length(p)     != LLL) p     <- rep_len(p,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(start) != LLL) start <- rep_len(start, LLL)
  ans <- p + shape + start

  bad0 <- !is.finite(shape) | !is.finite(start) |
          shape <= 0 | 1 <= shape |
          start != round(start) | start < 1
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(start, LLL) - 0.5
  approx.ans <- lo  # True at lhs
  hi <- 2 * lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate | p <= pzipfmb(hi, shape, start = start)
  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 2
  max.iter <- round(log2(1e300)) - 2
  while (!all(done) && iter < max.iter) {
    lo[!done] <- hi[!done]
    hi[!done] <- 2 * hi[!done] + 10.5  # Bug fixed
    done[!done] <- (p[!done] <= pzipfmb(hi[!done],
                    shape = shape[!done], start[!done]))
    iter <- iter + 1
  }

  foo <- function(q, shape, start, p)
    pzipfmb(q, shape = shape, start) - p

  lhs <- dont.iterate |
         p <= dzipfmb(start, shape = shape, start)
  
  approx.ans[!lhs] <-
    bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                    shape = shape[!lhs],
                    start = start[!lhs], p = p[!lhs])
  faa <- floor(approx.ans[!lhs])
  tmp <-
    ifelse(pzipfmb(faa, shape = shape[!lhs], start[!lhs]) < p[!lhs] &
           p[!lhs] <= pzipfmb(faa+1, shape = shape[!lhs], start[!lhs]),
           faa+1, faa)
  ans[!lhs] <- tmp

  vecTF <- !bad0 & !is.na(p) &
           p <= dzipfmb(start, shape = shape, start)
  ans[vecTF] <- start[vecTF]

  ans[!bad0 & !is.na(p) & p == 0] <- start[!bad0 & !is.na(p) & p == 0]
  ans[!bad0 & !is.na(p) & p == 1] <- Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qzipfmb



rzipfmb <- function(n, shape, start = 1) {
  qzipfmb(runif(n), shape, start = start)
}  # rzipfmb












dextlogF <- function(x, lambda, tau,
                   location = 0, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(lambda), length(tau),
             length(location), length(scale))
  if (length(x)        != LLL) x        <- rep_len(x,        LLL)
  if (length(lambda)   != LLL) lambda   <- rep_len(lambda,   LLL)
  if (length(tau)      != LLL) tau      <- rep_len(tau,      LLL)
  if (length(location) != LLL) location <- rep_len(location, LLL)
  if (length(scale)    != LLL) scale    <- rep_len(scale,    LLL)
  
  bad0 <- !is.finite(lambda) | !is.finite(tau) | 
          !is.finite(location) | !is.finite(scale) | 
          lambda <= 0 | tau < 0 | 1 < tau | scale <= 0
  bad <- bad0 | !is.finite(x)

  logpdf <- x + lambda + tau + location + scale
  
  if (any(!bad)) {
    zedd <- (x[!bad] - location[!bad]) / scale[!bad]
    logpdf[!bad] <- (1 - tau[!bad]) * zedd -
      lambda[!bad] * log1p(exp(zedd / lambda[!bad])) -
      lbeta((1 - tau[!bad]) * lambda[!bad],
                 tau[!bad]  * lambda[!bad]) -
      log(lambda[!bad] * scale[!bad])
  }
  
  logpdf[!bad0 & is.infinite(x)] <- log(0)
  logpdf[ bad0] <- NaN

  if (log.arg) logpdf else exp(logpdf)
}  # dextlogf()









extlogF1.control <-
  function(stepsize = 0.5,
           maxit = 100,
           ...) {
  list(stepsize = stepsize,
       maxit    = maxit)
}



 extlogF1 <-
  function(tau = c(0.25, 0.5, 0.75),  # NULL,  # \in (0, 1)
           parallel = TRUE ~ 0,  # FALSE,
           seppar = 0,
           tol0 = -0.001,  # Negative means relative, + means absolute
           llocation = "identitylink",
           ilocation = NULL,
           lambda.arg = NULL,  # 0.1,  # NULL means an adaptive value
           scale.arg = 1,  # Best to leave this alone
           ishrinkage = 0.95,
           digt = 4,
           idf.mu = 3,
           imethod = 1) {


  apply.parint.locat <- FALSE


  if (!is.Numeric(seppar, length.arg = 1) || seppar < 0)
    stop("bad input for argument 'seppar'")
  if (!is.Numeric(tol0, length.arg = 1))  # || tol0 < 0
    stop("bad input for argument 'tol0'")

  if (!is.Numeric(tau, positive = TRUE))
    stop("bad input for argument 'tau'")
  if (any(1 <= tau))
    stop("argument 'tau' must have values in (0, 1)")
  if (length(tau) > 1 && any(diff(tau) <= 0))
    stop("argument 'tau' must be an increasing sequence")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 4)
    stop("argument 'imethod' must be 1, 2 or ... 4")


  llocation <- llocation
  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  ilocat <- ilocation


  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")
  if (!is.Numeric(scale.arg, positive = TRUE))
    stop("bad input for argument 'scale.arg'")




  new("vglmff",
  blurb = c("One-parameter extended log-F distribution\n\n",
            "Links:      ",
            namesof("location", llocat, earg = elocat),
            "\n", "\n",
            "Quantiles:  location(tau)"),

  constraints = eval(substitute(expression({

    onemat <- matrix(1, M, 1)
    constraints.orig <- constraints


    cm1.locat <- diag(M)
    cmk.locat <- onemat
    con.locat <- cm.VGAM(cmk.locat,
                         x = x, bool = .parallel ,
                         constraints = constraints.orig,
                         apply.int   = .apply.parint.locat ,
                         cm.default           = cm1.locat,
                         cm.intercept.default = cm1.locat)
    
    constraints <- con.locat

  }), list( .parallel = parallel, .seppar = seppar, .tol0 = tol0,
            .apply.parint.locat = apply.parint.locat ))),



  infos = eval(substitute(function(...) {
    list(M1 = NA,  # 1,
         Q1 = NA,  # 1,
         tau        = .tau ,
         scale.arg  = .scale.arg ,
         lambda.arg = .lambda.arg ,
         seppar     = .seppar ,
         tol0.arg   = as.vector( .tol0 ),  # Could be negative
         digt       = .digt ,
         expected   = TRUE,
         multipleResponses = TRUE,
         parameters.names  = "location")
  }, list( .lambda.arg = lambda.arg, .scale.arg = scale.arg,
           .digt = digt,
           .tau = tau, .seppar = seppar, .tol0 = tol0 ))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = if (length( .tau ) > 1) 1 else Inf,
              ncol.y.max = if (length( .tau ) > 1) 1 else Inf,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    extra$ncoly <- ncoly <- ncol(y)
    if ((ncoly > 1) && (length( .tau ) > 1 ||
        length( .scale.arg ) > 1))
      stop("response must be a vector if 'tau' or 'scale.arg' ",
           "has a length greater than one")



    if ((midspread <- diff(quantile(y, probs = c(0.25, 0.75)))) == 0)
      stop("could not work out an adaptive 'lambda'") else                
    lambda <- if (is.null( .lambda.arg )) {
    muxfactor <- 25
    extra$midspread <- midspread
         as.vector((midspread / .scale.arg) / muxfactor)
    } else {
      as.vector( .lambda.arg )
    }
    extra$lambda.arg <- lambda
    tol0 <- as.vector( .tol0 )  # A negative value means relative
    if (tol0 < 0)
      tol0 <- as.vector(abs(tol0) * midspread)
    extra$tol0 <- tol0  # A positive value means absolute

    extra$M <- M <- max(length( .scale.arg ),
                        ncoly,
                        length( .tau ))  # Recycle
    extra$scale.arg <- rep_len( .scale.arg , M)
    extra$tau       <- rep_len( .tau       , M)
    extra$n <- n


    extra$tau.names <- tau.names <-
      paste("(tau = ", round(extra$tau, digits = .digt), ")", sep = "")
    extra$Y.names <- Y.names <- if (ncoly > 1) dimnames(y)[[2]] else "y"
    if (is.null(Y.names) || any(Y.names == ""))
      extra$Y.names <- Y.names <- paste("y", 1:ncoly, sep = "")
    extra$y.names <- y.names <-
      if (ncoly > 1) paste(Y.names, tau.names, sep = "") else tau.names

    extra$individual <- FALSE

    mynames1 <- param.names("location", M, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .llocat , earg = .elocat , tag = FALSE))

    seppar <- as.vector( .seppar )
    if (seppar > 0 && M > 1) {
    }  # seppar

    locat.init <- matrix(0, n, M)
    if (!length(etastart)) {

      for (jay in 1:M) {
        y.use <- if (ncoly > 1) y[, jay] else y
        locat.init[, jay] <- if ( .imethod == 1) {
           quantile(y.use, probs = extra$tau[jay])
        } else if ( .imethod == 2) {
          weighted.mean(y.use, w[, min(jay, ncol(w))])
        } else if ( .imethod == 3) {
          median(y.use)
        } else if ( .imethod == 4) {
          Fit5 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                 y = y.use, w = w, df = .idf.mu )
          c(predict(Fit5, x = x[, min(ncol(x), 2)])$y)
        } else {
          use.this <- weighted.mean(y.use, w[, min(jay, ncol(w))])
          (1 - .ishrinkage ) * y.use + .ishrinkage * use.this
        }

        if (length( .ilocat )) {
          locat.init <- matrix( .ilocat , n, M, byrow = TRUE)
        }

        if ( .llocat == "loglink") locat.init <- abs(locat.init)
        etastart <-
          cbind(theta2eta(locat.init, .llocat , earg = .elocat ))
      }
    }
    }), list( .imethod = imethod, .seppar = seppar, .tol0 = tol0,
              .idf.mu = idf.mu, .lambda.arg = lambda.arg,
              .ishrinkage = ishrinkage, .digt = digt,
              .elocat = elocat, .scale.arg = scale.arg,
              .llocat = llocat, .tau = tau,
              .ilocat = ilocat ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    locat <- eta2theta(eta, .llocat , earg = .elocat )
    if (length(locat) > extra$n)
      dimnames(locat) <- list(dimnames(eta)[[1]], extra$y.names)
    locat
  }, list( .elocat = elocat,
           .llocat = llocat, .scale.arg = scale.arg,
           .tau = tau, .lambda.arg = lambda.arg ))),
  last = eval(substitute(expression({
    misc$link <- setNames(rep_len( .llocat , M), mynames1)

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M) {
      misc$earg[[ii]] <- ( .elocat )
    }

    extra$eCDF <- numeric(M)
    locat <- as.matrix(locat)
    for (ii in 1:M) {
      y.use <- if (ncoly > 1) y[, ii] else y
      extra$eCDF[ii] <-
        100 * weighted.mean(y.use <= locat[, ii], w[, min(ii, ncol(w))])
    }
    names(extra$eCDF) <- y.names


    extra$scale.arg <- ( .scale.arg )
    }), list( .elocat = elocat,
              .llocat = llocat, .scale.arg = scale.arg, .tau = tau ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    ymat <- matrix(y, extra$n, extra$M)
    locat <- eta2theta(eta, .llocat , earg = .elocat )
    taumat <- matrix(extra$tau,       extra$n, extra$M, byrow = TRUE)
    Scale  <- matrix(extra$scale.arg, extra$n, extra$M, byrow = TRUE)
    lambda  <- extra$lambda.arg

    Ans <- if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dextlogF(x = c(ymat), location = c(locat),
                                 scale = c(Scale), tau = c(taumat),
                                 lambda = lambda, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }

    n <- nrow(eta)
    M <- NCOL(eta)
    seppar <- as.vector( .seppar )
    tol0 <- extra$tol0  # Absolute
    if (seppar > 0 && M > 1) {
      negdiffmat <- as.matrix(locat[, -ncol(locat)] - locat[, -1])
      negdiffmat <- matrix(pmax(0, tol0 + negdiffmat)^3, n, M-1)
      Adjustment <- seppar * sum(negdiffmat)
      Ans <- Ans - Adjustment
    }  # seppar

    Ans 
  }, list( .elocat = elocat, .scale.arg = scale.arg, .tau = tau,
           .llocat = llocat, .seppar = seppar ))),
  vfamily = c("extlogF1"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    locat <- eta2theta(eta, .llocat , earg = .elocat )
    okay1 <- all(is.finite(locat))
    okay1
  }, list( .elocat = elocat, .scale.arg = scale.arg, .tau = tau,
           .llocat = llocat ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    extra  <- object@extra
    locat  <- eta2theta(eta, .llocat , .elocat )
    Scale  <- matrix(extra$scale.arg, extra$n, extra$M, byrow = TRUE)
    taumat <- matrix(extra$tau,       extra$n, extra$M, byrow = TRUE)
    lambda <- extra$lambda.arg
    relogF(nsim * length(Scale), location = c(locat),
           lambda = lambda, scale = c(Scale), tau = c(taumat))
  }, list( .elocat = elocat, .scale.arg = scale.arg, .tau = tau,
           .llocat = llocat ))),

  deriv = eval(substitute(expression({
    seppar <- as.vector( .seppar )
    tol0 <- extra$tol0  # Absolute
    locat <- eta2theta(eta, .llocat , earg = .elocat )

    ymat <- matrix(y, n, M)
    Scale <- matrix(extra$scale.arg, n, M, byrow = TRUE)

    lambda <- extra$lambda.arg

    taumat <- matrix(extra$tau, n, M, byrow = TRUE)
    zedd <- (ymat - locat) / Scale

    dl.dlocat <- (taumat - 1 + plogis(zedd / lambda)) / Scale
    dlocat.deta <- dtheta.deta(locat, .llocat , earg = .elocat )

    sep.penalize <- seppar > 0 && M > 1
    if (sep.penalize) {
      dl.dlocat[, 1] <- dl.dlocat[, 1] - 3 * seppar *
                        pmax(0, tol0 + locat[, 1] - locat[, 2])^2
      dl.dlocat[, M] <- dl.dlocat[, M] + 3 * seppar *
                        pmax(0, tol0 + locat[, M-1] - locat[, M])^2
      if (M > 2)
      for (jay in 2:(M-1))
        dl.dlocat[, jay] <- dl.dlocat[, jay] - 3 * seppar * (
          pmax(0, tol0 + locat[, jay  ] - locat[, jay+1])^2 -
          pmax(0, tol0 + locat[, jay-1] - locat[, jay  ])^2)
    }  # seppar

    c(w) * cbind(dl.dlocat * dlocat.deta)
  }), list( .elocat = elocat, .seppar = seppar,
            .llocat = llocat, .scale.arg = scale.arg, .tau = tau ))),

  weight = eval(substitute(expression({
    ned2l.dlocat2 <- taumat * (1 - taumat) / ((1 + lambda) * Scale^2)

    if (sep.penalize) {
      ned2l.dlocat2[, 1] <- ned2l.dlocat2[, 1] + 6 * seppar *
                            pmax(0, tol0 + locat[,   1] - locat[, 2])
      ned2l.dlocat2[, M] <- ned2l.dlocat2[, M] + 6 * seppar *
                            pmax(0, tol0 + locat[, M-1] - locat[, M])
      if (M > 2)
        for (jay in 2:(M-1))
          ned2l.dlocat2[, jay] <-
          ned2l.dlocat2[, jay] + 6 * seppar * (
            pmax(0, tol0 + locat[, jay  ] - locat[, jay+1]) +
            pmax(0, tol0 + locat[, jay-1] - locat[, jay  ]))
      rhs.mat <- -6 * seppar * pmax(0, tol0 + locat[, -M] - locat[, -1])
      rhs.mat <- matrix(c(rhs.mat), n, M-1)  # Right dimension
    }  # seppar

    rhs.mat <- if (sep.penalize)
      rhs.mat * dlocat.deta[, -M] * dlocat.deta[, -1] else NULL
    wz <- cbind(ned2l.dlocat2 * dlocat.deta^2,
                rhs.mat)
    c(w) * wz
  }), list( .elocat = elocat, .scale.arg = scale.arg,
            .llocat = llocat ))))
}  # extlogF1()








 setClass("extlogF1",                 contains = "vglmff")



setMethod("showvglmS4VGAM",
          signature(VGAMff = "extlogF1"),
  function(object,
           VGAMff,
           ...) {
  cat("\nQuantiles:\n")
  print(eCDF(as(object, "vglm")))
  invisible(object)
  })




setMethod("showsummaryvglmS4VGAM",  signature(VGAMff = "extlogF1"),
  function(object,
           VGAMff,
           ...) {
  cat("Quantiles:\n")
  print(eCDF(as(object, "vglm"), all = TRUE))
  cat("\n")

  invisible(object)
})








is.crossing.vglm <- function(object, ...) {
  if (!is(object, "vglm"))
    stop("the object is not a VGLM")
  if (any(object@family@vfamily == "lms.bcn"))
    return(FALSE)
  if (!any(object@family@vfamily == "extlogF1"))
    stop("the object does not have a 'extlogF1' family function")
  locat <- fitted(object)
  if (ncol(locat) == 1) {
    TRUE
  } else {
    diffmat <- locat[, -1] - locat[, -ncol(locat)]
    crossq <- any(diffmat < 0)
    crossq
  }
}  # is.crossing.vglm



if (!isGeneric("is.crossing"))
  setGeneric("is.crossing", function(object, ...)
             standardGeneric("is.crossing"),
             package = "VGAM")


setMethod("is.crossing",  "vglm", function(object, ...)
    is.crossing.vglm(object, ...))







 fix.crossing.vglm <-
     function(object, maxit = 100, trace = FALSE,  # etastart = NULL,
              ...) {
  if (!is.crossing(object) ||
      any(object@family@vfamily == "lms.bcn"))
    return(object)

  if (!any(object@family@vfamily == "extlogF1"))
    stop("the object does not have a 'extlogF1' family function")

  object.save <- object
  M <- npred(object.save)
  if (M > 1 &
      !all(is.parallel(object.save)) &
      nrow(coef(object.save, matrix = TRUE)) > 1) {
    if (!has.intercept(object.save))
      stop("the object must have an intercept term")
    for (jay in 1:M) {  # M is an upperbound
      Hlist <- constraints(object, type = "term")
      locat <- fitted(object) 
      diffmat <- locat[,           -1, drop = FALSE] -
                 locat[, -ncol(locat), drop = FALSE]
      crossq <- any(diffmat < 0)
      if (crossq) {
        min.index <- which.min(apply(diffmat, 2, min))
        Hk <- Hlist[[2]]  # Next covariate past the intercept
        use.min.index <- 1 + which(cumsum(colSums(Hk)) == min.index)
        if (ncol(Hk) == 1) break
        Hk[, use.min.index - 1] <- Hk[, use.min.index - 1] +
                                   Hk[, use.min.index]
        Hk <- Hk[, -use.min.index, drop = FALSE]
        for (kay in 2:length(Hlist)) {  # Omit the intercept
          Hlist[[kay]] <- Hk
        }  # kay
      } else {  # crossq
        break
      }
      if (is.numeric(maxit))
        object@control$maxit <- maxit  # Overwrite this possibly
      if (is.logical(trace))
        object@control$trace <- trace  # Overwrite this possibly
      object <-
        vglm(formula = as.formula(formula(object)),
             family = extlogF1(tau = object.save@extra$tau,
                             lambda.arg = object.save@extra$lambda.arg,
                             scale.arg = object.save@extra$scale.arg,
                             llocation = linkfun(object.save)[1],
                             parallel = FALSE),
             data = get(object.save@misc$dataname),
             control = object@control,  # maxit, trace, etc.
             constraints = Hlist)  # New constraints; new object too
    }  # jay
  }  # M > 1 and Hk not all parallel

  object  # A new object without any crossing quantiles
}  # fix.crossing



if (!isGeneric("fix.crossing"))
  setGeneric("fix.crossing", function(object, ...)
             standardGeneric("fix.crossing"),
             package = "VGAM")


setMethod("fix.crossing",  "vglm", function(object, ...)
    fix.crossing.vglm(object, ...))







if (FALSE)
qtplot.extlogF1 <- function(lambda,
                          tau = c(0.25, 0.5, 0.75),
                          location = 0, scale = 1,
                          eta = NULL) {

    
  lp <- length(tau)
  lambda <- rep(lambda, length = lp)
  answer <- matrix(NA_real_, nrow(eta), lp,
                   dimnames = list(dimnames(eta)[[1]],
                                   as.character(tau)))
  for (ii in 1:lp) {
    answer[, ii] <- qtlogF(tau      = tau[ii],
                           lambda   = lambda[ii],
                           location = location,  # eta[, ii],
                           scale    = scale)
  }
  answer
}








 eCDF.vglm <-
     function(object, all = FALSE,
              ...) {

  okayfuns <- c("lms.bcn", "extlogF1")
  if (!any(object@family@vfamily %in% okayfuns))
    stop("the object does not have an empirical CDF defined on it")

  
  if (any(object@family@vfamily == "extlogF1")) {
    Ans <- object@extra$eCDF / 100
    if (all) {
      Ans <- cbind(Ans, object@extra$tau)
      rownames(Ans) <- rep(" ", NROW(Ans))  # NULL  #
      colnames(Ans) <- c("ecdf", "tau")
    }
  }  # "extlogF1"


  if (any(object@family@vfamily == "lms.bcn")) {
    Ans <- numeric(length(object@extra$percentiles))
    locat <- as.matrix(fitted(object))
    M <- npred(object)
    y <- depvar(object)
    w <- weights(object, type = "prior")
    for (ii in 1:M) {
      y.use <- if (ncol(y) > 1) y[, ii] else y
      Ans[ii] <- 1 *  # Was 100, but now 1
        weighted.mean(y.use <= locat[, ii], w[, min(ii, ncol(w))])
    }
    if (all) {
      Ans <- cbind(Ans, object@extra$percentiles / 100)  # tau, really
      rownames(Ans) <- rep(" ", NROW(Ans))  # NULL  #
      colnames(Ans) <- c("ecdf", "tau")
    } else {
      digt <- 4  # extlogF1() default
      tau.names <-
        paste("(tau = ", round(object@extra$percentiles / 100,
                               digits = digt), ")",
            sep = "")
      Y.names <- if (ncol(y) > 1) dimnames(y)[[2]] else "y"
      if (is.null(Y.names) || any(Y.names == ""))
        Y.names <- paste("y", 1:ncol(y), sep = "")
      y.names <- if (ncol(y) > 1)
        paste(Y.names, tau.names, sep = "") else tau.names
      names(Ans) <- y.names
    }
  }  # "lms.bcn"

  Ans
}  # eCDF.vglm



if (!isGeneric("eCDF"))
  setGeneric("eCDF", function(object, ...)
             standardGeneric("eCDF"),
             package = "VGAM")


setMethod("eCDF",  "vglm", function(object, ...)
    eCDF.vglm(object, ...))




         





























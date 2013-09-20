# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.











 double.cennormal <-
  function(r1 = 0, r2 = 0,
           lmu = "identity",
           lsd = "loge",
           imu = NULL, isd = NULL, zero = 2) {
  if (!is.Numeric(r1, allowable.length = 1, integer.valued = TRUE) ||
      r1 < 0)
    stop("bad input for 'r1'")
  if (!is.Numeric(r2, allowable.length = 1, integer.valued = TRUE) ||
      r2 < 0)
    stop("bad input for 'r2'")

  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")


  new("vglmff",
  blurb = c("Univariate Normal distribution with double censoring\n\n",
            "Links:    ",
            namesof("mu", lmu, earg = emu, tag = TRUE), ", ",
            namesof("sd", lsd, earg = esd, tag = TRUE),
            "\n",
            "Variance: sd^2"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }) , list( .zero = zero))),


  initialize = eval(substitute(expression({
    predictors.names <-
      c(namesof("mu", .lmu, earg =.emu, tag = FALSE),
        namesof("sd", .lsd, earg =.esd, tag = FALSE))

    if (ncol(y <- cbind(y)) != 1)
      stop("the response must be a vector or a one-column matrix")

    if (length(w) != n ||
      !is.Numeric(w, integer.valued = TRUE, positive = TRUE))
      stop("the argument 'weights' must be a vector ",
           "of positive integers")

    sumw <- sum(w)
    extra$bign <- sumw + .r1 + .r2  # Tot num of censored & uncensored obsns

    if (!length(etastart)) {
      yyyy.est <- if (length( .imu )) .imu else median(y)
      sd.y.est <- if (length( .isd )) .isd else {
        junk <- lm.wfit(x = x, y = c(y), w = c(w))
        1.25 * sqrt( sum(w * junk$resid^2) / junk$df.residual )
      }
      yyyy.est <- rep(yyyy.est , len = n)
      sd.y.est <- rep(sd.y.est , len = n)
      etastart <- cbind(mu = theta2eta(yyyy.est, .lmu , earg =.emu ),
                        sd = theta2eta(sd.y.est, .lsd , earg =.esd ))
    }
  }) , list( .lmu = lmu, .lsd = lsd,
             .emu = emu, .esd = esd,
             .imu = imu, .isd = isd,
             .r1 = r1, .r2 = r2 ))),
  linkinv = function(eta, extra = NULL) eta[, 1], 
  last = eval(substitute(expression({
    misc$link <-    c(mu = .lmu , sd = .lsd )

    misc$earg <- list(mu = .emu , sd = .esd )

    misc$multipleResponses <- FALSE
    misc$expected <- TRUE
    misc$r1 <- .r1
    misc$r2 <- .r2
  }) , list( .lmu = lmu, .lsd = lsd,
             .emu = emu, .esd = esd,
             .r1 = r1, .r2 = r2 ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    sd <- eta2theta(eta[, 2], .lsd, earg = .esd )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
      sum(w * dnorm(y, m = mu, sd = sd, log = TRUE)) +
      (if ( .r1 == 0) 0 else {
         z1 <- min((y - mu) / sd); Fz1 = pnorm(z1); .r1 * log(Fz1)}) +
      (if ( .r2 == 0) 0 else {
         z2 <- max((y - mu) / sd); Fz2 = pnorm(z2); .r2 * log1p(-Fz2)})
  } , list( .lmu = lmu, .lsd = lsd,
            .emu = emu, .esd = esd,
            .r1 = r1, .r2 = r2 ))),
  vfamily = c("double.cennormal"),
  deriv = eval(substitute(expression({
    sd <- eta2theta(eta[, 2], .lsd, earg =.esd)

    q1 <- .r1 / extra$bign
    q2 <- .r2 / extra$bign
    pee <- 1 - q1 - q2  # 1 if r1==r2==0
    z1 <- if ( .r1 == 0) - 100 else min((y - mu) / sd)  # 100==Inf
    z2 <- if ( .r2 == 0) + 100 else max((y - mu) / sd)  # 100==Inf
    fz1 <- if ( .r1 == 0) 0 else dnorm(z1)
    fz2 <- if ( .r2 == 0) 0 else dnorm(z2)
    Fz1 <- if ( .r1 == 0) 0.02 else pnorm(z1)  # 0/0 undefined
    Fz2 <- if ( .r2 == 0) 0.99 else pnorm(z2)

    dl.dmu <- (y - mu) / sd^2 +
             ((- .r1 * fz1/Fz1 + .r2 * fz2/(1-Fz2)) / sd) / (n*w)
    dl.dsd <- -1/sd + (y-mu)^2 / sd^3 +
             ((- .r1 * z1*fz1/Fz1 + .r2 * z2*fz2/(1-Fz2)) / sd) / (n*w)

    dmu.deta <- dtheta.deta(mu, .lmu, earg =.emu)
    dsd.deta <- dtheta.deta(sd, .lsd, earg =.esd)

    c(w) * cbind(dl.dmu * dmu.deta, dl.dsd * dsd.deta)
  }) , list( .lmu = lmu, .lsd = lsd,
             .emu = emu, .esd = esd,
             .r1 = r1, .r2 = r2 ))),
  weight=expression({
    wz <- matrix(as.numeric(NA), n, dimm(M))

    Q1 <- ifelse(q1 == 0, 1, q1)  # Saves division by 0 below; not elegant
    Q2 <- ifelse(q2 == 0, 1, q2)  # Saves division by 0 below; not elegant

    ed2l.dmu2 <- 1 / (sd^2) + 
                 ((fz1*(z1+fz1/Q1) - fz2*(z2-fz2/Q2)) / sd^2) / (pee*w)
    ed2l.dmusd <- ((fz1-fz2 + z1*fz1*(z1+fz1/Q1) -
                  z2*fz2*(z2-fz2/Q2)) / sd^2) / (pee*w)
    ed2l.dsd2 <- 2 / (sd^2) + 
                 ((z1*fz1-z2*fz2 + z1^2 *fz1 *(z1+fz1/Q1) -
                 z2^2 *fz2*(z2-fz2/Q2)) / sd^2) / (pee*w)

    wz[,iam(1,1,M)] <- w * ed2l.dmu2 * dmu.deta^2
    wz[,iam(2,2,M)] <- w * ed2l.dsd2 * dsd.deta^2
    wz[,iam(1,2,M)] <- w * ed2l.dmusd * dsd.deta * dmu.deta
    wz
  }))
}






dbisa <- function(x, shape, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  L <- max(length(x), length(shape), length(scale))
  if (length(x)     != L) x     <- rep(x,     len = L)
  if (length(shape) != L) shape <- rep(shape, len = L)
  if (length(scale) != L) scale <- rep(scale, len = L)

  logdensity <- rep(log(0), len = L)

  xok <- (x > 0)
  xifun <- function(x) {
    temp <- sqrt(x)
    temp - 1/temp
  }
  logdensity[xok] <-
    dnorm(xifun(x[xok] / scale[xok]) / shape[xok], log = TRUE) +
    log1p(scale[xok]/x[xok]) - log(2) - log(shape[xok]) -
    0.5 * log(x[xok]) - 0.5 * log(scale[xok])
  logdensity[scale <= 0] <- NaN
  logdensity[shape <= 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}


pbisa <- function(q, shape, scale = 1) {
  if (!is.Numeric(q))
    stop("bad input for argument 'q'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")
  if (!is.Numeric(scale, positive = TRUE))
    stop("bad input for argument 'scale'")
  ans <- pnorm(((temp <- sqrt(q/scale)) - 1/temp) / shape)
  ans[scale < 0 | shape < 0] <- NA
  ans[q <= 0] <- 0
  ans
}


qbisa <- function(p, shape, scale = 1) {
  if (!is.Numeric(p, positive = TRUE) || any(p >= 1))
      stop("argument 'p' must have values inside the interval (0,1)")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")
  if (!is.Numeric(scale, positive = TRUE))
    stop("bad input for argument 'scale'")
  A <- qnorm(p)
  temp1 <- A * shape * sqrt(4 + A^2 * shape^2)
  ans1 <- (2 + A^2 * shape^2 + temp1) * scale / 2
  ans2 <- (2 + A^2 * shape^2 - temp1) * scale / 2
  ifelse(p < 0.5, pmin(ans1, ans2), pmax(ans1, ans2))
}


rbisa <- function(n, shape, scale = 1) {

  A <- rnorm(n)
  temp1 <- A * shape
  temp1 <- temp1 * sqrt(4 + temp1^2)
  ans1 <- (2 + A^2 * shape^2 + temp1) * scale / 2
  ans2 <- (2 + A^2 * shape^2 - temp1) * scale / 2



  ans <- ifelse(A < 0, pmin(ans1, ans2), pmax(ans1, ans2))
  ans[shape <= 0] <- NaN
  ans[scale <= 0] <- NaN
  ans
}








 bisa <- function(lshape = "loge", lscale = "loge",
                  ishape = NULL,   iscale = 1,
                  imethod = 1, zero = NULL) {
  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")


  if (length(ishape) && !is.Numeric(ishape, positive = TRUE))
      stop("bad input for argument 'ishape'")
  if (!is.Numeric(iscale, positive = TRUE))
      stop("bad input for argument 'iscale'")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
      stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Birnbaum-Saunders distribution\n\n",
            "Links:    ",
            namesof("shape", lshape, earg = eshape, tag = TRUE), "; ",
            namesof("scale", lscale, earg = escale, tag = TRUE)),
  constraints = eval(substitute(expression({
      constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }) , list( .zero = zero))),
  initialize = eval(substitute(expression({
    if (ncol(y <- cbind(y)) != 1)
      stop("the response must be a vector or a one-column matrix")

    predictors.names <-
      c(namesof("shape", .lshape, earg = .eshape, tag = FALSE),
        namesof("scale", .lscale, earg = .escale, tag = FALSE))

    if (!length(etastart)) {
      scale.init <- rep( .iscale, len = n)
      shape.init <- if (is.Numeric( .ishape)) rep( .ishape, len = n) else {
      if ( .imethod == 1) {
        ybar <- rep(weighted.mean(y, w), len = n)
        ybarr <- rep(1 / weighted.mean(1/y, w), len = n)  # Reqrs y > 0
        sqrt(ybar / scale.init + scale.init / ybarr - 2)
      } else if ( .imethod == 2) {
        sqrt(2*( pmax(y, scale.init+0.1) / scale.init - 1))
      } else {
        ybar <- rep(weighted.mean(y, w), len = n)
        sqrt(2*(pmax(ybar, scale.init + 0.1) / scale.init - 1))
      }
    }
      etastart <- cbind(theta2eta(shape.init, .lshape, earg = .eshape),
                        theta2eta(scale.init, .lscale, earg = .escale))
    }
  }) , list( .lshape = lshape, .lscale = lscale,
             .ishape = ishape, .iscale = iscale,
             .eshape = eshape, .escale = escale,
             .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    sh <- eta2theta(eta[, 1], .lshape, earg = .eshape)
    sc <- eta2theta(eta[, 2], .lscale, earg = .escale)
    sc * (1 + sh^2 / 2)
  }, list( .lshape = lshape, .lscale = lscale,
           .eshape = eshape, .escale = escale ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape = .lshape, scale = .lscale)

    misc$earg <- list(shape = .eshape, scale = .escale)

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }) , list( .lshape = lshape, .lscale = lscale,
             .eshape = eshape, .escale = escale ))),
  loglikelihood = eval(substitute(
    function(mu,y,w,residuals = FALSE,eta, extra = NULL) {
    sh <- eta2theta(eta[, 1], .lshape , earg = .eshape )
    sc <- eta2theta(eta[, 2], .lscale , earg = .escale )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(w * dbisa(x = y, shape = sh, scale = sc, log = TRUE))
    }
  }, list( .lshape = lshape, .lscale = lscale,
           .eshape = eshape, .escale = escale ))),
  vfamily = c("bisa"),
  deriv = eval(substitute(expression({
    sh <- eta2theta(eta[, 1], .lshape, earg = .eshape)
    sc <- eta2theta(eta[, 2], .lscale, earg = .escale)

    dl.dsh <- ((y/sc - 2 + sc/y) / sh^2 - 1) / sh 
    dl.dsc <- -0.5 / sc + 1/(y+sc) + sqrt(y) * ((y+sc)/y) *
             (sqrt(y/sc) - sqrt(sc/y)) / (2 * sh^2 * sc^1.5)

    dsh.deta <- dtheta.deta(sh, .lshape, earg = .eshape)
    dsc.deta <- dtheta.deta(sc, .lscale, earg = .escale)

    c(w) * cbind(dl.dsh * dsh.deta,
                 dl.dsc * dsc.deta)
  }) , list( .lshape = lshape, .lscale = lscale,
             .eshape = eshape, .escale = escale ))),
  weight = eval(substitute(expression({
    wz <- matrix(as.numeric(NA), n, M)  # Diagonal!!
    wz[,iam(1,1,M)] <- 2 * dsh.deta^2 / sh^2
    hfunction <- function(alpha)
      alpha * sqrt(pi/2) - pi * exp(2/alpha^2) *
                           pnorm(2/alpha, lower.tail = FALSE)
    wz[,iam(2,2,M)] <- dsc.deta^2 *
                       (sh * hfunction(sh)  / sqrt(2*pi) + 1) / (sh*sc)^2
    c(w) * wz
  }), list( .zero = zero ))))
}



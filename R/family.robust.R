# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.















edhuber <- function(x, k = 0.862, mu = 0, sigma = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  zedd <- (x - mu) / sigma
  fk <- dnorm(k)
   eps <- 1 - 1 / (pnorm(k) - pnorm(-k) + 2 * fk /k)
  ceps <-     1 / (pnorm(k) - pnorm(-k) + 2 * fk / k)

  if (log.arg) {
    val <-  log(ceps) + dnorm(zedd, log = TRUE)
    val[zedd < (-k)] <- (log(ceps) + log(fk) +  ( k * (zedd+k)))[zedd < (-k)]
    val[zedd > (+k)] <- (log(ceps) + log(fk) +  (-k * (zedd-k)))[zedd > (+k)]
  } else {
    val <-  (ceps) * dnorm(zedd)
    val[zedd < (-k)] <- ((ceps) * fk * exp( k * (zedd + k)))[zedd < (-k)]
    val[zedd > (+k)] <- ((ceps) * fk * exp(-k * (zedd - k)))[zedd > (+k)]
  }
  list(val = if (log.arg) val - log(sigma) else val / sigma,
       eps = eps)
}


dhuber <- function(x, k = 0.862, mu = 0, sigma = 1, log = FALSE)
  edhuber(x, k, mu, sigma, log = log)$val





rhuber <- function(n, k = 0.862, mu = 0, sigma = 1) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  myl <- rep(0.0, len = use.n)

  lowlim <- 1
  upplim <- 0
  chunksize <- 2 * use.n
  while (lowlim <= use.n) {
    x <- rexp(chunksize)
    s <- sample(c(-1, 1), size = chunksize, replace = TRUE)
    y <- s*x/k
    u <- runif(chunksize)

    yok <- (abs(y) >= k | u <= exp(k * abs(y) - (k * k + y * y) / 2))
    sumyok <- sum(yok)
    if (sumyok > 0) {
      upplim <- upplim + sumyok

      if (upplim > use.n)
        myl <- rep(myl, len = upplim)

      myl[lowlim:upplim] <- y[yok]
      lowlim <- lowlim + sumyok
    }
  }
  myl <- rep(myl, len = use.n)  # Prune to right length

  rep(mu + sigma * myl, len = use.n)
}











qhuber <- function (p, k = 0.862, mu = 0, sigma = 1) {
  if(min(sigma) <= 0)
    stop("argument 'sigma' must be positive")
  if(min(k)     <= 0)
    stop("argument 'k' must be positive")

  cnorm <- sqrt(2 * pi) * ((2 * pnorm(k) - 1) + 2 * dnorm(k) / k)
  x <- pmin(p, 1 - p)
  q <- ifelse(x <= sqrt(2 * pi) * dnorm(k) / ( k * cnorm),
              log(k * cnorm * x) / k - k / 2,
              qnorm(abs(1 - pnorm(k) + x * cnorm / sqrt(2 * pi) -
                    dnorm(k) / k)))
  ifelse(p < 0.5, mu + q * sigma,
                  mu - q * sigma)
}




phuber <- function(q, k = 0.862, mu = 0, sigma = 1) {
  if (any(sigma <= 0))
    stop("argument 'sigma' must be positive")

  A1  <- (2 * dnorm(k) / k - 2 * pnorm(-k))
  eps <- A1 / (1 + A1)
  zedd <- (q - mu) / sigma
  x <- -abs(zedd)
  p <- ifelse(x <= -k ,
              exp(k^2 / 2) / k * exp(k * x) / sqrt(2 * pi),
              dnorm(k) / k + pnorm(x) - pnorm(-k))
  p <- p * (1 - eps)
  ifelse(zedd <= 0, p, 1 - p)
}





 huber2 <- function(llocation = "identity", lscale = "loge",
                    k = 0.862, imethod = 1, zero = 2) {


  A1 <- (2 * dnorm(k) / k - 2 * pnorm(-k))
  eps <- A1 / (1 + A1)

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 4)
       stop("argument 'imethod' must be 1 or 2 or 3 or 4")

  if (!is.Numeric(k, allowable.length = 1, positive = TRUE))
    stop("bad input for argument 'k'")

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  new("vglmff",
  blurb = c("Huber least favorable distribution\n\n",
            "Links: ",
            namesof("location",  llocat,  earg = elocat), ", ",
            namesof("scale",     lscale,  earg = escale), "\n\n",
            "Mean: location"),
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
      scale.y.est <- sqrt( sum(c(w) * junk$resid^2) / junk$df.residual )
      location.init <- if ( .llocat == "loge") pmax(1/1024, y) else {
        if ( .imethod == 3) {
          rep(weighted.mean(y, w), len = n)
        } else if ( .imethod == 2) {
          rep(median(rep(y, w)), len = n)
        } else if ( .imethod == 1) {
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
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat, earg = .elocat)
  }, list( .llocat = llocat,
           .elocat = elocat, .escale = escale ))),
  last = eval(substitute(expression({
    misc$link <-    c("location" = .llocat , "scale" = .lscale )

    misc$earg <- list("location" = .elocat , "scale" = .escale )

    misc$expected <- TRUE
    misc$k.huber <- .k
    misc$imethod <- .imethod
    misc$multipleResponses <- FALSE
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .k      = k,         .imethod = imethod ))),
 loglikelihood = eval(substitute(
   function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
   location <- eta2theta(eta[, 1], .llocat, earg = .elocat)
   myscale  <- eta2theta(eta[, 2], .lscale, earg = .escale)
   kay      <- .k
   if (residuals) stop("loglikelihood residuals not ",
                       "implemented yet") else {
     sum(c(w) * dhuber(y, k = kay, mu = location,  sigma = myscale,
                    log = TRUE))
   }
 }, list( .llocat = llocat, .lscale = lscale,
          .elocat = elocat, .escale = escale,
          .k      = k ))),
  vfamily = c("huber2"),
  deriv = eval(substitute(expression({
    mylocat <- eta2theta(eta[, 1], .llocat,  earg = .elocat)
    myscale <- eta2theta(eta[, 2], .lscale,  earg = .escale)
    myk     <- .k

    zedd <- (y - mylocat) / myscale
    cond2 <- (abs(zedd) <=  myk)
    cond3 <-     (zedd  >   myk)

    dl.dlocat        <- -myk + 0 * zedd  # cond1
    dl.dlocat[cond2] <- zedd[cond2]
    dl.dlocat[cond3] <-  myk  # myk is a scalar
    dl.dlocat <- dl.dlocat / myscale


    dl.dscale        <- (-myk * zedd)
    dl.dscale[cond2] <-      (zedd^2)[cond2]
    dl.dscale[cond3] <- ( myk * zedd)[cond3]
    dl.dscale <- (-1 + dl.dscale) / myscale

    dlocat.deta <- dtheta.deta(mylocat, .llocat, earg = .elocat)
    dscale.deta <- dtheta.deta(myscale, .lscale, earg = .escale)
    ans <- c(w) * cbind(dl.dlocat * dlocat.deta,
                        dl.dscale * dscale.deta)
    ans
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .eps    = eps,       .k      = k ))),
  weight = eval(substitute(expression({
    wz   <- matrix(as.numeric(NA), n, 2) # diag matrix; y is one-col too




    temp4 <- erf(myk / sqrt(2))
    ned2l.dlocat2 <- temp4 * (1 - .eps) / myscale^2

    ned2l.dscale2 <- (dnorm(myk) * (1 - myk^2) + temp4) *
                     2 * (1 - .eps) / (myk * myscale^2)

    wz[, iam(1, 1, M)] <- ned2l.dlocat2 * dlocat.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dscale2 * dscale.deta^2
    ans
    c(w) * wz
  }), list( .eps = eps ))))
}




 huber1 <- function(llocation = "identity",
                    k = 0.862,
                    imethod = 1) {


  A1 <- (2 * dnorm(k) / k - 2 * pnorm(-k))
  eps <- A1 / (1 + A1)

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 4)
    stop("argument 'imethod' must be 1 or 2 or 3 or 4")

  if (!is.Numeric(k, allowable.length = 1, positive = TRUE))
    stop("bad input for argument 'k'")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")


  new("vglmff",
  blurb = c("Huber least favorable distribution\n\n",
            "Links: ",
            namesof("location",  llocat,  earg = elocat), "\n\n",
            "Mean: location"),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <-
       c(namesof("location", .llocat, earg = .elocat, tag = FALSE))


    if (!length(etastart)) {
      junk <- lm.wfit(x = x, y = c(y), w = c(w))
      location.init <- if ( .llocat == "loge") pmax(1/1024, y) else {
        if ( .imethod == 3) {
          rep(weighted.mean(y, w), len = n)
        } else if ( .imethod == 2) {
          rep(median(rep(y, w)), len = n)
        } else if ( .imethod == 1) {
          junk$fitted
        } else {
          y
        }
      }
      etastart <- cbind(
           theta2eta(location.init,  .llocat, earg = .elocat))
    }
  }), list( .llocat = llocat,
            .elocat = elocat,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta, .llocat, earg = .elocat)
  }, list( .llocat = llocat,
           .elocat = elocat ))),
  last = eval(substitute(expression({
    misc$link <-    c("location" = .llocat )

    misc$earg <- list("location" = .elocat )

    misc$expected <- TRUE
    misc$k.huber <- .k
    misc$imethod <- .imethod
    misc$multipleResponses <- FALSE
  }), list( .llocat = llocat,
            .elocat = elocat,
            .k      = k,         .imethod = imethod ))),
 loglikelihood = eval(substitute(
   function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
   location <- eta2theta(eta, .llocat, earg = .elocat)
   kay      <- .k
   if (residuals) stop("loglikelihood residuals not ",
                       "implemented yet") else {
     sum(c(w) * dhuber(y, k = kay, mu = location,  sigma = 1,
                    log = TRUE))
   }
 }, list( .llocat = llocat,
          .elocat = elocat,
          .k      = k ))),
  vfamily = c("huber1"),
  deriv = eval(substitute(expression({
    mylocat <- eta2theta(eta, .llocat,  earg = .elocat)
    myk     <- .k

    zedd <- (y - mylocat) # / myscale
    cond2 <- (abs(zedd) <=  myk)
    cond3 <-     (zedd  >   myk)

    dl.dlocat        <- -myk + 0 * zedd # cond1
    dl.dlocat[cond2] <- zedd[cond2]
    dl.dlocat[cond3] <-  myk  # myk is a scalar
    dl.dlocat <- dl.dlocat # / myscale


    if (FALSE) {
      dl.dscale        <- (-myk * zedd)
      dl.dscale[cond2] <-      (zedd^2)[cond2]
      dl.dscale[cond3] <- ( myk * zedd)[cond3]
      dl.dscale <- (-1 + dl.dscale) / myscale
    }

    dlocat.deta <- dtheta.deta(mylocat, .llocat, earg = .elocat)
    ans <- c(w) * cbind(dl.dlocat * dlocat.deta)
    ans
  }), list( .llocat = llocat,
            .elocat = elocat,
            .eps    = eps,       .k      = k ))),
  weight = eval(substitute(expression({
    wz   <- matrix(as.numeric(NA), n, 1) # diag matrix; y is one-col too




    temp4 <- erf(myk / sqrt(2))
    ned2l.dlocat2 <- temp4 * (1 - .eps) # / myscale^2


    wz[, iam(1,1,M)] <- ned2l.dlocat2 * dlocat.deta^2
    ans
    c(w) * wz
  }), list( .eps = eps ))))
}




# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.















 Init.mu <-
  function(y, x = cbind("(Intercept)" = rep_len(1, nrow(as.matrix(y)))),
           w = x, imethod = 1, imu = NULL,
           ishrinkage = 0.95,
           pos.only = FALSE,
           probs.y = 0.35) {
    if (!is.matrix(x)) x <- as.matrix(x)
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(w)) w <- as.matrix(w)
    if (ncol(w) != ncol(y))
      w <- matrix(w, nrow = nrow(y), ncol = ncol(y))

    if (length(imu)) {
      MU.INIT <- matrix(imu, nrow(y), ncol(y), byrow = TRUE)
      return(MU.INIT)
    }


    if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 || ishrinkage > 1)
     warning("bad input for argument 'ishrinkage'; ",
             "using the value 0.95 instead")


    if (imethod > 6) {
      warning("argument 'imethod' should be 1 or 2 or... 6; ",
              "using the value 1")
      imethod <- 1
    }
    mu.init <- y
    for (jay in 1:ncol(y)) {
      TFvec <- if (pos.only) y[, jay] > 0 else TRUE
      locn.est <- if ( imethod %in% c(1, 4)) {
        weighted.mean(y[TFvec, jay], w[TFvec, jay]) + 1/16
      } else if ( imethod %in% c(3, 6)) {
        c(quantile(y[TFvec, jay], probs = probs.y ) + 1/16)
      } else {
        median(y[TFvec, jay]) + 1/16
      }

      if (imethod <= 3) {
        mu.init[, jay] <-      ishrinkage   * locn.est +
                          (1 - ishrinkage ) * y[, jay]
      } else {
        medabsres <- median(abs(y[, jay] - locn.est)) + 1/32
        allowfun <- function(z, maxtol = 1)
          sign(z) * pmin(abs(z), maxtol)
        mu.init[, jay] <- locn.est + (1 - ishrinkage ) *
                          allowfun(y[, jay] - locn.est, maxtol = medabsres)

        mu.init[, jay] <- abs(mu.init[, jay]) + 1 / 1024
      }
    }  # of for (jay)

    mu.init
  }








EIM.NB.specialp <-
  function(mu, size,
           y.max = NULL,  # Must be an integer
           cutoff.prob = 0.995,
           intercept.only = FALSE,
           extra.bit = TRUE) {


  if (intercept.only) {
    mu <- mu[1]
    size <- size[1]
  }

  y.min <- 0  # A fixed constant really

  if (!is.numeric(y.max)) {
    eff.p <- sort(c(cutoff.prob, 1 - cutoff.prob))
    y.max <- max(round(qnbinom(p = eff.p[2], mu = mu, size = size) * 1.1)) + 30
  }

  Y.mat <- if (intercept.only) y.min:y.max else
           matrix(y.min:y.max, length(mu), y.max-y.min+1, byrow = TRUE)
  neff.row <- ifelse(intercept.only, 1, nrow(Y.mat))
  neff.col <- ifelse(intercept.only, length(Y.mat), ncol(Y.mat))

  if (FALSE) {
  trigg.term <- if (intercept.only) {
    check2 <-
     sum(pnbinom(Y.mat, size = size, mu = mu, lower.tail = FALSE)
         / (Y.mat + size)^2)
    check2
  } else {
  check2 <-
    rowSums(pnbinom(Y.mat, size = size, mu = mu, lower.tail = FALSE)
            / (Y.mat + size)^2)
  check2
  }
  }  # FALSE


  if (TRUE) {
    answerC <- .C("eimpnbinomspecialp",
      as.integer(intercept.only),
      as.double(neff.row), as.double(neff.col),
      as.double(size),
      as.double(pnbinom(Y.mat, size = size, mu = mu, lower.tail = FALSE)),
      rowsums = double(neff.row))
    trigg.term <- answerC$rowsums
  }  # TRUE

  ned2l.dk2 <- trigg.term
  if (extra.bit)
    ned2l.dk2 <- ned2l.dk2 - 1 / size + 1 / (size + mu)
  ned2l.dk2
}  # EIM.NB.specialp()







EIM.NB.speciald <-
  function(mu, size,
           y.min = 0,  # 20160201; must be an integer
           y.max = NULL,  # Must be an integer
           cutoff.prob = 0.995,
           intercept.only = FALSE,
           extra.bit = TRUE) {





  if (intercept.only) {
    mu <- mu[1]
    size <- size[1]
  }

  if (!is.numeric(y.max)) {
    eff.p <- sort(c(cutoff.prob, 1 - cutoff.prob))
    y.max <- max(round(qnbinom(p = eff.p[2], mu = mu, size = size) * 1.1)) + 30
  }

  Y.mat <- if (intercept.only) y.min:y.max else
           matrix(y.min:y.max, length(mu), y.max-y.min+1, byrow = TRUE)
  trigg.term <- if (intercept.only) {
     dnbinom(Y.mat, size = size, mu = mu) %*% trigamma(Y.mat + size)
  } else {
     rowSums(dnbinom(Y.mat, size = size, mu = mu) *
             trigamma(Y.mat + size))
  }
  ned2l.dk2 <- trigamma(size) - trigg.term
  if (extra.bit)
    ned2l.dk2 <- ned2l.dk2 - 1 / size + 1 / (size + mu)
  ned2l.dk2
}  # end of EIM.NB.speciald()



NBD.Loglikfun2 <- function(munbval, sizeval,
                           y, x, w, extraargs) {
  sum(c(w) * dnbinom(x = y, mu = munbval,
                     size = sizeval, log = TRUE))
}





negbinomial.initialize.yj <-
  function(yvec, wvec = rep(1, length(yvec)),
           gprobs.y = ppoints(9),
           wm.yj = weighted.mean(yvec, w = wvec)) {
  try.mu <- c(quantile(yvec, probs = gprobs.y) + 1/16,
              wm.yj)
  if (median(try.mu) < 1) {
    y.pos <- yvec[yvec > 0]
    try.mu <- c(min(try.mu),  # 0.25,
                wm.yj,
                summary.default(y.pos)[c(1:3, 5)],
                quantile(y.pos, probs = gprobs.y) - 1/16)

  }
  unique(sort(try.mu))
}



negbinomial.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}



 negbinomial <-
  function(
           zero = "size",
           parallel = FALSE,
           deviance.arg = FALSE,
           type.fitted = c("mean", "quantiles"),
           percentiles = c(25, 50, 75),
           mds.min = 1e-3,
           nsimEIM = 500, cutoff.prob = 0.999,  # Maxiter = 5000,
           eps.trig = 1e-7,
           max.support = 4000,
           max.chunk.MB = 30,  # max.memory = Inf is allowed
           lmu = "loge", lsize = "loge",
           imethod = 1,
           imu = NULL,
           iprobs.y = NULL,  # 0.35,
           gprobs.y = ppoints(6),
           isize = NULL,
           gsize.mux = exp(c(-30, -20, -15, -10, -6:3))) {







  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
      stop("argument 'imethod' must be 1 or 2")


  if (!is.logical( deviance.arg ) || length( deviance.arg ) != 1)
    stop("argument 'deviance.arg' must be TRUE or FALSE")



  type.fitted <- match.arg(type.fitted,
                           c("mean", "quantiles"))[1]


  lmunb <- as.list(substitute(lmu))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  imunb <- imu

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")


  if (!is.Numeric(eps.trig, length.arg = 1,
                  positive = TRUE) || eps.trig > 1e-5)
    stop("argument 'eps.trig' must be positive and smaller in value")

  if (length(imunb) && !is.Numeric(imunb, positive = TRUE))
    stop("bad input for argument 'imu'")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("bad input for argument 'isize'")

  if (!is.Numeric(cutoff.prob, length.arg = 1) ||
    cutoff.prob < 0.95 ||
    cutoff.prob >= 1)
    stop("range error in the argument 'cutoff.prob'; ",
         "a value in [0.95, 1) is needed")

    if (!is.Numeric(nsimEIM, length.arg = 1, integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
    if (nsimEIM <= 10)
      warning("argument 'nsimEIM' should be an integer ",
               "greater than 10, say")


    if (is.logical( parallel ) && parallel  && length(zero))
      stop("need to set 'zero = NULL' when parallel = TRUE")



  ans <-
  new("vglmff",


  blurb = c("Negative binomial distribution\n\n",
            "Links:    ",
            namesof("mu",   lmunb, earg = emunb), ", ",
            namesof("size", lsize, earg = esize), "\n",
            "Mean:     mu\n",
            "Variance: mu * (1 + mu / size) for NB-2"),

  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints = constraints)

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .parallel = parallel, .zero = zero ))),



  infos = eval(substitute(function(...) {
    list(M1    = 2,
         Q1    = 1,
         expected = TRUE,
         imethod = .imethod ,
         mds.min = .mds.min ,
         multipleResponses = TRUE,
         parameters.names = c("mu", "size"),
         type.fitted = .type.fitted ,
         percentiles = .percentiles ,
         lmu   = .lmunb ,
         lsize = .lsize ,
         nsimEIM = .nsimEIM ,
         eps.trig = .eps.trig ,
         zero  = .zero ,
         max.chunk.MB = .max.chunk.MB ,
         cutoff.prob = .cutoff.prob
         )
  }, list( .zero = zero, .lsize = lsize, .lmunb = lmunb,
           .type.fitted = type.fitted,
           .percentiles = percentiles ,
           .eps.trig = eps.trig,
           .imethod = imethod,
           .mds.min = mds.min,
           .cutoff.prob = cutoff.prob,
           .max.chunk.MB = max.chunk.MB,
           .nsimEIM = nsimEIM ))),

  initialize = eval(substitute(expression({
    M1 <- 2

    temp12 <-
      w.y.check(w = w, y = y,
                Is.nonnegative.y = TRUE,
                Is.integer.y = TRUE,
                ncol.w.max = Inf,
                ncol.y.max = Inf,
                out.wy = TRUE,
                colsyperw = 1, maximize = TRUE)
    w <- temp12$w
    y <- temp12$y


    assign("CQO.FastAlgorithm",
          ( .lmunb == "loge") && ( .lsize == "loge"),
           envir = VGAMenv)

    if (any(function.name == c("cqo", "cao")) &&
        ((is.Numeric( .zero , length.arg = 1) && .zero != -2) ||
         (is.character( .zero ) && .zero != "size")))
        stop("argument zero = 'size' or zero = -2 is required")


    extra$type.fitted <- .type.fitted
    extra$percentiles <- .percentiles
    extra$colnames.y  <- colnames(y)
    M <- M1 * ncol(y)
    NOS <- ncoly <- ncol(y)  # Number of species
    predictors.names <-
     c(namesof(param.names("mu",   NOS),
                .lmunb , earg = .emunb , tag = FALSE),
       namesof(param.names("size", NOS),
                .lsize , earg = .esize , tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]

    gprobs.y <- .gprobs.y
    imunb <- .imunb  # Default is NULL
    if (length(imunb))
      imunb <- matrix(imunb, n, NOS, byrow = TRUE)

    if (!length(etastart)) {
      munb.init <-
      size.init <- matrix(NA_real_, n, NOS)
      if (length( .iprobs.y ))
        gprobs.y <-  .iprobs.y
      gsize.mux <- .gsize.mux  # gsize.mux is on a relative scale

      for (jay in 1:NOS) {  # For each response 'y_jay'... do:
        wm.yj <- weighted.mean(y[, jay], w = w[, jay])
        munb.init.jay <- if ( .imethod == 1 ) {
          negbinomial.initialize.yj(y[, jay], w[, jay],
                                    gprobs.y = gprobs.y,
                                    wm.yj = wm.yj)
        } else {
          wm.yj
        }
        if (length(imunb))
          munb.init.jay <- imunb[, jay]


        gsize <- gsize.mux * 0.5 * (mean(munb.init.jay) +
                                    wm.yj)
        if (length( .isize ))
          gsize <- .isize  # isize is on an absolute scale


        try.this <-
          grid.search2(munb.init.jay, gsize,
                       objfun = NBD.Loglikfun2,
                       y = y[, jay], w = w[, jay],
                       ret.objfun = TRUE)  # Last value is the loglik

        munb.init[, jay] <- try.this["Value1"]
        size.init[, jay] <- try.this["Value2"]
      }  # for (jay ...)



      newemu <- .emunb
      if ( .lmunb == "nbcanlink") {
        newemu$size <- size.init
        testing1 <- log(munb.init / (munb.init + size.init))
        testing2 <- theta2eta(munb.init, link = .lmunb , earg = newemu )
      }


      etastart <-
        cbind(theta2eta(munb.init, link = .lmunb , earg = newemu ),
              theta2eta(size.init, link = .lsize , earg = .esize ))
      etastart <-
        etastart[, interleave.VGAM(M, M1 = M1), drop = FALSE]
    }
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize,
            .imunb = imunb,
            .gprobs.y = gprobs.y, .gsize.mux = gsize.mux,
            .deviance.arg = deviance.arg,
            .isize = isize, .iprobs.y = iprobs.y,
            .nsimEIM = nsimEIM,
            .zero = zero, .imethod = imethod,
            .type.fitted = type.fitted,
            .percentiles = percentiles ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- ncol(eta) / c(M1 = 2)
    kmat <- NULL

    munb <- if ( .lmunb == "nbcanlink") {

      newemu <- .emunb
      kmat <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                        .lsize , earg = .esize )
      newemu$size <- kmat
      check.munb <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                .lmunb , earg = newemu )


      munb <- kmat / expm1(-eta[, c(TRUE, FALSE), drop = FALSE])
    } else {
      eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                .lmunb , earg = .emunb )
    }

   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }
    type.fitted <- match.arg(type.fitted,
                     c("mean", "quantiles"))[1]
    if (type.fitted == "mean") {
      return(label.cols.y(munb, colnames.y = extra$colnames.y,
                          NOS = NOS))
    }


    if (is.null(kmat))
      kmat <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                        .lsize , earg = .esize )
    percvec <- extra$percentiles
    lenperc <- length(percvec)
    jvec <- lenperc * (0:(NOS - 1))
    ans <- matrix(0, nrow(eta), lenperc * NOS)
    for (kay in 1:lenperc)
      ans[, jvec + kay] <-
        qnbinom(0.01 * percvec[kay], mu = munb, size = kmat)

    rownames(ans) <- rownames(eta)


    label.cols.y(ans, colnames.y = extra$colnames.y,
                 NOS = NOS, percentiles = percvec,
                 one.on.one = FALSE)
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize))),

  last = eval(substitute(expression({
    if (exists("CQO.FastAlgorithm", envir = VGAMenv))
        rm("CQO.FastAlgorithm", envir = VGAMenv)


    if (function.name == "cao")
      ind2 <- FALSE


    save.weights <- control$save.weights <- !all(ind2)


    temp0303 <- c(rep_len( .lmunb , NOS),
                  rep_len( .lsize , NOS))
    names(temp0303) <- c(param.names("mu",   NOS),
                         param.names("size", NOS))
    misc$link <- temp0303[interleave.VGAM(M, M1 = M1)] # Already named

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- newemu
      misc$earg[[M1*ii  ]] <- .esize
    }
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize ))),

  linkfun = eval(substitute(function(mu, extra = NULL) {
    M1 <- 2

    newemu <- .emunb

    eta.temp <- theta2eta(mu, .lmunb , earg = newemu)
    eta.size <- theta2eta(if (is.numeric( .isize )) .isize else 1.0,
                          .lsize , earg = .esize )
    eta.size <- 0 * eta.temp + eta.size  # Right dimension now.



    if ( .lmunb == "nbcanlink") {
      newemu$size <- eta2theta(eta.size, .lsize , earg = .esize )
    }



    eta.temp <- cbind(eta.temp, eta.size)
    eta.temp[, interleave.VGAM(ncol(eta.temp), M1 = M1), drop = FALSE]
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize,
                           .isize = isize ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE,  eta,
             extra = NULL, summation = TRUE) {
    munb <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                      .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                      .lsize , earg = .esize )



    newemu <- .emunb
    if ( .lmunb == "nbcanlink") {
      newemu$size <- kmat
    }

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dnbinom(x = y, mu = munb, size = kmat, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize))),

  vfamily = c("negbinomial"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    vecTF <- c(TRUE, FALSE)
    munb <- cbind(eta2theta(eta[,  vecTF], .lmunb , earg = .emunb ))
    size <- cbind(eta2theta(eta[, !vecTF], .lsize , earg = .esize ))
    rnbinom(nsim * length(munb), mu = munb, size = size)
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize ))),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    munb <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                      .lmunb , earg = .emunb )
    size <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                      .lsize , earg = .esize )

    smallval <- .mds.min  # .munb.div.size
    okay1 <- all(is.finite(munb)) && all(0 < munb) &&
             all(is.finite(size)) && all(0 < size)


    okay0 <- if ( .lmunb == "nbcanlink") all(eta < 0) else TRUE


    overdispersion <- if (okay1) all(smallval < munb / size) else FALSE
    if (!overdispersion)
      warning("parameter 'size' has very large values relative ",
              "to 'mu'; ",
              "try fitting a quasi-Poisson ",
              "model instead.")
    okay1 && overdispersion && okay0
  }, list( .lmunb = lmunb, .emunb = emunb,
           .lsize = lsize, .esize = esize,
           .mds.min = mds.min))),



  deriv = eval(substitute(expression({




  odd.iter <- 1   # iter %% 2
  even.iter <- 1  # 1 - odd.iter

  if ( iter == 1 && .deviance.arg ) {
    if (control$criterion != "coefficients" &&
        control$half.step)
      warning("Argument 'criterion' should be 'coefficients' ",
               "or 'half.step' should be 'FALSE' when ",
              "'deviance.arg = TRUE'")



    low.index <- ifelse(names(constraints)[1] == "(Intercept)", 2, 1)
    if (low.index <= length(constraints))
    for (iii in low.index:length(constraints)) {
      conmat <- constraints[[iii]]
      if (any(conmat[c(FALSE, TRUE), ] != 0))
        stop("argument 'deviance.arg' should only be TRUE for NB-2 ",
             "models; ",
             "non-zero elements detected for the 'size' parameter." )
    }
  }






    M1 <- 2
    NOS <- ncol(eta) / M1
    munb <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                      .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                      .lsize , earg = .esize )


    smallval <- .mds.min  # Something like this is needed
    if (any(big.size <- (munb / kmat < smallval))) {
      kmat[big.size] <- munb[big.size] / smallval
    }


    newemu <- .emunb
    if ( .lmunb == "nbcanlink") {
      newemu$size <- kmat
    }


    dl.dmunb <- y / munb - (1 + y/kmat) / (1 + munb/kmat)
    dl.dsize <- digamma(y + kmat) - digamma(kmat) -
                (y - munb) / (munb + kmat) + log1p(-munb / (kmat + munb))
    if (any(big.size)) {
      dl.dsize[big.size] <- 1e-8  # A small number
    }


    dsize.deta <- dtheta.deta(kmat, .lsize , earg = .esize )


    myderiv <- if ( .lmunb == "nbcanlink") {
      dmunb.deta1 <- 1 / nbcanlink(munb, size = kmat, wrt.param = 1,
                                   deriv = 1)

      dsize.deta1 <- 1 / nbcanlink(munb, size = kmat, wrt.param = 2,
                                   deriv = 1)
      c(w) * cbind(dl.dmunb * dmunb.deta1 *  odd.iter +
                   dl.dsize * dsize.deta1 * 1 * even.iter,
                   dl.dsize * dsize.deta  * even.iter)
    } else {
      dmunb.deta <- dtheta.deta(munb, .lmunb , earg = .emunb )
      c(w) * cbind(dl.dmunb * dmunb.deta,
                   dl.dsize * dsize.deta)
    }


    myderiv <- myderiv[, interleave.VGAM(M, M1 = M1)]
    myderiv
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize,
            .deviance.arg = deviance.arg,
            .mds.min = mds.min ))),



  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, n, M)


    max.support <- .max.support
    max.chunk.MB <- .max.chunk.MB


    ind2 <- matrix(FALSE, n, NOS)  # Used for SFS
    for (jay in 1:NOS) {
      eff.p <- sort(c( .cutoff.prob , 1 - .cutoff.prob ))
      Q.mins <- 0
      Q.maxs <- round(qnbinom(p = eff.p[2],
                              mu = munb[, jay],
                              size = kmat[, jay]) * 1.1) + 30


      eps.trig <- .eps.trig
      Q.MAXS <- if ( .lsize == "loge")
        pmax(10, ceiling(kmat[, jay] / sqrt(eps.trig))) else Inf
      Q.maxs <- pmin(Q.maxs, Q.MAXS)



      ind1 <- if (max.chunk.MB > 0)
        (Q.maxs - Q.mins < max.support) else FALSE
      if ((NN <- sum(ind1)) > 0) {
        Object.Size <- NN * 8 * max(Q.maxs - Q.mins) / (2^20)
        n.chunks <- if (intercept.only) 1 else
                    max(1, ceiling( Object.Size / max.chunk.MB))
        chunk.rows <- ceiling(NN / n.chunks)
        ind2[, jay] <- ind1  # Save this
        wind2 <- which(ind1)


        upr.ptr <- 0
        lwr.ptr <- upr.ptr + 1
        while (lwr.ptr <= NN) {
          upr.ptr <- min(upr.ptr + chunk.rows, NN)
          sind2 <- wind2[lwr.ptr:upr.ptr]
          wz[sind2, M1*jay] <-
            EIM.NB.specialp(mu          =   munb[sind2, jay],
                            size        = kmat[sind2, jay],
                            y.max = max(Q.maxs[sind2]),
                            cutoff.prob = .cutoff.prob ,
                            intercept.only = intercept.only)


          if (any(eim.kk.TF <- wz[sind2, M1*jay] <= 0)) {
            ind2[sind2[eim.kk.TF], jay] <- FALSE
          }


          lwr.ptr <- upr.ptr + 1
        }  # while
      }  # if
    }  # end of for (jay in 1:NOS)









    for (jay in 1:NOS) {
      run.varcov <- 0
      ii.TF <- !ind2[, jay]  # Not assigned above
      if (any(ii.TF)) {
        kkvec <- kmat[ii.TF, jay]
        muvec <-   munb[ii.TF, jay]
        for (ii in 1:( .nsimEIM )) {
          ysim <- rnbinom(sum(ii.TF), mu = muvec, size = kkvec)
          dl.dsize <- digamma(ysim + kkvec) - digamma(kkvec) -
                      (ysim - muvec) / (muvec + kkvec) +
                      log1p( -muvec / (kkvec + muvec))
          run.varcov <- run.varcov + dl.dsize^2
        }  # end of for loop

        run.varcov <- c(run.varcov / .nsimEIM )
        ned2l.dsize2 <- if (intercept.only)
          mean(run.varcov) else run.varcov

        wz[ii.TF, M1*jay] <- ned2l.dsize2
      }
    }



    save.weights <- !all(ind2)



    ned2l.dmunb2 <- 1 / munb - 1 / (munb + kmat)
    ned2l.dsize2 <- wz[, M1*(1:NOS), drop = FALSE]


    if ( .lmunb == "nbcanlink") {
      wz <- cbind(wz, matrix(0, n, M-1))  # Make it tridiagonal

      wz[,     M1*(1:NOS) - 1] <-
        (ned2l.dmunb2 * (munb/kmat)^2 * odd.iter +
         ned2l.dsize2 * even.iter * 1) *
          (munb + kmat)^2



      wz[, M + M1*(1:NOS) - 1] <-
        -(munb + kmat) * ned2l.dsize2 * dsize.deta * even.iter
    } else {
      wz[, c(TRUE, FALSE)] <- ned2l.dmunb2 * dmunb.deta^2
    }


    wz[, M1*(1:NOS)] <- wz[, M1*(1:NOS)] * dsize.deta^2


    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }), list( .cutoff.prob = cutoff.prob,
            .max.support = max.support,
            .max.chunk.MB = max.chunk.MB,
            .lmunb = lmunb, .lsize = lsize,
            .eps.trig = eps.trig,
            .nsimEIM = nsimEIM ))))




  if (deviance.arg) {
    ans@deviance <- eval(substitute(
      function(mu, y, w, residuals = FALSE, eta, extra = NULL,
               summation = TRUE) {






    size <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                      .lsize , earg = .esize )

    if (residuals) {
      stop("this part of the function has not been written yet.")
    } else {
      dev.elts <- 2 * c(w) *
                  (y * log(pmax(1, y) / mu) -
                  (y + size) * log((y + size) / (mu + size)))
      if (summation) {
        sum(dev.elts)
      } else {
        dev.elts
      }
    }
  }, list( .lsize = lsize, .esize = esize,
           .lmunb = lmunb, .emunb = emunb )))





  }





  ans
}  # negbinomial()








polya.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}



 polya <-
  function(
           zero = "size",
           type.fitted = c("mean", "prob"),
           mds.min = 1e-3,
           nsimEIM = 500,  cutoff.prob = 0.999,  # Maxiter = 5000,
           eps.trig = 1e-7,
           max.support = 4000,
           max.chunk.MB = 30,  # max.memory = Inf is allowed
           lprob = "logit", lsize = "loge",
           imethod = 1,
           iprob = NULL,
           iprobs.y = NULL,
           gprobs.y = ppoints(6),
           isize = NULL,
           gsize.mux = exp(c(-30, -20, -15, -10, -6:3)),
           imunb = NULL) {


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
      stop("argument 'imethod' must be 1 or 2")


  deviance.arg <- FALSE  # 20131212; for now

  type.fitted <- match.arg(type.fitted,
                           c("mean", "prob"))[1]



  if (length(iprob) && !is.Numeric(iprob, positive = TRUE))
    stop("bad input for argument 'iprob'")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("bad input for argument 'isize'")

  if (!is.Numeric(eps.trig, length.arg = 1,
                  positive = TRUE) || eps.trig > 0.001)
    stop("argument 'eps.trig' must be positive and smaller in value")

  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE))
    stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 10)
    warning("argument 'nsimEIM' should be an integer ",
            "greater than 10, say")


  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")



  ans <-
  new("vglmff",
  blurb = c("Polya (negative-binomial) distribution\n\n",
            "Links:    ",
            namesof("prob", lprob, earg = eprob), ", ",
            namesof("size", lsize, earg = esize), "\n",
            "Mean:     size * (1 - prob) / prob\n",
            "Variance: mean / prob"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         mds.min = .mds.min ,
         type.fitted  = .type.fitted ,
         eps.trig = .eps.trig ,
         parameters.names = c("prob", "size"),
         zero = .zero)
  }, list( .zero = zero, .eps.trig = eps.trig,
           .type.fitted = type.fitted,
           .mds.min = mds.min))),

  initialize = eval(substitute(expression({
    M1 <- 2
    if (any(function.name == c("cqo", "cao")))
      stop("polya() does not work with cqo() or cao(). ",
           "Try negbinomial()")


    temp12 <- w.y.check(w = w, y = y,
              Is.integer.y = TRUE,
              Is.nonnegative = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1, maximize = TRUE)
    w <- temp12$w
    y <- temp12$y


    M <- M1 * ncol(y)
    NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    predictors.names <-
      c(namesof(param.names("prob", NOS), .lprob , earg = .eprob ,
                tag = FALSE),
        namesof(param.names("size", NOS), .lsize , earg = .esize ,
                tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]

    if (is.null( .nsimEIM )) {
       save.weights <- control$save.weights <- FALSE
    }


    gprobs.y <- .gprobs.y
    imunb <- .imunb  # Default in NULL
    if (length(imunb))
      imunb <- matrix(imunb, n, NOS, byrow = TRUE)



    if (!length(etastart)) {
      munb.init <-
      size.init <- matrix(NA_real_, n, NOS)
      gprobs.y  <- .gprobs.y
      if (length( .iprobs.y ))
        gprobs.y <-  .iprobs.y
      gsize.mux <- .gsize.mux  # gsize.mux is on a relative scale

      for (jay in 1:NOS) {  # For each response 'y_jay'... do:
        munb.init.jay <- if ( .imethod == 1 ) {
          quantile(y[, jay], probs = gprobs.y) + 1/16
        } else {
          weighted.mean(y[, jay], w = w[, jay])
        }
        if (length(imunb))
          munb.init.jay <- imunb[, jay]


        gsize <- gsize.mux * 0.5 * (mean(munb.init.jay) +
                                    weighted.mean(y[, jay], w = w[, jay]))
        if (length( .isize ))
          gsize <- .isize  # isize is on an absolute scale


        try.this <-
          grid.search2(munb.init.jay, gsize,
                       objfun = NBD.Loglikfun2,
                       y = y[, jay], w = w[, jay],
                       ret.objfun = TRUE)  # Last value is the loglik

        munb.init[, jay] <- try.this["Value1"]
        size.init[, jay] <- try.this["Value2"]
      }  # for (jay ...)






      prob.init <- if (length( .iprob ))
                   matrix( .iprob , nrow(y), ncol(y), byrow = TRUE) else
                   size.init / (size.init + munb.init)


      etastart <-
        cbind(theta2eta(prob.init, .lprob , earg = .eprob),
              theta2eta(size.init, .lsize , earg = .esize))
      etastart <-
        etastart[, interleave.VGAM(M, M1 = M1), drop = FALSE]
      }
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize,
            .iprob = iprob, .isize = isize,
            .pinit = iprob,
            .gprobs.y = gprobs.y, .gsize.mux = gsize.mux,
            .iprobs.y = iprobs.y,
            .nsimEIM = nsimEIM, .zero = zero,
            .imethod = imethod , .imunb = imunb,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- ncol(eta) / c(M1 = 2)
    pmat <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                      .lprob , earg = .eprob )
    kmat <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                      .lsize , earg = .esize )

   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob"))[1]

    ans <- switch(type.fitted,
                  "mean"      = kmat * (1 - pmat) / pmat,
                  "prob"      = pmat)


    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lprob = lprob, .eprob = eprob,
           .lsize = lsize, .esize = esize))),
  last = eval(substitute(expression({
    temp0303 <- c(rep_len( .lprob , NOS),
                  rep_len( .lsize , NOS))
    names(temp0303) <- c(param.names("prob", NOS),
                         param.names("size", NOS))
    temp0303 <- temp0303[interleave.VGAM(M, M1 = M1)]
    misc$link <- temp0303  # Already named

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .eprob
      misc$earg[[M1*ii  ]] <- .esize
    }

    misc$isize <- .isize
    misc$imethod <- .imethod
    misc$nsimEIM <- .nsimEIM
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize,
            .isize = isize,
            .nsimEIM = nsimEIM,
            .imethod = imethod ))),


  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pmat  <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                       .lprob , earg = .eprob)
    temp300 <-         eta[, c(FALSE, TRUE), drop = FALSE]
    if ( .lsize == "loge") {
      bigval <- 68
      temp300[temp300 >  bigval] <-  bigval
      temp300[temp300 < -bigval] <- -bigval
    }
    kmat <- eta2theta(temp300, .lsize , earg = .esize )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dnbinom(y, prob = pmat, size = kmat, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lsize = lsize, .lprob = lprob,
           .esize = esize, .eprob = eprob ))),
  vfamily = c("polya"),



  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pmat <- eta2theta(eta[, c(TRUE, FALSE)], .lprob , .eprob )
    kmat <- eta2theta(eta[, c(FALSE, TRUE)], .lsize , .esize )
    rnbinom(nsim * length(pmat), prob = pmat, size = kmat)
  }, list( .lprob = lprob, .lsize = lsize,
           .eprob = eprob, .esize = esize ))),



  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pmat <- eta2theta(eta[, c(TRUE, FALSE)], .lprob , .eprob )
    size <- eta2theta(eta[, c(FALSE, TRUE)], .lsize , .esize )
    munb <- size * (1 / pmat - 1)

    smallval <- .mds.min  # .munb.div.size
    okay1 <- all(is.finite(munb)) && all(0 < munb) &&
             all(is.finite(size)) && all(0 < size) &&
             all(is.finite(pmat)) && all(0 < pmat & pmat < 1)
    overdispersion <- if (okay1) all(munb / size > smallval) else FALSE
    if (!overdispersion)
      warning("parameter 'size' has very large values; ",
              "try fitting a quasi-Poisson ",
              "model instead.")
    okay1 && overdispersion
  }, list( .lprob = lprob, .eprob = eprob,
           .lsize = lsize, .esize = esize,
           .mds.min = mds.min))),


  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1

    pmat  <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                       .lprob , earg = .eprob )
    temp3 <-           eta[, c(FALSE, TRUE), drop = FALSE]
    if ( .lsize == "loge") {
      bigval <- 68
      temp3[temp3 >  bigval] <-  bigval  # pmin() collapses matrices
      temp3[temp3 < -bigval] <- -bigval
     }
    kmat <- as.matrix(eta2theta(temp3, .lsize , earg = .esize ))

    dl.dprob <- kmat / pmat - y / (1.0 - pmat)
    dl.dkayy <- digamma(y + kmat) - digamma(kmat) + log(pmat)

    dprob.deta <- dtheta.deta(pmat, .lprob , earg = .eprob )
    dkayy.deta <- dtheta.deta(kmat, .lsize , earg = .esize )

    myderiv <- c(w) * cbind(dl.dprob * dprob.deta,
                            dl.dkayy * dkayy.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, M + M - 1)  # wz is 'tridiagonal'




    max.support <- .max.support
    max.chunk.MB <- .max.chunk.MB


    ind2 <- matrix(FALSE, n, NOS)  # Used for SFS
    for (jay in 1:NOS) {
      eff.p <- sort(c( .cutoff.prob , 1 - .cutoff.prob ))
      Q.mins <- 0
      Q.maxs <- round(qnbinom(p = eff.p[2],
                              mu = mu[, jay],
                              size = kmat[, jay]) * 1.1) + 30



      eps.trig <- .eps.trig
      Q.MAXS <-      pmax(10, ceiling(1 / sqrt(eps.trig)))
      Q.maxs <- pmin(Q.maxs, Q.MAXS)


      ind1 <- if (max.chunk.MB > 0)
        (Q.maxs - Q.mins < max.support) else FALSE
      if ((NN <- sum(ind1)) > 0) {
        Object.Size <- NN * 8 * max(Q.maxs - Q.mins) / (2^20)
        n.chunks <- if (intercept.only) 1 else
                    max(1, ceiling( Object.Size / max.chunk.MB))
        chunk.rows <- ceiling(NN / n.chunks)
        ind2[, jay] <- ind1  # Save this
        wind2 <- which(ind1)


        upr.ptr <- 0
        lwr.ptr <- upr.ptr + 1
        while (lwr.ptr <= NN) {
          upr.ptr <- min(upr.ptr + chunk.rows, NN)
          sind2 <- wind2[lwr.ptr:upr.ptr]
          wz[sind2, M1*jay] <-
            EIM.NB.specialp(mu          =   mu[sind2, jay],
                            size        = kmat[sind2, jay],
                            y.max = max(Q.maxs[sind2]),
                            cutoff.prob = .cutoff.prob ,
                            intercept.only = intercept.only,
                            extra.bit = FALSE)
          lwr.ptr <- upr.ptr + 1
        }  # while
      }  # if
    }  # end of for (jay in 1:NOS)









    for (jay in 1:NOS) {
      run.varcov <- 0
      ii.TF <- !ind2[, jay]  # Not assigned above
      if (any(ii.TF)) {
        ppvec <- pmat[ii.TF, jay]
        kkvec <- kmat[ii.TF, jay]
        muvec <-   mu[ii.TF, jay]
        for (ii in 1:( .nsimEIM )) {
          ysim <- rnbinom(sum(ii.TF), mu = muvec, size = kkvec)
          dl.dk <- digamma(ysim + kkvec) - digamma(kkvec) + log(ppvec)
          run.varcov <- run.varcov + dl.dk^2
        }  # end of for loop

        run.varcov <- c(run.varcov / .nsimEIM )
        ned2l.dk2 <- if (intercept.only) mean(run.varcov) else run.varcov

        wz[ii.TF, M1*jay] <- ned2l.dk2  # * (dk.deta2[ii.TF, jay])^2
      }
    }


    wz[,     M1*(1:NOS)    ] <- wz[,      M1 * (1:NOS)] * dkayy.deta^2


    save.weights <- !all(ind2)


    ned2l.dprob2 <- kmat / ((1 - pmat) * pmat^2)
    wz[,     M1*(1:NOS) - 1] <- ned2l.dprob2 * dprob.deta^2

    ned2l.dkayyprob <- -1 / pmat
    wz[, M + M1*(1:NOS) - 1] <- ned2l.dkayyprob * dkayy.deta * dprob.deta




    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }), list( .cutoff.prob = cutoff.prob, .eps.trig = eps.trig,
            .max.support = max.support,
            .max.chunk.MB = max.chunk.MB,
            .nsimEIM = nsimEIM ))))




  if (deviance.arg)
  ans@deviance <- eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    temp300 <-  eta[, c(FALSE, TRUE), drop = FALSE]


    if (ncol(as.matrix(y)) > 1 && ncol(as.matrix(w)) > 1)
      stop("cannot handle matrix 'w' yet")



    if ( .lsize == "loge") {
      bigval <- 68
      temp300[temp300 >  bigval] <-  bigval
      temp300[temp300 < -bigval] <- -bigval
    } else {
      stop("can only handle the 'loge' link")
    }
    kayy <-  eta2theta(temp300, .lsize , earg = .esize)
    devi <- 2 * (y * log(ifelse(y < 1, 1, y) / mu) +
                (y + kayy) * log((mu + kayy) / (kayy + y)))
    if (residuals) {
      sign(y - mu) * sqrt(abs(devi) * w)
    } else {
      dev.elts <- sum(c(w) * devi)
      if (summation) {
        sum(dev.elts)
      } else {
        dev.elts
      }
    }
  }, list( .lsize = lsize, .eprob = eprob,
           .esize = esize )))

  ans
}  # End of polya()











polyaR.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}



 polyaR <-
  function(
           zero = "size",
           type.fitted = c("mean", "prob"),
           mds.min = 1e-3,
           nsimEIM = 500,  cutoff.prob = 0.999,  # Maxiter = 5000,
           eps.trig = 1e-7,
           max.support = 4000,
           max.chunk.MB = 30,  # max.memory = Inf is allowed
           lsize = "loge", lprob = "logit",
           imethod = 1,
           iprob = NULL,
           iprobs.y = NULL,
           gprobs.y = ppoints(6),
           isize = NULL,
           gsize.mux = exp(c(-30, -20, -15, -10, -6:3)),
           imunb = NULL) {


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
      stop("argument 'imethod' must be 1 or 2")


  deviance.arg <- FALSE  # 20131212; for now


  type.fitted <- match.arg(type.fitted,
                           c("mean", "prob"))[1]


  if (!is.Numeric(eps.trig, length.arg = 1,
                  positive = TRUE) || eps.trig > 0.001)
    stop("argument 'eps.trig' must be positive and smaller in value")


  if (length(iprob) && !is.Numeric(iprob, positive = TRUE))
    stop("bad input for argument 'iprob'")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("bad input for argument 'isize'")

  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE))
    stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 10)
    warning("argument 'nsimEIM' should be an integer ",
            "greater than 10, say")


  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")



  ans <-
  new("vglmff",
  blurb = c("Polya (negative-binomial) distribution\n\n",
            "Links:    ",
            namesof("size", lsize, earg = esize), ", ",
            namesof("prob", lprob, earg = eprob), "\n",
            "Mean:     size * (1 - prob) / prob\n",
            "Variance: mean / prob"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         mds.min = .mds.min ,
         multipleResponses = TRUE,
         type.fitted  = .type.fitted ,
         parameters.names = c("size", "prob"),
         eps.trig = .eps.trig ,
         zero = .zero )
  }, list( .zero = zero, .eps.trig = eps.trig,
           .type.fitted = type.fitted,
           .mds.min = mds.min))),

  initialize = eval(substitute(expression({
    M1 <- 2
    if (any(function.name == c("cqo", "cao")))
      stop("polyaR() does not work with cqo() or cao(). ",
           "Try negbinomial()")


    temp12 <- w.y.check(w = w, y = y,
              Is.integer.y = TRUE,
              Is.nonnegative = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1, maximize = TRUE)
    w <- temp12$w
    y <- temp12$y


    M <- M1 * ncol(y)
    NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    predictors.names <-
      c(namesof(param.names("size", NOS),
                .lsize , earg = .esize , tag = FALSE),
        namesof(param.names("prob", NOS),
                .lprob , earg = .eprob , tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]

    if (is.null( .nsimEIM )) {
       save.weights <- control$save.weights <- FALSE
    }



    gprobs.y <- .gprobs.y
    imunb <- .imunb  # Default in NULL
    if (length(imunb))
      imunb <- matrix(imunb, n, NOS, byrow = TRUE)



    if (!length(etastart)) {
      munb.init <-
      size.init <- matrix(NA_real_, n, NOS)
      gprobs.y  <- .gprobs.y
      if (length( .iprobs.y ))
        gprobs.y <-  .iprobs.y
      gsize.mux <- .gsize.mux  # gsize.mux is on a relative scale

      for (jay in 1:NOS) {  # For each response 'y_jay'... do:
        munb.init.jay <- if ( .imethod == 1 ) {
          quantile(y[, jay], probs = gprobs.y) + 1/16
        } else {
          weighted.mean(y[, jay], w = w[, jay])
        }
        if (length(imunb))
          munb.init.jay <- imunb[, jay]


        gsize <- gsize.mux * 0.5 * (mean(munb.init.jay) +
                                    weighted.mean(y[, jay], w = w[, jay]))
        if (length( .isize ))
          gsize <- .isize  # isize is on an absolute scale


        try.this <-
          grid.search2(munb.init.jay, gsize,
                       objfun = NBD.Loglikfun2,
                       y = y[, jay], w = w[, jay],
                       ret.objfun = TRUE)  # Last value is the loglik

        munb.init[, jay] <- try.this["Value1"]
        size.init[, jay] <- try.this["Value2"]
      }  # for (jay ...)






      prob.init <- if (length( .iprob ))
                   matrix( .iprob , nrow(y), ncol(y), byrow = TRUE) else
                   size.init / (size.init + munb.init)


      etastart <-
        cbind(theta2eta(size.init, .lsize , earg = .esize ),
              theta2eta(prob.init, .lprob , earg = .eprob ))
      etastart <-
        etastart[, interleave.VGAM(M, M1 = M1), drop = FALSE]
      }
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize,
            .iprob = iprob, .isize = isize,
            .pinit = iprob,
            .gprobs.y = gprobs.y, .gsize.mux = gsize.mux,
            .iprobs.y = iprobs.y,
            .nsimEIM = nsimEIM, .zero = zero,
            .imethod = imethod , .imunb = imunb,
            .type.fitted = type.fitted ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- ncol(eta) / c(M1 = 2)
    kmat <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                      .lsize , earg = .esize )
    pmat <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                      .lprob , earg = .eprob )

   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob"))[1]

    ans <- switch(type.fitted,
                  "mean"      = kmat * (1 - pmat) / pmat,
                  "prob"      = pmat)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lprob = lprob, .eprob = eprob,
           .lsize = lsize, .esize = esize))),
  last = eval(substitute(expression({
    temp0303 <- c(rep_len( .lprob , NOS),
                  rep_len( .lsize , NOS))
    names(temp0303) <- c(param.names("size", NOS),
                         param.names("prob", NOS))
    temp0303 <- temp0303[interleave.VGAM(M, M1 = M1)]
    misc$link <- temp0303  # Already named

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .esize
      misc$earg[[M1*ii  ]] <- .eprob
    }

    misc$isize <- .isize
    misc$imethod <- .imethod
    misc$nsimEIM <- .nsimEIM
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize,
            .isize = isize,
            .nsimEIM = nsimEIM,
            .imethod = imethod ))),


  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pmat  <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                       .lprob , earg = .eprob)
    temp300 <-         eta[, c(TRUE, FALSE), drop = FALSE]
    if ( .lsize == "loge") {
      bigval <- 68
      temp300[temp300 >  bigval] <-  bigval
      temp300[temp300 < -bigval] <- -bigval
    }
    kmat <- eta2theta(temp300, .lsize , earg = .esize )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dnbinom(y, prob = pmat, size = kmat, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lsize = lsize, .lprob = lprob,
           .esize = esize, .eprob = eprob ))),
  vfamily = c("polyaR"),



  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    kmat <- eta2theta(eta[, c(TRUE, FALSE)], .lsize , .esize )
    pmat <- eta2theta(eta[, c(FALSE, TRUE)], .lprob , .eprob )
    rnbinom(nsim * length(pmat), prob = pmat, size = kmat)
  }, list( .lprob = lprob, .lsize = lsize,
           .eprob = eprob, .esize = esize ))),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    size <- eta2theta(eta[, c(TRUE, FALSE)], .lsize , .esize )
    pmat <- eta2theta(eta[, c(FALSE, TRUE)], .lprob , .eprob )
    munb <- size * (1 / pmat - 1)

    smallval <- .mds.min  # .munb.div.size
    overdispersion <- all(munb / size > smallval)
    ans <- all(is.finite(munb)) && all(0 < munb) &&
           all(is.finite(size)) && all(0 < size) &&
           all(is.finite(pmat)) && all(0 < pmat & pmat < 1) &&
           overdispersion
    if (!overdispersion)
      warning("parameter 'size' has very large values; ",
              "try fitting a quasi-Poisson ",
              "model instead.")
    ans
  }, list( .lprob = lprob, .eprob = eprob,
           .lsize = lsize, .esize = esize,
           .mds.min = mds.min))),


  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1

    pmat  <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                       .lprob , earg = .eprob)
    temp3 <-           eta[, c(TRUE, FALSE), drop = FALSE]
    if ( .lsize == "loge") {
      bigval <- 68
      temp3[temp3 >  bigval] <-  bigval  # pmin() collapses matrices
      temp3[temp3 < -bigval] <- -bigval
     }
    kmat <- as.matrix(eta2theta(temp3, .lsize , earg = .esize ))

    dl.dprob <- kmat / pmat - y / (1.0 - pmat)
    dl.dkayy <- digamma(y + kmat) - digamma(kmat) + log(pmat)

    dprob.deta <- dtheta.deta(pmat, .lprob , earg = .eprob)
    dkayy.deta <- dtheta.deta(kmat, .lsize , earg = .esize)

    myderiv <- c(w) * cbind(dl.dkayy * dkayy.deta,
                            dl.dprob * dprob.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize))),
  weight = eval(substitute(expression({
    wz <- matrix(0.0, n, M + M - 1)  # wz is 'tridiagonal'




    max.support <- .max.support
    max.chunk.MB <- .max.chunk.MB


    ind2 <- matrix(FALSE, n, NOS)  # Used for SFS
    for (jay in 1:NOS) {
      eff.p <- sort(c( .cutoff.prob , 1 - .cutoff.prob ))
      Q.mins <- 0
      Q.maxs <- round(qnbinom(p = eff.p[2],
                              mu = mu[, jay],
                              size = kmat[, jay]) * 1.1) + 30



      eps.trig <- .eps.trig
      Q.MAXS <-      pmax(10, ceiling(1 / sqrt(eps.trig) - kmat[, jay]))
      Q.maxs <- pmin(Q.maxs, Q.MAXS)



      ind1 <- if (max.chunk.MB > 0)
        (Q.maxs - Q.mins < max.support) else FALSE
      if ((NN <- sum(ind1)) > 0) {
        Object.Size <- NN * 8 * max(Q.maxs - Q.mins) / (2^20)
        n.chunks <- if (intercept.only) 1 else
                    max(1, ceiling( Object.Size / max.chunk.MB))
        chunk.rows <- ceiling(NN / n.chunks)
        ind2[, jay] <- ind1  # Save this
        wind2 <- which(ind1)


        upr.ptr <- 0
        lwr.ptr <- upr.ptr + 1
        while (lwr.ptr <= NN) {
          upr.ptr <- min(upr.ptr + chunk.rows, NN)
          sind2 <- wind2[lwr.ptr:upr.ptr]
          wz[sind2, M1*jay - 1] <-
            EIM.NB.specialp(mu          =   mu[sind2, jay],
                            size        = kmat[sind2, jay],
                            y.max = max(Q.maxs[sind2]),
                            cutoff.prob = .cutoff.prob ,
                            intercept.only = intercept.only,
                            extra.bit = FALSE)
          lwr.ptr <- upr.ptr + 1
        }  # while
      }  # if
    }  # end of for (jay in 1:NOS)









    for (jay in 1:NOS) {
      run.varcov <- 0
      ii.TF <- !ind2[, jay]  # Not assigned above
      if (any(ii.TF)) {
        ppvec <- pmat[ii.TF, jay]
        kkvec <- kmat[ii.TF, jay]
        muvec <-   mu[ii.TF, jay]
        for (ii in 1:( .nsimEIM )) {
          ysim <- rnbinom(sum(ii.TF), mu = muvec, size = kkvec)
          dl.dk <- digamma(ysim + kkvec) - digamma(kkvec) + log(ppvec)
          run.varcov <- run.varcov + dl.dk^2
        }  # end of for loop

        run.varcov <- c(run.varcov / .nsimEIM )
        ned2l.dk2 <- if (intercept.only) mean(run.varcov) else run.varcov

        wz[ii.TF, M1*jay - 1] <- ned2l.dk2  # * (dk.deta2[ii.TF, jay])^2
      }
    }


    wz[, M1*(1:NOS) - 1] <- wz[, M1*(1:NOS) - 1] * dkayy.deta^2


    save.weights <- !all(ind2)


    ned2l.dprob2 <- kmat / ((1 - pmat) * pmat^2)
    wz[,     M1*(1:NOS)    ] <- ned2l.dprob2 * dprob.deta^2

    ned2l.dkayyprob <- -1 / pmat
    wz[, M + M1*(1:NOS) - 1] <- ned2l.dkayyprob * dkayy.deta * dprob.deta




    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }), list( .cutoff.prob = cutoff.prob, .eps.trig = eps.trig,
            .max.support = max.support,
            .max.chunk.MB = max.chunk.MB,
            .nsimEIM = nsimEIM ))))




  if (deviance.arg)
  ans@deviance <- eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    temp300 <-  eta[, c(FALSE, TRUE), drop = FALSE]


    if (ncol(as.matrix(y)) > 1 && ncol(as.matrix(w)) > 1)
      stop("cannot handle matrix 'w' yet")



    if ( .lsize == "loge") {
      bigval <- 68
      temp300[temp300 >  bigval] <-  bigval
      temp300[temp300 < -bigval] <- -bigval
    } else {
      stop("can only handle the 'loge' link")
    }
    kayy <-  eta2theta(temp300, .lsize , earg = .esize)
    devi <- 2 * (y * log(ifelse(y < 1, 1, y) / mu) +
                (y + kayy) * log((mu + kayy) / (kayy + y)))
    if (residuals) {
      sign(y - mu) * sqrt(abs(devi) * w)
    } else {
      dev.elts <- sum(c(w) * devi)
      if (summation) {
        sum(dev.elts)
      } else {
        dev.elts
      }
    }
  }, list( .lsize = lsize, .eprob = eprob,
           .esize = esize )))

  ans
}  # End of polyaR()










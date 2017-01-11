# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.















dlog <- function(x, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(shape))
  if (length(x)     != N) x     <- rep_len(x,     N)
  if (length(shape) != N) shape <- rep_len(shape, N)
  ox <- !is.finite(x)
  zero <- ox | round(x) != x | x < 1
  ans <- rep_len(0.0, length(x))
  if (log.arg) {
    ans[ zero] <- log(0.0)
    ans[!zero] <- x[!zero] * log(shape[!zero]) - log(x[!zero]) -
                  log(-log1p(-shape[!zero]))
    ans[ox] <- log(0)  # 20141212 KaiH
  } else {
    ans[!zero] <- -(shape[!zero]^(x[!zero])) / (x[!zero] *
                   log1p(-shape[!zero]))
    ans[ox] <- 0.0  # 20141212 KaiH
  }
  ans[shape < 0 | 1 < shape] <- NaN
  ans
}



plog  <- function(q, shape, log.p = FALSE) {


  if (any(is.na(q))) stop("NAs not allowed for argument 'q'")
  if (any(is.na(shape)))
    stop("NAs not allowed for argument 'shape'")


  N <- max(length(q), length(shape))
  if (length(q)     != N) q     <- rep_len(q,     N)
  if (length(shape) != N) shape <- rep_len(shape, N)


  bigno <- 10
  owen1965 <- (q * (1 - shape) > bigno)
  if (specialCase <- any(owen1965)) {
    qqq <- q[owen1965]
    ppp <- shape[owen1965]
    pqp <- qqq * (1 - ppp)
    bigans <- (ppp^(1+qqq) / (1-ppp)) * (1/qqq -
              1 / (            pqp * (qqq-1)) +
              2 / ((1-ppp)   * pqp * (qqq-1) * (qqq-2)) -
              6 / ((1-ppp)^2 * pqp * (qqq-1) * (qqq-2) * (qqq-3)) +
          24 / ((1-ppp)^3 * pqp * (qqq-1) * (qqq-2) * (qqq-3) * (qqq-4)))
      bigans <- 1 + bigans / log1p(-ppp)
  }

  floorq <- pmax(1, floor(q))  # Ensures at least one element per q value
  floorq[owen1965] <- 1
  seqq <- sequence(floorq)
  seqp <- rep(shape, floorq)
  onevector <- (seqp^seqq / seqq) / (-log1p(-seqp))
  rlist <-  .C("tyee_C_cum8sum",
                as.double(onevector), answer = double(N),
                as.integer(N), as.double(seqq),
                as.integer(length(onevector)), notok=integer(1))
  if (rlist$notok != 0)
    stop("error in C function 'cum8sum'")
  ans <- if (log.p) log(rlist$answer) else rlist$answer
  if (specialCase)
    ans[owen1965] <- if (log.p) log(bigans) else bigans
  ans[q < 1] <- if (log.p) log(0.0) else 0.0
  ans[shape < 0 | 1 < shape] <- NaN
  ans
}





 qlog <- function(p, shape) {

  LLL <- max(length(p), length(shape))
  if (length(p)     != LLL) p     <- rep_len(p,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  ans <- rep_len(0, LLL)

  lo <- rep_len(1, LLL)
  approx.ans <- lo  # True at lhs
  hi <- 2 * lo + 10
  dont.iterate <- p == 1 | shape <= 0
  done <- p <= plog(hi, shape) | dont.iterate
  while (!all(done)) {
    hi.save <- hi[!done]
    hi[!done] <- 2 * lo[!done] + 10
    lo[!done] <- hi.save
    done[!done] <- (p[!done] <= plog(hi[!done], shape[!done]))
  }

  foo <- function(q, shape, p)
    plog(q, shape) - p

  lhs <- (p <= dlog(1, shape)) | dont.iterate

  approx.ans[!lhs] <-
    bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                    shape = shape[!lhs], p = p[!lhs])
  faa <- floor(approx.ans)
  ans <- ifelse(plog(faa, shape) < p & p <= plog(faa+1, shape),
                faa+1, faa)

  ans[p == 1] <- Inf
  ans[shape <= 0] <- NaN

  ans
}  # qlog




rlog <- function(n, shape) {
  qlog(runif(n), shape)
}







 logff <- function(lshape = "logit", gshape = ppoints(8), zero = NULL) {

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  new("vglmff",
  blurb = c("Logarithmic distribution f(y) = a * shape^y / y, ",
             "y = 1, 2, 3,...,\n",
             "            0 < shape < 1, a = -1 / log(1-shape)  \n\n",
             "Link:    ", namesof("shape", lshape, earg = eshape),
             "\n", "\n",
             "Mean:    a * shape / (1 - shape)", "\n"),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 1
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = "shape",
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    M <- M1 * ncoly


    mynames1  <- param.names("shape", ncoly)
    predictors.names <- namesof(mynames1, .lshape , earg = .eshape ,
                                tag = FALSE)


    if (!length(etastart)) {
      logff.Loglikfun <- function(shapeval, y, x, w, extraargs) {
        sum(c(w) * dlog(x = y, shape = shapeval, log = TRUE))
      }
      Init.shape <- matrix(0, n, M)
      shape.grid <- .gshape

      for (ilocal in 1:ncoly) {
        Init.shape[, ilocal] <- grid.search(shape.grid,
                                            objfun = logff.Loglikfun,
                                            y = y[, ilocal],  # x = x,
                                            w = w[, ilocal])
      }  # for
      etastart <- theta2eta(Init.shape, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .eshape = eshape, .gshape = gshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    aa <- -1 / log1p(-shape)
    aa * shape / (1 - shape)
  }, list( .lshape = lshape, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link <- c(rep_len( .lshape , ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .eshape
    }
  }), list( .lshape = lshape, .eshape = eshape ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dlog(x = y, shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("logff"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    okay0 <- if ( .lshape == "logfflink") all(0 < eta) else TRUE
    okay1 <- if (okay0) {
      shape <- eta2theta(eta, .lshape , earg = .eshape )
      all(is.finite(shape)) && all(0 < shape & shape < 1)
    } else {
      FALSE
    }
    okay0 && okay1
  }, list( .lshape = lshape, .eshape = eshape ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    rlog(nsim * length(shape), shape = shape)
  }, list( .lshape = lshape, .eshape = eshape ))),



  deriv = eval(substitute(expression({
    M1 <- 1
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    aa <- -1 / log1p(-shape)
    dl.dshape <- -aa / (1 - shape) + y / shape
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = eval(substitute(expression({
    ned2l.dshape2 <- aa * (1 - aa * shape) / (shape * (1-shape)^2)
    wz <- c(w) * ned2l.dshape2 * dshape.deta^2
    wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}









deflat.limit.oilog  <- function(shape) {
  if (any(shape <= 0 | 1 <= shape ))
    stop("argument 'shape' must be in (0, 1)")
  ans <- 1 / (1 - 1 / dlog(1, shape))
  ans
}



doilog <- function(x, shape, pstr1 = 0, log = FALSE) {



  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(pstr1))
  if (length(x)     != LLL) x     <- rep_len(x,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)

  ans <- rep(NA_real_, LLL)
  index1 <- (x == 1)
  if (log.arg) {
    ans[ index1] <- log(pstr1[ index1] + (1 - pstr1[ index1]) *
                       dlog(x[ index1], shape[ index1]))
    ans[!index1] <- log1p(-pstr1[!index1]) +
                          dlog(x[!index1], shape[!index1], log = TRUE)
  } else {
    ans[ index1] <-      pstr1[ index1] + (1 - pstr1[ index1]) *
                        dlog(x[ index1], shape[ index1])
    ans[!index1] <- (1 - pstr1[!index1]) *
                        dlog(x[!index1], shape[!index1])
  }


  ans[pstr1 < deflat.limit.oilog(shape) | 1 < pstr1] <- NaN
  ans[shape <= 0 | 1 <= shape] <- NaN
  ans
}  # doilog




poilog <- function(q, shape, pstr1 = 0) {

  LLL <- max(length(q), length(shape), length(pstr1))
  if (length(q)     != LLL) q     <- rep_len(q,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oilog(shape)

  ans <- plog(q, shape)  #, lower.tail = lower.tail, log.p = log.p
  ans <- ifelse(q < 1, 0, pstr1 + (1 - pstr1) * ans)

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN
  ans[shape <= 0] <- NaN
  ans[1 <= shape] <- NaN

  ans
}  # poilog





qoilog <- function(p, shape, pstr1 = 0) {

  LLL <- max(length(p), length(shape), length(pstr1))
  if (length(p)     != LLL) p     <- rep_len(p,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oilog(shape)

  ans[p <= pstr1] <- 1
  pindex <- (deflat.limit <= pstr1) & (pstr1 < p)
  ans[pindex] <-
    qlog((p[pindex] - pstr1[pindex]) / (1 - pstr1[pindex]),
         shape = shape[pindex])

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN

  ans[p < 0] <- NaN
  ans[1 < p] <- NaN
  ans[shape <= 0] <- NaN
  ans[1 <= shape] <- NaN

  ans
}  # qoilog



roilog <- function(n, shape, pstr1 = 0) {
  qoilog(runif(n), shape, pstr1 = pstr1)
}





 oilog <-
  function(lpstr1 = "logit", lshape = "logit",
         type.fitted = c("mean", "shape", "pobs1", "pstr1", "onempstr1"),
           ishape = NULL,
           gpstr1 = ppoints(8),
           gshape = ppoints(8),
           zero = NULL) {

  lpstr1 <- as.list(substitute(lpstr1))
  epstr1 <- link2list(lpstr1)
  lpstr1 <- attr(epstr1, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")



  type.fitted <- match.arg(type.fitted,
                   c("mean", "shape", "pobs1", "pstr1", "onempstr1"))[1]


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' values must be positive")


  new("vglmff",
  blurb = c("One-inflated logarithmic distribution\n\n",
            "Links:    ",
            namesof("pstr1", lpstr1, earg = epstr1 ), ", ",
            namesof("shape", lshape, earg = eshape ), "\n",
            "Mean:     pstr1 + (1 - pstr1) * a * shape / (1 - shape), ",
                       "a = -1 / log(1-shape), 0 < shape < 1"),

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
         parameters.names = c("pstr1", "shape"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),
  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    NOS <- ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("pstr1", ncoly)
    mynames2 <- param.names("shape", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lpstr1 , earg = .epstr1 , tag = FALSE),
          namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]


    if (!length(etastart)) {

      shape.init <-
      pstr1.init <- matrix(NA_real_, n, NOS)
      gpstr1 <- .gpstr1
      gshape <- .gshape

      oilog.Loglikfun <- function(pstr1, shape, y, x, w, extraargs) {
        sum(c(w) * doilog(x = y, pstr1 = pstr1,
                          shape = shape, log = TRUE))
      }


      for (jay in 1:NOS) {  # For each response 'y_jay'... do:


        try.this <-
          grid.search2(gpstr1, gshape,
                       objfun = oilog.Loglikfun,
                       y = y[, jay],  # x = x[TFvec, , drop = FALSE],
                       w = w[, jay],
                       ret.objfun = TRUE)  # Last value is the loglik
        pstr1.init[, jay] <-  try.this["Value1"]
        shape.init[, jay] <-  try.this["Value2"]
      }  # for (jay ...)

      etastart <- cbind(theta2eta(pstr1.init, .lpstr1 , earg = .epstr1 ),
                        theta2eta(shape.init, .lshape , earg = .eshape ))[,
                        interleave.VGAM(M, M1 = M1)]
      mustart <- NULL  # Since etastart has been computed.
    }  # End of !length(etastart)
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape,
                              .ishape = ishape,
            .gpstr1 = gpstr1,
            .gshape  = gshape,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- ncol(eta) / c(M1 = 2)
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "shape", "pobs1", "pstr1", "onempstr1"))[1]

    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )

    Meanfun <- function(shape) {
      aa <- -1 / log1p(-shape)
      Mean <- aa * shape / (1 - shape)
      Mean[shape <= 0 | 1 <= shape] <- NaN
      Mean
    }

    ans <-
      switch(type.fitted,
             "mean"      = pstr1 + (1 - pstr1) * Meanfun(shape),
             "shape"     = shape,
             "pobs1" = doizeta(1, shape = shape, pstr1 = pstr1),  # P(Y=1)
             "pstr1"     =     pstr1,
             "onempstr1" = 1 - pstr1)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$link <-
      c(rep_len( .lpstr1 , NOS),
        rep_len( .lshape , NOS))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .epstr1
      misc$earg[[M1*ii  ]] <- .eshape
    }
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doilog(x = y, pstr1 = pstr1, shape = shape,
                               log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),
  vfamily = c("oilog"),



  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    roilog(nsim * length(shape), shape = shape, pstr1 = pstr1)
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                       earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lshape  ,
                       earg = .eshape )
    okay1 <- all(is.finite(shape )) && all(0 < shape & shape < 1) &&
             all(is.finite(pstr1)) && all(pstr1 < 1)
    deflat.limit <- deflat.limit.oizeta(shape)
    okay2.deflat <- TRUE
    if (okay1 && !(okay2.deflat <- all(deflat.limit < pstr1)))
      warning("parameter 'pstr1' is too negative even allowing for ",
              "1-deflation.")
    okay1 && okay2.deflat
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),






  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- M / M1
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                       earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lshape ,
                       earg = .eshape )

    pmf1 <- dlog(1, shape)
    onempmf1 <- 1 - pmf1  # dozeta(1, shape = shape, pstr1 = pstr1)
    pobs1 <- pstr1 + (1 - pstr1) * pmf1
    index1 <- as.matrix(y == 1)

    mraa <- log1p(-shape)
    aaaa <- -1 / mraa

    dl.dpstr1 <- onempmf1 / pobs1
    dl.dpstr1[!index1] <- -1 / (1 - pstr1[!index1])




    dpmf1.dshape <- -1 / mraa - shape / ((1 - shape) * mraa^2)
    d2pmf1.dshape2 <- -2 / ((1 - shape) * mraa^2) -
                      shape * (2 + mraa) / ((1 - shape)^2 * mraa^3)

    dl.dshape <- (1 - pstr1) * dpmf1.dshape / pobs1  #
    dl.dshape[!index1] <- y[!index1] / shape[!index1] +
                         1 / ((1 - shape[!index1]) * mraa[!index1])

    dpstr1.deta <- dtheta.deta(pstr1, .lpstr1 , earg = .epstr1 )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    myderiv <- c(w) * cbind(dl.dpstr1 * dpstr1.deta,
                            dl.dshape * dshape.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape ))),
  weight = eval(substitute(expression({

    EY.y.gt.1 <- aaaa * shape^2 / ((1 - shape) * (1 - aaaa * shape))
    LHS <- ((1 - pstr1) / pobs1) * dpmf1.dshape^2 - d2pmf1.dshape2
    RHS <- EY.y.gt.1 / shape^2 - (1 + mraa) / ((1 - shape) * mraa)^2
    ned2l.dpstr12 <- onempmf1 / ((1 - pstr1) * pobs1)  #
    ned2l.dpstr1shape <- dpmf1.dshape / pobs1  #
    ned2l.dshape2 <- (1 - pstr1) * (LHS + (1 - pmf1) * RHS)

    wz <- array(c(c(w) * ned2l.dpstr12 * dpstr1.deta^2,
                  c(w) * ned2l.dshape2 * dshape.deta^2,
                  c(w) * ned2l.dpstr1shape * dpstr1.deta * dshape.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}  # oilog






dotlog <- function(x, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (log.arg) {
    ans <- dlog(x, shape, log = log.arg) - log1p(-dlog(1, shape))
    ans[x == 1] <- log(0)
  } else {
    ans <- dlog(x, shape) / (1 - dlog(1, shape))
    ans[x == 1] <- 0
  }
  ans
}  # dotlog



potlog  <- function(q, shape, log.p = FALSE) {
  if (log.p) log(plog(q, shape) - dlog(1, shape)) -
      log1p(-dlog(1, shape)) else
    (plog(q, shape) - dlog(1, shape)) / (1 - dlog(1, shape))
}



 qotlog <- function(p, shape) {

  ans <- qlog((1 - dlog(1, shape)) * p + dlog(1, shape), shape = shape)

  ans[p == 1] <- Inf
  ans[p < 0] <- NaN
  ans[1 < p] <- NaN

  ans[shape < 0 | 1 < shape] <- NaN
  ans
}  # qotlog



rotlog <- function(n, shape) {
  qotlog(runif(n), shape)
}



 otlog <- function(lshape = "logit", gshape = ppoints(8), zero = NULL) {

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")



  new("vglmff",
  blurb = c("One-truncated logarithmic distribution ",
            "f(y) = shape^y / ((-shape - log1p(-shape)) * y), ",
             "y = 2, 3,...,\n",
             "            0 < shape < 1,\n\n",
             "Link:    ", namesof("shape", lshape, earg = eshape)),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 1
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = "shape",
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (any(y <= 1))
      stop("cannot have any 1s in the response")


    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    M <- M1 * ncoly


    mynames1  <- param.names("shape", ncoly)
    predictors.names <- namesof(mynames1, .lshape , earg = .eshape ,
                                tag = FALSE)


    if (!length(etastart)) {
      dotlog.Loglikfun <- function(shapeval, y, x, w, extraargs) {
        sum(c(w) * dotlog(x = y, shape = shapeval, log = TRUE))
      }
      Init.shape <- matrix(0, n, M)
      shape.grid <- .gshape

      for (ilocal in 1:ncoly) {
        Init.shape[, ilocal] <- grid.search(shape.grid,
                                            objfun = dotlog.Loglikfun,
                                            y = y[, ilocal],  # x = x,
                                            w = w[, ilocal])
      }  # for
      etastart <- theta2eta(Init.shape, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .eshape = eshape, .gshape = gshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    aa <- -1 / log1p(-shape)
    ((aa * shape / (1 - shape)) - dlog(1, shape)) / (1 - dlog(1, shape))
  }, list( .lshape = lshape, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link <- c(rep_len( .lshape , ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .eshape
    }
  }), list( .lshape = lshape, .eshape = eshape ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dotlog(x = y, shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("otlog"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape & shape < 1)
    okay1
  }, list( .lshape = lshape, .eshape = eshape ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    rotlog(nsim * length(shape), shape = shape)
  }, list( .lshape = lshape, .eshape = eshape ))),



  deriv = eval(substitute(expression({
    M1 <- 1
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    aa <- -1 / log1p(-shape)
    dl.dshape <- y / shape +
                shape / ((1 - shape) * (shape + log1p(-shape)))
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = eval(substitute(expression({
    EY.logff <-  aa * shape / (1 - shape)

    d3 <- deriv3( ~ shape / ((1 - shape) * (shape + log(1 - shape))),
                  c("shape"), hessian = FALSE)
    eval.d3 <- eval(d3)
    d2pmf1.dshape2 <- c(attr(eval.d3, "gradient"))

    ned2l.dshape2 <-
      (EY.logff - dlog(1, shape)) / ((1 - dlog(1, shape)) * shape^2) -
      d2pmf1.dshape2
    wz <- c(w) * ned2l.dshape2 * dshape.deta^2
    wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}  # otlog






dotpospois <- function(x, lambda, log = FALSE) {
  if (!is.logical(larg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (larg) {
    ans <- dpospois(x, lambda, log = larg) - log1p(-dpospois(1, lambda))
    ans[x == 1] <- log(0)
  } else {
    ans <- dpospois(x, lambda) / (1 - dpospois(1, lambda))
    ans[x == 1] <- 0
  }
  ans
}  # dotpospois



potpospois  <- function(q, lambda, log.p = FALSE) {
  if (log.p) log(ppospois(q, lambda) - dpospois(1, lambda)) -
      log1p(-dpospois(1, lambda)) else
    (ppospois(q, lambda) - dpospois(1, lambda)) / (1-dpospois(1, lambda))
}



 qotpospois <- function(p, lambda) {
  ans <- qpospois((1 - dpospois(1, lambda)) * p +
                  dpospois(1, lambda), lambda = lambda)

  ans[p == 1 & 0 < lambda] <- Inf
  ans[p < 0] <- NaN
  ans[1 < p] <- NaN

  ans[lambda < 0] <- NaN
  ans
}  # qotpospois



rotpospois <- function(n, lambda) {
  qotpospois(runif(n), lambda)
}




 otpospoisson <-
    function(llambda = "loge",
             type.fitted = c("mean", "lambda", "prob0", "prob1"),
             ilambda = NULL, imethod = 1, zero = NULL) {

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")


  if (length( ilambda) && !is.Numeric(ilambda, positive = TRUE))
    stop("bad input for argument 'ilambda'")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "prob0", "prob1"))[1]


  new("vglmff",
  blurb = c("One-truncated Positive-Poisson distribution\n\n",
            "Links:    ",
            namesof("lambda", llambda, earg = elambda, tag = FALSE)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("lambda"),
         type.fitted  = .type.fitted ,
         llambda = .llambda ,
         elambda = .elambda )
  }, list( .llambda = llambda, .elambda = elambda,
           .type.fitted = type.fitted ))),

  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (any(y < 2))
      stop("response values must be 2 or more")

    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    M <- M1 * ncoly
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)


    mynames1 <- param.names("lambda", ncoly)
    predictors.names <- namesof(mynames1, .llambda , earg = .elambda ,
                                tag = FALSE)

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                             imu = .ilambda )

      etastart <- theta2eta(lambda.init, .llambda , earg = .elambda)
    }
  }), list( .llambda = llambda, .elambda = elambda,
            .ilambda = ilambda, .imethod = imethod,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- NCOL(eta) / c(M1 = 1)
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "prob0", "prob1"))[1]

    lambda <- eta2theta(eta, .llambda , earg = .elambda )
    ans <- switch(type.fitted,
                  "mean"      = lambda / ppois(1, lambda, lower = FALSE),
                  "lambda"    = lambda,
                  "prob0"     = ppois(0, lambda),  # P(Y=0) as it were
                  "prob1"     = ppois(1, lambda))  # P(Y=1) as it were
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .llambda = llambda, .elambda = elambda ))),
  last = eval(substitute(expression({
    misc$link <- rep_len( .llambda , M)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:M)
      misc$earg[[ii]] <- .elambda
  }), list( .llambda = llambda, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    lambda <- eta2theta(eta, .llambda , earg = .elambda )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dotpospois(x = y, lambda = lambda, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llambda = llambda, .elambda = elambda ))),
  vfamily = c("otpospoisson"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    lambda <- eta2theta(eta, .llambda , earg = .elambda )
    okay1 <- all(is.finite(lambda)) && all(0 < lambda)
    okay1
  }, list( .llambda = llambda, .elambda = elambda ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    lambda <- eta2theta(eta, .llambda , earg = .elambda )
    rotpospois(nsim * length(lambda), lambda)
  }, list( .llambda = llambda, .elambda = elambda ))),




  deriv = eval(substitute(expression({
    M1 <- 1
    lambda <- eta2theta(eta, .llambda , earg = .elambda )

    EY.cond <- 1 / ppois(1, lambda, lower.tail = FALSE)
    temp1 <- expm1(lambda)
    temp0 <- lambda * exp(-lambda)
    prob.geq.2 <- -expm1(-lambda) - temp0
    dl.dlambda <- y / lambda - 1 - temp0 / prob.geq.2

    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )

    c(w) * dl.dlambda * dlambda.deta
  }), list( .llambda = llambda, .elambda = elambda ))),
  weight = eval(substitute(expression({
    ned2l.dlambda2 <- EY.cond / lambda +
        ((1 - lambda) * exp(-lambda) - temp0^2 / prob.geq.2) / prob.geq.2
    wz <-  ned2l.dlambda2 * dlambda.deta^2
    c(w) * wz
  }), list( .llambda = llambda, .elambda = elambda ))))
}  # otpospoisson






doalog <- function(x, shape, pobs1 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(pobs1))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(shape)  != LLL) shape  <- rep_len(shape,  LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)
  ans <- rep_len(0.0, LLL)


  index1 <- (x == 1)

  if (log.arg) {
    ans[ index1] <- log(pobs1[index1])
    ans[!index1] <- log1p(-pobs1[!index1]) +
                    dotlog(x[!index1], shape[!index1], log = TRUE)
  } else {
    ans[ index1] <- pobs1[index1]
    ans[!index1] <- (1 - pobs1[!index1]) *
                    dotlog(x[!index1], shape[!index1])
  }
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans
}



poalog <- function(q, shape, pobs1 = 0) {
  LLL <- max(length(q), length(shape), length(pobs1))
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(shape)  != LLL) shape  <- rep_len(shape,  LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)
  ans <- rep_len(0.0, LLL)

  ans[q >  1] <-    pobs1[q > 1] +
                 (1-pobs1[q > 1]) * potlog(q[q > 1], shape[q > 1])
  ans[q <  1] <- 0
  ans[q == 1] <- pobs1[q == 1]

  ans <- pmax(0, ans)
  ans <- pmin(1, ans)
  ans[pobs1 < 0 | 1 < pobs1] <- NaN

  ans
}



qoalog <- function(p, shape, pobs1 = 0) {
  LLL <- max(length(p), length(shape), length(pobs1))
  if (length(p)      != LLL) p      <- rep_len(p,      LLL)
  if (length(shape)  != LLL) shape  <- rep_len(shape,  LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)

  ans <- rep_len(NaN, LLL)
  ind4 <- pobs1 < p
  ans[!ind4] <- 1
  ans[ ind4] <- qotlog((p[ind4] - pobs1[ind4]) / (1 - pobs1[ind4]),
                       shape = shape[ind4])
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[p < 0 | 1 < p] <- NaN
  ans
}



roalog <- function(n, shape, pobs1 = 0) {
  qoalog(runif(n), shape = shape, pobs1 = pobs1)
}






 oalog <-
  function(lpobs1 = "logit",
           lshape = "logit",
           type.fitted = c("mean", "shape", "pobs1", "onempobs1"),
           ipobs1 = NULL,
           gshape = ppoints(8),
           zero = NULL) {


  lpobs1 <- as.list(substitute(lpobs1))
  epobs1 <- link2list(lpobs1)
  lpobs1 <- attr(epobs1, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


 type.fitted <- match.arg(type.fitted,
                           c("mean", "shape", "pobs1", "onempobs1"))[1]


  new("vglmff",
  blurb = c("One-altered logarithmic distribution \n",
            "(Bernoulli and 1-truncated logarithmic distribution model)",
            "\n\n",
            "Links:    ",
            namesof("pobs1",  lpobs1, earg = epobs1, tag = FALSE), ", ",
            namesof("shape",  lshape, earg = eshape, tag = FALSE)),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("pobs1", "shape"),
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    M1 <- 2
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    extra$y1 <- y1 <- ifelse(y == 1, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y1), n, NOS)
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("pobs1", ncoly)
    mynames2 <- param.names("shape", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lpobs1 , earg = .epobs1 , tag = FALSE),
          namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly


    if (!length(etastart)) {
      dotlog.Loglikfun <- function(shapeval, y, x, w, extraargs) {
        sum(c(w) * dotlog(x = y, shape = shapeval, log = TRUE))
      }
      Init.shape <- matrix(0, n, ncoly)
      shape.grid <- .gshape

      for (jlocal in 1:ncoly) {
        index1 <- y[, jlocal] > 1
        Init.shape[, jlocal] <-
          grid.search(shape.grid,
                      objfun = dotlog.Loglikfun,
                      y = y[index1, jlocal],  # x = x,
                      w = w[index1, jlocal])
      }  # for
      etastart <-
        cbind(theta2eta(if (length( .ipobs1 )) .ipobs1 else
                        (0.5 + w * y1) / (1 + w),
                        .lpobs1 , earg = .epobs1 ),
              theta2eta(Init.shape, .lshape , earg = .eshape ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lshape = lshape, .eshape = eshape, .gshape = gshape,
            .lpobs1 = lpobs1, .epobs1 = epobs1,
            .ipobs1 = ipobs1,  # .ishape = ishape,
            .type.fitted = type.fitted
           ))),


  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "shape", "pobs1", "onempobs1"))[1]

    M1 <- 2
    NOS <- ncol(eta) / M1

    pobs1 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                             .lpobs1 , earg = .epobs1 ))
    shape <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                             .lshape , earg = .eshape ))

    aa <- -1 / log1p(-shape)
    otlog.mean <- ((aa * shape / (1 - shape)) -
                   dlog(1, shape)) / (1 - dlog(1, shape))

    ans <- switch(type.fitted,
                  "mean"      = pobs1 + (1 - pobs1) * otlog.mean,
                  "shape"     = shape,
                  "pobs1"     =      pobs1,  # P(Y=1)
                  "onempobs1" =  1 - pobs1)  # P(Y>1)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),
  last = eval(substitute(expression({
    temp.names <- c(rep_len( .lpobs1 , NOS),
                    rep_len( .lshape , NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epobs1
      misc$earg[[M1*ii  ]] <- .eshape
    }
  }), list( .lpobs1 = lpobs1, .lshape = lshape,
            .epobs1 = epobs1, .eshape = eshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    pobs1 <- cbind(eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                             .lpobs1, earg = .epobs1 ))
    shape <- cbind(eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                             .lshape, earg = .eshape ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doalog(x = y, pobs1 = pobs1, shape = shape,
                               log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),
  vfamily = c("oalog"),



  validparams = eval(substitute(function(eta, y, extra = NULL) {
    TFvec <- c(TRUE, FALSE)
    pobs1 <- eta2theta(eta[,  TFvec, drop = FALSE],
                       .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, !TFvec, drop = FALSE],
                       .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape & shape < 1) &&
             all(is.finite(pobs1)) && all(0 < pobs1 & pobs1 < 1)
    okay1
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    roalog(nsim * length(shape), shape = shape, pobs1 = pobs1)
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1  # extra$NOS
    y1 <- extra$y1
    skip <- extra$skip.these

    TFvec <- c(TRUE, FALSE)
    pobs1 <- eta2theta(eta[,  TFvec, drop = FALSE],
                       .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, !TFvec, drop = FALSE],
                       .lshape , earg = .eshape )

    aa <- -1 / log1p(-shape)
    dl.dshape <- y / shape +
                 shape / ((1 - shape) * (shape + log1p(-shape)))

    dl.dpobs1 <- -1 / (1 - pobs1)  # For y > 1 obsns

    for (spp. in 1:NOS) {
      dl.dpobs1[skip[, spp.], spp.] <- 1 / pobs1[skip[, spp.], spp.]
      dl.dshape[skip[, spp.], spp.] <- 0
    }
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    mu.phi1 <- pobs1

    temp3 <- if ( .lpobs1 == "logit") {
      c(w) * (y1 - mu.phi1)
    } else {
      c(w) * dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 ) *
             dl.dpobs1
    }

    ans <- cbind(temp3,
                 c(w) * dl.dshape * dshape.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lpobs1 = lpobs1, .lshape = lshape,
            .epobs1 = epobs1, .eshape = eshape ))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, M1 * NOS)  # EIM is diagonal


    EY.logff <-  aa * shape / (1 - shape)
    d3 <- deriv3( ~ shape / ((1 - shape) * (shape + log(1 - shape))),
                  c("shape"), hessian = FALSE)
    eval.d3 <- eval(d3)
    d2pmf1.dshape2 <- c(attr(eval.d3, "gradient"))

    ned2l.dshape2 <-
      (EY.logff - dlog(1, shape)) / ((1 - dlog(1, shape)) * shape^2) -
      d2pmf1.dshape2



    ned2l.dshape2 <- (1-pobs1) * ned2l.dshape2  #+stop("another quantity")
    wz[, NOS+(1:NOS)] <- c(w) * ned2l.dshape2 * dshape.deta^2


    tmp100 <- mu.phi1 * (1 - mu.phi1)
    tmp200 <- if ( .lpobs1 == "logit" && is.empty.list( .epobs1 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 )^2)
    }
    wz[, 1:NOS] <-  tmp200

    wz <- wz[, interleave.VGAM(ncol(wz), M1 = M1)]

    wz
  }), list( .lpobs1 = lpobs1,
            .epobs1 = epobs1 ))))
}  # End of oalog









doapospois <- function(x, lambda, pobs1 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(lambda), length(pobs1))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)
  ans <- rep_len(0.0, LLL)


  index1 <- (x == 1)

  if (log.arg) {
    ans[ index1] <- log(pobs1[index1])
    ans[!index1] <- log1p(-pobs1[!index1]) +
                    dotpospois(x[!index1], lambda[!index1], log = TRUE)
  } else {
    ans[ index1] <- pobs1[index1]
    ans[!index1] <- (1 - pobs1[!index1]) *
                    dotpospois(x[!index1], lambda[!index1])
  }
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[lambda < 0] <- NaN
  ans
}



poapospois <- function(q, lambda, pobs1 = 0) {
  LLL <- max(length(q), length(lambda), length(pobs1))
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)
  ans <- rep_len(0.0, LLL)

  ans[q >  1] <-    pobs1[q > 1] +
                 (1-pobs1[q > 1]) * potpospois(q[q > 1], lambda[q > 1])
  ans[q <  1] <- 0
  ans[q == 1] <- pobs1[q == 1]

  ans <- pmax(0, ans)
  ans <- pmin(1, ans)
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[lambda < 0] <- NaN

  ans
}



qoapospois <- function(p, lambda, pobs1 = 0) {
  LLL <- max(length(p), length(lambda), length(pobs1))
  if (length(p)      != LLL) p      <- rep_len(p,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)

  ans <- rep_len(NaN, LLL)
  ind4 <- pobs1 < p
  ans[!ind4] <- 1
  ans[ ind4] <- qotpospois((p[ind4] - pobs1[ind4]) / (1 - pobs1[ind4]),
                           lambda = lambda[ind4])
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[p < 0 | 1 < p] <- NaN
  ans[lambda < 0] <- NaN
  ans
}



roapospois <- function(n, lambda, pobs1 = 0) {
  qoapospois(runif(n), lambda = lambda, pobs1 = pobs1)
}






 oapospoisson <-
  function(lpobs1 = "logit",
           llambda = "loge",
           type.fitted = c("mean", "lambda", "pobs1", "onempobs1"),
           ipobs1 = NULL,
           zero = NULL) {


  lpobs1 <- as.list(substitute(lpobs1))
  epobs1 <- link2list(lpobs1)
  lpobs1 <- attr(epobs1, "function.name")

  llambd <- as.list(substitute(llambda))
  elambd <- link2list(llambd)
  llambd <- attr(elambd, "function.name")


 type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "pobs1", "onempobs1"))[1]


  new("vglmff",
  blurb = c("One-altered positive-Poisson distribution \n",
            "(Bernoulli and 1-truncated positive-Poisson ",
            "distribution model)\n\n",
            "Links:    ",
            namesof("pobs1",  lpobs1, earg = epobs1, tag = FALSE), ", ",
            namesof("lambda", llambd, earg = elambd, tag = FALSE)),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("pobs1", "lambda"),
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    M1 <- 2
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    extra$y1 <- y1 <- ifelse(y == 1, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y1), n, NOS)
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("pobs1",  ncoly)
    mynames2 <- param.names("lambda", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lpobs1 , earg = .epobs1 , tag = FALSE),
          namesof(mynames2, .llambd , earg = .elambd , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly


    if (!length(etastart)) {
      Init.lambda <- y - 0.25
      etastart <-
        cbind(theta2eta(if (length( .ipobs1 )) .ipobs1 else
                        (0.5 + w * y1) / (1 + w),
                        .lpobs1 , earg = .epobs1 ),
              theta2eta(Init.lambda, .llambd , earg = .elambd ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .llambd = llambd, .elambd = elambd,
            .lpobs1 = lpobs1, .epobs1 = epobs1,
            .ipobs1 = ipobs1,  # .ilambd = ilambd,
            .type.fitted = type.fitted
           ))),


  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "pobs1", "onempobs1"))[1]

    M1 <- 2
    NOS <- ncol(eta) / M1

    pobs1 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                             .lpobs1 , earg = .epobs1 ))
    lambd <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                             .llambd , earg = .elambd ))


    ans <- switch(type.fitted,
                  "mean"      = pobs1 + (1 - pobs1) *
                                lambd / ppois(1, lambd, lower = FALSE),
                  "lambda"      = lambd,
                  "pobs1"     =      pobs1,  # P(Y=1)
                  "onempobs1" =  1 - pobs1)  # P(Y>1)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpobs1 = lpobs1, .llambd = llambd,
           .epobs1 = epobs1, .elambd = elambd ))),
  last = eval(substitute(expression({
    temp.names <- c(rep_len( .lpobs1 , NOS),
                    rep_len( .llambd , NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epobs1
      misc$earg[[M1*ii  ]] <- .elambd
    }
  }), list( .lpobs1 = lpobs1, .llambd = llambd,
            .epobs1 = epobs1, .elambd = elambd ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pobs1  <- cbind(eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                              .lpobs1, earg = .epobs1))
    lambd <- cbind(eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                              .llambd, earg = .elambd ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doapospois(x = y, pobs1 = pobs1, lambda = lambd,
                                   log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs1 = lpobs1, .llambd = llambd,
           .epobs1 = epobs1, .elambd = elambd ))),
  vfamily = c("oapospoisson"),



  validparams = eval(substitute(function(eta, y, extra = NULL) {
    TFvec <- c(TRUE, FALSE)
    pobs1 <- eta2theta(eta[,  TFvec, drop = FALSE],
                       .lpobs1 , earg = .epobs1 )
    lambd <- eta2theta(eta[, !TFvec, drop = FALSE],
                       .llambd , earg = .elambd )
    okay1 <- all(is.finite(lambd)) && all(0 < lambd) &&
             all(is.finite(pobs1)) && all(0 < pobs1 & pobs1 < 1)
    okay1
  }, list( .lpobs1 = lpobs1, .llambd = llambd,
           .epobs1 = epobs1, .elambd = elambd ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpobs1 , earg = .epobs1 )
    lambd <- eta2theta(eta[, c(FALSE, TRUE)], .llambd , earg = .elambd )
    roapospois(nsim * length(lambd), lambd = lambd, pobs1 = pobs1)
  }, list( .lpobs1 = lpobs1, .llambd = llambd,
           .epobs1 = epobs1, .elambd = elambd ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1  # extra$NOS
    y1 <- extra$y1
    skip <- extra$skip.these

    TFvec <- c(TRUE, FALSE)
    pobs1  <- eta2theta(eta[,  TFvec, drop = FALSE],
                        .lpobs1 , earg = .epobs1 )
    lambda <- eta2theta(eta[, !TFvec, drop = FALSE],
                        .llambd , earg = .elambd )

    EY.cond <- 1 / ppois(1, lambda, lower.tail = FALSE)
    temp1 <- expm1(lambda)
    temp0 <- lambda * exp(-lambda)
    shape.geq.2 <- -expm1(-lambda) - temp0
    dl.dlambd <- y / lambda - 1 - temp0 / shape.geq.2


    dl.dpobs1 <- -1 / (1 - pobs1)  # For y > 1 obsns

    for (spp. in 1:NOS) {
      dl.dpobs1[skip[, spp.], spp.] <- 1 / pobs1[skip[, spp.], spp.]
      dl.dlambd[skip[, spp.], spp.] <- 0
    }
    dlambd.deta <- dtheta.deta(lambda, .llambd , earg = .elambd )
    mu.phi1 <- pobs1

    temp3 <- if ( .lpobs1 == "logit") {
      c(w) * (y1 - mu.phi1)
    } else {
      c(w) * dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 ) *
             dl.dpobs1
    }

    ans <- cbind(temp3,
                 c(w) * dl.dlambd * dlambd.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lpobs1 = lpobs1, .llambd = llambd,
            .epobs1 = epobs1, .elambd = elambd ))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, M1 * NOS)  # EIM is diagonal

    ned2l.dlambd2 <- EY.cond / lambda +
        ((1 - lambda) *
         exp(-lambda) - temp0^2 / shape.geq.2) / shape.geq.2

    ned2l.dlambd2 <- (1 - pobs1) * ned2l.dlambd2
    wz[, NOS+(1:NOS)] <- c(w) * ned2l.dlambd2 * dlambd.deta^2


    tmp100 <- mu.phi1 * (1 - mu.phi1)
    tmp200 <- if ( .lpobs1 == "logit" && is.empty.list( .epobs1 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 )^2)
    }
    wz[, 1:NOS] <-  tmp200

    wz <- wz[, interleave.VGAM(ncol(wz), M1 = M1)]

    wz
  }), list( .lpobs1 = lpobs1,
            .epobs1 = epobs1 ))))
}  # End of oapospoisson








doazeta <- function(x, shape, pobs1 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(pobs1))
  if (length(x)     != LLL) x     <- rep_len(x,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pobs1) != LLL) pobs1 <- rep_len(pobs1, LLL)
  ans <- rep_len(0.0, LLL)


  index1 <- (x == 1)

  if (log.arg) {
    ans[ index1] <- log(pobs1[index1])
    ans[!index1] <- log1p(-pobs1[!index1]) +
                    dotzeta(x[!index1], shape[!index1], log = TRUE)
  } else {
    ans[ index1] <- pobs1[index1]
    ans[!index1] <- (1 - pobs1[!index1]) *
                    dotzeta(x[!index1], shape[!index1])
  }
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[shape <= 0] <- NaN
  ans
}



poazeta <- function(q, shape, pobs1 = 0) {
  LLL <- max(length(q), length(shape), length(pobs1))
  if (length(q)     != LLL) q     <- rep_len(q,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pobs1) != LLL) pobs1 <- rep_len(pobs1, LLL)
  ans <- rep_len(0.0, LLL)

  ans[q >  1] <-    pobs1[q > 1] +
                 (1-pobs1[q > 1]) * potzeta(q[q > 1], shape[q > 1])
  ans[q <  1] <- 0
  ans[q == 1] <- pobs1[q == 1]

  ans <- pmax(0, ans)
  ans <- pmin(1, ans)
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[shape <= 0] <- NaN

  ans
}



qoazeta <- function(p, shape, pobs1 = 0) {
  LLL <- max(length(p), length(shape), length(pobs1))
  if (length(p)     != LLL) p     <- rep_len(p,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pobs1) != LLL) pobs1 <- rep_len(pobs1, LLL)

  ans <- rep_len(NaN, LLL)
  ind4 <- pobs1 < p
  ans[!ind4] <- 1
  ans[ ind4] <- qotzeta((p[ind4] - pobs1[ind4]) / (1 - pobs1[ind4]),
                        shape = shape[ind4])
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[p < 0 | 1 < p] <- NaN
  ans[shape <= 0] <- NaN
  ans
}



roazeta <- function(n, shape, pobs1 = 0) {
  qoazeta(runif(n), shape = shape, pobs1 = pobs1)
}






 oazeta <-
  function(lpobs1 = "logit",
           lshape = "loge",
           type.fitted = c("mean", "shape", "pobs1", "onempobs1"),
           gshape = exp((-4:3)/4),
           ishape = NULL,
           ipobs1 = NULL,
           zero = NULL) {


  lpobs1 <- as.list(substitute(lpobs1))
  epobs1 <- link2list(lpobs1)
  lpobs1 <- attr(epobs1, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


 type.fitted <- match.arg(type.fitted,
                          c("mean", "shape", "pobs1", "onempobs1"))[1]


  new("vglmff",
  blurb = c("One-altered zeta distribution \n",
            "(Bernoulli and 1-truncated zeta distribution model)\n\n",
            "Links:    ",
            namesof("pobs1", lpobs1, earg = epobs1, tag = FALSE), ", ",
            namesof("shape", lshape, earg = eshape, tag = FALSE)),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("pobs1", "shape"),
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    M1 <- 2
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    extra$y1 <- y1 <- ifelse(y == 1, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y1), n, NOS)
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("pobs1", ncoly)
    mynames2 <- param.names("shape", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lpobs1 , earg = .epobs1 , tag = FALSE),
          namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly


    if (!length(etastart)) {
      otzetaff.Loglikfun <- function(shape, y, x, w, extraargs) {
        sum(c(w) * dotzeta(x = y, shape, log = TRUE))
      }

      gshape <- .gshape
      if (!length( .ishape )) {
        shape.init <- matrix(NA_real_, n, M/M1, byrow = TRUE)
        for (jay in 1:ncoly) {
          index1 <- y[, jay] > 1
          shape.init[, jay] <-
            grid.search(gshape, objfun = otzetaff.Loglikfun,  # x = x,
                        y = y[index1, jay], w = w[index1, jay])
        }
      } else {
        shape.init <- matrix( .ishape , n, M, byrow = TRUE)
      }

      etastart <-
        cbind(theta2eta(if (length( .ipobs1 )) .ipobs1 else
                        (0.5 + w * y1) / (1 + w),
                        .lpobs1 , earg = .epobs1 ),
              theta2eta(shape.init, .lshape , earg = .eshape ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lshape = lshape, .eshape = eshape,
            .lpobs1 = lpobs1, .epobs1 = epobs1,
            .ipobs1 = ipobs1, .ishape = ishape,
                              .gshape = gshape,
            .type.fitted = type.fitted
           ))),


  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "shape", "pobs1", "onempobs1"))[1]

    M1 <- 2
    NOS <- ncol(eta) / M1

    pobs1 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                             .lpobs1 , earg = .epobs1 ))
    shape <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                             .lshape , earg = .eshape ))
    if (type.fitted == "mean") {
      ans <- shape
      ans[shape > 1] <- zeta(shape[shape > 1])/zeta(shape[shape > 1] + 1)
      ans[shape <= 1] <- NA
      pmf.1 <- dzeta(1, shape)
      mean.otzeta <- (ans - pmf.1) / (1 - pmf.1)
    }

    ans <- switch(type.fitted,
                  "mean"      = pobs1 + (1 - pobs1) * mean.otzeta,
                  "shape"     = shape,
                  "pobs1"     =      pobs1,  # P(Y=1)
                  "onempobs1" =  1 - pobs1)  # P(Y>1)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),
  last = eval(substitute(expression({
    temp.names <- c(rep_len( .lpobs1 , NOS),
                    rep_len( .lshape , NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epobs1
      misc$earg[[M1*ii  ]] <- .eshape
    }
  }), list( .lpobs1 = lpobs1, .lshape = lshape,
            .epobs1 = epobs1, .eshape = eshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pobs1 <- cbind(eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                             .lpobs1, earg = .epobs1 ))
    shape <- cbind(eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                             .lshape, earg = .eshape ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doazeta(x = y, pobs1 = pobs1, shape = shape,
                                log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),
  vfamily = c("oazeta"),



  validparams = eval(substitute(function(eta, y, extra = NULL) {
    TFvec <- c(TRUE, FALSE)
    pobs1 <- eta2theta(eta[,  TFvec, drop = FALSE],
                       .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, !TFvec, drop = FALSE],
                       .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape) &&
             all(is.finite(pobs1)) && all(0 < pobs1 & pobs1 < 1)
    okay1
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    roazeta(nsim * length(shape), shape = shape, pobs1 = pobs1)
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1  # extra$NOS
    y1 <- extra$y1
    skip <- extra$skip.these

    TFvec <- c(TRUE, FALSE)
    pobs1 <- eta2theta(eta[,  TFvec, drop = FALSE],
                       .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, !TFvec, drop = FALSE],
                       .lshape , earg = .eshape )

    BBBB  <- zeta(shape + 1) - 1
    fred1 <- zeta(shape + 1, deriv = 1)
    dl.dshape <- -log(y) - fred1 / BBBB

    dl.dpobs1 <- -1 / (1 - pobs1)  # For y > 1 obsns

    for (spp. in 1:NOS) {
      dl.dpobs1[skip[, spp.], spp.] <- 1 / pobs1[skip[, spp.], spp.]
      dl.dshape[skip[, spp.], spp.] <- 0
    }
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    mu.phi1 <- pobs1

    temp3 <- if ( .lpobs1 == "logit") {
      c(w) * (y1 - mu.phi1)
    } else {
      c(w) * dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 ) *
             dl.dpobs1
    }

    ans <- cbind(temp3,
                 c(w) * dl.dshape * dshape.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lpobs1 = lpobs1, .lshape = lshape,
            .epobs1 = epobs1, .eshape = eshape ))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, M1 * NOS)  # EIM is diagonal

    ned2l.dshape2 <- (zeta(shape + 1, deriv = 2) - fred1^2 / BBBB) / BBBB

    ned2l.dshape2 <- (1 - pobs1) * ned2l.dshape2
    wz[, NOS+(1:NOS)] <- c(w) * ned2l.dshape2 * dshape.deta^2


    tmp100 <- mu.phi1 * (1 - mu.phi1)
    tmp200 <- if ( .lpobs1 == "logit" && is.empty.list( .epobs1 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 )^2)
    }
    wz[, 1:NOS] <-  tmp200

    wz <- wz[, interleave.VGAM(ncol(wz), M1 = M1)]

    wz
  }), list( .lpobs1 = lpobs1,
            .epobs1 = epobs1 ))))
}  # End of oazeta













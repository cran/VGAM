# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.














hzeta.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 hzeta <- function(lshape = "loglog", ishape = NULL, nsimEIM = 100) {

  stopifnot(ishape > 0)
  stopifnot(nsimEIM > 10,
            length(nsimEIM) == 1,
            nsimEIM == round(nsimEIM))



  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  new("vglmff",
  blurb = c("Haight's Zeta distribution f(y) = (2y-1)^(-shape) - ",
            "(2y+1)^(-shape),\n",
            "    shape>0, y = 1, 2,....\n\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape), "\n\n",
            "Mean:     (1-2^(-shape)) * zeta(shape) if shape>1",
            "\n",
            "Variance: (1-2^(1-shape)) * zeta(shape-1) - mean^2 if shape>2"),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = FALSE,
         multipleResponses = FALSE,
         parameters.names = c("shape"),
         lshape = .lshape ,
         nsimEIM = .nsimEIM )
  }, list( .nsimEIM = nsimEIM, .lshape = lshape ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              Is.integer.y = TRUE,
              Is.positive.y = TRUE)


    predictors.names <-
      namesof("shape", .lshape , earg = .eshape , tag = FALSE)

    if (!length(etastart)) {
      a.init <- if (length( .ishape)) .ishape else {
        if ((meany <- weighted.mean(y, w)) < 1.5) 3.0 else
        if (meany < 2.5) 1.4 else 1.1
      }
      a.init <- rep_len(a.init, n)
      etastart <- theta2eta(a.init, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .eshape = eshape, .ishape = ishape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    mu <- (1-2^(-shape)) * zeta(shape)
    mu[shape <= 1] <- Inf
    mu
  }, list( .lshape = lshape, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape = .lshape)

    misc$earg <- list(shape = .eshape )

    misc$nsimEIM <- .nsimEIM

  }), list( .lshape = lshape, .eshape = eshape, .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dhzeta(x = y, shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("hzeta"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape)
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
    rhzeta(nsim * length(shape), shape = shape)
  }, list( .lshape = lshape, .eshape = eshape ))),



  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .eshape )

    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    d3 <- deriv3(~ log((2*y-1)^(-shape) - (2*y+1)^(-shape)),
                 "shape", hessian = FALSE)
    eval.d3 <- eval(d3)

    dl.dshape <-  attr(eval.d3, "gradient")

    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = eval(substitute(expression({
    sd3 <- deriv3(~ log((2*ysim-1)^(-shape) - (2*ysim+1)^(-shape)),
                  "shape", hessian = FALSE)
    run.var <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- rhzeta(n, shape = shape)
      eval.sd3 <- eval(sd3)
      dl.dshape <-  attr(eval.d3, "gradient")
      rm(ysim)
      temp3 <- dl.dshape
      run.var <- ((ii-1) * run.var + temp3^2) / ii
    }
    wz <- if (intercept.only)
        matrix(colMeans(cbind(run.var)),
               n, dimm(M), byrow = TRUE) else cbind(run.var)

    wz <- wz * dshape.deta^2
    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))
}




dhzeta <- function(x, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(shape, positive = TRUE))
    stop("'shape' must be numeric and have positive values")

  nn <- max(length(x), length(shape))
  if (length(x)     != nn) x     <- rep_len(x,     nn)
  if (length(shape) != nn) shape <- rep_len(shape, nn)

  ox <- !is.finite(x)
  zero <- ox | round(x) != x | x < 1
  ans <- rep_len(0, nn)
  ans[!zero] <- (2*x[!zero]-1)^(-shape[!zero]) -
                (2*x[!zero]+1)^(-shape[!zero])
  if (log.arg) log(ans) else ans
}



phzeta <- function(q, shape, log.p = FALSE) {


  nn <- max(length(q), length(shape))
  q     <- rep_len(q,     nn)
  shape <- rep_len(shape, nn)
  oq <- !is.finite(q)
  zero <- oq | q < 1
  q <- floor(q)
  ans <- 0 * q
  ans[!zero] <- 1 - (2*q[!zero]+1)^(-shape[!zero])

  ans[q == -Inf] <- 0  # 20141215 KaiH
  ans[q ==  Inf] <- 1  # 20141215 KaiH

  ans[shape <= 0] <- NaN
  if (log.p) log(ans) else ans
}



qhzeta <- function(p, shape) {

  if (!is.Numeric(p, positive = TRUE) ||
      any(p >= 1))
    stop("argument 'p' must have values inside the interval (0,1)")

  nn <- max(length(p), length(shape))
  p     <- rep_len(p,     nn)
  shape <- rep_len(shape, nn)
  ans <- (((1 - p)^(-1/shape) - 1) / 2)  # p is in (0,1)
  ans[shape <= 0] <- NaN
  floor(ans + 1)
}


rhzeta <- function(n, shape) {


  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  shape <- rep_len(shape, use.n)
  ans <- (runif(use.n)^(-1/shape) - 1) / 2
  ans[shape <= 0] <- NaN
  floor(ans + 1)
}












dkumar <- function(x, shape1, shape2, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  N <- max(length(x), length(shape1), length(shape2))
  if (length(x)      != N) x      <- rep_len(x,      N)
  if (length(shape1) != N) shape1 <- rep_len(shape1, N)
  if (length(shape2) != N) shape2 <- rep_len(shape2, N)

  logdensity <- rep_len(log(0), N)
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



qkumar <- function(p, shape1, shape2,
                   lower.tail = TRUE, log.p = FALSE) {



  if (!is.logical(lower.tail) || length(lower.tail) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- (-expm1((1/shape2) * log(-expm1(ln.p))))^(1/shape1)
      ans[ln.p > 0] <- NaN
    } else {
      ans <- (-expm1((1/shape2) * log1p(-p)))^(1/shape1)
      ans[p < 0] <- NaN
      ans[p == 0] <- 0
      ans[p == 1] <- 1
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- (-expm1(ln.p / shape2))^(1/shape1)
      ans[ln.p > 0] <- NaN
      ans
    } else {
      ans <- (-expm1((1/shape2) * log(p)))^(1/shape1)
      ans[p < 0] <- NaN
      ans[p == 0] <- 1
      ans[p == 1] <- 0
      ans[p > 1] <- NaN
    }
  }
  ans[(shape1 <= 0) | (shape2 <= 0)] = NaN
  ans
}



pkumar <- function(q, shape1, shape2,
                   lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ans <- log(-expm1(shape2 * log1p(-q^shape1)))
      ans[q <= 0 ] <- -Inf
      ans[q >= 1] <- 0
    } else {
      ans <- -expm1(shape2 * log1p(-q^shape1))
      ans[q <= 0] <- 0
      ans[q >= 1] <- 1
    }
  } else {
    if (log.p) {
      ans <- shape2 * log1p(-q^shape1)
      ans[q <= 0] <- 0
      ans[q >= 1] <- -Inf
    } else {
      ans <- exp(shape2 * log1p(-q^shape1))
      ans[q <= 0] <- 1
      ans[q >= 1] <- 0
    }
  }

  ans[(shape1 <= 0) | (shape2 <= 0)] <- NaN
  ans
}










 kumar <-
  function(lshape1 = "loge", lshape2 = "loge",
           ishape1 = NULL,   ishape2 = NULL,
           gshape1 = exp(2*ppoints(5) - 1), tol12 = 1.0e-4, zero = NULL) {
  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")
  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")

  if (length(ishape1) &&
     (!is.Numeric(ishape1, length.arg = 1, positive = TRUE)))
    stop("bad input for argument 'ishape1'")
  if (length(ishape2) && !is.Numeric(ishape2))
    stop("bad input for argument 'ishape2'")

  if (!is.Numeric(tol12, length.arg = 1, positive = TRUE))
    stop("bad input for argument 'tol12'")
  if (!is.Numeric(gshape1, positive = TRUE))
    stop("bad input for argument 'gshape1'")


  new("vglmff",
  blurb = c("Kumaraswamy distribution\n\n",
            "Links:    ",
              namesof("shape1", lshape1, eshape1, tag = FALSE), ", ",
              namesof("shape2", lshape2, eshape2, tag = FALSE), "\n",
            "Mean:     shape2 * beta(1 + 1 / shape1, shape2)"),
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
         parameters.names = c("shape1", "shape2"),
         lshape1 = .lshape1 ,
         lshape2 = .lshape2 ,
         zero = .zero )
  }, list( .zero = zero, .lshape1 = lshape1, .lshape2 = lshape2 ))),

  initialize = eval(substitute(expression({
    checklist <- w.y.check(w = w, y = y, Is.positive.y = TRUE,
                           ncol.w.max = Inf, ncol.y.max = Inf,
                           out.wy = TRUE, colsyperw = 1, maximize = TRUE)
    w <- checklist$w
    y <- checklist$y  # Now 'w' and 'y' have the same dimension.
    if (any((y <= 0) | (y >= 1)))
      stop("the response must be in (0, 1)")

    extra$ncoly <- ncoly <- ncol(y)
    extra$M1 <- M1 <- 2
    M <- M1 * ncoly
    mynames1 <- param.names("shape1", ncoly)
    mynames2 <- param.names("shape2", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lshape1 , earg = .eshape1 , tag = FALSE),
          namesof(mynames2, .lshape2 , earg = .eshape2 , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]

    if (!length(etastart)) {
      kumar.Loglikfun <- function(shape1, y, x, w, extraargs) {
        mediany <- colSums(y * w) / colSums(w)
        shape2 <- log(0.5) / log1p(-(mediany^shape1))
        sum(c(w) * dkumar(y, shape1 = shape1, shape2 = shape2, log = TRUE))
      }

      shape1.grid <- .gshape1
      shape1.init <- if (length( .ishape1 )) .ishape1 else
        grid.search(shape1.grid, objfun = kumar.Loglikfun,
                    y = y,  x = x, w = w)
      shape1.init <- matrix(shape1.init, n, ncoly, byrow = TRUE)

      mediany <- colSums(y * w) / colSums(w)
      shape2.init <- if (length( .ishape2 )) .ishape2 else
        log(0.5) / log1p(-(mediany^shape1.init))
      shape2.init <- matrix(shape2.init, n, ncoly, byrow = TRUE)

      etastart <- cbind(theta2eta(shape1.init, .lshape1 , earg = .eshape1 ),
                        theta2eta(shape2.init, .lshape2 , earg = .eshape2 ))[,
                  interleave.VGAM(M, M1 = M1)]
    }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .ishape1 = ishape1, .ishape2 = ishape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .gshape1 = gshape1 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2 , earg = .eshape2 )
    shape2 * (base::beta(1 + 1/shape1, shape2))
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  last = eval(substitute(expression({
    misc$link <- c(rep_len( .lshape1 , ncoly),
                   rep_len( .lshape2 , ncoly))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .eshape1
      misc$earg[[M1*ii  ]] <- .eshape2
    }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  loglikelihood = eval(substitute(
  function(mu, y, w, residuals = FALSE, eta, extra = NULL, summation = TRUE) {
    shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2 , earg = .eshape2 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dkumar(x = y, shape1, shape2, log = TRUE)
      if (summation) sum(ll.elts) else ll.elts
    }
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  vfamily = c("kumar"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2 , earg = .eshape2 )
    okay1 <- all(is.finite(shape1)) && all(0 < shape1) &&
             all(is.finite(shape2)) && all(0 < shape2)
    okay1
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  simslot = eval(substitute(
  function(object, nsim) {
    eta <- predict(object)
    shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2 , earg = .eshape2 )
    rkumar(nsim * length(shape1), shape1 = shape1, shape2 = shape2)
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  deriv = eval(substitute(expression({
    shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2 , earg = .eshape2 )
    dshape1.deta <- dtheta.deta(shape1, link = .lshape1 , earg = .eshape1 )
    dshape2.deta <- dtheta.deta(shape2, link = .lshape2 , earg = .eshape2 )
    dl.dshape1 <- 1 / shape1 + log(y) - (shape2 - 1) * log(y) *
                  (y^shape1) / (1 - y^shape1)
    dl.dshape2 <- 1 / shape2 + log1p(-y^shape1)
    dl.deta <- c(w) * cbind(dl.dshape1 * dshape1.deta,
                            dl.dshape2 * dshape2.deta)
    dl.deta[, interleave.VGAM(M, M1 = M1)]
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  weight = eval(substitute(expression({
    ned2l.dshape11 <- (1 + (shape2 / (shape2 - 2)) *
      ((digamma(shape2) -  digamma(2))^2 -
      (trigamma(shape2) - trigamma(2)))) / shape1^2
    ned2l.dshape22 <- 1 / shape2^2
    ned2l.dshape12 <-
       (digamma(2) - digamma(1 + shape2)) / ((shape2 - 1) * shape1)

    index1 <- (abs(shape2 - 1) < .tol12 )  # Fix up singular point at shape2 == 1
    ned2l.dshape12[index1] <- -trigamma(2) / shape1[index1]
    index2 <- (abs(shape2 - 2) < .tol12 )  # Fix up singular point at shape2 == 2
    ned2l.dshape11[index2] <- (1 - 2 * psigamma(2, deriv = 2)) / shape1[index2]^2

    wz <- array(c(c(w) * ned2l.dshape11 * dshape1.deta^2,
                  c(w) * ned2l.dshape22 * dshape2.deta^2,
                  c(w) * ned2l.dshape12 * dshape1.deta * dshape2.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2, .tol12 = tol12 ))))
}






drice <- function(x, sigma, vee, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)




  N <- max(length(x), length(vee), length(sigma))
  if (length(x)      != N) x      <- rep_len(x,      N)
  if (length(vee)    != N) vee    <- rep_len(vee   , N)
  if (length(sigma ) != N) sigma  <- rep_len(sigma , N)

  logdensity <- rep_len(log(0), N)
  xok <- (x > 0)
  x.abs <- abs(x[xok] * vee[xok] / sigma[xok]^2)
  logdensity[xok] <- log(x[xok]) - 2 * log(sigma[xok]) +
                     (-(x[xok]^2+vee[xok]^2)/(2*sigma[xok]^2)) +
                     log(besselI(x.abs, nu = 0, expon.scaled = TRUE)) +
                     x.abs
  logdensity[sigma <= 0] <- NaN
  logdensity[vee < 0] <- NaN

  logdensity[is.infinite(x)] <- -Inf  # 20141209 KaiH

  if (log.arg) logdensity else exp(logdensity)
}



rrice <- function(n, sigma, vee) {
  theta <- 1  # any number
  X <- rnorm(n, mean = vee * cos(theta), sd = sigma)
  Y <- rnorm(n, mean = vee * sin(theta), sd = sigma)
  sqrt(X^2 + Y^2)
}




marcumQ <- function(a, b, m = 1,
                    lower.tail = TRUE, log.p = FALSE, ... ) {
  pchisq(b^2, df = 2*m, ncp = a^2,
         lower.tail = lower.tail, log.p = log.p, ... )
}



price <- function(q, sigma, vee,
                  lower.tail = TRUE, log.p = FALSE, ...) {
  ans <- marcumQ(vee/sigma, q/sigma, m = 1,
                 lower.tail = lower.tail, log.p = log.p, ... )
  ans
}



qrice <- function(p, sigma, vee,
                  lower.tail = TRUE, log.p = FALSE, ... ) {
  sqrt(qchisq(p, df = 2, ncp = (vee/sigma)^2,
              lower.tail = lower.tail, log.p = log.p, ... )) * sigma
}









riceff.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}



 riceff <- function(lsigma = "loge", lvee = "loge",
                    isigma = NULL, ivee = NULL,
                    nsimEIM = 100, zero = NULL, nowarning = FALSE) {



  lvee     <- as.list(substitute(lvee))
  evee     <- link2list(lvee)
  lvee     <- attr(evee, "function.name")


  lsigma <- as.list(substitute(lsigma))
  esigma <- link2list(lsigma)
  lsigma <- attr(esigma, "function.name")



  if (length(ivee) && !is.Numeric(ivee, positive = TRUE))
    stop("bad input for argument 'ivee'")
  if (length(isigma) && !is.Numeric(isigma, positive = TRUE))
    stop("bad input for argument 'isigma'")

  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
    stop("'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Rice distribution\n\n",
            "Links:    ",
            namesof("sigma", lsigma, earg = esigma, tag = FALSE), ", ",
            namesof("vee",   lvee,   earg = evee,   tag = FALSE), "\n",
            "Mean:     ",
            "sigma*sqrt(pi/2)*exp(z/2)*((1-z)*",
            "besselI(-z/2, nu = 0) - z * besselI(-z/2, nu = 1)) ",
            "where z=-vee^2/(2*sigma^2)"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = FALSE,
         multipleResponses = FALSE,
         parameters.names = c("sigma", "vee"),
         nsimEIM = .nsimEIM,
         lsigma = .lsigma ,
         lvee = .lvee ,
         zero = .zero )
  }, list( .zero = zero, .lsigma = lsigma, .lvee = lvee,
           .nsimEIM = nsimEIM ))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      c(namesof("sigma", .lsigma , earg = .esigma , tag = FALSE),
        namesof("vee",   .lvee   , earg = .evee   , tag = FALSE))




    if (!length(etastart)) {
      riceff.Loglikfun <- function(vee, y, x, w, extraargs) {
        sigma.init <- sd(rep(y, w))
        sum(c(w) * (log(y) - 2*log(sigma.init) +
                    log(besselI(y*vee/sigma.init^2, nu = 0)) -
                   (y^2 + vee^2) / (2*sigma.init^2)))
      }
    vee.grid <-
      seq(quantile(rep(y, w), probs = seq(0, 1, 0.2))["20%"],
          quantile(rep(y, w), probs = seq(0, 1, 0.2))["80%"], len = 11)
    vee.init <- if (length( .ivee )) .ivee else
      grid.search(vee.grid, objfun = riceff.Loglikfun, y = y,  x = x, w = w)
      vee.init <- rep_len(vee.init, length(y))
      sigma.init <- if (length( .isigma )) .isigma else
          sqrt(max((weighted.mean(y^2, w) - vee.init^2)/2, 0.001))
      sigma.init <- rep_len(sigma.init, length(y))

      etastart <-
        cbind(theta2eta(sigma.init, .lsigma , earg = .esigma ),
              theta2eta(vee.init,   .lvee ,   earg = .evee ))
    }
  }), list( .lvee = lvee, .lsigma = lsigma,
            .ivee = ivee, .isigma = isigma,
            .evee = evee, .esigma = esigma ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    vee   <- eta2theta(eta[, 1], link = .lvee ,   earg = .evee )
    sigma <- eta2theta(eta[, 2], link = .lsigma , earg = .esigma )
    temp9 <- -vee^2 / (2*sigma^2)


      sigma * sqrt(pi/2) *
      ((1-temp9) * besselI(-temp9/2, nu = 0, expon = TRUE) -
          temp9  * besselI(-temp9/2, nu = 1, expon = TRUE))
  }, list( .lvee = lvee, .lsigma = lsigma,
           .evee = evee, .esigma = esigma ))),
  last = eval(substitute(expression({
    misc$link <-    c("sigma" = .lsigma , "vee" = .lvee )

    misc$earg <- list("sigma" = .esigma , "vee" = .evee )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lvee = lvee, .lsigma = lsigma,
            .evee = evee, .esigma = esigma, .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    sigma <- eta2theta(eta[, 1], link = .lsigma , earg = .esigma )
    vee   <- eta2theta(eta[, 2], link = .lvee   , earg = .evee )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * drice(x = y, sigma = sigma, vee = vee, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lvee = lvee, .lsigma = lsigma,
           .evee = evee, .esigma = esigma ))),
  vfamily = c("riceff"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    sigma <- eta2theta(eta[, 1], link = .lsigma , earg = .esigma )
    vee   <- eta2theta(eta[, 2], link = .lvee   , earg = .evee )
    okay1 <- all(is.finite(sigma)) && all(0 < sigma) &&
             all(is.finite(vee  )) && all(0 < vee  )
    okay1
  }, list( .lvee = lvee, .lsigma = lsigma,
           .evee = evee, .esigma = esigma ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    sigma <- eta2theta(eta[, 1], link = .lsigma , earg = .esigma )
    vee   <- eta2theta(eta[, 2], link = .lvee   , earg = .evee )
    rrice(nsim * length(vee),
          vee = vee, sigma = sigma)
  }, list( .lvee = lvee, .lsigma = lsigma,
           .evee = evee, .esigma = esigma ))),



  deriv = eval(substitute(expression({
    sigma <- eta2theta(eta[, 1], link = .lsigma , earg = .esigma )
    vee   <- eta2theta(eta[, 2], link = .lvee   , earg = .evee )

    dvee.deta <- dtheta.deta(vee, link = .lvee , earg = .evee )
    dsigma.deta <- dtheta.deta(sigma, link = .lsigma , earg = .esigma )

    temp8 <- y * vee / sigma^2
    dl.dvee <- -vee/sigma^2 + (y/sigma^2) *
               besselI(temp8, nu = 1) / besselI(temp8, nu = 0)
    dl.dsigma <- -2/sigma + (y^2 + vee^2)/(sigma^3) -
                 (2 * temp8 / sigma) *
                 besselI(temp8, nu = 1) / besselI(temp8, nu = 0)

    c(w) * cbind(dl.dsigma * dsigma.deta,
                 dl.dvee   * dvee.deta)

  }), list( .lvee = lvee, .lsigma = lsigma,
            .evee = evee, .esigma = esigma, .nsimEIM = nsimEIM ))),
  weight = eval(substitute(expression({
    run.var <- run.cov <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- rrice(n, vee = vee, sigma = sigma)
      temp8 <- ysim * vee / sigma^2
      dl.dvee <- -vee/sigma^2 + (ysim/sigma^2) *
                 besselI(temp8, nu = 1) / besselI(temp8, nu = 0)
      dl.dsigma <- -2/sigma + (ysim^2 + vee^2)/(sigma^3) -
                   (2 * temp8 / sigma) *
                   besselI(temp8, nu = 1) / besselI(temp8, nu = 0)

      rm(ysim)
      temp3 <- cbind(dl.dsigma, dl.dvee)
      run.var <- ((ii-1) * run.var + temp3^2) / ii
      run.cov <- ((ii-1) * run.cov + temp3[, 1] * temp3[, 2]) / ii
    }
    wz <- if (intercept.only)
        matrix(colMeans(cbind(run.var, run.cov)),
               n, dimm(M), byrow = TRUE) else cbind(run.var, run.cov)

    dtheta.detas <- cbind(dsigma.deta, dvee.deta)
    index0 <- iam(NA_real_, NA_real_, M = M, both = TRUE, diag = TRUE)
    wz <- wz * dtheta.detas[, index0$row] *
               dtheta.detas[, index0$col]
    c(w) * wz
  }), list( .lvee = lvee, .lsigma = lsigma,
            .evee = evee, .esigma = esigma, .nsimEIM = nsimEIM ))))
}





dskellam <- function(x, mu1, mu2, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  L <- max(length(x), length(mu1), length(mu2))
  if (length(x)      != L) x      <- rep_len(x,      L)
  if (length(mu1)    != L) mu1    <- rep_len(mu1,    L)
  if (length(mu2)    != L) mu2    <- rep_len(mu2,    L)

  ok2 <- is.finite(mu1) & is.finite(mu2) & (mu1 >= 0) & (mu2 >= 0)
  ok3 <- (mu1 == 0) & (mu2 >  0)
  ok4 <- (mu1 >  0) & (mu2 == 0)
  ok5 <- (mu1 == 0) & (mu2 == 0)
    if (log.arg) {
      ans <- -mu1 - mu2 + 2 * sqrt(mu1*mu2) +
             0.5 * x * log(mu1) - 0.5 * x * log(mu2) +
             log(besselI(2 * sqrt(mu1*mu2),

                         nu = abs(x),

                         expon.scaled = TRUE))
      ans[ok3] <- dpois(x = -x[ok3], lambda = mu2[ok3], log = TRUE)
      ans[ok4] <- dpois(x = -x[ok4], lambda = mu1[ok4], log = TRUE)
      ans[ok5] <- dpois(x =  x[ok5], lambda = 0.0,      log = TRUE)
      ans[x != round(x)] = log(0.0)
    } else {
      ans <- (mu1/mu2)^(x/2) * exp(-mu1-mu2 + 2 * sqrt(mu1*mu2)) *
             besselI(2 * sqrt(mu1*mu2),

                     nu = abs(x),

                     expon.scaled = TRUE)
      ans[ok3] <- dpois(x = -x[ok3], lambda = mu2[ok3])
      ans[ok4] <- dpois(x = -x[ok4], lambda = mu1[ok4])
      ans[ok5] <- dpois(x =  x[ok5], lambda = 0.0)
      ans[x != round(x)] <- 0.0
    }
    ans[!ok2] <- NaN
    ans
}






rskellam <- function(n, mu1, mu2) {
  rpois(n, mu1) - rpois(n, mu2)
}



skellam.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 skellam <- function(lmu1 = "loge", lmu2 = "loge",
                     imu1 = NULL,   imu2 = NULL,
                     nsimEIM = 100, parallel = FALSE, zero = NULL) {

  lmu1 <- as.list(substitute(lmu1))
  emu1 <- link2list(lmu1)
  lmu1 <- attr(emu1, "function.name")

  lmu2 <- as.list(substitute(lmu2))
  emu2 <- link2list(lmu2)
  lmu2 <- attr(emu2, "function.name")


  if (length(imu1) &&
      !is.Numeric(imu1, positive = TRUE))
    stop("bad input for argument 'imu1'")
  if (length(imu2) &&
      !is.Numeric(imu2, positive = TRUE))
    stop("bad input for argument 'imu2'")



  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
    stop("argument 'nsimEIM' should be an integer greater than 50")

  new("vglmff",
  blurb = c("Skellam distribution\n\n",
         "Links:    ",
         namesof("mu1", lmu1, earg = emu1, tag = FALSE), ", ",
         namesof("mu2", lmu2, earg = emu2, tag = FALSE), "\n",
         "Mean:     mu1-mu2", "\n",
         "Variance: mu1+mu2"),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints = constraints,
                           apply.int = TRUE)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .parallel = parallel, .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = FALSE,
         multipleResponses = FALSE,
         parameters.names = c("mu1", "mu2"),
         nsimEIM = .nsimEIM,
         lmu1 = .lmu1 ,
         lmu2 = .lmu2 ,
         zero = .zero )
  }, list( .zero = zero, .lmu1 = lmu1, .lmu2 = lmu2,
           .nsimEIM = nsimEIM ))),
  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <- c(
      namesof("mu1", .lmu1, earg = .emu1, tag = FALSE),
      namesof("mu2", .lmu2, earg = .emu2, tag = FALSE))


    if (!length(etastart)) {
      junk <- lm.wfit(x = x, y = c(y), w = c(w))
      var.y.est <- sum(c(w) * junk$resid^2) / junk$df.residual
      mean.init <- weighted.mean(y, w)

      mu1.init <- max((var.y.est + mean.init) / 2, 0.01)
      mu2.init <- max((var.y.est - mean.init) / 2, 0.01)
      mu1.init <- rep_len(if (length( .imu1 )) .imu1 else mu1.init, n)
      mu2.init <- rep_len(if (length( .imu2 )) .imu2 else mu2.init, n)

      etastart <- cbind(theta2eta(mu1.init, .lmu1, earg = .emu1 ),
                        theta2eta(mu2.init, .lmu2, earg = .emu2 ))
      }
  }), list( .lmu1 = lmu1, .lmu2 = lmu2,
            .imu1 = imu1, .imu2 = imu2,
            .emu1 = emu1, .emu2 = emu2 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
      mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1 )
      mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2 )
      mu1 - mu2
  }, list( .lmu1 = lmu1, .lmu2 = lmu2,
           .emu1 = emu1, .emu2 = emu2 ))),
  last = eval(substitute(expression({
      misc$link <-    c("mu1" = .lmu1, "mu2" = .lmu2)

      misc$earg <- list("mu1" = .emu1, "mu2" = .emu2 )

      misc$expected <- TRUE
      misc$nsimEIM <- .nsimEIM
  }), list( .lmu1 = lmu1, .lmu2 = lmu2,
            .emu1 = emu1, .emu2 = emu2,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1 )
    mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {

      ll.elts <-
        if ( is.logical( .parallel ) &&
             length( .parallel ) == 1 &&
             .parallel )
          c(w) * log(besselI(2*mu1, nu = y, expon = TRUE)) else
          c(w) * (-mu1 - mu2 +
                  0.5 * y * log(mu1) -
                  0.5 * y * log(mu2) +
                  2 * sqrt(mu1*mu2) +  # Use this when expon = TRUE
                  log(besselI(2 * sqrt(mu1*mu2), nu = y, expon = TRUE)))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmu1 = lmu1, .lmu2 = lmu2,
           .emu1 = emu1, .emu2 = emu2,
           .parallel = parallel ))),
  vfamily = c("skellam"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1 )
    mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2 )
    okay1 <- all(is.finite(mu1)) && all(0 < mu1) &&
             all(is.finite(mu2)) && all(0 < mu2)
    okay1
  }, list( .lmu1 = lmu1, .lmu2 = lmu2,
           .emu1 = emu1, .emu2 = emu2 ))),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1 )
    mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2 )
    rskellam(nsim * length(mu1), mu1, mu2)
  }, list( .lmu1 = lmu1, .lmu2 = lmu2,
           .emu1 = emu1, .emu2 = emu2,
           .parallel = parallel ))),





  deriv = eval(substitute(expression({
    mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1 )
    mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2 )

    dmu1.deta <- dtheta.deta(mu1, link = .lmu1, earg = .emu1 )
    dmu2.deta <- dtheta.deta(mu2, link = .lmu2, earg = .emu2 )

    temp8 <- 2 * sqrt(mu1*mu2)
    temp9 <-  besselI(temp8, nu = y  , expon = TRUE)
    temp7 <- (besselI(temp8, nu = y-1, expon = TRUE) +
              besselI(temp8, nu = y+1, expon = TRUE)) / 2
    temp6 <- temp7 / temp9

    dl.dmu1 <- -1 + 0.5 * y / mu1 + sqrt(mu2/mu1) * temp6
    dl.dmu2 <- -1 - 0.5 * y / mu2 + sqrt(mu1/mu2) * temp6

    c(w) * cbind(dl.dmu1 * dmu1.deta,
                 dl.dmu2 * dmu2.deta)
  }), list( .lmu1 = lmu1, .lmu2 = lmu2,
            .emu1 = emu1, .emu2 = emu2,
            .nsimEIM = nsimEIM ))),
  weight = eval(substitute(expression({
    run.var <- run.cov <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- rskellam(n, mu1=mu1, mu2=mu2)
      temp9 <-  besselI(temp8, nu = ysim,   expon = TRUE)
      temp7 <- (besselI(temp8, nu = ysim-1, expon = TRUE) +
                besselI(temp8, nu = ysim+1, expon = TRUE)) / 2
      temp6 <- temp7 / temp9
      dl.dmu1 <- -1 + 0.5 * ysim/mu1 + sqrt(mu2/mu1) * temp6
      dl.dmu2 <- -1 - 0.5 * ysim/mu2 + sqrt(mu1/mu2) * temp6
      rm(ysim)
      temp3 <- cbind(dl.dmu1, dl.dmu2)
      run.var <- ((ii-1) * run.var + temp3^2) / ii
      run.cov <- ((ii-1) * run.cov + temp3[, 1] * temp3[, 2]) / ii
    }
    wz <- if (intercept.only)
          matrix(colMeans(cbind(run.var, run.cov)),
                 n, dimm(M), byrow = TRUE) else
          cbind(run.var, run.cov)

    dtheta.detas <- cbind(dmu1.deta, dmu2.deta)
    index0 <- iam(NA_real_, NA_real_, M = M, both = TRUE, diag = TRUE)
    wz <- wz * dtheta.detas[, index0$row] *
               dtheta.detas[, index0$col]
    c(w) * wz
  }), list( .lmu1 = lmu1, .lmu2 = lmu2,
            .emu1 = emu1, .emu2 = emu2,
            .nsimEIM = nsimEIM ))))
}





dyules <- function(x, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if ( log.arg ) {
    ans <- log(shape) + lbeta(abs(x), shape+1)
    ans[(x != round(x)) | (x < 1)] <- log(0)
  } else {
    ans <- shape * beta(x, shape+1)
    ans[(x != round(x)) | (x < 1)] <- 0
  }
  ans[!is.finite(shape) | (shape <= 0)] <- NaN
  ans
}




pyules <- function(q, shape, lower.tail = TRUE, log.p = FALSE) {


  tq <- trunc(q)

  if (lower.tail) {
    ans <- 1 - tq * beta(abs(tq), shape+1)
    ans[q < 1] <- 0
    ans[is.infinite(q) & 0 < q] <- 1  # 20141215 KaiH
  } else {
    ans <-     tq * beta(abs(tq), shape+1)
    ans[q < 1] <- 1
    ans[is.infinite(q) & 0 < q] <- 0  # 20160713
  }

  ans[shape <= 0] <- NaN
  if (log.p) log(ans) else ans
  ans
}



 qyules <- function(p, shape) {

  LLL <- max(length(p), length(shape))
  if (length(p)     != LLL) p     <- rep_len(p,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  ans <- rep_len(0, LLL)

  lo <- rep_len(1, LLL)
  approx.ans <- lo  # True at lhs
  hi <- 2 * lo + 10
  dont.iterate <- p == 1 | shape <= 0
  done <- p <= pyules(hi, shape) | dont.iterate
  while (!all(done)) {
    hi.save <- hi[!done]
    hi[!done] <- 2 * lo[!done] + 10
    lo[!done] <- hi.save
    done[!done] <- (p[!done] <= pyules(hi[!done], shape[!done]))
  }

  foo <- function(q, shape, p)
    pyules(q, shape) - p

  lhs <- (p <= dyules(1, shape)) | dont.iterate

  approx.ans[!lhs] <- bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                                      shape = shape[!lhs], p = p[!lhs])
  faa <- floor(approx.ans)
  ans <- ifelse(pyules(faa, shape) < p & p <= pyules(faa+1, shape), faa+1, faa)

  ans[p == 1] <- Inf
  ans[shape <= 0] <- NaN

  ans
}  # qyules



ryules <- function(n, shape) {

  rgeom(n, prob = exp(-rexp(n, rate = shape))) + 1
}





yulesimon.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 yulesimon <- function(lshape = "loge",
                       ishape = NULL, nsimEIM = 200,
                       zero = NULL) {

  if (length(ishape) &&
      !is.Numeric(ishape, positive = TRUE))
    stop("argument 'ishape' must be > 0")



  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
    stop("argument 'nsimEIM' should be an integer greater than 50")



  new("vglmff",
  blurb = c("Yule-Simon distribution f(y) = shape * beta(y, shape + 1), ",
            "shape > 0, y = 1, 2,..\n\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape), "\n\n",
            "Mean:     shape / (shape - 1), provided shape > 1\n",
            "Variance: shape^2 / ((shape - 1)^2 * (shape - 2)), ",
            "provided shape > 2"),
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
         nsimEIM = .nsimEIM ,
         parameters.names = c("shape"),
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM ))),

  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    ncoly <- ncol(y)

    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1  <- param.names("shape", ncoly)
    predictors.names <-
      namesof(mynames1, .lshape , earg = .eshape , tag = FALSE)

    if (!length(etastart)) {
      wmeany <- colSums(y * w) / colSums(w) + 1/8

      shape.init <- wmeany / (wmeany - 1)
      shape.init <- matrix(if (length( .ishape )) .ishape else
                           shape.init, n, M, byrow = TRUE)
      etastart <- theta2eta(shape.init, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .eshape = eshape, .ishape = ishape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    ans <- shape <- eta2theta(eta, .lshape , earg = .eshape )
    ans[shape >  1] <- shape / (shape - 1)
    ans[shape <= 1] <- NA
    ans
  }, list( .lshape = lshape, .eshape = eshape ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <- c(rep_len( .lshape , ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .eshape
    }

    misc$M1 <- M1
    misc$ishape <- .ishape
    misc$nsimEIM <- .nsimEIM
  }), list( .lshape = lshape, .eshape = eshape, .nsimEIM = nsimEIM,
            .ishape = ishape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dyules(x = y, shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("yulesimon"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape)
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
    ryules(nsim * length(shape), shape = shape)
  }, list( .lshape = lshape, .eshape = eshape ))),







  deriv = eval(substitute(expression({
    M1 <- 1
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    dl.dshape <- 1/shape + digamma(1+shape) - digamma(1+shape+y)
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = eval(substitute(expression({

    run.var <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- ryules(n, shape <- shape)
      dl.dshape <- 1/shape + digamma(1+shape) - digamma(1+shape+ysim)
      rm(ysim)
      temp3 <- dl.dshape
      run.var <- ((ii-1) * run.var + temp3^2) / ii
    }
    wz <- if (intercept.only)
        matrix(colMeans(cbind(run.var)),
               n, M, byrow = TRUE) else cbind(run.var)

    wz <- wz * dshape.deta^2


    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))
}  # yule.simon()







dlind <- function(x, theta, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if ( log.arg ) {
    ans <- 2 * log(theta) + log1p(x) - theta * x - log1p(theta)
    ans[x < 0 | is.infinite(x)] <- log(0)  # 20141209 KaiH
  } else {
    ans <- theta^2 * (1 + x) * exp(-theta * x) / (1 + theta)
    ans[x < 0 | is.infinite(x)] <- 0  # 20141209 KaiH
  }
  ans[theta <= 0] <- NaN
  ans
}



plind <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {


  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ans <- log(-expm1(-theta * q + log1p(q / (1 + 1/theta))))
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- -expm1(-theta * q + log1p(q / (1 + 1/theta)))
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- -theta * q + log1p(q / (1 + 1/theta))
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      ans <- exp(-theta * q + log1p(q / (1 + 1/theta)))
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  }
  ans[theta <= 0] <- NaN
  ans
}






rlind <- function(n, theta) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
             stop("bad input for argument 'n'") else n



  ifelse(runif(use.n) < rep_len(1 / (1 + 1/theta), use.n),
         rexp(use.n, theta),
         rgamma(use.n, shape = 2, scale = 1 / theta))
}



 lindley <- function(link = "loge",
                     itheta = NULL, zero = NULL) {


  if (length(itheta) &&
      !is.Numeric(itheta, positive = TRUE))
    stop("argument 'itheta' must be > 0")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")




  new("vglmff",
  blurb = c("Lindley distribution f(y) = ",
            "theta^2 * (1 + y) * exp(-theta * y) / (1 + theta), ",
            "theta > 0, y > 0,\n\n",
            "Link:    ",
            namesof("theta", link, earg = earg), "\n\n",
            "Mean:     (theta + 2) / (theta * (theta + 1))\n",
            "Variance: (theta^2 + 4 * theta + 2) / (theta * (theta + 1))^2"),

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
         parameters.names = c("theta"),
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
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
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1  <- param.names("theta", ncoly)
    predictors.names <- namesof(mynames1, .link , earg = .earg , tag = FALSE)

    if (!length(etastart)) {
      wmeany <- colSums(y * w) / colSums(w) + 1/8


      theta.init <- 1 / (wmeany + 1)
      theta.init <- matrix(if (length( .itheta )) .itheta else
                           theta.init, n, M, byrow = TRUE)
      etastart <- theta2eta(theta.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .itheta = itheta ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    theta <- eta2theta(eta, .link , earg = .earg )
    (theta + 2) / (theta * (theta + 1))
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <- c(rep_len( .link , ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$M1 <- M1
    misc$itheta <- .itheta
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .link = link, .earg = earg,
            .itheta = itheta ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    theta <- eta2theta(eta, .link , earg = .earg )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dlind(x = y, theta = theta, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("lindley"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    theta <- eta2theta(eta, .link , earg = .earg )
    okay1 <- all(is.finite(theta)) && all(0 < theta)
    okay1
  }, list( .link = link, .earg = earg ))),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    theta <- eta2theta(eta, .link , earg = .earg )
    rlind(nsim * length(theta), theta = theta)
  }, list( .link = link, .earg = earg ))),



  deriv = eval(substitute(expression({
    M1 <- 1
    theta <- eta2theta(eta, .link , earg = .earg )

    dl.dtheta <- 2 / theta - 1 / (1 + theta) - y

    DTHETA.DETA <- dtheta.deta(theta, .link , earg = .earg )

    c(w) * dl.dtheta * DTHETA.DETA
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({

    ned2l.dtheta2 <- (theta^2 + 4 * theta + 2) / (theta * (1 + theta))^2

    c(w) * ned2l.dtheta2 * DTHETA.DETA^2
  }), list( .zero = zero ))))
}






dpoislindley <- function(x, theta, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if ( log.arg ) {
    ans <- 2 * log(theta) + log(theta + 2 + x) -
           (x+3) * log1p(theta)
    ans[(x != round(x)) | (x < 0)] <- log(0)
  } else {
    ans <- theta^2 * (theta + 2 + x) / (theta + 1)^(x+3)
    ans[(x != round(x)) | (x < 0)] <- 0
  }
  ans[ # !is.finite(theta) |
     (theta <= 0)] <- NA
  ans
}







dslash <- function(x, mu = 0, sigma = 1, log = FALSE,
                   smallno = .Machine$double.eps * 1000) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(sigma) || any(sigma <= 0))
    stop("argument 'sigma' must be positive")
  L <- max(length(x), length(mu), length(sigma))
  if (length(x)     != L) x     <- rep_len(x,     L)
  if (length(mu)    != L) mu    <- rep_len(mu,    L)
  if (length(sigma) != L) sigma <- rep_len(sigma, L)

  zedd <- (x-mu)/sigma
  if (log.arg) {
    ifelse(abs(zedd) < smallno,
           -log(2*sigma*sqrt(2*pi)),
           log1p(-exp(-zedd^2/2)) - log(sqrt(2*pi)*sigma*zedd^2))
  } else {
    ifelse(abs(zedd) < smallno,
           1/(2*sigma*sqrt(2*pi)),
           -expm1(-zedd^2/2)/(sqrt(2*pi)*sigma*zedd^2))
  }
}




pslash <- function(q, mu = 0, sigma = 1, very.negative = -10000,
                   lower.tail = TRUE, log.p = FALSE) {
  if (anyNA(q))
    stop("argument 'q' must have non-missing values")
  if (!is.Numeric(mu))
    stop("argument 'mu' must have finite and non-missing values")
  if (!is.Numeric(sigma, positive = TRUE))
    stop("argument 'sigma' must have positive finite non-missing values")
  if (!is.Numeric(very.negative, length.arg = 1) ||
     (very.negative >= 0))
    stop("argument 'very.negative' must be quite negative")

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  L <- max(length(q), length(mu), length(sigma))
  if (length(q)     != L) q     <- rep_len(q,     L)
  if (length(mu)    != L) mu    <- rep_len(mu,    L)
  if (length(sigma) != L) sigma <- rep_len(sigma, L)

  zedd <- (q - mu)/sigma
  ans <- as.numeric(q * NA)
  extreme.q <- FALSE
  for (ii in 1:L) {
    use.trick <- (-abs(zedd[ii]) <= very.negative)
    if (use.trick) {
      ans[ii] <- ifelse(zedd[ii] < 0, 0.0, 1.0)
      extreme.q <- TRUE
    } else
    if ((zedd[ii] >= very.negative) &&
         zedd[ii] <= 0.0) {
      temp2 <- integrate(dslash, lower = q[ii], upper = mu[ii],
                         mu = mu[ii], sigma = sigma[ii])
      if (temp2$message != "OK")
        warning("integrate() failed on 'temp2'")
      ans[ii] <- 0.5 - temp2$value
    } else {
      temp1 <- integrate(dslash, lower = mu[ii], upper =  q[ii],
                         mu = mu[ii], sigma = sigma[ii])
      if (temp1$message != "OK")
        warning("integrate() failed")
      ans[ii] <- 0.5 + temp1$value
    }
  }
  if (extreme.q)
    warning("returning 0 or 1 values for extreme values of argument 'q'")

  if (lower.tail) {
    if (log.p) log(ans) else ans
  } else {
    if (log.p) log1p(-ans) else -expm1(log(ans))
  }
}




rslash <- function (n, mu = 0, sigma = 1) {
  rnorm(n = n, mean = mu, sd = sigma) / runif(n = n)
}



slash.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 slash <- function(lmu = "identitylink", lsigma = "loge",
                   imu = NULL, isigma = NULL,
                   gprobs.y = ppoints(8),
                   nsimEIM = 250, zero = NULL,
                   smallno = .Machine$double.eps * 1000) {

  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lsigma <- as.list(substitute(lsigma))
  esigma <- link2list(lsigma)
  lsigma <- attr(esigma, "function.name")


  if (length(isigma) &&
      !is.Numeric(isigma, positive = TRUE))
    stop("argument 'isigma' must be > 0")



  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
    stop("argument 'nsimEIM' should be an integer greater than 50")

  if (!is.Numeric(gprobs.y, positive = TRUE) ||
      max(gprobs.y) >= 1)
    stop("bad input for argument 'gprobs.y'")
  if (!is.Numeric(smallno, positive = TRUE) ||
      smallno > 0.1)
    stop("bad input for argument 'smallno'")


  new("vglmff",
  blurb = c("Slash distribution\n\n",
         "Links:    ",
         namesof("mu",    lmu,    earg = emu,    tag = FALSE), ", ",
         namesof("sigma", lsigma, earg = esigma, tag = FALSE), "\n",
         paste(
         "1-exp(-(((y-mu)/sigma)^2)/2))/(sqrt(2*pi)*",
         "sigma*((y-mu)/sigma)^2)",
         "\ty!=mu",
         "\n1/(2*sigma*sqrt(2*pi))",
         "\t\t\t\t\t\t\ty=mu\n")),

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
         parameters.names = c("mu", "sigma"),
         lmu    = .lmu ,
         lsigma = .lsigma ,
         zero = .zero )
  }, list( .zero = zero, .lmu = lmu, .lsigma = lsigma ))),


  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <- c(
        namesof("mu",    .lmu ,    earg = .emu,    tag = FALSE),
        namesof("sigma", .lsigma , earg = .esigma, tag = FALSE))


    if (!length(etastart)) {

      slash.Loglikfun <- function(mu, y, x, w, extraargs) {
          sigma <- if (is.Numeric(.isigma)) .isigma else
            max(0.01,
               ((quantile(rep(y, w), prob = 0.75)/2)-mu)/qnorm(0.75))
          zedd <- (y-mu)/sigma
          sum(c(w) * ifelse(abs(zedd)<.smallno,
                         -log(2*sigma*sqrt(2*pi)),
                         log1p(-exp(-zedd^2/2)) -
                         log(sqrt(2*pi) * sigma * zedd^2)))
      }
      gprobs.y <- .gprobs.y
      mu.grid <- quantile(rep(y, w), probs = gprobs.y)
      mu.grid <- seq(mu.grid[1], mu.grid[2], length=100)
      mu.init <- if (length( .imu )) .imu else
                 grid.search(mu.grid, objfun = slash.Loglikfun,
                             y = y,  x = x, w = w)
      sigma.init <- if (is.Numeric(.isigma)) .isigma else
        max(0.01,
           ((quantile(rep(y, w), prob = 0.75)/2) -
                      mu.init) / qnorm(0.75))
      mu.init <- rep_len(mu.init, length(y))
      etastart <- matrix(0, n, 2)
      etastart[, 1] <- theta2eta(mu.init, .lmu , earg = .emu )
      etastart[, 2] <- theta2eta(sigma.init, .lsigma , earg = .esigma )
    }
  }), list( .lmu = lmu, .lsigma = lsigma,
            .imu = imu, .isigma = isigma,
            .emu = emu, .esigma = esigma,
            .gprobs.y = gprobs.y, .smallno = smallno))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
      NA * eta2theta(eta[, 1], link = .lmu , earg = .emu )
  }, list( .lmu = lmu, .emu = emu ))),
  last = eval(substitute(expression({
    misc$link <-    c("mu" = .lmu , "sigma" = .lsigma )

    misc$earg <- list("mu" = .emu , "sigma" = .esigma )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma, .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mu    <- eta2theta(eta[, 1], link = .lmu    , earg = .emu )
    sigma <- eta2theta(eta[, 2], link = .lsigma , earg = .esigma )
    zedd <- (y - mu) / sigma
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dslash(x = y, mu = mu, sigma = sigma, log = TRUE,
                               smallno = .smallno)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmu = lmu, .lsigma = lsigma,
           .emu = emu, .esigma = esigma, .smallno = smallno ))),
  vfamily = c("slash"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    mu    <- eta2theta(eta[, 1], link = .lmu    , earg = .emu    )
    sigma <- eta2theta(eta[, 2], link = .lsigma , earg = .esigma )

    okay1 <- all(is.finite(mu))    &&
             all(is.finite(sigma)) && all(0 < sigma)
    okay1
  }, list( .lmu = lmu, .lsigma = lsigma,
           .emu = emu, .esigma = esigma ))),






  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    mu    <- eta2theta(eta[, 1], link = .lmu    , earg = .emu )
    sigma <- eta2theta(eta[, 2], link = .lsigma , earg = .esigma )
    rslash(nsim * length(sigma), mu = mu, sigma = sigma)
  }, list( .lmu = lmu, .lsigma = lsigma,
           .emu = emu, .esigma = esigma, .smallno = smallno ))),




  deriv = eval(substitute(expression({
    mu    <- eta2theta(eta[, 1], link = .lmu    , earg = .emu    )
    sigma <- eta2theta(eta[, 2], link = .lsigma , earg = .esigma )

    dmu.deta    <- dtheta.deta(mu,    link = .lmu    , earg = .emu    )
    dsigma.deta <- dtheta.deta(sigma, link = .lsigma , earg = .esigma )

    zedd <- (y - mu) / sigma
    d3 <- deriv3(~ w * log(1 - exp(-(((y - mu) / sigma)^2) / 2)) -
                 log(sqrt(2 * pi) * sigma * ((y - mu) / sigma)^2),
                 c("mu", "sigma"))
    eval.d3 <- eval(d3)
    dl.dthetas <-  attr(eval.d3, "gradient")
    dl.dmu    <- dl.dthetas[, 1]
    dl.dsigma <- dl.dthetas[, 2]
    ind0 <- (abs(zedd) < .smallno)
    dl.dmu[ind0] <- 0
    dl.dsigma[ind0] <- -1 / sigma[ind0]
    c(w) * cbind(dl.dmu * dmu.deta, dl.dsigma * dsigma.deta)
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma, .smallno = smallno ))),
  weight = eval(substitute(expression({
    run.varcov <- 0
    ind1 <- iam(NA_real_, NA_real_, M = M, both = TRUE, diag = TRUE)
    sd3 <- deriv3(~ w * log(1 - exp(-(((ysim - mu) / sigma)^2) / 2))-
                  log(sqrt(2 * pi) * sigma * ((ysim - mu) / sigma)^2),
                  c("mu", "sigma"))
    for (ii in 1:( .nsimEIM )) {
      ysim <- rslash(n, mu = mu, sigma = sigma)
      seval.d3 <- eval(sd3)

      dl.dthetas <-  attr(seval.d3, "gradient")
      dl.dmu    <- dl.dthetas[, 1]
      dl.dsigma <- dl.dthetas[, 2]

      temp3 <- cbind(dl.dmu, dl.dsigma)
      run.varcov <- run.varcov + temp3[, ind1$row] * temp3[, ind1$col]
    }
    run.varcov <- run.varcov / .nsimEIM
    wz <- if (intercept.only)
        matrix(colMeans(run.varcov, na.rm = FALSE),
               n, ncol(run.varcov), byrow = TRUE) else run.varcov
    dthetas.detas <- cbind(dmu.deta, dsigma.deta)
    wz <- wz * dthetas.detas[, ind1$row] * dthetas.detas[, ind1$col]
    c(w) * wz
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma,
            .nsimEIM = nsimEIM, .smallno = smallno ))))
}




dnefghs <- function(x, tau, log = FALSE) {



  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(tau))
  if (length(x)   != N) x   <- rep_len(x,   N)
  if (length(tau) != N) tau <- rep_len(tau, N)

  logdensity <- log(sin(pi*tau)) + (1-tau)*x - log(pi) - log1pexp(x)
  logdensity[tau < 0] <- NaN
  logdensity[tau > 1] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}



 nefghs <- function(link = "logit",
                    itau = NULL, imethod = 1) {

  if (length(itau) &&
      !is.Numeric(itau, positive = TRUE) ||
      any(itau >= 1))
    stop("argument 'itau' must be in (0, 1)")

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
       imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  new("vglmff",
  blurb = c("Natural exponential family generalized hyperbolic ",
            "secant distribution\n",
            "f(y) = sin(pi*tau)*exp((1-tau)*y)/(pi*(1+exp(y))\n\n",
            "Link:    ",
            namesof("tau", link, earg = earg), "\n\n",
            "Mean:     pi / tan(pi * tau)\n"),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("tau"),
         ltau = .link )
  }, list( .link = link ))),

  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      namesof("tau", .link , earg = .earg , tag = FALSE)


    if (!length(etastart)) {
      wmeany <- if ( .imethod == 1) weighted.mean(y, w) else
                median(rep(y, w))
      if (abs(wmeany) < 0.01)
        wmeany <- 0.01
      tau.init <- atan(pi / wmeany) / pi + 0.5
      tau.init[tau.init < 0.03] <- 0.03
      tau.init[tau.init > 0.97] <- 0.97
      tau.init <- rep_len(if (length( .itau )) .itau else tau.init, n)
      etastart <- theta2eta(tau.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg,
            .itau = itau,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    tau <- eta2theta(eta, .link , earg = .earg )
    pi / tan(pi * tau)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c(tau = .link )

    misc$earg <- list(tau = .earg )

    misc$expected <- TRUE
    misc$imethod <- .imethod
  }), list( .link = link, .earg = earg,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    tau <- eta2theta(eta, .link , earg = .earg )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dnefghs(x = y, tau = tau, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("nefghs"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    tau <- eta2theta(eta, .link , earg = .earg )
    okay1 <- all(is.finite(tau)) && all(0 < tau)
    okay1
  }, list( .link = link, .earg = earg ))),


  deriv = eval(substitute(expression({
    tau <- eta2theta(eta, .link , earg = .earg )
    dl.dtau <- pi / tan(pi * tau) - y
    dtau.deta <- dtheta.deta(tau, .link , earg = .earg )
    w * dl.dtau * dtau.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    ned2l.dtau2 <- (pi / sin(pi * tau))^2
    wz <- ned2l.dtau2 * dtau.deta^2
    c(w) * wz
  }), list( .link = link ))))
}




dlogF <- function(x, shape1, shape2, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)




  logdensity <- shape1*x - lbeta(shape1, shape2) -
                (shape1 + shape2) * log1pexp(x)

  logdensity[is.infinite(x)] <- -Inf  # 20141209 KaiH

  if (log.arg) logdensity else exp(logdensity)
}




 logF <- function(lshape1 = "loge", lshape2 = "loge",
                  ishape1 = NULL, ishape2 = 1,
                  imethod = 1) {

  if (length(ishape1) &&
      !is.Numeric(ishape1, positive = TRUE))
    stop("argument 'ishape1' must be positive")
  if ( # length(ishape2) &&
      !is.Numeric(ishape2, positive = TRUE))
    stop("argument 'ishape2' must be positive")


  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")


  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
      stop("argument 'imethod' must be 1 or 2")

  new("vglmff",
  blurb = c("log F distribution\n",
            "f(y) = exp(-shape2 * y) / (beta(shape1, shape2) * ",
            "(1 + exp(-y))^(shape1 + shape2))\n\n",
            "Link:    ",
            namesof("shape1", lshape1, earg = eshape1), ", ",
            namesof("shape2", lshape2, earg = eshape2), "\n\n",
            "Mean:     digamma(shape1) - digamma(shape2)"),
  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("shape1", "shape2"),
         lshape1 = .lshape1 ,
         lshape2 = .lshape2 ,
         imethod = .imethod )
  }, list( .imethod = imethod, .lshape1 = lshape1, .lshape2 = lshape2 ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <- c(
      namesof("shape1", .lshape1 , earg = .eshape1 , tag = FALSE),
      namesof("shape2", .lshape2 , earg = .eshape2 , tag = FALSE))


    if (!length(etastart)) {
      wmeany <- if ( .imethod == 1) weighted.mean(y, w) else
                median(rep(y, w))


      shape1.init <- shape2.init <- rep_len( .ishape2 , n)
      shape1.init <- if (length( .ishape1)) rep_len( .ishape1, n) else {
                index1 <- (y > wmeany)
                shape1.init[ index1] <- shape2.init[ index1] + 1/1
                shape1.init[!index1] <- shape2.init[!index1] - 1/1
                shape1.init <- pmax(shape1.init, 1/8)
                shape1.init
              }
      etastart <-
          cbind(theta2eta(shape1.init, .lshape1 , earg = .eshape1 ),
                theta2eta(shape2.init, .lshape2 , earg = .eshape2 ))
    }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .ishape1 = ishape1, .ishape2 = ishape2,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .eshape2 )
    digamma(shape1) - digamma(shape2)
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1 = .lshape1 , shape2 = .lshape2 )

    misc$earg <- list(shape1 = .eshape1 , shape2 = .eshape2 )

    misc$expected <- TRUE
    misc$imethod <- .imethod
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .eshape2 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dlogF(x = y, shape1 = shape1,
                              shape2 = shape2, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  vfamily = c("logF"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .eshape2 )
    okay1 <- all(is.finite(shape1)) && all(0 < shape1) &&
             all(is.finite(shape2)) && all(0 < shape2)
    okay1
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),













  deriv = eval(substitute(expression({
    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .eshape2 )

    tmp888 <- digamma(shape1 + shape2) - log1pexp(-y)
    dl.dshape1 <- tmp888 - digamma(shape1)
    dl.dshape2 <- tmp888 - digamma(shape2) - y

    dshape1.deta <- dtheta.deta(shape1, .lshape1 , earg = .eshape1 )
    dshape2.deta <- dtheta.deta(shape2, .lshape2 , earg = .eshape2 )

    c(w) * cbind(dl.dshape1 * dshape1.deta,
                 dl.dshape2 * dshape2.deta)
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  weight = eval(substitute(expression({
    tmp888 <- trigamma(shape1 + shape2)
    ned2l.dshape12 <- trigamma(shape1) - tmp888
    ned2l.dshape22 <- trigamma(shape2) - tmp888
    ned2l.dshape1shape2 <- -tmp888

    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M = M)] <- ned2l.dshape12 * dshape1.deta^2
    wz[, iam(2, 2, M = M)] <- ned2l.dshape22 * dshape2.deta^2
    wz[, iam(1, 2, M = M)] <- ned2l.dshape1shape2 * dshape1.deta * dshape2.deta

    c(w) * wz
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))))
}







dbenf <- function(x, ndigits = 1, log = FALSE) {
  if (!is.Numeric(ndigits, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


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
  if (!is.Numeric(ndigits, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
             stop("bad input for argument 'n'") else n
  myrunif <- runif(use.n)

  ans <- rep_len(lowerlimit, use.n)
  for (ii in (lowerlimit+1):upperlimit) {
      indexTF <- (pbenf(ii-1, ndigits = ndigits) < myrunif) &
                 (myrunif <= pbenf(ii, ndigits = ndigits))
      ans[indexTF] <- ii
  }
  ans
}



pbenf <- function(q, ndigits = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (!is.Numeric(ndigits, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)

  ans <- q * NA
  floorq <- floor(q)
  indexTF <- is.finite(q) & (floorq >= lowerlimit)

  if (ndigits == 1) {
    if (lower.tail) {
      if (log.p) {
        ans[indexTF] <- log(log10(1 + floorq[indexTF]))
        ans[q <  lowerlimit ] <- -Inf
        ans[q >= upperlimit] <- 0
      } else {
        ans[indexTF] <- log10(1 + floorq[indexTF])
        ans[q <  lowerlimit] <- 0
        ans[q >= upperlimit] <- 1
      }
    } else {
      if (log.p) {
        ans[indexTF] <- log1p(-log10(1 + floorq[indexTF]))
        ans[q <  lowerlimit] <- 0
        ans[q >= upperlimit] <- -Inf
      } else {
        ans[indexTF] <- log10(10 / (1 + floorq[indexTF]))
        ans[q <  lowerlimit] <- 1
        ans[q >= upperlimit] <- 0
      }
    }
  } else {
    if (lower.tail) {
      if (log.p) {
        ans[indexTF] <- log(log10((1 + floorq[indexTF])/10))
        ans[q <  lowerlimit ] <- -Inf
        ans[q >= upperlimit] <- 0
      } else {
        ans[indexTF] <- log10((1 + floorq[indexTF])/10)
        ans[q <  lowerlimit] <- 0
        ans[q >= upperlimit] <- 1
     }
    } else {
      if (log.p) {
        ans[indexTF] <- log(log10(100/(1 + floorq[indexTF])))
        ans[q <  lowerlimit] <- 0
        ans[q >= upperlimit] <- -Inf
      } else {
        ans[indexTF] <- log10(100/(1 + floorq[indexTF]))
        ans[q <  lowerlimit] <- 1
        ans[q >= upperlimit] <- 0
      }
    }
  }
  ans
}



if (FALSE)
qbenf <- function(p, ndigits = 1) {

  if (!is.Numeric(ndigits, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  bad <- !is.na(p) & !is.nan(p) & ((p < 0) | (p > 1))
  if (any(bad))
    stop("bad input for argument 'p'")

  ans <- rep_len(lowerlimit, length(p))
  for (ii in (lowerlimit+1):upperlimit) {
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




qbenf <- function(p, ndigits = 1,
                  lower.tail = TRUE, log.p = FALSE) {
  if (!is.Numeric(ndigits, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (log.p) {
    bad <- ((p > 0) | is.na(p) | is.nan(p))
  } else {
    bad <- ((p < 0) | (p > 1) | is.na(p) | is.nan(p))
  }
  if (any(bad))
    stop("bad input for argument 'p'")

  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  ans <- rep_len(lowerlimit, length(p))

  if (lower.tail) {
    for (ii in (lowerlimit+1):upperlimit) {
      indexTF <- is.finite(p) &
                 (pbenf(ii-1, ndigits = ndigits,
                        lower.tail = lower.tail, log.p = log.p) < p) &
              (p <= pbenf(ii, ndigits = ndigits,
                          lower.tail = lower.tail, log.p = log.p))
      ans[indexTF] <- ii
    }
  } else {  ## when lower.tail = F, pbenf(ii-1) >= p & pben(ii) < p
    for (ii in (lowerlimit+1):upperlimit) {
      indexTF <- is.finite(p) &
                 (pbenf(ii-1, ndigits = ndigits,
                        lower.tail = lower.tail, log.p = log.p) >= p) &
                 (p > pbenf(ii, ndigits = ndigits,
                            lower.tail = lower.tail, log.p = log.p))
      ans[indexTF] <- ii
    }
  }

  if (lower.tail) {
    if (log.p) {
      ans[p > 0] <- NaN
      ans[p == -Inf] <- lowerlimit
    } else {
      ans[p < 0] <- NaN
      ans[p == 0] <- lowerlimit
      ans[p == 1] <- upperlimit
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ans[p > 0] <- NaN
      ans[p == -Inf] <- upperlimit
    } else {
      ans[p < 0] <- NaN
      ans[p == 0] <- upperlimit
      ans[p == 1] <- lowerlimit
      ans[p > 1] <- NaN
    }
  }
  ans
}


















 truncgeometric <-
  function(upper.limit = Inf,  # lower.limit = 1,  # Inclusive
           link = "logit", expected = TRUE,
           imethod = 1, iprob = NULL, zero = NULL) {

  if (is.finite(upper.limit) &&
      !is.Numeric(upper.limit, integer.valued = TRUE,
                  positive = TRUE))
    stop("bad input for argument 'upper.limit'")

  if (any(upper.limit < 0))
    stop("bad input for argument 'upper.limit'")



  if (!is.logical(expected) || length(expected) != 1)
    stop("bad input for argument 'expected'")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")



  uu.ll <- min(upper.limit)


  new("vglmff",
  blurb = c("Truncated geometric distribution ",
            "(P[Y=y] =\n",
            "     ",
            "prob * (1 - prob)^y / [1-(1-prob)^",
             uu.ll+1, "], y = 0,1,...,",
             uu.ll, ")\n",
            "Link:     ",
            namesof("prob", link, earg = earg), "\n",
            "Mean:     mu = 1 / prob - 1 ",
            ifelse(is.finite(upper.limit),
                   paste("- (", upper.limit+1, ") * (1 - prob)^",
                         upper.limit+1, " / (1 - ",
                         "(1 - prob)^", upper.limit+1, ")", sep = ""),
                         "")),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = .expected ,
         imethod = .imethod ,
         multipleResponses = TRUE,
         parameters.names = c("prob"),
         upper.limit = .upper.limit ,
         zero = .zero )
  }, list( .zero = zero,
           .expected = expected,
           .imethod = imethod,
           .upper.limit = upper.limit ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
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
    extra$M1 <- M1
    M <- M1 * ncoly
    extra$upper.limit <- matrix( .upper.limit , n, ncoly, byrow = TRUE)

    if (any(y > extra$upper.limit))
      stop("some response values greater than argument 'upper.limit'")


    mynames1 <- param.names("prob", ncoly)
    predictors.names <- namesof(mynames1, .link , earg = .earg , tag = FALSE)


    if (!length(etastart)) {
      prob.init <- if ( .imethod == 2)
                      1 / (1 + y + 1/16) else
                  if ( .imethod == 3)
                      1 / (1 + apply(y, 2, median) + 1/16) else
                      1 / (1 + colSums(y * w) / colSums(w) + 1/16)

      if (!is.matrix(prob.init))
        prob.init <- matrix(prob.init, n, M, byrow = TRUE)


      if (length( .iprob ))
        prob.init <- matrix( .iprob , n, M, byrow = TRUE)


        etastart <- theta2eta(prob.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg,
            .upper.limit = upper.limit,
            .imethod = imethod, .iprob = iprob ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob <- eta2theta(eta, .link , earg = .earg )
    QQQ <- 1 - prob
    upper.limit <- extra$upper.limit
    tmp1 <- QQQ^(upper.limit+1)
    answer <- 1 / prob - 1 - (upper.limit+1) * tmp1 / (1 - tmp1)
    answer[!is.finite(answer)] <- 1 / prob[!is.finite(answer)] - 1
    answer
  }, list( .link = link, .earg = earg ))),

  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <- c(rep_len( .link , ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$M1 <- M1
    misc$multipleResponses <- TRUE
    misc$expected <- .expected
    misc$imethod <- .imethod
    misc$iprob <- .iprob
  }), list( .link = link, .earg = earg,
            .iprob = iprob,
            .upper.limit = upper.limit,
            .expected = expected, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    prob <- eta2theta(eta, .link , earg = .earg )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      upper.limit <- extra$upper.limit
      ll.elts <- c(w) * (dgeom(x = y, prob = prob, log = TRUE) -
                         log1p(-(1.0 - prob)^(1 + upper.limit)))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("truncgeometric"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    prob <- eta2theta(eta, .link , earg = .earg )
    okay1 <- all(is.finite(prob)) && all(0 < prob & prob < 1)
    okay1
  }, list( .link = link, .earg = earg ))),
  deriv = eval(substitute(expression({
    prob <- eta2theta(eta, .link , earg = .earg )
    sss <- upper.limit <- extra$upper.limit  # Is a matrix

    QQQ <- 1 - prob
    tmp1 <- QQQ^(upper.limit + 1)
    dl.dprob <- 1 / prob  + (0 - y) / (1 - prob) -
                (1+upper.limit) * QQQ^(upper.limit - 0) / (1 - tmp1)
    dl.dprob[!is.finite(upper.limit)] <-  1 / prob[!is.finite(upper.limit)] +
      (0 - y[!is.finite(upper.limit)]) / (1 - prob[!is.finite(upper.limit)])


    dprobdeta <- dtheta.deta(prob, .link , earg = .earg )
    c(w) * cbind(dl.dprob * dprobdeta)
  }), list( .link = link, .earg = earg,
            .upper.limit = upper.limit,
            .expected = expected ))),
  weight = eval(substitute(expression({

    eim.oim.fun <- function(mu.y, sss)
      ifelse(is.finite(sss),
             1/prob^2 + (0 + mu.y) / QQQ^2 - (1+sss) *
             ((sss-0) * QQQ^(sss-1) / (1 - tmp1) +
             (1+sss) * QQQ^(2*sss) / (1 - tmp1)^2),
             1 / (prob^2 * (1 - prob)))


    ned2l.dprob2 <- if ( .expected ) {
      eim.oim.fun(mu, sss)
    } else {
      eim.oim.fun(y, sss)
    }
    wz <- ned2l.dprob2 * dprobdeta^2
    if ( !( .expected ))
      wz <- wz - dl.dprob * d2theta.deta2(prob, .link , earg = .earg )
    c(w) * wz
  }), list( .link = link, .earg = earg,
            .expected = expected ))))
}









 betaff <-
  function(A = 0, B = 1,
           lmu = "logit",
           lphi = "loge",
           imu = NULL, iphi = NULL,  # imethod = 1,
           gprobs.y = ppoints(8),  # (1:9)/10,
           gphi  = exp(-3:5)/4,
           zero = NULL) {



  if (!is.Numeric(A, length.arg = 1) ||
      !is.Numeric(B, length.arg = 1) || A >= B)
    stop("A must be < B, and both must be of length one")

  stdbeta <- (A == 0 && B == 1)


  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")



  lphi <- as.list(substitute(lphi))
  ephi <- link2list(lphi)
  lphi <- attr(ephi, "function.name")



  if (length(imu) && (!is.Numeric(imu, positive = TRUE) ||
     any(imu <= A) || any(imu >= B)))
    stop("bad input for argument 'imu'")
  if (length(iphi) && !is.Numeric(iphi, positive = TRUE))
    stop("bad input for argument 'iphi'")


  new("vglmff",
  blurb = c("Beta distribution parameterized by mu and a ",
            "precision parameter\n",
            if (stdbeta) paste("f(y) = y^(mu*phi-1) * (1-y)^((1-mu)*phi-1)",
            "/ beta(mu*phi,(1-mu)*phi),\n",
            "      0<y<1, 0<mu<1, phi>0\n\n") else
            paste("f(y) = (y-",A,")^(mu1*phi-1) * (",B,
            "-y)^(((1-mu1)*phi)-1) / \n(beta(mu1*phi,(1-mu1)*phi) * (",
            B, "-", A, ")^(phi-1)),\n",
            A," < y < ",B, ", ", A," < mu < ",B,
            ", mu = ", A, " + ", (B-A), " * mu1",
            ", phi > 0\n\n", sep = ""),
            "Links:    ",
            namesof("mu",  lmu,  earg = emu),  ", ",
            namesof("phi", lphi, earg = ephi)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M)
  }), list( .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("mu", "phi"),
         A = .A ,
         B = .B ,
         zero = .zero )
  }, list( .zero = zero,
           .A = A, .B = B ))),

  initialize = eval(substitute(expression({
    if (min(y) <= .A || max(y) >= .B)
      stop("data not within (A, B)")


    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    extra$A <- .A  # Needed for @validparams
    extra$B <- .B


    predictors.names <- c(namesof("mu",  .lmu ,  .emu , short = TRUE),
                          namesof("phi", .lphi , .ephi, short = TRUE))
    if (!length(etastart)) {
      NOS <- 1
      muu.init <-
      phi.init <- matrix(NA_real_, n, NOS)
      gprobs.y <- .gprobs.y
      gphi <- if (length( .iphi )) .iphi else .gphi

      betaff.Loglikfun <- function(muu, phi, y, x, w, extraargs) {
        zedd <- (y   - extraargs$A) / ( extraargs$B - extraargs$A)
        m1u  <- (muu - extraargs$A) / ( extraargs$B - extraargs$A)
        shape1 <- phi * m1u
        shape2 <- (1 - m1u) * phi
        sum(c(w) * (dbeta(x = zedd, shape1, shape2, log = TRUE) -
                    log(abs( extraargs$B - extraargs$A ))))
      }


      for (jay in 1:NOS) {  # For each response 'y_jay'... do:
        gmuu <- if (length( .imu )) .imu else quantile(y[, jay], probs = gprobs.y)


        try.this <-
          grid.search2(gmuu, gphi,
                       objfun = betaff.Loglikfun,
                       y = y[, jay],
                       w = w[, jay],
                       extraargs = list(A = .A , B = .B ),
                       ret.objfun = TRUE)  # Last value is the loglik
        muu.init[, jay] <-  try.this["Value1"]
        phi.init[, jay] <-  try.this["Value2"]
      }  # for (jay ...)


if (FALSE) {
      mu.init <- if (is.Numeric( .imu )) .imu else {
                   if ( .imethod == 1) weighted.mean(y, w) else
                      (y + weighted.mean(y, w)) / 2
                 }
      mu1.init <- (mu.init - .A ) / ( .B - .A )  # In (0,1)
      phi.init <- if (is.Numeric( .iphi )) .iphi else
         max(0.01, -1 + ( .B - .A )^2 * mu1.init*(1-mu1.init)/var(y))
  }



      etastart <- matrix(0, n, 2)
      etastart[, 1] <- theta2eta(muu.init, .lmu  , earg = .emu  )
      etastart[, 2] <- theta2eta(phi.init, .lphi , earg = .ephi )
    }
  }), list( .lmu = lmu, .lphi = lphi, .imu = imu, .iphi = iphi,
            .A = A, .B = B, .emu = emu, .ephi = ephi,
            .gprobs.y = gprobs.y, .gphi = gphi  ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
     mu <- eta2theta(eta[, 1], .lmu , .emu )
     mu
  }, list( .lmu = lmu, .emu = emu, .A = A, .B = B))),
  last = eval(substitute(expression({
    misc$link <-    c(mu = .lmu , phi = .lphi )
    misc$earg <- list(mu = .emu , phi = .ephi )
    misc$limits <- c( .A , .B )
    misc$stdbeta <- .stdbeta
  }), list( .lmu = lmu, .lphi = lphi, .A = A, .B = B,
            .emu = emu, .ephi = ephi,
            .stdbeta = stdbeta ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mu  <- eta2theta(eta[, 1], .lmu  , earg = .emu  )
    phi <- eta2theta(eta[, 2], .lphi , earg = .ephi )
    m1u <- if ( .stdbeta ) mu else (mu - .A ) / ( .B - .A )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      shape1 <- phi * m1u
      shape2 <- (1 - m1u) * phi
      zedd <- (y - .A) / ( .B - .A)
      ll.elts <-
        c(w) * (dbeta(x = zedd, shape1 = shape1, shape2 = shape2,
                      log = TRUE) -
                log( abs( .B - .A )))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmu = lmu, .lphi = lphi, .A = A, .B = B,
           .emu = emu, .ephi = ephi,
           .stdbeta = stdbeta ))),
  vfamily = "betaff",
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    mu  <- eta2theta(eta[, 1], .lmu  , .emu  )
    phi <- eta2theta(eta[, 2], .lphi , .ephi )
    okay1 <- all(is.finite(mu )) && all(extra$A < mu & mu < extra$B) &&
             all(is.finite(phi)) && all(0 < phi)
    okay1
  }, list( .lmu = lmu, .lphi = lphi, .A = A, .B = B,
           .emu = emu, .ephi = ephi ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")

    eta <- predict(object)
    mu  <- eta2theta(eta[, 1], .lmu  , earg = .emu  )
    phi <- eta2theta(eta[, 2], .lphi , earg = .ephi )
    m1u <- if ( .stdbeta ) mu else (mu - .A ) / ( .B - .A )
    shape1 <- phi * m1u
    shape2 <- (1 - m1u) * phi
    .A + ( .B - .A ) *
    rbeta(nsim * length(shape1), shape1 = shape1, shape2 = shape2)
  }, list( .lmu = lmu, .lphi = lphi, .A = A, .B = B,
           .emu = emu, .ephi = ephi,
           .stdbeta = stdbeta ))),





  deriv = eval(substitute(expression({
    mu  <- eta2theta(eta[, 1], .lmu  , .emu  )
    phi <- eta2theta(eta[, 2], .lphi , .ephi )
    m1u <- if ( .stdbeta ) mu else (mu - .A) / ( .B - .A)
    dmu.deta <- dtheta.deta(mu, .lmu , .emu )
    dmu1.dmu <- 1 / ( .B - .A )
    dphi.deta <- dtheta.deta(phi, .lphi , .ephi )
    temp1 <- m1u*phi
    temp2 <- (1-m1u)*phi
    if ( .stdbeta ) {
      dl.dmu1 <- phi*(digamma(temp2) - digamma(temp1) + log(y) - log1p(-y))
      dl.dphi <- digamma(phi) - mu*digamma(temp1) - (1-mu)*digamma(temp2) +
          mu*log(y) + (1-mu)*log1p(-y)
    } else {
      dl.dmu1 <- phi*(digamma(temp2) - digamma(temp1) +
                     log(y-.A) - log( .B-y))
      dl.dphi <- digamma(phi) - m1u*digamma(temp1) -
                (1-m1u)*digamma(temp2) +
                m1u*log(y-.A) + (1-m1u)*log( .B-y) - log( .B -.A)
    }
      c(w) * cbind(dl.dmu1 * dmu1.dmu * dmu.deta,
                   dl.dphi * dphi.deta)
  }), list( .lmu = lmu, .lphi = lphi,
            .emu = emu, .ephi = ephi,
            .A = A, .B = B,
            .stdbeta = stdbeta ))),
  weight = eval(substitute(expression({
    ned2l.dmu12 <- (trigamma(temp1) + trigamma(temp2)) * phi^2
    ned2l.dphi2 <- -trigamma(phi) + trigamma(temp1) * m1u^2 +
                    trigamma(temp2) * (1-m1u)^2
    ned2l.dmu1phi <- temp1 * trigamma(temp1) - temp2 * trigamma(temp2)
    wz <- matrix(NA_real_, n, dimm(M))
    wz[, iam(1, 1, M)] <- ned2l.dmu12 * dmu1.dmu^2 * dmu.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dphi2 * dphi.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dmu1phi * dmu1.dmu * dmu.deta * dphi.deta
    c(w) * wz
  }), list( .A = A, .B = B ))))
}





 betaR <-
  function(lshape1 = "loge", lshape2 = "loge",
           i1 = NULL, i2 = NULL, trim = 0.05,
           A = 0, B = 1, parallel = FALSE, zero = NULL) {

  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")

  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")


  if (length( i1 ) && !is.Numeric( i1, positive = TRUE))
    stop("bad input for argument 'i1'")
  if (length( i2 ) && !is.Numeric( i2, positive = TRUE))
    stop("bad input for argument 'i2'")

  if (!is.Numeric(A, length.arg = 1) ||
     !is.Numeric(B, length.arg = 1) ||
     A >= B)
    stop("A must be < B, and both must be of length one")

  stdbeta <- (A == 0 && B == 1)  # stdbeta == T iff standard beta distn



  new("vglmff",
  blurb = c("Two-parameter Beta distribution ",
            "(shape parameters parameterization)\n",
            if (stdbeta)
            paste("y^(shape1-1) * (1-y)^(shape2-1) / B(shape1,shape2),",
            "0 <= y <= 1, shape1>0, shape2>0\n\n") else
            paste("(y-",A,")^(shape1-1) * (",B,
            "-y)^(shape2-1) / [B(shape1,shape2) * (",
            B, "-", A, ")^(shape1+shape2-1)], ",
             A," <= y <= ",B," shape1>0, shape2>0\n\n", sep = ""),
            "Links:    ",
            namesof("shape1", lshape1, earg = eshape1),  ", ",
            namesof("shape2", lshape2, earg = eshape2)),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints, apply.int  = TRUE)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M)
  }), list( .parallel = parallel, .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         A  = .A,
         B  = .B,
         multipleResponses = FALSE,
         zero = .zero )
  }, list( .A = A, .B = B,
           .zero = zero ))),
  initialize = eval(substitute(expression({
    if (min(y) <= .A || max(y) >= .B)
      stop("data not within (A, B)")

    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")


    w.y.check(w = w, y = y)


    predictors.names <-
        c(namesof("shape1", .lshape1 , earg = .eshape1 , short = TRUE),
          namesof("shape2", .lshape2 , earg = .eshape2 , short = TRUE))

    if (!length(etastart)) {
      mu1d <- mean(y, trim = .trim )
      uu <- (mu1d - .A) / ( .B - .A)
      DD <- ( .B - .A)^2
      pinit <- max(0.01, uu^2 * (1 - uu) * DD / var(y) - uu)
      qinit <- max(0.01, pinit * (1 - uu) / uu)
      etastart <- matrix(0, n, 2)
      etastart[, 1] <- theta2eta( pinit, .lshape1 , earg = .eshape1 )
      etastart[, 2] <- theta2eta( qinit, .lshape2 , earg = .eshape2 )
    }
    if (is.Numeric( .i1 ))
      etastart[, 1] <- theta2eta( .i1 , .lshape1 , earg = .eshape1 )
    if (is.Numeric( .i2 ))
      etastart[, 2] <- theta2eta( .i2 , .lshape2 , earg = .eshape2 )
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .i1 = i1, .i2 = i2, .trim = trim, .A = A, .B = B,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shapes <- cbind(eta2theta(eta[, 1], .lshape1 , earg = .eshape1 ),
                    eta2theta(eta[, 2], .lshape2 , earg = .eshape2 ))
    .A + ( .B - .A ) * shapes[, 1] / (shapes[, 1] + shapes[, 2])
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1 = .lshape1 , shape2 = .lshape2 )
    misc$earg <- list(shape1 = .eshape1 , shape2 = .eshape2 )
    misc$limits <- c( .A , .B )
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .A = A, .B = B,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    shapes <- cbind(eta2theta(eta[, 1], .lshape1 , earg = .eshape1 ),
                    eta2theta(eta[, 2], .lshape2 , earg = .eshape2 ))
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      zedd <- (y - .A ) / ( .B - .A )
      ll.elts <-
        c(w) * (dbeta(x = zedd, shape1 = shapes[, 1],
                                shape2 = shapes[, 2],
                      log = TRUE) - log( abs( .B - .A )))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  vfamily = "betaR",
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shapes <- cbind(eta2theta(eta[, 1], .lshape1 , earg = .eshape1 ),
                    eta2theta(eta[, 2], .lshape2 , earg = .eshape2 ))
    okay1 <- all(is.finite(shapes)) && all(0 < shapes)
    okay1
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")

    eta <- predict(object)
    shapes <- cbind(eta2theta(eta[, 1], .lshape1 , earg = .eshape1 ),
                    eta2theta(eta[, 2], .lshape2 , earg = .eshape2 ))
    .A + ( .B - .A ) *
    rbeta(nsim * length(shapes[, 1]),
          shape1 = shapes[, 1], shape2 = shapes[, 2])
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),



  deriv = eval(substitute(expression({
    shapes <- cbind(eta2theta(eta[, 1], .lshape1 , earg = .eshape1 ),
                    eta2theta(eta[, 2], .lshape2 , earg = .eshape2 ))

    dshapes.deta <-
      cbind(dtheta.deta(shapes[, 1], .lshape1 , earg = .eshape1),
            dtheta.deta(shapes[, 2], .lshape2 , earg = .eshape2))

    dl.dshapes <- cbind(log(y - .A ), log( .B - y)) -
                  digamma(shapes) +
                  digamma(shapes[, 1] + shapes[, 2]) - log( .B - .A )

    c(w) * dl.dshapes * dshapes.deta
  }), list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  weight = expression({
    trig.sum <- trigamma(shapes[, 1] + shapes[, 2])
    ned2l.dshape12 <- trigamma(shapes[, 1]) - trig.sum
    ned2l.dshape22 <- trigamma(shapes[, 2]) - trig.sum
    ned2l.dshape1shape2 <- -trig.sum
    wz <- matrix(NA_real_, n, dimm(M))  # dimm(M) == 3
    wz[, iam(1, 1, M)] <- ned2l.dshape12      * dshapes.deta[, 1]^2
    wz[, iam(2, 2, M)] <- ned2l.dshape22      * dshapes.deta[, 2]^2
    wz[, iam(1, 2, M)] <- ned2l.dshape1shape2 * dshapes.deta[, 1] *
                                                dshapes.deta[, 2]
    c(w) * wz
  }))
}





 betaprime <-
  function(lshape = "loge", ishape1 = 2, ishape2 = NULL, zero = NULL) {

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  new("vglmff",
  blurb = c("Beta-prime distribution\n",
            "y^(shape1-1) * (1+y)^(-shape1-shape2) / Beta(shape1,shape2),",
            " y>0, shape1>0, shape2>0\n\n",
            "Links:    ",
            namesof("shape1", lshape, earg = eshape),  ", ",
            namesof("shape2", lshape, earg = eshape), "\n",
            "Mean:     shape1/(shape2-1) provided shape2>1"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M)
  }), list( .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("shape1", "shape2"),
         lshape1 = .lshape ,
         lshape2 = .lshape ,
         zero = .zero )
  }, list( .zero = zero, .lshape = lshape ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
      c(namesof("shape1", .lshape , earg = .eshape , short = TRUE),
        namesof("shape2", .lshape , earg = .eshape , short = TRUE))
    if (is.numeric( .ishape1) && is.numeric( .ishape2 )) {
      vec <- c( .ishape1, .ishape2 )
      vec <- c(theta2eta(vec[1], .lshape , earg = .eshape ),
               theta2eta(vec[2], .lshape , earg = .eshape ))
      etastart <- matrix(vec, n, 2, byrow = TRUE)
    }
    if (!length(etastart)) {
      init1 <- if (length( .ishape1 ))
        rep_len( .ishape1 , n) else rep_len(1, n)
      init2 <- if (length( .ishape2 ))
        rep_len( .ishape2 , n) else 1 + init1 / (y + 0.1)
      etastart <-
        matrix(theta2eta(c(init1, init2), .lshape , earg = .eshape ),
               n, 2, byrow = TRUE)
    }
  }), list( .lshape = lshape, .eshape = eshape,
            .ishape1 = ishape1, .ishape2 = ishape2 ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    shapes <- eta2theta(eta, .lshape , earg = .eshape )
    ifelse(shapes[, 2] > 1, shapes[, 1] / (shapes[, 2] - 1), NA)
  }, list( .lshape = lshape, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1 = .lshape , shape2 = .lshape )
    misc$earg <- list(shape1 = .eshape , shape2 = .eshape )
  }), list( .lshape = lshape, .eshape = eshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shapes <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * ((shapes[, 1]-1) * log(y) -
                 lbeta(shapes[, 1], shapes[, 2]) -
                (shapes[, 2] + shapes[, 1]) * log1p(y))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = "betaprime",
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shapes <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shapes)) && all(0 < shapes)
    okay1
  }, list( .lshape = lshape, .eshape = eshape ))),
  deriv = eval(substitute(expression({
    shapes <- eta2theta(eta, .lshape , earg = .eshape )
    dshapes.deta <- dtheta.deta(shapes, .lshape , earg = .eshape )
    dl.dshapes <- cbind(log(y) - log1p(y) - digamma(shapes[, 1]) +
                        digamma(shapes[, 1] + shapes[, 2]),
                        - log1p(y) - digamma(shapes[, 2]) +
                        digamma(shapes[, 1] + shapes[, 2]))
    c(w) * dl.dshapes * dshapes.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = expression({
    temp2 <- trigamma(shapes[, 1] + shapes[, 2])
    ned2l.dshape12 <- trigamma(shapes[, 1]) - temp2
    ned2l.dshape22 <- trigamma(shapes[, 2]) - temp2
    ned2l.dshape1shape2 <- -temp2

    wz <- matrix(NA_real_, n, dimm(M))  #3=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.dshape12 * dshapes.deta[, 1]^2
    wz[, iam(2, 2, M)] <- ned2l.dshape22 * dshapes.deta[, 2]^2
    wz[, iam(1, 2, M)] <- ned2l.dshape1shape2 *
                          dshapes.deta[, 1] * dshapes.deta[, 2]

    c(w) * wz
  }))
}  # betaprime










 zoabetaR <-
  function(lshape1 = "loge", lshape2 = "loge",
           lpobs0 = "logit", lpobs1 = "logit",
           ishape1 = NULL, ishape2 = NULL, trim = 0.05,
           type.fitted = c("mean", "pobs0", "pobs1", "beta.mean"),
           parallel.shape = FALSE,
           parallel.pobs = FALSE,
           zero = NULL) {

  A <- 0
  B <- 1

  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")

  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")

  lprobb0 <- as.list(substitute(lpobs0))
  eprobb0 <- link2list(lprobb0)
  lprobb0 <- attr(eprobb0, "function.name")

  lprobb1 <- as.list(substitute(lpobs1))
  eprobb1 <- link2list(lprobb1)
  lprobb1 <- attr(eprobb1, "function.name")

  if (length( ishape1 ) && !is.Numeric( ishape1, positive = TRUE))
    stop("bad input for argument 'ishape1'")
  if (length( ishape2 ) && !is.Numeric( ishape2, positive = TRUE))
    stop("bad input for argument 'ishape2'")

  if (!is.Numeric(A, length.arg = 1) ||
      !is.Numeric(B, length.arg = 1) ||
     A >= B)
    stop("A must be < B, and both must be of length one")

  stdbeta <- (A == 0 && B == 1)  # stdbeta == TRUE iff standard beta distn



  type.fitted <- match.arg(type.fitted,
                   c("mean", "pobs0", "pobs1", "beta.mean"))[1]


  new("vglmff",
  blurb = c("Standard Beta distribution with 0- and \n",
            "1-inflation ",
            "(shape parameters parameterization)\n",
            if (stdbeta)
            paste("y^(shape1-1) * (1-y)^(shape2-1) / beta(shape1,shape2),",
            "0 <= y <= 1, shape1>0, shape2>0\n\n") else
            paste("(y-",A,")^(shape1-1) * (",B,
            "-y)^(shape2-1) / [beta(shape1,shape2) * (",
            B, "-", A, ")^(shape1+shape2-1)], ",
             A," <= y <= ",B," shape1>0, shape2>0, ",
            "0 < pobs0 < 1, 0 < pobs1 < 1 \n\n", sep = ""),
            "Links:    ",
            namesof("shape1", lshape1, earg = eshape1),  ", ",
            namesof("shape2", lshape2, earg = eshape2),  ", ",
            namesof("pobs0",  lprobb0, earg = eprobb0),  ", ",
            namesof("pobs1",  lprobb1, earg = eshape1)),



  constraints = eval(substitute(expression({


    constraints.orig <- constraints


    if (is.logical( .parallel.probb ) && .parallel.probb &&
        (cind0[1] + cind1[1] <= 1))
      warning("argument 'parallel.pobs' specified when there is only ",
              "one of 'pobs0' and 'pobs1'")




    cmk.s <- kronecker(matrix(1, NOS, 1), rbind(1, 1, 0, 0))
    cmk.S <- kronecker(diag(NOS), rbind(diag(2), 0*diag(2)))
    con.s <- cm.VGAM(cmk.s, x = x,
                     bool = .parallel.shape ,  # Same as .parallel.b
                     constraints = constraints.orig,
                     apply.int = TRUE,
                     cm.default           = cmk.S,
                     cm.intercept.default = cmk.S)



    cmk.p <- kronecker(matrix(1, NOS, 1), rbind(0, 0, 1, 1))
    cmk.P <- kronecker(diag(NOS), rbind(0*diag(2), diag(2)))
    con.p <- cm.VGAM(cmk.p,
                     x = x,
                     bool = .parallel.probb ,  #
                     constraints = constraints.orig,
                     apply.int = TRUE,
                     cm.default           = cmk.P,
                     cm.intercept.default = cmk.P)

    con.use <- con.s
    for (klocal in seq_along(con.s)) {
      con.use[[klocal]] <-
        cbind(con.s[[klocal]], con.p[[klocal]])
 # Delete rows that are not needed:
      if (!cind0[1]) {
        con.use[[klocal]] <- (con.use[[klocal]])[c(TRUE, TRUE, FALSE, TRUE), ]
      }
      if (!cind1[1]) {
        con.use[[klocal]] <- (con.use[[klocal]])[c(TRUE, TRUE, TRUE, FALSE), ]
      }
      col.delete <- apply(con.use[[klocal]], 2, function(HkCol) all(HkCol == 0))
      con.use[[klocal]] <- (con.use[[klocal]])[, !col.delete]
    }


    constraints <- con.use



   constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = M1)

  }), list( .parallel.shape = parallel.shape,
            .parallel.probb = parallel.pobs,
            .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(M1 = NA,  # Either 3 or 4, data-dependent
         Q1 = 1,
         A  = .A ,
         B  = .B ,
         expected = TRUE,
         multipleResponses = TRUE,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .A = A, .B = B,
           .type.fitted = type.fitted,
           .zero = zero ))),
  initialize = eval(substitute(expression({
    if (min(y) <  .A || max(y)  > .B)
      stop("data not within [A, B]")

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    ncoly <- NOS <- ncol(y)
    if (ncoly > 1 && !( .stdbeta ))
      stop("can only input multiple responses with the standard beta")

    cind0 <- colSums(ind0 <- y == 0) > 0
    cind1 <- colSums(ind1 <- y == 1) > 0
    if (!any(cind0 | cind1))
      stop("no 0s or 1s in the responses to perform 0- and/or ",
           "1-inflation! ",
           "Try using betaff() or betaR() instead.")

    if (ncoly > 1 && !all(cind0 == cind0[1]) &&  # FALSE &&
                     !all(cind0 == cind0[1]))
      stop("with multiple responses, cannot have 0-inflation in ",
           "some responses and 1-inflation in other responses")
    M1 <- 2 + cind0[1] + cind1[1]  # 4 when there is both 0 & 1-inflation
    M <- M1 * NOS

    mynames1 <- param.names("shape1", ncoly)
    mynames2 <- param.names("shape2", ncoly)
    mynames3 <- param.names("pobs0",  ncoly)
    mynames4 <- param.names("pobs1",  ncoly)
    predictors.names <-
        c(namesof(mynames1, .lshape1 , earg = .eshape1 , short = TRUE),
          namesof(mynames2, .lshape2 , earg = .eshape2 , short = TRUE),
          if (cind0[1])
            namesof(mynames3, .lprobb0 , earg = .eprobb0 , short = TRUE)
          else NULL,
          if (cind1[1])
            namesof(mynames4, .lprobb1 , earg = .eprobb1 , short = TRUE)
          else NULL)[
          interleave.VGAM(M, M1 = M1)]

    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)
    extra$M1          <- M1  # Determined from the data
    extra$cind0       <- cind0
    extra$cind1       <- cind1

    if (!length(etastart)) {

      p0init <- matrix(colMeans(ind0), n, ncoly, byrow = TRUE)
      p1init <- matrix(colMeans(ind1), n, ncoly, byrow = TRUE)


      mu1d <- matrix(NA_real_, n, NOS)
      for (jay in 1:ncoly) {
        yy <- y[, jay]
        yy <- yy[ .A < yy & yy < .B ]
        mu1d[, jay] <- weighted.mean(yy, trim = .trim )
      }
      uu <- (mu1d - .A ) / ( .B - .A )
      DD <- ( .B - .A )^2
      p.init <- if (is.Numeric( .ishape1 ))
        matrix( .ishape1 , n, ncoly, byrow = TRUE) else
        uu^2 * (1 - uu) * DD / var(yy) - uu
      p.init[p.init < 0.01] <- 0.01
      q.init <- if (is.Numeric( .ishape2 ))
        matrix( .ishape2 , n, ncoly, byrow = TRUE) else
        p.init * (1 - uu) / uu
      q.init[q.init < 0.01] <- 0.01
      etastart <- cbind(
        theta2eta(p.init, .lshape1 , earg = .eshape1 ),
        theta2eta(q.init, .lshape2 , earg = .eshape2 ),
        if (cind0[1])
        theta2eta(p0init, .lprobb0 , earg = .eprobb0 )
        else NULL,
        if (cind1[1])
        theta2eta(p1init, .lprobb1 , earg = .eprobb1 )
        else NULL)[,
          interleave.VGAM(M, M1 = M1)]
    }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .lprobb0 = lprobb0, .lprobb1 = lprobb1,
            .eprobb0 = eprobb0, .eprobb1 = eprobb1,
            .ishape1 = ishape1, .ishape2 = ishape2,
            .trim = trim, .A = A, .B = B,
            .type.fitted = type.fitted,
            .stdbeta = stdbeta ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    M1 <- extra$M1
    cind0 <- extra$cind0
    cind1 <- extra$cind1
    NOS <- ncol(eta) / M1
    shape1 <- eta2theta(eta[, c(TRUE, rep(FALSE, M1 - 1)), drop = FALSE],
                        .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE,
                                rep(FALSE, M1 - 2)), drop = FALSE],
                        .lshape2 , earg = .eshape2 )
    probb0 <- if (cind0[1])
      eta2theta(eta[, c(FALSE, FALSE, TRUE,
                        if (cind1[1]) FALSE else NULL), drop = FALSE],
                .lprobb0 , earg = .eprobb0 ) else 0
    probb1 <- if (cind1[1])
      eta2theta(eta[, c(FALSE, FALSE,
                        if (cind0[1]) FALSE else NULL, TRUE),
                    drop = FALSE],
                .lprobb1 , earg = .eprobb1 ) else 0

    type.fitted <- match.arg(extra$type.fitted,
                     c("mean", "pobs0", "pobs1", "beta.mean"))[1]

    ans <-
      switch(type.fitted,
             "mean"      = (1 - probb0) * shape1 / (shape1 + shape2) +
                                probb1  * shape2 / (shape1 + shape2),
      "beta.mean" = shape1/(shape1+shape2),  # zz Mux by (1-pobs0-pobs1)??
             "pobs0"     = probb0,
             "pobs1"     = probb1)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B,
           .eshape1 = eshape1, .eshape2 = eshape2,
           .lprobb0 = lprobb0, .lprobb1 = lprobb1,
           .eprobb0 = eprobb0, .eprobb1 = eprobb1 ))),
  last = eval(substitute(expression({
    misc$link <- rep_len( c( .lshape1 , .lshape2 ,
                             if (cind0[1]) .lprobb0 else NULL,
                             if (cind1[1]) .lprobb1 else NULL), M)
    names(misc$link) <- c(mynames1, mynames2,
                          if (cind0[1]) mynames3 else NULL,
                          if (cind1[1]) mynames4 else NULL)[
                        interleave.VGAM(M, M1 = M1)]

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    jay <- 1
    while (jay <= M) {
      misc$earg[[jay]] <- .eshape1
      jay <- jay + 1
      misc$earg[[jay]] <- .eshape2
      jay <- jay + 1
      if (cind0[1]) {
        misc$earg[[jay]] <- .eprobb0
        jay <- jay + 1
      }
      if (cind1[1]) {
        misc$earg[[jay]] <- .eprobb1
        jay <- jay + 1
      }
    }

    misc$supportlimits <- c( .A , .B )
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .lprobb0 = lprobb0, .lprobb1 = lprobb1,
            .eprobb0 = eprobb0, .eprobb1 = eprobb1,
            .A = A, .B = B ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    M1 <- 4
    M1 <- extra$M1
    cind0 <- extra$cind0
    cind1 <- extra$cind1
    NOS <- ncol(eta) / M1
    shape1 <- eta2theta(eta[, c(TRUE, rep(FALSE, M1 - 1)), drop = FALSE],
                        .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE,
                                rep(FALSE, M1 - 2)), drop = FALSE],
                        .lshape2 , earg = .eshape2 )
    probb0 <- if (cind0[1])
      eta2theta(eta[, c(FALSE, FALSE, TRUE,
                        if (cind1[1]) FALSE else NULL), drop = FALSE],
                .lprobb0 , earg = .eprobb0 ) else 0
    probb1 <- if (cind1[1])
      eta2theta(eta[, c(FALSE, FALSE,
                        if (cind0[1]) FALSE else NULL, TRUE),
                    drop = FALSE],
                .lprobb1 , earg = .eprobb1 ) else 0

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      zedd <- (y - .A ) / ( .B - .A )
      ll.elts <-
        c(w) * (dzoabeta(x = zedd, shape1 = shape1, shape2 = shape2,
                         pobs0 = probb0, pobs1 = probb1,
                      log = TRUE) - log( abs( .B - .A )))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B,
           .eshape1 = eshape1, .eshape2 = eshape2,
           .lprobb0 = lprobb0, .lprobb1 = lprobb1,
           .eprobb0 = eprobb0, .eprobb1 = eprobb1 ))),
  vfamily = "zoabetaR",









  validparams = eval(substitute(function(eta, y, extra = NULL) {
    M1 <- 4
    M1 <- extra$M1
    cind0 <- extra$cind0
    cind1 <- extra$cind1
    NOS <- ncol(eta) / M1
    shape1 <- eta2theta(eta[, c(TRUE, rep(FALSE, M1 - 1)), drop = FALSE],
                        .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE,
                                rep(FALSE, M1 - 2)), drop = FALSE],
                        .lshape2 , earg = .eshape2 )
    probb0 <- if (cind0[1])
      eta2theta(eta[, c(FALSE, FALSE, TRUE,
                        if (cind1[1]) FALSE else NULL), drop = FALSE],
                .lprobb0 , earg = .eprobb0 ) else 0.5
    probb1 <- if (cind1[1])
      eta2theta(eta[, c(FALSE, FALSE,
                        if (cind0[1]) FALSE else NULL, TRUE), drop = FALSE],
                .lprobb1 , earg = .eprobb1 ) else 0.5



    okay1 <- all(is.finite(shape1)) && all(0 < shape1) &&
             all(is.finite(shape2)) && all(0 < shape2) &&
             all(is.finite(probb0)) && all(0 < probb0 & probb0 < 1) &&
             all(is.finite(probb1)) && all(0 < probb1 & probb1 < 1)
    okay1
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B,
           .eshape1 = eshape1, .eshape2 = eshape2,
           .lprobb0 = lprobb0, .lprobb1 = lprobb1,
           .eprobb0 = eprobb0, .eprobb1 = eprobb1 ))),



  deriv = eval(substitute(expression({
    M1 <- 4
    M1 <- extra$M1
    cind0 <- extra$cind0
    cind1 <- extra$cind1
    NOS <- ncol(eta) / M1
    shape1 <- eta2theta(eta[, c(TRUE, rep(FALSE, M1 - 1)), drop = FALSE],
                        .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE,
                                rep(FALSE, M1 - 2)), drop = FALSE],
                        .lshape2 , earg = .eshape2 )
    probb0 <- if (cind0[1])
      eta2theta(eta[, c(FALSE, FALSE, TRUE,
                        if (cind1[1]) FALSE else NULL), drop = FALSE],
                .lprobb0 , earg = .eprobb0 ) else 0
    probb1 <- if (cind1[1])
      eta2theta(eta[, c(FALSE, FALSE,
                        if (cind0[1]) FALSE else NULL, TRUE), drop = FALSE],
                .lprobb1 , earg = .eprobb1 ) else 0


    dshape1.deta <- dtheta.deta(shape1, .lshape1 , earg = .eshape1 )
    dshape2.deta <- dtheta.deta(shape2, .lshape2 , earg = .eshape2 )
    dprobb0.deta <- dtheta.deta(probb0, .lprobb0 , earg = .eprobb0 )
    dprobb1.deta <- dtheta.deta(probb1, .lprobb1 , earg = .eprobb1 )

    index0 <- y == 0
    index1 <- y == 1
    indexi <- !index0 & !index1  # In the interior, i.e., (0, 1)
    dig.sum <- digamma(shape1 + shape2)
    QQ <- 1 - probb0 - probb1

    if (cind0[1]) {
      dl.dprobb0 <- -1 / QQ
      dl.dprobb0[index0] <- 1 / probb0[index0]
      dl.dprobb0[index1] <- 0
    }

    if (cind1[1]) {
      dl.dprobb1 <- -1 / QQ
      dl.dprobb1[index0] <- 0
      dl.dprobb1[index1] <- 1 / probb1[index1]
    }

    dl.dshape1 <-    log(y) - digamma(shape1) + dig.sum
    dl.dshape2 <- log1p(-y) - digamma(shape2) + dig.sum
    dl.dshape1[!indexi] <- 0
    dl.dshape2[!indexi] <- 0

    myderiv <- c(w) *
      cbind(dl.dshape1 * dshape1.deta,
            dl.dshape2 * dshape2.deta,
            if (cind0[1]) dl.dprobb0 * dprobb0.deta else NULL,
            if (cind1[1]) dl.dprobb1 * dprobb1.deta else NULL)
    colnames(myderiv) <- NULL
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .lprobb0 = lprobb0, .lprobb1 = lprobb1,
            .eprobb0 = eprobb0, .eprobb1 = eprobb1 ))),
  weight = expression({
    trig.sum <- trigamma(shape1 + shape2)


ned2l.dshape12 <- (trigamma(shape1) - trig.sum) * QQ
ned2l.dshape22 <- (trigamma(shape2) - trig.sum) * QQ
ned2l.dprobb02 <- (1 - probb1) / (probb0 * QQ)
ned2l.dprobb12 <- (1 - probb0) / (probb1 * QQ)
ned2l.dshape1shape2 <- -trig.sum * QQ  # (1 - probb0 - probb0) zz
ned2l.dshape2probb0 <- 0
ned2l.dprobb0probb1 <- 1 / QQ
ned2l.dshape1probb0 <- 0
ned2l.dshape2probb1 <- 0
ned2l.dshape1probb1 <- 0

ned2l.dshape1probb0 <- 0
    wz <- array(c(c(w) * ned2l.dshape12 * dshape1.deta^2,
                  c(w) * ned2l.dshape22 * dshape2.deta^2,
   if (cind0[1])  c(w) * ned2l.dprobb02 * dprobb0.deta^2 else NULL,
   if (cind1[1])  c(w) * ned2l.dprobb12 * dprobb1.deta^2 else NULL,
                  c(w) * ned2l.dshape1shape2 * dshape1.deta * dshape2.deta,
   if (cind0[1])  c(w) * ned2l.dshape2probb0 * dshape2.deta * dprobb0.deta,
                  c(w) * ned2l.dprobb0probb1 * dprobb0.deta * dprobb1.deta,
   if (cind0[1])  c(w) * ned2l.dshape1probb0 * dshape1.deta * dprobb0.deta,
   if (cind1[1])  c(w) * ned2l.dshape2probb1 * dshape2.deta * dprobb1.deta,
   if (cind1[1])  c(w) * ned2l.dshape1probb1 * dshape1.deta * dprobb1.deta),
                dim = c(n, M / M1, M1*(M1+1)/2))

    wz <- arwz2wz(wz, M = M, M1 = M1) # tridiagonal but unexploited here
    wz
  }))
}  # zoabetaR






dtopple <- function(x, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  L <- max(length(x), length(shape))
  if (length(x)     != L) x     <- rep_len(x,     L)
  if (length(shape) != L) shape <- rep_len(shape, L)

  logdensity <- rep_len(log(0), L)
  xok <- (0 <= x) & (x <= 1)
  logdensity[xok] <-
    log(2) + log(shape[xok]) + log1p(-x[xok]) +
    (shape[xok] - 1) * (log(x[xok]) + log(2) + log1p(-x[xok]/2))
  logdensity[shape >= 1] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}



ptopple <- function(q, shape, lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  if (lower.tail) {
    if (log.p) {
      ans <- shape * (log(q) + log(2) + log1p(-q/2))
      ans[q <= 0 ] <- -Inf
      ans[q >= 1] <- 0
    } else {
      ans <- (q * (2 - q))^shape
      ans[q <= 0] <- 0
      ans[q >= 1] <- 1
    }
  } else {
    if (log.p) {
      ans <- log1p(-(q * (2 - q))^shape)
      ans[q <= 0] <- 0
      ans[q >= 1] <- -Inf
    } else {
      ans <- exp(log1p(-(q * (2 - q))^shape))
      ans[q <= 0] <- 1
      ans[q >= 1] <- 0
    }
  }
  ans[shape <= 0] <- NaN
  ans[shape >= 1] <- NaN
  ans
}



qtopple <- function(p, shape) {
  ans <- -expm1(0.5 * log1p(-p^(1/shape)))
  ans[shape <= 0] <- NaN
  ans[shape >= 1] <- NaN
  ans
}



rtopple <- function(n, shape) {
  qtopple(runif(n), shape)
}





 topple <- function(lshape = "logit", zero = NULL,
                    gshape = ppoints(8)) {


  lshape <- as.list(substitute(lshape))  # orig
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")



  new("vglmff",
  blurb = c("Topp-Leone distribution F(y;shape) = (y * (2 - y))^shape, ",
            "0 < y < 1, 0 < shape < 1\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M)
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
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (any(y >= 1))
      stop("response must be in (0, 1)")


    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1  <- param.names("shape", ncoly)
    predictors.names <-
      namesof(mynames1, .lshape , earg = .eshape , tag = FALSE)


    if (!length(etastart)) {
      shape.init <- matrix(0, nrow(x), ncoly)
      gshape <- .gshape
      topple.Loglikfun <- function(shape, y, x = NULL, w, extraargs = NULL) {
        sum(c(w) * dtopple(x = y, shape = shape, log = TRUE))
      }

      for (jay in 1:ncoly) {
        shape.init[, jay] <- grid.search(gshape, objfun = topple.Loglikfun,
                                         y = y[, jay], w = w[, jay])
      }
      etastart <- theta2eta(shape.init, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .gshape = gshape,
            .eshape = eshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    1 - (gamma(1 + shape))^2 * 4^shape / gamma(2 * (1 + shape))
  }, list( .lshape = lshape,
           .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ilocal in 1:ncoly) {
      misc$earg[[ilocal]] <- .eshape
    }

    misc$link <- rep_len( .lshape , ncoly)
    names(misc$link) <- mynames1
  }), list( .lshape = lshape, .eshape = eshape ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dtopple(x = y, shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("topple"),
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
    rtopple(nsim * length(shape), shape = c(shape))
  }, list( .lshape = lshape,
           .eshape = eshape ))),


  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    dl.dshape <- 1 / shape + log(y) + log(2) + log1p(-y/2)
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = eval(substitute(expression({
    ned2l.dshape2 <- 1 / shape^2
    wz <- c(w) * ned2l.dshape2 * dshape.deta^2
    wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}




dzeta <- function(x, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  LLL <- max(length(shape), length(x))
  if (length(x)     != LLL) x     <- rep_len(x,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)

  ox <- !is.finite(x)
  zero <- ox | round(x) != x | x < 1
  ans <- rep_len(if (log.arg) log(0) else 0, LLL)
  if (any(!zero)) {
      if (log.arg) {
          ans[!zero] <- (-shape[!zero]-1)*log(x[!zero]) -
                        log(zeta(shape[!zero]+1))
      } else {
          ans[!zero] <- x[!zero]^(-shape[!zero]-1) / zeta(shape[!zero]+1)
      }
  }
  if (any(ox))
    ans[ox] <- if (log.arg) log(0) else 0
  ans[shape <= 0] <- NaN  # Added 20160617
  ans
}



 pzeta <- function(q, shape, lower.tail = TRUE) {


  LLL <- max(lenq <- length(q), lens <- length(shape))
  if (length(q)     != LLL) q     <- rep_len(q,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  ans <- rep_len(0, LLL)

  aa <- 12  # Same as Zeta.aux()
  qfloor <- floor(q)
  for (nn in 1:(aa-1))
    ans <- ans + as.numeric(nn <= qfloor) / nn^(shape+1)

  vecTF <- (aa-1 <= qfloor)
  if (lower.tail) {
    if (any(vecTF))
      ans[vecTF] <- zeta(shape[vecTF]+1) -
                    Zeta.aux(shape[vecTF]+1, qfloor[vecTF]+1)
  } else {
    ans <- zeta(shape+1) - ans
    if (any(vecTF))
      ans[vecTF] <- Zeta.aux(shape[vecTF]+1, qfloor[vecTF]+1)
  }
  ans / zeta(shape+1)
}  # pzeta




 qzeta <- function(p, shape) {

  LLL <- max(lenp <- length(p), lens <- length(shape))
  if (length(p)     != LLL) p     <- rep_len(p,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  ans <- rep_len(0, LLL)

  lo <- rep_len(1, LLL)
  approx.ans <- lo  # True at lhs
  hi <- 2 * lo + 10
  dont.iterate <- p == 1 | shape <= 0
  done <- p <= pzeta(hi, shape) | dont.iterate
  while (!all(done)) {
    hi.save <- hi[!done]
    hi[!done] <- 2 * lo[!done] + 10
    lo[!done] <- hi.save
    done[!done] <- (p[!done] <= pzeta(hi[!done], shape[!done]))
  }

  foo <- function(q, shape, p)
    pzeta(q, shape) - p

  lhs <- (p <= dzeta(1, shape)) | dont.iterate

  approx.ans[!lhs] <- bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                                      shape = shape[!lhs], p = p[!lhs])
  faa <- floor(approx.ans)
  ans <- ifelse(pzeta(faa, shape) < p & p <= pzeta(faa+1, shape), faa+1, faa)

  ans[p == 1] <- Inf
  ans[shape <= 0] <- NaN

  ans
}  # qzeta



rzeta <- function(n, shape) {
  qzeta(runif(n), shape)
}





 zetaff <-
    function(lshape = "loge",
             ishape = NULL,
             gshape = exp(-3:4)/4,
             zero = NULL) {


  if (length(ishape) && !is.Numeric(ishape, positive = TRUE))
    stop("argument 'ishape' must be > 0")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  new("vglmff",
  blurb = c("Zeta distribution ",
            "f(y) = 1/(y^(shape+1) zeta(shape+1)), shape>0, y = 1, 2,..\n\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape), "\n\n",
            "Mean:     zeta(shape) / zeta(shape+1), provided shape>1\n",
            "Variance: zeta(shape-1) / zeta(shape+1) - mean^2, provided shape>2"),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         multipleResponses = TRUE,
         parameters.names = "shape",
         zero = .zero ,
         lshape = .lshape )
  }, list( .lshape = lshape,
           .zero = zero ))),
  initialize = eval(substitute(expression({

   temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              Is.positive.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    ncoly <- ncol(y)
    mynames1 <- param.names("shape", ncoly)
    predictors.names <-
      namesof(mynames1, .lshape , earg = .eshape , tag = FALSE)

    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    if (!length(etastart)) {
      zetaff.Loglikfun <- function(shape, y, x, w, extraargs) {
        sum(c(w) * dzeta(x = y, shape, log = TRUE))
      }


      gshape <- .gshape
      if (!length( .ishape )) {
        shape.init <- matrix(NA_real_, n, M, byrow = TRUE)
        for (jay in 1:ncoly) {
          shape.init[, jay] <- grid.search(gshape, objfun = zetaff.Loglikfun,
                                           y = y[, jay], x = x, w = w[, jay])
        }
      } else {
        shape.init <- matrix( .ishape , n, M, byrow = TRUE)
      }
      etastart <- theta2eta(shape.init, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .eshape = eshape,
            .ishape = ishape, .gshape = gshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    ans <- pp <- eta2theta(eta, .lshape , earg = .eshape )
    ans[pp > 1] <- zeta(pp[pp > 1]) / zeta(pp[pp > 1] + 1)
    ans[pp <= 1] <- NA
    ans
  }, list( .lshape = lshape, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$link <- rep_len( .lshape , ncoly)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (jay in 1:ncoly) {
      misc$earg[[jay]] <- .eshape
    }

  }), list( .lshape = lshape, .eshape = eshape ))),
  loglikelihood = eval(substitute( function(mu, y, w, residuals = FALSE,
             eta, extra = NULL, summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dzeta(x = y, shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("zetaff"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape)
    okay1
  }, list( .lshape = lshape, .eshape = eshape ))),
  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .eshape )

    fred0 <- zeta(shape+1)
    fred1 <- zeta(shape+1, deriv = 1)
    dl.dshape <- -log(y) - fred1 / fred0

    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = expression({
    NOS <- NCOL(y)
    nd2l.dshape2 <- zeta(shape + 1, deriv = 2) / fred0 - (fred1/fred0)^2
    wz <- nd2l.dshape2 * dshape.deta^2
    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }))
}



gharmonic2 <- function(n, shape = 1) {



  if (!is.Numeric(n, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'n'")

  LLL <- max(length(n), length(shape))
  if (length(n)     != LLL) n     <- rep_len(n,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)

  aa <- 12
  ans <- rep_len(0, LLL)
  for (ii in 1:aa)
     ans <- ans + as.numeric(ii <= n) / ii^shape

  vecTF <- (aa < n)
  if (any(vecTF))
    ans[vecTF] <- zeta(shape[vecTF]) - Zeta.aux(shape[vecTF], 1 + n[vecTF])
  ans
}



gharmonic <- function(n, shape = 1, deriv = 0) {





  if (!is.Numeric(n, integer.valued = TRUE, positive = TRUE))
      stop("bad input for argument 'n'")
  if (!is.Numeric(deriv, length.arg = 1, integer.valued = TRUE) ||
      deriv < 0)
    stop("bad input for argument 'deriv'")

  lognexponent <- deriv
  sign <- ifelse(deriv %% 2 == 0, 1, -1)

  ans <-
  if (length(n) == 1 && length(shape) == 1) {
    if (lognexponent != 0) sum(log(1:n)^lognexponent * (1:n)^(-shape)) else
      sum((1:n)^(-shape))
  } else {
    LEN <- max(length(n), length(shape))
    n <- rep_len(n, LEN)
    ans <- shape <- rep_len(shape, LEN)
    if (lognexponent != 0) {
      for (ii in 1:LEN)
        ans[ii] <- sum(log(1:n[ii])^lognexponent * (1:n[ii])^(-shape[ii]))
    } else {
      for (ii in 1:LEN)
        ans[ii] <- sum((1:n[ii])^(-shape[ii]))
    }
    ans
  }
  sign * ans
}




dzipf <- function(x, N, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if (!is.Numeric(x))
    stop("bad input for argument 'x'")
  if (!is.Numeric(N, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'N'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")
  nn <- max(length(x), length(N), length(shape))
  if (length(x)     != nn) x     <- rep_len(x,     nn)
  if (length(N)     != nn) N     <- rep_len(N,     nn)
  if (length(shape) != nn) shape <- rep_len(shape, nn)

  ox <- !is.finite(x)
  zero <- ox | round(x) != x | x < 1 | x > N
  ans <- (if (log.arg) log(0) else 0) * x
  if (any(!zero))
    if (log.arg) {
      ans[!zero] <- (-shape[!zero]) * log(x[!zero]) -
                    log(gharmonic2(N[!zero], shape[!zero]))
    } else {
      ans[!zero] <- x[!zero]^(-shape[!zero]) / gharmonic2(N[!zero],
                                                          shape[!zero])
    }
  ans
}





pzipf <- function(q, N, shape, log.p = FALSE) {

  if (!is.Numeric(N, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'N'")

  nn <- max(length(q), length(N), length(shape))
  if (length(q)     != nn) q     <- rep_len(q,     nn)
  if (length(N)     != nn) N     <- rep_len(N,     nn)
  if (length(shape) != nn) shape <- rep_len(shape, nn)
    oq <- !is.finite(q)
  dont.iterate <- shape <= 0
    zeroOR1 <- oq | q < 1 | N <= q | dont.iterate
    floorq <- floor(q)
    ans <- 0 * floorq
    ans[oq | q >= N] <- 1
    if (any(!zeroOR1))
      ans[!zeroOR1] <- gharmonic2(floorq[!zeroOR1], shape[!zeroOR1]) /
                       gharmonic2(     N[!zeroOR1], shape[!zeroOR1])

    ans[shape <= 0] <- NaN

    if (log.p) log(ans) else ans
}






qzipf <- function(p, N, shape) {
  if (!is.Numeric(p))
    stop("bad input for argument 'p'")
  if (!is.Numeric(N, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'N'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")

  nn <- max(length(p), length(N), length(shape))
  if (length(p)     != nn) p     <- rep_len(p,     nn)
  if (length(N)     != nn) N     <- rep_len(N,     nn)
  if (length(shape) != nn) shape <- rep_len(shape, nn)

  a <- rep_len(1, nn)
  b <- rep_len(N, nn)
  approx.ans <- a  # True at lhs

  foo <- function(q, N, shape, p)
    pzipf(q, N, shape) - p

  dont.iterate <- p == 1 | shape <= 0
  lhs <- (p <= dzipf(1, N, shape)) | dont.iterate

  approx.ans[!lhs] <-
    bisection.basic(foo, a[!lhs], b[!lhs], shape = shape[!lhs], tol = 1/16,
                    p = p[!lhs], N = N[!lhs])
  faa <- floor(approx.ans)
  ans <- ifelse(pzipf(faa, N, shape) < p & p <= pzipf(faa+1, N, shape),
                faa+1, faa)

  ans[shape <= 0] <- NaN
  ans[p == 1] <- N

  ans
}  # qzipf



rzipf <- function(n, N, shape) {
  qzipf(runif(n), N, shape)
}







 zipf <- function(N = NULL, lshape = "loge", ishape = NULL) {

  if (length(N) &&
    (!is.Numeric(N, positive = TRUE,
                 integer.valued = TRUE, length.arg = 1) ||
      N <= 1))
    stop("bad input for argument 'N'")
  enteredN <- length(N)

  if (length(ishape) && !is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' must be > 0")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  new("vglmff",
  blurb = c("Zipf distribution f(y;s) = y^(-s) / sum((1:N)^(-s)),",
            " s > 0, y = 1, 2,...,N",
            ifelse(enteredN, paste(" = ", N, sep = ""), ""),
            "\n\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape),
            "\n\n",
            "Mean:    gharmonic(N, shape-1) / gharmonic(N, shape)"),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         multipleResponses = FALSE,
         parameters.names = "shape",
         N = enteredN,
         lshape = .lshape )
  }, list( .lshape = lshape,
           .enteredN = enteredN
         ))),
  initialize = eval(substitute(expression({


    w.y.check(w = w, y = y,
              Is.integer.y = TRUE)


    predictors.names <- namesof("shape", .lshape , earg = .eshape ,
                                tag = FALSE)

    NN <- .N
    if (!is.Numeric(NN, length.arg = 1,
                    positive = TRUE, integer.valued = TRUE))
        NN <- max(y)
    if (max(y) > NN)
      stop("maximum of the response is greater than argument 'N'")
    if (any(y < 1))
      stop("all response values must be in 1, 2, 3,...,N( = ", NN,")")
    extra$N <- NN

    if (!length(etastart)) {
      llfun <- function(shape, y, N, w) {
        sum(c(w) * dzipf(x = y, N = extra$N, shape = shape, log = TRUE))
      }
      shape.init <- if (length( .ishape )) .ishape else
        getInitVals(gvals = seq(0.1, 3, length.out = 19),
                    llfun = llfun,
                    y = y, N = extra$N, w = w)
      shape.init <- rep_len(shape.init, length(y))
      if ( .lshape == "loglog") shape.init[shape.init <= 1] <- 1.2
      etastart <- theta2eta(shape.init, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .eshape = eshape, .ishape = ishape, .N = N ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    gharmonic2(extra$N, shape = shape - 1) / gharmonic2(extra$N, shape = shape)
  }, list( .lshape = lshape, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$expected <- FALSE
    misc$link <-    c(shape = .lshape)
    misc$earg <- list(shape = .eshape )
    misc$N <- extra$N
  }), list( .lshape = lshape, .eshape = eshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dzipf(x = y, N = extra$N, shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("zipf"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape)
    okay1
  }, list( .lshape = lshape, .eshape = eshape ))),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    extra <- object@extra
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    rzipf(nsim * length(shape), N = extra$N, shape = shape)
  }, list( .lshape = lshape, .eshape = eshape ))),



  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    fred1 <- gharmonic(extra$N, shape, deriv = 1)

    fred0 <- gharmonic2(extra$N, shape)

    dl.dshape <- -log(y) - fred1 / fred0
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    d2shape.deta2 <- d2theta.deta2(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = expression({
    d2l.dshape <- gharmonic(extra$N, shape, deriv = 2) / fred0 -
                  (fred1/fred0)^2
    wz <- c(w) * (dshape.deta^2 * d2l.dshape - d2shape.deta2 * dl.dshape)
    wz
  }))
}







deflat.limit.oizeta  <- function(shape) {
  if (any(shape <= 0))
    stop("argument 'shape' must be positive")
  ans <- -dzeta(1, shape) / pzeta(1, shape, lower.tail = FALSE)
  ans
}



doizeta <- function(x, shape, pstr1 = 0, log = FALSE) {



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
                        dzeta(x[ index1], shape[ index1]))
    ans[!index1] <- log1p(-pstr1[!index1]) +
                         dzeta(x[!index1], shape[!index1], log = TRUE)
  } else {
    ans[ index1] <-      pstr1[ index1] + (1 - pstr1[ index1]) *
                       dzeta(x[ index1], shape[ index1])
    ans[!index1] <- (1 - pstr1[!index1]) *
                       dzeta(x[!index1], shape[!index1])
  }


  ans[pstr1 < deflat.limit.oizeta(shape)] <- NaN
  ans[pstr1 > 1] <- NaN

  ans
}  # doizeta




poizeta <- function(q, shape, pstr1 = 0) {

  LLL <- max(length(q), length(shape), length(pstr1))
  if (length(q)     != LLL) q     <- rep_len(q,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oizeta(shape)

  ans <- pzeta(q, shape)  #, lower.tail = lower.tail, log.p = log.p
  ans <- ifelse(q < 1, 0, pstr1 + (1 - pstr1) * ans)

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN
  ans[shape <= 0] <- NaN

  ans
}  # poizeta





qoizeta <- function(p, shape, pstr1 = 0) {

  LLL <- max(length(p), length(shape), length(pstr1))
  if (length(p)     != LLL) p     <- rep_len(p,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oizeta(shape)

  ans[p <= pstr1] <- 1
  pindex <- (deflat.limit <= pstr1) & (pstr1 < p)
  ans[pindex] <-
    qzeta((p[pindex] - pstr1[pindex]) / (1 - pstr1[pindex]),
          shape = shape[pindex])

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN

  ans[p < 0] <- NaN
  ans[1 < p] <- NaN
  ans[shape <= 0] <- NaN

  ans
}  # qoizeta



roizeta <- function(n, shape, pstr1 = 0) {
  qoizeta(runif(n), shape, pstr1 = pstr1)
}





 oizeta <-
  function(lpstr1 = "logit", lshape = "loge",
           type.fitted = c("mean", "shape", "pobs1", "pstr1", "onempstr1"),
           ishape = NULL,
           gpstr1 = ppoints(8),
           gshape = exp((-3:3) / 4), # grid for finding shape.init
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
  blurb = c("One-inflated zeta regression\n\n",
            "Links:    ",
            namesof("pstr1",  lpstr1, earg = epstr1 ), ", ",
            namesof("shape", lshape, earg = eshape ), "\n",
            "Mean:     pstr1 + (1 - pstr1) * zeta(shape) / ",
                       "zeta(1 + shape), if shape > 1"),

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
              Is.nonnegative.y = TRUE,
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

    mynames1 <- param.names("pstr1",  ncoly)
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

      oizeta.Loglikfun <- function(pstr1, shape, y, x, w, extraargs) {
        sum(c(w) * doizeta(x = y, pstr1 = pstr1,
                           shape = shape, log = TRUE))
      }


      for (jay in 1:NOS) {  # For each response 'y_jay'... do:


        try.this <-
          grid.search2(gpstr1, gshape,
                       objfun = oizeta.Loglikfun,
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
            .gshape = gshape,
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
      Mean <- shape
      Mean[shape > 1] <-
        zeta(shape[shape > 1]) / zeta(1 + shape[shape > 1])
      Mean[shape <= 1] <- NA
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
      ll.elts <- c(w) * doizeta(x = y, pstr1 = pstr1, shape = shape,
                                log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),
  vfamily = c("oizeta"),



  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    roizeta(nsim * length(shape), shape = shape, pstr1 = pstr1)
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                        earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lshape ,
                        earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape) &&
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

    pmf1 <- dzeta(1, shape)
    onempmf1 <- 1 - pmf1  # dozeta(1, shape = shape, pstr1 = pstr1)
    pobs1 <- pstr1 + (1 - pstr1) * pmf1
    index1 <- as.matrix(y == 1)

    zeta0 <- zeta(shape + 1)
    zeta1 <- zeta(shape + 1, deriv = 1)
    zeta2 <- zeta(shape + 1, deriv = 2)

    dl.dpstr1 <- onempmf1 / pobs1
    dl.dpstr1[!index1] <- -1 / (1 - pstr1[!index1])

    dpmf1.dshape <- -zeta1 / zeta0^2

   d2pmf1.dshape2 <- (2 * zeta1^2 / zeta0 - zeta2) / zeta0^2

    dl.dshape <- (1 - pstr1) * dpmf1.dshape / pobs1  #
    dl.dshape[!index1] <- -log(y[!index1]) - zeta1[!index1] / zeta0[!index1]

    dpstr1.deta <- dtheta.deta(pstr1, .lpstr1 , earg = .epstr1 )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    myderiv <- c(w) * cbind(dl.dpstr1 * dpstr1.deta,
                            dl.dshape * dshape.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape ))),
  weight = eval(substitute(expression({

    LHS <- ((1 - pstr1) / pobs1) * dpmf1.dshape^2 - d2pmf1.dshape2
    RHS <- (zeta2 - zeta1^2 / zeta0) / zeta0
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
}  # oizeta








deflat.limit.oizipf  <- function(N, shape) {
  if (any(shape <= 0))
    stop("argument 'shape' must be positive")
  ans <- 1 / (1 - 1 / dzipf(1, N, shape))
  ans
}





doizipf <- function(x, N, shape, pstr1 = 0, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if (!is.Numeric(x))
    stop("bad input for argument 'x'")
  if (!is.Numeric(N, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'N'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")
  nn <- max(length(x), length(N), length(shape), length(pstr1))
  if (length(x)    != nn) x     <- rep_len(x,     nn)
  if (length(N)    != nn) N     <- rep_len(N,     nn)
  if (length(shape)!= nn) shape <- rep_len(shape, nn)
  if (length(pstr1)!= nn) pstr1 <- rep_len(pstr1, nn)

  ans <- rep(NA_real_, nn)
  index1 <- (x == 1)
  if (log.arg) {
    ans[ index1] <- log(pstr1[ index1] + (1 - pstr1[ index1]) *
                      dzipf(x[ index1], N[ index1], shape[ index1]))
    ans[!index1] <- log1p(-pstr1[!index1]) +
                      dzipf(x[!index1], N[!index1], shape[!index1], log = TRUE)
  } else {
    ans[ index1] <-      pstr1[ index1] + (1 - pstr1[ index1]) *
                       dzipf(x[ index1], N[ index1], shape[ index1])
    ans[!index1] <- (1 - pstr1[!index1]) *
                       dzipf(x[!index1], N[!index1], shape[!index1])
  }


  deflat.limit <- deflat.limit.oizipf(N, shape)
  ans[pstr1 < deflat.limit] <- NaN
  ans[pstr1 > 1] <- NaN

  ans
}





poizipf <- function(q, N, shape, pstr1 = 0) {

  LLL <- max(length(q), length(N), length(shape), length(pstr1))
  if (length(q)     != LLL) q     <- rep_len(q,     LLL)
  if (length(N)     != LLL) N     <- rep_len(N,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oizipf(N, shape)

  ans <- pzipf(q, N, shape)  #, lower.tail = lower.tail, log.p = log.p
  ans <- ifelse(q < 1, 0, pstr1 + (1 - pstr1) * ans)

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN
  ans[s <= 0] <- NaN

  ans
}






qoizipf <- function(p, N, shape, pstr1 = 0) {

  if (!is.Numeric(p))
    stop("bad input for argument 'p'")
  if (!is.Numeric(N, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'N'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")

  nn <- max(length(p), length(N), length(s), length(pstr1))
  if (length(p)     != nn) p     <- rep_len(p,     nn)
  if (length(N)     != nn) N     <- rep_len(N,     nn)
  if (length(shape) != nn) shape <- rep_len(shape, nn)
  if (length(pstr1) != nn) pstr1 <- rep_len(pstr1, nn)


  ans    <- rep_len(NA_real_, nn)
  deflat.limit <- deflat.limit.oizipf(N, shape)

  dont.iterate <- 1 < p
  ans[p <= pstr1] <- 1
  pindex <- (pstr1 < p) & (deflat.limit <= pstr1) & !dont.iterate
  if (any(pindex))
  ans[pindex] <-
    qzipf((p[pindex] - pstr1[pindex]) / (1 - pstr1[pindex]),
          N = N[pindex], shape = shape[pindex])

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN

  ans[shape < 0] <- NaN
  ans[p < 0] <- NaN
  ans[1 < p] <- NaN

  ans
}



roizipf <- function(n, N, shape, pstr1 = 0) {
  qoizipf(runif(n), N, shape, pstr1 = pstr1)
}





 oizipf <-
  function(N = NULL, lpstr1 = "logit", lshape = "loge",
           type.fitted = c("mean", "shape", "pobs1", "pstr1", "onempstr1"),
           ishape = NULL,
           gpstr1 = ppoints(8),
           gshape = exp((-3:3) / 4), # grid for finding shape.init
           zero = NULL) {

  if (length(N) &&
     (!is.Numeric(N, positive = TRUE,
                 integer.valued = TRUE, length.arg = 1) ||
      N <= 1))
    stop("bad input for argument 'N'")
  enteredN <- length(N)

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
  blurb = c("One-inflated Zipf distribution f(y; pstr1, shape) = pstr1 + ",
            "(1 - pstr1) * y^(-shape) / sum((1:N)^(-shape)),",
            " 0 < shape, y = 1, 2,...,N",
            ifelse(enteredN, paste(" = ", N, sep = ""), ""),
            "\n\n",
            "Links:    ",
            namesof("pstr1", lpstr1, earg = epstr1 ), ", ",
            namesof("shape", lshape, earg = eshape ), "\n",
            "Mean:     pstr1 + (1 - pstr1) * ",
            "gharmonic(N, shape-1) / gharmonic(N, shape)"),

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



    NN <- .N
    if (!is.Numeric(NN, length.arg = 1,
                    positive = TRUE, integer.valued = TRUE))
      NN <- max(y)
    if (max(y) > NN)
      stop("maximum of the response is greater than argument 'N'")
    extra$N <- NN



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

      oizipf.Loglikfun <- function(pstr1, shape, y, x, w, extraargs) {
        sum(c(w) * doizipf(x = y, pstr1 = pstr1, N = extraargs$N,
                           s = shape, log = TRUE))
      }


      for (jay in 1:NOS) {  # For each response 'y_jay'... do:


        try.this <-
          grid.search2(gpstr1, gshape,
                       objfun = oizipf.Loglikfun,
                       y = y[, jay],  # x = x[TFvec, , drop = FALSE],
                       w = w[, jay],
                       extraargs = list(N = extra$N),
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
            .gshape = gshape,
            .type.fitted = type.fitted,
            .N = N ))),
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

    Meanfun <- function(shape, extra) {
      Mean <- shape
      Mean <- ( gharmonic2(extra$N, shape = shape - 1)
              / gharmonic2(extra$N, shape = shape))
      Mean[shape <= 0] <- NaN
      Mean
    }

    ans <-
      switch(type.fitted,
             "mean"      = pstr1 + (1 - pstr1) * Meanfun(shape, extra),
             "shape"     = shape,
             "pobs1"     = doizipf(1, N = extra$N, s = shape, pstr1 = pstr1),
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
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doizipf(x = y, pstr1 = pstr1, s = shape,
                                N = extra$N, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),
  vfamily = c("oizipf"),



  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    roizipf(nsim * length(shape), s = shape, pstr1 = pstr1,
            N = object@extra$N)
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                       earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lshape ,
                       earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape) &&
             all(is.finite(pstr1)) && all(pstr1 < 1)
    deflat.limit <- deflat.limit.oizipf(N = extra$N, s = shape)
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

    pmf1 <- dzipf(1, N = extra$N, shape = shape)
    onempmf1 <- 1 - pmf1  # dozeta(1, shape = shape, pstr1 = pstr1)
    pobs1 <- pstr1 + (1 - pstr1) * pmf1
    index1 <- as.matrix(y == 1)

    ghar0 <-  gharmonic2(extra$N, shape)
    ghar1 <-   gharmonic(extra$N, shape, deriv = 1)
    ghar2 <-   gharmonic(extra$N, shape, deriv = 2)

    dl.dpstr1 <- onempmf1 / pobs1
    dl.dpstr1[!index1] <- -1 / (1 - pstr1[!index1])

    dpmf1.dshape <- -ghar1 / ghar0^2

    d2pmf1.dshape2 <- (2 * ghar1^2 / ghar0 - ghar2) / ghar0^2

    dl.dshape <- (1 - pstr1) * dpmf1.dshape / pobs1  #
    dl.dshape[!index1] <- -log(y[!index1]) - ghar1[!index1] / ghar0[!index1]

    dpstr1.deta <- dtheta.deta(pstr1, .lpstr1 , earg = .epstr1 )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    myderiv <- c(w) * cbind(dl.dpstr1 * dpstr1.deta,
                            dl.dshape * dshape.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape ))),
  weight = eval(substitute(expression({

    LHS <- ((1 - pstr1) / pobs1) * dpmf1.dshape^2 - d2pmf1.dshape2
    RHS <- (ghar2 - ghar1^2 / ghar0) / ghar0
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
}  # oizipf





dotzeta <- function(x, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (log.arg) {
    ans <- dzeta(x, shape, log = log.arg) - log1p(-dzeta(1, shape))
    ans[x == 1] <- log(0)
  } else {
    ans <- dzeta(x, shape) / (1 - dzeta(1, shape))
    ans[x == 1] <- 0
  }
  ans[shape < 0] <- NaN
  ans
}  # dotzeta



potzeta <- function(q, shape, log.p = FALSE) {
  if (log.p) log(pzeta(q, shape) - dzeta(1, shape)) -
      log1p(-dzeta(1, shape)) else
    (pzeta(q, shape) - dzeta(1, shape)) / (1 - dzeta(1, shape))
}



 qotzeta <- function(p, shape) {
  ans <- qzeta((1 - dzeta(1, shape)) * p + dzeta(1, shape), shape = shape)

  ans[p == 1] <- Inf
  ans[p < 0 | 1 < p] <- NaN

  ans[shape < 0] <- NaN
  ans
}  # qotzeta



rotzeta <- function(n, shape) {
  qotzeta(runif(n), shape)
}





 otzeta <-
    function(lshape = "loge",
             ishape = NULL,
             gshape = exp((-4:3)/4),
             zero = NULL) {

  if (length(ishape) && !is.Numeric(ishape, positive = TRUE))
    stop("argument 'ishape' must be > 0")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  new("vglmff",
  blurb = c("One-truncated Zeta distribution ",
            "f(y) = 1/(y^(shape+1) * (zeta(shape+1) - 1 - 1/2^(shape+1)))",
            " 0<shape, y = 2, 3,...\n\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape)),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = FALSE,  # NR == FS
         multipleResponses = TRUE,
         parameters.names = "shape",
         zero = .zero ,
         lshape = .lshape )
  }, list( .lshape = lshape,
           .zero = zero ))),
  initialize = eval(substitute(expression({

   temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              Is.positive.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (any(y <= 1))
      stop("no 1s in the response allowed!")


    ncoly <- ncol(y)
    mynames1 <- param.names("shape", ncoly)
    predictors.names <-
      namesof(mynames1, .lshape , earg = .eshape , tag = FALSE)

    M1 <- 1
    M <- M1 * ncoly


    if (!length(etastart)) {
      otzetaff.Loglikfun <- function(shape, y, x, w, extraargs) {
        sum(c(w) * dotzeta(x = y, shape, log = TRUE))
      }


      gshape <- .gshape
      if (!length( .ishape )) {
        shape.init <- matrix(NA_real_, n, M, byrow = TRUE)
        for (jay in 1:ncoly) {
          shape.init[, jay] <- grid.search(gshape, objfun = otzetaff.Loglikfun,
                                           y = y[, jay], x = x, w = w[, jay])
        }
      } else {
        shape.init <- matrix( .ishape , n, M, byrow = TRUE)
      }
      etastart <- theta2eta(shape.init, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .eshape = eshape,
            .ishape = ishape, .gshape = gshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    ans <- pp <- eta2theta(eta, .lshape , earg = .eshape )
    ans[pp > 1] <- zeta(pp[pp > 1]) / zeta(pp[pp > 1] + 1)
    ans[pp <= 1] <- NA
    pmf.1 <- dzeta(1, pp)
    (ans - pmf.1) / (1 - pmf.1)
  }, list( .lshape = lshape, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$link <- rep_len( .lshape , ncoly)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (jay in 1:ncoly) {
      misc$earg[[jay]] <- .eshape
    }

  }), list( .lshape = lshape, .eshape = eshape ))),
  loglikelihood = eval(substitute( function(mu, y, w, residuals = FALSE,
             eta, extra = NULL, summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dotzeta(x = y, shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("otzeta"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape)
    okay1
  }, list( .lshape = lshape, .eshape = eshape ))),
  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    BBBB  <- zeta(shape + 1) - 1
    fred1 <- zeta(shape + 1, deriv = 1)
    dl.dshape <- -log(y) - fred1 / BBBB
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = expression({
    ned2l.dshape2 <- (zeta(shape + 1, deriv = 2) - fred1^2 / BBBB) / BBBB
    wz <- ned2l.dshape2 * dshape.deta^2
    c(w) * wz
  }))
}







ddiffzeta <- function(x, shape, start = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(shape), length(x), length(start))
  if (length(x)     != LLL) x     <- rep_len(x,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(start) != LLL) start <- rep_len(start, LLL)

  ox <- !is.finite(x)
  zero <- ox | round(x) != x | x < start
  ans <- rep_len(if (log.arg) log(0) else 0, LLL)
  if (any(!zero)) {
    ans[!zero] <- (start[!zero] /      x[!zero]) ^(shape[!zero]) -
                  (start[!zero] / (1 + x[!zero]))^(shape[!zero])
    if (log.arg)
      ans[!zero] <- log(ans[!zero])
  }
  if (any(ox))
    ans[ox] <- if (log.arg) log(0) else 0
  ans[shape <= 0] <- NaN
  ans[start != round(start) | start < 1] <- NaN
  ans
}



 pdiffzeta <- function(q, shape, start = 1, lower.tail = TRUE) {

  LLL <- max(length(shape), length(q), length(start))
  if (length(q)     != LLL) q     <- rep_len(q,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(start) != LLL) start <- rep_len(start, LLL)

  if (lower.tail) {
    ans <- 1 - (start / floor(1 + q))^shape
  } else {
    ans <-     (start / floor(1 + q))^shape
  }
  ans[q < start] <- if (lower.tail) 0 else 1
  ans[shape <= 0] <- NaN
  ans[start != round(start) | start < 1] <- NaN
  ans
}  # pdiffzeta




 qdiffzeta <- function(p, shape, start = 1) {

  LLL <- max(length(p), length(shape), length(start))
  if (length(p)     != LLL) p     <- rep_len(p,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(start) != LLL) start <- rep_len(start, LLL)

  lo <- rep_len(start, LLL)
  approx.ans <- lo  # True at lhs
  hi <- 2 * lo + 10
  dont.iterate <- p == 1 | shape <= 0 | start != round(start) | start < 1
  done <- p <= pdiffzeta(hi, shape, start = start) | dont.iterate
  max.iter <- 100
  iter <- 0
  while (!all(done) && iter < max.iter) {
    hi.save <- hi[!done]
    hi[!done] <- 2 * lo[!done] + 10
    lo[!done] <- hi.save
    done[!done] <- is.infinite(hi[!done]) |
                   (p[!done] <= pdiffzeta(hi[!done], shape[!done], start[!done]))
    iter <- iter + 1
  }

  foo <- function(q, shape, start, p)
    pdiffzeta(q, shape, start) - p

  lhs <- (p <= ddiffzeta(start, shape, start = start)) | dont.iterate

  approx.ans[!lhs] <- bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                                      shape = shape[!lhs],
                                      start = start[!lhs], p = p[!lhs])
  faa <- floor(approx.ans)
  ans <- ifelse(pdiffzeta(faa  , shape, start = start) <  p &
                p <= pdiffzeta(faa+1, shape, start = start), faa+1, faa)

  ans[p == 1] <- Inf
  ans[shape <= 0] <- NaN
  ans[start != round(start) | start < 1] <- NaN

  ans
}  # qdiffzeta




rdiffzeta <- function(n, shape, start = 1) {
  rr <- runif(n)
  qdiffzeta(rr, shape, start = start)
}



 diffzeta <- function(start = 1, lshape = "loge", ishape = NULL) {


  if (!is.Numeric(start, positive = TRUE,
                  integer.valued = TRUE, length.arg = 1))
    stop("bad input for argument 'start'")
  enteredstart <- length(start)

  if (length(ishape) && !is.Numeric(ishape, positive = TRUE))
    stop("argument 'ishape' must be > 0")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  new("vglmff",
  blurb = c("Difference in 2 Zipf distributions ",
            "f(y;s) = y^(-shape) / sum((1:start)^(-shape)), ",
            "shape > 0, start, start+1,...",
            ifelse(enteredstart, paste("start = ", start, sep = ""), ""),
            "\n\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape),
            "\n\n",
            "Mean:    gharmonic(start, shape-1) / gharmonic(start, shape)"),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         start = .start ,
         parameters.names = "shape")
  }, list( .start = start ))),
  initialize = eval(substitute(expression({
  start <- .start
  temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              Is.positive.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (any(y < start))
      stop("some response values less than 'start'")


    predictors.names <- namesof("shape", .lshape , earg = .eshape , tag = FALSE)

    extra$start <- start
    if (!length(etastart)) {
      llfun <- function(shape, y, start, w) {
        sum(c(w) * ddiffzeta(x = y, start = extra$start, shape = shape, log = TRUE))
      }
      shape.init <- if (length( .ishape )) .ishape else
        getInitVals(gvals = seq(0.1, 3.0, length.out = 19),
                    llfun = llfun,
                    y = y, start = extra$start, w = w)
      shape.init <- rep_len(shape.init, length(y))
      if ( .lshape == "loglog") shape.init[shape.init <= 1] <- 1.2
      etastart <- theta2eta(shape.init, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape,
            .eshape = eshape, .ishape = ishape, .start = start ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    aa <- extra$start
    if (length(aa) != 1 || aa < 1 || round(aa) != aa)
      stop("the 'start' variable must be of unit length")
    if (aa == 1)
      return(zeta(shape))
    mymat <- matrix(1:aa, NROW(eta), aa, byrow = TRUE)
    temp1 <- rowSums(1 / mymat^shape)
    (aa^shape) * (zeta(shape) - temp1 + 1 / aa^(shape-1))
  }, list( .lshape = lshape, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$expected <- FALSE
    misc$link <-    c(shape = .lshape )
    misc$earg <- list(shape = .eshape )
    misc$start <- extra$start
  }), list( .lshape = lshape, .eshape = eshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * ddiffzeta(x = y, start = extra$start,
                                  shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("diffzeta"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape)
    okay1
  }, list( .lshape = lshape, .eshape = eshape ))),






  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    temp1 <- extra$start /  y
    temp2 <- extra$start / (y+1)
    AA <- temp1^shape - temp2^shape
    Aprime <- log(temp1) * temp1^shape -
              log(temp2) * temp2^shape

    dl.dshape <- Aprime / AA
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = expression({
    ned2l.dshape <- (Aprime / AA)^2  # Not quite FS. Half FS.
    wz <- c(w) * ned2l.dshape * dshape.deta^2
    wz
  }))
}
















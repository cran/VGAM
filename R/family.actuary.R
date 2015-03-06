# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.













dgumbelII <- function(x, scale = 1, shape, log = FALSE) {


  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(scale))
  if (length(x)       != LLL) x       <- rep(x,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)


  ans <- x
  index0 <- (x < 0) & is.finite(x) & !is.na(x)

  ans[!index0] <- log(shape[!index0] / scale[!index0]) +
            (shape[!index0] + 1) * log(scale[!index0] / x[!index0]) -
             (x[!index0] / scale[!index0])^(-shape[!index0])
  ans[index0] <- log(0)
  ans[x == Inf] <- log(0)

  if (log.arg) {
  } else {
    ans <- exp(ans)
    ans[index0] <- 0
    ans[x == Inf] <- 0
  }
  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}



pgumbelII <- function(q, scale = 1, shape,
                      lower.tail = TRUE, log.p = FALSE) {

  # 20150121 KaiH
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  # 20150121 KaiH
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  LLL <- max(length(q), length(shape), length(scale)) 
  if (length(q)       != LLL) q       <- rep(q,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  
  # 20150121 KaiH 
  if (lower.tail) {
    if (log.p) { 
      ans <- -(q / scale)^(-shape)
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else { 
      ans <- exp(-(q / scale)^(-shape))
      ans[q <= 0] <- 0 
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- log(-expm1(-(q / scale)^(-shape)))
      ans[q <= 0] <- 0 
      ans[q == Inf] <- -Inf
    } else { 
      ans <- -expm1(-(q / scale)^(-shape))
      ans[q <= 0] <- 1 
      ans[q == Inf] <- 0
    }
  } 
  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}










qgumbelII <- function(p, scale = 1, shape,
                      lower.tail = TRUE, log.p = FALSE) {



  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")



  LLL <- max(length(p), length(shape), length(scale))
  if (length(p)       != LLL) p       <- rep(p,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)


  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- scale * (-ln.p)^(-1 / shape)
      ans[ln.p > 0] <- NaN
    } else { # Default
      ans <- scale * (-log(p))^(-1 / shape)
      ans[p < 0] <- NaN
      ans[p == 0] <- 0
      ans[p == 1] <- Inf
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- scale * (-log(-expm1(ln.p)))^(-1 / shape)
      ans[ln.p > 0] <- NaN
    } else {
      ans <- scale * (-log1p(-p))^(-1 / shape)
      ans[p < 0] <- NaN
      ans[p == 0] <- Inf
      ans[p == 1] <- 0
      ans[p > 1] <- NaN
    }
  }

  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}


rgumbelII <- function(n, scale = 1, shape) {
  qgumbelII(runif(n), shape = shape, scale = scale)
}









 gumbelII <-
  function(lscale = "loge", lshape = "loge",
           iscale = NULL,   ishape = NULL,
           probs.y = c(0.2, 0.5, 0.8),
           perc.out = NULL,  # 50,
           imethod = 1, zero = -1, nowarning = FALSE) {




  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE))
    stop("bad input for argument 'zero'")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")
  if (!is.Numeric(probs.y, positive  = TRUE) ||
      length(probs.y) < 2 ||
      max(probs.y) >= 1)
    stop("bad input for argument 'probs.y'")
  if (length(perc.out))
    if (!is.Numeric(perc.out, positive  = TRUE) ||
        max(probs.y) >= 100)
    stop("bad input for argument 'perc.out'")


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' values must be positive")
  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
      stop("argument 'iscale' values must be positive")


  new("vglmff",
  blurb = c("Gumbel Type II distribution\n\n",
            "Links:    ",
            namesof("scale", lscale, escale), ", ",
            namesof("shape", lshape, eshape), "\n",
            "Mean:     scale^(1/shape) * gamma(1 - 1 / shape)\n",
            "Variance: scale^(2/shape) * (gamma(1 - 2/shape) - ",
                      "gamma(1 + 1/shape)^2)"),
 constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         perc.out = .perc.out ,
         zero = .zero )
  }, list( .zero = zero,
           .perc.out = perc.out
         ))),

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
    M1 <- 2
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1 <- paste("scale", if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("shape", if (ncoly > 1) 1:ncoly else "", sep = "")


    predictors.names <-
        c(namesof(mynames1, .lscale , .escale , tag = FALSE),
          namesof(mynames2, .lshape , .eshape , tag = FALSE))[
          interleave.VGAM(M, M = M1)]


    Shape.init <- matrix(if (length( .ishape )) .ishape else 0 + NA,
                         n, ncoly, byrow = TRUE)
    Scale.init <- matrix(if (length( .iscale )) .iscale else 0 + NA,
                         n, ncoly, byrow = TRUE)

    if (!length(etastart)) {
      if (!length( .ishape ) ||
          !length( .iscale )) {
        for (ilocal in 1:ncoly) {

          anyc <- FALSE # extra$leftcensored | extra$rightcensored
          i11 <- if ( .imethod == 1) anyc else FALSE # can be all data
          probs.y <- .probs.y
          xvec <- log(-log(probs.y))
          fit0 <- lsfit(y  = xvec,
                        x  = log(quantile(y[!i11, ilocal],
                                          probs = probs.y )))


          if (!is.Numeric(Shape.init[, ilocal]))
            Shape.init[, ilocal] <- -fit0$coef["X"]
          if (!is.Numeric(Scale.init[, ilocal]))
            Scale.init[, ilocal] <-
              exp(fit0$coef["Intercept"] / Shape.init[, ilocal])
        }  # ilocal

        etastart <-
          cbind(theta2eta(Scale.init, .lscale , .escale ),
                theta2eta(Shape.init, .lshape , .eshape ))[,
                interleave.VGAM(M, M = M1)]
      }
    }
  }), list(
            .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .iscale = iscale, .ishape = ishape,
            .probs.y = probs.y,
            .imethod = imethod ) )),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , .escale )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , .eshape )
    Shape <- as.matrix(Shape)

    if (length( .perc.out ) > 1 && ncol(Shape) > 1)
      stop("argument 'perc.out' should be of length one since ",
           "there are multiple responses")

    if (!length( .perc.out )) {
      return(Scale * gamma(1 - 1 / Shape))
    }

    ans <- if (length( .perc.out ) > 1) {
      qgumbelII(p = matrix( .perc.out / 100, length(Shape),
                           length( .perc.out ), byrow = TRUE),
                shape = Shape, scale = Scale)
    } else {
      qgumbelII(p = .perc.out / 100, shape = Shape, scale = Scale)
    }
    colnames(ans) <- paste(as.character( .perc.out ), "%", sep = "")
    ans
  }, list(
           .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape,
           .perc.out = perc.out ) )),
  last = eval(substitute(expression({


    M1 <- extra$M1
    misc$link <-
      c(rep( .lscale , length = ncoly),
        rep( .lshape , length = ncoly))[interleave.VGAM(M, M = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .escale
      misc$earg[[M1*ii  ]] <- .eshape
    }

    misc$M1 <- M1
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$perc.out <- .perc.out
    misc$true.mu <- FALSE # @fitted is not a true mu


  }), list(
            .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .perc.out = perc.out,
            .imethod = imethod ) )),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    Scale <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , .escale )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dgumbelII(x = y, shape = Shape,
                                  scale = Scale, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape
         ) )),
  vfamily = c("gumbelII"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    Scale <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , .escale )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , .eshape )
    rgumbelII(nsim * length(Scale), shape = Shape, scale = Scale)
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape
         ) )),



  deriv = eval(substitute(expression({
    M1 <- 2
    Scale <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , .escale )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , .eshape )

    dl.dshape <- 1 / Shape + log(Scale / y) -
                 log(Scale / y) * (Scale / y)^Shape
    dl.dscale <- Shape / Scale - (Shape / y) * (Scale / y)^(Shape - 1)


    dscale.deta <- dtheta.deta(Scale, .lscale , .escale )
    dshape.deta <- dtheta.deta(Shape, .lshape , .eshape )

    myderiv <- c(w) * cbind(dl.dscale, dl.dshape) *
                      cbind(dscale.deta, dshape.deta)
    myderiv[, interleave.VGAM(M, M = M1)]
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape
          ) )),
  weight = eval(substitute(expression({
    EulerM <- -digamma(1.0)


    ned2l.dshape2 <- (1 + trigamma(2) + digamma(2)^2) / Shape^2
    ned2l.dscale2 <-  (Shape / Scale)^2
    ned2l.dshapescale <- digamma(2) / Scale


    wz <- array(c(c(w) * ned2l.dscale2 * dscale.deta^2,
                  c(w) * ned2l.dshape2 * dshape.deta^2,
                  c(w) * ned2l.dshapescale * dscale.deta * dshape.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }), list( .lscale = lscale, .lshape = lshape ))))
}





dmbeard <- function(x, shape, scale = 1, rho, epsilon, log = FALSE) {


  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(scale),
             length(rho), length(epsilon))
  if (length(x)       != LLL) x       <- rep(x,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  if (length(rho)     != LLL) rho     <- rep(rho,     length.out = LLL)
  if (length(epsilon) != LLL) epsilon <- rep(epsilon, length.out = LLL)


  index0 <- (x < 0)

  ans <- log(epsilon * exp(-x * scale) + shape) +
            (-epsilon * x -
            ((rho * epsilon - 1) / (rho * scale)) *
            (log1p(rho * shape) -
             log(exp(-x * scale) + rho * shape) - scale * x)) - 
            log(exp(-x * scale) + shape * rho)

  ans[index0] <- log(0)
  ans[x == Inf] <- log(0)

  if (log.arg) {
  } else {
    ans <- exp(ans)
    ans[index0] <- 0
    ans[x == Inf] <- 0
  }
  ans[shape <= 0 | scale <= 0 | rho <= 0 | epsilon <= 0] <- NaN
  ans
}


pmbeard <- function(q, shape, scale = 1, rho, epsilon) {

  LLL <- max(length(q), length(shape), length(scale),
             length(rho), length(epsilon))
  if (length(q)       != LLL) q       <- rep(q,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  if (length(rho)     != LLL) rho     <- rep(rho,     length.out = LLL)
  if (length(epsilon) != LLL) epsilon <- rep(epsilon, length.out = LLL)


  ans <- -expm1(-epsilon * q -
               ((rho * epsilon - 1) / (rho * scale)) *
               (log1p(rho * shape) -
                log(exp(-scale * q) + rho * shape) - scale * q))
  ans[(q <= 0)] <- 0
  ans[shape <= 0 | scale <= 0 | rho <= 0 | epsilon <= 0] <- NaN
  ans[q == Inf] <- 1
  ans
}







dmperks <- function(x, scale = 1, shape, epsilon, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(scale), length(epsilon))
  if (length(x)       != LLL) x       <- rep(x,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  if (length(epsilon) != LLL) epsilon <- rep(epsilon, length.out = LLL)


  index0 <- (x < 0)
  ans <- log(epsilon * exp(-x * scale) + shape) +
            (-epsilon * x -
            ((epsilon - 1) / scale) *
            (log1p(shape) -
             log(shape + exp(-x * scale)) -x * scale)) - 
            log(exp(-x * scale) + shape)

  ans[index0] <- log(0)
  ans[x == Inf] <- log(0)
  if (log.arg) {
  } else {
    ans <- exp(ans)
    ans[index0] <- 0
    ans[x == Inf] <- 0
  }
  ans[shape <= 0 | scale <= 0 | epsilon <= 0] <- NaN
  ans
}



pmperks <- function(q, scale = 1, shape, epsilon) {

  LLL <- max(length(q), length(shape), length(scale))
  if (length(q)       != LLL) q       <- rep(q,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)


  ans <- -expm1(-epsilon * q -
               ((epsilon - 1) / scale) *
               (log1p(shape) -
                log(shape + exp(-q * scale)) - q * scale))

  ans[(q <= 0)] <- 0
  ans[shape <= 0 | scale <= 0] <- NaN
  ans[q == Inf] <- 1
  ans
}












dbeard <- function(x, shape, scale = 1, rho, log = FALSE) {

 warning("does not integrate to unity")

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(scale), length(rho))
  if (length(x)       != LLL) x       <- rep(x,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  if (length(rho)     != LLL) rho     <- rep(rho,     length.out = LLL)

  index0 <- (x < 0)
    ans <- log(shape) - x * scale * (rho^(-1 / scale)) +
           log(rho) + log(scale) +
           (rho^(-1 / scale)) * log1p(shape * rho) -
           (1 + rho^(-1 / scale)) *
           log(shape * rho + exp(-x * scale))
    ans[index0] <- log(0)
    ans[x == Inf] <- log(0)


  if (log.arg) {
  } else {
    ans <- exp(ans)
    ans[index0] <- 0
    ans[x == Inf] <- 0
  }
  ans[shape <= 0 | scale <= 0 | rho <= 0] <- NaN
  ans
}






dbeard <- function(x, shape, scale = 1, rho, log = FALSE) {
  alpha <- shape
  beta  <- scale

 warning("does not integrate to unity")

  ret <- ifelse(x <= 0 | beta <= 0, NaN,
                exp(alpha+beta*x)*(1+exp(alpha+rho))**(exp(-rho/beta))/
                (1+exp(alpha+rho+beta*x))**(1+exp(-rho/beta)))
  ret
}



qbeard <- function(x, u = 0.5, alpha = 1, beta = 1,rho = 1) {
  ret <-
    ifelse(x <= 0 | u <= 0 | u >= 1 |
           length(x) != length(u) | beta <= 0,
           NaN,
           (1/beta) * (log((u**(-beta*exp(rho))) *
           (1+exp(alpha+rho+beta*x))-1)-alpha-rho)-x)

  return(ret)
}










dperks <- function(x, scale = 1, shape, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(scale))
  if (length(x)     != LLL) x     <- rep(x,     length.out = LLL)
  if (length(shape) != LLL) shape <- rep(shape, length.out = LLL)
  if (length(scale) != LLL) scale <- rep(scale, length.out = LLL)

  index0 <- (x < 0)
    ans <- log(shape) - x +
           log1p(shape) / scale -
           (1 + 1 / scale) * log(shape + exp(-x * scale))
    ans[index0] <- log(0)
    ans[x == Inf] <- log(0)

  if (log.arg) {
  } else {
    ans <- exp(ans)
    ans[index0] <- 0
    ans[x == Inf] <- 0
  }
  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}



pperks <- function(q, scale = 1, shape,
                   lower.tail = TRUE, log.p = FALSE) {


  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  LLL <- max(length(q), length(shape), length(scale))
  if (length(q)       != LLL) q       <- rep(q,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)

  logS <- -q + (log1p(shape) -
          log(shape + exp(-q * scale))) / scale


  if (lower.tail) {
    if (log.p) {
      ans <- log(-expm1(logS))
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- -expm1(logS)
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- logS
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      ans <- exp(logS)
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  }

  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}


qperks <- function(p, scale = 1, shape, lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  LLL <- max(length(p), length(shape), length(scale))
  if (length(p)       != LLL) p       <- rep(p,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)


  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      tmp <- scale * log(-expm1(ln.p))
      onemFb <- exp(tmp)
      ans <- (log1p(shape - onemFb) - log(shape) - tmp) / scale
      ans[ln.p > 0] <- NaN
    } else {
      tmp <- scale * log1p(-p)
      onemFb <- exp(tmp)
      ans <- (log1p(shape - onemFb) - log(shape) - tmp) / scale
      ans[p < 0] <- NaN
      ans[p == 0] <- 0
      ans[p == 1] <- Inf
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      tmp <- scale * ln.p
      onemFb <- exp(tmp)
      ans <- (log1p(shape - onemFb) - log(shape) - tmp) / scale
      ans[ln.p > 0] <- NaN
    } else {
      tmp <- scale * log(p)
      onemFb <- exp(tmp)
      ans <- (log1p(shape - onemFb) - log(shape) - tmp) / scale
      ans[p < 0] <- NaN
      ans[p == 0] <- Inf
      ans[p == 1] <- 0
      ans[p > 1] <- NaN
    }
  }

  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}



rperks <- function(n, scale = 1, shape) {
  qperks(runif(n), scale = scale, shape = shape)
}





perks.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 perks <-
  function(lscale = "loge",    lshape = "loge",
           iscale = NULL,      ishape = NULL,
           gscale = exp(-5:5), gshape = exp(-5:5),
           nsimEIM = 500,
           oim.mean = FALSE,
           zero = NULL, nowarning = FALSE) {



  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")


  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE))
    stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
    warning("argument 'nsimEIM' should be an integer ",
            "greater than 50, say")


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' values must be positive")
  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
      stop("argument 'iscale' values must be positive")




    if (!is.logical(oim.mean) || length(oim.mean) != 1)
      stop("bad input for argument 'oim.mean'")



  new("vglmff",
  blurb = c("Perks' distribution\n\n",
            "Links:    ",
            namesof("scale", lscale, escale), ", ",
            namesof("shape", lshape, eshape), "\n",
            "Median:   qperks(p = 0.5, scale = scale, shape = shape)"),

  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         nsimEIM = .nsimEIM,
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM ))),
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
    M1 <- 2
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1 <- paste("scale", if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("shape", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lscale , .escale , tag = FALSE),
          namesof(mynames2, .lshape , .eshape , tag = FALSE))[
          interleave.VGAM(M, M = M1)]



    if (!length(etastart)) {

      matH <- matrix(if (length( .ishape )) .ishape else 0 + NA,
                     n, ncoly, byrow = TRUE)
      matC <- matrix(if (length( .iscale )) .iscale else 0 + NA,
                     n, ncoly, byrow = TRUE)

      shape.grid <- .gshape
      scale.grid <- .gscale

      for (spp. in 1:ncoly) {
        yvec <- y[, spp.]
        wvec <- w[, spp.]

        perks.Loglikfun <- function(scaleval, y, x, w, extraargs) {
          ans <-
          sum(c(w) * dperks(x = y, shape = extraargs$Shape,
                            scale = scaleval, log = TRUE))
          ans
        }


        mymat <- matrix(-1, length(shape.grid), 2)
        for (jlocal in 1:length(shape.grid)) {
          mymat[jlocal, ] <-
            grid.search(scale.grid,
                        objfun = perks.Loglikfun,
                        y = yvec, x = x, w = wvec,
                        ret.objfun = TRUE,
                        extraargs = list(Shape = shape.grid[jlocal]))
        }
        index.shape <- which(mymat[, 2] == max(mymat[, 2]))[1]

        if (!length( .ishape ))
          matH[, spp.] <- shape.grid[index.shape]
        if (!length( .iscale ))
          matC[, spp.] <- mymat[index.shape, 1]
      }  # spp.

      etastart <-
          cbind(theta2eta(matC, .lscale , .escale ),
                theta2eta(matH, .lshape , .eshape ))[,
                interleave.VGAM(M, M = M1)]
    }  # End of !length(etastart)
  }), list( .lscale = lscale, .lshape = lshape,
            .eshape = eshape, .escale = escale,
            .gshape = gshape, .gscale = gscale,
            .ishape = ishape, .iscale = iscale
            ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , .escale )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , .eshape )

    qperks(p = 0.5, shape = Shape, scale = Scale)
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape ))),
  last = eval(substitute(expression({

    misc$link <-
      c(rep( .lscale , length = ncoly),
        rep( .lshape , length = ncoly))[interleave.VGAM(M, M = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .escale
      misc$earg[[M1*ii  ]] <- .eshape
    }


    misc$M1 <- M1
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    Scale <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , .escale )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dperks(x = y, shape = Shape,
                               scale = Scale, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape ))),
  vfamily = c("perks"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    Scale <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , .escale )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , .eshape )
    rperks(nsim * length(Scale), shape = Shape, scale = Scale)
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape ))),






 
  deriv = eval(substitute(expression({
    M1 <- 2
    scale <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                       .lscale , .escale )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                       .lshape , .eshape )


    temp2 <- exp(y * scale)
    temp3 <- 1 + shape * temp2
    dl.dshape <- 1 / shape + 1 / (scale * (1 + shape)) -
                 (1 + 1 / scale) * temp2 / temp3
    dl.dscale <- y - log1p(shape) / scale^2 +
                 log1p(shape * temp2) / scale^2 -
                 (1 + 1 / scale) * shape * y * temp2 / temp3

    dshape.deta <- dtheta.deta(shape, .lshape , .eshape )
    dscale.deta <- dtheta.deta(scale, .lscale , .escale )

    dthetas.detas <- cbind(dscale.deta, dshape.deta)
    myderiv <- c(w) * cbind(dl.dscale, dl.dshape) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = M1)]
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape ))),


  weight = eval(substitute(expression({

    NOS <- M / M1
    dThetas.detas <- dthetas.detas[, interleave.VGAM(M, M = M1)]

    wz <- matrix(0.0, n, M + M - 1)  # wz is 'tridiagonal' 

    ind1 <- iam(NA, NA, M = M1, both = TRUE, diag = TRUE)


    for (spp. in 1:NOS) {
      run.varcov <- 0
      Scale <- scale[, spp.]
      Shape <- shape[, spp.]




      if (FALSE && intercept.only && .oim.mean ) {

 stop("this is wrong")
      temp8 <- (1 + Shape * exp(Scale * y[, spp.]))^2
      nd2l.dadb <- 2 * y[, spp.] * exp(Scale * y[, spp.]) / temp8

      nd2l.dada <- 1 / Shape^2 + 1 / (1 + Shape)^2 -
        2 * exp(2 * Scale * y[, spp.]) / temp8

      nd2l.dbdb <- 2 * Shape * y[, spp.]^2 * exp(Scale * y[, spp.]) / temp8


      ave.oim11 <- weighted.mean(nd2l.dada, w[, spp.])
      ave.oim12 <- weighted.mean(nd2l.dadb, w[, spp.])
      ave.oim22 <- weighted.mean(nd2l.dbdb, w[, spp.])
      run.varcov <- cbind(ave.oim11, ave.oim22, ave.oim12)
    } else {

      for (ii in 1:( .nsimEIM )) {
        ysim <- rperks(n = n, shape = Shape, scale = Scale)
if (ii < 3) {
}

        temp2 <- exp(ysim * Scale)
        temp3 <- 1 + Shape * temp2
        dl.dshape <- 1 / Shape + 1 / (Scale * (1 + Shape)) -
                     (1 + 1 / Scale) * temp2 / temp3
        dl.dscale <- ysim - log1p(Shape) / Scale^2 +
                     log1p(Shape * temp2) / Scale^2 -
                     (1 + 1 / Scale) * Shape * ysim * temp2 / temp3


        temp7 <- cbind(dl.dscale, dl.dshape)
if (ii < 3) {
}
        run.varcov <- run.varcov +
                      temp7[, ind1$row.index] *
                      temp7[, ind1$col.index]
      }
      run.varcov <- cbind(run.varcov / .nsimEIM )

    }



      wz1 <- if (intercept.only)
          matrix(colMeans(run.varcov),
                 nrow = n, ncol = ncol(run.varcov), byrow = TRUE) else
          run.varcov

      wz1 <- wz1 * dThetas.detas[, M1 * (spp. - 1) + ind1$row] *
                   dThetas.detas[, M1 * (spp. - 1) + ind1$col]


      for (jay in 1:M1)
        for (kay in jay:M1) {
          cptr <- iam((spp. - 1) * M1 + jay,
                      (spp. - 1) * M1 + kay,
                      M = M)
          wz[, cptr] <- wz1[, iam(jay, kay, M = M1)]
        }
    }  # End of for (spp.) loop



    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / M1)
  }), list( .lscale = lscale,
            .escale = escale,
            .nsimEIM = nsimEIM, .oim.mean = oim.mean ))))
}  # perks()








dmakeham <- function(x, scale = 1, shape, epsilon = 0, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(scale), length(epsilon))
  if (length(x)       != LLL) x       <- rep(x,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  if (length(epsilon) != LLL) epsilon <- rep(epsilon, length.out = LLL)

  index0 <- (x < 0)
  ans <- log(epsilon * exp(-x * scale) + shape) +
         x * (scale - epsilon) -
         (shape / scale) * expm1(x * scale)
  ans[index0] <- log(0)
  ans[x == Inf] <- log(0)
  if (log.arg) {
  } else {
    ans <- exp(ans)
    ans[index0] <- 0
    ans[x == Inf] <- 0
  }
  ans[shape <= 0 | scale <= 0 | epsilon < 0] <- NaN
  ans
}



pmakeham <- function(q, scale = 1, shape, epsilon = 0,
                     lower.tail = TRUE, log.p = FALSE) {


  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  LLL <- max(length(q), length(shape), length(scale), length(epsilon))
  if (length(q)       != LLL) q       <- rep(q,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  if (length(epsilon) != LLL) epsilon <- rep(epsilon, length.out = LLL)

  if (lower.tail) {
    if (log.p) {
      ans <- log(-expm1(-q * epsilon - (shape / scale) * expm1(scale * q)))
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- -expm1(-q * epsilon - (shape / scale) * expm1(scale * q))
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- -q * epsilon - (shape / scale) * expm1(scale * q)
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      ans <- exp(-q * epsilon - (shape / scale) * expm1(scale * q))
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  }

  ans[shape <= 0 | scale <= 0 | epsilon < 0] <- NaN
  ans
}



qmakeham <- function(p, scale = 1, shape, epsilon = 0,
                     lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  LLL <- max(length(p), length(shape), length(scale), length(epsilon))
  if (length(p)       != LLL) p       <- rep(p,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  if (length(epsilon) != LLL) epsilon <- rep(epsilon, length.out = LLL)


  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- shape / (scale * epsilon) - log(-expm1(ln.p)) / epsilon -
        lambertW((shape / epsilon) * exp(shape / epsilon) *
                   exp(log(-expm1(ln.p)) * (-scale / epsilon))) / scale
      ans[ln.p == 0] <- Inf
      ans[ln.p > 0] <- NaN
    } else {
      ans <- shape / (scale * epsilon) - log1p(-p) / epsilon -
        lambertW((shape / epsilon) * exp(shape / epsilon) *
                   exp( (-scale / epsilon) * log1p(-p) )) / scale
      ans[p < 0] <- NaN
      ans[p == 0] <- 0
      ans[p == 1] <- Inf
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <-  shape / (scale * epsilon) - ln.p / epsilon -
        lambertW((shape / epsilon) * exp(shape / epsilon) *
                  exp(ln.p * (-scale / epsilon))) / scale
      ans[ln.p == -Inf] <- Inf
      ans[ln.p > 0] <- NaN
    } else {
      ans <- shape / (scale * epsilon) - log(p) / epsilon -
        lambertW((shape / epsilon) * exp(shape / epsilon) *
                  p^(-scale / epsilon)) / scale
      ans[p < 0] <- NaN
      ans[p == 0] <- Inf
      ans[p == 1] <- 0
      ans[p > 1] <- NaN
    }
  }

  ans[epsilon == 0] <-
    qgompertz(p     =     p[epsilon == 0],
              shape = shape[epsilon == 0],
              scale = scale[epsilon == 0],
              lower.tail = lower.tail,
              log.p = log.p)

  ans[shape <= 0 | scale <= 0 | epsilon < 0] <- NaN
  ans
}



rmakeham <- function(n, scale = 1, shape, epsilon = 0) {
  qmakeham(runif(n), scale = scale, shape = shape, epsilon = epsilon)
}




makeham.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 makeham <-
  function(lscale = "loge", lshape = "loge", lepsilon = "loge",
           iscale = NULL,   ishape = NULL,   iepsilon = NULL,  # 0.3,
           gscale = exp(-5:5),
           gshape = exp(-5:5),
           gepsilon = exp(-4:1),
           nsimEIM = 500,
           oim.mean = TRUE,
           zero = NULL, nowarning = FALSE) {







  lepsil <- lepsilon
  iepsil <- iepsilon


  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lepsil <- as.list(substitute(lepsil))
  eepsil <- link2list(lepsil)
  lepsil <- attr(eepsil, "function.name")

  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE))
    stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
    warning("argument 'nsimEIM' should be an integer ",
            "greater than 50, say")


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' values must be positive")
  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
      stop("argument 'iscale' values must be positive")
  if (length(iepsil))
    if (!is.Numeric(iepsil, positive = TRUE))
      stop("argument 'iepsil' values must be positive")





    if (!is.logical(oim.mean) || length(oim.mean) != 1)
      stop("bad input for argument 'oim.mean'")




  new("vglmff",
  blurb = c("Makeham distribution\n\n",
            "Links:    ",
            namesof("scale",   lscale, escale), ", ",
            namesof("shape",   lshape, eshape), ", ",
            namesof("epsilon", lepsil, eepsil), "\n",
            "Median:   qmakeham(p = 0.5, scale, shape, epsilon)"),

  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 3
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         nsimEIM = .nsimEIM,
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM ))),
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

    M1 <- 3
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1 <- paste("scale",   if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("shape",   if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames3 <- paste("epsilon", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lscale , .escale , tag = FALSE),
          namesof(mynames2, .lshape , .eshape , tag = FALSE),
          namesof(mynames3, .lepsil , .eepsil , tag = FALSE))[
          interleave.VGAM(M, M = M1)]


    if (!length(etastart)) {

      matC <- matrix(if (length( .iscale )) .iscale else 0 + NA,
                     n, ncoly, byrow = TRUE)
      matH <- matrix(if (length( .ishape )) .ishape else 0 + NA,
                     n, ncoly, byrow = TRUE)

      matE <- matrix(if (length( .iepsil )) .iepsil else 0.3,
                     n, ncoly, byrow = TRUE)


      shape.grid <- unique(sort(c( .gshape )))
      scale.grid <- unique(sort(c( .gscale )))



      for (spp. in 1:ncoly) {
        yvec <- y[, spp.]
        wvec <- w[, spp.]

        makeham.Loglikfun <- function(scaleval, y, x, w, extraargs) {
          ans <-
          sum(c(w) * dmakeham(x = y, shape = extraargs$Shape,
                              epsilon = extraargs$Epsil,
                              scale = scaleval, log = TRUE))
          ans
        }

        mymat <- matrix(-1, length(shape.grid), 2)
        for (jlocal in 1:length(shape.grid)) {
          mymat[jlocal, ] <-
            grid.search(scale.grid,
                        objfun = makeham.Loglikfun,
                        y = yvec, x = x, w = wvec,
                        ret.objfun = TRUE,
                        extraargs = list(Shape = shape.grid[jlocal],
                                         Epsil = matE[1, spp.]))
        }
        index.shape <- which(mymat[, 2] == max(mymat[, 2]))[1]

        if (!length( .ishape ))
          matH[, spp.] <- shape.grid[index.shape]
        if (!length( .iscale ))
          matC[, spp.] <- mymat[index.shape, 1]
      }  # spp.





      epsil.grid <- c( .gepsil )
      for (spp. in 1:ncoly) {
        yvec <- y[, spp.]
        wvec <- w[, spp.]

        makeham.Loglikfun2 <- function(epsilval, y, x, w, extraargs) {
          ans <-
          sum(c(w) * dmakeham(x = y, shape = extraargs$Shape,
                              epsilon = epsilval, 
                              scale = extraargs$Scale, log = TRUE))
          ans
        }
        Init.epsil <-
            grid.search(epsil.grid,
                        objfun = makeham.Loglikfun2,
                        y = yvec, x = x, w = wvec,
                        extraargs = list(Shape = matH[1, spp.],
                                         Scale = matC[1, spp.]))

        matE[, spp.] <- Init.epsil
      }  # spp.


      etastart <- cbind(theta2eta(matC, .lscale , .escale ),
                        theta2eta(matH, .lshape , .eshape ),
                        theta2eta(matE, .lepsil , .eepsil ))[,
                        interleave.VGAM(M, M = M1)]
    }  # End of !length(etastart)
  }), list(
            .lshape = lshape, .lscale = lscale, .lepsil = lepsil,
            .eshape = eshape, .escale = escale, .eepsil = eepsil,
            .gshape = gshape, .gscale = gscale, .gepsil = gepsilon,
            .ishape = ishape, .iscale = iscale, .iepsil = iepsil
          ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    scale <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lscale , .escale )
    shape <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lshape , .eshape )
    epsil <- eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lepsil , .eepsil )
    qmakeham(p = 0.5, scale = scale, shape = shape, epsil = epsil)
  }, list(
            .lshape = lshape, .lscale = lscale, .lepsil = lepsil,
            .eshape = eshape, .escale = escale, .eepsil = eepsil
         ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <-
      c(rep( .lscale , length = ncoly),
        rep( .lshape , length = ncoly),
        rep( .lepsil , length = ncoly))[interleave.VGAM(M, M = M1)]
    temp.names <- c(mynames1, mynames2, mynames3)[
                    interleave.VGAM(M, M = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-2]] <- .escale
      misc$earg[[M1*ii-1]] <- .eshape
      misc$earg[[M1*ii  ]] <- .eepsil
    }

    misc$M1 <- M1
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list(
            .lshape = lshape, .lscale = lscale, .lepsil = lepsil,
            .eshape = eshape, .escale = escale, .eepsil = eepsil,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    scale <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lscale , .escale )
    shape <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lshape , .eshape )
    epsil <- eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lepsil , .eepsil )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dmakeham(x = y, scale = scale, shape = shape,
                                 epsil = epsil, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .lscale = lscale, .lepsil = lepsil,
           .eshape = eshape, .escale = escale, .eepsil = eepsil ))),
  vfamily = c("makeham"),
 



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    Scale <- eta2theta(eta[, c(TRUE, FALSE, FALSE), drop = FALSE],
                       .lscale , .escale )
    shape <- eta2theta(eta[, c(FALSE, TRUE, FALSE), drop = FALSE],
                       .lshape , .eshape )
    epsil <- eta2theta(eta[, c(FALSE, FALSE, TRUE), drop = FALSE],
                       .lepsil , .eepsil )
    rmakeham(nsim * length(Scale),
             scale = c(Scale), shape = c(shape), epsilon = c(epsil))
  }, list( .lshape = lshape, .lscale = lscale, .lepsil = lepsil,
           .eshape = eshape, .escale = escale, .eepsil = eepsil ))),




  deriv = eval(substitute(expression({
    scale <- eta2theta(eta[, c(TRUE, FALSE, FALSE), drop = FALSE],
                       .lscale , .escale )
    shape <- eta2theta(eta[, c(FALSE, TRUE, FALSE), drop = FALSE],
                       .lshape , .eshape )
    epsil <- eta2theta(eta[, c(FALSE, FALSE, TRUE), drop = FALSE],
                       .lepsil , .eepsil )

    temp2 <- exp(y * scale)
    temp3 <- epsil + shape * temp2
    dl.dshape <- temp2 / temp3 - expm1(y * scale) / scale
    dl.dscale <- shape * y * temp2 / temp3 +
                 shape * expm1(y * scale) / scale^2 -
                 shape * y * temp2 / scale
    dl.depsil <- 1 / temp3 - y

    dshape.deta <- dtheta.deta(shape, .lshape , .eshape )
    dscale.deta <- dtheta.deta(scale, .lscale , .escale )
    depsil.deta <- dtheta.deta(epsil, .lepsil , .eepsil )

    dthetas.detas <- cbind(dscale.deta, dshape.deta, depsil.deta)
    myderiv <- c(w) * cbind(dl.dscale,
                            dl.dshape,
                            dl.depsil) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = M1)]
  }), list( .lshape = lshape, .lscale = lscale, .lepsil = lepsil,
            .eshape = eshape, .escale = escale, .eepsil = eepsil ))),

  weight = eval(substitute(expression({
    NOS <- M / M1
    dThetas.detas <- dthetas.detas[, interleave.VGAM(M, M = M1)]
    wz <- matrix(0.0, n, M + M - 1 + M - 2)  # wz has half-bandwidth 3

    ind1 <- iam(NA, NA, M = M1, both = TRUE, diag = TRUE)  # Use simulated EIM

    for (spp. in 1:NOS) {
      run.varcov <- 0
      Shape <- shape[, spp.]
      Scale <- scale[, spp.]
      Epsil <- epsil[, spp.]

      for (ii in 1:( .nsimEIM )) {
        ysim <- rmakeham(n = n, scale = Scale, shape = Shape, epsil = Epsil)
        temp2 <- exp(ysim * Scale)
        temp3 <- Epsil + Shape * temp2
        dl.dshape <- temp2 / temp3 - expm1(ysim * Scale) / Scale
        dl.dscale <- Shape * ysim * temp2 / temp3 +
                     Shape * expm1(ysim * Scale) / Scale^2 -
                     Shape * ysim * temp2 / Scale
        dl.depsil <- 1 / temp3 - ysim

        temp7 <- cbind(dl.dscale, dl.dshape, dl.depsil)
        run.varcov <- run.varcov + temp7[, ind1$row.index] *
                                   temp7[, ind1$col.index]
      }
      run.varcov <- cbind(run.varcov / .nsimEIM )


      wz1 <- if (intercept.only)
        matrix(colMeans(run.varcov, na.rm = TRUE),
               nrow = n, ncol = ncol(run.varcov), byrow = TRUE) else run.varcov

      wz1 <- wz1 * dThetas.detas[, M1 * (spp. - 1) + ind1$row] *
                   dThetas.detas[, M1 * (spp. - 1) + ind1$col]
      for (jay in 1:M1)
        for (kay in jay:M1) {  # Now copy wz1 into wz
          cptr <- iam((spp. - 1) * M1 + jay,
                      (spp. - 1) * M1 + kay, M = M)
          wz[, cptr] <- wz1[, iam(jay, kay, M = M1)]
        }
    }  # End of for (spp.) loop
    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / M1)
  }), list( .lshape = lshape, .lscale = lscale, .lepsil = lepsil,
            .eshape = eshape, .escale = escale, .eepsil = eepsil,
            .nsimEIM = nsimEIM, .oim.mean = oim.mean ))))
}  # makeham()








dgompertz <- function(x, scale = 1, shape, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(scale))
  if (length(x)     != LLL) x     <- rep(x,     length.out = LLL)
  if (length(shape) != LLL) shape <- rep(shape, length.out = LLL)
  if (length(scale) != LLL) scale <- rep(scale, length.out = LLL)


  index0 <- (x < 0)
  index1 <- abs(x * scale) < 0.1 & is.finite(x * scale)
  ans <- log(shape) + x * scale - (shape / scale) * (exp(x * scale) - 1)
  ans[index1] <- log(shape[index1]) + x[index1] * scale[index1] -
                 (shape[index1] / scale[index1]) *
                 expm1(x[index1] * scale[index1])
  ans[index0] <- log(0)
  ans[x == Inf] <- log(0)
  if (log.arg) {
  } else {
    ans <- exp(ans)
    ans[index0] <- 0
    ans[x == Inf] <- 0
  }
  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}



pgompertz <- function(q, scale = 1, shape,
                      lower.tail = TRUE, log.p = FALSE) {


  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  LLL <- max(length(q), length(shape), length(scale))
  if (length(q)       != LLL) q       <- rep(q,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)


  if (lower.tail) {
    if (log.p) {
      ans <- log1p(-exp((-shape / scale) * expm1(scale * q)))
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- -expm1((-shape / scale) * expm1(scale * q))
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- (-shape / scale) * expm1(scale * q)
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      ans <- exp((-shape / scale) * expm1(scale * q))
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  }
  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}



qgompertz <- function(p, scale = 1, shape,
                      lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  LLL <- max(length(p), length(shape), length(scale))
  if (length(p)       != LLL) p       <- rep(p,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)

  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- log1p((-scale / shape) * log(-expm1(ln.p))) / scale
      ans[ln.p > 0] <- NaN
    } else {
      ans <- log1p((-scale / shape) * log1p(-p)) / scale
      ans[p < 0] <- NaN
      ans[p == 0] <- 0
      ans[p == 1] <- Inf
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- log1p((-scale / shape) * ln.p) / scale
      ans[ln.p > 0] <- NaN
    } else {
      ans <- log1p((-scale / shape) * log(p)) / scale
      ans[p < 0] <- NaN
      ans[p == 0] <- Inf
      ans[p == 1] <- 0
      ans[p > 1] <- NaN
    }
  }
  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}





rgompertz <- function(n, scale = 1, shape) {
  qgompertz(runif(n), scale = scale, shape = shape)
}







gompertz.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 gompertz <-
  function(lscale = "loge", lshape = "loge",
           iscale = NULL,   ishape = NULL,
           nsimEIM = 500,
           zero = NULL, nowarning = FALSE) {





  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE))
    stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
    warning("argument 'nsimEIM' should be an integer ",
            "greater than 50, say")


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' values must be positive")
  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
      stop("argument 'iscale' values must be positive")





  new("vglmff",
  blurb = c("Gompertz distribution\n\n",
            "Links:    ",
            namesof("scale", lscale, escale ), ", ",
            namesof("shape", lshape, eshape ), "\n",
            "Median:     scale * log(2 - 1 / shape)"),

  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         nsimEIM = .nsimEIM,
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM ))),
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
    M1 <- 2
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1 <- paste("scale", if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("shape", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lscale , .escale , tag = FALSE),
          namesof(mynames2, .lshape , .eshape , tag = FALSE))[
          interleave.VGAM(M, M = M1)]



    if (!length(etastart)) {

      matH <- matrix(if (length( .ishape )) .ishape else 0 + NA,
                     n, ncoly, byrow = TRUE)
      matC <- matrix(if (length( .iscale )) .iscale else 0 + NA,
                     n, ncoly, byrow = TRUE)

      shape.grid <- c(exp(-seq(4, 0.1, len = 07)), 1,
                      exp( seq(0.1, 4, len = 07)))
      scale.grid <- c(exp(-seq(4, 0.1, len = 07)), 1,
                      exp( seq(0.1, 4, len = 07)))

      for (spp. in 1:ncoly) {
        yvec <- y[, spp.]
        wvec <- w[, spp.]


        gompertz.Loglikfun <- function(scaleval, y, x, w, extraargs) {
          ans <-
          sum(c(w) * dgompertz(x = y, shape = extraargs$Shape,
                               scale = scaleval, log = TRUE))
          ans 
        }

        mymat <- matrix(-1, length(shape.grid), 2)
        for (jlocal in 1:length(shape.grid)) {
          mymat[jlocal, ] <-
            grid.search(scale.grid,
                        objfun = gompertz.Loglikfun,
                        y = yvec, x = x, w = wvec,
                        ret.objfun = TRUE,
                        extraargs = list(Shape = shape.grid[jlocal]))
        }
        index.shape <- which(mymat[, 2] == max(mymat[, 2]))[1]

        if (!length( .ishape ))
          matH[, spp.] <- shape.grid[index.shape]
        if (!length( .iscale ))
          matC[, spp.] <- mymat[index.shape, 1]
      }  # spp.

      etastart <- cbind(theta2eta(matC, .lscale , .escale ),
                        theta2eta(matH, .lshape , .eshape ))[,
                        interleave.VGAM(M, M = M1)]
    }  # End of !length(etastart)
  }), list( .lshape = lshape, .lscale = lscale,
            .eshape = eshape, .escale = escale,
            .ishape = ishape, .iscale = iscale
          ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    scale <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , .escale )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , .eshape )
    log1p((scale / shape) * log(2)) / scale
  }, list( .lshape = lshape, .lscale = lscale,
           .eshape = eshape, .escale = escale ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <-
      c(rep( .lscale , length = ncoly),
        rep( .lshape , length = ncoly))[interleave.VGAM(M, M = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .escale
      misc$earg[[M1*ii  ]] <- .eshape
    }

    misc$M1 <- M1
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list( .lshape = lshape, .lscale = lscale,
            .eshape = eshape, .escale = escale,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    scale <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , .escale )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dgompertz(x = y, scale = scale,
                                  shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
    }, list( .lshape = lshape, .lscale = lscale,
             .eshape = eshape, .escale = escale ))),
  vfamily = c("gompertz"),
 



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    Scale <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , .escale )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , .eshape )
    rgompertz(nsim * length(Scale),
              shape = c(Shape), scale = c(Scale))
    }, list( .lshape = lshape, .lscale = lscale,
             .eshape = eshape, .escale = escale ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    scale <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lscale ,
                       .escale )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lshape ,
                       .eshape )


    temp2 <- exp(y * scale)
    temp4 <- -expm1(y * scale)
    dl.dshape <- 1 / shape + temp4 / scale
    dl.dscale <- y * (1 - shape * temp2 / scale) -
                 shape * temp4 / scale^2

    dscale.deta <- dtheta.deta(scale, .lscale , .escale )
    dshape.deta <- dtheta.deta(shape, .lshape , .eshape )

    dthetas.detas <- cbind(dscale.deta, dshape.deta)
    myderiv <- c(w) * cbind(dl.dscale, dl.dshape) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = M1)]
  }), list( .lshape = lshape, .lscale = lscale,
            .eshape = eshape, .escale = escale ))),


  weight = eval(substitute(expression({

    NOS <- M / M1
    dThetas.detas <- dthetas.detas[, interleave.VGAM(M, M = M1)]

    wz <- matrix(0.0, n, M + M - 1)  # wz is 'tridiagonal' 

    ind1 <- iam(NA, NA, M = M1, both = TRUE, diag = TRUE)


    for (spp. in 1:NOS) {
      run.varcov <- 0
      Shape <- shape[, spp.]
      Scale <- scale[, spp.]

      for (ii in 1:( .nsimEIM )) {
        ysim <- rgompertz(n = n, shape = Shape, scale = Scale)
if (ii < 3) {
}

        temp2 <- exp(ysim * scale)
        temp4 <- -expm1(ysim * scale)
        dl.dshape <- 1 / shape + temp4 / scale
        dl.dscale <- ysim * (1 - shape * temp2 / scale) -
                     shape * temp4 / scale^2


        temp7 <- cbind(dl.dscale, dl.dshape)
        run.varcov <- run.varcov +
                      temp7[, ind1$row.index] *
                      temp7[, ind1$col.index]
      }
      run.varcov <- cbind(run.varcov / .nsimEIM )

      wz1 <- if (intercept.only)
          matrix(colMeans(run.varcov),
                 nrow = n, ncol = ncol(run.varcov), byrow = TRUE) else
          run.varcov

      wz1 <- wz1 * dThetas.detas[, M1 * (spp. - 1) + ind1$row] *
                   dThetas.detas[, M1 * (spp. - 1) + ind1$col]


      for (jay in 1:M1)
        for (kay in jay:M1) {
          cptr <- iam((spp. - 1) * M1 + jay,
                      (spp. - 1) * M1 + kay,
                      M = M)
          wz[, cptr] <- wz1[, iam(jay, kay, M = M1)]
        }
    }  # End of for (spp.) loop



    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / M1)
  }), list( .lscale = lscale,
            .escale = escale,
            .nsimEIM = nsimEIM ))))
}  # gompertz()






 dmoe <- function (x, alpha = 1, lambda = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x),
             length(alpha),
             length(lambda))
  if (length(x)      != LLL) x      <- rep(x,      length.out = LLL)
  if (length(alpha)  != LLL) alpha  <- rep(alpha,  length.out = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, length.out = LLL)

  index0 <- (x < 0)
  if (log.arg) {
    ans <- log(lambda) + (lambda * x) -
           2 * log(expm1(lambda * x) + alpha)
    ans[index0] <- log(0)
  } else {
    ans <- lambda * exp(lambda * x) / (expm1(lambda * x) + alpha)^2
    ans[index0] <- 0
  }
  ans[alpha <= 0 | lambda <= 0] <- NaN
  ans
}



 pmoe <- function (q, alpha = 1, lambda = 1) {
  ret <- ifelse(alpha <= 0 | lambda <= 0, NaN,
                1 - 1 / (expm1(lambda * q) + alpha))
  ret[q < log(2 - alpha) / lambda] <- 0
  ret
}



qmoe <- function (p, alpha = 1, lambda = 1) {
  ifelse(p < 0 | p > 1 | alpha <= 0 | lambda <= 0, NaN,
        log1p(-alpha + 1 / (1 - p)) / lambda)
}



rmoe <- function (n, alpha = 1, lambda = 1) {
  qmoe(p = runif(n), alpha = alpha, lambda = lambda)
}




exponential.mo.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}




 exponential.mo <-
  function(lalpha = "loge", llambda = "loge",
           ealpha = list(), elambda = list(),
           ialpha = 1,      ilambda = NULL,
           imethod = 1,
           nsimEIM = 200,
           zero = NULL) {

  stop("fundamentally unable to estimate the parameters as ",
       "the support of the density depends on the parameters")


  lalpha <- as.list(substitute(lalpha))
  ealpha <- link2list(lalpha)
  lalpha <- attr(ealpha, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  lalpha0 <- lalpha
  ealpha0 <- ealpha
  ialpha0 <- ialpha



  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE))
    stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
    warning("argument 'nsimEIM' should be an integer ",
            "greater than 50, say")

  if (length(ialpha0))
    if (!is.Numeric(ialpha0, positive = TRUE))
      stop("argument 'ialpha' values must be positive")
  if (length(ilambda))
    if (!is.Numeric(ilambda, positive = TRUE))
      stop("argument 'ilambda' values must be positive")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")



  new("vglmff",
  blurb = c("Marshall-Olkin exponential distribution\n\n",
            "Links:    ",
            namesof("alpha",  lalpha0, ealpha0 ), ", ",
            namesof("lambda", llambda, elambda ), "\n",
            "Median:     log(3 - alpha) / lambda"),

  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         nsimEIM = .nsimEIM,
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM ))),
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

    M1 <- 2
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1 <- paste("alpha",   if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("lambda",  if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lalpha0 , .ealpha0 , tag = FALSE),
          namesof(mynames2, .llambda , .elambda , tag = FALSE))[
          interleave.VGAM(M, M = M1)]



    if (!length(etastart)) {

      matL <- matrix(if (length( .ilambda )) .ilambda else 0,
                     n, ncoly, byrow = TRUE)
      matA <- matrix(if (length( .ialpha0 )) .ialpha0 else 0,
                     n, ncoly, byrow = TRUE)


      for (spp. in 1:ncoly) {
        yvec <- y[, spp.]

        moexpon.Loglikfun <- function(lambdaval, y, x, w, extraargs) {
          ans <-
          sum(c(w) * log(dmoe(x = y, alpha = extraargs$alpha,
                              lambda = lambdaval)))
          ans
        }
        Alpha.init <- .ialpha0
        lambda.grid <- seq(0.1, 10.0, len = 21)
        Lambda.init <- grid.search(lambda.grid,
                                   objfun = moexpon.Loglikfun,
                                   y = y, x = x, w = w,
                                   extraargs = list(alpha = Alpha.init))

        if (length(mustart)) {
          Lambda.init <- Lambda.init / (1 - Phimat.init)
        }

        if (!length( .ialpha0 ))
          matA[, spp.] <- Alpha0.init
        if (!length( .ilambda ))
          matL[, spp.] <- Lambda.init
      }  # spp.

      etastart <- cbind(theta2eta(matA, .lalpha0, .ealpha0 ),
                        theta2eta(matL, .llambda, .elambda ))[,
                        interleave.VGAM(M, M = M1)]
      mustart <- NULL # Since etastart has been computed.
    }  # End of !length(etastart)
  }), list( .lalpha0 = lalpha0, .llambda = llambda,
            .ealpha0 = ealpha0, .elambda = elambda,
            .ialpha0 = ialpha0, .ilambda = ilambda,
            .imethod = imethod
          ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha0 <- eta2theta(eta[, c(TRUE, FALSE)], .lalpha0 , .ealpha0 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , .elambda )
    log(3 - alpha0) / lambda
  }, list( .lalpha0 = lalpha0, .llambda = llambda,
           .ealpha0 = ealpha0, .elambda = elambda ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <-
      c(rep( .lalpha0 , length = ncoly),
        rep( .llambda , length = ncoly))[interleave.VGAM(M, M = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .ealpha0
      misc$earg[[M1*ii  ]] <- .elambda
    }

    misc$M1 <- M1
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list( .lalpha0 = lalpha0, .llambda = llambda,
            .ealpha0 = ealpha0, .elambda = elambda,
            .nsimEIM = nsimEIM,
            .imethod = imethod ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    alpha0 <- eta2theta(eta[, c(TRUE, FALSE)], .lalpha0 , .ealpha0 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , .elambda )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * log(dmoe(x = y, alpha = alpha0,
                                 lambda = lambda))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
    }, list( .lalpha0 = lalpha0, .llambda = llambda,
             .ealpha0 = ealpha0, .elambda = elambda ))),
  vfamily = c("exponential.mo"),
 
  deriv = eval(substitute(expression({
    M1 <- 2
    alpha0 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lalpha0 ,
                       .ealpha0 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .llambda ,
                       .elambda )

    temp2 <- (expm1(lambda * y) + alpha0)
    dl.dalpha0 <- -2 / temp2
    dl.dlambda <- 1 / lambda + y - 2 * y * exp(lambda * y) / temp2

    dalpha0.deta <- dtheta.deta(alpha0, .lalpha0 , .ealpha0 )
    dlambda.deta <- dtheta.deta(lambda, .llambda , .elambda )

    dthetas.detas <- cbind(dalpha0.deta,
                           dlambda.deta)
    myderiv <- c(w) * cbind(dl.dalpha0, dl.dlambda) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = M1)]
  }), list( .lalpha0 = lalpha0, .llambda = llambda,
            .ealpha0 = ealpha0, .elambda = elambda ))),


  weight = eval(substitute(expression({

    NOS <- M / M1
    dThetas.detas <- dthetas.detas[, interleave.VGAM(M, M = M1)]

    wz <- matrix(0.0, n, M + M - 1)  # wz is 'tridiagonal' 

    ind1 <- iam(NA, NA, M = M1, both = TRUE, diag = TRUE)


    for (spp. in 1:NOS) {
      run.varcov <- 0
      Alph <- alpha0[, spp.]
      Lamb <- lambda[, spp.]

      for (ii in 1:( .nsimEIM )) {
        ysim <- rmoe(n = n, alpha = Alph, lambda = Lamb)
if (ii < 3) {
}

        temp2 <- (expm1(lambda * ysim) + alpha0)
        dl.dalpha0 <- -2 / temp2
        dl.dlambda <- 1 / lambda + ysim -
                      2 * ysim * exp(lambda * ysim) / temp2


        temp3 <- cbind(dl.dalpha0, dl.dlambda)
        run.varcov <- run.varcov +
                      temp3[, ind1$row.index] *
                      temp3[, ind1$col.index]
      }
      run.varcov <- cbind(run.varcov / .nsimEIM)

      wz1 <- if (intercept.only)
          matrix(colMeans(run.varcov),
                 nrow = n, ncol = ncol(run.varcov), byrow = TRUE) else
          run.varcov

      wz1 <- wz1 * dThetas.detas[, M1 * (spp. - 1) + ind1$row] *
                   dThetas.detas[, M1 * (spp. - 1) + ind1$col]


      for (jay in 1:M1)
        for (kay in jay:M1) {
          cptr <- iam((spp. - 1) * M1 + jay,
                      (spp. - 1) * M1 + kay,
                      M = M)
          wz[, cptr] <- wz1[, iam(jay, kay, M = M1)]
        }
    }  # End of for (spp.) loop




    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / M1)
  }), list( .llambda = llambda,
            .elambda = elambda,
            .nsimEIM = nsimEIM ))))
}  # exponential.mo()









 genbetaII <- function(lss,
                       lshape1.a = "loge",
                       lscale = "loge",
                       lshape2.p = "loge",
                       lshape3.q = "loge",
                       ishape1.a = NULL,
                       iscale = NULL,
                       ishape2.p = 1.0,
                       ishape3.q = 1.0,
                       zero = NULL) {

  if (!is.logical(lss) || lss)
    stop("argument 'lss' not specified correctly. ",
         "See online help for important information")


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  lshape1.a <- as.list(substitute(lshape1.a))
  eshape1.a <- link2list(lshape1.a)
  lshape1.a <- attr(eshape1.a, "function.name")

  lshape2.p <- as.list(substitute(lshape2.p))
  eshape2.p <- link2list(lshape2.p)
  lshape2.p <- attr(eshape2.p, "function.name")

  lshape3.q <- as.list(substitute(lshape3.q))
  eshape3.q <- link2list(lshape3.q)
  lshape3.q <- attr(eshape3.q, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  new("vglmff",
  blurb = c("Generalized Beta II distribution\n\n",
            "Links:    ",
            namesof("shape1.a", lshape1.a, earg = eshape1.a), ", ", 
            namesof("scale",    lscale,    earg = escale), ", ", 
            namesof("shape2.p", lshape2.p, earg = eshape2.p), ", ", 
            namesof("shape3.q", lshape3.q, earg = eshape3.q), "\n", 
            "Mean:     scale * gamma(shape2.p + 1/shape1.a) * ",
                      "gamma(shape3.q - 1/shape1.a) / ",
                      "(gamma(shape2.p) * gamma(shape3.q))"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1)


    predictors.names <-
      c(namesof("shape1.a", .lshape1.a , earg = .eshape1.a , tag = FALSE),
        namesof("scale"   , .lscale    , earg = .escale    , tag = FALSE),
        namesof("shape2.p", .lshape2.p , earg = .eshape2.p , tag = FALSE),
        namesof("shape3.q", .lshape3.q , earg = .eshape3.q , tag = FALSE))

    if (!length( .ishape1.a ) || !length( .iscale )) {
      qvec <- c( .25, .5, .75)  # Arbitrary; could be made an argument
      ishape3.q <- if (length( .ishape3.q)) .ishape3.q else 1
      xvec <- log( (1-qvec)^(-1/ ishape3.q ) - 1 )
      fit0 <- lsfit(x = xvec, y = log(quantile(y, qvec )))
    }

    if (!length(etastart)) {
      aa      <- rep(if (length( .ishape1.a )) .ishape1.a else
                     abs(1 / fit0$coef[2]),
                     length.out = n)
      scale   <- rep(if (length( .iscale )) .iscale else
                     exp(fit0$coef[1]), length = n)
      qq      <- rep(if (length( .ishape3.q )) .ishape3.q else 1.0,
                     length.out = n)
      parg    <- rep(if (length( .ishape2.p )) .ishape2.p else 1.0,
                     length.out = n)


      outOfRange <- (qq - 1/aa <= 0)
      qq[outOfRange] <- 1 / aa[outOfRange] + 1
      outOfRange <- (parg + 1/aa <= 0)
      parg[outOfRange] <- 1 / aa[outOfRange] + 1
    

      etastart <-
        cbind(theta2eta(aa,    .lshape1.a , earg = .eshape1.a ),
              theta2eta(scale, .lscale ,    earg = .escale ),
              theta2eta(parg,  .lshape2.p , earg = .eshape2.p ),
              theta2eta(qq,    .lshape3.q , earg = .eshape3.q ))
    }
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .lshape2.p = lshape2.p, .lshape3.q = lshape3.q,
            .eshape1.a = eshape1.a, .escale = escale, 
            .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
            .ishape1.a = ishape1.a, .iscale = iscale, 
            .ishape2.p = ishape2.p, .ishape3.q = ishape3.q ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    aa     <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale  <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    parg   <- eta2theta(eta[, 3], .lshape2.p , earg = .eshape2.p )
    qq     <- eta2theta(eta[, 4], .lshape3.q , earg = .eshape3.q )
    ans <- Scale * exp(lgamma(parg + 1/aa) +
                       lgamma(qq   - 1/aa) - lgamma(parg) - lgamma(qq))
    ans[parg + 1/aa <= 0] <- NA
    ans[qq   - 1/aa <= 0] <- NA
    ans[aa          <= 0] <- NA
    ans[Scale       <= 0] <- NA
    ans[parg        <= 0] <- NA
    ans[qq          <= 0] <- NA
    ans
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale, 
           .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
           .lshape2.p = lshape2.p, .lshape3.q = lshape3.q ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1.a = .lshape1.a , scale = .lscale ,
                      shape2.p = .lshape2.p , shape3.q = .lshape3.q )

    misc$earg <- list(shape1.a = .eshape1.a , scale = .escale ,
                      shape2.p = .eshape2.p , shape3.q = .eshape3.q )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale, 
            .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
            .lshape2.p = lshape2.p, .lshape3.q = lshape3.q ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    aa     <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    scale  <- eta2theta(eta[, 2], .lscale , earg = .escale )
    parg   <- eta2theta(eta[, 3], .lshape2.p , earg = .eshape2.p )
    qq     <- eta2theta(eta[, 4], .lshape3.q , earg = .eshape3.q )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * (log(aa) + (aa * parg - 1) * log(y) -
                aa * parg * log(scale) +
              - lbeta(parg, qq) - (parg + qq) * log1p((y/scale)^aa))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale, 
           .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
           .lshape2.p = lshape2.p, .lshape3.q = lshape3.q ))),
  vfamily = c("genbetaII"),
  deriv = eval(substitute(expression({
    aa     <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    scale  <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    parg   <- eta2theta(eta[, 3], .lshape2.p , earg = .eshape2.p )
    qq     <- eta2theta(eta[, 4], .lshape3.q , earg = .eshape3.q )

    temp1 <- log(y/scale)
    temp2 <- (y/scale)^aa
    temp3 <- digamma(parg + qq)
    temp3a <- digamma(parg)
    temp3b <- digamma(qq)
    temp4 <- log1p(temp2)

    dl.da <- 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
    dl.dscale <- (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
    dl.dp <- aa * temp1 + temp3 - temp3a - temp4
    dl.dq <- temp3 - temp3b - temp4

    da.deta     <- dtheta.deta(aa,    .lshape1.a , earg = .eshape1.a )
    dscale.deta <- dtheta.deta(scale, .lscale ,    earg = .escale )
    dp.deta     <- dtheta.deta(parg,  .lshape2.p , earg = .eshape2.p )
    dq.deta     <- dtheta.deta(qq,    .lshape3.q , earg = .eshape3.q )

    c(w) * cbind( dl.da * da.deta,
                  dl.dscale * dscale.deta,
                  dl.dp * dp.deta,
                  dl.dq * dq.deta )
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale, 
            .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
            .lshape2.p = lshape2.p, .lshape3.q = lshape3.q ))),
  weight = eval(substitute(expression({
    temp5  <- trigamma(parg + qq)
    temp5a <- trigamma(parg)
    temp5b <- trigamma(qq)

    ned2l.da <- (1 + parg+qq + parg * qq * (temp5a + temp5b +
                (temp3b - temp3a + (parg-qq)/(parg*qq))^2 -
                (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1+parg+qq))
    ned2l.dscale <- aa^2 * parg * qq / (scale^2 * (1+parg+qq))
    ned2l.dp <- temp5a - temp5
    ned2l.dq <- temp5b - temp5
    ned2l.dascale <- (parg - qq - parg * qq *
                     (temp3a -temp3b)) / (scale*(1 + parg+qq))
    ned2l.dap <- -(qq   * (temp3a -temp3b) -1) / (aa*(parg+qq))
    ned2l.daq <- -(parg * (temp3b -temp3a) -1) / (aa*(parg+qq))
    ned2l.dscalep <-  aa * qq   / (scale*(parg+qq))
    ned2l.dscaleq <- -aa * parg / (scale*(parg+qq))
    ned2l.dpq <- -temp5

    wz <- matrix(as.numeric(NA), n, dimm(M))  # M==4 means 10=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.da * da.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dscale * dscale.deta^2
    wz[, iam(3, 3, M)] <- ned2l.dp * dp.deta^2
    wz[, iam(4, 4, M)] <- ned2l.dq * dq.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dascale * da.deta * dscale.deta
    wz[, iam(1, 3, M)] <- ned2l.dap * da.deta * dp.deta
    wz[, iam(1, 4, M)] <- ned2l.daq * da.deta * dq.deta
    wz[, iam(2, 3, M)] <- ned2l.dscalep * dscale.deta * dp.deta
    wz[, iam(2, 4, M)] <- ned2l.dscaleq * dscale.deta * dq.deta
    wz[, iam(3, 4, M)] <- ned2l.dpq * dp.deta * dq.deta
    wz <- c(w) * wz
    wz
  }), list( .lshape1.a <- lshape1.a, .lscale = lscale,
            .eshape1.a <- eshape1.a, .escale = escale, 
            .eshape2.p <- eshape2.p, .eshape3.q = eshape3.q,
            .lshape2.p <- lshape2.p, .lshape3.q = lshape3.q ))))
}


rsinmad <- function(n, scale = 1, shape1.a, shape3.q)
  qsinmad(runif(n), shape1.a = shape1.a, scale = scale,
          shape3.q = shape3.q)


rlomax <- function(n, scale = 1, shape3.q)
  rsinmad(n, scale = scale, shape1.a = 1, shape3.q = shape3.q)


rfisk <- function(n, scale = 1, shape1.a)
  rsinmad(n, scale = scale, shape1.a = shape1.a, shape3.q = 1)


rparalogistic <- function(n, scale = 1, shape1.a)
  rsinmad(n, scale = scale, shape1.a = shape1.a, shape3.q = shape1.a)


rdagum <- function(n, scale = 1, shape1.a, shape2.p)
  qdagum(runif(n), scale = scale, shape1.a = shape1.a,
         shape2.p = shape2.p)


rinv.lomax <- function(n, scale = 1, shape2.p)
  rdagum(n, scale = scale, shape1.a = 1, shape2.p = shape2.p)


rinv.paralogistic <- function(n, scale = 1, shape1.a)
  rdagum(n, scale = scale, shape1.a = shape1.a, shape2.p = shape1.a)




qsinmad <- function(p, scale = 1, shape1.a, shape3.q,
                    lower.tail = TRUE,
                    log.p = FALSE) {



  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")



  LLL <- max(length(p), length(shape1.a), length(scale),
                        length(shape3.q))
  if (length(p) != LLL)
    p         <- rep(p,         length.out = LLL)
  if (length(shape1.a) != LLL)
    shape1.a  <- rep(shape1.a,  length.out = LLL)
  if (length(scale) != LLL)
    scale     <- rep(scale,     length.out = LLL)
  if (length(shape3.q) != LLL)
    shape3.q  <- rep(shape3.q,  length.out = LLL)


  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- scale * expm1((-1/shape3.q) * log(-expm1(ln.p)))^(1/shape1.a)
    } else {
      ans <- scale * expm1((-1/shape3.q) * log1p(-p))^(1/shape1.a)
      ans[p == 0] <- 0
      ans[p == 1] <- Inf
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- scale * expm1(-ln.p / shape3.q)^(1/shape1.a)
    } else {
      ans <- scale * expm1(-log(p) / shape3.q)^(1/shape1.a)
      ans[p == 0] <- Inf
      ans[p == 1] <- 0
    }
  }

  ans[scale    <= 0 | shape1.a <= 0 | shape3.q <= 0] <- NaN
  ans
}




qlomax <- function(p, scale = 1, shape3.q,
                   lower.tail = TRUE, log.p = FALSE)
  qsinmad(p, shape1.a = 1, scale = scale, shape3.q = shape3.q,
          lower.tail = lower.tail, log.p = log.p)

qfisk <- function(p, scale = 1, shape1.a,
                  lower.tail = TRUE, log.p = FALSE)
  qsinmad(p, shape1.a = shape1.a, scale = scale, shape3.q = 1,
          lower.tail = lower.tail, log.p = log.p)

qparalogistic <- function(p, scale = 1, shape1.a,
                          lower.tail = TRUE, log.p = FALSE)
  qsinmad(p, shape1.a = shape1.a, scale = scale,
          shape3.q = shape1.a,  ## 20150121 KaiH; add shape3.q = shape1.a
          lower.tail = lower.tail, log.p = log.p)





qdagum <- function(p, scale = 1, shape1.a, shape2.p,
                   lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  LLL <- max(length(p), length(shape1.a), length(scale),
                        length(shape2.p))
  if (length(p) != LLL)
    p         <- rep(p,         length.out = LLL)
  if (length(shape1.a) != LLL)
    shape1.a  <- rep(shape1.a,  length.out = LLL)
  if (length(scale) != LLL)
    scale     <- rep(scale,     length.out = LLL)
  if (length(shape2.p) != LLL)
    shape2.p  <- rep(shape2.p,  length.out = LLL)

  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- scale * (expm1(-ln.p/shape2.p))^(-1/shape1.a)
      ans[ln.p > 0] <- NaN
    } else {
      ans <- scale * (expm1(-log(p)/shape2.p))^(-1/shape1.a)
      ans[p < 0] <- NaN
      ans[p == 0] <- 0
      ans[p == 1] <- Inf
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- scale * (expm1(-log(-expm1(ln.p))/shape2.p))^(-1/shape1.a)
      ans[ln.p > 0] <- NaN
    } else {
      ans <- scale * (expm1(-log1p(-p)/shape2.p))^(-1/shape1.a)
      ans[p < 0] <- NaN
      ans[p == 0] <- Inf
      ans[p == 1] <- 0
      ans[p > 1] <- NaN
    }
  }

  ans[scale <= 0 | shape1.a <= 0 | shape2.p <= 0] <- NaN
  ans
}




qinv.lomax <- function(p, scale = 1, shape2.p,
                       lower.tail = TRUE, log.p = FALSE)
  qdagum(p, scale = scale, shape1.a = 1, shape2.p = shape2.p,
         lower.tail = lower.tail, log.p = log.p)


qinv.paralogistic <- function(p, scale = 1, shape1.a,
                              lower.tail = TRUE, log.p = FALSE)
  qdagum(p, scale = scale, shape1.a = shape1.a,
         shape2.p = shape1.a,   ##  20150121 Kai; add shape2.p = shape1.a
         lower.tail = lower.tail, log.p = log.p)




psinmad <- function(q, scale = 1, shape1.a, shape3.q,
                    lower.tail = TRUE, log.p = FALSE) {


  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  LLL <- max(length(q), length(shape1.a), length(scale),
                        length(shape3.q))
  if (length(q) != LLL)
    q         <- rep(q,         length.out = LLL)
  if (length(shape1.a) != LLL)
    shape1.a  <- rep(shape1.a,  length.out = LLL)
  if (length(scale) != LLL)
    scale     <- rep(scale,     length.out = LLL)
  if (length(shape3.q) != LLL)
    shape3.q  <- rep(shape3.q,  length.out = LLL)

  # 20150121 KaiH
  if (lower.tail) {
    if (log.p) {
      ans <- log1p(-(1 + (q / scale)^shape1.a)^(-shape3.q))
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- exp(log1p(-(1 + (q / scale)^shape1.a)^(-shape3.q)))
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- (-shape3.q) * log1p((q / scale)^shape1.a)
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      ans <- (1 + (q / scale)^shape1.a)^(-shape3.q)
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  }
  ans[scale <= 0 | shape1.a <= 0 | shape3.q <= 0] <- NaN
  ans
}






plomax <- function(q, scale = 1, shape3.q,   # Change the order
                   lower.tail = TRUE, log.p = FALSE)
  psinmad(q, shape1.a = 1, scale = scale, shape3.q = shape3.q,
          lower.tail = lower.tail, log.p = log.p)


pfisk <- function(q, scale = 1, shape1.a,
                  lower.tail = TRUE, log.p = FALSE)
  psinmad(q, shape1.a = shape1.a, scale = scale, shape3.q = 1,
          lower.tail = lower.tail, log.p = log.p)


pparalogistic <- function(q, scale = 1, shape1.a,   # Change the order
                          lower.tail = TRUE, log.p = FALSE)
  psinmad(q, shape1.a = shape1.a, scale = scale,
          shape3.q = shape1.a,  # Add shape3.q = shape1.a
          lower.tail = lower.tail, log.p = log.p)






pdagum <- function(q, scale = 1, shape1.a, shape2.p,
                   lower.tail = TRUE, log.p = FALSE) {



  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  LLL <- max(length(q), length(shape1.a), length(scale),
                        length(shape2.p))
  if (length(q) != LLL)
    q         <- rep(q,         length.out = LLL)
  if (length(shape1.a) != LLL)
    shape1.a  <- rep(shape1.a,  length.out = LLL)
  if (length(scale) != LLL)
    scale     <- rep(scale,     length.out = LLL)
  if (length(shape2.p) != LLL)
    shape2.p  <- rep(shape2.p,  length.out = LLL)


  if (lower.tail) {
    if (log.p) {
      ans <- (-shape2.p) * log1p((q/scale)^(-shape1.a))
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- exp( (-shape2.p) * log1p((q/scale)^(-shape1.a)) )
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- log1p(-(1 + (q/scale)^(-shape1.a))^(-shape2.p))
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      stop("unfinished")
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  }

  ans[shape1.a <= 0 | scale <= 0 | shape2.p <= 0] <- NaN
  ans
}





pinv.lomax <- function(q, scale = 1, shape2.p,
                       lower.tail = TRUE, log.p = FALSE)
  pdagum(q, scale = scale, shape1.a = 1, shape2.p = shape2.p,
         lower.tail = lower.tail, log.p = log.p)


pinv.paralogistic <- function(q, scale = 1, shape1.a,
                              lower.tail = TRUE, log.p = FALSE)
  pdagum(q, scale = scale, shape1.a = shape1.a,
         shape2.p = shape1.a,
         lower.tail = lower.tail, log.p = log.p)



dsinmad <- function(x, scale = 1, shape1.a, shape3.q, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL      <- max(length(x),     length(shape1.a),
                  length(scale), length(shape3.q))
  x        <- rep(x,         length.out = LLL);
  shape1.a <- rep(shape1.a,  length.out = LLL)
  scale    <- rep(scale,     length.out = LLL);
  shape3.q <- rep(shape3.q,  length.out = LLL)

  Loglik <- rep(log(0), length.out = LLL)
  xok <- (x > 0) & !is.na(x)  # Avoids log(x) if x<0, and handles NAs
  Loglik[xok] <- log(shape1.a[xok]) + log(shape3.q[xok]) +
                 (shape1.a[xok]-1) * log(x[xok]) -
                shape1.a[xok] * log(scale[xok]) -
           (1 + shape3.q[xok]) * log1p((x[xok]/scale[xok])^shape1.a[xok])
  x.eq.0 <- (x == 0) & !is.na(x)
  Loglik[x.eq.0] <- log(shape1.a[x.eq.0]) + log(shape3.q[x.eq.0]) -
                    shape1.a[x.eq.0] * log(scale[x.eq.0])
  Loglik[is.na(x)]  <- NA
  Loglik[is.nan(x)] <- NaN
  Loglik[x == Inf]  <- log(0)

  if (log.arg) Loglik else exp(Loglik)
}



dlomax <- function(x, scale = 1, shape3.q, log = FALSE)
  dsinmad(x, scale = scale, shape1.a = 1, shape3.q = shape3.q, log = log)



dfisk <- function(x, scale = 1, shape1.a, log = FALSE)
  dsinmad(x, scale = scale, shape1.a = shape1.a, shape3.q = 1, log = log)



dparalogistic <- function(x, scale = 1, shape1.a, log = FALSE)
  dsinmad(x, scale = scale, shape1.a = shape1.a, shape3.q = shape1.a,
          log = log)



ddagum <- function(x, scale = 1, shape1.a, shape2.p, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x),
             length(shape1.a),
             length(scale),
             length(shape2.p))
  x        <- rep(x,        length.out = LLL)
  shape1.a <- rep(shape1.a, length.out = LLL)
  scale    <- rep(scale,    length.out = LLL)
  shape2.p <- rep(shape2.p, length.out = LLL)

  Loglik <- rep(log(0), length.out = LLL)
  xok <- (x > 0) & !is.na(x)  # Avoids log(x) if x<0, and handles NAs
  Loglik[xok] <- log(shape1.a[xok]) +
                 log(shape2.p[xok]) +
                 (shape1.a[xok] * shape2.p[xok]-1) * log(    x[xok]) -
                  shape1.a[xok] * shape2.p[xok]    * log(scale[xok]) -
                 (1 + shape2.p[xok]) *
                 log1p((x[xok]/scale[xok])^shape1.a[xok])
  Loglik[shape2.p <= 0] <- NaN

  x.eq.0 <- (x == 0) & !is.na(x)
  Loglik[x.eq.0] <- log(shape1.a[x.eq.0]) +
                    log(shape2.p[x.eq.0]) -
                    shape1.a[x.eq.0] * shape2.p[x.eq.0] *
                    log(scale[x.eq.0])
  Loglik[is.na(x)]  <- NA
  Loglik[is.nan(x)] <- NaN
  Loglik[x == Inf]  <- log(0)

  if (log.arg) Loglik else exp(Loglik)
}



dinv.lomax <- function(x, scale = 1, shape2.p, log = FALSE)
  ddagum(x, scale = scale, shape1.a = 1, shape2.p = shape2.p,
         log = log)


dinv.paralogistic <- function(x, scale = 1, shape1.a, log = FALSE)
  ddagum(x, scale = scale, shape1.a = shape1.a, shape2.p = shape1.a,
         log = log)



 sinmad <- function(lss,
                    lshape1.a = "loge",
                    lscale = "loge",
                    lshape3.q = "loge",
                    ishape1.a = NULL, 
                    iscale = NULL,
                    ishape3.q = 1.0, 
                    zero = NULL) {

  if (!is.logical(lss) || lss)
    stop("argument 'lss' not specified correctly. ",
         "See online help for important information")




  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  lshape1.a <- as.list(substitute(lshape1.a))
  eshape1.a <- link2list(lshape1.a)
  lshape1.a <- attr(eshape1.a, "function.name")

  lshape3.q <- as.list(substitute(lshape3.q))
  eshape3.q <- link2list(lshape3.q)
  lshape3.q <- attr(eshape3.q, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  new("vglmff",
  blurb = c("Singh-Maddala distribution\n\n",
          "Links:    ",
          namesof("shape1.a", lshape1.a, earg = eshape1.a), ", ", 
          namesof("scale",    lscale,    earg = escale),    ", ", 
          namesof("shape3.q", lshape3.q, earg = eshape3.q), "\n", 
          "Mean:     scale * gamma(1 + 1/shape1.a) * ",
                    "gamma(shape3.q - 1/shape1.a) / ",
                    "gamma(shape3.q)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1)


    predictors.names <-
        c(namesof("shape1.a", .lshape1.a , earg = .eshape1.a , tag = FALSE),
          namesof("scale",    .lscale ,    earg = .escale ,    tag = FALSE),
          namesof("shape3.q", .lshape3.q , earg = .eshape3.q , tag = FALSE))
    parg <- 1

    if (!length( .ishape1.a) || !length( .iscale )) {
        qvec <- c( .25, .5, .75)  # Arbitrary; could be made an argument
        ishape3.q <- if (length( .ishape3.q)) .ishape3.q else 1
        xvec <- log( (1-qvec)^(-1/ ishape3.q ) - 1 )
        fit0 <- lsfit(x = xvec, y = log(quantile(y, qvec )))
    }



    if (!length(etastart)) {
      aa    <- rep(if (length( .ishape1.a)) .ishape1.a else 1/fit0$coef[2],
                   length.out = n)
      scale <- rep(if (length( .iscale )) .iscale else exp(fit0$coef[1]),
                   length.out = n)
      qq    <- rep(if (length( .ishape3.q)) .ishape3.q else 1.0,
                   length.out = n)


      outOfRange <- (aa * qq <= 1)
      qq[outOfRange] <- 1 / aa[outOfRange] + 1 


        etastart <-
          cbind(theta2eta(aa,    .lshape1.a , earg = .eshape1.a ),
                theta2eta(scale, .lscale    , earg = .escale ),
                theta2eta(qq,    .lshape3.q , earg = .eshape3.q ))
    }
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .lshape3.q = lshape3.q,
            .eshape1.a = eshape1.a, .escale = escale, 
            .eshape3.q = eshape3.q,
            .ishape1.a = ishape1.a, .iscale = iscale, 
            .ishape3.q = ishape3.q ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    aa     <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale  <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    parg   <- 1
    qq     <- eta2theta(eta[, 3], .lshape3.q , earg = .eshape3.q )

    ans <- Scale * exp(lgamma(parg + 1/aa) +
                      lgamma(qq   - 1/aa) - lgamma(parg) - lgamma(qq))
    ans[parg + 1/aa <= 0] <- NA
    ans[qq   - 1/aa <= 0] <- NA
    ans[aa          <= 0] <- NA
    ans[Scale       <= 0] <- NA
    ans[qq          <= 0] <- NA
    ans
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale, 
           .eshape3.q = eshape3.q,
           .lshape3.q = lshape3.q ))),

  last = eval(substitute(expression({
    misc$link <-
         c(shape1.a = .lshape1.a ,
           scale    = .lscale ,
           shape3.q = .lshape3.q )

    misc$earg <-
      list(shape1.a = .eshape1.a ,
           scale    = .escale ,
           shape3.q = .eshape3.q )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,

            .eshape1.a = eshape1.a, .escale = escale, 
            .eshape3.q = eshape3.q,
            .lshape3.q = lshape3.q ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    scale <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    parg <- 1
    qq    <- eta2theta(eta[, 3], .lshape3.q, earg = .eshape3.q )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dsinmad(x = y, shape1.a = aa, scale = scale,
                                shape3.q = qq, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .lshape3.q = lshape3.q,
           .eshape1.a = eshape1.a, .escale = escale,
           .eshape3.q = eshape3.q ))),
  vfamily = c("sinmad"),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    scale <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    qq    <- eta2theta(eta[, 3], .lshape3.q,  earg = .eshape3.q )
    rsinmad(nsim * length(qq), shape1.a = aa, scale = scale,
            shape3.q = qq)
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .lshape3.q = lshape3.q,
           .eshape1.a = eshape1.a, .escale = escale,
           .eshape3.q = eshape3.q ))),


  deriv = eval(substitute(expression({
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    scale <- eta2theta(eta[, 2], .lscale    , earg = .escale )
    parg <- 1
    qq <-    eta2theta(eta[, 3], .lshape3.q , earg = .eshape3.q )

    temp1 <- log(y/scale)
    temp2 <- (y/scale)^aa
    temp3a <- digamma(parg)
    temp3b <- digamma(qq)

    dl.da <- 1 / aa + parg * temp1 - (parg + qq) * temp1 / (1 + 1 / temp2)
    dl.dscale <- (aa / scale) * (-parg + (parg + qq) / (1 + 1 / temp2))
    dl.dq <- digamma(parg + qq) - temp3b - log1p(temp2)

    da.deta     <- dtheta.deta(aa,    .lshape1.a , earg = .eshape1.a )
    dscale.deta <- dtheta.deta(scale, .lscale ,    earg = .escale )
    dq.deta     <- dtheta.deta(qq,    .lshape3.q , earg = .eshape3.q )

    c(w) * cbind(dl.da     * da.deta,
                 dl.dscale * dscale.deta,
                 dl.dq     * dq.deta )
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale, 
            .eshape3.q = eshape3.q,
            .lshape3.q = lshape3.q ))),

  weight = eval(substitute(expression({
    ned2l.da <- (1 + parg+qq + parg * qq * (trigamma(parg) + trigamma(qq) +
                (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1 + parg + qq))
    ned2l.dscale <- aa^2 * parg * qq / (scale^2 * (1+parg+qq))
    ned2l.dq <- 1/qq^2
    ned2l.dascale <- (parg - qq - parg*qq *
                     (temp3a -temp3b)) / (scale*(1 + parg+qq))
    ned2l.daq <- -(parg * (temp3b -temp3a) -1) / (aa*(parg+qq))
    ned2l.dscaleq <- -aa * parg / (scale*(parg+qq))

    wz <- matrix(as.numeric(NA), n, dimm(M))  #M==3 means 6=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.da * da.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dscale * dscale.deta^2
    wz[, iam(3, 3, M)] <- ned2l.dq * dq.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dascale * da.deta * dscale.deta
    wz[, iam(1, 3, M)] <- ned2l.daq * da.deta * dq.deta
    wz[, iam(2, 3, M)] <- ned2l.dscaleq * dscale.deta * dq.deta
    wz <- c(w) * wz
    wz
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale, 
            .eshape3.q = eshape3.q,
            .lshape3.q = lshape3.q ))))
}


 dagum <- function(lss,
                   lshape1.a = "loge",
                   lscale = "loge",
                   lshape2.p = "loge",
                   ishape1.a = NULL, 
                   iscale = NULL,
                   ishape2.p = 1.0, 
                   zero = NULL) {

  if (!is.logical(lss) || lss)
    stop("argument 'lss' not specified correctly. ",
         "See online help for important information")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
      stop("bad input for argument 'zero'")


  lshape1.a <- as.list(substitute(lshape1.a))
  eshape1.a <- link2list(lshape1.a)
  lshape1.a <- attr(eshape1.a, "function.name")

  lshape2.p <- as.list(substitute(lshape2.p))
  eshape2.p <- link2list(lshape2.p)
  lshape2.p <- attr(eshape2.p, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  new("vglmff",
  blurb = c("Dagum distribution\n\n",
          "Links:    ",
          namesof("shape1.a",     lshape1.a,     earg = eshape1.a), ", ", 
          namesof("scale",        lscale,        earg = escale), ", ", 
          namesof("shape2.p",     lshape2.p,     earg = eshape2.p), "\n", 
        "Mean:     scale * gamma(shape2.p + 1/shape1.a) * ",
                  "gamma(1 - 1/shape1.a) / ",
                  "gamma(shape2.p)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1)

    predictors.names <-
      c(namesof("shape1.a", .lshape1.a , earg = .eshape1.a , tag = FALSE),
        namesof("scale",    .lscale ,    earg = .escale ,    tag = FALSE),
        namesof("shape2.p", .lshape2.p , earg = .eshape2.p , tag = FALSE))

    if (!length( .ishape1.a) || !length( .iscale )) {
        qvec <- c( .25, .5, .75)  # Arbitrary; could be made an argument
        ishape2.p <- if (length( .ishape2.p)) .ishape2.p else 1
        xvec <- log( qvec^(-1/ ishape2.p ) - 1 )
        fit0 <- lsfit(x = xvec, y = log(quantile(y, qvec )))
    }

        if (!length(etastart)) {
          parg <- rep(if (length( .ishape2.p )) .ishape2.p else 1.0,
                      length.out = n)
          aa <- rep(if (length( .ishape1.a )) .ishape1.a else
                    -1/fit0$coef[2],
                   length.out = n)
          scale <- rep(if (length( .iscale )) .iscale else
                       exp(fit0$coef[1]),
                       length.out = n)


      outOfRange <- (parg + 1/aa <= 0)
      parg[outOfRange] <- 1 / aa[outOfRange] + 1
    




      etastart <-
        cbind(theta2eta(aa,    .lshape1.a ,  earg = .eshape1.a ),
              theta2eta(scale, .lscale ,     earg = .escale ),
              theta2eta(parg,  .lshape2.p ,  earg = .eshape2.p ))
    }
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .lshape2.p = lshape2.p,
            .eshape1.a = eshape1.a, .escale = escale, 
            .eshape2.p = eshape2.p,
            .ishape1.a = ishape1.a, .iscale = iscale, 
            .ishape2.p = ishape2.p ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    aa     <- eta2theta(eta[, 1], .lshape1.a ,  earg = .eshape1.a )
    Scale  <- eta2theta(eta[, 2], .lscale ,     earg = .escale )
    parg   <- eta2theta(eta[, 3], .lshape2.p ,  earg = .eshape2.p )
    qq     <- 1

    ans <- Scale * exp(lgamma(parg + 1/aa) +
                       lgamma(qq   - 1/aa) - lgamma(parg) - lgamma(qq))
    ans[parg + 1/aa <= 0] <- NA
    ans[qq   - 1/aa <= 0] <- NA
    ans[aa          <= 0] <- NA
    ans[Scale       <= 0] <- NA
    ans[parg        <= 0] <- NA
    ans
    }, list( .lshape1.a = lshape1.a, .lscale = lscale,
             .eshape1.a = eshape1.a, .escale = escale, 
             .eshape2.p = eshape2.p,
             .lshape2.p = lshape2.p ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1.a = .lshape1.a , scale = .lscale ,
                      p = .lshape2.p )

    misc$earg <- list(shape1.a = .eshape1.a , scale = .escale ,
                      p = .eshape2.p )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .lshape2.p = lshape2.p,
            .eshape1.a = eshape1.a, .escale = escale,
            .eshape2.p = eshape2.p ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    aa    <- eta2theta(eta[, 1], .lshape1.a ,  earg = .eshape1.a )
    Scale <- eta2theta(eta[, 2], .lscale ,     earg = .escale )
    parg  <- eta2theta(eta[, 3], .lshape2.p ,  earg = .eshape2.p )
    qq <- 1
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * ddagum(x = y, shape1.a = aa, scale = Scale,
                               shape2.p = parg, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .lshape2.p = lshape2.p, 
           .eshape1.a = eshape1.a, .escale = escale,
           .eshape2.p = eshape2.p ))),
  vfamily = c("dagum"),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    parg  <- eta2theta(eta[, 3], .lshape2.p,  earg = .eshape2.p )
    rdagum(nsim * length(parg), shape1.a = aa, scale = Scale,
           shape2.p = parg)
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .lshape2.p = lshape2.p, 
           .eshape1.a = eshape1.a, .escale = escale,
           .eshape2.p = eshape2.p ))),


  deriv = eval(substitute(expression({
    aa    <- eta2theta(eta[, 1], .lshape1.a ,  earg = .eshape1.a )
    Scale <- eta2theta(eta[, 2], .lscale ,     earg = .escale )
    parg  <- eta2theta(eta[, 3], .lshape2.p ,  earg = .eshape2.p )
    qq <- 1

    temp1 <- log(y / Scale)
    temp2 <- (y / Scale)^aa
    temp3a <- digamma(parg)
    temp3b <- digamma(qq)

    dl.da <- 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
    dl.dscale <- (aa/Scale) * (-parg + (parg+qq) / (1+1/temp2))
    dl.dp <- aa * temp1 + digamma(parg + qq) - temp3a - log1p(temp2)

    da.deta     <- dtheta.deta(aa,    .lshape1.a , earg = .eshape1.a )
    dscale.deta <- dtheta.deta(Scale, .lscale ,    earg = .escale )
    dp.deta     <- dtheta.deta(parg,  .lshape2.p , earg = .eshape2.p )

    c(w) * cbind( dl.da     * da.deta,
                  dl.dscale * dscale.deta,
                  dl.dp     * dp.deta )
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale, 
            .eshape2.p = eshape2.p,
            .lshape2.p = lshape2.p ))),
  weight = eval(substitute(expression({
    ned2l.da <- (1 + parg+qq + parg * qq * (trigamma(parg) + trigamma(qq) +
                (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1+parg+qq))
    ned2l.dscale <- aa^2 * parg * qq / (Scale^2 * (1+parg+qq))
    ned2l.dp <- 1 / parg^2 
    ned2l.dascale <- (parg - qq - parg * qq *(temp3a -temp3b)
                   ) / (Scale * (1 + parg+qq))
    ned2l.dap= -(qq   * (temp3a -temp3b) -1) / (aa*(parg+qq))
    ned2l.dscalep <-  aa * qq   / (Scale * (parg + qq))
    wz <- matrix(as.numeric(NA), n, dimm(M))  #M==3 means 6=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.da     * da.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dscale * dscale.deta^2
    wz[, iam(3, 3, M)] <- ned2l.dp     * dp.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dascale * da.deta * dscale.deta
    wz[, iam(1, 3, M)] <- ned2l.dap * da.deta * dp.deta
    wz[, iam(2, 3, M)] <- ned2l.dscalep * dscale.deta * dp.deta
    wz <- c(w) * wz
    wz
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale, 
            .eshape2.p = eshape2.p,
            .lshape2.p = lshape2.p ))))
}



 betaII <-
   function(  # lss,
            lscale = "loge", lshape2.p = "loge", lshape3.q = "loge",
            iscale = NULL,   ishape2.p = 2,      ishape3.q = 2,
            zero = NULL) {



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")



  lshape2.p <- as.list(substitute(lshape2.p))
  eshape2.p <- link2list(lshape2.p)
  lshape2.p <- attr(eshape2.p, "function.name")

  lshape3.q <- as.list(substitute(lshape3.q))
  eshape3.q <- link2list(lshape3.q)
  lshape3.q <- attr(eshape3.q, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  new("vglmff",
  blurb = c("Beta II distribution\n\n",
          "Links:    ",
          namesof("scale",    lscale,    earg = escale), ", ", 
          namesof("shape2.p", lshape2.p, earg = eshape2.p), ", ", 
          namesof("shape3.q", lshape3.q, earg = eshape3.q), "\n", 
        "Mean:     scale * gamma(shape2.p + 1) * ",
                  "gamma(shape3.q - 1) / ",
                  "(gamma(shape2.p) * gamma(shape3.q))"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1)


    predictors.names <-
    c(namesof("scale",    .lscale ,   earg = .escale ,   tag = FALSE),
      namesof("shape2.p", .lshape2.p, earg = .eshape2.p, tag = FALSE),
      namesof("shape3.q", .lshape3.q, earg = .eshape3.q, tag = FALSE))

    if (!length( .iscale )) {
      qvec <- c(0.25, .5, .75)  # Arbitrary; could be made an argument
      ishape3.q <- if (length( .ishape3.q)) .ishape3.q else 1
      xvec <- log( (1-qvec)^(-1/ ishape3.q ) - 1 )
      fit0 <- lsfit(x = xvec, y = log(quantile(y, qvec )))
    }

    if (!length(etastart)) {
      scale <- rep(if (length( .iscale )) .iscale else
                   exp(fit0$coef[1]),
                   length.out = n)
      qq    <- rep(if (length( .ishape3.q)) .ishape3.q else 1.0,
                   length.out = n)
      parg  <- rep(if (length( .ishape2.p)) .ishape2.p else 1.0,
                   length.out = n)



      aa     <- 1
      outOfRange <- (parg + 1/aa <= 0)
      parg[outOfRange] <- 1 / aa[outOfRange] + 1
      outOfRange <- (qq   - 1/aa <= 0)
      qq[outOfRange] <- 1 / aa + 1


      etastart <-
        cbind(theta2eta(scale, .lscale ,   earg = .escale ),
              theta2eta(parg,  .lshape2.p, earg = .eshape2.p ),
              theta2eta(qq,    .lshape3.q, earg = .eshape3.q))
    }
  }), list( .lscale = lscale,
            .escale = escale, 
            .lshape2.p = lshape2.p, .lshape3.q = lshape3.q,
            .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
            .iscale = iscale, 
            .ishape2.p = ishape2.p,
            .ishape3.q = ishape3.q ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    aa     <- 1
    Scale  <- eta2theta(eta[, 1], .lscale , earg = .escale )
    parg   <- eta2theta(eta[, 2], .lshape2.p, earg = .eshape2.p )
    qq     <- eta2theta(eta[, 3], .lshape3.q, earg = .eshape3.q)

    ans <- Scale * exp(lgamma(parg + 1/aa) +
                       lgamma(qq   - 1/aa) - lgamma(parg) - lgamma(qq))
    ans[parg + 1/aa <= 0] <- NA
    ans[qq   - 1/aa <= 0] <- NA
    ans[Scale       <= 0] <- NA
    ans[parg        <= 0] <- NA
    ans[qq          <= 0] <- NA
    ans
    }, list( .lscale = lscale,
             .escale = escale, 
             .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
             .lshape2.p = lshape2.p, .lshape3.q = lshape3.q ))),
  last = eval(substitute(expression({
    misc$link <-    c(scale    = .lscale ,
                      shape2.p = .lshape2.p,
                      shape3.q = .lshape3.q)

    misc$earg <- list(scale    = .escale ,
                      shape2.p = .eshape2.p,
                      shape3.q = .eshape3.q)

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list(
    .lscale = lscale, .lshape2.p = lshape2.p, .lshape3.q = lshape3.q,
    .escale = escale, .eshape2.p = eshape2.p, .eshape3.q = eshape3.q ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    aa <- 1
    scale <- eta2theta(eta[, 1], .lscale ,   earg = .escale )
    parg  <- eta2theta(eta[, 2], .lshape2.p, earg = .eshape2.p )
    qq    <- eta2theta(eta[, 3], .lshape3.q, earg = .eshape3.q)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * (log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
               (-lbeta(parg, qq)) - (parg+qq)*log1p((y/scale)^aa))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lscale = lscale,
           .escale = escale, 
           .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
           .lshape2.p = lshape2.p, .lshape3.q = lshape3.q ))),
  vfamily = c("betaII"),
  deriv = eval(substitute(expression({
    aa <- 1
    scale <- eta2theta(eta[, 1], .lscale   , earg = .escale )
    parg  <- eta2theta(eta[, 2], .lshape2.p, earg = .eshape2.p )
    qq    <- eta2theta(eta[, 3], .lshape3.q, earg = .eshape3.q)

    temp1 <- log(y/scale)
    temp2 <- (y/scale)^aa
    temp3 <- digamma(parg + qq)
    temp3a <- digamma(parg)
    temp3b <- digamma(qq)
    temp4 <- log1p(temp2)

    dl.dscale <- (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
    dl.dp <- aa * temp1 + temp3 - temp3a - temp4
    dl.dq <- temp3 - temp3b - temp4

    dscale.deta <- dtheta.deta(scale, .lscale , earg = .escale )
    dp.deta <- dtheta.deta(parg, .lshape2.p, earg = .eshape2.p )
    dq.deta <- dtheta.deta(qq, .lshape3.q, earg = .eshape3.q)

    c(w) * cbind( dl.dscale * dscale.deta,
                  dl.dp * dp.deta,
                  dl.dq * dq.deta )
  }), list( .lscale = lscale,
            .escale = escale, 
            .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
            .lshape2.p = lshape2.p, .lshape3.q = lshape3.q ))),
  weight = eval(substitute(expression({
    temp5  <- trigamma(parg + qq)
    ned2l.dscale <- aa^2 * parg * qq / (scale^2 * (1+parg+qq))
    ned2l.dp <- trigamma(parg) - temp5
    ned2l.dq <- trigamma(qq) - temp5
    ned2l.dscalep <-  aa * qq   / (scale*(parg+qq))
    ned2l.dscaleq <- -aa * parg / (scale*(parg+qq))
    ned2l.dpq <- -temp5

    wz <- matrix(as.numeric(NA), n, dimm(M))  #M==3 means 6=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.dscale * dscale.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dp * dp.deta^2
    wz[, iam(3, 3, M)] <- ned2l.dq * dq.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dscalep * dscale.deta * dp.deta
    wz[, iam(1, 3, M)] <- ned2l.dscaleq * dscale.deta * dq.deta
    wz[, iam(2, 3, M)] <- ned2l.dpq * dp.deta * dq.deta
    wz <- c(w) * wz
    wz
  }), list(
    .lscale = lscale, .lshape2.p = lshape2.p, .lshape3.q = lshape3.q,
    .escale = escale, .eshape2.p = eshape2.p, .eshape3.q = eshape3.q ))))
}



 lomax <- function(  # lss,
                   lscale = "loge",    lshape3.q = "loge",
                   iscale = NULL,      ishape3.q = NULL,  # 2.0, 
                                       gshape3.q = exp(-5:5),
                   zero = NULL) {


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")



  lshape3.q <- as.list(substitute(lshape3.q))
  eshape3.q <- link2list(lshape3.q)
  lshape3.q <- attr(eshape3.q, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  new("vglmff",
  blurb = c("Lomax distribution\n\n",
          "Links:    ",
          namesof("scale",    lscale,    earg = escale), ", ", 
          namesof("shape3.q", lshape3.q, earg = eshape3.q), "\n", 
          "Mean:     scale / (shape3.q - 1)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
        stop("response must be a vector or a one-column matrix")

    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1)



    mynames1 <- "scale"
    mynames2 <- "shape3.q"
    predictors.names <-
      c(namesof(mynames1,  .lscale    , earg = .escale ,    tag = FALSE),
        namesof(mynames2,  .lshape3.q , earg = .eshape3.q , tag = FALSE))

    aa <- parg <- 1



                   
    lomax.Loglikfun <- function(shape3.q, y, x, w, extraargs) {
      qvec <- c(0.25, 0.5, 0.75)  # Arbitrary; could be made an argument
      xvec <- log( (1-qvec)^(-1/ shape3.q ) - 1 )
      fit0 <- lsfit(x = xvec, y = log(quantile(y, qvec )))
      init.scale <- exp(fit0$coef[1])

      ans <- sum(c(w) * dlomax(x = y, shape3.q = shape3.q,
                               scale = init.scale, log = TRUE))
      ans
    }

    shape3.q.grid <- .gshape3.q
    yvec <- y
    wvec <- w
    Init.shape3.q <-
        grid.search(shape3.q.grid,
                    objfun = lomax.Loglikfun,
                    y = yvec, x = x, w = wvec,
                    extraargs = NULL)




    if (!length( .iscale )) {
      qvec <- c(0.25, 0.5, 0.75)  # Arbitrary; could be made an argument
      ishape3.q <- if (length( .ishape3.q )) .ishape3.q else Init.shape3.q
      xvec <- log( (1-qvec)^(-1/ ishape3.q ) - 1 )
      fit0 <- lsfit(x = xvec, y = log(quantile(y, qvec )))
    }

    if (!length(etastart)) {
      qq    <- rep(if (length( .ishape3.q)) .ishape3.q else Init.shape3.q,
                   length.out = n)
      scale <- rep(if (length( .iscale )) .iscale else
                   exp(fit0$coef[1]),
                   length.out = n)


      aa     <- 1
      outOfRange <- (qq   - 1/aa <= 0)
      qq[outOfRange] <- 1 / aa + 1



      etastart <-
        cbind(theta2eta(scale, .lscale   , earg = .escale ),
              theta2eta(qq,    .lshape3.q, earg = .eshape3.q))
    }
  }), list( .lscale = lscale, .lshape3.q = lshape3.q,
            .escale = escale, .eshape3.q = eshape3.q,
                              .gshape3.q = gshape3.q,
            .iscale = iscale, .ishape3.q = ishape3.q ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    aa     <- 1
    Scale  <- eta2theta(eta[, 1], .lscale    , earg = .escale )
    parg   <- 1
    qq     <- eta2theta(eta[, 2], .lshape3.q , earg = .eshape3.q )





    ans <- Scale * exp(lgamma(parg + 1/aa) +
                       lgamma(qq   - 1/aa) - lgamma(parg) - lgamma(qq))
    ans[parg + 1/aa <= 0] <- NA
    ans[qq   - 1/aa <= 0] <- NA
    ans[Scale       <= 0] <- NA
    ans[qq          <= 0] <- NA
    ans
  }, list( .lscale = lscale, .lshape3.q = lshape3.q,
           .escale = escale, .eshape3.q = eshape3.q ))),
  last = eval(substitute(expression({
    misc$link <-    c(scale = .lscale , shape3.q = .lshape3.q )

    misc$earg <- list(scale = .escale , shape3.q = .eshape3.q )


    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lscale = lscale, .lshape3.q = lshape3.q,
            .escale = escale, .eshape3.q = eshape3.q ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    aa <- 1
    scale <- eta2theta(eta[, 1], .lscale ,    earg = .escale )
    parg <- 1
    qq <-    eta2theta(eta[, 2], .lshape3.q,  earg = .eshape3.q )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dlomax(x = y, scale = scale,
                               shape3.q = qq, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lscale = lscale, .lshape3.q = lshape3.q,
           .escale = escale, .eshape3.q = eshape3.q ))),
  vfamily = c("lomax"),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    scale <- eta2theta(eta[, 1], .lscale ,    earg = .escale )
    qq <-    eta2theta(eta[, 2], .lshape3.q,  earg = .eshape3.q )
    rlomax(nsim * length(qq), scale = scale, shape3.q = qq)
  }, list( .lscale = lscale, .lshape3.q = lshape3.q,
           .escale = escale, .eshape3.q = eshape3.q ))),


  deriv = eval(substitute(expression({
    aa <- 1
    scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    parg <- 1
    qq <- eta2theta(eta[, 2], .lshape3.q, earg = .eshape3.q)
    temp2 <- (y/scale)^aa

    dl.dscale <- (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
    dl.dq <- digamma(parg + qq) - digamma(qq) - log1p(temp2)

    dscale.deta <- dtheta.deta(scale, .lscale ,   earg = .escale )
    dq.deta     <- dtheta.deta(qq,    .lshape3.q, earg = .eshape3.q)

    c(w) * cbind( dl.dscale * dscale.deta,
                  dl.dq * dq.deta )
  }), list( .lscale = lscale, .lshape3.q = lshape3.q,
            .escale = escale, .eshape3.q = eshape3.q ))),
  weight = eval(substitute(expression({
    ned2l.dscale <- aa^2 * parg * qq / (scale^2 * (1+parg+qq))
    ned2l.dq <- 1/qq^2 
    ned2l.dscaleq <- -aa * parg / (scale*(parg+qq))

    wz <- matrix(as.numeric(NA), n, dimm(M))  #M == 2 means 3=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.dscale * dscale.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dq * dq.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dscaleq * dscale.deta * dq.deta
    wz <- c(w) * wz
    wz
  }), list( .lscale = lscale, .lshape3.q = lshape3.q,
            .escale = escale, .eshape3.q = eshape3.q ))))
}



 fisk <- function(lss,
                  lshape1.a = "loge", lscale = "loge",
                  ishape1.a = NULL,   iscale = NULL,
                  zero = NULL) {

  if (!is.logical(lss) || lss)
    stop("argument 'lss' not specified correctly. ",
         "See online help for important information")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  lshape1.a <- as.list(substitute(lshape1.a))
  eshape1.a <- link2list(lshape1.a)
  lshape1.a <- attr(eshape1.a, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  new("vglmff",
  blurb = c("Fisk distribution\n\n",
            "Links:    ",
            namesof("shape1.a", lshape1.a, earg = eshape1.a), ", ", 
            namesof("scale",    lscale,    earg = escale), "\n", 
            "Mean:     scale * gamma(1 + 1/shape1.a) * ",
                      "gamma(1 - 1/shape1.a)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1)


    predictors.names <-
      c(namesof("shape1.a", .lshape1.a , earg = .eshape1.a , tag = FALSE),
        namesof("scale",    .lscale    , earg = .escale    , tag = FALSE))

    qq <- parg <- 1

    if (!length( .iscale )) {
      qvec <- c( 0.25, 0.5, 0.75)  # Arbitrary; could be made an argument
      xvec <- log( 1 / qvec - 1 )
      fit0 <- lsfit(x = xvec, y = log(quantile(y, qvec )))
    }

    if (!length(etastart)) {
      aa <- rep(if (length( .ishape1.a)) .ishape1.a else
                abs(-1 / fit0$coef[2]),
                length.out = n)
      scale <- rep(if (length( .iscale )) .iscale else
                   exp(fit0$coef[1]),
                   length.out = n)


      parg   <- 1
      qq     <- 1
      outOfRange <- (parg + 1/aa <= 0)
        parg[outOfRange] <- 1 / aa[outOfRange] + 1
        outOfRange <- (qq   - 1/aa <= 0)
      qq[outOfRange] <- 1 / aa + 1



      etastart <-
        cbind(theta2eta(aa,    .lshape1.a , earg = .eshape1.a ),
              theta2eta(scale, .lscale    , earg = .escale ))
    }
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale, 
            .ishape1.a = ishape1.a, .iscale = iscale ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    aa     <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale  <- eta2theta(eta[, 2], .lscale    , earg = .escale )
    parg   <- 1
    qq     <- 1

    ans <- Scale * exp(lgamma(parg + 1/aa) +
                       lgamma(qq   - 1/aa) - lgamma(parg) - lgamma(qq))
    ans[parg + 1/aa <= 0] <- NA
    ans[qq   - 1/aa <= 0] <- NA
    ans[aa          <= 0] <- NA
    ans[Scale       <= 0] <- NA
    ans
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1.a = .lshape1.a , scale = .lscale )

    misc$earg <- list(shape1.a = .eshape1.a , scale = .escale )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale <- eta2theta(eta[, 2], .lscale    , earg = .escale )
    parg <- qq <- 1
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dfisk(x = y, shape1.a = aa, scale = Scale, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale))),
  vfamily = c("fisk"),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    rfisk(nsim * length(aa), shape1.a = aa, scale = Scale)
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale ))),


  deriv = eval(substitute(expression({
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    scale <- eta2theta(eta[, 2], .lscale    , earg = .escale )
    parg <- qq <- 1

    temp1 <- log(y/scale)
    temp2 <- (y/scale)^aa
    temp3a <- digamma(parg)
    temp3b <- digamma(qq)

    dl.da <- 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
    dl.dscale <- (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))

    da.deta     <- dtheta.deta(aa,   .lshape1.a , earg = .eshape1.a )
    dscale.deta <- dtheta.deta(scale, .lscale   , earg = .escale )

    c(w) * cbind( dl.da * da.deta,
                  dl.dscale * dscale.deta )
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale))),
    weight = eval(substitute(expression({
    ned2l.da <- (1 + parg+qq + parg * qq * (trigamma(parg) + trigamma(qq) +
                (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1+parg+qq))
    ned2l.dscale <- aa^2 * parg * qq / (scale^2 * (1 + parg + qq))
    ned2l.dascale <- (parg - qq - parg*qq*(temp3a - temp3b)) / (
                      scale * (1 + parg + qq))

    wz <- matrix(as.numeric(NA), n, dimm(M))  #M == 2 means 3=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.da * da.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dscale * dscale.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dascale * da.deta * dscale.deta

    wz <- c(w) * wz
    wz
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale))))
}


 inv.lomax <- function(  # lss,
                       lscale = "loge",
                       lshape2.p = "loge",
                       iscale = NULL,
                       ishape2.p = 1.0, 
                       zero = NULL) {

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  lshape2.p <- as.list(substitute(lshape2.p))
  eshape2.p <- link2list(lshape2.p)
  lshape2.p <- attr(eshape2.p, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")




  new("vglmff",
  blurb = c("Inverse Lomax distribution\n\n",
          "Links:    ",
          namesof("scale",    lscale,    earg = escale), ", ", 
          namesof("shape2.p", lshape2.p, earg = eshape2.p), "\n", 
          "Mean:     does not exist"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1)



    predictors.names <-
      c(namesof("scale",    .lscale ,   earg = .escale    , tag = FALSE),
        namesof("shape2.p", .lshape2.p, earg = .eshape2.p , tag = FALSE))

    qq <- aa <- 1

    if (!length( .iscale )) {
      qvec <- c(0.25, .5, .75)  # Arbitrary; could be made an argument
      ishape2.p <- if (length( .ishape2.p)) .ishape2.p else 1
      xvec <- log( qvec^(-1/ ishape2.p ) - 1 )
      fit0 <- lsfit(x = xvec, y = log(quantile(y, qvec )))
    }
        if (!length(etastart)) {
          scale <- rep(if (length( .iscale )) .iscale else
                       exp(fit0$coef[1]),
                       length.out = n)
          parg <- rep(if (length( .ishape2.p)) .ishape2.p else 1.0,
                      length.out = n)




          etastart <-
            cbind(theta2eta(scale, .lscale ,   earg = .escale ),
                  theta2eta(parg,  .lshape2.p, earg = .eshape2.p ))
        }
    }), list( .lscale = lscale, .lshape2.p = lshape2.p,
              .escale = escale, .eshape2.p = eshape2.p,
              .iscale = iscale, .ishape2.p = ishape2.p ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale  <- eta2theta(eta[, 1], .lscale ,   earg = .escale )
    parg   <- eta2theta(eta[, 2], .lshape2.p, earg = .eshape2.p )



    qinv.lomax(p = 0.5, scale = Scale, shape2.p = parg)
    }, list( .lscale = lscale,
             .escale = escale, 
             .eshape2.p = eshape2.p,
             .lshape2.p = lshape2.p ))),
  last = eval(substitute(expression({
    misc$link <-    c(scale = .lscale , shape2.p = .lshape2.p )

    misc$earg <- list(scale = .escale , shape2.p = .eshape2.p )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lscale = lscale,
            .escale = escale, 
            .eshape2.p = eshape2.p,
            .lshape2.p = lshape2.p ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    aa <- 1
    Scale <- eta2theta(eta[, 1], .lscale    , earg = .escale    )
    parg  <- eta2theta(eta[, 2], .lshape2.p , earg = .eshape2.p )
    qq <- 1
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dinv.lomax(x = y, scale = Scale,
                                  shape2.p = parg, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lscale = lscale, .lshape2.p = lshape2.p,
           .escale = escale, .eshape2.p = eshape2.p ))),
  vfamily = c("inv.lomax"),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    Scale <- eta2theta(eta[, 1], .lscale    , earg = .escale    )
    parg  <- eta2theta(eta[, 2], .lshape2.p , earg = .eshape2.p )
    rinv.lomax(nsim * length(Scale), scale = Scale, shape2.p = parg)
  }, list( .lscale = lscale, .lshape2.p = lshape2.p,
           .escale = escale, .eshape2.p = eshape2.p ))),


  deriv = eval(substitute(expression({
    aa <- qq <- 1 
    scale <- eta2theta(eta[, 1], .lscale    , earg = .escale    )
    parg  <- eta2theta(eta[, 2], .lshape2.p , earg = .eshape2.p )

    temp1 <- log(y/scale)
    temp2 <- (y/scale)^aa

    dl.dscale <- (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
    dl.dp <- aa * temp1 + digamma(parg + qq) - digamma(parg) - log1p(temp2)

    dscale.deta <- dtheta.deta(scale, .lscale ,    earg = .escale )
    dp.deta     <- dtheta.deta(parg,  .lshape2.p , earg = .eshape2.p )

    c(w) * cbind( dl.dscale * dscale.deta,
                  dl.dp     * dp.deta )
  }), list( .lscale = lscale, .lshape2.p = lshape2.p,
            .escale = escale, .eshape2.p = eshape2.p ))),
  weight = eval(substitute(expression({
    ned2l.dscale <- aa^2 * parg * qq / (scale^2 * (1 + parg + qq))
    ned2l.dp <- 1 / parg^2 
    ned2l.dscalep <-  aa * qq   / (scale * (parg + qq))

    wz <- matrix(as.numeric(NA), n, dimm(M))  #M == 2 means 3=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.dscale * dscale.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dp * dp.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dscalep * dscale.deta * dp.deta
    wz <- c(w) * wz
    wz
  }), list( .lscale = lscale, .lshape2.p = lshape2.p,
            .escale = escale, .eshape2.p = eshape2.p ))))
}


 paralogistic <- function(lss,
                          lshape1.a = "loge",
                          lscale = "loge",
                          ishape1.a = 2,
                          iscale = NULL,
                          zero = NULL) {

  if (!is.logical(lss) || lss)
    stop("argument 'lss' not specified correctly. ",
         "See online help for important information")




  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  lshape1.a <- as.list(substitute(lshape1.a))
  eshape1.a <- link2list(lshape1.a)
  lshape1.a <- attr(eshape1.a, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  new("vglmff",
  blurb = c("Paralogistic distribution\n\n",
          "Links:    ",
          namesof("shape1.a", lshape1.a, earg = eshape1.a), ", ", 
          namesof("scale",    lscale,    earg = escale), "\n", 
          "Mean:     scale * gamma(1 + 1/shape1.a) * ",
                    "gamma(shape1.a - 1/shape1.a) / gamma(shape1.a)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1)



    predictors.names <-
      c(namesof("shape1.a", .lshape1.a , earg = .eshape1.a , tag = FALSE),
        namesof("scale",    .lscale ,    earg = .escale ,    tag = FALSE))

    parg <- 1

    if (!length( .ishape1.a) || !length( .iscale )) {
      qvec <- c( .25, .5, .75)  # Arbitrary; could be made an argument
      ishape1.a <- if (length( .ishape1.a)) .ishape1.a else 1
      xvec <- log( (1-qvec)^(-1/ ishape1.a ) - 1 )
      fit0 <- lsfit(x = xvec, y = log(quantile(y, qvec )))
    }

    if (!length(etastart)) {
      aa <- rep(if (length( .ishape1.a)) .ishape1.a else
                abs(1/fit0$coef[2]),
                length.out = n)
      scale <- rep(if (length( .iscale )) .iscale else
                   exp(fit0$coef[1]),
                   length.out = n)




      parg   <- 1
      qq     <- aa
      outOfRange <- (parg + 1/aa <= 0)
      parg[outOfRange] <- 1 / aa[outOfRange] + 1
      outOfRange <- (qq   - 1/aa <= 0)
      aa[outOfRange] <-
      qq[outOfRange] <- 2 # Need aa > 1, where aa == qq


        etastart <-
          cbind(theta2eta(aa,    .lshape1.a , earg = .eshape1.a ),
                theta2eta(scale, .lscale ,    earg = .escale ))
      }
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale, 
            .ishape1.a = ishape1.a, .iscale = iscale
            ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    aa     <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale  <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    parg   <- 1
    qq     <- aa

    ans <- Scale * exp(lgamma(parg + 1/aa) +
                       lgamma(qq   - 1/aa) - lgamma(parg) - lgamma(qq))
    ans[parg + 1/aa <= 0] <- NA
    ans[qq   - 1/aa <= 0] <- NA
    ans[aa          <= 0] <- NA
    ans[Scale       <= 0] <- NA
    ans
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1.a = .lshape1.a , scale = .lscale )

    misc$earg <- list(shape1.a = .eshape1.a , scale = .escale )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    parg <- 1
    qq <- aa
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dparalogistic(x = y, shape1.a = aa,
                                      scale = Scale, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale))),
  vfamily = c("paralogistic"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale <- eta2theta(eta[, 2], .lscale  ,   earg = .escale )
    rparalogistic(nsim * length(Scale), shape1.a = aa, scale = Scale)
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale))),



  deriv = eval(substitute(expression({
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    scale <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    parg <- 1
    qq <- aa

    temp1 <- log(y/scale)
    temp2 <- (y/scale)^aa
    temp3a <- digamma(parg)
    temp3b <- digamma(qq)

    dl.da     <- 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
    dl.dscale <- (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))

    da.deta     <- dtheta.deta(aa,    .lshape1.a , earg = .eshape1.a )
    dscale.deta <- dtheta.deta(scale, .lscale ,    earg = .escale )

    c(w) * cbind( dl.da * da.deta,
               dl.dscale * dscale.deta)
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale))),
  weight = eval(substitute(expression({
    ned2l.da <- (1 + parg + qq +
                 parg * qq * (trigamma(parg) + trigamma(qq) +
                (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1+parg+qq))
    ned2l.dscale <- aa^2 * parg * qq / (scale^2 * (1+parg+qq))
    ned2l.dascale <- (parg - qq - parg*qq*(temp3a -temp3b)) /(
                      scale*(1 + parg+qq))

    wz <- matrix(as.numeric(NA), n, dimm(M))  #M == 2 means 3=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.da * da.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dscale * dscale.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dascale * da.deta * dscale.deta
    wz <- c(w) * wz
    wz
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale))))
}



 inv.paralogistic <- function(lss,
                              lshape1.a = "loge", lscale = "loge",
                              ishape1.a = 2,      iscale = NULL,
                              zero = NULL) {

  if (!is.logical(lss) || lss)
    stop("argument 'lss' not specified correctly. ",
         "See online help for important information")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  lshape1.a <- as.list(substitute(lshape1.a))
  eshape1.a <- link2list(lshape1.a)
  lshape1.a <- attr(eshape1.a, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")


  new("vglmff",
  blurb = c("Inverse paralogistic distribution\n\n",
            "Links:    ",
            namesof("shape1.a", lshape1.a, earg = eshape1.a), ", ", 
            namesof("scale",    lscale,    earg = escale), "\n", 
            "Mean:     scale * gamma(shape1.a + 1/shape1.a) * ",
                      "gamma(1 - 1/shape1.a)/gamma(shape1.a)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1)


    predictors.names <-
      c(namesof("shape1.a", .lshape1.a ,  earg = .eshape1.a , tag = FALSE),
        namesof("scale",    .lscale    ,  earg = .escale ,    tag = FALSE))

    if (!length( .ishape1.a) || !length( .iscale )) {
      qvec <- c( .25, .5, .75)  # Arbitrary; could be made an argument
      ishape2.p <- if (length( .ishape1.a )) .ishape1.a else 1
      xvec <- log( qvec^(-1/ ishape2.p ) - 1 )
      fit0 <- lsfit(x = xvec, y = log(quantile(y, qvec )))
    }

    qq <- 1
    if (!length(etastart)) {
      aa <- rep(if (length( .ishape1.a)) .ishape1.a else -1/fit0$coef[2],
                length = n)
      scale <- rep(if (length( .iscale )) .iscale else
                  exp(fit0$coef[1]), length = n)





    parg <- aa
    qq <- 1
      outOfRange <- (parg + 1/aa <= 0)
      parg[outOfRange] <-
        aa[outOfRange] <- 2
      outOfRange <- (qq   - 1/aa <= 0)
      qq[outOfRange] <- 1 / aa[outOfRange] + 1



      etastart <-
        cbind(theta2eta(aa,    .lshape1.a , earg = .eshape1.a ),
              theta2eta(scale, .lscale ,    earg = .escale ))
    }
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale,
            .ishape1.a = ishape1.a, .iscale = iscale ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    aa     <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale  <- eta2theta(eta[, 2], .lscale    , earg = .escale )
    parg <- aa
    qq <- 1

    ans <- Scale * exp(lgamma(parg + 1/aa) +
                       lgamma(qq   - 1/aa) - lgamma(parg) - lgamma(qq))
    ans[parg + 1/aa <= 0] <- NA
    ans[qq   - 1/aa <= 0] <- NA
    ans[aa          <= 0] <- NA
    ans[Scale       <= 0] <- NA
    ans
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1.a = .lshape1.a , scale = .lscale )

    misc$earg <- list(shape1.a = .eshape1.a , scale = .escale )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    parg <- aa
    qq <- 1
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dinv.paralogistic(x = y, shape1.a = aa,
                                         scale = Scale, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale))),
  vfamily = c("inv.paralogistic"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    Scale <- eta2theta(eta[, 2], .lscale  ,   earg = .escale )
    rinv.paralogistic(nsim * length(Scale), shape1.a = aa, scale = Scale)
  }, list( .lshape1.a = lshape1.a, .lscale = lscale,
           .eshape1.a = eshape1.a, .escale = escale))),



  deriv = eval(substitute(expression({
    aa    <- eta2theta(eta[, 1], .lshape1.a , earg = .eshape1.a )
    scale <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
    parg <- aa 
    qq <- 1

    temp1 <- log(y/scale)
    temp2 <- (y/scale)^aa
    temp3a <- digamma(parg)
    temp3b <- digamma(qq)

    dl.da <- 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
    dl.dscale <- (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))

    da.deta     <- dtheta.deta(aa,    .lshape1.a , earg = .eshape1.a )
    dscale.deta <- dtheta.deta(scale, .lscale ,    earg = .escale )

    c(w) * cbind( dl.da     * da.deta,
                  dl.dscale * dscale.deta )
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale))),

  weight = eval(substitute(expression({
    ned2l.da <- (1 + parg + qq +
                 parg * qq * (trigamma(parg) + trigamma(qq) +
                (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1 + parg + qq))
    ned2l.dscale <- aa^2 * parg * qq / (scale^2 * (1 + parg + qq))
    ned2l.dascale <- (parg - qq -
                      parg * qq * (temp3a -temp3b)) / (scale *
                     (1 + parg + qq))

    wz <- matrix(as.numeric(NA), n, dimm(M))  # M==3 means 6=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.da * da.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dscale * dscale.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dascale * da.deta * dscale.deta
    wz <- c(w) * wz
    wz
  }), list( .lshape1.a = lshape1.a, .lscale = lscale,
            .eshape1.a = eshape1.a, .escale = escale))))
}









 if (FALSE)
 genlognormal <- function(link.sigma = "loge", link.r = "loge",
                          init.sigma = 1, init.r = 1, zero = NULL) {
warning("20040402; does not work, possibly because first derivs are ",
        "not continuous (sign() is used). Certainly, the derivs wrt ",
        "mymu are problematic (run with maxit=4:9 and look at weight ",
        "matrices). Possibly fundamentally cannot be estimated by IRLS. ",
        "Pooling does not seem to help")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")



  link.sigma <- as.list(substitute(link.sigma))
  esigma <- link2list(link.sigma)
  link.sigma <- attr(esigma, "function.name")

  link.r <- as.list(substitute(link.r))
  er <- link2list(link.r)
  link.r <- attr(er, "function.name")


  new("vglmff",
  blurb = c("Three-parameter generalized lognormal distribution\n\n",
          "Links:    ",
          "loc; ",
          namesof("sigma", link.sigma, earg = esigma, tag = TRUE), ", ",
          namesof("r",     link.r,     earg = er,     tag = TRUE)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")
    predictors.names <-
    c(namesof("loc", "identitylink", earg = list(),     tag = FALSE),
      namesof("sigma", .link.sigma, earg = .esigma, tag = FALSE),
      namesof("r",     .link.r,     earg = .er,     tag = FALSE))

    if (!length( .init.sigma) || !length( .init.r)) {
      init.r <- if (length( .init.r)) .init.r else 1
      sigma.init <- (0.5 * sum(abs(log(y) -
                               mean(log(y )))^init.r))^(1/init.r)
    }
    if (any(y <= 0))
      stop("y must be positive")

    if (!length(etastart)) {
      sigma.init <- rep(if (length( .init.sigma)) .init.sigma else
                        sigma.init, length.out = n)
      r.init <- if (length( .init.r)) .init.r else init.r
      etastart <-
        cbind(mu = rep(log(median(y)), length.out = n),
              sigma = sigma.init,
              r = r.init)
    }
  }), list( .link.sigma = link.sigma, .link.r = link.r,
            .init.sigma = init.sigma, .init.r = init.r ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    mymu  <- eta2theta(eta[, 1], "identitylink" , earg = list())
    sigma <- eta2theta(eta[, 2], .link.sigma , earg = .esigma )
    rrrrr <- eta2theta(eta[, 3], .link.r , earg = .er )
    rrrrr
  }, list( .link.sigma = link.sigma, .link.r = link.r ))),

  last = eval(substitute(expression({
    misc$link = c(loc = "identitylink",
                  "sigma" = .link.sigma,
                  r = .link.r )
    misc$expected = TRUE
  }), list( .link.sigma = link.sigma, .link.r = link.r ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mymu <- eta2theta(eta[, 1], "identitylink", earg = list())
    sigma <- eta2theta(eta[, 2], .link.sigma, earg = .esigma)
    r <- eta2theta(eta[, 3], .link.r, earg = .er)
    temp89 <- (abs(log(y)-mymu)/sigma)^r
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * (-log(r^(1/r) * sigma) - lgamma(1+1/r) - temp89/r)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link.sigma = link.sigma, .link.r = link.r ))),
  vfamily = c("genlognormal3"),
  deriv = eval(substitute(expression({
    mymu  <- eta2theta(eta[, 1], "identitylink", earg = list())
    sigma <- eta2theta(eta[, 2], .link.sigma, earg = .esigma)

    r <- eta2theta(eta[, 3], .link.r, earg = .er)
    ss <- 1 + 1/r
    temp33 <- (abs(log(y)-mymu)/sigma)
    temp33r1 <- temp33^(r-1)

    dl.dmymu <- temp33r1 * sign(log(y)-mymu) / sigma
    dl.dsigma <- (temp33*temp33r1 - 1) / sigma
    dl.dr <- (log(r) - 1 + digamma(ss) + temp33*temp33r1)/r^2 -
             temp33r1 * log(temp33r1) / r

    dmymu.deta <- dtheta.deta(mymu, "identitylink", earg = list())
    dsigma.deta <- dtheta.deta(sigma, .link.sigma, earg = .esigma)
    dr.deta <- dtheta.deta(r, .link.r, earg = .er)

    c(w) * cbind(dl.dmymu * dmymu.deta, 
                 dl.dsigma * dsigma.deta, 
                 dl.dr * dr.deta)
  }), list( .link.sigma = link.sigma, .link.r = link.r ))),
  weight = expression({
    wz <- matrix(0, n, 6)  # 5 will have small savings of 1 column

    B <- log(r) + digamma(ss)
    ned2l.dmymu2 <- (r-1) * gamma(1-1/r) / (sigma^2 * r^(2/r) * gamma(ss))
    ned2l.dsigma2 <- r / sigma^2
    ned2l.dr2 <- (ss * trigamma(ss) + B^2 - 1) / r^3 
    ned2l.dsigmar <- -B / (r * sigma)

    wz[, iam(1, 1, M)] <- ned2l.dmymu2 * dmymu.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dsigma2 * dsigma.deta^2
    wz[, iam(3, 3, M)] <- ned2l.dr2 * dr.deta^2
    wz[, iam(2, 3, M)] <- ned2l.dsigmar * dsigma.deta * dr.deta
    wz = c(w) * wz
    wz
  }))
}








# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.













dgumbelII <- function(x, shape, scale = 1, log = FALSE) {


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


pgumbelII <- function(q, shape, scale = 1) {

  LLL <- max(length(q), length(shape), length(scale))
  if (length(q)       != LLL) q       <- rep(q,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)


  ans <- exp(-(q / scale)^(-shape))
  ans[(q <= 0)] <- 0
  ans[shape <= 0 | scale <= 0] <- NaN
  ans[q == Inf] <- 1
  ans
}



qgumbelII <- function(p, shape, scale = 1) {

  LLL <- max(length(p), length(shape), length(scale))
  if (length(p)       != LLL) p       <- rep(p,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)


  ans <- scale * (-log(p))^(-1 / shape)
  ans[p < 0] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans[p > 1] <- NaN
  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}


rgumbelII <- function(n, shape, scale = 1) {
  qgumbelII(runif(n), shape = shape, scale = scale)
}









 gumbelII <-
  function(lshape = "loge", lscale = "loge",
           ishape = NULL,   iscale = NULL,
           probs.y = c(0.2, 0.5, 0.8),
           perc.out = NULL, # 50,
           imethod = 1, zero = -2)
{


  lshape <- as.list(substitute(lshape))
  e.shape <- link2list(lshape)
  l.shape <- attr(e.shape, "function.name")

  lscale <- as.list(substitute(lscale))
  e.scale <- link2list(lscale)
  l.scale <- attr(e.scale, "function.name")


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE))
    stop("bad input for argument 'zero'")
  if (!is.Numeric(imethod, allowable.length = 1,
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
            namesof("shape", l.shape, e.shape), ", ",
            namesof("scale", l.scale, e.scale), "\n",
            "Mean:     scale^(1/shape) * gamma(1 - 1 / shape)\n",
            "Variance: scale^(2/shape) * (gamma(1 - 2/shape) - ",
                      "gamma(1 + 1/shape)^2)"),
 constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
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
    Musual <- 2
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    mynames1 <- paste("shape", if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("scale", if (ncoly > 1) 1:ncoly else "", sep = "")


    predictors.names <-
        c(namesof(mynames1, .l.shape , .e.shape , tag = FALSE),
          namesof(mynames2, .l.scale , .e.scale , tag = FALSE))[
          interleave.VGAM(M, M = Musual)]


    Shape.init <- matrix(if(length( .ishape )) .ishape else 0 + NA,
                         n, ncoly, byrow = TRUE)
    Scale.init <- matrix(if(length( .iscale )) .iscale else 0 + NA,
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
        } # ilocal

        etastart <-
          cbind(theta2eta(Shape.init, .l.shape , .e.shape ),
                theta2eta(Scale.init, .l.scale , .e.scale ))[,
                interleave.VGAM(M, M = Musual)]
      }
    }
  }), list(
            .l.scale = l.scale, .l.shape = l.shape,
            .e.scale = e.scale, .e.shape = e.shape,
            .iscale = iscale, .ishape = ishape,
            .probs.y = probs.y,
            .imethod = imethod ) )),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Shape <- eta2theta(eta[, c(TRUE, FALSE)], .l.shape , .e.shape )
    Scale <- eta2theta(eta[, c(FALSE, TRUE)], .l.scale , .e.scale )
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
           .l.scale = l.scale, .l.shape = l.shape,
           .e.scale = e.scale, .e.shape = e.shape,
           .perc.out = perc.out ) )),
  last = eval(substitute(expression({


    Musual <- extra$Musual
    misc$link <-
      c(rep( .l.shape , length = ncoly),
        rep( .l.scale , length = ncoly))[interleave.VGAM(M, M = Musual)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = Musual)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for(ii in 1:ncoly) {
      misc$earg[[Musual*ii-1]] <- .e.shape
      misc$earg[[Musual*ii  ]] <- .e.scale
    }

    misc$Musual <- Musual
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$perc.out <- .perc.out
    misc$true.mu <- FALSE # @fitted is not a true mu


  }), list(
            .l.scale = l.scale, .l.shape = l.shape,
            .e.scale = e.scale, .e.shape = e.shape,
            .perc.out = perc.out,
            .imethod = imethod ) )),
  loglikelihood = eval(substitute(
   function(mu, y, w, residuals = FALSE,eta, extra = NULL) {
    Shape <- eta2theta(eta[, c(TRUE, FALSE)], .l.shape , .e.shape )
    Scale <- eta2theta(eta[, c(FALSE, TRUE)], .l.scale , .e.scale )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
      sum(c(w) * dgumbelII(x = y, shape = Shape,
                           scale = Scale, log = TRUE))
  }, list( .l.scale = l.scale, .l.shape = l.shape,
           .e.scale = e.scale, .e.shape = e.shape
         ) )),
  vfamily = c("gumbelII"),
  deriv = eval(substitute(expression({
    Musual <- 2
    Shape <- eta2theta(eta[, c(TRUE, FALSE)], .l.shape , .e.shape )
    Scale <- eta2theta(eta[, c(FALSE, TRUE)], .l.scale , .e.scale )

    dl.dshape <- 1 / Shape + log(Scale / y) -
                 log(Scale / y) * (Scale / y)^Shape
    dl.dscale <- Shape / Scale - (Shape / y) * (Scale / y)^(Shape - 1)


    dshape.deta <- dtheta.deta(Shape, .l.shape , .e.shape )
    dscale.deta <- dtheta.deta(Scale, .l.scale , .e.scale )

    myderiv <- c(w) * cbind(dl.dshape, dl.dscale) *
                      cbind(dshape.deta, dscale.deta)
    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list( .l.scale = l.scale, .l.shape = l.shape,
            .e.scale = e.scale, .e.shape = e.shape
          ) )),
  weight = eval(substitute(expression({
    EulerM <- -digamma(1.0)


    ned2l.dshape2 <- (1 + trigamma(2) + digamma(2)^2) / Shape^2
    ned2l.dscale2 <-  (Shape / Scale)^2
    ned2l.dshapescale <- digamma(2) / Scale

    wz <- matrix(0.0, n, M + M - 1) # wz is tridiagonal

    ind11 <- ind22 <- ind12 <- NULL
    for (ii in 1:(M / Musual)) {
      ind11 <- c(ind11, iam(Musual*ii - 1, Musual*ii - 1, M))
      ind22 <- c(ind22, iam(Musual*ii - 0, Musual*ii - 0, M))
      ind12 <- c(ind12, iam(Musual*ii - 1, Musual*ii - 0, M))
    }
    wz[, ind11] <- ned2l.dshape2 * dshape.deta^2
    wz[, ind22] <- ned2l.dscale2 * dscale.deta^2
    wz[, ind12] <- ned2l.dshapescale * dscale.deta * dshape.deta

    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / Musual)
  }), list( .l.scale = l.scale, .l.shape = l.shape ))))
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


  index0 = (x < 0)

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







dmperks <- function(x, shape, scale = 1, epsilon, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(scale), length(epsilon))
  if (length(x)       != LLL) x       <- rep(x,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  if (length(epsilon) != LLL) epsilon <- rep(epsilon, length.out = LLL)


  index0 = (x < 0)
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



pmperks <- function(q, shape, scale = 1, epsilon) {

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

  index0 = (x < 0)
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
alpha=shape;  beta=scale;

 warning("does not integrate to unity")

  ret=ifelse(x<=0 | beta<=0,NaN,
      exp(alpha+beta*x)*(1+exp(alpha+rho))**(exp(-rho/beta))/
      (1+exp(alpha+rho+beta*x))**(1+exp(-rho/beta)))
  ret
}



qbeard=function(x,u=0.5,alpha=1,beta=1,rho=1) {
  ret     = ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u) | beta<=0,
  NaN, (1/beta)*
  (log((u**(-beta*exp(rho)))*
  (1+exp(alpha+rho+beta*x))-1)-alpha-rho)-x)

  return(ret)
}










dperks <- function(x, shape, scale = 1, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(scale))
  if (length(x)     != LLL) x     <- rep(x,     length.out = LLL)
  if (length(shape) != LLL) shape <- rep(shape, length.out = LLL)
  if (length(scale) != LLL) scale <- rep(scale, length.out = LLL)

  index0 = (x < 0)
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



pperks <- function(q, shape, scale = 1) {

  LLL <- max(length(q), length(shape), length(scale))
  if (length(q)       != LLL) q       <- rep(q,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)

  logS <- -q + (log1p(shape) -
          log(shape + exp(-q * scale))) / scale
  ans <- -expm1(logS)

  ans[(q <= 0)] <- 0
  ans[shape <= 0 | scale <= 0] <- NaN
  ans[q == Inf] <- 1
  ans
}


qperks <- function(p, shape, scale = 1) {

  LLL <- max(length(p), length(shape), length(scale))
  if (length(p)       != LLL) p       <- rep(p,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)

  tmp <- scale * log1p(-p)
  onemFb <- exp(tmp)
  ans <- (log1p(shape - onemFb) - log(shape) - tmp) / scale
  ans[p < 0] <- NaN
  ans[p == 0] <- 0
  ans[p > 1] <- NaN
  ans[p == 1] <- Inf
  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}


rperks <- function(n, shape, scale = 1) {
  qperks(runif(n), shape = shape, scale = scale)
}





perks.control <- function(save.weight = TRUE, ...)
{
  list(save.weight = save.weight)
}


 perks <-
  function(lshape = "loge", lscale = "loge",
           ishape = NULL,   iscale = NULL,
           nsimEIM = 500,
           oim.mean = FALSE,
           zero = NULL)
{

  lshape <- as.list(substitute(lshape))
  e.shape <- link2list(lshape)
  l.shape <- attr(e.shape, "function.name")

  lscale <- as.list(substitute(lscale))
  e.scale <- link2list(lscale)
  l.scale <- attr(e.scale, "function.name")


  if (!is.Numeric(nsimEIM, allowable.length = 1,
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
            namesof("shape", l.shape, e.shape), ", ",
            namesof("scale", l.scale, e.scale), "\n",
            "Median:     qperks(p = 0.5, shape, scale)"),

  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
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
    Musual <- 2
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    mynames1 <- paste("shape", if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("scale", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .l.shape , .e.shape , tag = FALSE),
          namesof(mynames2, .l.scale , .e.scale , tag = FALSE))[
          interleave.VGAM(M, M = Musual)]



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

        perks.Loglikfun <- function(scaleval, y, x, w, extraargs) {
          ans <-
          sum(c(w) * dperks(x = y, shape = extraargs$Shape,
                            scale = scaleval, log = TRUE))
          ans
        }

        mymat <- matrix(-1, length(shape.grid), 2)
        for (jlocal in 1:length(shape.grid)) {
          mymat[jlocal, ] <-
            getMaxMin(scale.grid,
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
      } # spp.

      etastart <-
          cbind(theta2eta(matH, .l.shape , .e.shape ),
                theta2eta(matC, .l.scale , .e.scale ))[,
                interleave.VGAM(M, M = Musual)]
    } # End of !length(etastart)
  }), list( .l.scale = l.scale, .l.shape = l.shape,
            .e.scale = e.scale, .e.shape = e.shape,
            .ishape = ishape, .iscale = iscale
            ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Shape <- eta2theta(eta[, c(TRUE, FALSE)], .l.shape , .e.shape )
    Scale <- eta2theta(eta[, c(FALSE, TRUE)], .l.scale , .e.scale )

    qperks(p = 0.5, shape = Shape, scale = Scale)
  }, list( .l.scale = l.scale, .l.shape = l.shape,
           .e.scale = e.scale, .e.shape = e.shape ))),
  last = eval(substitute(expression({

    misc$link <-
      c(rep( .l.shape , length = ncoly),
        rep( .l.scale , length = ncoly))[interleave.VGAM(M, M = Musual)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = Musual)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for(ii in 1:ncoly) {
      misc$earg[[Musual*ii-1]] <- .e.shape
      misc$earg[[Musual*ii  ]] <- .e.scale
    }


    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list( .l.scale = l.scale, .l.shape = l.shape,
            .e.scale = e.scale, .e.shape = e.shape,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Shape <- eta2theta(eta[, c(TRUE, FALSE)], .l.shape , .e.shape )
    Scale <- eta2theta(eta[, c(FALSE, TRUE)], .l.scale , .e.scale )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dperks(x = y, shape = Shape,
                        scale = Scale, log = TRUE))
    }
  }, list( .l.scale = l.scale, .l.shape = l.shape,
           .e.scale = e.scale, .e.shape = e.shape ))),
  vfamily = c("perks"),
 
  deriv = eval(substitute(expression({
    Musual <- 2
    shape <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                       .l.shape , .e.shape )
    scale <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                       .l.scale , .e.scale )


    temp2 <- exp(y * scale)
    temp3 <- 1 + shape * temp2
    dl.dshape <- 1 / shape + 1 / (scale * (1 + shape)) -
                 (1 + 1 / scale) * temp2 / temp3
    dl.dscale <- y - log1p(shape) / scale^2 +
                 log1p(shape * temp2) / scale^2 -
                 (1 + 1 / scale) * shape * y * temp2 / temp3

    dshape.deta <- dtheta.deta(shape, .l.shape , .e.shape )
    dscale.deta <- dtheta.deta(scale, .l.scale , .e.scale )

    dthetas.detas <- cbind(dshape.deta, dscale.deta)
    myderiv <- c(w) * cbind(dl.dshape, dl.dscale) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list( .l.scale = l.scale, .l.shape = l.shape,
            .e.scale = e.scale, .e.shape = e.shape ))),


  weight = eval(substitute(expression({

    NOS <- M / Musual
    dThetas.detas <- dthetas.detas[, interleave.VGAM(M, M = Musual)]

    wz <- matrix(0.0, n, M + M - 1) # wz is 'tridiagonal' 

    ind1 <- iam(NA, NA, M = Musual, both = TRUE, diag = TRUE)


    for(spp. in 1:NOS) {
      run.varcov <- 0
      Shape <- shape[, spp.]
      Scale <- scale[, spp.]




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

      for(ii in 1:( .nsimEIM )) {
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


        temp7 <- cbind(dl.dshape, dl.dscale)
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

      wz1 <- wz1 * dThetas.detas[, Musual * (spp. - 1) + ind1$row] *
                   dThetas.detas[, Musual * (spp. - 1) + ind1$col]


      for(jay in 1:Musual)
        for(kay in jay:Musual) {
          cptr <- iam((spp. - 1) * Musual + jay,
                      (spp. - 1) * Musual + kay,
                      M = M)
          wz[, cptr] <- wz1[, iam(jay, kay, M = Musual)]
        }
    } # End of for(spp.) loop



    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / Musual)
  }), list( .l.scale = l.scale,
            .e.scale = e.scale,
            .nsimEIM = nsimEIM, .oim.mean = oim.mean ))))
} # perks()








dmakeham <- function(x, shape, scale = 1, epsilon = 0, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(scale), length(epsilon))
  if (length(x)       != LLL) x       <- rep(x,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  if (length(epsilon) != LLL) epsilon <- rep(epsilon, length.out = LLL)

  index0 = (x < 0)
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



pmakeham <- function(q, shape, scale = 1, epsilon = 0) {

  LLL <- max(length(q), length(shape), length(scale), length(epsilon))
  if (length(q)       != LLL) q       <- rep(q,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  if (length(epsilon) != LLL) epsilon <- rep(epsilon, length.out = LLL)


  ans <- -expm1(-q * epsilon - (shape / scale) * expm1(scale * q))
  ans[(q <= 0)] <- 0
  ans[shape <= 0 | scale <= 0 | epsilon < 0] <- NaN
  ans[q == Inf] <- 1
  ans
}



qmakeham <- function(p, shape, scale = 1, epsilon = 0) {

  LLL <- max(length(p), length(shape), length(scale), length(epsilon))
  if (length(p)       != LLL) p       <- rep(p,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)
  if (length(epsilon) != LLL) epsilon <- rep(epsilon, length.out = LLL)


  ans <- shape / (scale * epsilon) - log1p(-p) / epsilon -
  lambertW((shape / epsilon) * exp(shape / epsilon) *
          (1 - p)^(-(scale / epsilon))) / scale
  ans[epsilon == 0] <-
    qgompertz(p     =     p[epsilon == 0],
              shape = shape[epsilon == 0],
              scale = scale[epsilon == 0])
  ans[p < 0] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans[p > 1] <- NaN
  ans[shape <= 0 | scale <= 0 | epsilon < 0] <- NaN
  ans
}


rmakeham <- function(n, shape, scale = 1, epsilon = 0) {
  qmakeham(runif(n), shape = shape, scale = scale, epsilon = epsilon)
}




makeham.control <- function(save.weight = TRUE, ...)
{
  list(save.weight = save.weight)
}


 makeham <-
  function(lshape = "loge", lscale = "loge", lepsilon = "loge",
           ishape = NULL,   iscale = NULL,   iepsilon = 0.3,
           nsimEIM = 500,
           oim.mean = TRUE,
           zero = NULL)
{





  lepsil <- lepsilon
  iepsil <- iepsilon


  lshape <- as.list(substitute(lshape))
  e.shape <- link2list(lshape)
  l.shape <- attr(e.shape, "function.name")

  lscale <- as.list(substitute(lscale))
  e.scale <- link2list(lscale)
  l.scale <- attr(e.scale, "function.name")

  lepsil <- as.list(substitute(lepsil))
  e.epsil <- link2list(lepsil)
  l.epsil <- attr(e.epsil, "function.name")

  if (!is.Numeric(nsimEIM, allowable.length = 1,
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
            namesof("shape",   l.shape, e.shape), ", ",
            namesof("scale",   l.scale, e.scale), ", ",
            namesof("epsilon", l.epsil, e.epsil), "\n",
            "Median:   qmakeham(p = 0.5, shape, scale, epsilon)"),

  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 3
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 3,
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

    Musual <- 3
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    mynames1 <- paste("shape",   if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("scale",   if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames3 <- paste("epsilon", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .l.shape , .e.shape , tag = FALSE),
          namesof(mynames2, .l.scale , .e.scale , tag = FALSE),
          namesof(mynames3, .l.epsil , .e.epsil , tag = FALSE))[
          interleave.VGAM(M, M = Musual)]


    if (!length(etastart)) {

      matC <- matrix(if (length( .iscale )) .iscale else 0 + NA,
                     n, ncoly, byrow = TRUE)
      matH <- matrix(if (length( .ishape )) .ishape else 0 + NA,
                     n, ncoly, byrow = TRUE)

      matE <- matrix(if (length( .iepsil )) .iepsil else 0.3,
                     n, ncoly, byrow = TRUE)


      shape.grid <- c(exp(-seq(4, 0.1, len = 05)), 1,
                      exp( seq(0.1, 4, len = 05)))
      scale.grid <- c(exp(-seq(4, 0.1, len = 05)), 1,
                      exp( seq(0.1, 4, len = 05)))



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
            getMaxMin(scale.grid,
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
      } # spp.





      epsil.grid <- c(exp(-seq(4, 0.1, len = 05)), 1,
                      exp( seq(0.1, 1, len = 05)))
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
            getMaxMin(epsil.grid,
                      objfun = makeham.Loglikfun2,
                      y = yvec, x = x, w = wvec,
                      extraargs = list(Shape = matH[1, spp.],
                                       Scale = matC[1, spp.]))

        matE[, spp.] <- Init.epsil
      } # spp.


      etastart <- cbind(theta2eta(matH, .l.shape , .e.shape ),
                        theta2eta(matC, .l.scale , .e.scale ),
                        theta2eta(matE, .l.epsil , .e.epsil ))[,
                        interleave.VGAM(M, M = Musual)]
    } # End of !length(etastart)
  }), list(
            .l.shape = l.shape, .l.scale = l.scale, .l.epsil = l.epsil,
            .e.shape = e.shape, .e.scale = e.scale, .e.epsil = e.epsil,
            .ishape = ishape, .iscale = iscale, .iepsil = iepsil
          ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .l.shape , .e.shape )
    scale <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .l.scale , .e.scale )
    epsil <- eta2theta(eta[, c(FALSE, FALSE, TRUE)], .l.epsil , .e.epsil )
    qmakeham(p = 0.5, shape = shape, scale = scale, epsil = epsil)
  }, list(
            .l.shape = l.shape, .l.scale = l.scale, .l.epsil = l.epsil,
            .e.shape = e.shape, .e.scale = e.scale, .e.epsil = e.epsil
         ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <-
      c(rep( .l.shape , length = ncoly),
        rep( .l.scale , length = ncoly),
        rep( .l.epsil , length = ncoly))[interleave.VGAM(M, M = Musual)]
    temp.names <- c(mynames1, mynames2, mynames3)[
                    interleave.VGAM(M, M = Musual)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for(ii in 1:ncoly) {
      misc$earg[[Musual*ii-2]] <- .e.shape
      misc$earg[[Musual*ii-1]] <- .e.scale
      misc$earg[[Musual*ii  ]] <- .e.epsil
    }

    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list(
            .l.shape = l.shape, .l.scale = l.scale, .l.epsil = l.epsil,
            .e.shape = e.shape, .e.scale = e.scale, .e.epsil = e.epsil,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    shape <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .l.shape , .e.shape )
    scale <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .l.scale , .e.scale )
    epsil <- eta2theta(eta[, c(FALSE, FALSE, TRUE)], .l.epsil , .e.epsil )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dmakeham(x = y, shape = shape, scale = scale,
                          epsil = epsil, log = TRUE))
    }
  }, list(
            .l.shape = l.shape, .l.scale = l.scale, .l.epsil = l.epsil,
            .e.shape = e.shape, .e.scale = e.scale, .e.epsil = e.epsil
         ))),
  vfamily = c("makeham"),
 
  deriv = eval(substitute(expression({
    Musual <- 3
    shape <- eta2theta(eta[, c(TRUE, FALSE, FALSE), drop = FALSE],
                       .l.shape , .e.shape )
    scale <- eta2theta(eta[, c(FALSE, TRUE, FALSE), drop = FALSE],
                       .l.scale , .e.scale )
    epsil <- eta2theta(eta[, c(FALSE, FALSE, TRUE), drop = FALSE],
                       .l.epsil , .e.epsil )


    temp2 <- exp(y * scale)
    temp3 <- epsil + shape * temp2
    dl.dshape <- temp2 / temp3 - expm1(y * scale) / scale
    dl.dscale <- shape * y * temp2 / temp3 +
                 shape * expm1(y * scale) / scale^2 -
                 shape * y * temp2 / scale

    dl.depsil <- 1 / temp3 - y

    dshape.deta <- dtheta.deta(shape, .l.shape , .e.shape )
    dscale.deta <- dtheta.deta(scale, .l.scale , .e.scale )
    depsil.deta <- dtheta.deta(epsil, .l.epsil , .e.epsil )

    dthetas.detas <- cbind(dshape.deta, dscale.deta, depsil.deta)
    myderiv <- c(w) * cbind(dl.dshape,
                            dl.dscale,
                            dl.depsil) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list(
            .l.shape = l.shape, .l.scale = l.scale, .l.epsil = l.epsil,
            .e.shape = e.shape, .e.scale = e.scale, .e.epsil = e.epsil
          ))),


  weight = eval(substitute(expression({

    NOS <- M / Musual
    dThetas.detas <- dthetas.detas[, interleave.VGAM(M, M = Musual)]

    wz <- matrix(0.0, n, M + M - 1 + M - 2) # wz has half-bw 3

    ind1 <- iam(NA, NA, M = Musual, both = TRUE, diag = TRUE)


    for(spp. in 1:NOS) {
      run.varcov <- 0
      Shape <- shape[, spp.]
      Scale <- scale[, spp.]
      Epsil <- epsil[, spp.]




      if (FALSE && intercept.only && .oim.mean ) {

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

      for(ii in 1:( .nsimEIM )) {
        ysim <- rmakeham(n = n, shape = Shape, scale = Scale,
                         epsil = Epsil)
if (ii < 3) {
}

        temp2 <- exp(ysim * Scale)
        temp3 <- Epsil + Shape * temp2
 if (!is.Numeric(temp2))
  stop("temp2 is not Numeric")
 if (!is.Numeric(temp3))
  stop("temp3 is not Numeric")
        dl.dshape <- temp2 / temp3 - expm1(ysim * Scale) / Scale
        dl.dscale <- Shape * ysim * temp2 / temp3 +
                     Shape * expm1(ysim * Scale) / Scale^2 -
                     Shape * ysim * temp2 / Scale
        dl.depsil <- 1 / temp3 - ysim



        temp7 <- cbind(dl.dshape, dl.dscale, dl.depsil)
if (ii < 3) {
}
        run.varcov <- run.varcov +
                      temp7[, ind1$row.index] *
                      temp7[, ind1$col.index]
      }
      run.varcov <- cbind(run.varcov / .nsimEIM )

    }



      for (ilocal in 1:ncol(run.varcov)) {
        indexInf <- is.finite(run.varcov[, ilocal])
        run.varcov[!indexInf, ilocal] <-
          mean(run.varcov[indexInf, ilocal])
      }



      wz1 <- if (intercept.only)
          matrix(colMeans(run.varcov, na.rm = TRUE),
                 nrow = n, ncol = ncol(run.varcov), byrow = TRUE) else
          run.varcov


      wz1 <- wz1 * dThetas.detas[, Musual * (spp. - 1) + ind1$row] *
                   dThetas.detas[, Musual * (spp. - 1) + ind1$col]


      for(jay in 1:Musual)
        for(kay in jay:Musual) {
          cptr <- iam((spp. - 1) * Musual + jay,
                      (spp. - 1) * Musual + kay,
                      M = M)
          wz[, cptr] <- wz1[, iam(jay, kay, M = Musual)]
        }
    } # End of for(spp.) loop



    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / Musual)
  }), list(
            .l.shape = l.shape, .l.scale = l.scale, .l.epsil = l.epsil,
            .e.shape = e.shape, .e.scale = e.scale, .e.epsil = e.epsil,
            .nsimEIM = nsimEIM, .oim.mean = oim.mean ))))
} # makeham()








dgompertz <- function(x, shape, scale = 1, log = FALSE) {

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



pgompertz <- function(q, shape, scale = 1) {

  LLL <- max(length(q), length(shape), length(scale))
  if (length(q)       != LLL) q       <- rep(q,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)

  ans <- -expm1((-shape / scale) * expm1(scale * q))
  ans[(q <= 0)] <- 0
  ans[shape <= 0 | scale <= 0] <- NaN
  ans[q == Inf] <- 1
  ans
}


qgompertz <- function(p, shape, scale = 1) {

  LLL <- max(length(p), length(shape), length(scale))
  if (length(p)       != LLL) p       <- rep(p,       length.out = LLL)
  if (length(shape)   != LLL) shape   <- rep(shape,   length.out = LLL)
  if (length(scale)   != LLL) scale   <- rep(scale,   length.out = LLL)

  ans <- log1p((-scale / shape) * log1p(-p)) / scale
  ans[p < 0] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans[p > 1] <- NaN
  ans[shape <= 0 | scale <= 0] <- NaN
  ans
}


rgompertz <- function(n, shape, scale = 1) {
  qgompertz(runif(n), shape = shape, scale = scale)
}







gompertz.control <- function(save.weight = TRUE, ...)
{
  list(save.weight = save.weight)
}


 gompertz <-
  function(lshape = "loge", lscale = "loge",
           ishape = NULL,   iscale = NULL,
           nsimEIM = 500,
           zero = NULL)
{



  lshape <- as.list(substitute(lshape))
  e.shape <- link2list(lshape)
  l.shape <- attr(e.shape, "function.name")

  lscale <- as.list(substitute(lscale))
  e.scale <- link2list(lscale)
  l.scale <- attr(e.scale, "function.name")



  if (!is.Numeric(nsimEIM, allowable.length = 1,
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
            namesof("shape", l.shape, e.shape ), ", ",
            namesof("scale", l.scale, e.scale ), "\n",
            "Median:     scale * log(2 - 1 / shape)"),

  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
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
    Musual <- 2
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    mynames1 <- paste("shape", if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("scale", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .l.shape , .e.shape , tag = FALSE),
          namesof(mynames2, .l.scale , .e.scale , tag = FALSE))[
          interleave.VGAM(M, M = Musual)]



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
            getMaxMin(scale.grid,
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
      } # spp.

      etastart <- cbind(theta2eta(matH, .l.shape , .e.shape ),
                        theta2eta(matC, .l.scale , .e.scale ))[,
                        interleave.VGAM(M, M = Musual)]
    } # End of !length(etastart)
  }), list( .l.shape = l.shape, .l.scale = l.scale,
            .e.shape = e.shape, .e.scale = e.scale,
            .ishape = ishape, .iscale = iscale
          ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta[, c(TRUE, FALSE)], .l.shape , .e.shape )
    scale <- eta2theta(eta[, c(FALSE, TRUE)], .l.scale , .e.scale )
    log1p((scale / shape) * log(2)) / scale
  }, list( .l.shape = l.shape, .l.scale = l.scale,
           .e.shape = e.shape, .e.scale = e.scale ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <-
      c(rep( .l.shape , length = ncoly),
        rep( .l.scale , length = ncoly))[interleave.VGAM(M, M = Musual)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = Musual)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for(ii in 1:ncoly) {
      misc$earg[[Musual*ii-1]] <- .e.shape
      misc$earg[[Musual*ii  ]] <- .e.scale
    }

    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list( .l.shape = l.shape, .l.scale = l.scale,
            .e.shape = e.shape, .e.scale = e.scale,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    shape <- eta2theta(eta[, c(TRUE, FALSE)], .l.shape , .e.shape )
    scale <- eta2theta(eta[, c(FALSE, TRUE)], .l.scale , .e.scale )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dgompertz(x = y, shape = shape,
                           scale = scale, log = TRUE))
    }
    }, list( .l.shape = l.shape, .l.scale = l.scale,
             .e.shape = e.shape, .e.scale = e.scale ))),
  vfamily = c("gompertz"),
 
  deriv = eval(substitute(expression({
    Musual <- 2
    shape <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .l.shape ,
                       .e.shape )
    scale <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .l.scale ,
                       .e.scale )


    temp2 <- exp(y * scale)
    temp4 <- -expm1(y * scale)
    dl.dshape <- 1 / shape + temp4 / scale
    dl.dscale <- y * (1 - shape * temp2 / scale) -
                 shape * temp4 / scale^2

    dshape.deta <- dtheta.deta(shape, .l.shape , .e.shape )
    dscale.deta <- dtheta.deta(scale, .l.scale , .e.scale )

    dthetas.detas <- cbind(dshape.deta, dscale.deta)
    myderiv <- c(w) * cbind(dl.dshape, dl.dscale) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list( .l.shape = l.shape, .l.scale = l.scale,
            .e.shape = e.shape, .e.scale = e.scale ))),


  weight = eval(substitute(expression({

    NOS <- M / Musual
    dThetas.detas <- dthetas.detas[, interleave.VGAM(M, M = Musual)]

    wz <- matrix(0.0, n, M + M - 1) # wz is 'tridiagonal' 

    ind1 <- iam(NA, NA, M = Musual, both = TRUE, diag = TRUE)


    for(spp. in 1:NOS) {
      run.varcov <- 0
      Shape <- shape[, spp.]
      Scale <- scale[, spp.]

      for(ii in 1:( .nsimEIM )) {
        ysim <- rgompertz(n = n, shape = Shape, scale = Scale)
if (ii < 3) {
}

        temp2 <- exp(ysim * scale)
        temp4 <- -expm1(ysim * scale)
        dl.dshape <- 1 / shape + temp4 / scale
        dl.dscale <- ysim * (1 - shape * temp2 / scale) -
                     shape * temp4 / scale^2


        temp7 <- cbind(dl.dshape, dl.dscale)
        run.varcov <- run.varcov +
                      temp7[, ind1$row.index] *
                      temp7[, ind1$col.index]
      }
      run.varcov <- cbind(run.varcov / .nsimEIM )

      wz1 <- if (intercept.only)
          matrix(colMeans(run.varcov),
                 nrow = n, ncol = ncol(run.varcov), byrow = TRUE) else
          run.varcov

      wz1 <- wz1 * dThetas.detas[, Musual * (spp. - 1) + ind1$row] *
                   dThetas.detas[, Musual * (spp. - 1) + ind1$col]


      for(jay in 1:Musual)
        for(kay in jay:Musual) {
          cptr <- iam((spp. - 1) * Musual + jay,
                      (spp. - 1) * Musual + kay,
                      M = M)
          wz[, cptr] <- wz1[, iam(jay, kay, M = Musual)]
        }
    } # End of for(spp.) loop



    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / Musual)
  }), list( .l.scale = l.scale,
            .e.scale = e.scale,
            .nsimEIM = nsimEIM ))))
} # gompertz()






dmoe <- function (x, alpha = 1, lambda = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(alpha), length(lambda))
  if (length(x)      != LLL) x      <- rep(x,      length.out = LLL)
  if (length(alpha)  != LLL) alpha  <- rep(alpha,  length.out = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, length.out = LLL)

  index0 = (x < 0)
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



rmoe <- function (n, alpha = 1, lambda = 1)
{

  qmoe(p = runif(n), alpha = alpha, lambda = lambda)
}




exponential.mo.control <- function(save.weight = TRUE, ...)
{
    list(save.weight = save.weight)
}




 exponential.mo <-
  function(lalpha = "loge", llambda = "loge",
           ealpha = list(), elambda = list(),
           ialpha = 1,      ilambda = NULL,
           imethod = 1,
           nsimEIM = 200,
           zero = NULL)
{

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



  if (!is.Numeric(nsimEIM, allowable.length = 1,
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


  if (!is.Numeric(imethod, allowable.length = 1,
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
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
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

    Musual <- 2
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    mynames1 <- paste("alpha",   if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("lambda",  if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lalpha0 , .ealpha0 , tag = FALSE),
          namesof(mynames2, .llambda , .elambda , tag = FALSE))[
          interleave.VGAM(M, M = Musual)]



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
        Lambda.init <- getMaxMin(lambda.grid,
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
      } # spp.

      etastart <- cbind(theta2eta(matA, .lalpha0, .ealpha0 ),
                        theta2eta(matL, .llambda, .elambda ))[,
                        interleave.VGAM(M, M = Musual)]
      mustart <- NULL # Since etastart has been computed.
    } # End of !length(etastart)
  }), list( .lalpha0 = lalpha0, .llambda = llambda,
            .ealpha0 = ealpha0, .elambda = elambda,
            .ialpha0 = ialpha0, .ilambda = ilambda,
            .imethod = imethod
          ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha0 = eta2theta(eta[, c(TRUE, FALSE)], .lalpha0 , .ealpha0 )
    lambda = eta2theta(eta[, c(FALSE, TRUE)], .llambda , .elambda )
    log(3 - alpha0) / lambda
  }, list( .lalpha0 = lalpha0, .llambda = llambda,
           .ealpha0 = ealpha0, .elambda = elambda ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <-
      c(rep( .lalpha0 , length = ncoly),
        rep( .llambda , length = ncoly))[interleave.VGAM(M, M = Musual)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = Musual)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for(ii in 1:ncoly) {
      misc$earg[[Musual*ii-1]] <- .ealpha0
      misc$earg[[Musual*ii  ]] <- .elambda
    }

    misc$Musual <- Musual
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$nsimEIM = .nsimEIM
  }), list( .lalpha0 = lalpha0, .llambda = llambda,
            .ealpha0 = ealpha0, .elambda = elambda,
            .nsimEIM = nsimEIM,
            .imethod = imethod ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    alpha0 = eta2theta(eta[, c(TRUE, FALSE)], .lalpha0 , .ealpha0 )
    lambda = eta2theta(eta[, c(FALSE, TRUE)], .llambda , .elambda )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * log(dmoe(x = y, alpha = alpha0,
                          lambda = lambda)))
    }
    }, list( .lalpha0 = lalpha0, .llambda = llambda,
             .ealpha0 = ealpha0, .elambda = elambda ))),
  vfamily = c("exponential.mo"),
 
  deriv = eval(substitute(expression({
    Musual <- 2
    alpha0 = eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lalpha0 ,
                       .ealpha0 )
    lambda = eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .llambda ,
                       .elambda )

    temp2 = (expm1(lambda * y) + alpha0)
    dl.dalpha0 = -2 / temp2
    dl.dlambda = 1 / lambda + y - 2 * y * exp(lambda * y) / temp2

    dalpha0.deta = dtheta.deta(alpha0, .lalpha0 , .ealpha0 )
    dlambda.deta = dtheta.deta(lambda, .llambda , .elambda )

    dthetas.detas = cbind(dalpha0.deta,
                          dlambda.deta)
    myderiv = c(w) * cbind(dl.dalpha0, dl.dlambda) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list( .lalpha0 = lalpha0, .llambda = llambda,
            .ealpha0 = ealpha0, .elambda = elambda ))),


  weight = eval(substitute(expression({

    NOS = M / Musual
    dThetas.detas = dthetas.detas[, interleave.VGAM(M, M = Musual)]

    wz = matrix(0.0, n, M + M - 1) # wz is 'tridiagonal' 

    ind1 = iam(NA, NA, M = Musual, both = TRUE, diag = TRUE)


    for(spp. in 1:NOS) {
      run.varcov = 0
      Alph = alpha0[, spp.]
      Lamb = lambda[, spp.]

      for(ii in 1:( .nsimEIM )) {
        ysim = rmoe(n = n, alpha = Alph, lambda = Lamb)
if (ii < 3) {
}

        temp2 = (expm1(lambda * ysim) + alpha0)
        dl.dalpha0 = -2 / temp2
        dl.dlambda = 1 / lambda + ysim -
                     2 * ysim * exp(lambda * ysim) / temp2


        temp3 = cbind(dl.dalpha0, dl.dlambda)
        run.varcov = run.varcov +
                     temp3[, ind1$row.index] *
                     temp3[, ind1$col.index]
      }
      run.varcov = cbind(run.varcov / .nsimEIM)

      wz1 = if (intercept.only)
          matrix(colMeans(run.varcov),
                 nrow = n, ncol = ncol(run.varcov), byrow = TRUE) else
          run.varcov

      wz1 = wz1 * dThetas.detas[, Musual * (spp. - 1) + ind1$row] *
                  dThetas.detas[, Musual * (spp. - 1) + ind1$col]


      for(jay in 1:Musual)
        for(kay in jay:Musual) {
          cptr = iam((spp. - 1) * Musual + jay,
                     (spp. - 1) * Musual + kay,
                     M = M)
          wz[, cptr] = wz1[, iam(jay, kay, M = Musual)]
        }
    } # End of for(spp.) loop




    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / Musual)
  }), list( .llambda = llambda,
            .elambda = elambda,
            .nsimEIM = nsimEIM ))))
} # exponential.mo()






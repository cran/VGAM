# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.



























 getMaxMin <- function(vov, objfun, y, x, w, extraargs = NULL,
                       maximize = TRUE, abs.arg = FALSE,
                       ret.objfun = FALSE) {
  if (!is.vector(vov))
    stop("'vov' must be a vector")
  objvals <- vov
  for(ii in 1:length(vov))
    objvals[ii] <- objfun(vov[ii], y = y, x = x, w = w,
                          extraargs = extraargs)
  try.this <- if (abs.arg) {
               if (maximize) vov[abs(objvals) == max(abs(objvals))] else
               vov[abs(objvals) == min(abs(objvals))]
             } else {
               if (maximize) vov[objvals == max(objvals)] else
               vov[objvals == min(objvals)]
             }
  if (!length(try.this))
    stop("something has gone wrong!")
  ans <- if (length(try.this) == 1)
    try.this else sample(try.this, size = 1)
  if (ret.objfun) c(ans, objvals[ans == vov]) else ans
}




 mccullagh89 <- function(ltheta = "rhobit",
                         lnu = logoff(offset = 0.5),
                         itheta = NULL, inu = NULL,
                         zero = NULL) {



  ltheta <- as.list(substitute(ltheta))
  etheta <- link2list(ltheta)
  ltheta <- attr(etheta, "function.name")

  lnu <- as.list(substitute(lnu))
  enu <- link2list(lnu)
  lnu <- attr(enu, "function.name")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  new("vglmff",
  blurb = c("McCullagh (1989)'s distribution \n",
            "f(y) = (1-2*theta*y+theta^2)^(-nu) * [1 - y^2]^(nu-1/2) /\n",
            "       Beta[nu+1/2, 1/2], ",
            "  -1 < y < 1, -1 < theta < 1, nu > -1/2\n",
            "Links:     ",
            namesof("theta", ltheta, earg = etheta), ", ",
            namesof("nu",    lnu,    earg = enu),
            "\n", "\n",
            "Mean:     nu*theta/(1+nu)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    w.y.check(w, y)

    y <- as.numeric(y)
    if (any(y <= -1 | y >= 1))
      stop("all y values must be in (-1, 1)")

    predictors.names <-
      c(namesof("theta", .ltheta , earg = .etheta, tag = FALSE),
        namesof("nu",    .lnu ,    earg = .enu,    tag = FALSE))

    if (!length(etastart)) {
      theta.init <- if (length( .itheta ))
          rep( .itheta, length = n) else {
          mccullagh89.aux <- function(thetaval, y, x, w, extraargs)
          mean((y-thetaval)*(thetaval^2-1)/(1-2*thetaval*y+thetaval^2))
          theta.grid <- seq(-0.9, 0.9, by=0.05)
          try.this <- getMaxMin(theta.grid, objfun = mccullagh89.aux,
                               y = y,  x = x, w = w, maximize = FALSE,
                               abs.arg = TRUE)
          try.this <- rep(try.this, length.out = n)
          try.this
      }
      tmp <- y / (theta.init - y)
      tmp[tmp < -0.4] <- -0.4
      tmp[tmp > 10.0] <- 10.0
      nu.init <- rep(if (length( .inu)) .inu else tmp, length = n)
      nu.init[!is.finite(nu.init)] <- 0.4
      etastart <-
        cbind(theta2eta(theta.init, .ltheta , earg = .etheta ),
              theta2eta(nu.init,    .lnu,     earg = .enu ))
    }
  }), list( .ltheta = ltheta, .lnu = lnu, .inu = inu, .itheta = itheta,
            .etheta = etheta, .enu = enu ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Theta <- eta2theta(eta[, 1], .ltheta , earg = .etheta )
    nu <- eta2theta(eta[, 2], .lnu, earg = .enu )
    nu * Theta / (1 + nu)
  }, list( .ltheta = ltheta, .lnu = lnu,
           .etheta = etheta, .enu = enu ))),
  last = eval(substitute(expression({
    misc$link <-    c("theta" = .ltheta , "nu" = .lnu )

    misc$earg <- list("theta" = .etheta , "nu" = .enu )

  }), list( .ltheta = ltheta, .lnu = lnu, .etheta = etheta, .enu = enu ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Theta <- eta2theta(eta[, 1], .ltheta , earg = .etheta )
    nu <- eta2theta(eta[, 2], .lnu, earg = .enu )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
      sum(c(w) * ((nu-0.5)*log1p(-y^2) -
                nu * log1p(-2*Theta*y + Theta^2) -
              lbeta(nu + 0.5, 0.5)))
  }, list( .ltheta = ltheta, .lnu = lnu,
           .etheta = etheta, .enu = enu ))),
  vfamily = c("mccullagh89"),
  deriv = eval(substitute(expression({
    Theta <- eta2theta(eta[, 1], .ltheta , earg = .etheta )
    nu    <- eta2theta(eta[, 2], .lnu, earg = .enu )

    dTheta.deta <- dtheta.deta(Theta, .ltheta , earg = .etheta )
    dnu.deta <- dtheta.deta(nu, .lnu, earg = .enu )

    dl.dTheta <- 2 * nu * (y-Theta) / (1 -2*Theta*y + Theta^2)
    dl.dnu <- log1p(-y^2) - log1p(-2*Theta*y + Theta^2) -
             digamma(nu + 0.5) + digamma(nu + 1)

    c(w) * cbind(dl.dTheta * dTheta.deta,
                 dl.dnu * dnu.deta)
  }), list( .ltheta = ltheta, .lnu = lnu,
            .etheta = etheta, .enu = enu ))),
  weight = eval(substitute(expression({
    d2l.dTheta2 <- (2 * nu^2 / (1+nu)) / (1-Theta^2)
    d2l.dnu2 <- trigamma(nu+0.5) - trigamma(nu+1)

    wz <- matrix(as.numeric(NA), n, M) # diagonal matrix
    wz[, iam(1, 1, M)] <- d2l.dTheta2 * dTheta.deta^2
    wz[, iam(2, 2, M)] <- d2l.dnu2 * dnu.deta^2

    c(w) * wz
  }), list( .ltheta = ltheta, .lnu = lnu ))))
}




hzeta.control <- function(save.weight = TRUE, ...) {
  list(save.weight = save.weight)
}



 hzeta <- function(link = "loglog", ialpha = NULL, nsimEIM = 100) {

  stopifnot(ialpha > 0)
  stopifnot(nsimEIM > 10,
            length(nsimEIM) == 1,
            nsimEIM == round(nsimEIM))



  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Haight's Zeta distribution f(y) = (2y-1)^(-alpha) - ",
            "(2y+1)^(-alpha),\n",
            "    alpha>0, y = 1, 2,....\n\n",
            "Link:    ",
            namesof("alpha", link, earg = earg), "\n\n",
            "Mean:     (1-2^(-alpha)) * zeta(alpha) if alpha>1",
            "\n",
            "Variance: (1-2^(1-alpha)) * zeta(alpha-1) - mean^2 if alpha>2"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y)

    if (any(y < 1))
      stop("all y values must be in 1, 2, 3,....")

    predictors.names <-
      namesof("alpha", .link , earg = .earg , tag = FALSE)

    if (!length(etastart)) {
      a.init <- if (length( .ialpha)) .ialpha else {
        if ((meany <- weighted.mean(y, w)) < 1.5) 3.0 else
        if (meany < 2.5) 1.4 else 1.1 
      }
      a.init <- rep(a.init, length = n) 
      etastart <- theta2eta(a.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .ialpha = ialpha ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha <- eta2theta(eta, .link , earg = .earg )
    mu <- (1-2^(-alpha)) * zeta(alpha)
    mu[alpha <= 1] <- Inf
    mu
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c(alpha = .link)

    misc$earg <- list(alpha = .earg )

    misc$nsimEIM <- .nsimEIM

  }), list( .link = link, .earg = earg, .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    alpha <- eta2theta(eta, .link , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dhzeta(x = y, alpha = alpha, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("hzeta"),
  deriv = eval(substitute(expression({
    alpha <- eta2theta(eta, .link , earg = .earg ) 

    dalpha.deta <- dtheta.deta(alpha, .link , earg = .earg )

    d3 <- deriv3(~ log((2*y-1)^(-alpha) - (2*y+1)^(-alpha)),
                "alpha", hessian = FALSE)
    eval.d3 <- eval(d3)

    dl.dalpha <-  attr(eval.d3, "gradient")

    c(w) * dl.dalpha * dalpha.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    sd3 <- deriv3(~ log((2*ysim-1)^(-alpha) - (2*ysim+1)^(-alpha)),
                 "alpha", hessian = FALSE)
    run.var <- 0
    for(ii in 1:( .nsimEIM )) {
      ysim <- rhzeta(n, alpha=alpha)
      eval.sd3 <- eval(sd3)
      dl.dalpha <-  attr(eval.d3, "gradient")
      rm(ysim)
      temp3 <- dl.dalpha
      run.var <- ((ii-1) * run.var + temp3^2) / ii
    }
    wz <- if (intercept.only)
        matrix(colMeans(cbind(run.var)),
               n, dimm(M), byrow = TRUE) else cbind(run.var)

    wz <- wz * dalpha.deta^2
    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))
}




dhzeta <- function(x, alpha, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(alpha, positive = TRUE))
    stop("'alpha' must be numeric and have positive values")

  nn <- max(length(x), length(alpha))
  x <- rep(x, length.out = nn);
  alpha <- rep(alpha, length.out = nn)

  ox <- !is.finite(x)
  zero <- ox | round(x) != x | x < 1
  ans <- rep(0, length.out = nn)
  ans[!zero] <- (2*x[!zero]-1)^(-alpha[!zero]) -
                (2*x[!zero]+1)^(-alpha[!zero])
  if (log.arg) log(ans) else ans
}


phzeta <- function(q, alpha) {


  nn <- max(length(q), length(alpha))
  q <- rep(q, length.out = nn)
  alpha <- rep(alpha, length.out = nn)
  oq <- !is.finite(q)
  zero <- oq | q < 1
  q <- floor(q)
  ans <- 0 * q
  ans[!zero] <- 1 - (2*q[!zero]+1)^(-alpha[!zero])

  ans[alpha <= 0] <- NaN

  ans
}


qhzeta <- function(p, alpha) {

  if (!is.Numeric(p, positive = TRUE) ||
      any(p >= 1))
    stop("argument 'p' must have values inside the interval (0,1)")

  nn <- max(length(p), length(alpha))
  p <- rep(p, length.out = nn)
  alpha <- rep(alpha, length.out = nn)
  ans <- (((1 - p)^(-1/alpha) - 1) / 2) # p is in (0,1)
  ans[alpha <= 0] <- NaN
  floor(ans + 1)
}


rhzeta <- function(n, alpha) {


  ans <- (runif(n)^(-1/alpha) - 1) / 2
  ans[alpha <= 0] <- NaN
  floor(ans + 1)
}






 dirmultinomial <- function(lphi = "logit",
                            iphi = 0.10, parallel = FALSE, zero = "M") {




  lphi <- as.list(substitute(lphi))
  ephi <- link2list(lphi)
  lphi <- attr(ephi, "function.name")


  if (length(zero) && 
      !(is.Numeric(zero, integer.valued = TRUE, positive = TRUE) ||
      is.character(zero )))
    stop("bad input for argument 'zero'")

  if (!is.Numeric(iphi, positive = TRUE) ||
      max(iphi) >= 1.0)
    stop("bad input for argument 'iphi'")




  new("vglmff",
  blurb = c("Dirichlet-multinomial distribution\n\n",
            "Links:    ",
            "log(prob[1]/prob[M]), ..., log(prob[M-1]/prob[M]), ",
            namesof("phi", lphi, earg = ephi), "\n", "\n",
            "Mean:     shape_j / sum_j(shape_j)"),
  constraints = eval(substitute(expression({
    .ZERO <- .zero
    if (is.character( .ZERO)) .ZERO <- eval(parse(text = .ZERO))
    .PARALLEL <- .parallel
    if (is.logical( .PARALLEL) && .PARALLEL) {
      mycmatrix <- if (length( .ZERO ))
          stop("can only handle parallel = TRUE when zero = NULL") else
          cbind(rbind(matrix(1, M - 1, 1), 0),
                rbind(matrix(0, M - 1, 1), 1))
    } else {
      mycmatrix <- if (M == 1) diag(1) else diag(M)
    }
    constraints <- cm.vgam(mycmatrix, x, .PARALLEL ,
                           constraints, apply.int = TRUE)
    constraints <- cm.zero.vgam(constraints, x, .ZERO , M)
  }), list( .parallel = parallel, .zero = zero ))),
  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- TRUE
    eval(process.categorical.data.vgam)

    if (length(mustart.orig))
      mustart <- mustart.orig

    y <- as.matrix(y)
    ycount <- as.matrix(y * c(w))
    M <- ncol(y)

    if (max(abs(ycount - round(ycount))) > 1.0e-6)
      warning("there appears to be non-integer responses")

    if (min(ycount) < 0)
      stop("all values of the response (matrix) must be non-negative")

    predictors.names <-
      c(paste("log(prob[,", 1:(M-1), "]/prob[,", M, "])", sep = ""),
        namesof("phi", .lphi , short = TRUE))

    extra$n2 <- w # aka omega, must be integer # as.vector(apply(y, 1, sum))

    if (!length(etastart)) {
      if (length(mustart.orig)) {
        prob.init <- mustart
      } else {
        prob.init <- colSums(ycount)
        prob.init <- prob.init / sum(prob.init)
        prob.init <- matrix(prob.init, n, M, byrow = TRUE)
      }

      phi.init <- rep( .iphi , length.out = n)
      etastart <-
        cbind(log(prob.init[, -M] / prob.init[, M]),
              theta2eta(phi.init, .lphi , earg = .ephi ))
    }

    mustart <- NULL # Since etastart has been computed.

  }), list( .lphi = lphi, .ephi = ephi, .iphi = iphi ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    M <- if (is.matrix(eta)) ncol(eta) else 1
    temp <- cbind(exp(eta[, -M, drop = FALSE]), 1)
    prop.table(temp, 1)
  }, list( .ephi = ephi, .lphi = lphi ))),
  last = eval(substitute(expression({

    misc$link <- c(rep("loge", length = M-1), .lphi )
    names(misc$link) <- c(
      paste("prob[,", 1:(M-1), "]/prob[,", M, "])", sep = ""),
      "phi")

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:(M-1))
      misc$earg[[ii]] <- list()
    misc$earg[[M]] <- .ephi

    misc$expected <- TRUE

    if (intercept.only) {
      misc$shape <- probs[1,] * (1/phi[1]-1) # phi & probs computed in @deriv
    }
  }), list( .ephi = ephi, .lphi = lphi ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    M <- if (is.matrix(eta)) ncol(eta) else 1
    probs <- cbind(exp(eta[, -M]), 1)
    probs <- prop.table(probs, 1)
    phi <- eta2theta(eta[, M], .lphi , earg = .ephi )
    n <- length(phi)
    ycount <- as.matrix(y * c(w))

    ycount <- round(ycount)

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      ans <- rep(0.0, length.out = n)
      omega <- extra$n2
      for(jay in 1:M) {
        maxyj <- max(ycount[, jay])
        loopOveri <- (n < maxyj)
        if (loopOveri) {
          for(iii in 1:n) {
              rrr <- 1:ycount[iii, jay] # a vector
              if (ycount[iii, jay] > 0)
                ans[iii] <- ans[iii] + sum(log((1-phi[iii]) *
                            probs[iii, jay] + (rrr-1)*phi[iii]))
          }
        } else {
          for(rrr in 1:maxyj) {
              index <- (rrr <= ycount[, jay]) & (ycount[, jay] > 0)
              if (any(index))
                  ans[index] <- ans[index] + log((1-phi[index]) *
                                probs[index, jay] + (rrr-1) * phi[index])
          }
        }
      } # end of jay loop

      maxomega <- max(omega)
      loopOveri <- n < maxomega
      if (loopOveri) {
        for(iii in 1:n) {
          rrr <- 1:omega[iii]
          ans[iii]<- ans[iii] - sum(log1p(-phi[iii] + (rrr-1) * phi[iii]))
        }
      } else {
        for(rrr in 1:maxomega) {
          ind8 <- rrr <= omega
          ans[ind8] <- ans[ind8] - log1p(-phi[ind8] + (rrr-1) * phi[ind8])
        }
      }
      sum(ans)
    }
  }, list( .ephi = ephi, .lphi = lphi ))),
  vfamily = c("dirmultinomial"),
  deriv = eval(substitute(expression({
    probs <- cbind(exp(eta[, -M]), 1)
    probs <- prop.table(probs, 1)

    phi <- eta2theta(eta[, M], .lphi , earg = .ephi )

    dl.dprobs <- matrix(0.0, n, M-1)
    dl.dphi <- rep(0.0, length.out = n)

    omega <- extra$n2
    ycount <- as.matrix(y * c(w))

    ycount <- round(ycount)

    for(jay in 1:M) {
      maxyj <- max(ycount[, jay])
      loopOveri <- n < maxyj
      if (loopOveri) {
        for(iii in 1:n) {
          rrr <- 1:ycount[iii, jay]
          if (ycount[iii, jay] > 0) {
            PHI <- phi[iii]
            dl.dphi[iii] <- dl.dphi[iii] +
 sum((rrr-1-probs[iii, jay]) / ((1-PHI)*probs[iii, jay] + (rrr-1)*PHI))

            tmp9 <- (1-PHI) / ((1-PHI)*probs[iii, jay] + (rrr-1)*PHI)
            if (jay < M) {
              dl.dprobs[iii, jay] <- dl.dprobs[iii, jay] + sum(tmp9)
            } else {
              for(jay2 in 1:(M-1))
                dl.dprobs[iii, jay2]<-dl.dprobs[iii, jay2]-sum(tmp9)
            }
          }
        }
      } else {
        for(rrr in 1:maxyj) {
          index <- (rrr <= ycount[, jay]) & (ycount[, jay] > 0)
          PHI <- phi[index]
          dl.dphi[index] <- dl.dphi[index] +
            (rrr-1-probs[index, jay]) / ((1-PHI)*probs[index, jay] +
            (rrr-1)*PHI)
          tmp9 <- (1-PHI) / ((1-PHI)*probs[index, jay] + (rrr-1)*PHI)
          if (jay < M) {
              dl.dprobs[index, jay] <- dl.dprobs[index, jay] + tmp9
          } else {
              for(jay2 in 1:(M-1))
                  dl.dprobs[index, jay2] <- dl.dprobs[index, jay2] - tmp9
          }
        }
      }
    } # end of jay loop
    maxomega <- max(omega)
    loopOveri <- n < maxomega
    if (loopOveri) {
      for(iii in 1:n) {
        rrr <- 1:omega[iii]
        dl.dphi[iii]<-dl.dphi[iii] - sum((rrr-2)/(1 + (rrr-2)*phi[iii]))
      }
    } else {
      for(rrr in 1:maxomega) {
        index <- rrr <= omega
        dl.dphi[index] <-
        dl.dphi[index] - (rrr-2)/(1 + (rrr-2)*phi[index])
      }
    }

    dprobs.deta <- probs[, -M] * (1 - probs[, -M]) # n x (M-1)
    dphi.deta <- dtheta.deta(phi, .lphi , earg = .ephi )

    ans <- cbind(dl.dprobs * dprobs.deta,
                 dl.dphi   * dphi.deta)
    ans
  }), list( .ephi = ephi, .lphi = lphi ))),
    weight = eval(substitute(expression({
      wz <- matrix(0, n, dimm(M))
      loopOveri <- (n < maxomega)
      if (loopOveri) {
          for(iii in 1:n) {
              rrr <- 1:omega[iii]  # A vector
              PHI <- phi[iii]
              pYiM.ge.rrr <- 1 - pbetabinom.ab(q = rrr-1,
                                               size = omega[iii],
                  shape1<-probs[iii, M]*(1/PHI-1),
                  shape2<-(1-probs[iii, M])*(1/PHI-1))  # A vector
              denomM <- ((1-PHI)*probs[iii, M] + (rrr-1)*PHI)^2  # A vector
              wz[iii, iam(M, M, M)] <- wz[iii, iam(M, M, M)] +
                      sum(probs[iii, M]^2 * pYiM.ge.rrr / denomM) -
                      sum(1 / (1 + (rrr-2)*PHI)^2)
              for(jay in 1:(M-1)) {
                  denomj <- ((1-PHI)*probs[iii, jay] + (rrr-1)*PHI)^2
                  pYij.ge.rrr <- 1 - pbetabinom.ab(q = rrr-1,
                                                   size = omega[iii],
                      shape1<-probs[iii, jay]*(1/PHI-1),
                      shape2<-(1-probs[iii, jay])*(1/PHI-1))
                  wz[iii, iam(jay, jay, M)] <- wz[iii, iam(jay, jay, M)] + 
                      sum(pYij.ge.rrr / denomj) + 
                      sum(pYiM.ge.rrr / denomM)
                  for(kay in jay:(M-1)) if (kay > jay) {
                    wz[iii, iam(jay, kay, M)] <- wz[iii, iam(jay, kay, M)] +
                        sum(pYiM.ge.rrr / denomM)
                  }
                  wz[iii, iam(jay, M, M)] <- wz[iii, iam(jay, M, M)] +
                          sum(probs[iii, jay] * pYij.ge.rrr / denomj) -
                          sum(probs[iii, M]   * pYiM.ge.rrr / denomM)
                  wz[iii, iam(M, M, M)] <- wz[iii, iam(M, M, M)] +
                          sum(probs[iii, jay]^2 * pYij.ge.rrr / denomj)
              } # end of jay loop
          } # end of iii loop
      } else {
          for(rrr in 1:maxomega) {
              ind5 <- rrr <= omega
              PHI <- phi[ind5]
              pYiM.ge.rrr <- 1 - pbetabinom.ab(q = rrr-1,
                                               size = omega[ind5],
                  shape1<-probs[ind5, M]*(1/PHI-1),
                  shape2<-(1-probs[ind5, M])*(1/PHI-1))
              denomM <- ((1-PHI)*probs[ind5, M] + (rrr-1)*PHI)^2
              wz[ind5, iam(M, M, M)] <- wz[ind5, iam(M, M, M)] +
                      probs[ind5, M]^2 * pYiM.ge.rrr / denomM -
                      1 / (1 + (rrr-2)*PHI)^2
              for(jay in 1:(M-1)) {
                  denomj <- ((1-PHI)*probs[ind5, jay] + (rrr-1)*PHI)^2
                  pYij.ge.rrr <- 1 - pbetabinom.ab(q = rrr-1,
                                                   size = omega[ind5],
                      shape1<-probs[ind5, jay]*(1/PHI-1),
                      shape2<-(1-probs[ind5, jay])*(1/PHI-1))
                  wz[ind5, iam(jay, jay, M)] <- wz[ind5, iam(jay, jay, M)] + 
                      pYij.ge.rrr / denomj + pYiM.ge.rrr / denomM 
                  for(kay in jay:(M-1)) if (kay > jay) {
                    wz[ind5, iam(jay, kay, M)] <- wz[ind5, iam(jay, kay, M)] +
                        pYiM.ge.rrr / denomM 
                  }
                  wz[ind5, iam(jay, M, M)] <- wz[ind5, iam(jay, M, M)] +
                      probs[ind5, jay] * pYij.ge.rrr / denomj -
                      probs[ind5, M]   * pYiM.ge.rrr / denomM
                  wz[ind5, iam(M, M, M)] <- wz[ind5, iam(M, M, M)] +
                      probs[ind5, jay]^2 * pYij.ge.rrr / denomj
              } # end of jay loop
          } # end of rrr loop
      }

      for(jay in 1:(M-1))
        for(kay in jay:(M-1))
          wz[, iam(jay, kay, M)] <- wz[, iam(jay, kay, M)] * (1-phi)^2
      for(jay in 1:(M-1))
        wz[, iam(jay, M, M)] <- wz[, iam(jay, M, M)] * (phi-1) / phi
      wz[, iam(M, M, M)] <- wz[, iam(M, M, M)] / phi^2

      d1Thetas.deta <- cbind(dprobs.deta,
                             dphi.deta)
      index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
      wz <- wz * d1Thetas.deta[, index$row] * d1Thetas.deta[, index$col]
      wz
  }), list( .ephi = ephi, .lphi = lphi ))))
}





dirmul.old <- function(link = "loge", init.alpha = 0.01,
                       parallel = FALSE, zero = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(zero) &&
    !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")

  if (!is.Numeric(init.alpha, positive = TRUE))
    stop("'init.alpha' must contain positive values only")


  new("vglmff",
  blurb = c("Dirichlet-Multinomial distribution\n\n",
            "Links:     ",
            namesof("shape1", link, earg = earg), ", ..., ",
            namesof("shapeM", link, earg = earg), "\n\n",
            "Posterior mean:    (n_j + shape_j)/(2*sum(n_j) + ",
                                "sum(shape_j))\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(1, M, 1), x, .parallel ,
                           constraints, apply.int = TRUE)
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .parallel = parallel, .zero = zero ))),
  initialize = eval(substitute(expression({
    y <- as.matrix(y)
    M <- ncol(y)
      if (any(y != round(y )))
        stop("all y values must be integer-valued")

      predictors.names <- namesof(paste("shape", 1:M, sep = ""),
                                  .link , earg = .earg , short = TRUE)

      extra$n2 <- rowSums(y) # Nb. don't multiply by 2
      extra$y  <- y

      if (!length(etastart)) {
        yy <- if (is.numeric( .init.alpha))
            matrix( .init.alpha, n, M, byrow = TRUE) else
            matrix(runif(n*M), n, M)
        etastart <- theta2eta(yy, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .init.alpha = init.alpha ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .link , earg = .earg )
    M <- if (is.matrix(eta)) ncol(eta) else 1
    sumshape <- as.vector(shape %*% rep(1, length.out = M))
    (extra$y + shape) / (extra$n2 + sumshape)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <- rep( .link , length = M)
    names(misc$link) <- paste("shape", 1:M, sep = "")

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$pooled.weight <- pooled.weight
  }), list( .link = link, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    shape <- eta2theta(eta, .link , earg = .earg )
    M <- if (is.matrix(eta)) ncol(eta) else 1
    sumshape <- as.vector(shape %*% rep(1, length.out = M))
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
    sum(c(w) * (lgamma(sumshape) - lgamma(extra$n2 + sumshape ))) +
    sum(c(w) * (lgamma(y + shape) - lgamma(shape )))
  }, list( .link = link, .earg = earg ))),
  vfamily = c("dirmul.old"),
  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .link , earg = .earg )

    sumshape <- as.vector(shape %*% rep(1, length.out = M))
    dl.dsh <- digamma(sumshape) - digamma(extra$n2 + sumshape) +
             digamma(y + shape) - digamma(shape)

    dsh.deta <- dtheta.deta(shape, .link , earg = .earg )

    c(w) * dl.dsh * dsh.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
    wz <- matrix(trigamma(sumshape) - trigamma(extra$n2 + sumshape),
                nrow = n, ncol = dimm(M))
    wz[, 1:M] <- wz[, 1:M] + trigamma(y + shape) - trigamma(shape)
    wz <- -wz * dsh.deta[, index$row] * dsh.deta[, index$col]


    if (TRUE && intercept.only) {
      sumw <- sum(w)
      for(ii in 1:ncol(wz))
        wz[, ii] <- sum(wz[, ii]) / sumw
      pooled.weight <- TRUE
      wz <- c(w) * wz # Put back the weights
    } else
        pooled.weight <- FALSE

    wz
  }), list( .link = link, .earg = earg ))))
}






rdiric <- function(n, shape, dimension = NULL) {

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (!is.numeric(dimension))
    dimension <- length(shape)
  shape <- rep(shape, length.out = dimension)

  ans <- rgamma(use.n * dimension,
                rep(shape, rep(use.n, dimension)))
  dim(ans) <- c(use.n, dimension) 


  ans <- ans / rowSums(ans)
  ans
}




 dirichlet <- function(link = "loge", parallel = FALSE, zero = NULL) {


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(zero) &&
    !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")



  new("vglmff",
  blurb = c("Dirichlet distribution\n\n",
            "Links:     ",
            namesof("shapej", link, earg = earg), "\n\n",
            "Mean:     shape_j/(1 + sum(shape_j)), j = 1,..,ncol(y)"),
  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(1, M, 1), x, .parallel ,
                           constraints, apply.int = TRUE)
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .parallel = parallel, .zero = zero ))),
  initialize = eval(substitute(expression({
    y <- as.matrix(y)
    M <- ncol(y)

    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = Inf,
              out.wy = FALSE,
              colsyperw = NULL,
              maximize = FALSE)

    if (any(y <= 0) || any(y >= 1))
      stop("all y values must be > 0 and < 1")

    mynames1 <- paste("shape", 1:M, sep = "")
    predictors.names <-
      namesof(mynames1, .link , earg = .earg , short = TRUE)
    if (!length(etastart)) {
      yy <- matrix(t(y) %*% rep(1 / nrow(y), nrow(y)), nrow(y), M,
                   byrow = TRUE)
      etastart <- theta2eta(yy, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .link , earg = .earg )
    prop.table(shape, 1)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <- rep( .link , length.out = M)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for(ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$expected <- TRUE
  }), list( .link = link, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    shape <- eta2theta(eta, .link , earg = .earg )
    M <- if (is.matrix(eta)) ncol(eta) else 1
    sumshape <- as.vector(shape %*% rep(1, length.out = M))
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * lgamma(sumshape)) -
      sum(c(w) * lgamma(shape)) +
      sum(c(w) * (shape-1) * log(y))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("dirichlet"),
  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .link , earg = .earg )

    sumshape <- as.vector(shape %*% rep(1, length.out = M))
    dl.dsh <- digamma(sumshape) - digamma(shape) + log(y)

    dsh.deta <- dtheta.deta(shape, .link , earg = .earg )

    c(w) * dl.dsh * dsh.deta
  }), list( .link = link, .earg = earg ))),
  weight = expression({
    index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
    wz <- matrix(trigamma(sumshape), nrow = n, ncol = dimm(M))
    wz[, 1:M] <- wz[, 1:M] - trigamma(shape)
    wz <- -c(w) * wz * dsh.deta[, index$row] * dsh.deta[, index$col]
    wz
  }))
}




 zeta <- function(x, deriv = 0) {



  deriv.arg <- deriv
  rm(deriv)
  if (!is.Numeric(deriv.arg, allowable.length = 1,
                  integer.valued = TRUE))
    stop("'deriv' must be a single non-negative integer")
  if (deriv.arg < 0 || deriv.arg > 2)
    stop("'deriv' must be 0, 1, or 2")


  if (deriv.arg > 0)
    return(Zeta.derivative(x, deriv.arg = deriv.arg))



  if (any(special <- Re(x) <= 1)) {
    ans <- x
    ans[special] <- Inf   # For Re(x) == 1

    special3 <- Re(x) < 1
    ans[special3] <- NA # For 0 < Re(x) < 1

    special4 <- (0 < Re(x)) & (Re(x) < 1) & (Im(x) == 0)
    ans[special4] <- Zeta.derivative(x[special4], deriv.arg = deriv.arg)


    special2 <- Re(x) < 0
    if (any(special2)) {
      x2 <- x[special2]
      cx <- 1-x2
      ans[special2] <- 2^(x2) * pi^(x2-1) * sin(pi*x2/2) *
                      gamma(cx) * Recall(cx)
    }

    if (any(!special)) {
      ans[!special] <- Recall(x[!special])
    }
    return(ans)
  }

  a <- 12; k <- 8
  B <- c(1/6, -1/30,1/42,-1/30,5/66,-691/2730,7/6,-3617/510)
  ans <- 0
  for(ii in 1:(a-1))
     ans <- ans + 1.0 / ii^x
  ans <- ans + 1.0 / ((x-1.0)* a^(x-1.0)) + 1.0 / (2.0 * a^x)

  term <- (x/2) / a^(x+1)
  ans <- ans + term * B[1]

  for(mm in 2:k) {
    term <- term * (x+2*mm-2) * (x+2*mm-3) / (a * a * 2 * mm * (2*mm-1))
    ans <- ans + term * B[mm]
  }
  ans
}



 Zeta.derivative <- function(x, deriv.arg = 0) {


    if (!is.Numeric(deriv.arg, allowable.length = 1,
                    integer.valued = TRUE))
        stop("'deriv.arg' must be a single non-negative integer")
    if (deriv.arg < 0 || deriv.arg > 2)
        stop("'deriv.arg' must be 0, 1, or 2")

    if (any(Im(x) != 0))
        stop("Sorry, currently can only handle x real, not complex")
    if (any(x < 0))
        stop("Sorry, currently cannot handle x < 0")

    ok <- is.finite(x) & x > 0 & x != 1   # Handles NAs
    ans <- rep(as.numeric(NA), length(x))
    nn <- sum(ok)  # Effective length (excludes x < 0 and x = 1 values)
    if (nn)
        ans[ok] <- dotC(name = "vzetawr", as.double(x[ok]), ans = double(nn),
                  as.integer(deriv.arg), as.integer(nn))$ans



    if (deriv.arg == 0)
        ans[is.finite(x) & abs(x) < 1.0e-12] <- -0.5

    ans
}



dzeta <- function(x, p, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


    if (!is.Numeric(p, positive = TRUE)) # || min(p) <= 1
        stop("'p' must be numeric and > 0")
    LLL <- max(length(p), length(x))
    x <- rep(x, length.out = LLL);
    p <- rep(p, length.out = LLL)

    ox <- !is.finite(x)
    zero <- ox | round(x) != x | x < 1
    if (any(zero)) warning("non-integer x and/or x < 1 or NAs")
    ans <- rep(if (log.arg) log(0) else 0, length.out = LLL)
    if (any(!zero)) {
        if (log.arg) {
            ans[!zero] <- (-p[!zero]-1)*log(x[!zero]) - log(zeta(p[!zero]+1))
        } else {
            ans[!zero] <- x[!zero]^(-p[!zero]-1) / zeta(p[!zero]+1)
        }
    }
    if (any(ox)) ans[ox] <- NA
    ans
}


 zetaff <- function(link = "loge", init.p = NULL, zero = NULL) {


  if (length(init.p) && !is.Numeric(init.p, positive = TRUE))
    stop("argument 'init.p' must be > 0")

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  new("vglmff",
  blurb = c("Zeta distribution ",
            "f(y) = 1/(y^(p+1) zeta(p+1)), p>0, y = 1, 2,..\n\n",
            "Link:    ",
            namesof("p", link, earg = earg), "\n\n",
            "Mean:     zeta(p) / zeta(p+1), provided p>1\n",
            "Variance: zeta(p-1) / zeta(p+1) - mean^2, provided p>2"),
  infos = eval(substitute(function(...) {
    list(Musual = 1,
         multipleResponses = TRUE,
         zero = .zero ,
         link = .link)
  }, list( .link = link,
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

    mynames1 <- paste("p", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
      namesof(mynames1, .link , earg = .earg , tag = FALSE)

    Musual <- 1
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    if (!length(etastart)) {
      zetaff.Loglikfun <- function(pp, y, x, w, extraargs) {
        sum(c(w) * dzeta(x = y, p = pp, log = TRUE))
      }


      p.grid <- seq(0.1, 3.0, length.out = 19)
      pp.init <- matrix( if (length( .init.p )) .init.p else -1,
                       n, M, byrow = TRUE)
      if (!length( .init.p ))
      for (spp. in 1:ncoly) {
        pp.init[, spp.] <- getMaxMin(p.grid, objfun = zetaff.Loglikfun,
                                     y = y[, spp.], x = x, w = w[, spp.])
        if ( .link == "loglog")
          pp.init[pp.init <= 1, spp.] <- 1.2
      }

      etastart <- theta2eta(pp.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .init.p = init.p ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    ans <- pp <- eta2theta(eta, .link , earg = .earg )
    ans[pp > 1] <- zeta(pp[pp > 1]) / zeta(pp[pp > 1] + 1)
    ans[pp <= 1] <- NA
    ans
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual

    misc$link <- rep( .link , length = ncoly)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for(ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$multipleResponses <- TRUE
  }), list( .link = link, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    pp <- eta2theta(eta, .link , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dzeta(x = y, p = pp, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("zetaff"),
  deriv = eval(substitute(expression({
    pp <- eta2theta(eta, .link , earg = .earg )

    fred1 <- zeta(pp+1)
    fred2 <- zeta(pp+1, deriv = 1)
    dl.dpp <- -log(y) - fred2 / fred1

    dpp.deta <- dtheta.deta(pp, .link , earg = .earg )

    c(w) * dl.dpp * dpp.deta
  }), list( .link = link, .earg = earg ))),
  weight = expression({
    NOS <- ncol(y)
    nd2l.dpp2 <- zeta(pp + 1, deriv = 2) / fred1 - (fred2/fred1)^2
    wz <- nd2l.dpp2 * dpp.deta^2
    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }))
}



gharmonic <- function(n, s = 1, lognexponent = 0) {

    if (!is.Numeric(n, integer.valued = TRUE, positive = TRUE))
        stop("bad input for argument 'n'")
    if (!is.Numeric(lognexponent, allowable.length = 1))
        stop("bad input for argument 'lognexponent'")
    if (length(n) == 1 && length(s) == 1) {
        if (lognexponent != 0) sum(log(1:n)^lognexponent * (1:n)^(-s)) else
            sum((1:n)^(-s))
    } else {
        LEN <- max(length(n), length(s))
        n <- rep(n, length.out = LEN)
        ans <- s <- rep(s, length.out = LEN)
        if (lognexponent != 0) {
            for(ii in 1:LEN)
                ans[ii] <- sum(log(1:n[ii])^lognexponent * (1:n[ii])^(-s[ii]))
        } else
            for(ii in 1:LEN)
                ans[ii] <- sum((1:n[ii])^(-s[ii]))
        ans
    }
}



rzipf <- function(n, N, s) {
 r <- runif(n)
 sapply(r, function(x) {min(which(pzipf(1:N, N, s) > x))})
}






dzipf <- function(x, N, s, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


    if (!is.Numeric(x))
      stop("bad input for argument 'x'")
    if (!is.Numeric(N, integer.valued = TRUE, positive = TRUE))
      stop("bad input for argument 'N'")
    if (!is.Numeric(s, positive = TRUE))
      stop("bad input for argument 's'")
    nn <- max(length(x), length(N), length(s))
    x <- rep(x, length.out = nn);
    N <- rep(N, length.out = nn);
    s <- rep(s, length.out = nn);
    ox <- !is.finite(x)
    zero <- ox | round(x) != x | x < 1 | x > N
    ans <- (if (log.arg) log(0) else 0) * x
    if (any(!zero))
        if (log.arg) {
          ans[!zero] <- (-s[!zero]) * log(x[!zero]) -
                       log(gharmonic(N[!zero], s[!zero]))
        } else {
          ans[!zero] <- x[!zero]^(-s[!zero]) / gharmonic(N[!zero], s[!zero])
        }
    ans
}



pzipf <- function(q, N, s) {
    if (!is.Numeric(q))
        stop("bad input for argument 'q'")
    if (!is.Numeric(N, integer.valued = TRUE, positive = TRUE))
        stop("bad input for argument 'N'")
    if (!is.Numeric(s, positive = TRUE))
        stop("bad input for argument 's'")

    nn <- max(length(q), length(N), length(s))
    q <- rep(q, length.out = nn);
    N <- rep(N, length.out = nn);
    s <- rep(s, length.out = nn);
    oq <- !is.finite(q)
    zeroOR1 <- oq | q < 1 | q >= N
    floorq <- floor(q)
    ans <- 0 * floorq
    ans[oq | q >= N] <- 1
    if (any(!zeroOR1))
        ans[!zeroOR1] <- gharmonic(floorq[!zeroOR1], s[!zeroOR1]) /
                        gharmonic(N[!zeroOR1], s[!zeroOR1])
    ans
}


 zipf <- function(N = NULL, link = "loge", init.s = NULL) {
  if (length(N) &&
    (!is.Numeric(N, positive = TRUE,
                 integer.valued = TRUE, allowable.length = 1) ||
      N <= 1))
    stop("bad input for argument 'N'")
  enteredN <- length(N)
  if (length(init.s) && !is.Numeric(init.s, positive = TRUE))
      stop("argument 'init.s' must be > 0")

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Zipf distribution f(y;s) = y^(-s) / sum((1:N)^(-s)),",
            " s > 0, y = 1, 2,...,N",
            ifelse(enteredN, paste(" = ",N,sep = ""), ""),
            "\n\n",
            "Link:    ",
            namesof("s", link, earg = earg),
            "\n\n",
            "Mean:    gharmonic(N,s-1) / gharmonic(N,s)"),
  initialize = eval(substitute(expression({


    w.y.check(w = w, y = y,
              Is.integer.y = TRUE)


    predictors.names <- namesof("s", .link , earg = .earg , tag = FALSE)

    NN <- .N
    if (!is.Numeric(NN, allowable.length = 1,
                    positive = TRUE, integer.valued = TRUE))
        NN <- max(y)
    if (max(y) > NN)
        stop("maximum of the response is greater than argument 'N'")
    if (any(y < 1))
        stop("all response values must be in 1, 2, 3,...,N( = ", NN,")")
    extra$N <- NN
    if (!length(etastart)) {
        llfun <- function(ss, y, N, w) {
            sum(c(w) * dzipf(x = y, N=extra$N, s=ss, log = TRUE))
        }
        ss.init <- if (length( .init.s )) .init.s else
            getInitVals(gvals = seq(0.1, 3.0, length.out = 19),
                        llfun=llfun,
                        y = y, N=extra$N, w = w)
        ss.init <- rep(ss.init, length = length(y))
        if ( .link == "loglog") ss.init[ss.init <= 1] = 1.2
        etastart <- theta2eta(ss.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .init.s = init.s, .N = N ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    ss <- eta2theta(eta, .link , earg = .earg )
    gharmonic(extra$N, s=ss - 1) / gharmonic(extra$N, s=ss)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$expected <- FALSE
    misc$link <-    c(s = .link)
    misc$earg <- list(s = .earg )
    misc$N <- extra$N
  }), list( .link = link, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    ss <- eta2theta(eta, .link , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dzipf(x = y, N=extra$N, s=ss, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("zipf"),
  deriv = eval(substitute(expression({
    ss <- eta2theta(eta, .link , earg = .earg )
    fred1 <- gharmonic(extra$N, ss)
    fred2 <- gharmonic(extra$N, ss, lognexp = 1)
    dl.dss <- -log(y) + fred2 / fred1
    dss.deta <- dtheta.deta(ss, .link , earg = .earg )
    d2ss.deta2 <- d2theta.deta2(ss, .link , earg = .earg )
    c(w) * dl.dss * dss.deta
  }), list( .link = link, .earg = earg ))),
  weight = expression({
    d2l.dss <- gharmonic(extra$N, ss, lognexp = 2) / fred1 - (fred2/fred1)^2
    wz <- c(w) * (dss.deta^2 * d2l.dss - d2ss.deta2 * dl.dss)
    wz
  }))
}



cauchy.control <- function(save.weight = TRUE, ...) {
    list(save.weight = save.weight)
}


 cauchy <- function(llocation = "identity", lscale = "loge",
                    ilocation = NULL, iscale = NULL,
                    iprobs = seq(0.2, 0.8, by=0.2),
                    imethod = 1, nsimEIM = NULL, zero = 2) {

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")


  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
      stop("bad input for argument 'zero'")
  if (length(nsimEIM) &&
     (!is.Numeric(nsimEIM, allowable.length = 1, integer.valued = TRUE) ||
      nsimEIM <= 50))
    stop("argument 'nsimEIM' should be an integer greater than 50")
  if (length(iscale) && !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")
  if (!is.Numeric(iprobs, positive = TRUE) || max(iprobs) >= 1)
    stop("bad input for argument 'iprobs'")



  new("vglmff",
  blurb = c("Two parameter Cauchy distribution ",
            "(location & scale unknown)\n\n",
            "Link:    ",
            namesof("location", llocat, earg = elocat), "\n",
            namesof("scale",    lscale,    earg = escale), "\n\n",
            "Mean:     NA\n",
            "Variance: NA"),
  constraints = eval(substitute(expression({
      constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    predictors.names <- c(
      namesof("location", .llocat , earg = .elocat , tag = FALSE),
      namesof("scale",    .lscale ,    earg = .escale ,    tag = FALSE))



    w.y.check(w = w, y = y)



    if (!length(etastart)) {
      loc.init <- if (length( .ilocat)) .ilocat else {
        if ( .imethod == 2) median(rep(y, w)) else 
        if ( .imethod == 3) y else {
            cauchy2.Loglikfun <- function(loc, y, x, w, extraargs) {
                 iprobs <- .iprobs
                 qy <- quantile(rep(y, w), probs = iprobs)
                 ztry <- tan(pi*(iprobs-0.5))
                 btry <- (qy - loc) / ztry
                 scal <- median(btry, na.rm = TRUE)
                 if (scal <= 0) scal <- 0.1
                 sum(c(w) * dcauchy(x = y, loc = loc, scale = scal,
                                    log = TRUE))
             }
             loc.grid <- c(quantile(y, probs = seq(0.1, 0.9, by=0.05)))
             try.this <- getMaxMin(loc.grid, objfun = cauchy2.Loglikfun,
                                  y = y,  x = x, w = w)
                try.this <- rep(c(try.this), length.out = n)
                try.this
            }
        }
        loc.init <- rep(c(loc.init), length.out = n)


            sca.init <- if (length( .iscale )) .iscale else {
                iprobs <- .iprobs
                qy <- quantile(rep(y, w), probs = iprobs)
                ztry <- tan(pi*(iprobs-0.5))
                btry <- (qy - loc.init[1]) / ztry
                sca.init <- median(btry, na.rm = TRUE)
                if (sca.init <= 0) sca.init <- 0.01
                sca.init
            }

            sca.init <- rep(c(sca.init), length.out = n)
            if ( .llocat == "loge") loc.init <- abs(loc.init)+0.01
            etastart <-
              cbind(theta2eta(loc.init, .llocat , earg = .elocat ),
                    theta2eta(sca.init, .lscale ,    earg = .escale ))
        }
    }), list( .ilocat = ilocat,
              .elocat = elocat, .llocat = llocat,
              .iscale = iscale, .escale = escale, .lscale = lscale,
              .iprobs = iprobs, .imethod = imethod ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta[, 1], .llocat , earg = .elocat )
    }, list( .llocat = llocat,
             .elocat = elocat ))),
    last = eval(substitute(expression({
        misc$expected <- TRUE
        misc$link <-    c("location" = .llocat , "scale" =.lscale)
        misc$earg <- list("location" = .elocat , "scale" = .escale )
        misc$imethod <- .imethod
    }), list( .escale = escale, .elocat = elocat,
              .imethod = imethod,
              .llocat = llocat, .lscale = lscale ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
        myscale  <- eta2theta(eta[, 2], .lscale ,    earg = .escale )
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(c(w) * dcauchy(x = y, loc=locat, sc=myscale, log = TRUE))
        }
    }, list( .escale = escale, .lscale = lscale,
             .elocat = elocat, .llocat = llocat ))),
    vfamily = c("cauchy"),
    deriv = eval(substitute(expression({
        location <- eta2theta(eta[, 1], .llocat , earg = .elocat )
        myscale  <- eta2theta(eta[, 2], .lscale , earg = .escale )
        dlocation.deta <- dtheta.deta(location, .llocat , earg = .elocat )
        dscale.deta    <- dtheta.deta(myscale, .lscale , earg = .escale )
        Z <- (y-location) / myscale
        dl.dlocation <- 2 * Z / ((1 + Z^2) * myscale)
        dl.dscale <- (Z^2 - 1) / ((1 + Z^2) * myscale)
        c(w) * cbind(dl.dlocation * dlocation.deta,
                     dl.dscale * dscale.deta)
    }), list( .escale = escale, .lscale = lscale,
              .elocat = elocat, .llocat = llocat ))),
    weight = eval(substitute(expression({
        run.varcov <- 0
        ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
        dthetas.detas = cbind(dlocation.deta, dscale.deta)
        if (length( .nsimEIM )) {
            for(ii in 1:( .nsimEIM )) {
                ysim <- rcauchy(n, loc = location, scale = myscale)
                Z <- (ysim-location) / myscale
                dl.dlocation <- 2 * Z / ((1 + Z^2) * myscale)
                dl.dscale <- (Z^2 - 1) / ((1 + Z^2) * myscale)
                rm(ysim)
                temp3 <- matrix(c(dl.dlocation, dl.dscale), n, 2)
                run.varcov <- ((ii-1) * run.varcov +
                           temp3[, ind1$row.index] *
                           temp3[, ind1$col.index]) / ii
            }
            wz <- if (intercept.only)
                matrix(colMeans(run.varcov),
                       n, ncol(run.varcov), byrow = TRUE) else run.varcov

            wz <- wz * dthetas.detas[, ind1$row] *
                      dthetas.detas[, ind1$col]
            wz <- c(w) * matrix(wz, n, dimm(M))
        } else {
            wz <- cbind(matrix(0.5 / myscale^2,n,2), matrix(0,n,1)) *
                 dthetas.detas[, ind1$row] * dthetas.detas[, ind1$col]
            wz <- c(w) * wz[, 1:M]  # diagonal wz
        }

        wz
    }), list( .escale = escale, .lscale = lscale, .nsimEIM = nsimEIM,
              .elocat = elocat, .llocat = llocat ))))
}







 cauchy1 <- function(scale.arg = 1, llocation = "identity",
                     ilocation = NULL, imethod = 1) {


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation



  if (!is.Numeric(scale.arg, positive = TRUE))
    stop("bad input for 'scale.arg'")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")



  new("vglmff",
  blurb = c("One-parameter Cauchy distribution ",
            "(location unknown, scale known)\n\n",
            "Link:    ",
            namesof("location", llocat, earg = elocat), "\n\n",
            "Mean:     NA\n",
            "Variance: NA"),
  initialize = eval(substitute(expression({
    predictors.names <- namesof("location", .llocat ,
                                earg = .elocat , tag = FALSE)


    w.y.check(w = w, y = y)



        if (!length(etastart)) {
          loc.init <- if (length( .ilocat)) .ilocat else {
            if ( .imethod == 2) median(rep(y, w)) else 
            if ( .imethod == 3) y else {
              cauchy1.Loglikfun <- function(loc, y, x, w, extraargs) {
                 scal <- extraargs
                 sum(c(w) * dcauchy(x = y, loc = loc, scale = scal,
                                    log = TRUE))
               }
               loc.grid <- quantile(y, probs = seq(0.1, 0.9,
                                                  by = 0.05))
                 try.this <- getMaxMin(loc.grid,
                                    objfun = cauchy1.Loglikfun,
                                    y = y,  x = x, w = w,
                                    extraargs = .scale.arg )
              try.this <- rep(try.this, length.out = n)
              try.this
            }
          }
          loc.init <- rep(loc.init, length.out = n)
          if ( .llocat == "loge") loc.init = abs(loc.init)+0.01
          etastart <-
            theta2eta(loc.init, .llocat , earg = .elocat )
        }
    }), list( .scale.arg = scale.arg, .ilocat = ilocat,
              .elocat = elocat, .llocat = llocat,
              .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta, .llocat , earg = .elocat )
  }, list( .llocat = llocat,
           .elocat = elocat ))),
  last = eval(substitute(expression({
    misc$link <-    c("location" = .llocat)
    misc$earg <- list("location" = .elocat )

    misc$expected <- TRUE
    misc$scale.arg <- .scale.arg 
  }), list( .scale.arg = scale.arg, .elocat = elocat,
           .llocat = llocat ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    locat <- eta2theta(eta, .llocat , earg = .elocat )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dcauchy(x = y, loc=locat, scale = .scale.arg,
          log = TRUE))
    }
  }, list( .scale.arg = scale.arg, .elocat = elocat,
           .llocat = llocat ))),
  vfamily = c("cauchy1"),
  deriv = eval(substitute(expression({
    locat <- eta2theta(eta, .llocat , earg = .elocat )
    temp <- (y-locat)/.scale.arg
    dl.dlocat <- 2 * temp / ((1 + temp^2) * .scale.arg)

    dlocation.deta <- dtheta.deta(locat, .llocat , earg = .elocat )

    c(w) * dl.dlocat * dlocation.deta
  }), list( .scale.arg = scale.arg, .elocat = elocat,
            .llocat = llocat ))),
  weight = eval(substitute(expression({
    wz <- c(w) * dlocation.deta^2 / ( .scale.arg^2 * 2)
    wz
  }), list( .scale.arg = scale.arg, .elocat = elocat,
            .llocat = llocat ))))
}






 logistic1 <- function(llocation = "identity",
                       scale.arg = 1, imethod = 1) {
  if (!is.Numeric(scale.arg, allowable.length = 1, positive = TRUE))
    stop("'scale.arg' must be a single positive number")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")



  new("vglmff",
  blurb = c("One-parameter logistic distribution ",
            "(location unknown, scale known)\n\n",
            "Link:    ",
            namesof("location", llocat, earg = elocat), "\n\n",
            "Mean:     location", "\n",
            "Variance: (pi*scale)^2 / 3"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y)


    predictors.names <- namesof("location", .llocat , 
                                earg = .elocat , tag = FALSE)


    if (!length(etastart)) {
      locat.init <- if ( .imethod == 1) y else median(rep(y, w))
      locat.init <- rep(locat.init, length.out = n)
      if ( .llocat == "loge")
        locat.init <- abs(locat.init) + 0.001
      etastart <-
        theta2eta(locat.init, .llocat , earg = .elocat )
    }
  }), list( .imethod = imethod, .llocat = llocat,
            .elocat = elocat ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta, .llocat , earg = .elocat )
  }, list( .llocat = llocat,
           .elocat = elocat ))),
  last = eval(substitute(expression({
    misc$expected <- TRUE
    misc$link <-    c(location = .llocat)
    misc$earg <- list(location = .elocat )
    misc$scale.arg <- .scale.arg 
  }), list( .llocat = llocat, 
            .elocat = elocat, .scale.arg = scale.arg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    locat <- eta2theta(eta, .llocat , earg = .elocat )
    zedd <- (y-locat)/.scale.arg
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dlogis(x = y, locat = locat,
                     scale = .scale.arg, log = TRUE))
    }
  }, list( .llocat = llocat,
           .elocat = elocat, .scale.arg = scale.arg ))),
  vfamily = c("logistic1"),
  deriv = eval(substitute(expression({
    locat <- eta2theta(eta, .llocat , earg = .elocat )

    ezedd <- exp(-(y-locat) / .scale.arg )
    dl.dlocat <- (1 - ezedd) / ((1 + ezedd) * .scale.arg)
    dlocat.deta <- dtheta.deta(locat, .llocat ,
                                 earg = .elocat )

    c(w) * dl.dlocat * dlocat.deta
  }), list( .llocat = llocat,
            .elocat = elocat, .scale.arg = scale.arg ))),
  weight = eval(substitute(expression({
    wz <- c(w) * dlocat.deta^2 / ( .scale.arg^2 * 3) 
    wz
  }), list( .scale.arg = scale.arg ))))
}




 erlang <-
  function(shape.arg, link = "loge",
           imethod = 1, zero = NULL)
{

  if (!is.Numeric(shape.arg, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE))
      stop("'shape' must be a positive integer")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
      stop("argument 'imethod' must be 1 or 2 or 3")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  new("vglmff",
  blurb = c("Erlang distribution\n\n",
            "Link:    ", namesof("scale", link, earg = earg), "\n", "\n",
            "Mean:     shape * scale", "\n",
            "Variance: shape * scale^2"),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 1
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 1,
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    ncoly <- ncol(y)
    Musual <- 1
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    mynames1  <- paste("scale", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
      namesof(mynames1, .link , earg = .earg , tag = FALSE)


    if (!length(etastart)) {
        if ( .imethod == 1) {
          sc.init <- y / .shape.arg
        }
        if ( .imethod == 2) {
          sc.init <- (colSums(y * w) / colSums(w))/ .shape.arg
        }
        if ( .imethod == 3) {
          sc.init <- median(y) / .shape.arg
        }

        if ( !is.matrix(sc.init))
          sc.init <- matrix(sc.init, n, M, byrow = TRUE)


        etastart <-
          theta2eta(sc.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg,
            .shape.arg = shape.arg, .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    sc <- eta2theta(eta, .link , earg = .earg )
    .shape.arg * sc 
  }, list( .link = link, .earg = earg, .shape.arg = shape.arg ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <- c(rep( .link , length = ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for(ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$shape.arg <- .shape.arg 
  }), list( .link = link, .earg = earg, .shape.arg = shape.arg ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    sc <- eta2theta(eta, .link , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * (( .shape.arg - 1) * log(y) - y / sc -
                 .shape.arg * log(sc) - lgamma( .shape.arg )))
    }
  }, list( .link = link, .earg = earg, .shape.arg = shape.arg ))),
  vfamily = c("erlang"),
  deriv = eval(substitute(expression({
    sc <- eta2theta(eta, .link , earg = .earg )
    dl.dsc <- (y / sc - .shape.arg) / sc
    dsc.deta <- dtheta.deta(sc, .link , earg = .earg )
    c(w) * dl.dsc * dsc.deta
  }), list( .link = link, .earg = earg, .shape.arg = shape.arg ))),
  weight = eval(substitute(expression({
    ned2l.dsc2 <- .shape.arg / sc^2
    wz <- c(w) * dsc.deta^2 * ned2l.dsc2
    wz
  }), list( .earg = earg, .shape.arg = shape.arg ))))
}





dbort <- function(x, Qsize = 1, a = 0.5, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(x))
    stop("bad input for argument 'x'")
  if (!is.Numeric(Qsize, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'Qsize'")
  if (!is.Numeric(a, positive = TRUE) || max(a) >= 1)
    stop("bad input for argument 'a'")
  N <- max(length(x), length(Qsize), length(a))
  x <- rep(x, length.out = N);
  Qsize <- rep(Qsize, length.out = N);
  a <- rep(a, length.out = N);

  xok <- (x >= Qsize) & (x == round(x)) & (a > 0) & (a < 1)
  ans <- rep(if (log.arg) log(0) else 0, length.out = N) # loglikelihood
  ans[xok] <- lgamma(1 + Qsize[xok]) - lgamma(x[xok] + 1 - Qsize[xok]) +
             (x[xok] - 1 - Qsize[xok]) * log(x[xok]) +
             (x[xok] - Qsize[xok]) * log(a[xok]) - a[xok] * x[xok]
  if (!log.arg) {
    ans[xok] <- exp(ans[xok])
  }
  ans
}


rbort <- function(n, Qsize = 1, a = 0.5) {

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n
  if (!is.Numeric(Qsize, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'Qsize'")
  if (!is.Numeric(a, positive = TRUE) ||
      max(a) >= 1)
    stop("bad input for argument 'a'")

  N <- use.n
  qsize <- rep(Qsize, length.out = N);
  a <- rep(a, length.out = N)
  totqsize <- qsize
  fini <- (qsize < 1)
  while(any(!fini)) {
    additions <- rpois(sum(!fini), a[!fini])
    qsize[!fini] <- qsize[!fini] + additions
    totqsize[!fini] <- totqsize[!fini] + additions
    qsize <- qsize - 1
    fini <- fini | (qsize < 1)
  }
  totqsize
}



 borel.tanner <- function(Qsize = 1, link = "logit",
                          imethod = 1) {
  if (!is.Numeric(Qsize, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'Qsize'")

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 4)
    stop("argument 'imethod' must be 1 or 2, 3 or 4")




  new("vglmff",
  blurb = c("Borel-Tanner distribution\n\n",
            "Link:    ",
            namesof("a", link, earg = earg), "\n\n",
            "Mean:     Qsize/(1-a)",
            "\n",
            "Variance: Qsize*a / (1-a)^3"),
  initialize = eval(substitute(expression({
    if (any(y < .Qsize ))
      stop("all y values must be >= ", .Qsize )


    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE)


    predictors.names <- namesof("a", .link , earg = .earg , tag = FALSE)

    if (!length(etastart)) {
      a.init <- switch(as.character( .imethod ),
              "1" = 1 - .Qsize / (y+1/8),
              "2" = rep(1 - .Qsize / weighted.mean(y, w), length.out = n),
              "3" = rep(1 - .Qsize / median(y), length.out = n),
              "4" = rep(0.5, length.out = n))
      etastart <-
          theta2eta(a.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .Qsize = Qsize,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    aa <- eta2theta(eta, .link , earg = .earg )
    .Qsize / (1 - aa)
  }, list( .link = link, .earg = earg, .Qsize = Qsize ))),
  last = eval(substitute(expression({
    misc$link <-    c(a = .link)

    misc$earg <- list(a = .earg )

    misc$expected <- TRUE
    misc$Qsize <- .Qsize 
  }), list( .link = link, .earg = earg, .Qsize = Qsize ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    aa <- eta2theta(eta, .link , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dbort(x = y, Qsize =  .Qsize, a = aa, log = TRUE))
    }
  }, list( .link = link, .earg = earg, .Qsize = Qsize ))),
  vfamily = c("borel.tanner"),
  deriv = eval(substitute(expression({
    aa <- eta2theta(eta, .link , earg = .earg )
    dl.da <- (y - .Qsize) / aa - y 
    da.deta <- dtheta.deta(aa, .link , earg = .earg )
    c(w) * dl.da * da.deta
  }), list( .link = link, .earg = earg, .Qsize = Qsize ))),
  weight = eval(substitute(expression({
    ned2l.da2 <- .Qsize / (aa * (1 - aa))
    wz <- c(w) * ned2l.da2 * da.deta^2
    wz
  }), list( .Qsize = Qsize ))))
}





dfelix <- function(x, a = 0.25, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(x))
    stop("bad input for argument 'x'")
  if (!is.Numeric(a, positive = TRUE))
    stop("bad input for argument 'a'")
  N <- max(length(x), length(a))
  x <- rep(x, length.out = N);
  a <- rep(a, length.out = N);

  xok <- (x %% 2 == 1) & (x == round(x)) & (x >= 1) & (a > 0) & (a < 0.5)
  ans <- rep(if (log.arg) log(0) else 0, length.out = N) # loglikelihood
  ans[xok] <- ((x[xok]-3)/2) * log(x[xok]) + ((x[xok]-1)/2) * log(a[xok]) -
             lgamma(x[xok]/2 + 0.5) - a[xok] * x[xok]
  if (!log.arg) {
    ans[xok] <- exp(ans[xok])
  }
  ans
}



 felix <- function(link = elogit(min = 0, max = 0.5), imethod = 1) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
      stop("argument 'imethod' must be 1 or 2, 3 or 4")


  new("vglmff",
  blurb = c("Felix distribution\n\n",
            "Link:    ",
            namesof("a", link, earg = earg), "\n\n",
            "Mean:     1/(1-2*a)"),
  initialize = eval(substitute(expression({
      if (any(y < 1) ||
          any((y+1)/2 != round((y+1)/2)))
        warning("response should be positive, odd and integer-valued")

    w.y.check(w = w, y = y)



      predictors.names <-
        namesof("a", .link , earg = .earg , tag = FALSE)

      if (!length(etastart)) {
          wymean <- weighted.mean(y, w)
          a.init <- switch(as.character( .imethod ),
              "1" = (y - 1 + 1/8) / (2 * (y + 1/8) + 1/8),
              "2" = rep((wymean-1+1/8) / (2*(wymean+1/8)+1/8),
                         length.out = n),
              "3" = rep((median(y)-1+1/8) / (2*(median(y)+1/8)+1/8),
                         length.out = n),
              "4" = rep(0.25,
                         length.out = n))
          etastart <-
            theta2eta(a.init, .link , earg = .earg )
      }
  }), list( .link = link, .earg = earg,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    aa <- eta2theta(eta, .link , earg = .earg )
    1 / (1 - 2 * aa)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$expected <- TRUE

    misc$link <-    c(a = .link)

    misc$earg <- list(a = .earg )
  }), list( .link = link, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    aa <- eta2theta(eta, .link , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
       sum(c(w) * dfelix(x = y, a = aa, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("felix"),
  deriv = eval(substitute(expression({
    aa <- eta2theta(eta, .link , earg = .earg )
    dl.da <- (y - 1) / (2 * aa) - y 
    da.deta <- dtheta.deta(aa, .link , earg = .earg )
    c(w) * dl.da * da.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    ned2l.da2 <- 1 / (aa * (1 - 2 * aa))
    wz <- c(w) * da.deta^2 * ned2l.da2
    wz
  }), list( .link = link ))))
}





 betaff <-
  function(A = 0, B = 1,
           lmu = "logit",
           lphi = "loge",
           imu = NULL, iphi = NULL, imethod = 1, zero = NULL)
{


  stdbeta <- (A == 0 && B == 1)


  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")



  lphi <- as.list(substitute(lphi))
  ephi <- link2list(lphi)
  lphi <- attr(ephi, "function.name")


  if (!is.Numeric(A, allowable.length = 1) ||
      !is.Numeric(B, allowable.length = 1) || A >= B)
    stop("A must be < B, and both must be of length one")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
      stop("bad input for argument 'zero'")

  if (length(imu) && (!is.Numeric(imu, positive = TRUE) ||
     any(imu <= A) || any(imu >= B)))
      stop("bad input for argument 'imu'")
  if (length(iphi) && !is.Numeric(iphi, positive = TRUE))
      stop("bad input for argument 'iphi'")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
      stop("argument 'imethod' must be 1 or 2")


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
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    if (min(y) <= .A || max(y) >= .B)
      stop("data not within (A, B)")


    w.y.check(w = w, y = y)


      predictors.names <- c(namesof("mu",  .lmu ,  .emu , short = TRUE),
                            namesof("phi", .lphi , .ephi, short = TRUE))
      if (!length(etastart)) {
        mu.init <- if (is.Numeric( .imu )) .imu else
                  {if ( .imethod == 1) weighted.mean(y, w) else
                   median(rep(y, w))}
        mu1.init <- (mu.init - .A) / ( .B - .A)  # In (0,1)
        phi.init <- if (is.Numeric( .iphi)) .iphi else
           max(0.01, -1 + ( .B - .A)^2 * mu1.init*(1-mu1.init)/var(y))
        etastart <- matrix(0, n, 2)
        etastart[, 1] <- theta2eta(mu.init, .lmu , earg = .emu )
        etastart[, 2] <- theta2eta(phi.init, .lphi , earg = .ephi )
    }
  }), list( .lmu = lmu, .lphi = lphi, .imu = imu, .iphi = iphi,
            .A = A, .B = B, .emu = emu, .ephi = ephi,
            .imethod = imethod ))),

    linkinv = eval(substitute(function(eta, extra = NULL) {
       mu <- eta2theta(eta[, 1], .lmu , .emu )
       mu
    }, list( .lmu = lmu, .emu = emu, .A = A, .B = B))),
    last = eval(substitute(expression({
        misc$link <-    c(mu = .lmu , phi = .lphi)
        misc$earg <- list(mu = .emu , phi = .ephi)
        misc$limits <- c( .A, .B)
        misc$stdbeta <- .stdbeta
    }), list( .lmu = lmu, .lphi = lphi, .A = A, .B = B,
              .emu = emu, .ephi = ephi,
              .stdbeta = stdbeta ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL){
        mu <- eta2theta(eta[, 1], .lmu , .emu )
        m1u <- if ( .stdbeta ) mu else (mu - .A) / ( .B - .A)
        phi <- eta2theta(eta[, 2], .lphi , .ephi )
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            shape1 <- phi * m1u
            shape2 <- (1 - m1u) * phi
            zedd <- (y - .A) / ( .B - .A)
            sum(c(w) * (dbeta(x = zedd, shape1 = shape1, shape2 = shape2,
                           log = TRUE) -
                     log( abs( .B - .A ))))
        }
    }, list( .lmu = lmu, .lphi = lphi, .A = A, .B = B,
             .emu = emu, .ephi = ephi,
             .stdbeta = stdbeta ))),
    vfamily = "betaff",
    deriv = eval(substitute(expression({
        mu <- eta2theta(eta[, 1], .lmu , .emu )
        phi <- eta2theta(eta[, 2], .lphi , .ephi )
        m1u <- if ( .stdbeta ) mu else (mu - .A) / ( .B - .A)
        dmu.deta <- dtheta.deta(mu, .lmu , .emu )
        dmu1.dmu <- 1 / ( .B - .A)
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
        d2l.dmu12 <- phi^2 * (trigamma(temp1) + trigamma(temp2))
        d2l.dphi2 <- -trigamma(phi) + trigamma(temp1) * m1u^2 +
            trigamma(temp2) * (1-m1u)^2
        d2l.dmu1phi <- temp1*trigamma(temp1) - temp2*trigamma(temp2)
        wz <- matrix(as.numeric(NA), n, dimm(M))
        wz[, iam(1, 1, M)] <- d2l.dmu12 * dmu1.dmu^2 * dmu.deta^2
        wz[, iam(2, 2, M)] <- d2l.dphi2 * dphi.deta^2
        wz[, iam(1, 2, M)] <- d2l.dmu1phi * dmu1.dmu * dmu.deta * dphi.deta
        c(w) * wz
    }), list( .A = A, .B = B ))))
}





 beta.ab <- function(lshape1 = "loge", lshape2 = "loge",
                     i1 = NULL, i2 = NULL, trim = 0.05,
                     A = 0, B = 1, parallel = FALSE, zero = NULL) {

  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")

  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")

  if (length( i1 ) && !is.Numeric( i1, positive = TRUE))
    stop("bad input for argument 'i1'")
  if (length( i2 ) && !is.Numeric( i2, positive = TRUE))
    stop("bad input for argument 'i2'")

  if (!is.Numeric(A, allowable.length = 1) ||
     !is.Numeric(B, allowable.length = 1) ||
     A >= B)
    stop("A must be < B, and both must be of length one")

  stdbeta <- (A == 0 && B == 1) # stdbeta == T iff standard beta distn



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
    constraints <- cm.vgam(matrix(1, M, 1), x, .parallel ,
                           constraints, apply.int  = TRUE)
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .parallel = parallel, .zero = zero ))),
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
      mu1d <- mean(y, trim = .trim)
      uu <- (mu1d - .A) / ( .B - .A) 
      DD <- ( .B - .A)^2 
      pinit <- max(0.01, uu^2 * (1 - uu) * DD / var(y) - uu)
      qinit <- max(0.01, pinit * (1 - uu) / uu)
      etastart <- matrix(0, n, 2)
      etastart[, 1] <- theta2eta( pinit, .lshape1 , earg = .eshape1 )
      etastart[, 2] <- theta2eta( qinit, .lshape2 , earg = .eshape2 )
    }
      if (is.Numeric( .i1 ))
        etastart[, 1] <- theta2eta( .i1, .lshape1 , earg = .eshape1 )
      if (is.Numeric( .i2 ))
        etastart[, 2] <- theta2eta( .i2, .lshape2 , earg = .eshape2 )
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .i1 = i1, .i2 = i2, .trim = trim, .A = A, .B = B,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shapes <- cbind(eta2theta(eta[, 1], .lshape1 , earg = .eshape1 ),
                   eta2theta(eta[, 2], .lshape2 , earg = .eshape2 ))
    .A + ( .B-.A) * shapes[, 1] / (shapes[, 1] + shapes[, 2])
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B, 
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1 = .lshape1 , shape2 = .lshape2)
    misc$earg <- list(shape1 = .eshape1 , shape2 = .eshape2)
    misc$limits <- c( .A, .B)
  }), list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B, 
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL){
    shapes <- cbind(eta2theta(eta[, 1], .lshape1 , earg = .eshape1 ),
                   eta2theta(eta[, 2], .lshape2 , earg = .eshape2 ))
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      zedd <- (y - .A) / ( .B - .A)
      sum(c(w) * (dbeta(x = zedd, shape1 = shapes[, 1],
                        shape2 = shapes[, 2],
                        log = TRUE) - log( abs( .B - .A ))))
    }
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B, 
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  vfamily = "beta.ab",
  deriv = eval(substitute(expression({
    shapes <- cbind(eta2theta(eta[, 1], .lshape1 , earg = .eshape1 ),
                   eta2theta(eta[, 2], .lshape2 , earg = .eshape2 ))

    dshapes.deta <-
      cbind(dtheta.deta(shapes[, 1], .lshape1 , earg = .eshape1),
            dtheta.deta(shapes[, 2], .lshape2 , earg = .eshape2))

    dl.dshapes <- cbind(log(y-.A), log( .B-y)) - digamma(shapes) +
                 digamma(shapes[, 1] + shapes[, 2]) - log( .B - .A)

    c(w) * dl.dshapes * dshapes.deta
  }), list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B, 
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  weight = expression({
    temp2 <- trigamma(shapes[, 1]+shapes[, 2])
    d2l.dshape12 <- temp2 - trigamma(shapes[, 1])
    d2l.dshape22 <- temp2 - trigamma(shapes[, 2])
    d2l.dshape1shape2 <- temp2

    wz <- matrix(as.numeric(NA), n, dimm(M))   #3=dimm(M)
    wz[, iam(1, 1, M)] <- d2l.dshape12 * dshapes.deta[, 1]^2
    wz[, iam(2, 2, M)] <- d2l.dshape22 * dshapes.deta[, 2]^2
    wz[, iam(1, 2, M)] <- d2l.dshape1shape2 * dshapes.deta[, 1] *
                                             dshapes.deta[, 2]

    -c(w) * wz
  }))
}




 simple.exponential <- function() {
  new("vglmff",
  blurb = c("Simple exponential distribution\n",
            "Link:    log(rate)\n"),
  deviance = function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    devy <- -log(y) - 1
    devmu <- -log(mu) - y/mu
    devi <- 2 * (devy - devmu)
    if (residuals) sign(y - mu) * sqrt(abs(devi) * w) else sum(w * devi)
  },
  initialize = expression({
    predictors.names <- "log(rate)"
    mustart <- y + (y == 0) / 8
  }),
  linkinv = function(eta, extra = NULL)
    exp(-eta),
  linkfun = function(mu, extra = NULL)
    -log(mu),
  vfamily = "simple.exponential",
  deriv = expression({
    rate <- 1 / mu
    dl.drate <- mu - y
    drate.deta <- dtheta.deta(rate, "loge")
    c(w) * dl.drate * drate.deta
  }),
  weight = expression({
    ned2l.drate2 <- 1 / rate^2
    wz <- c(w) * drate.deta^2 * ned2l.drate2
    wz
  }))
}




 exponential <- function(link = "loge",
                         location = 0, expected = TRUE,
                         shrinkage.init = 0.95,
                         zero = NULL) {
  if (!is.Numeric(location, allowable.length = 1))
    stop("bad input for argument 'location'")

  if (!is.logical(expected) || length(expected) != 1)
    stop("bad input for argument 'expected'")

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")

  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
      shrinkage.init < 0 ||
      shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")


  new("vglmff",
  blurb = c("Exponential distribution\n\n",
            "Link:     ",
            namesof("rate", link, earg, tag = TRUE), "\n",
            "Mean:     ", "mu = ", 
            if (location == 0) "1 / rate" else
            paste(location, "+ 1 / rate"), "\n",
            "Variance: ",
            if (location == 0) "Exponential: mu^2" else
            paste("(mu - ", location, ")^2", sep = "")),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(Musual = 1,
         zero = .zero )
  }, list( .zero = zero ))),

  deviance = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    devy <- -log(y - .location) - 1
    devmu <- -log(mu - .location) - (y - .location) / (mu - .location)
    devi <- 2 * (devy - devmu)
    if (residuals) {
      sign(y - mu) * sqrt(abs(devi) * w)
    } else  {
      sum(w * devi)
    }
  }, list( .location = location, .earg = earg ))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    ncoly <- ncol(y)
    Musual <- 1
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    extra$Loc <- matrix( .location , n, ncoly, byrow = TRUE)


    if (any(y <= extra$Loc))
      stop("all responses must be greater than ", extra$Loc)

    mynames1 <- if (M == 1) "rate" else paste("rate", 1:M, sep = "")
    predictors.names <-
      namesof(mynames1, .link , earg = .earg , short = TRUE)

    if (length(mustart) + length(etastart) == 0)
      mustart <- matrix(colSums(y * w) / colSums(w),
                        n, M, byrow = TRUE) * .sinit +
                 (1 - .sinit) * y +
                 1 / 8


    if (!length(etastart))
      etastart <- theta2eta(1 / (mustart - extra$Loc),
                            .link , earg = .earg )
  }), list( .location = location,
            .link = link, .earg = earg,
            .sinit = shrinkage.init ))),
  linkinv = eval(substitute(function(eta, extra = NULL)
    extra$Loc + 1 / eta2theta(eta, .link , earg = .earg ),
  list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <- rep( .link , length = M)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for(ii in 1:M) {
      misc$earg[[ii]] <- .earg
    }

    misc$location <- .location
    misc$expected <- .expected
    misc$multipleResponses <- TRUE
    misc$Musual <- Musual
  }), list( .link = link, .earg = earg,
            .expected = expected, .location = location ))),
  linkfun = eval(substitute(function(mu, extra = NULL) 
    theta2eta(1 / (mu - extra$Loc), .link , earg = .earg ),
  list( .link = link, .earg = earg ))),
  vfamily = c("exponential"),
  deriv = eval(substitute(expression({
    rate <- 1 / (mu - extra$Loc)
    dl.drate <- mu - y
    drate.deta <- dtheta.deta(rate, .link , earg = .earg )
    c(w) * dl.drate * drate.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    ned2l.drate2 <- (mu - extra$Loc)^2
    wz <- ned2l.drate2 * drate.deta^2
    if (! .expected ) {
      d2rate.deta2 <- d2theta.deta2(rate, .link , earg = .earg )
      wz <- wz - dl.drate * d2rate.deta2
    }
    c(w) * wz
  }), list( .link = link, .expected = expected, .earg = earg ))))
}




 gamma1 <- function(link = "loge", zero = NULL) {


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")



  new("vglmff",
  blurb = c("1-parameter Gamma distribution\n",
            "Link:     ",
            namesof("shape", link, earg = earg, tag = TRUE), "\n", 
            "Mean:       mu (=shape)\n",
            "Variance:   mu (=shape)"),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 1
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 1,
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


    M <- if (is.matrix(y)) ncol(y) else 1
    Musual <- 1

    mynames1 <- if (M == 1) "shape" else paste("shape", 1:M, sep = "")
    predictors.names <-
      namesof(mynames1, .link , earg = .earg , short = TRUE)

    if (!length(etastart))
      etastart <- cbind(theta2eta(y + 1/8, .link , earg = .earg ))
  }), list( .link = link, .earg = earg ))), 
  linkinv = eval(substitute(function(eta, extra = NULL)
    eta2theta(eta, .link , earg = .earg )),
  list( .link = link, .earg = earg )),
  last = eval(substitute(expression({
    misc$link <- rep( .link , length = M)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$Musual <- Musual
  }), list( .link = link, .earg = earg ))),
  linkfun = eval(substitute(function(mu, extra = NULL)
    theta2eta(mu, .link , earg = .earg )),
  list( .link = link, .earg = earg )),
  loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
    sum(c(w) * dgamma(x = y, shape = mu, scale = 1, log = TRUE))
  },
  vfamily = c("gamma1"),
  deriv = eval(substitute(expression({
    shape <- mu
    dl.dshape <- log(y) - digamma(shape)
    dshape.deta <- dtheta.deta(shape, .link , earg = .earg )
    ans <- c(w) * dl.dshape * dshape.deta
    ans
    c(w) * dl.dshape * dshape.deta
  }), list( .link = link, .earg = earg ))),
  weight = expression({
    ned2l.dshape <- trigamma(shape)
    wz <- ned2l.dshape * dshape.deta^2
    c(w) * wz
  }))
}




 gamma2.ab <-
  function(lrate = "loge", lshape = "loge",
           irate = NULL, ishape = NULL, expected = TRUE, zero = 2)
{

  lrate <- as.list(substitute(lrate))
  erate <- link2list(lrate)
  lrate <- attr(erate, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  if (length( irate) && !is.Numeric(irate, positive = TRUE))
    stop("bad input for argument 'irate'")
  if (length( ishape) && !is.Numeric(ishape, positive = TRUE))
    stop("bad input for argument 'ishape'")

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")

  if (!is.logical(expected) || length(expected) != 1)
    stop("bad input for argument 'expected'")




  new("vglmff",
  blurb = c("2-parameter Gamma distribution\n",
            "Links:    ",
            namesof("rate",  lrate,  earg = erate), ", ", 
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     mu = shape/rate\n",
            "Variance: (mu^2)/shape = shape/rate^2"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              Is.positive.y = TRUE)


      predictors.names <-
        c(namesof("rate",  .lrate  , earg = .erate  , tag = FALSE),
          namesof("shape", .lshape , earg = .eshape , tag = FALSE))

      if (!length(etastart)) {
        mymu = y + 0.167 * (y == 0)
        junk <- lsfit(x, y, wt = w, intercept = FALSE)
        var.y.est <- sum(c(w) * junk$resid^2) / (nrow(x) - length(junk$coef))
        init.shape <-  if (length( .ishape )) .ishape else mymu^2 / var.y.est
        init.rate <-  if (length( .irate)) .irate else init.shape / mymu
        init.rate <- rep(init.rate, length.out = n)
        init.shape <- rep(init.shape, length.out = n)
        if ( .lshape == "loglog")
          init.shape[init.shape <= 1] <- 3.1 # Hopefully value is big enough
          etastart <-
            cbind(theta2eta(init.rate,  .lrate  , earg = .erate ),
                  theta2eta(init.shape, .lshape , earg = .eshape ))
      }
  }), list( .lrate = lrate, .lshape = lshape,
            .irate = irate, .ishape = ishape,
            .erate = erate, .eshape = eshape))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 2], .lshape , earg = .eshape ) / (
    eta2theta(eta[, 1], .lrate  , earg = .erate ))
  }, list( .lrate = lrate, .lshape = lshape,
           .erate = erate, .eshape = eshape))),
  last = eval(substitute(expression({
    misc$link <-    c(rate = .lrate , shape = .lshape)
    misc$earg <- list(rate = .erate, shape = .eshape )
  }), list( .lrate = lrate, .lshape = lshape,
            .erate = erate, .eshape = eshape))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    rate  <- eta2theta(eta[, 1], .lrate  , earg = .erate )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
            sum(c(w) * dgamma(x = y, shape = shape, rate=rate, log = TRUE))
      }
  }, list( .lrate = lrate, .lshape = lshape,
           .erate = erate, .eshape = eshape))),
  vfamily = c("gamma2.ab"),
  deriv = eval(substitute(expression({
      rate <- eta2theta(eta[, 1], .lrate , earg = .erate )
      shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
      dl.drate <- mu - y
      dl.dshape <- log(y*rate) - digamma(shape)
      dratedeta <- dtheta.deta(rate, .lrate , earg = .erate )
      dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
      c(w) * cbind(dl.drate * dratedeta,
                   dl.dshape * dshape.deta)
  }), list( .lrate = lrate, .lshape = lshape,
            .erate = erate, .eshape = eshape))),
  weight = eval(substitute(expression({
    d2l.dshape2 <- -trigamma(shape)
    d2l.drate2 <- -shape/(rate^2)
    d2l.drateshape <- 1/rate
    wz <- matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
    wz[, iam(1, 1, M)] <- -d2l.drate2 * dratedeta^2
    wz[, iam(2, 2, M)] <- -d2l.dshape2 * dshape.deta^2
    wz[, iam(1, 2, M)] <- -d2l.drateshape * dratedeta * dshape.deta
    if (! .expected) {
      d2ratedeta2 <- d2theta.deta2(rate, .lrate , earg = .erate )
      d2shapedeta2 <- d2theta.deta2(shape, .lshape , earg = .eshape )
      wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] - dl.drate * d2ratedeta2
      wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] - dl.dshape * d2shapedeta2
    }
    c(w) * wz
  }), list( .lrate = lrate, .lshape = lshape,
            .erate = erate, .eshape = eshape, .expected = expected ))))
}



 gamma2 <-
  function(lmu = "loge", lshape = "loge",
           imethod = 1,  ishape = NULL,
           parallel = FALSE, apply.parint = FALSE,
           deviance.arg = FALSE, zero = -2)
{




  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  if (length(zero) && !is.Numeric(zero, integer.valued = TRUE))
    stop("bad input for argument 'zero'")

  if (length( ishape) && !is.Numeric(ishape, positive = TRUE))
    stop("bad input for argument 'ishape'")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  if (!is.logical(apply.parint) ||
      length(apply.parint) != 1)
    stop("argument 'apply.parint' must be a single logical")


  if (is.logical(parallel) && parallel && length(zero))
    stop("set 'zero = NULL' if 'parallel = TRUE'")


    ans <- 
    new("vglmff",
    blurb = c("2-parameter Gamma distribution",
              " (McCullagh and Nelder 1989 parameterization)\n",
              "Links:    ",
              namesof("mu",    lmu,    earg = emu), ", ", 
              namesof("shape", lshape, earg = eshape), "\n",
              "Mean:     mu\n",
              "Variance: (mu^2)/shape"),
    constraints = eval(substitute(expression({

    constraints <- cm.vgam(matrix(1, M, 1), x, .parallel , constraints,
                           apply.int = .apply.parint )

        dotzero <- .zero
        Musual <- 2
        eval(negzero.expression)
        constraints <- cm.zero.vgam(constraints, x, z_Index, M)
  }), list( .zero = zero,
            .parallel = parallel, .apply.parint = apply.parint ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    Musual <- 2

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


      assign("CQO.FastAlgorithm", ( .lmu == "loge" && .lshape == "loge"),
             envir = VGAM:::VGAMenv)
      if (any(function.name == c("cqo","cao")) &&
         is.Numeric( .zero , allowable.length = 1) && .zero != -2)
        stop("argument zero = -2 is required")

      M <- Musual * ncol(y)
      NOS <- ncoly <- ncol(y)  # Number of species


      temp1.names =
        if (NOS == 1) "mu"    else paste("mu",    1:NOS, sep = "")
      temp2.names =
        if (NOS == 1) "shape" else paste("shape", 1:NOS, sep = "")
      predictors.names <-
          c(namesof(temp1.names, .lmu ,    earg = .emu ,    tag = FALSE),
            namesof(temp2.names, .lshape , earg = .eshape, tag = FALSE))
      predictors.names <- predictors.names[interleave.VGAM(M, M = Musual)]




    if (is.logical( .parallel ) & .parallel & ncoly > 1)
      warning("the constraint matrices may not be correct with ",
              "multiple responses")



      if (!length(etastart)) {
        init.shape <- matrix(1.0, n, NOS)
        mymu <- y # + 0.167 * (y == 0)  # imethod == 1 (the default)
        if ( .imethod == 2) {
                for(ii in 1:ncol(y)) {
                    mymu[, ii] <- weighted.mean(y[, ii], w = w[, ii])
                }
        }
        for(spp in 1:NOS) {
          junk <- lsfit(x, y[, spp], wt = w[, spp], intercept = FALSE)
          var.y.est <- sum(w[, spp] * junk$resid^2) / (n - length(junk$coef))
          init.shape[, spp] <- if (length( .ishape )) .ishape else
              mymu[, spp]^2 / var.y.est
          if ( .lshape == "loglog")
              init.shape[init.shape[, spp] <= 1,spp] <- 3.1
        }
        etastart <-
              cbind(theta2eta(mymu, .lmu , earg = .emu ),
                    theta2eta(init.shape, .lshape , earg = .eshape ))
        etastart <-
            etastart[, interleave.VGAM(M, M = Musual), drop = FALSE]
      }
  }), list( .lmu = lmu, .lshape = lshape, .ishape = ishape,
            .emu = emu, .eshape = eshape,
            .parallel = parallel, .apply.parint = apply.parint,
            .zero = zero, .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Musual <- 2
    NOS <- ncol(eta) / Musual
    eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
              .lmu , earg = .emu )
  }, list( .lmu = lmu, .emu = emu ))),
  last = eval(substitute(expression({
    if (exists("CQO.FastAlgorithm", envir = VGAM:::VGAMenv))
        rm("CQO.FastAlgorithm", envir = VGAM:::VGAMenv)

    tmp34 <- c(rep( .lmu ,    length = NOS),
              rep( .lshape , length = NOS))
    names(tmp34) =
       c(if (NOS == 1) "mu"    else paste("mu",    1:NOS, sep = ""), 
         if (NOS == 1) "shape" else paste("shape", 1:NOS, sep = ""))
    tmp34 <- tmp34[interleave.VGAM(M, M = 2)]
    misc$link <- tmp34 # Already named

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:NOS) {
      misc$earg[[Musual*ii-1]] <- .emu
      misc$earg[[Musual*ii  ]] <- .eshape
    }

    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$parallel <- .parallel
    misc$apply.parint <- .apply.parint
  }), list( .lmu = lmu, .lshape = lshape,
            .emu = emu, .eshape = eshape,
            .parallel = parallel, .apply.parint = apply.parint ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    temp <- theta2eta(mu, .lmu , earg = .emu )
    temp <- cbind(temp, NA * temp)
    temp[, interleave.VGAM(ncol(temp), M = 2), drop = FALSE]
  }, list( .lmu = lmu, .emu = emu ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Musual <- 2
    NOS <- ncol(eta) / Musual
    mymu <- mu  # eta2theta(eta[, 2*(1:NOS)-1], .lmu , earg = .emu )
    shapemat <- eta2theta(eta[, Musual * (1:NOS), drop = FALSE],
                         .lshape , earg = .eshape )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
        sum(c(w) * dgamma(x = y,
                       shape = c(shapemat),
                       scale = c(mymu/shapemat),
                       log = TRUE))
    }
  }, list( .lmu = lmu, .lshape = lshape,
           .emu = emu, .eshape = eshape))),
  vfamily = c("gamma2"),
  deriv = eval(substitute(expression({
    Musual <- 2
    NOS <- ncol(eta) / Musual

    mymu  <- eta2theta(eta[, Musual * (1:NOS) - 1],
                      .lmu ,    earg = .emu    )
    shape <- eta2theta(eta[, Musual * (1:NOS)],
                      .lshape , earg = .eshape )

    dl.dmu <- shape * (y / mymu - 1) / mymu
    dl.dshape <- log(y) + log(shape) - log(mymu) + 1 - digamma(shape) -
                y / mymu

    dmu.deta    <- dtheta.deta(mymu,  .lmu ,    earg = .emu )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    myderiv <- c(w) * cbind(dl.dmu    * dmu.deta,
                           dl.dshape * dshape.deta)
    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list( .lmu = lmu, .lshape = lshape,
            .emu = emu, .eshape = eshape))),
  weight = eval(substitute(expression({
    ned2l.dmu2 <- shape / (mymu^2)
    ned2l.dshape2 <- trigamma(shape) - 1 / shape
    wz <- matrix(as.numeric(NA), n, M) # 2 = M; diagonal!

    wz[, Musual*(1:NOS)-1] <- ned2l.dmu2 * dmu.deta^2
    wz[, Musual*(1:NOS)  ] <- ned2l.dshape2 * dshape.deta^2


    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)

  }), list( .lmu = lmu ))))



  if (deviance.arg) ans@deviance <- eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {


    if (ncol(as.matrix(y)) > 1 && ncol(as.matrix(w)) > 1)
      stop("cannot handle matrix 'w' yet")


    Musual <- 2
    NOS <- ncol(eta) / 2
    temp300 <-  eta[, 2*(1:NOS), drop = FALSE]
    shape <-  eta2theta(temp300, .lshape , earg = .eshape )
    devi <- -2 * (log(y/mu) - y/mu + 1)
    if (residuals) {
       warning("not 100% sure about these deviance residuals!")
       sign(y - mu) * sqrt(abs(devi) * w)
    } else
       sum(c(w) * devi)
  }, list( .lshape = lshape )))
  ans
}



 geometric <- function(link = "logit", expected = TRUE,
                       imethod = 1, iprob = NULL, zero = NULL) {

  if (!is.logical(expected) || length(expected) != 1)
    stop("bad input for argument 'expected'")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


 if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")



  new("vglmff",
  blurb = c("Geometric distribution ",
            "(P[Y=y] = prob * (1 - prob)^y, y = 0, 1, 2,...)\n",
            "Link:     ",
            namesof("prob", link, earg = earg), "\n",
            "Mean:     mu = (1 - prob) / prob\n",
            "Variance: mu * (1 + mu) = (1 - prob) / prob^2"),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 1
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 1,
         zero = .zero )
  }, list( .zero = zero ))),


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
    Musual <- 1
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    mynames1  <- paste("prob", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
      namesof(mynames1, .link , earg = .earg , tag = FALSE)


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
            .imethod = imethod, .iprob = iprob ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob <- eta2theta(eta, .link , earg = .earg )
    (1 - prob) / prob 
  }, list( .link = link, .earg = earg ))),

  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <- c(rep( .link , length = ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for(ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$expected <- .expected
    misc$imethod <- .imethod
    misc$iprob <- .iprob
  }), list( .link = link, .earg = earg,
            .iprob = iprob,
            .expected = expected, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    prob <- eta2theta(eta, .link , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
        sum(c(w) * dgeom(x = y, prob = prob, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("geometric"),
  deriv = eval(substitute(expression({
    prob <- eta2theta(eta, .link , earg = .earg )

    dl.dprob <- -y / (1 - prob) + 1 / prob 

    dprobdeta <- dtheta.deta(prob, .link , earg = .earg )
    c(w) * cbind(dl.dprob * dprobdeta)
  }), list( .link = link, .earg = earg, .expected = expected ))),
  weight = eval(substitute(expression({
    ned2l.dprob2 <- if ( .expected ) {
      1 / (prob^2 * (1 - prob))
    } else {
      y / (1 - prob)^2 + 1 / prob^2
    }
    wz <- ned2l.dprob2 * dprobdeta^2
    if ( !( .expected ))
      wz <- wz - dl.dprob * d2theta.deta2(prob, .link , earg = .earg )
    c(w) * wz
  }), list( .link = link, .earg = earg,
            .expected = expected ))))
}




dbetageom <- function(x, shape1, shape2, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(x))
    stop("bad input for argument 'x'")
  if (!is.Numeric(shape1, positive = TRUE))
    stop("bad input for argument 'shape1'")
  if (!is.Numeric(shape2, positive = TRUE))
    stop("bad input for argument 'shape2'")
  N <- max(length(x), length(shape1), length(shape2))
  x <- rep(x, length.out = N);
  shape1 <- rep(shape1, length.out = N);
  shape2 <- rep(shape2, length.out = N)

  loglik <- lbeta(1+shape1, shape2+abs(x)) - lbeta(shape1, shape2)
  xok <- (x == round(x) & x >= 0)
  loglik[!xok] <- log(0)
  if (log.arg) {
    loglik
  } else {
    exp(loglik)
  }
}


pbetageom <- function(q, shape1, shape2, log.p = FALSE) {
    if (!is.Numeric(q))
      stop("bad input for argument 'q'")
    if (!is.Numeric(shape1, positive = TRUE))
      stop("bad input for argument 'shape1'")
    if (!is.Numeric(shape2, positive = TRUE))
      stop("bad input for argument 'shape2'")
    N <- max(length(q), length(shape1), length(shape2))
    q <- rep(q, length.out = N);
    shape1 <- rep(shape1, length.out = N);
    shape2 <- rep(shape2, length.out = N)
    ans <- q * 0  # Retains names(q)
    if (max(abs(shape1-shape1[1])) < 1.0e-08 &&
       max(abs(shape2-shape2[1])) < 1.0e-08) {
        qstar <- floor(q)
        temp <- if (max(qstar) >= 0) dbetageom(x = 0:max(qstar), 
               shape1 = shape1[1], shape2 = shape2[1]) else 0*qstar
        unq <- unique(qstar)
        for(ii in unq) {
            index <- (qstar == ii)
            ans[index] <- if (ii >= 0) sum(temp[1:(1+ii)]) else 0
        }
    } else
    for(ii in 1:N) {
        qstar <- floor(q[ii])
        ans[ii] <- if (qstar >= 0) sum(dbetageom(x = 0:qstar,
                 shape1 = shape1[ii], shape2 = shape2[ii])) else 0
    }
    if (log.p) log(ans) else ans
}


rbetageom <- function(n, shape1, shape2) {
  rgeom(n = n, prob = rbeta(n = n, shape1 = shape1, shape2 = shape2))
}








negbinomial.control <- function(save.weight = TRUE, ...) {
    list(save.weight = save.weight)
}



 negbinomial <- function(lmu = "loge", lsize = "loge",
                         imu = NULL,   isize = NULL,
                         probs.y = 0.75,
                         nsimEIM = 100, cutoff = 0.995, Maxiter = 5000,
                         deviance.arg = FALSE, imethod = 1,
                         parallel = FALSE,
                         shrinkage.init = 0.95, zero = -2) {













  lmuuu <- as.list(substitute(lmu))
  emuuu <- link2list(lmuuu)
  lmuuu <- attr(emuuu, "function.name")
  imuuu <- imu

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")


  if (length(imuuu) && !is.Numeric(imuuu, positive = TRUE))
    stop("bad input for argument 'imu'")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("bad input for argument 'isize'")

  if (!is.Numeric(cutoff, allowable.length = 1) ||
    cutoff < 0.8 ||
    cutoff >= 1)
    stop("range error in the argument 'cutoff'")
  if (!is.Numeric(Maxiter, integer.valued = TRUE, allowable.length = 1) ||
    Maxiter < 100)
    stop("bad input for argument 'Maxiter'")
  if (!is.Numeric(imethod, allowable.length = 1,
    integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")
  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
    shrinkage.init < 0 ||
     shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")

  if (!is.null(nsimEIM)) {
    if (!is.Numeric(nsimEIM, allowable.length = 1, integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
    if (nsimEIM <= 10)
      warning("argument 'nsimEIM' should be an integer ",
               "greater than 10, say")
  }


    if (!is.logical( parallel ) || length( parallel ) != 1)
      stop("argument 'parallel' must be TRUE or FALSE")
    if ( parallel  && length(zero))
      stop("need to set 'zero = NULL' when parallel = TRUE")



  ans <- 
  new("vglmff",




  blurb = c("Negative-binomial distribution\n\n",
            "Links:    ",
            namesof("mu",   lmuuu, earg = emuuu), ", ",
            namesof("size", lsize, earg = esize), "\n",
            "Mean:     mu\n",
            "Variance: mu * (1 + mu / size) for NB-2"),

  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)

    if ( .parallel && ncol(cbind(y)) > 1)
      stop("univariate responses needed if 'parallel = TRUE'")
    constraints <- cm.vgam(matrix(1, M, 1), x, .parallel , constraints)
  }), list( .parallel = parallel, .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero)
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({
    Musual <- 2

    temp5 <- w.y.check(w = w, y = y,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1, maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    assign("CQO.FastAlgorithm",
          ( .lmuuu == "loge") && ( .lsize == "loge"),
           envir = VGAM:::VGAMenv)

    if (any(function.name == c("cqo", "cao")) &&
        is.Numeric( .zero , allowable.length = 1) &&
        .zero != -2)
        stop("argument zero = -2 is required")


    if (any(y < 0))
      stop("negative values not allowed for the 'negbinomial' family")
    if (any(round(y) != y))
      stop("integer-values only allowed for the 'negbinomial' family")
    if (ncol(w) > ncol(y))
      stop("number of columns of prior-'weights' is greater than ",
           "the number of responses")

    M <- Musual * ncol(y) 
    NOS <- ncoly <- ncol(y) # Number of species
    predictors.names <-
     c(namesof(if (NOS == 1) "mu"   else paste("mu",   1:NOS, sep = ""),
                .lmuuu, earg = .emuuu, tag = FALSE),
       namesof(if (NOS == 1) "size" else paste("size", 1:NOS, sep = ""),
                .lsize , earg = .esize , tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M = Musual)]

    if (is.null( .nsimEIM )) {
       save.weight <- control$save.weight <- FALSE
    }

    if (is.numeric( .mu.init ))
      MU.INIT <- matrix( .mu.init , nrow(y), ncol(y), byrow = TRUE)


    if (!length(etastart)) {
      mu.init <- y
      for(iii in 1:ncol(y)) {
        use.this <- if ( .imethod == 1) {
          weighted.mean(y[, iii], w[, iii]) + 1/16
        } else if ( .imethod == 3) {
          c(quantile(y[, iii], probs = .probs.y ) + 1/16)
        } else {
          median(y[, iii]) + 1/16
        }

        if (is.numeric( .mu.init )) {
          mu.init[, iii] <- MU.INIT[, iii]
        } else {
          medabsres <- median(abs(y[, iii] - use.this)) + 1/32
          allowfun <- function(z, maxtol = 1) sign(z) * pmin(abs(z), maxtol)
          mu.init[, iii] <- use.this + (1 - .sinit ) *
                           allowfun(y[, iii] - use.this, maxtol = medabsres)

          mu.init[, iii] <- abs(mu.init[, iii]) + 1 / 1024
        }
      } # of for(iii)

      if ( is.Numeric( .k.init )) {
        kay.init <- matrix( .k.init, nrow = n, ncol = NOS, byrow = TRUE)
      } else {
        negbinomial.Loglikfun <- function(kmat, y, x, w, extraargs) {
            mu <- extraargs
            sum(c(w) * dnbinom(x = y, mu = mu, size = kmat, log = TRUE))
        }
        k.grid <- 2^((-7):7)
        k.grid <- 2^(seq(-8, 8, length = 40))
        kay.init <- matrix(0, nrow = n, ncol = NOS)
        for(spp. in 1:NOS) {
          kay.init[, spp.] <- getMaxMin(k.grid,
                                       objfun = negbinomial.Loglikfun,
                                       y = y[, spp.], x = x, w = w[, spp.],
                                       extraargs = mu.init[, spp.])
        }
      }



    newemu <- .emuuu
    if ( .lmuuu == "nbcanlink") {
      newemu$size <- kay.init
    }




      etastart <-
        cbind(theta2eta(mu.init , link = .lmuuu , earg = newemu ),
              theta2eta(kay.init, link = .lsize , earg = .esize ))
      etastart <-
        etastart[, interleave.VGAM(M, M = Musual), drop = FALSE]
      }
  }), list( .lmuuu = lmuuu, .lsize = lsize,
            .emuuu = emuuu, .esize = esize,
            .mu.init = imu,
            .k.init = isize, .probs.y = probs.y,
            .sinit = shrinkage.init, .nsimEIM = nsimEIM,
            .zero = zero, .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Musual <- 2
    NOS <- ncol(eta) / Musual
    eta.k <- eta[, Musual * (1:NOS) , drop = FALSE]
    kmat <- eta2theta(eta.k, .lsize , earg = .esize )




    newemu <- .emuuu
    if ( .lmuuu == "nbcanlink") {
      newemu$size <- kmat
    }



    eta2theta(eta[, Musual * (1:NOS) - 1, drop = FALSE], .lmuuu ,
              earg = newemu)
  }, list( .lmuuu = lmuuu, .lsize = lsize,
           .emuuu = emuuu, .esize = esize))),

  last = eval(substitute(expression({
    if (exists("CQO.FastAlgorithm", envir = VGAM:::VGAMenv))
        rm("CQO.FastAlgorithm", envir = VGAM:::VGAMenv)

    temp0303 <- c(rep( .lmuuu, length = NOS),
                 rep( .lsize , length = NOS))
    names(temp0303) =
      c(if (NOS == 1) "mu"   else paste("mu",   1:NOS, sep = ""),
        if (NOS == 1) "size" else paste("size", 1:NOS, sep = ""))
    temp0303 <- temp0303[interleave.VGAM(M, M = 2)]
    misc$link <- temp0303 # Already named

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:NOS) {
      misc$earg[[Musual*ii-1]] <- newemu
      misc$earg[[Musual*ii  ]] <- .esize
    }

    misc$cutoff <- .cutoff 
    misc$imethod <- .imethod 
    misc$nsimEIM <- .nsimEIM
    misc$expected <- TRUE
    misc$shrinkage.init <- .sinit
    misc$multipleResponses <- TRUE
  }), list( .lmuuu = lmuuu, .lsize = lsize,
            .emuuu = emuuu, .esize = esize,
            .cutoff = cutoff,
            .nsimEIM = nsimEIM,
            .sinit = shrinkage.init,
            .imethod = imethod ))),

  linkfun = eval(substitute(function(mu, extra = NULL) {
    Musual <- 2

    newemu <- .emuuu

    eta.temp <- theta2eta(mu, .lmuuu , earg = newemu)
    eta.kayy <- theta2eta(if (is.numeric( .isize )) .isize else 1.0,
                     .lsize , earg = .esize )
    eta.kayy <- 0 * eta.temp + eta.kayy  # Right dimension now.





    if ( .lmuuu == "nbcanlink") {
      newemu$size <- eta2theta(eta.kayy, .lsize , earg = .esize )
    }



    eta.temp <- cbind(eta.temp, eta.kayy)
    eta.temp[, interleave.VGAM(ncol(eta.temp), M = Musual), drop = FALSE]
  }, list( .lmuuu = lmuuu, .lsize = lsize,
           .emuuu = emuuu, .esize = esize,
                           .isize = isize ))),

  loglikelihood = eval(substitute(
      function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      Musual <- 2
      NOS <- ncol(eta) / Musual

      eta.k <- eta[, Musual*(1:NOS), drop = FALSE]
      if ( .lsize == "loge") {
          bigval <- 68
          eta.k <- ifelse(eta.k >  bigval,  bigval, eta.k)
          eta.k <- ifelse(eta.k < -bigval, -bigval, eta.k)
      }
      kmat <- eta2theta(eta.k, .lsize , earg = .esize )





    newemu <- .emuuu
    if ( .lmuuu == "nbcanlink") {
      newemu$size <- kmat
    }





    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
      sum(c(w) * dnbinom(x = y, mu = mu, size = kmat, log = TRUE))
  }, list( .lsize = lsize,
           .lmuuu = lmuuu, .emuuu = emuuu, .esize = esize))),

  vfamily = c("negbinomial"),

  deriv = eval(substitute(expression({
    Musual <- 2
    NOS <- ncol(eta) / Musual
    M <- ncol(eta)
    eta.k <- eta[, Musual*(1:NOS)  , drop = FALSE]
    if ( .lsize == "loge") {
      bigval <- 68
      eta.k <- ifelse(eta.k >  bigval,  bigval, eta.k)
      eta.k <- ifelse(eta.k < -bigval, -bigval, eta.k)
    }
    kmat <- eta2theta(eta.k, .lsize , earg = .esize )




    newemu <- .emuuu
    if ( .lmuuu == "nbcanlink") {
      newemu$size <- kmat
    }



    dl.dmu <- y / mu - (y + kmat) / (mu + kmat)
    dl.dk <- digamma(y + kmat) - digamma(kmat) -
            (y + kmat) / (mu + kmat) + 1 + log(kmat / (kmat + mu))

    if ( .lmuuu == "nbcanlink")
      newemu$wrt.eta <- 1
    dmu.deta <- dtheta.deta(mu, .lmuuu , earg = newemu) # eta1

    if ( .lmuuu == "nbcanlink")
      newemu$wrt.eta <- 2
    dk.deta1 <- dtheta.deta(mu, .lmuuu , earg = newemu) # eta2

    dk.deta2 <- dtheta.deta(kmat, .lsize , earg = .esize)

    myderiv <- c(w) * cbind(dl.dmu * dmu.deta,
                           dl.dk  * dk.deta2)


    if ( .lmuuu == "nbcanlink") {
      myderiv[, 1:NOS] =
      myderiv[, 1:NOS] + c(w) * dl.dk * dk.deta1
    }




    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list( .lmuuu = lmuuu, .lsize = lsize,
            .emuuu = emuuu, .esize = esize))),

  weight = eval(substitute(expression({
    wz <- matrix(as.numeric(NA), n, M)


    if (is.null( .nsimEIM )) {
      fred2 <- dotFortran(name = "enbin9", ans = double(n*NOS),
                  as.double(kmat), as.double(mu), as.double( .cutoff ),
                  as.integer(n), ok = as.integer(1), as.integer(NOS),
                  sumpdf = double(1), as.double( .Machine$double.eps ),
                  as.integer( .Maxiter ))
      if (fred2$ok != 1)
        stop("error in Fortran subroutine exnbin9")
      dim(fred2$ans) <- c(n, NOS)
      ned2l.dk2 <- -fred2$ans - 1/kmat + 1/(kmat+mu)
      wz[, Musual*(1:NOS)] <- dk.deta2^2 * ned2l.dk2



    } else {

      run.varcov <- matrix(0, n, NOS)

      for(ii in 1:( .nsimEIM )) {
        ysim <- rnbinom(n = n*NOS, mu = c(mu), size = c(kmat))
        if (NOS > 1) dim(ysim) = c(n, NOS)
        dl.dk <- digamma(ysim + kmat) - digamma(kmat) -
                (ysim + kmat) / (mu + kmat) +
                1 + log(kmat / (kmat + mu))
        run.varcov <- run.varcov + dl.dk^2
      } # end of for loop

      run.varcov <- cbind(run.varcov / .nsimEIM )
      ned2l.dk2 <- if (intercept.only)
          matrix(colMeans(run.varcov),
                 n, ncol(run.varcov), byrow = TRUE) else run.varcov

      wz[, Musual*(1:NOS)] <- ned2l.dk2 * dk.deta2^2
    } # end of else


    ned2l.dmu2 <- 1 / mu - 1 / (mu + kmat)
    wz[, Musual*(1:NOS) - 1] <- ned2l.dmu2 * dmu.deta^2



    if ( .lmuuu == "nbcanlink") {
      wz[, Musual*(1:NOS)-1] <-
      wz[, Musual*(1:NOS)-1] + ned2l.dk2 * dk.deta1^2

      wz <- cbind(wz,
                 kronecker(ned2l.dk2 * dk.deta1 * dk.deta2,
                           if (NOS > 1) cbind(1, 0) else 1))
    }




    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }), list( .cutoff = cutoff,
            .Maxiter = Maxiter,
            .lmuuu = lmuuu,
            .nsimEIM = nsimEIM ))))



  if (deviance.arg) ans@deviance = eval(substitute(
      function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      Musual <- 2
      NOS <- ncol(eta) / Musual



    if (ncol(as.matrix(y)) > 1 && ncol(as.matrix(w)) > 1)
      stop("cannot handle matrix 'w' yet")


      temp300 <-  eta[, Musual*(1:NOS), drop = FALSE]
      if ( .lsize == "loge") {
          bigval <- 68
          temp300[temp300 >  bigval] <-  bigval
          temp300[temp300 < -bigval] <- -bigval
      } else stop("can only handle the 'loge' link")
      kmat <-  eta2theta(temp300, .lsize , earg = .esize )
      devi <- 2 * (y * log(ifelse(y < 1, 1, y) / mu) +
             (y + kmat) * log((mu + kmat) / (kmat + y)))
      if (residuals) {
         sign(y - mu) * sqrt(abs(devi) * w)
      } else {
         sum(c(w) * devi)
      }
  }, list( .lsize = lsize, .emuuu = emuuu,
           .esize = esize)))

  ans
}







polya.control <- function(save.weight = TRUE, ...) {
    list(save.weight = save.weight)
}



 polya <-
  function(lprob = "logit", lsize = "loge",
           iprob = NULL,    isize = NULL,
           probs.y = 0.75,
           nsimEIM = 100,
           deviance.arg = FALSE, imethod = 1,
           shrinkage.init = 0.95, zero = -2) {



  if (length(iprob) && !is.Numeric(iprob, positive = TRUE))
    stop("bad input for argument 'iprob'")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("bad input for argument 'isize'")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
     stop("argument 'imethod' must be 1 or 2 or 3")
  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
     shrinkage.init < 0 ||
     shrinkage.init > 1)
     stop("bad input for argument 'shrinkage.init'")

  if (!is.Numeric(nsimEIM, allowable.length = 1,
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



  ans =
  new("vglmff",
  blurb = c("Polya (negative-binomial) distribution\n\n",
            "Links:    ",
            namesof("prob", lprob, earg = eprob), ", ",
            namesof("size", lsize, earg = esize), "\n",
            "Mean:     size * (1 - prob) / prob\n",
            "Variance: mean / prob"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)

  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero)
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({
    Musual <- 2
    if (any(function.name == c("cqo", "cao")))
      stop("polya() does not work with cqo() or cao(). ",
           "Try negbinomial()")



    temp5 <- w.y.check(w = w, y = y,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1, maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    M <- Musual * ncol(y)
    NOS <- ncoly <- ncol(y)  # Number of species

    predictors.names <-
      c(namesof(if (NOS == 1) "prob" else
                paste("prob", 1:NOS, sep = ""),
               .lprob , earg = .eprob , tag = FALSE),
        namesof(if (NOS == 1) "size" else
                paste("size", 1:NOS, sep = ""),
               .lsize ,  earg = .esize ,  tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M = 2)]

    if (is.null( .nsimEIM )) {
       save.weight <- control$save.weight <- FALSE
    }

    
    PROB.INIT <- if (is.numeric( .pinit )) {
      matrix( .pinit, nrow(y), ncol(y), byrow = TRUE)
    } else {
      NULL
    }

    if (!length(etastart)) {
      mu.init <- y
      for(iii in 1:ncol(y)) {
        use.this <- if ( .imethod == 1) {
          weighted.mean(y[, iii], w[, iii]) + 1/16
        } else if ( .imethod == 3) {
          c(quantile(y[, iii], probs <- .probs.y) + 1/16)
        } else {
          median(y[, iii]) + 1/16
        }

        if (FALSE) {
          mu.init[, iii] <- MU.INIT[, iii]
        } else {
          medabsres <- median(abs(y[, iii] - use.this)) + 1/32
          allowfun <- function(z, maxtol = 1) sign(z) * pmin(abs(z), maxtol)
          mu.init[, iii] <- use.this + (1 - .sinit) * allowfun(y[, iii] -
                          use.this, maxtol = medabsres)

          mu.init[, iii] <- abs(mu.init[, iii]) + 1 / 1024
        }
      }



      if ( is.Numeric( .kinit )) {
        kayy.init <- matrix( .kinit, nrow = n, ncol = NOS, byrow = TRUE)
      } else {
        negbinomial.Loglikfun <- function(kmat, y, x, w, extraargs) {
            mu <- extraargs
            sum(c(w) * dnbinom(x = y, mu = mu, size = kmat, log = TRUE))
        }
        k.grid <- 2^((-7):7)
        k.grid <- 2^(seq(-8, 8, length = 40))
        kayy.init <- matrix(0, nrow = n, ncol = NOS)
        for(spp. in 1:NOS) {
          kayy.init[, spp.] <- getMaxMin(k.grid,
                             objfun = negbinomial.Loglikfun,
                             y = y[, spp.], x = x, w = w,
                             extraargs = mu.init[, spp.])
        }
      }

      prob.init <- if (length(PROB.INIT)) PROB.INIT else
                  kayy.init / (kayy.init + mu.init)


      etastart <-
        cbind(theta2eta(prob.init, .lprob , earg = .eprob),
              theta2eta(kayy.init, .lsize , earg = .esize))
      etastart <-
        etastart[, interleave.VGAM(M, M = Musual), drop = FALSE]
      }
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize,
            .pinit = iprob, .kinit = isize,
            .probs.y = probs.y,
            .sinit = shrinkage.init, .nsimEIM = nsimEIM, .zero = zero,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Musual <- 2
    NOS <- ncol(eta) / Musual
    pmat <- eta2theta(eta[, Musual*(1:NOS) - 1, drop = FALSE],
                     .lprob , earg = .eprob)
    kmat <- eta2theta(eta[, Musual*(1:NOS)-  0, drop = FALSE],
                     .lsize , earg = .esize)
    kmat / (kmat + pmat)
  }, list( .lprob = lprob, .eprob = eprob,
           .lsize = lsize, .esize = esize))),
  last = eval(substitute(expression({
    temp0303 <- c(rep( .lprob , length = NOS),
                 rep( .lsize , length = NOS))
    names(temp0303) =
      c(if (NOS == 1) "prob" else paste("prob", 1:NOS, sep = ""),
        if (NOS == 1) "size" else paste("size", 1:NOS, sep = ""))
    temp0303 <- temp0303[interleave.VGAM(M, M = 2)]
    misc$link <- temp0303 # Already named

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:NOS) {
      misc$earg[[Musual*ii-1]] <- .eprob
      misc$earg[[Musual*ii  ]] <- .esize
    }

    misc$isize <- .isize  
    misc$imethod <- .imethod 
    misc$nsimEIM <- .nsimEIM
    misc$expected <- TRUE
    misc$shrinkage.init <- .sinit
    misc$Musual <- 2
    misc$multipleResponses <- TRUE
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize,
            .isize = isize,
            .nsimEIM = nsimEIM,
            .sinit = shrinkage.init, .imethod = imethod ))),


  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Musual <- 2
    NOS <- ncol(eta) / Musual
    pmat  <- eta2theta(eta[, Musual*(1:NOS) - 1, drop = FALSE],
                      .lprob , earg = .eprob)
    temp300 <-         eta[, Musual*(1:NOS)    , drop = FALSE]
    if ( .lsize == "loge") {
      bigval <- 68
      temp300 <- ifelse(temp300 >  bigval,  bigval, temp300)
      temp300 <- ifelse(temp300 < -bigval, -bigval, temp300)
    }
    kmat <- eta2theta(temp300, .lsize , earg = .esize)
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
      sum(c(w) * dnbinom(x = y, prob = pmat, size = kmat, log = TRUE))
  }, list( .lsize = lsize, .lprob = lprob,
           .esize = esize, .eprob = eprob ))),
  vfamily = c("polya"),
  deriv = eval(substitute(expression({
    Musual <- 2
    NOS <- ncol(eta) / Musual
    M <- ncol(eta)

    pmat  <- eta2theta(eta[, Musual*(1:NOS) - 1, drop = FALSE],
                      .lprob , earg = .eprob)
    temp3 <-           eta[, Musual*(1:NOS)    , drop = FALSE]
    if ( .lsize == "loge") {
      bigval <- 68
      temp3 <- ifelse(temp3 >  bigval,  bigval, temp3)
      temp3 <- ifelse(temp3 < -bigval, -bigval, temp3)
    }
    kmat <- eta2theta(temp3, .lsize , earg = .esize)

    dl.dprob <- kmat / pmat - y / (1.0 - pmat)
    dl.dkayy <- digamma(y + kmat) - digamma(kmat) + log(pmat)

    dprob.deta <- dtheta.deta(pmat, .lprob , earg = .eprob)
    dkayy.deta <- dtheta.deta(kmat, .lsize , earg = .esize)
    dthetas.detas <- cbind(dprob.deta, dkayy.deta)
    dThetas.detas <- dthetas.detas[, interleave.VGAM(M, M = Musual)]
    myderiv <- c(w) * cbind(dl.dprob, dl.dkayy) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize))),
  weight = eval(substitute(expression({
    wz <- matrix(0.0, n, M + M - 1) # wz is 'tridiagonal' 

    ind1 <- iam(NA, NA, M = Musual, both = TRUE, diag = TRUE)
    mumat <- as.matrix(mu)


    for(spp. in 1:NOS) {
      run.varcov <- 0
      kvec <- kmat[, spp.]
      pvec <- pmat[, spp.]

      for(ii in 1:( .nsimEIM )) {
        ysim <- rnbinom(n = n, prob = pvec, size = kvec)

        dl.dprob <- kvec / pvec - ysim / (1.0 - pvec)
        dl.dkayy <- digamma(ysim + kvec) - digamma(kvec) + log(pvec)
        temp3 <- cbind(dl.dprob, dl.dkayy)
        run.varcov <- run.varcov +
                     temp3[, ind1$row.index] *
                     temp3[, ind1$col.index]
      }
      run.varcov <- cbind(run.varcov / .nsimEIM)

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



    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }), list( .nsimEIM = nsimEIM ))))




  if (deviance.arg) ans@deviance = eval(substitute(
      function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Musual <- 2
    NOS <- ncol(eta) / Musual
    temp300 <-  eta[, Musual*(1:NOS), drop = FALSE]




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
    if (residuals)
      sign(y - mu) * sqrt(abs(devi) * w) else
      sum(c(w) * devi)
    }, list( .lsize = lsize, .eprob = eprob,
             .esize = esize)))

    ans
} # End of polya()




 simple.poisson <- function() {
  new("vglmff",
  blurb = c("Poisson distribution\n\n",
            "Link:     log(lambda)",
            "\n",
            "Variance: lambda"),
  deviance = function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    nz <- y > 0
    devi <-  - (y - mu)
    devi[nz] <- devi[nz] + y[nz] * log(y[nz]/mu[nz])
    if (residuals) sign(y - mu) * sqrt(2 * abs(devi) * w) else
        2 * sum(c(w) * devi)
  },
  initialize = expression({
    if (ncol(cbind(w)) != 1)
      stop("prior weight must be a vector or a one-column matrix")

    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")

    predictors.names <- "log(lambda)"

    mu <- (weighted.mean(y, w) + y) / 2 + 1/8

    if (!length(etastart))
      etastart <- log(mu)
  }), 
  linkinv = function(eta, extra = NULL)
    exp(eta),
  last = expression({
    misc$link <-    c(lambda = "loge")
    misc$earg <- list(lambda = list())
  }),
  link = function(mu, extra = NULL)
    log(mu),
  vfamily = "simple.poisson",
  deriv = expression({
    lambda <- mu
    dl.dlambda <- -1 + y/lambda
    dlambda.deta <- dtheta.deta(theta = lambda, link = "loge")
    c(w) * dl.dlambda * dlambda.deta
  }),
  weight = expression({
    d2l.dlambda2 <- 1 / lambda
    c(w) * d2l.dlambda2 * dlambda.deta^2
  }))
}













 studentt <-  function(ldf = "loglog", idf = NULL,
                       tol1 = 0.1, imethod = 1)
{


  ldof <- as.list(substitute(ldf))
  edof <- link2list(ldof)
  ldof <- attr(edof, "function.name")
  idof <- idf


  if (length(idof))
    if (!is.Numeric(idof) || any(idof <= 1))
      stop("argument 'idf' should be > 1")

  if (!is.Numeric(tol1, positive  = TRUE))
    stop("argument 'tol1' should be positive")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
      stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Student t-distribution\n\n",
            "Link:     ",
            namesof("df", ldof, earg = edof), "\n",
            "Variance: df / (df - 2) if df > 2\n"),
  infos = eval(substitute(function(...) {
    list(Musual = 1,
         tol1 = .tol1 )
  }, list( .tol1 = tol1 ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y)


    predictors.names <- namesof("df", .ldof , earg = .edof , tag = FALSE)

    if (!length(etastart)) {

      init.df <- if (length( .idof )) .idof else {
        VarY <- var(y)
        MadY <- mad(y)
        if (VarY <= (1 + .tol1 )) VarY <- 1.12
        if ( .imethod == 1) {
          2 * VarY / (VarY - 1)
        } else if ( .imethod == 2) {
          ifelse(MadY < 1.05, 30, ifelse(MadY > 1.2, 2, 5))
        } else
          10
      }


      etastart <- rep(theta2eta(init.df, .ldof , earg = .edof ),
                      length.out = length(y))
    }
  }), list( .ldof = ldof, .edof = edof, .idof = idof,
            .tol1 = tol1, .imethod = imethod ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Dof <- eta2theta(eta, .ldof , earg = .edof )
    ans <- 0 * eta
    ans[Dof <= 1] <- NA
    ans
  }, list( .ldof = ldof, .edof = edof ))),
  last = eval(substitute(expression({
    misc$link <-    c(df = .ldof )
    misc$earg <- list(df = .edof )
    misc$imethod <- .imethod
    misc$expected = TRUE
  }), list( .ldof = ldof,
            .edof = edof, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Dof <-  eta2theta(eta, .ldof , earg = .edof )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
        sum(c(w) * dt(x = y, df = Dof, log = TRUE))
    }
  }, list( .ldof = ldof, .edof = edof ))), 
  vfamily = c("studentt"),
  deriv = eval(substitute(expression({
    Dof <- eta2theta(eta, .ldof , earg = .edof )
    ddf.deta <-  dtheta.deta(theta = Dof, .ldof , earg = .edof )

    DDS  <- function(df)          digamma((df + 1) / 2) -  digamma(df / 2)
    DDSp <- function(df)  0.5 * (trigamma((df + 1) / 2) - trigamma(df / 2))

    temp0 <- 1 / Dof
    temp1 <-  temp0 * y^2
    dl.ddf <- 0.5 * (-temp0 - log1p(temp1) +
              (Dof + 1) * y^2 / (Dof^2 * (1 + temp1)) + DDS(Dof))
    c(w) * dl.ddf * ddf.deta
  }), list( .ldof = ldof, .edof = edof ))),
  weight = eval(substitute(expression({

    const2 <- (Dof + 0) / (Dof + 3)
    const2[!is.finite(Dof)] <- 1  # Handles Inf

    tmp6 <- DDS(Dof)
    nedl2.dnu2 <- 0.5 * (tmp6 * (const2 * tmp6 - 2 / (Dof + 1)) - DDSp(Dof))
 
    wz <- c(w) * nedl2.dnu2 * ddf.deta^2
    wz
  }), list( .ldof = ldof, .edof = edof ))))
}






    Kayfun.studentt <- function(df, bigno = .Machine$double.eps^(-0.46)) {
      ind1 <- is.finite(df)

      const4 <- dnorm(0)
      ans <- df

      if (any(ind1))
        ans[ind1] <- exp(lgamma((df[ind1] + 1) / 2) -
                         lgamma( df[ind1]      / 2)) / sqrt(pi * df[ind1])
      ans[df <= 0] <- NaN
      ind2 <- (df >= bigno)
      if (any(ind2)) {
        dff <- df[ind2]
        ans[ind2] <- const4 # 1 / const3  # for handling df = Inf
      }
      ans[!ind1] <- const4 # 1 / const3  # for handling df = Inf

      ans
    }




 studentt3 <- function(llocation = "identity",
                       lscale    = "loge",
                       ldf       = "loglog",
                       ilocation = NULL, iscale = NULL, idf = NULL,
                       imethod = 1,
                       zero = -(2:3)) {



  lloc <- as.list(substitute(llocation))
  eloc <- link2list(lloc)
  lloc <- attr(eloc, "function.name")

  lsca <- as.list(substitute(lscale))
  esca <- link2list(lsca)
  lsca <- attr(esca, "function.name")

  ldof <- as.list(substitute(ldf))
  edof <- link2list(ldof)
  ldof <- attr(edof, "function.name")


  iloc <- ilocation
  isca <- iscale
  idof <- idf
 

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
      stop("argument 'imethod' must be 1 or 2 or 3")

  if (length(iloc))
    if (!is.Numeric(iloc))
      stop("bad input in argument 'ilocation'")
  if (length(isca))
    if (!is.Numeric(isca, positive = TRUE))
      stop("argument 'iscale' should be positive")
  if (length(idof))
    if (!is.Numeric(idof) || any(idof <= 1))
      stop("argument 'idf' should be > 1")



  new("vglmff",
  blurb = c("Student t-distribution\n\n",
            "Link:     ",
            namesof("location", lloc, earg = eloc), ", ",
            namesof("scale",    lsca, earg = esca), ", ",
            namesof("df",       ldof, earg = edof), "\n",
            "Variance: scale^2 * df / (df - 2) if df > 2\n"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 3
    eval(negzero.expression)
  }), list( .zero = zero ))),
    infos = eval(substitute(function(...) {
      list(Musual = 3,
           zero = .zero)
    }, list( .zero = zero ))),
  initialize = eval(substitute(expression({
    Musual <- 3



    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$Musual <- Musual
    M <- Musual * ncoly #

    mynames1 <- paste("location", if (NOS > 1) 1:NOS else "", sep = "")
    mynames2 <- paste("scale",    if (NOS > 1) 1:NOS else "", sep = "")
    mynames3 <- paste("df",       if (NOS > 1) 1:NOS else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lloc , earg = .eloc , tag = FALSE),
          namesof(mynames2, .lsca , earg = .esca , tag = FALSE),
          namesof(mynames3, .ldof , earg = .edof , tag = FALSE))
    predictors.names <-
      predictors.names[interleave.VGAM(Musual * NOS, M = Musual)]

    if (!length(etastart)) {
      init.loc <- if (length( .iloc )) .iloc else {
        if ( .imethod == 2) apply(y, 2, median) else
        if ( .imethod == 3) (colMeans(y) + t(y)) / 2 else {
           colSums(w * y) / colSums(w)
        }
      }

      sdvec <- apply(y, 2, sd)
      init.sca <- if (length( .isca )) .isca else
                  sdvec / 2.3

      sdvec    <- rep(sdvec,
                      length.out <- max(length(sdvec),
                                       length(init.sca)))
      init.sca <- rep(init.sca,
                      length.out <- max(length(sdvec),
                                       length(init.sca)))
      ind9 <- (sdvec / init.sca <= (1 + 0.12))
      sdvec[ind9] <- sqrt(1.12) * init.sca[ind9]
      init.dof <- if (length( .idof )) .idof else
                (2 * (sdvec / init.sca)^2) / ((sdvec / init.sca)^2  - 1)
      if (!is.Numeric(init.dof) || init.dof <= 1)
        init.dof <- rep(3, length.out = ncoly)

      mat1 <- matrix(theta2eta(init.loc, .lloc , earg = .eloc ), n, NOS,
                     byrow = TRUE)
      mat2 <- matrix(theta2eta(init.sca, .lsca , earg = .esca ), n, NOS,
                     byrow = TRUE)
      mat3 <- matrix(theta2eta(init.dof, .ldof , earg = .edof ), n, NOS,
                     byrow = TRUE)
      etastart <- cbind(mat1, mat2, mat3)
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = Musual)]
    }
  }), list( .lloc = lloc, .eloc = eloc, .iloc = iloc,
            .lsca = lsca, .esca = esca, .isca = isca,
            .ldof = ldof, .edof = edof, .idof = idof,
            .imethod = imethod ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS    <- extra$NOS
    Musual <- extra$Musual
    Loc <-  eta2theta(eta[, Musual*(1:NOS)-2], .lloc , earg = .eloc )
    Dof <-  eta2theta(eta[, Musual*(1:NOS)-0], .ldof , earg = .edof )
    Loc[Dof <= 1] <- NA
    Loc
  }, list( .lloc = lloc, .eloc = eloc,
           .lsca = lsca, .esca = esca,
           .ldof = ldof, .edof = edof ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <- c(rep( .lloc , length = NOS),
                   rep( .lsca , length = NOS),
                   rep( .ldof , length = NOS))
    misc$link <- misc$link[interleave.VGAM(Musual * NOS, M = Musual)]
    temp.names <- c(mynames1, mynames2, mynames3)
    temp.names <- temp.names[interleave.VGAM(Musual * NOS, M = Musual)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", Musual * NOS)
    names(misc$earg) <- temp.names
    for(ii in 1:NOS) {
      misc$earg[[Musual*ii-2]] <- .eloc
      misc$earg[[Musual*ii-1]] <- .esca
      misc$earg[[Musual*ii  ]] <- .edof
    }
 
    misc$Musual <- Musual
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .ldof = ldof, .edof = edof,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    NOS <- extra$NOS
    Musual <- extra$Musual
    Loc <-  eta2theta(eta[, Musual*(1:NOS)-2], .lloc , earg = .eloc )
    Sca <-  eta2theta(eta[, Musual*(1:NOS)-1], .lsca , earg = .esca )
    Dof <-  eta2theta(eta[, Musual*(1:NOS)-0], .ldof , earg = .edof )
    zedd <- (y - Loc) / Sca
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
        sum(c(w) * (dt(x = zedd, df = Dof, log = TRUE) - log(Sca)))
    }
  }, list(  .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .ldof = ldof, .edof = edof ))), 
  vfamily = c("studentt3"),
  deriv = eval(substitute(expression({
    Musual <- extra$Musual
    NOS <- extra$NOS
    Loc <- eta2theta(eta[, Musual*(1:NOS)-2], .lloc , earg = .eloc )
    Sca <- eta2theta(eta[, Musual*(1:NOS)-1], .lsca , earg = .esca )
    Dof <- eta2theta(eta[, Musual*(1:NOS)-0], .ldof , earg = .edof )

    dloc.deta <- cbind(dtheta.deta(theta = Loc, .lloc , earg = .eloc ))
    dsca.deta <- cbind(dtheta.deta(theta = Sca, .lsca , earg = .esca ))
    ddof.deta <- cbind(dtheta.deta(theta = Dof, .ldof , earg = .edof ))

    zedd  <- (y - Loc) / Sca
    temp0 <- 1 / Dof
    temp1 <- temp0 * zedd^2
    dl.dloc <- (Dof + 1) * zedd / (Sca * (Dof + zedd^2))
    dl.dsca <- zedd * dl.dloc - 1 / Sca
    dl.ddof <- 0.5 * (-temp0 - log1p(temp1) +
                     (Dof+1) * zedd^2 / (Dof^2 * (1 + temp1)) +
                     digamma((Dof+1)/2) - digamma(Dof/2))
 
    ans <- c(w) * cbind(dl.dloc * dloc.deta,
                        dl.dsca * dsca.deta,
                        dl.ddof * ddof.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M = Musual)]
    ans
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .ldof = ldof, .edof = edof ))),
  weight = eval(substitute(expression({

    const1 <- (Dof + 1) / (Dof + 3)
    const2 <- (Dof + 0) / (Dof + 3)
    const1[!is.finite(Dof)] <- 1  # Handles Inf
    const2[!is.finite(Dof)] <- 1  # Handles Inf

    const4 <- dnorm(0)
    ned2l.dlocat2 <-      const1 / (Sca * (Kayfun.studentt(Dof) / const4))^2
    ned2l.dscale2 <- 2  * const2 /  Sca^2

    DDS  <- function(df)          digamma((df + 1) / 2) -  digamma(df / 2)
    DDSp <- function(df)  0.5 * (trigamma((df + 1) / 2) - trigamma(df / 2))


    tmp6 <- DDS(Dof)
    edl2.dnu2 <- 0.5 * (tmp6 * (const2 * tmp6 - 2 / (Dof + 1)) - DDSp(Dof))
    ned2l.dshape2 <- cbind(edl2.dnu2)  # cosmetic name change

    ned2l.dshape.dlocat <- cbind(0 * Sca)
    ned2l.dshape.dscale <- cbind((-1 / (Dof + 1.0) + const2 * DDS(Dof)) / Sca)



    wz <- array(c(c(w) * ned2l.dlocat2 * dloc.deta^2,
                  c(w) * ned2l.dscale2 * dsca.deta^2,
                  c(w) * ned2l.dshape2 * ddof.deta^2,
                  c(w) * ned2l.dshape2 * 0,
                  c(w) * ned2l.dshape.dscale  * dsca.deta * ddof.deta,
                  c(w) * ned2l.dshape.dlocat * dloc.deta * ddof.deta),
                dim = c(n, M / Musual, 6))
    wz <- arwz2wz(wz, M = M, Musual = Musual)



 if (FALSE) {
    wz <- matrix(0.0, n, dimm(M))
    wz[, Musual*(1:NOS) - 2] <- ned2l.dlocat2 * dloc.deta^2
    wz[, Musual*(1:NOS) - 1] <- ned2l.dscale2 * dsca.deta^2
    wz[, Musual*(1:NOS) - 0] <- ned2l.dshape2 * ddof.deta^2

    for (ii in ((1:NOS) - 1)) {
      ind3 <- 1 + ii
      wz[, iam(ii*Musual + 1, ii*Musual + 3, M = M)] <-
           ned2l.dshape.dlocat[, ind3] *
           dloc.deta[, ind3] * ddof.deta[, ind3]
      wz[, iam(ii*Musual + 2, ii*Musual + 3, M = M)] <-
           ned2l.dshape.dscale[, ind3] *
           dsca.deta[, ind3] * ddof.deta[, ind3]
    }

  while (all(wz[, ncol(wz)] == 0))
    wz <- wz[, -ncol(wz)]
 }



    wz
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .ldof = ldof, .edof = edof ))))
}





 studentt2 <- function(df = Inf,
                       llocation = "identity",
                       lscale    = "loge",
                       ilocation = NULL, iscale = NULL,
                       imethod = 1,
                       zero = -2) {

  lloc <- as.list(substitute(llocation))
  eloc <- link2list(lloc)
  lloc <- attr(eloc, "function.name")

  lsca <- as.list(substitute(lscale))
  esca <- link2list(lsca)
  lsca <- attr(esca, "function.name")



  iloc <- ilocation; isca <- iscale
  doff <- df


  if (is.finite(doff))
    if (!is.Numeric(doff, positive = TRUE))
    stop("argument 'df' must be positive")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
      stop("argument 'imethod' must be 1 or 2 or 3")

  if (length(iloc))
    if (!is.Numeric(iloc))
      stop("bad input in argument 'ilocation'")
  if (length(isca))
    if (!is.Numeric(isca, positive = TRUE))
      stop("argument 'iscale' should be positive")


  new("vglmff",
  blurb = c("Student t-distribution\n\n",
            "Link:     ",
            namesof("location", lloc, earg = eloc), ", ",
            namesof("scale",    lsca, earg = esca), "\n",
            "Variance: scale^2 * df / (df - 2) if df > 2\n"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),
    infos = eval(substitute(function(...) {
      list(Musual = 2,
           zero = .zero)
    }, list( .zero = zero ))),
  initialize = eval(substitute(expression({
    Musual <- 2


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    extra$NOS <- NOS <- ncoly <- ncol(y) # Number of species
    extra$Musual <- Musual
    M <- Musual * ncoly #

    mynames1 <- paste("location", if (NOS > 1) 1:NOS else "", sep = "")
    mynames2 <- paste("scale",    if (NOS > 1) 1:NOS else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lloc , earg = .eloc , tag = FALSE),
          namesof(mynames2, .lsca , earg = .esca , tag = FALSE))
    predictors.names <-
      predictors.names[interleave.VGAM(Musual * NOS, M = Musual)]

    if (!length(etastart)) {

      init.loc <- if (length( .iloc )) .iloc else {
        if ( .imethod == 2) apply(y, 2, median) else
        if ( .imethod == 3) (colMeans(y) + t(y)) / 2 else {
           colSums(w * y) / colSums(w)
        }
      }

      sdvec <- apply(y, 2, sd)
      init.sca <- if (length( .isca )) .isca else
                  sdvec / 2.3

      mat1 <- matrix(theta2eta(init.loc, .lloc , earg = .eloc ), n, NOS,
                     byrow = TRUE)
      mat2 <- matrix(theta2eta(init.sca, .lsca , earg = .esca ), n, NOS,
                     byrow = TRUE)
      etastart <- cbind(mat1, mat2)
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = Musual)]
    }
  }), list( .lloc = lloc, .eloc = eloc, .iloc = iloc,
            .lsca = lsca, .esca = esca, .isca = isca,
            .doff = doff,
            .imethod = imethod ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- extra$NOS
    Musual <- extra$Musual
    Loc <-  eta2theta(eta[, Musual*(1:NOS) - 1], .lloc , earg = .eloc )
    Dof <- matrix( .doff , nrow(cbind(Loc)), NOS, byrow = TRUE)
    Loc[Dof <= 1] <- NA
    Loc
  }, list( .lloc = lloc, .eloc = eloc,
           .lsca = lsca, .esca = esca,
           .doff = doff ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <- c(rep( .lloc , length = NOS),
                   rep( .lsca , length = NOS))
    temp.names <- c(mynames1, mynames2)
    temp.names <- temp.names[interleave.VGAM(Musual * NOS, M = Musual)]
    names(misc$link) <- temp.names
    misc$earg <- vector("list", Musual * NOS)
    names(misc$earg) <- temp.names
    for(ii in 1:NOS) {
      misc$earg[[Musual*ii-1]] <- .eloc
      misc$earg[[Musual*ii-0]] <- .esca
    }
 
    misc$Musual <- Musual
    misc$simEIM <- TRUE
    misc$df <- .doff
    misc$imethod <- .imethod
    misc$expected = TRUE
    misc$multipleResponses <- TRUE
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .doff = doff,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    NOS <- extra$NOS
    Musual <- extra$Musual
    Loc <- eta2theta(eta[, Musual*(1:NOS)-1], .lloc , earg = .eloc )
    Sca <- eta2theta(eta[, Musual*(1:NOS)-0], .lsca , earg = .esca )
    Dof <- matrix( .doff , nrow(cbind(Loc)), NOS, byrow = TRUE)
    zedd <- (y - Loc) / Sca
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
        sum(c(w) * (dt(x = zedd, df = Dof, log = TRUE) - log(Sca)))
    }
  }, list(  .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .doff = doff ))), 
  vfamily = c("studentt2"),
  deriv = eval(substitute(expression({
    Musual <- extra$Musual
    NOS <- extra$NOS
    Loc <- eta2theta(eta[, Musual*(1:NOS)-1], .lloc , earg = .eloc )
    Sca <- eta2theta(eta[, Musual*(1:NOS)-0], .lsca , earg = .esca )
    Dof <- matrix( .doff , n, NOS, byrow = TRUE)

    dlocat.deta <- dtheta.deta(theta = Loc, .lloc , earg = .eloc )
    dscale.deta <- dtheta.deta(theta = Sca, .lsca , earg = .esca )

    zedd  <- (y - Loc) / Sca
    temp0 <- 1 / Dof
    temp1 <- temp0 * zedd^2
    dl.dlocat <- (Dof + 1) * zedd / (Sca * (Dof + zedd^2))
    dl.dlocat[!is.finite(Dof)] <- zedd / Sca  # Adjust for df=Inf
    dl.dscale <- zedd * dl.dlocat - 1 / Sca
 
    ans <- c(w) * cbind(dl.dlocat * dlocat.deta,
                        dl.dscale * dscale.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M = Musual)]
    ans
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .doff = doff ))),
  weight = eval(substitute(expression({

    const1 <- (Dof + 1) / (Dof + 3)
    const2 <- (Dof + 0) / (Dof + 3)
    const1[!is.finite( Dof )] <- 1  # Handles Inf
    const2[!is.finite( Dof )] <- 1  # Handles Inf

    const4 <- dnorm(0)
    ned2l.dlocat2 <-        const1 / (Sca * (Kayfun.studentt(Dof) / const4))^2

    ned2l.dscale2 <- 2.0  * const2 /  Sca^2                 # 2.0 seems to work

    wz <- matrix(as.numeric(NA), n, M)  #2=M; diagonal!
    wz[, Musual*(1:NOS) - 1] <- ned2l.dlocat2 * dlocat.deta^2
    wz[, Musual*(1:NOS)    ] <- ned2l.dscale2 * dscale.deta^2

    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .doff = doff  ))))
}





 
 chisq <- function(link = "loge", zero = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")



  new("vglmff",
  blurb = c("Chi-squared distribution\n\n",
            "Link:     ",
            namesof("df", link, earg = earg, tag = FALSE)),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 1
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 1,
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
    Musual <- 1
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    extra$ncoly <- NOS <- ncoly # Number of species
    mynames1 <- paste("df", if (NOS > 1) 1:NOS else "", sep = "")
    predictors.names <-
      namesof(mynames1, .link , earg = .earg , tag = FALSE)

    if (!length(mustart) && !length(etastart))
      mustart <- y + (1 / 8) * (y == 0)
  }), list( .link = link, .earg = earg ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta, .link , earg = .earg )
  }, list( .link = link, .earg = earg ))),

  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <- c(rep( .link , length = ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for(ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .link = link, .earg = earg ))),

  linkfun = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, .link , earg = .earg )
  }, list( .link = link, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    mydf <- eta2theta(eta, .link , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
        sum(c(w) * dchisq(x = y, df = mydf, ncp = 0, log = TRUE))
  }, list( .link = link, .earg = earg ))),
  vfamily = "chisq",
  deriv = eval(substitute(expression({
    mydf <- eta2theta(eta, .link , earg = .earg )
    dl.dv <- (log(y / 2) - digamma(mydf / 2)) / 2
    dv.deta <- dtheta.deta(mydf, .link , earg = .earg )
    c(w) * dl.dv * dv.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    ned2l.dv2 <- trigamma(mydf / 2) / 4
    wz <- ned2l.dv2 * dv.deta^2
    c(w) * wz
  }), list( .link = link, .earg = earg ))))
}







dsimplex <- function(x, mu = 0.5, dispersion = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)
  sigma <- dispersion 

  deeFun <- function(y, mu)
      (((y - mu) / (mu * (1 - mu)))^2) / (y * (1 - y))
  logpdf <- (-0.5 * log(2 * pi) - log(sigma) - 1.5 * log(x) -
            1.5 * log1p(-x) - 0.5 * deeFun(x, mu) / sigma^2)
  logpdf[x     <= 0.0] <- -Inf # log(0.0)
  logpdf[x     >= 1.0] <- -Inf # log(0.0)
  logpdf[mu    <= 0.0] <- NaN
  logpdf[mu    >= 1.0] <- NaN
  logpdf[sigma <= 0.0] <- NaN
  if (log.arg) logpdf else exp(logpdf)
}


rsimplex <- function(n, mu = 0.5, dispersion = 1) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
               stop("bad input for argument 'n'") else n

  oneval <- (length(mu) == 1 && length(dispersion) == 1)
  answer <- rep(0.0, length.out = use.n)
  mu <- rep(mu, length.out = use.n);
  dispersion <- rep(dispersion, length.out = use.n)
  Kay1 <- 3 * (dispersion * mu * (1-mu))^2

  if (oneval) {
    Kay1 <- Kay1[1] # Since oneval means there is only one unique value
    mymu <-   mu[1]
    myroots <- polyroot(c(-mymu^2, Kay1+2*mymu^2, -3*Kay1+1-2*mymu, 2*Kay1))
    myroots <- myroots[abs(Im(myroots)) < 0.00001]
    myroots <- Re(myroots)
    myroots <- myroots[myroots >= 0.0]
    myroots <- myroots[myroots <= 1.0]
    pdfmax <- dsimplex(myroots, mymu, dispersion[1])
    pdfmax <- rep(max(pdfmax), length.out = use.n) # For multiple peaks
  } else {
    pdfmax <- numeric(use.n)
    for (ii in 1:use.n) {
      myroots <- polyroot(c(-mu[ii]^2, Kay1[ii]+2*mu[ii]^2,
                           -3*Kay1[ii]+1-2*mu[ii], 2*Kay1[ii]))
      myroots <- myroots[abs(Im(myroots)) < 0.00001]
      myroots <- Re(myroots)
      myroots <- myroots[myroots >= 0.0]
      myroots <- myroots[myroots <= 1.0]
      pdfmax[ii] <- max(dsimplex(myroots, mu[ii], dispersion[ii]))
    }
  }

  index <- 1:use.n
  nleft <- length(index)
  while (nleft > 0) {
    xx <- runif(nleft) # , 0, 1
    yy <- runif(nleft, max = pdfmax[index])
    newindex <- (1:nleft)[yy < dsimplex(xx, mu[index], dispersion[index])]
    if (length(newindex)) {
      answer[index[newindex]] <- xx[newindex]
      index <- setdiff(index, index[newindex])
      nleft <- nleft - length(newindex)
    }
  }
  answer
}






 simplex <- function(lmu = "logit", lsigma = "loge",
                     imu = NULL, isigma = NULL,
                     imethod = 1, shrinkage.init = 0.95,
                     zero = 2) {







  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lsigma <- as.list(substitute(lsigma))
  esigma <- link2list(lsigma)
  lsigma <- attr(esigma, "function.name")


  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
       imethod > 3)
      stop("argument 'imethod' must be 1 or 2 or 3")
  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
      shrinkage.init < 0 ||
      shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  new("vglmff",
  blurb = c("Univariate simplex distribution\n\n",
            "f(y) = [2*pi*sigma^2*(y*(1-y))^3]^(-0.5) * \n",
            "       exp[-0.5*(y-mu)^2 / (sigma^2 * y * ",
            "(1-y) * mu^2 * (1-mu)^2)],\n",
            "   0 < y < 1, 0 < mu < 1, sigma > 0\n\n",
            "Links:     ",
            namesof("mu", lmu, earg = emu), ", ",
            namesof("sigma", lsigma, earg = esigma), "\n\n",
            "Mean:              mu\n",
            "Variance function: V(mu) = mu^3 * (1 - mu)^3"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
      if (any(y <= 0.0 | y >= 1.0))
        stop("all 'y' values must be in (0,1)")


    w.y.check(w = w, y = y,
              Is.positive.y = TRUE)


      predictors.names <- c(
          namesof("mu",    .lmu ,    earg = .emu ,    tag = FALSE),
          namesof("sigma", .lsigma, earg = .esigma, tag = FALSE))

      deeFun <- function(y, mu)
          (((y - mu) / (mu * (1 - mu)))^2) / (y * (1 - y))

      if (!length(etastart)) {

          use.this =
            if ( .imethod == 3) weighted.mean(y, w = w) else
            if ( .imethod == 1) median(y) else
                                mean(y, trim = 0.1)


          init.mu <- (1 - .sinit) * y + .sinit * use.this
          mu.init <- rep(if (length( .imu )) .imu else init.mu, length = n)
          sigma.init <- if (length( .isigma )) rep( .isigma, leng = n) else {
          use.this <- deeFun(y, mu=init.mu)
          rep(sqrt( if ( .imethod == 3) weighted.mean(use.this, w) else
                    if ( .imethod == 1) median(use.this) else
                                        mean(use.this, trim = 0.1)),
              length = n)
          }
          etastart <-
            cbind(theta2eta(mu.init,    .lmu ,    earg = .emu),
                  theta2eta(sigma.init, .lsigma, earg = .esigma))
      }
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma,
            .imu = imu, .isigma = isigma,
            .sinit = shrinkage.init, .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
      eta2theta(eta[, 1], .lmu , earg = .emu)
  }, list( .lmu = lmu, .emu = emu ))),
  last = eval(substitute(expression({
      misc$link <- c(mu = .lmu ,
                    sigma = .lsigma)
      misc$earg <- list(mu = .emu ,
                       sigma = .esigma)
      misc$imu   <- .imu
      misc$isigma <- .isigma
      misc$imethod <- .imethod
      misc$shrinkage.init <- .sinit
  }), list( .lmu = lmu, .lsigma = lsigma,
            .imu = imu, .isigma = isigma,
            .emu = emu, .esigma = esigma,
            .sinit = shrinkage.init, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      sigma <- eta2theta(eta[, 2], .lsigma, earg = .esigma)
      if (residuals) stop("loglikelihood residuals not ",
                          "implemented yet") else {
        sum(c(w) * dsimplex(x = y, mu = mu, dispersion = sigma, log = TRUE))
      }
  }, list( .lsigma = lsigma, .emu = emu,
           .esigma = esigma ))),
  vfamily = c("simplex"),
  deriv = eval(substitute(expression({
      deeFun <- function(y, mu)
          (((y - mu) / (mu * (1 - mu)))^2) / (y * (1 - y))
      sigma       <- eta2theta(eta[, 2], .lsigma, earg = .esigma)
      dmu.deta    <- dtheta.deta(mu,    .lmu ,    earg = .emu)
      dsigma.deta <- dtheta.deta(sigma, .lsigma, earg = .esigma)

      dl.dmu <- (y - mu) * (deeFun(y, mu) +
               1 / (mu * (1 - mu))^2) / (mu * (1 - mu) * sigma^2)

      dl.dsigma <- (deeFun(y, mu) / sigma^2 - 1) / sigma
      cbind(dl.dmu * dmu.deta,
            dl.dsigma * dsigma.deta)
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma ))),
  weight = eval(substitute(expression({
      wz <- matrix(0.0, n, M)  # Diagonal!!
      eim11 <- 3 / (mu * (1 - mu)) + 1 / (sigma^2 * (mu * (1 - mu))^3)
      wz[, iam(1, 1, M)] <- eim11 * dmu.deta^2
      wz[, iam(2, 2, M)] <- (2 / sigma^2) * dsigma.deta^2
      c(w) * wz
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma ))))
}








 rig <- function(lmu = "identity", llambda = "loge",
                 imu = NULL, ilambda = 1) {


  if (!is.Numeric(ilambda, positive = TRUE))
      stop("bad input for 'ilambda'")


  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")


  new("vglmff",
  blurb = c("Reciprocal inverse Gaussian distribution \n",
            "f(y) = [lambda/(2*pi*y)]^(0.5) * \n",
            "       exp[-0.5*(lambda/y) * (y-mu)^2], ",
            "  0 < y,\n",
            "Links:     ",
            namesof("mu",     lmu, earg = emu), ", ",
            namesof("lambda", llambda, earg = elambda), "\n\n",
            "Mean:     mu"),
  initialize = eval(substitute(expression({


    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <- 
      c(namesof("mu",     .lmu ,    earg = .emu ,     tag = FALSE),
        namesof("lambda", .llambda , earg = .elambda , tag = FALSE))
    if (!length(etastart)) {
      mu.init <- rep(if (length( .imu )) .imu else
                     median(y), length = n)
      lambda.init <- rep(if (length( .ilambda )) .ilambda else
                     sqrt(var(y)), length = n)
      etastart <-
        cbind(theta2eta(mu.init, .lmu , earg = .emu),
              theta2eta(lambda.init, .llambda , earg = .elambda ))
    }
  }), list( .lmu = lmu, .llambda = llambda,
            .emu = emu, .elambda = elambda,
            .imu = imu, .ilambda = ilambda ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .lmu , earg = .emu)
  }, list( .lmu = lmu,
           .emu = emu, .elambda = elambda ))),
  last = eval(substitute(expression({
    misc$d3 <- d3    # because save.weights = FALSE
    misc$link <-    c(mu = .lmu , lambda = .llambda )
    misc$earg <- list(mu = .emu , lambda = .elambda )
    misc$pooled.weight <- pooled.weight
  }), list( .lmu = lmu, .llambda = llambda,
            .emu = emu, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      lambda <- eta2theta(eta[, 2], .llambda , earg = .elambda )
      if (residuals) stop("loglikelihood residuals not ",
                          "implemented yet") else
      sum(c(w) * (-0.5*log(y) + 0.5*log(lambda) - (0.5*lambda/y) * (y-mu)^2))
  }, list( .llambda = llambda,
           .emu = emu, .elambda = elambda ))),
  vfamily = c("rig"),
  deriv = eval(substitute(expression({
      if (iter == 1) {
          d3 <- deriv3(~ w * 
               (-0.5*log(y) + 0.5*log(lambda) - (0.5*lambda/y) * (y-mu)^2),
                      c("mu", "lambda"), hessian= TRUE)
      }

      lambda <- eta2theta(eta[, 2], .llambda , earg = .elambda )

      eval.d3 <- eval(d3)
      dl.dthetas <-  attr(eval.d3, "gradient")

      dmu.deta <- dtheta.deta(mu, .lmu , earg = .emu)
      dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
      dtheta.detas <- cbind(dmu.deta, dlambda.deta)

      dl.dthetas * dtheta.detas
  }), list( .lmu = lmu, .llambda = llambda,
            .emu = emu, .elambda = elambda ))),
  weight = eval(substitute(expression({
      d2l.dthetas2 =  attr(eval.d3, "hessian")

      wz <- matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
      wz[, iam(1, 1, M)] <- -d2l.dthetas2[, 1, 1] * dtheta.detas[, 1]^2
      wz[, iam(2, 2, M)] <- -d2l.dthetas2[, 2, 2] * dtheta.detas[, 2]^2
      wz[, iam(1, 2, M)] <- -d2l.dthetas2[, 1, 2] * dtheta.detas[, 1] *
                                               dtheta.detas[, 2]
      if (!.expected) {
          d2mudeta2 <- d2theta.deta2(mu, .lmu , earg = .emu)
          d2lambda <- d2theta.deta2(lambda, .llambda , earg = .elambda )
          wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] - dl.dthetas[, 1] * d2mudeta2
          wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] - dl.dthetas[, 2] * d2lambda
      }

      if (intercept.only) {
          sumw <- sum(w)
          for(ii in 1:ncol(wz))
              wz[, ii] <- sum(wz[, ii]) / sumw
          pooled.weight <- TRUE
          wz <- c(w) * wz   # Put back the weights
      } else
          pooled.weight <- FALSE

      wz
  }), list( .lmu = lmu, .llambda = llambda, .expected = FALSE,
            .emu = emu, .elambda = elambda ))))
}



 hypersecant <- function(link.theta = elogit(min = -pi/2, max = pi/2),
                         init.theta = NULL) {


  link.theta <- as.list(substitute(link.theta))
  earg <- link2list(link.theta)
  link.theta <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Hyperbolic Secant distribution \n",
            "f(y) = exp(theta*y + log(cos(theta ))) / (2*cosh(pi*y/2))\n",
            "  for all y,\n",
            "Link:     ",
            namesof("theta", link.theta , earg = earg), "\n\n",
            "Mean:     tan(theta)"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)




    predictors.names <-
      namesof("theta", .link.theta , earg = .earg , tag = FALSE)
    if (!length(etastart)) {
      theta.init <- rep(if (length( .init.theta )) .init.theta else
                       median(y), length = n)
      etastart <-
        theta2eta(theta.init, .link.theta , earg = .earg )
    }
  }), list( .link.theta = link.theta , .earg = earg,
            .init.theta = init.theta ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    theta <- eta2theta(eta, .link.theta , earg = .earg )
    tan(theta)
  }, list( .link.theta = link.theta , .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <- c(theta = .link.theta )
    misc$earg <- list(theta = .earg )
    misc$expected <- TRUE
  }), list( .link.theta = link.theta , .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    theta <- eta2theta(eta, .link.theta , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
    sum(c(w) * (theta*y + log(cos(theta)) - log(cosh(pi*y/2 ))))
  }, list( .link.theta = link.theta , .earg = earg ))),
  vfamily = c("hypersecant"),
  deriv = eval(substitute(expression({
    theta <- eta2theta(eta, .link.theta , earg = .earg )
    dl.dthetas <-  y - tan(theta)
    dparam.deta <- dtheta.deta(theta, .link.theta , earg = .earg )
    c(w) * dl.dthetas * dparam.deta
  }), list( .link.theta = link.theta , .earg = earg ))),
  weight = expression({
    d2l.dthetas2 <-  1 / cos(theta)^2
    wz <- c(w) * d2l.dthetas2 * dparam.deta^2
    wz
  }))
}



 hypersecant.1 <- function(link.theta = elogit(min = -pi/2, max = pi/2),
                           init.theta = NULL) {


  link.theta <- as.list(substitute(link.theta))
  earg <- link2list(link.theta)
  link.theta <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Hyperbolic Secant distribution \n",
            "f(y) = (cos(theta)/pi) * y^(-0.5+theta/pi) * \n",
            "       (1-y)^(-0.5-theta/pi), ",
            "  0 < y < 1,\n",
            "Link:     ",
            namesof("theta", link.theta , earg = earg), "\n\n",
            "Mean:     0.5 + theta/pi", "\n",
            "Variance: (pi^2 - 4*theta^2) / (8*pi^2)"),
  initialize = eval(substitute(expression({
    if (any(y <= 0 | y >= 1))
      stop("all response 'y' values must be in (0,1)")


    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1)




    predictors.names <-
      namesof("theta", .link.theta , earg = .earg , tag = FALSE)

    if (!length(etastart)) {
    theta.init <- rep(if (length( .init.theta )) .init.theta else
                     median(y), length = n)

    etastart <-
        theta2eta(theta.init, .link.theta , earg = .earg )
    }
  }), list( .link.theta = link.theta , .earg = earg,
            .init.theta = init.theta ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    theta <- eta2theta(eta, .link.theta , earg = .earg )
    0.5 + theta / pi
  }, list( .link.theta = link.theta , .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <- c(theta = .link.theta )
    misc$earg <- list(theta = .earg )
    misc$expected <- TRUE
  }), list( .link.theta = link.theta , .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    theta <- eta2theta(eta, .link.theta , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
    sum(c(w) * (log(cos(theta)) + (-0.5+theta/pi)*log(y) +
            (-0.5-theta/pi)*log1p(-y )))
  }, list( .link.theta = link.theta , .earg = earg ))),
  vfamily = c("hypersecant.1"),
  deriv = eval(substitute(expression({
    theta <- eta2theta(eta, .link.theta , earg = .earg )
    dl.dthetas <-  -tan(theta) + log(y/(1-y)) / pi 
    dparam.deta <- dtheta.deta(theta, .link.theta , earg = .earg )
    c(w) * dl.dthetas * dparam.deta
  }), list( .link.theta = link.theta , .earg = earg ))),
  weight = expression({
    d2l.dthetas2 <-  1 / cos(theta)^2
    wz <- c(w) * d2l.dthetas2 * dparam.deta^2
    wz
  }))
}



 leipnik <- function(lmu = "logit", llambda = "loge",
                     imu = NULL,    ilambda = NULL) {




  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")


  if (is.Numeric(ilambda) && any(ilambda <= -1))
    stop("argument 'ilambda' must be > -1")



  new("vglmff",
  blurb = c("Leipnik's distribution \n",
            "f(y) = (y(1-y))^(-1/2) * [1 + (y-mu)^2 / (y*(1-y))]^(-lambda/2) /\n",
            "       Beta[(lambda+1)/2, 1/2], ",
            "  0 < y < 1,  lambda > -1\n",
            "Links:     ",
            namesof("mu", lmu, earg = emu), ", ",
            namesof("lambda", llambda, earg = elambda), "\n\n",
            "Mean:     mu\n",
            "Variance: mu*(1-mu)"),
  initialize = eval(substitute(expression({
      if (any(y <= 0 | y >= 1))
        stop("all response 'y' values must be in (0,1)")


    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1)




      predictors.names <-
        c(namesof("mu",     .lmu ,     earg = .emu ,     tag = FALSE),
          namesof("lambda", .llambda , earg = .elambda , tag = FALSE))

    if (!length(etastart)) {
      mu.init <- rep(if (length( .imu )) .imu else
                    (y), length = n)
      lambda.init <- rep(if (length( .ilambda )) .ilambda else
                     1/var(y), length = n)
      etastart <-
       cbind(theta2eta(mu.init,     .lmu ,     earg = .emu ),
             theta2eta(lambda.init, .llambda , earg = .elambda ))
    }
  }), list( .lmu = lmu, .llambda = llambda,
            .emu = emu, .elambda = elambda,
            .imu = imu, .ilambda = ilambda ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .lmu , earg = .emu)
  }, list( .lmu = lmu,
           .emu = emu, .elambda = elambda ))),
  last = eval(substitute(expression({
    misc$link <-    c(mu = .lmu , lambda = .llambda )
    misc$earg <- list(mu = .emu , lambda = .elambda )

    misc$pooled.weight <- pooled.weight
    misc$expected <- FALSE
  }), list( .lmu = lmu, .llambda = llambda,
            .emu = emu, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    lambda <- eta2theta(eta[, 2], .llambda , earg = .elambda )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
    sum(c(w) * (-0.5*log(y*(1-y)) - 0.5 * lambda *
            log1p((y-mu)^2 / (y*(1-y ))) - lgamma((lambda+1)/2) +
            lgamma(1+ lambda/2 )))
  }, list( .llambda = llambda,
           .emu = emu, .elambda = elambda ))),
  vfamily = c("leipnik"),
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta[, 2], .llambda , earg = .elambda )
    dl.dthetas =
      cbind(dl.dmu = lambda*(y-mu) / (y*(1-y)+(y-mu)^2),
            dl.dlambda= -0.5 * log1p((y-mu)^2 / (y*(1-y))) -
            0.5*digamma((lambda+1)/2) +
            0.5*digamma(1+lambda/2))

    dmu.deta <- dtheta.deta(mu, .lmu , earg = .emu)
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    dtheta.detas <- cbind(dmu.deta, dlambda.deta)

    c(w) * dl.dthetas * dtheta.detas
  }), list( .lmu = lmu, .llambda = llambda,
            .emu = emu, .elambda = elambda ))),
  weight = eval(substitute(expression({
    denominator <- y*(1-y) + (y-mu)^2
    d2l.dthetas2 <-  array(NA, c(n, 2, 2))
    d2l.dthetas2[, 1, 1] <- c(w) * lambda*(-y*(1-y)+(y-mu)^2)/denominator^2
    d2l.dthetas2[, 1, 2] <- 
    d2l.dthetas2[, 2, 1] <- c(w) * (y-mu) / denominator
    d2l.dthetas2[, 2, 2] <- c(w) * (-0.25*trigamma((lambda+1)/2) +
                                 0.25*trigamma(1+lambda/2))

    wz <- matrix(as.numeric(NA), n, dimm(M)) #3=dimm(M)
    wz[, iam(1, 1, M)] <- -d2l.dthetas2[, 1, 1] * dtheta.detas[, 1]^2
    wz[, iam(2, 2, M)] <- -d2l.dthetas2[, 2, 2] * dtheta.detas[, 2]^2
    wz[, iam(1, 2, M)] <- -d2l.dthetas2[, 1, 2] * dtheta.detas[, 1] *
                                                 dtheta.detas[, 2]
    if (!.expected) {
      d2mudeta2 <- d2theta.deta2(mu, .lmu , earg = .emu)
      d2lambda <- d2theta.deta2(lambda, .llambda , earg = .elambda )
      wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] - dl.dthetas[, 1] * d2mudeta2
      wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] - dl.dthetas[, 2] * d2lambda
    }

    if (intercept.only) {
        sumw <- sum(w)
        for(ii in 1:ncol(wz))
          wz[, ii] <- sum(wz[, ii]) / sumw
        pooled.weight <- TRUE
        wz <- c(w) * wz # Put back the weights
    } else
        pooled.weight <- FALSE

    wz
  }), list( .lmu = lmu, .llambda = llambda, .expected = FALSE,
            .emu = emu, .elambda = elambda ))))
}





 invbinomial <- function(lrho = elogit(min = 0.5, max = 1),
                         llambda = "loge",
                         irho = NULL,
                         ilambda = NULL,
                         zero = NULL) {



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")



  lrho <- as.list(substitute(lrho))
  erho <- link2list(lrho)
  lrho <- attr(erho, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")


  new("vglmff",
  blurb = c("Inverse binomial distribution\n\n",
            "Links:    ",
            namesof("rho", lrho, earg = erho), ", ", 
            namesof("lambda", llambda, earg = elambda), "\n", 
            "Mean:     lambda*(1-rho)/(2*rho-1)\n",
            "Variance: lambda*rho*(1-rho)/(2*rho-1)^3\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
    c(namesof("rho", .lrho, earg = .erho, tag = FALSE),
      namesof("lambda", .llambda , earg = .elambda , tag = FALSE))

    if (!length(etastart)) {
      covarn <- sd(c(y))^2 / weighted.mean(y, w)
      temp1 <- 0.5 + (1 + sqrt(1+8*covarn)) / (8*covarn)
      temp2 <- 0.5 + (1 - sqrt(1+8*covarn)) / (8*covarn)
      init.rho <- rep(if (length( .irho)) .irho else {
        ifelse(temp1 > 0.5 && temp1 < 1, temp1, temp2)
      }, length = n)
      init.lambda <- rep(if (length( .ilambda)) .ilambda else {
        (2*init.rho-1) * weighted.mean(y, w) / (1-init.rho)
      }, length = n)
      etastart <-
        cbind(theta2eta(init.rho, .lrho, earg = .erho),
              theta2eta(init.lambda, .llambda , earg = .elambda ))
    }
  }), list( .llambda = llambda, .lrho = lrho,
            .elambda = elambda, .erho = erho,
            .ilambda = ilambda, .irho = irho ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    rho <- eta2theta(eta[, 1], .lrho, earg = .erho)
    lambda <- eta2theta(eta[, 2], .llambda , earg = .elambda )
    ifelse(rho > 0.5, lambda*(1-rho)/(2*rho-1), NA)
  }, list( .llambda = llambda, .lrho = lrho,
           .elambda = elambda, .erho = erho ))),
  last = eval(substitute(expression({
    misc$link <- c(rho= .lrho, lambda = .llambda )
    misc$earg <- list(rho= .erho, lambda = .elambda )
    misc$pooled.weight <- pooled.weight
  }), list( .llambda = llambda, .lrho = lrho,
            .elambda = elambda, .erho = erho ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    rho <- eta2theta(eta[, 1], .lrho, earg = .erho)
    lambda <- eta2theta(eta[, 2], .llambda , earg = .elambda )

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
    sum(c(w) * (log(lambda) - lgamma(2*y+lambda) - lgamma(y+1) -
           lgamma(y+lambda+1) + y*log(rho) + y*log1p(-rho) +
           lambda*log(rho)))
  }, list( .llambda = llambda, .lrho = lrho,
           .elambda = elambda, .erho = erho ))),
  vfamily = c("invbinomial"),
  deriv = eval(substitute(expression({
    rho <- eta2theta(eta[, 1], .lrho, earg = .erho)
    lambda <- eta2theta(eta[, 2], .llambda , earg = .elambda )

    dl.drho <- (y + lambda)/rho - y/(1-rho)
    dl.dlambda <- 1/lambda - digamma(2*y+lambda) - digamma(y+lambda+1) +
                 log(rho)

    drho.deta <- dtheta.deta(rho, .lrho, earg = .erho)
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )

    c(w) * cbind(dl.drho * drho.deta,
                 dl.dlambda * dlambda.deta )
  }), list( .llambda = llambda, .lrho = lrho,
              .elambda = elambda, .erho = erho ))),
  weight = eval(substitute(expression({
    ned2l.drho2 <- (mu+lambda) / rho^2 + mu / (1-rho)^2
    d2l.dlambda2 <- 1/(lambda^2) + trigamma(2*y+lambda)+trigamma(y+lambda+1)
    ned2l.dlambdarho <- -1/rho

    wz <- matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.drho2 * drho.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dlambdarho * dlambda.deta * drho.deta
    wz[, iam(2, 2, M)] <-  d2l.dlambda2 * dlambda.deta^2

    d2rhodeta2 <- d2theta.deta2(rho, .lrho, earg = .erho)
    d2lambda.deta2 <- d2theta.deta2(lambda, .llambda , earg = .elambda )
    wz <- c(w) * wz

    if (intercept.only) {
      pooled.weight <- TRUE

      wz[, iam(2, 2, M)] <-  sum(wz[, iam(2, 2, M)]) / sum(w)

    } else {
      pooled.weight <- FALSE
    }

    wz
  }), list( .llambda = llambda, .lrho = lrho,
            .elambda = elambda, .erho = erho ))))
}



 genpoisson <- function(llambda = elogit(min = -1, max = 1),
                        ltheta = "loge",
                        ilambda = NULL, itheta = NULL,
                        use.approx = TRUE,
                        imethod = 1, zero = 1) {


  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  ltheta <- as.list(substitute(ltheta))
  etheta <- link2list(ltheta)
  ltheta <- attr(etheta, "function.name")


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.logical(use.approx) || length(use.approx) != 1)
    stop("'use.approx' must be logical value")




  new("vglmff",
  blurb = c("Generalized Poisson distribution\n\n",
            "Links:    ",
            namesof("lambda", llambda, earg = elambda), ", ", 
            namesof("theta", ltheta, earg = etheta), "\n", 
            "Mean:     theta / (1-lambda)\n",
            "Variance: theta / (1-lambda)^3"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)





    predictors.names <-
       c(namesof("lambda", .llambda , earg = .elambda , tag = FALSE),
         namesof("theta",  .ltheta ,  earg = .etheta,  tag = FALSE))
    init.lambda <- if ( .imethod == 1)
      1 - sqrt(weighted.mean(y, w) / var(y)) else 0.5
    init.theta  <- if ( .imethod == 1)
      sqrt((0.01 + weighted.mean(y, w)^3) / var(y)) else
      median(y) * (1-init.lambda)
    if (init.theta <= 0)
      init.theta <- 0.1
    cutpt <- if (init.lambda < 0) {
      mmm <- max(trunc(-init.theta / init.lambda), 4)
      max(-1, -init.theta /mmm)
    } else -1
    if (init.lambda <= cutpt)
      init.lambda <- cutpt + 0.1
    if (init.lambda >= 1)
      init.lambda <- 0.9
    if (!length(etastart)) {
      lambda <- rep(if (length( .ilambda)) .ilambda else
                   init.lambda, length = n)
      theta <- rep(if (length( .itheta)) .itheta else init.theta ,
                  length = n)
      etastart <-
        cbind(theta2eta(lambda, .llambda , earg = .elambda ),
              theta2eta(theta,  .ltheta ,  earg = .etheta ))
    }
  }), list( .ltheta = ltheta, .llambda = llambda,
            .etheta = etheta, .elambda = elambda,
            .imethod = imethod,
            .itheta = itheta, .ilambda = ilambda )) ),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    lambda <- eta2theta(eta[, 1], .llambda , earg = .elambda )
    theta  <- eta2theta(eta[, 2], .ltheta , earg = .etheta )
    theta / (1 - lambda)
  }, list( .ltheta = ltheta, .llambda = llambda,
           .etheta = etheta, .elambda = elambda ))),
  last = eval(substitute(expression({
    misc$link <-    c(lambda = .llambda , theta = .ltheta )
    misc$earg <- list(lambda = .elambda , theta = .etheta )
    if (! .use.approx )
      misc$pooled.weight <- pooled.weight
  }), list( .ltheta = ltheta, .llambda = llambda,
            .use.approx = use.approx,
            .etheta = etheta, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    lambda <- eta2theta(eta[, 1], .llambda , earg = .elambda )
    theta  <- eta2theta(eta[, 2], .ltheta  , earg = .etheta  )
    index <- (y == 0)
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
    sum(w[index] * (-theta[index])) +
    sum(w[!index] * (-y[!index]*lambda[!index]-theta[!index] +
        (y[!index]-1)*log(theta[!index]+y[!index]*lambda[!index]) +
        log(theta[!index]) - lgamma(y[!index]+1)) )
  }, list( .ltheta = ltheta, .llambda = llambda,
           .etheta = etheta, .elambda = elambda ))),
  vfamily = c("genpoisson"),
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta[, 1], .llambda , earg = .elambda )
    theta <- eta2theta(eta[, 2], .ltheta , earg = .etheta )
    dl.dlambda <- -y + y*(y-1)/(theta+y*lambda)
    dl.dtheta <- -1 + (y-1)/(theta+y*lambda) + 1/theta
    dTHETA.deta <- dtheta.deta(theta, .ltheta , earg = .etheta )
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    c(w) * cbind(dl.dlambda * dlambda.deta,
                 dl.dtheta * dTHETA.deta )
  }), list( .ltheta = ltheta, .llambda = llambda,
            .etheta = etheta, .elambda = elambda ))),
  weight = eval(substitute(expression({
    wz <- matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
    if ( .use.approx ) {
        BBB <- (theta+2)*(theta+2*lambda-theta*lambda)-(theta^2)*(1-lambda)
        d2l.dlambda2 <- 2 * theta * (theta+2) / ((1-lambda) * BBB)
        d2l.dtheta2 <- 2 * (1 + lambda * (2/theta - 1)) / BBB
        d2l.dthetalambda <-  2 * theta / BBB
        wz[, iam(1, 1, M)] <- d2l.dlambda2 * dlambda.deta^2
        wz[, iam(2, 2, M)] <- d2l.dtheta2 * dTHETA.deta^2
        wz[, iam(1, 2, M)] <- d2l.dthetalambda * dTHETA.deta * dlambda.deta
        wz <- c(w) * wz
    } else {
        d2l.dlambda2 <- -y^2 * (y-1) / (theta+y*lambda)^2
        d2l.dtheta2 <- -(y-1)/(theta+y*lambda)^2 - 1 / theta^2
        d2l.dthetalambda <-  -y * (y-1) / (theta+y*lambda)^2 
        wz[, iam(1, 1, M)] <- -d2l.dlambda2 * dlambda.deta^2
        wz[, iam(2, 2, M)] <- -d2l.dtheta2 * dTHETA.deta^2
        wz[, iam(1, 2, M)] <- -d2l.dthetalambda * dTHETA.deta * dlambda.deta

        d2THETA.deta2 <- d2theta.deta2(theta, .ltheta , earg = .etheta )
        d2lambdadeta2 <- d2theta.deta2(lambda, .llambda , earg = .elambda )
        wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] - dl.dlambda * d2lambdadeta2
        wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] - dl.dtheta * d2THETA.deta2
        wz <- c(w) * wz

        if (intercept.only) {
          sumw <- sum(w)
          for(ii in 1:ncol(wz))
            wz[, ii] <- sum(wz[, ii]) / sumw
          pooled.weight <- TRUE
          wz <- c(w) * wz   # Put back the weights
        } else
          pooled.weight <- FALSE
        }
    wz
  }), list( .ltheta = ltheta, .llambda = llambda,
            .use.approx = use.approx,
            .etheta = etheta, .elambda = elambda ))))
}






dlgamma <- function(x, location = 0, scale = 1, k = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(scale, positive = TRUE))
    stop("bad input for argument 'scale'")
  if (!is.Numeric(k, positive = TRUE))
    stop("bad input for argument 'k'")
  z <- (x-location) / scale
  if (log.arg) {
    k * z - exp(z) - log(scale) - lgamma(k)
  } else {
    exp(k * z - exp(z)) / (scale * gamma(k))
  }
}


plgamma <- function(q, location = 0, scale = 1, k = 1) {

  zedd <- (q - location) / scale
  ans <- pgamma(exp(zedd), k)
  ans[scale <  0] <- NaN
  ans
}


qlgamma <- function(p, location = 0, scale = 1, k = 1) {
  if (!is.Numeric(scale, positive = TRUE))
    stop("bad input for argument 'scale'")

  ans <- location + scale * log(qgamma(p, k))
  ans[scale <  0] <- NaN
  ans
}


rlgamma <- function(n, location = 0, scale = 1, k = 1) {
  ans <- location + scale * log(rgamma(n, k))
  ans[scale < 0] <- NaN
  ans
}



 lgammaff <- function(link = "loge", init.k = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Log-gamma distribution ",
            "f(y) = exp(ky - e^y)/gamma(k)), k>0\n\n",
            "Link:    ",
            namesof("k", link, earg = earg), "\n", "\n",
            "Mean:    digamma(k)", "\n"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)


    predictors.names <-
      namesof("k", .link , earg = .earg , tag = FALSE) 

    if (!length(etastart)) {
      k.init <- if (length( .init.k))
               rep( .init.k, length.out = length(y)) else {
               medy = median(y)
          if (medy < 2) 5 else if (medy < 4) 20 else exp(0.7 * medy)
        }
      etastart <- theta2eta(k.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .init.k = init.k ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    k <- eta2theta(eta, .link , earg = .earg )
    digamma(k)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <- c(k = .link )
    misc$earg <- list(k = .earg )
    misc$expected <- TRUE
  }), list( .link = link, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      kk <- eta2theta(eta, .link , earg = .earg )
      if (residuals) stop("loglikelihood residuals not ",
                          "implemented yet") else {
          sum(c(w) * dlgamma(x = y, location = 0, scale = 1,
                             k = kk, log = TRUE))
      }
    }, list( .link = link, .earg = earg ))),
  vfamily = c("lgammaff"),
  deriv = eval(substitute(expression({
    kk <- eta2theta(eta, .link , earg = .earg ) 
    dl.dk <- y - digamma(kk)
    dk.deta <- dtheta.deta(kk, .link , earg = .earg )
    c(w) * dl.dk * dk.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    ned2l.dk2 <- trigamma(kk)
    wz <- c(w) * dk.deta^2 * ned2l.dk2
    wz
  }), list( .link = link, .earg = earg ))))
}







 lgamma3ff <-
  function(llocation = "identity", lscale = "loge", lshape = "loge",
           ilocation = NULL, iscale = NULL, ishape = 1, zero = NULL) {

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")
  if (length(iscale) &&
      !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")




  new("vglmff",
  blurb = c("Log-gamma distribution",
            " f(y) = exp(k(y-a)/b - e^((y-a)/b))/(b*gamma(k)), ",
            "location=a, scale=b>0, shape=k>0\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n\n",
            "Mean:     a + b*digamma(k)", "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
      c(namesof("location", .llocat , earg = .elocat , tag = FALSE),
        namesof("scale", .lscale , earg = .escale , tag = FALSE),
        namesof("shape", .lshape , earg = .eshape, tag = FALSE))


    if (!length(etastart)) {
      k.init <- if (length( .ishape ))
               rep( .ishape, length.out = length(y)) else {
          rep(exp(median(y)), length.out = length(y))
      }
      scale.init <- if (length( .iscale ))
          rep( .iscale, length.out = length(y)) else {
          rep(sqrt(var(y) / trigamma(k.init)), length.out = length(y))
      }
      loc.init <- if (length( .ilocat ))
          rep( .ilocat, length.out = length(y)) else {
          rep(median(y) - scale.init * digamma(k.init),
              length.out = length(y))
      }
      etastart <-
        cbind(theta2eta(loc.init, .llocat , earg = .elocat ),
              theta2eta(scale.init, .lscale , earg = .escale ),
              theta2eta(k.init, .lshape , earg = .eshape ))
    }
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape,
            .ilocat = ilocat, .iscale = iscale, .ishape = ishape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat , earg = .elocat ) +
    eta2theta(eta[, 2], .lscale , earg = .escale ) *
    digamma(eta2theta(eta[, 3], .lshape , earg = .eshape ))
  }, list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
           .elocat = elocat, .escale = escale, .eshape = eshape))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat , scale = .lscale ,
                     shape = .lshape)

    misc$earg <- list(location = .elocat , scale = .escale ,
                     shape = .eshape )
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    aa <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    bb <- eta2theta(eta[, 2], .lscale , earg = .escale )
    kk <- eta2theta(eta[, 3], .lshape , earg = .eshape )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
            sum(c(w) * dlgamma(x = y, locat=aa, scale=bb, k=kk, log = TRUE))
      }
  }, list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
           .elocat = elocat, .escale = escale, .eshape = eshape))),
  vfamily = c("lgamma3ff"),
  deriv = eval(substitute(expression({
    a <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    b <- eta2theta(eta[, 2], .lscale , earg = .escale )
    k <- eta2theta(eta[, 3], .lshape , earg = .eshape )

    zedd <- (y-a)/b
    dl.da <- (exp(zedd) - k) / b
    dl.db <- (zedd * (exp(zedd) - k) - 1) / b
    dl.dk <- zedd - digamma(k)

    da.deta <- dtheta.deta(a, .llocat , earg = .elocat )
    db.deta <- dtheta.deta(b, .lscale , earg = .escale )
    dk.deta <- dtheta.deta(k, .lshape , earg = .eshape )

    c(w) * cbind(dl.da * da.deta,
                 dl.db * db.deta,
                 dl.dk * dk.deta)
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape))),
  weight = eval(substitute(expression({
    ned2l.da2 <- k / b^2
    ned2l.db2 <- (1 + k*(trigamma(k+1) + (digamma(k+1))^2)) / b^2
    ned2l.dk2 <- trigamma(k)
    ned2l.dadb <- (1 + k*digamma(k)) / b^2
    ned2l.dadk <- 1 / b
    ned2l.dbdk <- digamma(k) / b

    wz <- matrix(as.numeric(NA), n, dimm(M))
    wz[, iam(1, 1, M)] <- ned2l.da2 * da.deta^2
    wz[, iam(2, 2, M)] <- ned2l.db2 * db.deta^2
    wz[, iam(3, 3, M)] <- ned2l.dk2 * dk.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dadb * da.deta * db.deta
    wz[, iam(1, 3, M)] <- ned2l.dadk * da.deta * dk.deta
    wz[, iam(2, 3, M)] <- ned2l.dbdk * db.deta * dk.deta
    wz <- c(w) * wz
    wz
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape))))
}



 prentice74 <- function(
        llocation = "identity", lscale = "loge", lshape = "identity",
        ilocation = NULL, iscale = NULL, ishape = NULL, zero = 2:3) {

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")
  if (length(iscale) &&
     !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")



  new("vglmff",
  blurb = c("Log-gamma distribution (Prentice, 1974)",
            " f(y) = |q| * exp(w/q^2 - e^w) / (b*gamma(1/q^2)) ,\n",
            "w = (y-a)*q/b + digamma(1/q^2), ",
            "location = a, scale = b > 0, shape = q\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n", "\n",
            "Mean:     a", "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({


    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)


    predictors.names <-
    c(namesof("location", .llocat , earg = .elocat , tag = FALSE),
      namesof("scale", .lscale , earg = .escale , tag = FALSE),
      namesof("shape", .lshape , earg = .eshape, tag = FALSE))



    if (!length(etastart)) {
        sdy <- sqrt(var(y))
        k.init <- if (length( .ishape ))
            rep( .ishape, length.out = length(y)) else {
            skewness <- mean((y-mean(y))^3) / sdy^3 # <0 Left Skewed
            rep(-skewness, length.out = length(y))
        }
        scale.init <- if (length( .iscale ))
            rep( .iscale, length.out = length(y)) else {
            rep(sdy, length.out = length(y))
        }
        loc.init <- if (length( .iloc ))
                   rep( .iloc, length.out = length(y)) else {
              rep(median(y), length.out = length(y))
          }
          etastart <-
            cbind(theta2eta(loc.init, .llocat , earg = .elocat ),
                  theta2eta(scale.init, .lscale , earg = .escale ),
                  theta2eta(k.init, .lshape , earg = .eshape ))
      }
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape,
            .iloc = ilocat, .iscale = iscale, .ishape = ishape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat , earg = .elocat )
  }, list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
           .elocat = elocat, .escale = escale, .eshape = eshape))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat , scale = .lscale ,
                     shape = .lshape )
    misc$earg <- list(location = .elocat , scale = .escale ,
                     shape = .eshape )
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    a <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    b <- eta2theta(eta[, 2], .lscale , earg = .escale )
    k <- eta2theta(eta[, 3], .lshape , earg = .eshape )
    tmp55 <- k^(-2)
    doubw <- (y-a)*k/b + digamma(tmp55)
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
    sum(c(w)*(log(abs(k)) - log(b) - lgamma(tmp55) +
              doubw * tmp55 - exp(doubw )))
  }, list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
           .elocat = elocat, .escale = escale, .eshape = eshape))),
  vfamily = c("prentice74"),
  deriv = eval(substitute(expression({
    a <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    b <- eta2theta(eta[, 2], .lscale , earg = .escale )
    k <- eta2theta(eta[, 3], .lshape , earg = .eshape )

    tmp55 <- k^(-2)
    mustar <- digamma(tmp55)
    doubw <- (y-a)*k/b + mustar
    sigmastar2 <- trigamma(tmp55)

    dl.da <- k*(exp(doubw) - tmp55) / b
    dl.db <- ((doubw - mustar) * (exp(doubw) - tmp55) - 1) / b
    dl.dk <- 1/k - 2 * (doubw - mustar) / k^3 - (exp(doubw) - tmp55) *
            ((doubw - mustar) / k - 2 * sigmastar2 / k^3)

    da.deta <- dtheta.deta(a, .llocat , earg = .elocat )
    db.deta <- dtheta.deta(b, .lscale , earg = .escale )
    dk.deta <- dtheta.deta(k, .lshape , earg = .eshape )

    c(w) * cbind(dl.da * da.deta,
                 dl.db * db.deta,
                 dl.dk * dk.deta)
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape))),
  weight = eval(substitute(expression({
    ned2l.da2 <- 1 / b^2
    ned2l.db2 <- (1 + sigmastar2*tmp55) / b^2
    ned2l.dk2 <- tmp55 - 3*sigmastar2*tmp55^2 + 4*sigmastar2*tmp55^4 *
               (sigmastar2 - k^2)
    ned2l.dadb <- k / b^2
    ned2l.dadk <- (2*(sigmastar2*tmp55^2 - tmp55) - 1) / b
    ned2l.dbdk <- (sigmastar2*tmp55 - 1) / (b*k)

    wz <- matrix(as.numeric(NA), n, dimm(M))
    wz[, iam(1, 1, M)] <- ned2l.da2 * da.deta^2
    wz[, iam(2, 2, M)] <- ned2l.db2 * db.deta^2
    wz[, iam(3, 3, M)] <- ned2l.dk2 * dk.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dadb * da.deta * db.deta
    wz[, iam(1, 3, M)] <- ned2l.dadk * da.deta * dk.deta
    wz[, iam(2, 3, M)] <- ned2l.dbdk * db.deta * dk.deta
    wz <- c(w) * wz
    wz
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape))))
}



dgengamma <- function(x, scale = 1, d = 1, k = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(scale, positive = TRUE))
    stop("bad input for argument 'scale'")
  if (!is.Numeric(d, positive = TRUE))
    stop("bad input for argument 'd'")
  if (!is.Numeric(k, positive = TRUE))
    stop("bad input for argument 'k'")

  N <- max(length(x), length(scale), length(d), length(k))
  x <- rep(x, length.out = N);
  scale <- rep(scale, length.out = N);
  d <- rep(d, length.out = N);
  k <- rep(k, length.out = N); 

  Loglik <- rep(log(0), length.out = N)
  xok <- x > 0
  if (any(xok)) {
    zedd <- (x[xok]/scale[xok])^d[xok]
    Loglik[xok] <- log(d[xok]) + (-d[xok]*k[xok]) * log(scale[xok]) +
               (d[xok]*k[xok]-1) * log(x[xok]) - zedd - lgamma(k[xok])
  }
  if (log.arg) {
    Loglik
  } else {
    exp(Loglik)
  }
}




pgengamma <- function(q, scale = 1, d = 1, k = 1) {
  zedd <- (q / scale)^d
  ans <- pgamma(zedd, k)
  ans[scale <  0] <- NaN
  ans[d     <= 0] <- NaN
  ans
}


qgengamma <- function(p, scale = 1, d = 1, k = 1) {
  ans <- scale * qgamma(p, k)^(1/d)
  ans[scale <  0] <- NaN
  ans[d     <= 0] <- NaN
  ans
}


rgengamma <- function(n, scale = 1, d = 1, k = 1) {

  ans <- scale * rgamma(n, k)^(1/d)
  ans[scale <  0] <- NaN
  ans[d     <= 0] <- NaN
  ans
}


 gengamma <- function(lscale = "loge", ld = "loge", lk = "loge",
                      iscale = NULL, id = NULL, ik = NULL, zero = NULL) {

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  ld <- as.list(substitute(ld))
  ed <- link2list(ld)
  ld <- attr(ed, "function.name")

  lk <- as.list(substitute(lk))
  ek <- link2list(lk)
  lk <- attr(ek, "function.name")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")
  if (length(iscale) &&
      !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")




  new("vglmff",
  blurb = c("Generalized gamma distribution",
            " f(y) = d * b^(-d*k) * y^(d*k-1) * exp(-(y/b)^d) /  gamma(k),\n",
            "scale=b>0, d>0, k>0, y>0\n\n",
            "Links:    ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("d", ld, earg = ed), ", ",
            namesof("k", lk, earg = ek), "\n", "\n",
            "Mean:     b * gamma(k+1/d) / gamma(k)", "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
      c(namesof("scale", .lscale , earg = .escale , tag = FALSE),
        namesof("d", .ld , earg = .ed , tag = FALSE),
        namesof("k", .lk , earg = .ek , tag = FALSE))


    if (!length(etastart)) {
      b.init <- if (length( .iscale ))
          rep( .iscale, length.out = length(y)) else {
          rep(mean(y^2) / mean(y), length.out = length(y))
      }
      k.init <- if (length( .ik ))
          rep( .ik , length.out = length(y)) else {
          rep(mean(y) / b.init, length.out = length(y))
      }
      d.init <- if (length( .id ))
          rep( .id , length.out = length(y)) else {
          rep(digamma(k.init) / mean(log(y / b.init)),
              length.out = length(y))
      }
        etastart <-
          cbind(theta2eta(b.init, .lscale , earg = .escale ),
                theta2eta(d.init, .ld , earg = .ed ),
                theta2eta(k.init, .lk , earg = .ek ))
    }
  }), list( .lscale = lscale, .ld = ld, .lk = lk,
            .escale = escale, .ed = ed, .ek = ek,
            .iscale = iscale, .id = id, .ik = ik ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    b <- eta2theta(eta[, 1], .lscale , earg = .escale )
    d <- eta2theta(eta[, 2], .ld , earg = .ed )
    k <- eta2theta(eta[, 3], .lk , earg = .ek )
    b * gamma(k + 1 / d) / gamma(k)
  }, list( .lscale = lscale, .lk = lk, .ld = ld,
           .escale = escale, .ek = ek, .ed = ed ))),
  last = eval(substitute(expression({
    misc$link <-    c(scale = .lscale , d = .ld , k = .lk )
    misc$earg <- list(scale = .escale , d = .ed , k = .ek )
    misc$expected <- TRUE
  }), list( .lscale = lscale, .ld = ld, .lk = lk,
            .escale = escale, .ed = ed, .ek = ek ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    b <- eta2theta(eta[, 1], .lscale , earg = .escale )
    d <- eta2theta(eta[, 2], .ld , earg = .ed )
    k <- eta2theta(eta[, 3], .lk , earg = .ek )

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dgengamma(x = y, scale = b, d = d, k = k, log = TRUE))
    }
  }, list( .lscale = lscale, .ld = ld, .lk = lk,
           .escale = escale, .ed = ed, .ek = ek ))),
  vfamily = c("gengamma"),
  deriv = eval(substitute(expression({
    b <- eta2theta(eta[, 1], .lscale , earg = .escale )
    d <- eta2theta(eta[, 2], .ld , earg = .ed )
    k <- eta2theta(eta[, 3], .lk , earg = .ek )

    tmp22 <- (y/b)^d
    tmp33 <- log(y/b)
    dl.db <- d * (tmp22 - k) / b
    dl.dd <- 1/d + tmp33 * (k - tmp22)
    dl.dk <- d * tmp33 - digamma(k)

    db.deta <- dtheta.deta(b, .lscale , earg = .escale )
    dd.deta <- dtheta.deta(d, .ld , earg = .ed )
    dk.deta <- dtheta.deta(k, .lk , earg = .ek )

    c(w) * cbind(dl.db * db.deta,
                 dl.dd * dd.deta,
                 dl.dk * dk.deta)
  }), list( .lscale = lscale, .ld = ld, .lk = lk,
            .escale = escale, .ed = ed, .ek = ek ))),
  weight = eval(substitute(expression({
    ned2l.db2 <- k * (d/b)^2
    ned2l.dd2 <- (1 + k * (trigamma(k+1) + (digamma(k+1))^2)) / d^2 
    ned2l.dk2 <- trigamma(k)
    ned2l.dbdd <- -(1 + k*digamma(k)) / b
    ned2l.dbdk <- d / b
    ned2l.dddk <- -digamma(k) / d

    wz <- matrix(as.numeric(NA), n, dimm(M))
    wz[, iam(1, 1, M)] <- ned2l.db2 * db.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dd2 * dd.deta^2
    wz[, iam(3, 3, M)] <- ned2l.dk2 * dk.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dbdd * db.deta * dd.deta
    wz[, iam(1, 3, M)] <- ned2l.dbdk * db.deta * dk.deta
    wz[, iam(2, 3, M)] <- ned2l.dddk * dd.deta * dk.deta

    wz <- c(w) * wz
    wz
  }), list( .lscale = lscale, .ld = ld, .lk = lk,
            .escale = escale, .ed = ed, .ek = ek ))))
}




dlog <- function(x, prob, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

    if (!is.Numeric(prob, positive = TRUE) || max(prob) >= 1)
        stop("bad input for argument 'prob'")
    N <- max(length(x), length(prob))
    if (length(x) != N)
        x <- rep(x, length.out = N)
    if (length(prob) != N)
        prob <- rep(prob, length.out = N)
    ox <- !is.finite(x)
    zero <- ox | round(x) != x | x < 1
    ans <- rep(0.0, length.out = length(x))
        if (log.arg) {
            ans[ zero] <- log(0.0)
            ans[!zero] <- x[!zero] * log(prob[!zero]) - log(x[!zero]) -
                         log(-log1p(-prob[!zero]))
        } else {
            ans[!zero] <- -(prob[!zero]^(x[!zero])) / (x[!zero] *
                         log1p(-prob[!zero]))
        }
    if (any(ox))
        ans[ox] <- NA
    ans
}



plog  <- function(q, prob, log.p = FALSE) {
    if (!is.Numeric(q)) stop("bad input for argument 'q'")
    if (!is.Numeric(prob, positive = TRUE) || max(prob) >= 1)
        stop("bad input for argument 'prob'")
    N <- max(length(q), length(prob))
    q <- rep(q, length.out = N);
    prob <- rep(prob, length.out = N);

    bigno <- 10
    owen1965 <- (q * (1 - prob) > bigno)
    if (specialCase <- any(owen1965)) {
        qqq <- q[owen1965]
        ppp <- prob[owen1965]
        pqp <- qqq * (1 - ppp)
        bigans <- (ppp^(1+qqq) / (1-ppp)) * (1/qqq -
                 1 / (            pqp * (qqq-1)) +
                 2 / ((1-ppp)   * pqp * (qqq-1) * (qqq-2)) -
                 6 / ((1-ppp)^2 * pqp * (qqq-1) * (qqq-2) * (qqq-3)) +
                24 / ((1-ppp)^3 * pqp * (qqq-1) * (qqq-2) * (qqq-3) * (qqq-4)))
        bigans <- 1 + bigans / log1p(-ppp)
    }

    floorq <- pmax(1, floor(q)) # Ensures at least one element per q value
    floorq[owen1965] <- 1
    seqq <- sequence(floorq)
    seqp <- rep(prob, floorq)
    onevector <- (seqp^seqq / seqq) / (-log1p(-seqp))
    rlist <-  dotC(name = "tyee_C_cum8sum",
                  as.double(onevector), answer = double(N),
                  as.integer(N), as.double(seqq),
                  as.integer(length(onevector)), notok=integer(1))
    if (rlist$notok != 0) stop("error in 'cum8sum'")
    ans <- if (log.p) log(rlist$answer) else rlist$answer
    if (specialCase)
        ans[owen1965] <- if (log.p) log(bigans) else bigans
    ans[q < 1] <- if (log.p) log(0.0) else 0.0
    ans
}







rlog <- function(n, prob, Smallno = 1.0e-6) {

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
               stop("bad input for argument 'n'") else n

  if (!is.Numeric(prob, allowable.length = 1, positive = TRUE) ||
      max(prob) >= 1)
    stop("bad input for argument 'prob'")
  if (!is.Numeric(Smallno, positive = TRUE, allowable.length = 1) ||
      Smallno > 0.01 ||
     Smallno < 2 * .Machine$double.eps)
    stop("bad input for argument 'Smallno'")
  ans <- rep(0.0, length.out = use.n)

  ptr1 <- 1; ptr2 <- 0
  a <- -1 / log1p(-prob)
  mean <- a*prob/(1-prob)    # E(Y)
  sigma <- sqrt(a * prob * (1 - a * prob)) / (1 - prob) # sd(Y)
  ymax <- dlog(x = 1, prob)
  while(ptr2 < use.n) {
    Lower <- 0.5 # A continuity correction is used = 1 - 0.5.
    Upper <- mean + 5 * sigma
    while(plog(q = Upper, prob) < 1 - Smallno)
      Upper <- Upper + sigma
    Upper <- Upper + 0.5
    x <- round(runif(2 * use.n, min = Lower, max = Upper))
    index <- runif(2 * use.n, max = ymax) < dlog(x,prob)
    sindex <- sum(index)
    if (sindex) {
      ptr2 <- min(use.n, ptr1 + sindex - 1)
      ans[ptr1:ptr2] <- (x[index])[1:(1+ptr2-ptr1)]
      ptr1 <- ptr2 + 1
    }
  }
  ans
}








 logff <- function(link = "logit", init.c = NULL, zero = NULL) {
  if (length(init.c) &&
     (!is.Numeric(init.c, positive = TRUE) || max(init.c) >= 1))
      stop("init.c must be in (0,1)")

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")



  new("vglmff",
  blurb = c("Logarithmic distribution f(y) = a * c^y / y, ",
             "y = 1, 2, 3,...,\n",
             "            0 < c < 1, a = -1 / log(1-c)  \n\n",
             "Link:    ", namesof("c", link, earg = earg), "\n", "\n",
             "Mean:    a * c / (1 - c)", "\n"),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 1
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 1,
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
    Musual <- 1
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    mynames1  <- paste("c", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
      namesof(mynames1, .link , earg = .earg , tag = FALSE)


    if (!length(etastart)) {
      logff.Loglikfun <- function(probval, y, x, w, extraargs) {
        sum(c(w) * dlog(x = y, prob = probval, log = TRUE))
      }
      Init.c <- matrix(if (length( .init.c )) .init.c else 0,
                       n, M, byrow = TRUE)

      if (!length( .init.c ))
        for(ilocal in 1:ncoly) {
          prob.grid <- seq(0.05, 0.95, by = 0.05)
          Init.c[, ilocal] <- getMaxMin(prob.grid,
                                        objfun = logff.Loglikfun,
                                        y = y[, ilocal], x = x,
                                        w = w[, ilocal])

        }
      etastart <- theta2eta(Init.c, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .init.c = init.c ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    cc <- eta2theta(eta, .link , earg = .earg )
    aa <- -1 / log1p(-cc)
    aa * cc / (1 - cc)
  }, list( .link = link, .earg = earg ))),

  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <- c(rep( .link , length = ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for(ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .link = link, .earg = earg ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    cc <- eta2theta(eta, .link , earg = .earg )
    aa <- -1 / log1p(-cc)
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {
        sum(c(w) * dlog(x = y, prob = -expm1(-1/aa), log = TRUE))
      }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("logff"),
  deriv = eval(substitute(expression({
    Musual <- 1
    cc <- eta2theta(eta, .link , earg = .earg )
    aa <- -1 / log1p(-cc)
    dl.dc <- 1 / ((1 - cc) * log1p(-cc)) + y / cc
    dc.deta <- dtheta.deta(cc, .link , earg = .earg )
    c(w) * dl.dc * dc.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    ned2l.dc2 <- aa * (1 - aa * cc) / (cc * (1-cc)^2)
    wz <- c(w) * dc.deta^2 * ned2l.dc2
    wz
  }), list( .link = link, .earg = earg ))))
}




 levy <- function(delta = NULL, link.gamma = "loge",
                  idelta = NULL, igamma = NULL) {



  delta.known = is.Numeric(delta, allowable.length = 1)

  link.gamma <- as.list(substitute(link.gamma))
  earg <- link2list(link.gamma)
  link.gamma <- attr(earg, "function.name")



  new("vglmff",
  blurb = c("Levy distribution f(y) = sqrt(gamma/(2*pi)) * ",
            "(y-delta)^(-3/2) * \n",
            "          exp(-gamma / (2*(y-delta ))),\n",
            "          delta < y, gamma > 0",
            if (delta.known) paste(", delta = ", delta, ",", sep = ""),
            "\n\n",
            if (delta.known) "Link:    " else "Links:   ",
            namesof("gamma", link.gamma, earg = earg),
            if (! delta.known) 
                c(", ", namesof("delta", "identity", earg = list())),
            "\n\n",
            "Mean:    NA", 
            "\n"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)




    predictors.names <-
      c(namesof("gamma", .link.gamma, earg = .earg , tag = FALSE),
        if ( .delta.known) NULL else 
        namesof("delta", "identity", earg = list(), tag = FALSE))


    if (!length(etastart)) {
      delta.init <- if ( .delta.known) {
                     if (min(y,na.rm = TRUE) <= .delta)
                         stop("delta must be < min(y)")
                     .delta 
                   } else {
                     if (length( .idelta)) .idelta else
                         min(y,na.rm = TRUE) - 1.0e-4 *
                         diff(range(y,na.rm = TRUE))
                   }
      gamma.init <- if (length( .igamma)) .igamma else
                   median(y - delta.init) # = 1/median(1/(y-delta.init))
      gamma.init <- rep(gamma.init, length = length(y))
      etastart <-
        cbind(theta2eta(gamma.init, .link.gamma , earg = .earg ),
                        if ( .delta.known) NULL else delta.init)
                       
    }
  }), list( .link.gamma = link.gamma, .earg = earg,
            .delta.known = delta.known,
            .delta = delta,
            .idelta = idelta,
            .igamma = igamma ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta <- as.matrix(eta)
    mygamma <- eta2theta(eta[, 1], .link.gamma, earg = .earg )
    delta <- if ( .delta.known) .delta else eta[, 2]


    NA * mygamma
  }, list( .link.gamma = link.gamma, .earg = earg,
           .delta.known = delta.known,
           .delta = delta ))),
  last = eval(substitute(expression({
    misc$link <- if ( .delta.known) NULL else c(delta = "identity")
    misc$link <- c(gamma = .link.gamma, misc$link)
    misc$earg <- if ( .delta.known) list(gamma = .earg ) else
                list(gamma = .earg , delta = list())
    if ( .delta.known)
      misc$delta <- .delta
  }), list( .link.gamma = link.gamma, .earg = earg,
            .delta.known = delta.known,
            .delta = delta ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    eta <- as.matrix(eta)
    mygamma <- eta2theta(eta[, 1], .link.gamma, earg = .earg )
    delta <- if ( .delta.known) .delta else eta[, 2]
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else
    sum(c(w) * 0.5 * (log(mygamma) -3 * log(y - delta) -
                      mygamma / (y - delta)))
  }, list( .link.gamma = link.gamma, .earg = earg,
           .delta.known = delta.known,
           .delta = delta ))),
  vfamily = c("levy"),
  deriv = eval(substitute(expression({
    eta <- as.matrix(eta)
    mygamma <- eta2theta(eta[, 1], .link.gamma, earg = .earg )
    delta <- if ( .delta.known) .delta else eta[, 2]
    if (! .delta.known)
      dl.ddelta  <- (3 - mygamma / (y-delta)) / (2 * (y-delta))
    dl.dgamma <- 0.5 * (1 / mygamma - 1 / (y-delta))
    dgamma.deta <- dtheta.deta(mygamma, .link.gamma, earg = .earg )
    c(w) * cbind(dl.dgamma * dgamma.deta, 
                 if ( .delta.known) NULL else dl.ddelta)
  }), list( .link.gamma = link.gamma, .earg = earg,
            .delta.known = delta.known,
            .delta = delta ))),
  weight = eval(substitute(expression({
    wz <- matrix(as.numeric(NA), n, dimm(M)) # M = if (delta is known) 1 else 2
    wz[, iam(1, 1, M)] <- 1 * dgamma.deta^2 
    if (! .delta.known) {
      wz[, iam(1, 2, M)] <-  3 * dgamma.deta
      wz[, iam(2, 2, M)] <-  21
    }
    wz <- c(w) * wz / (2 * mygamma^2) 
    wz
  }), list( .link.gamma = link.gamma, .earg = earg,
           .delta.known = delta.known,
           .delta = delta ))))
}


        



dlino <- function(x, shape1, shape2, lambda = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  loglik <-  dbeta(x = x, shape1 = shape1, shape2 = shape2, log = TRUE) +
            shape1 * log(lambda) -
            (shape1+shape2) * log1p(-(1-lambda)*x)
  if (log.arg) loglik else exp(loglik)
}


plino <- function(q, shape1, shape2, lambda = 1) {
  ans <- pbeta(q = lambda * q / (1 - (1-lambda)*q),
              shape1 = shape1, shape2 = shape2)
  ans[lambda <= 0] <- NaN
  ans
}


qlino <- function(p, shape1, shape2, lambda = 1) {
  Y <- qbeta(p = p, shape1 = shape1, shape2 = shape2)
  ans <- Y / (lambda + (1-lambda)*Y)
  ans[lambda <= 0] <- NaN
  ans
}


rlino <- function(n, shape1, shape2, lambda = 1) {
  Y <- rbeta(n = n, shape1 = shape1, shape2 = shape2)
  ans <- Y / (lambda + (1 - lambda) * Y)
  ans[lambda <= 0] <- NaN
  ans
}



 lino <- function(lshape1 = "loge",
                  lshape2 = "loge",
                  llambda = "loge",
                  ishape1 = NULL, ishape2 = NULL, ilambda = 1,
                  zero = NULL) {

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
      stop("bad input for argument 'zero'")
  if (!is.Numeric(ilambda, positive = TRUE))
      stop("bad input for argument 'ilambda'")



  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")

  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")


  new("vglmff",
  blurb = c("Generalized Beta distribution (Libby and Novick, 1982)\n\n",
            "Links:    ",
            namesof("shape1", lshape1, earg = eshape1), ", ", 
            namesof("shape2", lshape2, earg = eshape2), ", ", 
            namesof("lambda", llambda, earg = elambda), "\n", 
            "Mean:     something complicated"),
  constraints = eval(substitute(expression({
      constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    if (min(y) <= 0 || max(y) >= 1)
      stop("values of the response must be between 0 and 1 (0,1)")

    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1)




    predictors.names <-
      c(namesof("shape1", .lshape1 , earg = .eshape1 , tag = FALSE),
        namesof("shape2", .lshape2 , earg = .eshape2 , tag = FALSE),
        namesof("lambda", .llambda , earg = .elambda , tag = FALSE))




      if (!length(etastart)) {
        lambda.init <- rep(if (length( .ilambda )) .ilambda else 1,
                          length = n)
        sh1.init <- if (length( .ishape1 ))
                     rep( .ishape1, length = n) else NULL
        sh2.init <- if (length( .ishape2 ))
                     rep( .ishape2, length = n) else NULL
            txY.init <- lambda.init * y / (1+lambda.init*y - y)
            mean1 <- mean(txY.init)
            mean2 <- mean(1/txY.init)
            if (!is.Numeric(sh1.init))
                sh1.init <- rep((mean2 - 1) / (mean2 - 1/mean1), length = n)
            if (!is.Numeric(sh2.init))
                sh2.init <- rep(sh1.init * (1-mean1) / mean1, length = n)
            etastart <-
              cbind(theta2eta(sh1.init, .lshape1 , earg = .eshape1),
                    theta2eta(sh2.init, .lshape2 , earg = .eshape2),
                    theta2eta(lambda.init, .llambda , earg = .elambda ))
        }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
            .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda,
            .ishape1=ishape1, .ishape2=ishape2, .ilambda = ilambda ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    sh1 <- eta2theta(eta[, 1], .lshape1 , earg = .eshape1)
    sh2 <- eta2theta(eta[, 2], .lshape2 , earg = .eshape2)
    lambda <- eta2theta(eta[, 3], .llambda , earg = .elambda )
    rep(as.numeric(NA), length = nrow(eta))
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
           .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1 = .lshape1 , shape2 = .lshape2 ,
                     lambda = .llambda )
    misc$earg <- list(shape1 = .eshape1 , shape2 = .eshape2 ,
                     lambda = .elambda )
  }), list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
            .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    sh1 <- eta2theta(eta[, 1], .lshape1 , earg = .eshape1)
    sh2 <- eta2theta(eta[, 2], .lshape2 , earg = .eshape2)
    lambda <- eta2theta(eta[, 3], .llambda , earg = .elambda )
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {
      sum(c(w) * dlino(y, shape1 = sh1, shape2 = sh2,
                       lambda = lambda, log = TRUE))
    }
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
           .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda ))),
  vfamily = c("lino"),
  deriv = eval(substitute(expression({
    sh1 <- eta2theta(eta[, 1], .lshape1 , earg = .eshape1)
    sh2 <- eta2theta(eta[, 2], .lshape2 , earg = .eshape2)
    lambda <- eta2theta(eta[, 3], .llambda , earg = .elambda )

    temp1 <- log1p(-(1-lambda) * y)
    temp2 <- digamma(sh1+sh2)

    dl.dsh1 <- log(lambda) + log(y) - digamma(sh1) + temp2 - temp1
    dl.dsh2 <- log1p(-y) - digamma(sh2) + temp2 - temp1
    dl.dlambda <- sh1/lambda - (sh1+sh2) * y / (1 - (1-lambda) * y)

    dsh1.deta <- dtheta.deta(sh1, .lshape1 , earg = .eshape1)
    dsh2.deta <- dtheta.deta(sh2, .lshape2 , earg = .eshape2)
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )

    c(w) * cbind( dl.dsh1    * dsh1.deta,
                  dl.dsh2    * dsh2.deta,
                  dl.dlambda * dlambda.deta)
  }), list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
            .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda ))),
  weight = eval(substitute(expression({
    temp3 <- trigamma(sh1+sh2)

    ned2l.dsh1 <- trigamma(sh1) - temp3
    ned2l.dsh2 <- trigamma(sh2) - temp3
    ned2l.dlambda2 <- sh1 * sh2 / (lambda^2 * (sh1+sh2+1))
    ned2l.dsh1sh2 <- -temp3
    ned2l.dsh1lambda <- -sh2 / ((sh1+sh2)*lambda)
    ned2l.dsh2lambda <-  sh1 / ((sh1+sh2)*lambda)

    wz <- matrix(as.numeric(NA), n, dimm(M))  #M==3 means 6=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.dsh1 * dsh1.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dsh2 * dsh2.deta^2
    wz[, iam(3, 3, M)] <- ned2l.dlambda2 * dlambda.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dsh1sh2 * dsh1.deta * dsh2.deta
    wz[, iam(1, 3, M)] <- ned2l.dsh1lambda * dsh1.deta * dlambda.deta
    wz[, iam(2, 3, M)] <- ned2l.dsh2lambda * dsh2.deta * dlambda.deta
    wz <- c(w) * wz
    wz
  }), list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
            .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda ))))
}






 betaprime <- function(link = "loge", i1 = 2, i2 = NULL, zero = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Beta-prime distribution\n",
            "y^(shape1-1) * (1+y)^(-shape1-shape2) / Beta(shape1,shape2),",
            " y>0, shape1>0, shape2>0\n\n",
            "Links:    ",
            namesof("shape1", link, earg = earg),  ", ",
            namesof("shape2", link, earg = earg), "\n",
            "Mean:     shape1/(shape2-1) provided shape2>1"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
      c(namesof("shape1", .link , earg = .earg , short = TRUE),
        namesof("shape2", .link , earg = .earg , short = TRUE))
    if (is.numeric( .i1) && is.numeric( .i2)) {
      vec <- c( .i1, .i2)
      vec <- c(theta2eta(vec[1], .link , earg = .earg ),
              theta2eta(vec[2], .link , earg = .earg ))
      etastart <- matrix(vec, n, 2, byrow = TRUE)
    }
    if (!length(etastart)) {
      init1 <- if (length( .i1)) 
        rep( .i1, length.out = n) else rep(1, length.out = n)
      init2 <- if (length( .i2))
        rep( .i2, length.out = n) else 1 + init1 / (y + 0.1)
      etastart <-
        matrix(theta2eta(c(init1, init2), .link , earg = .earg ),
               n, 2, byrow = TRUE)
    }
  }), list( .link = link, .earg = earg, .i1 = i1, .i2 = i2 ))), 

  linkinv = eval(substitute(function(eta, extra = NULL) {
      shapes <- eta2theta(eta, .link , earg = .earg )
    ifelse(shapes[, 2] > 1, shapes[, 1] / (shapes[, 2]-1), NA)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <- c(shape1 = .link , shape2 = .link)
    misc$earg <- list(shape1 = .earg , shape2 = .earg )
  }), list( .link = link, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL){
    shapes <- eta2theta(eta, .link , earg = .earg )
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {
        sum(c(w) *((shapes[, 1]-1) * log(y) -
                 lbeta(shapes[, 1], shapes[, 2]) -
                (shapes[, 2]+shapes[, 1]) * log1p(y)))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = "betaprime",
  deriv = eval(substitute(expression({
    shapes <- eta2theta(eta, .link , earg = .earg )
    dshapes.deta <- dtheta.deta(shapes, .link , earg = .earg )
    dl.dshapes <- cbind(log(y) - log1p(y) - digamma(shapes[, 1]) + 
                       digamma(shapes[, 1]+shapes[, 2]),
                       - log1p(y) - digamma(shapes[, 2]) + 
                       digamma(shapes[, 1]+shapes[, 2]))
    c(w) * dl.dshapes * dshapes.deta
  }), list( .link = link, .earg = earg ))),
  weight = expression({
    temp2 <- trigamma(shapes[, 1] + shapes[, 2])
    d2l.dshape12 <- temp2 - trigamma(shapes[, 1])
    d2l.dshape22 <- temp2 - trigamma(shapes[, 2])
    d2l.dshape1shape2 <- temp2

    wz <- matrix(as.numeric(NA), n, dimm(M))   #3=dimm(M)
    wz[, iam(1, 1, M)] <- d2l.dshape12 * dshapes.deta[, 1]^2
    wz[, iam(2, 2, M)] <- d2l.dshape22 * dshapes.deta[, 2]^2
    wz[, iam(1, 2, M)] <- d2l.dshape1shape2 *
                      dshapes.deta[, 1] * dshapes.deta[, 2]

    -c(w) * wz
  }))
}






dmaxwell <- function(x, a, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  L <- max(length(x), length(a))
  x <- rep(x, length.out = L); a = rep(a, length.out = L);
  logdensity <- rep(log(0), length.out = L)
  xok <- (x > 0)
  logdensity[xok] <- 0.5 * log(2/pi) + 1.5 * log(a[xok]) +
                     2 * log(x[xok]) - 0.5 * a[xok] * x[xok]^2
  logdensity[a <= 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}


pmaxwell <- function(q, a) {
  L <- max(length(q), length(a))
  q <- rep(q, length.out = L);
  a <- rep(a, length.out = L); 
  ans <- ifelse(q > 0,
                erf(q*sqrt(a/2)) - q*exp(-0.5*a*q^2) * sqrt(2*a/pi),
                0)
  ans[a <= 0] <- NaN
  ans
}


rmaxwell <- function(n, a) {

  sqrt(2 * rgamma(n = n, 1.5) / a)
}


qmaxwell <- function(p, a) {
  if (!is.Numeric(p, positive = TRUE) || any(p >= 1)) 
    stop("bad input for argument 'p'")
  if (any(a <= 0))
    stop("argument 'a' must be positive")

  N <- max(length(p), length(a));
  p <- rep(p, length.out = N);
  a <- rep(a, length.out = N)
  sqrt(2 * qgamma(p = p, 1.5) / a)
}




 maxwell <- function(link = "loge", zero = NULL) {


  link <- as.list(substitute(link)) # orig
  earg <- link2list(link)
  link <- attr(earg, "function.name")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  new("vglmff",
  blurb = c("Maxwell distribution f(y;a) = sqrt(2/pi) * a^(3/2) * y^2 *",
            " exp(-0.5*a*y^2), y>0, a>0\n",
            "Link:    ",
            namesof("a", link, earg = earg),
            "\n", "\n",
            "Mean:    sqrt(8 / (a * pi))"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 1,
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
    Musual <- 1
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    mynames1  <- paste("a", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
      namesof(mynames1, .link , earg = .earg )


    if (!length(etastart)) {
      a.init <- 8 / (pi * (y + 0.1)^2)
      etastart <- theta2eta(a.init, .link , earg = .earg )
    }
  }), list( .link = link,
            .earg = earg ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    aa <- eta2theta(eta, .link , earg = .earg )
    sqrt(8 / (aa * pi))
  }, list( .link = link,
           .earg = earg ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for(ilocal in 1:ncoly) {
      misc$earg[[ilocal]] <- .earg
    }

    misc$link <- rep( .link , length = ncoly)
    names(misc$link) <- mynames1

    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .link = link, .earg = earg ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    aa <- eta2theta(eta, .link , earg = .earg )
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else
      sum(c(w) * dmaxwell(x = y, a = aa, log = TRUE))
  }, list( .link = link,
           .earg = earg ))),
  vfamily = c("maxwell"),
  deriv = eval(substitute(expression({
    aa <- eta2theta(eta, .link , earg = .earg )

    dl.da <- 1.5 / aa - 0.5 * y^2

    da.deta <- dtheta.deta(aa, .link , earg = .earg )

    c(w) * dl.da * da.deta
  }), list( .link = link,
            .earg = earg ))),
  weight = eval(substitute(expression({
    ned2l.da2 <- 1.5 / aa^2
    wz <- c(w) * ned2l.da2 * da.deta^2
    wz
  }), list( .link = link, .earg = earg ))))
}







dnaka <- function(x, shape, scale = 1, log = FALSE) {
    if (!is.logical(log.arg <- log) || length(log) != 1)
      stop("bad input for argument 'log'")
    rm(log)

    L <- max(length(x), length(shape), length(scale))
    x     <- rep(x,     length.out = L)
    shape <- rep(shape, length.out = L)
    scale <- rep(scale, length.out = L);

    logdensity <- rep(log(0), length.out = L)
    xok <- (x > 0)
    logdensity[xok] <- dgamma(x = x[xok]^2, shape = shape[xok],
                             scale = scale[xok]/shape[xok], log = TRUE) +
                      log(2) + log(x[xok])
    if (log.arg) logdensity else exp(logdensity)
}


pnaka <- function(q, shape, scale = 1) {
    if (!is.Numeric(q))
        stop("bad input for argument 'q'")
    if (!is.Numeric(shape, positive = TRUE))
        stop("bad input for argument 'shape'")
    if (!is.Numeric(scale, positive = TRUE))
        stop("bad input for argument 'scale'")
    L <- max(length(q), length(shape), length(scale))
    q     <- rep(q,     length.out = L)
    shape <- rep(shape, length.out = L)
    scale <- rep(scale, length.out = L);
    ifelse(q <= 0, 0, pgamma(shape * q^2 / scale, shape))
}


qnaka <- function(p, shape, scale = 1, ...) {
  if (!is.Numeric(p, positive = TRUE) || max(p) >= 1)
    stop("bad input for argument 'p'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")
  if (!is.Numeric(scale, positive = TRUE))
    stop("bad input for argument 'scale'")

  L <- max(length(p), length(shape), length(scale))
  p <- rep(p, length.out = L); shape = rep(shape, length.out = L);
  scale <- rep(scale, length.out = L);
  ans <- rep(0.0, length.out = L)

  myfun <- function(x, shape, scale = 1, p)
    pnaka(q = x, shape = shape, scale = scale) - p
  for(ii in 1:L) {
    EY <- sqrt(scale[ii]/shape[ii]) *
          gamma(shape[ii] + 0.5) / gamma(shape[ii])
    Upper <- 5 * EY
    while(pnaka(q = Upper, shape = shape[ii],
                           scale = scale[ii]) < p[ii])
        Upper <- Upper + scale[ii]
    ans[ii] <- uniroot(f = myfun, lower = 0, upper = Upper,
                       shape = shape[ii], scale = scale[ii],
                       p = p[ii], ...)$root
  }
  ans
}


rnaka <- function(n, shape, scale = 1, Smallno = 1.0e-6) {

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (!is.Numeric(scale, positive = TRUE, allowable.length = 1))
    stop("bad input for argument 'scale'")
  if (!is.Numeric(shape, positive = TRUE, allowable.length = 1))
    stop("bad input for argument 'shape'")
  if (!is.Numeric(Smallno, positive = TRUE, allowable.length = 1) ||
      Smallno > 0.01 ||
      Smallno < 2 * .Machine$double.eps)
    stop("bad input for argument 'Smallno'")
  ans <- rep(0.0, length.out = use.n)

  ptr1 <- 1; ptr2 <- 0
  ymax <- dnaka(x = sqrt(scale * (1 - 0.5 / shape)),
               shape = shape, scale = scale)
  while(ptr2 < use.n) {
    EY <- sqrt(scale / shape) * gamma(shape + 0.5) / gamma(shape)
    Upper <- EY + 5 * scale
    while(pnaka(q = Upper, shape = shape, scale = scale) < 1 - Smallno)
      Upper <- Upper + scale
    x <- runif(2*use.n, min = 0, max = Upper)
    index <- runif(2*use.n, max = ymax) < dnaka(x, shape = shape,
                                               scale = scale)
    sindex <- sum(index)
    if (sindex) {
      ptr2 <- min(use.n, ptr1 + sindex - 1)
      ans[ptr1:ptr2] <- (x[index])[1:(1+ptr2-ptr1)]
      ptr1 <- ptr2 + 1
    }
  }
  ans
}






 nakagami <- function(lshape = "loge", lscale = "loge",
                      ishape = NULL, iscale = 1) {

  if (!is.null(iscale) && !is.Numeric(iscale, positive = TRUE))
    stop("argument 'iscale' must be a positive number or NULL")


  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")


  new("vglmff",
  blurb = c("Nakagami distribution f(y) = 2 * (shape/scale)^shape *\n",
            "                             ",
            "y^(2*shape-1) * exp(-shape*y^2/scale) / gamma(shape),\n",
            "                             ",
            "y>0, shape>0, scale>0\n",
            "Links:    ",
            namesof("shape", lshape, earg = eshape), ", ",
            namesof("scale", lscale, earg = escale),
            "\n",
            "\n",
            "Mean:    sqrt(scale/shape) * gamma(shape+0.5) / gamma(shape)"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
      c(namesof("shape", .lshape , earg = .eshape, tag = FALSE),
        namesof("scale", .lscale , earg = .escale , tag = FALSE))


    if (!length(etastart)) {
        init2 <- if (is.Numeric( .iscale, positive = TRUE))
                    rep( .iscale, length.out = n) else
                    rep(1, length.out = n)
        init1 <- if (is.Numeric( .ishape, positive = TRUE))
                    rep( .ishape, length.out = n) else
                rep(init2 / (y+1/8)^2, length.out = n)
        etastart <-
          cbind(theta2eta(init1, .lshape , earg = .eshape ),
                theta2eta(init2, .lscale , earg = .escale ))
    }
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .ishape = ishape, .iscale = iscale ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        shape <- eta2theta(eta[, 1], .lshape , earg = .eshape )
        scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
        sqrt(scale/shape) * gamma(shape+0.5) / gamma(shape)
    }, list( .lscale = lscale, .lshape = lshape,
             .escale = escale, .eshape = eshape))),
    last = eval(substitute(expression({
        misc$link <-    c(shape = .lshape , scale = .lscale)
        misc$earg <- list(shape = .eshape, scale = .escale )
        misc$expected = TRUE
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        shape <- eta2theta(eta[, 1], .lshape , earg = .eshape )
        scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else
            sum(c(w) * dnaka(x = y, shape = shape, scale = scale, log = TRUE))
    }, list( .lscale = lscale, .lshape = lshape,
             .escale = escale, .eshape = eshape))),
    vfamily = c("nakagami"),
    deriv = eval(substitute(expression({
        shape <- eta2theta(eta[, 1], .lshape , earg = .eshape )
        Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
        dl.dshape <- 1 + log(shape/Scale) - digamma(shape) +
                    2 * log(y) - y^2 / Scale
        dl.dscale <- -shape/Scale + shape * (y/Scale)^2
        dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
        dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
        c(w) * cbind(dl.dshape * dshape.deta,
                     dl.dscale * dscale.deta)
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape))),
    weight = eval(substitute(expression({
        d2l.dshape2 <- trigamma(shape) - 1/shape
        d2l.dscale2 <- shape / Scale^2
        wz <- matrix(as.numeric(NA), n, M)  # diagonal
        wz[, iam(1, 1, M)] <- d2l.dshape2 * dshape.deta^2
        wz[, iam(2, 2, M)] <- d2l.dscale2 * dscale.deta^2
        c(w) * wz
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape))))
}



drayleigh <- function(x, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  L <- max(length(x), length(scale))
  x <- rep(x, length.out = L); scale = rep(scale, length.out = L);
  logdensity <- rep(log(0), length.out = L)
  xok <- (x > 0)
  logdensity[xok] <- log(x[xok]) - 0.5 * (x[xok]/scale[xok])^2 -
                    2 * log(scale[xok])
  if (log.arg) logdensity else exp(logdensity)
}


prayleigh <- function(q, scale = 1) {
  if (any(scale <= 0))
    stop("argument 'scale' must be positive")
  L <- max(length(q), length(scale)) 
  q <- rep(q, length.out = L); scale = rep(scale, length.out = L);
  ifelse(q > 0,  -expm1(-0.5*(q/scale)^2), 0)
}


qrayleigh <- function(p, scale = 1) {
  if (any(p <= 0) || any(p >= 1))
    stop("argument 'p' must be between 0 and 1")
  ans <- scale * sqrt(-2 * log1p(-p))
  ans[scale <= 0] <- NaN
  ans
}


rrayleigh <- function(n, scale = 1) {
  ans <- scale * sqrt(-2 * log(runif(n)))
  ans[scale <= 0] <- NaN
  ans
}



 rayleigh <- function(lscale = "loge",
                      nrfs = 1 / 3 + 0.01,
                      oim.mean = TRUE, zero = NULL) {
  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")


  if (!is.Numeric(nrfs, allowable.length = 1) ||
      nrfs < 0 ||
      nrfs > 1)
    stop("bad input for 'nrfs'")

  if (!is.logical(oim.mean) || length(oim.mean) != 1)
    stop("bad input for argument 'oim.mean'")


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")



  new("vglmff",
  blurb = c("Rayleigh distribution\n\n",
            "f(y) = y*exp(-0.5*(y/scale)^2)/scale^2, y>0, scale>0\n\n",
            "Link:    ",
            namesof("scale", lscale, earg = escale), "\n\n",
            "Mean:    scale * sqrt(pi / 2)"),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 1
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 1,
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
    Musual <- 1
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    mynames1  <- paste("scale", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
      namesof(mynames1, .lscale , earg = .escale , tag = FALSE)


    if (!length(etastart)) {
      Ymat <- matrix(colSums(y) / colSums(w), n, ncoly, byrow = TRUE)
      b.init <- (Ymat + 1/8) / sqrt(pi/2)
      etastart <- theta2eta(b.init, .lscale , earg = .escale )
    }
  }), list( .lscale = lscale, .escale = escale))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale <- eta2theta(eta, .lscale , earg = .escale )
    Scale * sqrt(pi / 2)
  }, list( .lscale = lscale, .escale = escale))),

  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <- c(rep( .lscale , length = ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for(ii in 1:ncoly) {
      misc$earg[[ii]] <- .escale
    }

    misc$Musual <- Musual
    misc$multipleResponses <- TRUE
    misc$nrfs <- .nrfs
  }), list( .lscale = lscale,
            .escale = escale, .nrfs = nrfs  ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Scale <- eta2theta(eta, .lscale , earg = .escale )
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {
        sum(c(w) * drayleigh(x = y, scale = Scale, log = TRUE))
    }
  }, list( .lscale = lscale, .escale = escale))),

  vfamily = c("rayleigh"),
  deriv = eval(substitute(expression({
    Scale <- eta2theta(eta, .lscale , earg = .escale )

    dl.dScale <- ((y/Scale)^2 - 2) / Scale

    dScale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )

    c(w) * dl.dScale * dScale.deta
  }), list( .lscale = lscale, .escale = escale))),

  weight = eval(substitute(expression({
    d2l.dScale2 <- (3 * (y/Scale)^2 - 2) / Scale^2
    ned2l.dScale2 <- 4 / Scale^2
    wz <- c(w) * dScale.deta^2 *
         ((1 - .nrfs) * d2l.dScale2 + .nrfs * ned2l.dScale2)


    if (intercept.only && .oim.mean ) {
      ave.oim <- weighted.mean(d2l.dScale2, w)
      if (ave.oim > 0) {
        wz <- c(w) * dScale.deta^2 * ave.oim
      }
    }

    wz
  }), list( .lscale = lscale,
            .escale = escale,
            .nrfs = nrfs, .oim.mean = oim.mean ))))
}





dparetoIV <- function(x, location = 0, scale = 1, inequality = 1,
                      shape = 1, log = FALSE) {
    if (!is.logical(log.arg <- log) || length(log) != 1)
      stop("bad input for argument 'log'")
    rm(log)

    N <- max(length(x), length(location), length(scale),
            length(inequality), length(shape))
    x <- rep(x, length.out = N);
    location <- rep(location, length.out = N)
    scale <- rep(scale, length.out = N);
    inequality <- rep(inequality, length.out = N)
    shape <- rep(shape, length.out = N)

    logdensity <- rep(log(0), length.out = N)
    xok <- (x > location)
    zedd <- (x - location) / scale
    logdensity[xok] <- log(shape[xok]) -
                      log(scale[xok]) -  log(inequality[xok]) +
                      (1/inequality[xok]-1) * log(zedd[xok]) - 
                      (shape[xok]+1) *
                      log1p(zedd[xok]^(1/inequality[xok]))
    if (log.arg) logdensity else exp(logdensity)
}


pparetoIV <-
  function(q, location = 0, scale = 1, inequality = 1, shape = 1) {
  if (!is.Numeric(q))
    stop("bad input for argument 'q'")
  if (!is.Numeric(scale, positive = TRUE)) 
    stop("bad input for argument 'scale'")
  if (!is.Numeric(inequality, positive = TRUE)) 
    stop("bad input for argument 'inequality'")
  if (!is.Numeric(shape, positive = TRUE)) 
    stop("bad input for argument 'shape'")

  N <- max(length(q), length(location), length(scale),
          length(inequality), length(shape))
  q <- rep(q, length.out = N);
  location <- rep(location, length.out = N)
  scale <- rep(scale, length.out = N);
  inequality <- rep(inequality, length.out = N)
  shape <- rep(shape, length.out = N)
  answer <- q * 0
  ii <- q > location
  zedd <- (q[ii] - location[ii]) / scale[ii]
  answer[ii] <- 1 - (1 + zedd^(1/inequality[ii]))^(-shape[ii])
  answer
}


qparetoIV <-
  function(p, location = 0, scale = 1, inequality = 1, shape = 1) {
  if (!is.Numeric(p, positive = TRUE) || any(p >= 1)) 
    stop("bad input for argument 'p'")
  if (!is.Numeric(inequality, positive = TRUE)) 
    stop("bad input for argument 'inequality'")
  if (!is.Numeric(shape, positive = TRUE)) 
    stop("bad input for argument 'shape'")
  ans <- location + scale * (-1 + (1-p)^(-1/shape))^inequality
  ans[scale <= 0] <- NaN
  ans[shape <= 0] <- NaN
  ans
}


rparetoIV =
  function(n, location = 0, scale = 1, inequality = 1, shape = 1) {
  if (!is.Numeric(inequality, positive = TRUE)) 
    stop("bad input for argument 'inequality'")
  ans <- location + scale * (-1 + runif(n)^(-1/shape))^inequality
  ans[scale <= 0] <- NaN
  ans[shape <= 0] <- NaN
  ans
}


dparetoIII <- function(x, location = 0, scale = 1, inequality = 1,
                       log = FALSE)
  dparetoIV(x = x, location = location, scale = scale,
            inequality = inequality, shape = 1, log = log)

pparetoIII <- function(q, location = 0, scale = 1, inequality = 1)
  pparetoIV(q = q, location = location, scale = scale,
            inequality = inequality, shape = 1)

qparetoIII <- function(p, location = 0, scale = 1, inequality = 1)
  qparetoIV(p = p, location = location, scale = scale,
            inequality = inequality, shape = 1)

rparetoIII <- function(n, location = 0, scale = 1, inequality = 1)
  rparetoIV(n = n, location= location, scale = scale,
            inequality = inequality, shape = 1)



dparetoII <- function(x, location = 0, scale = 1, shape = 1, log = FALSE)
  dparetoIV(x = x, location = location, scale = scale,
            inequality = 1, shape = shape,
            log = log)

pparetoII <- function(q, location = 0, scale = 1, shape = 1)
  pparetoIV(q = q, location = location, scale = scale,
            inequality = 1, shape = shape)

qparetoII <- function(p, location = 0, scale = 1, shape = 1)
  qparetoIV(p = p, location = location, scale = scale,
            inequality = 1, shape = shape)

rparetoII <- function(n, location = 0, scale = 1, shape = 1)
  rparetoIV(n = n, location = location, scale = scale,
            inequality = 1, shape = shape)


dparetoI <- function(x, scale = 1, shape = 1)
  dparetoIV(x = x, location = scale, scale = scale, inequality = 1,
            shape = shape)

pparetoI <- function(q, scale = 1, shape = 1)
  pparetoIV(q = q, location = scale, scale = scale, inequality = 1,
            shape = shape)

qparetoI <- function(p, scale = 1, shape = 1)
  qparetoIV(p = p, location = scale, scale = scale, inequality = 1,
            shape = shape)

rparetoI <- function(n, scale = 1, shape = 1)
  rparetoIV(n = n, location = scale, scale = scale, inequality = 1,
            shape = shape)



 paretoIV <- function(location = 0,
                      lscale = "loge",
                      linequality = "loge",
                      lshape = "loge",
                      iscale = 1, iinequality = 1, ishape = NULL,
                      imethod = 1) {

  if (!is.Numeric(location))
    stop("argument 'location' must be numeric")
  if (is.Numeric(iscale) && any(iscale <= 0))
    stop("argument 'iscale' must be positive")
  if (is.Numeric(iinequality) && any(iinequality <= 0))
    stop("argument 'iinequality' must be positive")
  if (is.Numeric(ishape) && any(ishape <= 0))
    stop("argument 'ishape' must be positive")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE) ||
      imethod > 2)
    stop("bad input for argument 'imethod'")

  if (linequality == "nloge" && location != 0)
      warning("The Burr distribution has 'location = 0' and ",
              "'linequality = nloge'")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  linequ <- as.list(substitute(linequality))
  einequ <- link2list(linequ)
  linequ <- attr(einequ, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  iinequ <- iinequality


  new("vglmff",
  blurb = c("Pareto(IV) distribution F(y)=1-[1+((y - ", location,
            ")/scale)^(1/inequality)]^(-shape),",
            "\n",
            "         y > ",
            location,
            ", scale > 0, inequality > 0, shape > 0,\n",
            "Links:    ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("inequality", linequ, earg = einequ ),
            ", ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:    location + scale * NA"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
     c(namesof("scale", .lscale , earg = .escale , tag = FALSE),
       namesof("inequality", .linequ ,
               earg = .einequ , tag = FALSE),
       namesof("shape", .lshape , earg = .eshape, tag = FALSE))



    extra$location <- location <- .location
    if (any(y <= location))
      stop("the response must have values > than the 'location' argument")

    if (!length(etastart)) {
      inequ.init <- if (length( .iinequ )) .iinequ else  1
      scale.init <- if (length( .iscale )) .iscale else 1
      shape.init <- if (length( .ishape )) .ishape else NULL

      if (!length(shape.init)) {
        zedd <- (y - location) / scale.init
        if ( .imethod == 1) {
          A1 <- weighted.mean(1/(1 + zedd^(1/inequ.init)), w = w)
          A2 <- weighted.mean(1/(1 + zedd^(1/inequ.init))^2, w = w)
        } else {
          A1 <- median(1/(1 + zedd^(1/inequ.init )))
          A2 <- median(1/(1 + zedd^(1/inequ.init))^2)
        }
        shape.init <- max(0.01, (2*A2-A1)/(A1-A2))
      }

          etastart <- cbind(
            theta2eta(rep(scale.init, length.out = n),
                      .lscale , earg = .escale ),
            theta2eta(rep(inequ.init, length.out = n),
                      .linequ, earg = .einequ),
            theta2eta(rep(shape.init, length.out = n),
                      .lshape , earg = .eshape ))
      }
  }), list( .location = location, .lscale = lscale,
      .linequ = linequ, .lshape = lshape, .imethod = imethod,
      .escale = escale, .einequ = einequ, .eshape = eshape,
      .iscale = iscale, .iinequ = iinequ, .ishape = ishape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
      location <- extra$location
      Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
      inequ <- eta2theta(eta[, 2], .linequ, earg = .einequ)
      shape <- eta2theta(eta[, 3], .lshape , earg = .eshape )
      location + Scale * NA
  }, list( .lscale = lscale, .linequ = linequ, .lshape = lshape,
           .escale = escale, .einequ = einequ, .eshape = eshape))),
    last = eval(substitute(expression({
        misc$link <-    c("scale"      = .lscale ,
                          "inequality" = .linequ,
                          "shape"      = .lshape)
        misc$earg <- list("scale"      = .escale ,
                          "inequality" = .einequ,
                          "shape"      = .eshape )
        misc$location = extra$location # Use this for prediction
    }), list( .lscale = lscale, .linequ = linequ,
              .escale = escale, .einequ = einequ,
              .lshape = lshape,
              .eshape = eshape))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        location <- extra$location
        Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
        inequ <- eta2theta(eta[, 2], .linequ, earg = .einequ)
        shape <- eta2theta(eta[, 3], .lshape , earg = .eshape )
        zedd <- (y - location) / Scale
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(c(w) * dparetoIV(x = y, location = location, scale = Scale,
                                 inequ = inequ, shape = shape,
                                 log = TRUE))
        }
    }, list( .lscale = lscale, .linequ = linequ,
             .escale = escale, .einequ = einequ,
             .lshape = lshape,
             .eshape = eshape))),
    vfamily = c("paretoIV"),
    deriv = eval(substitute(expression({
        location = extra$location
        Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
        inequ <- eta2theta(eta[, 2], .linequ, earg = .einequ)
        shape <- eta2theta(eta[, 3], .lshape , earg = .eshape )
        zedd <- (y - location) / Scale
        temp100 <- 1 + zedd^(1/inequ)
        dl.dscale <- (shape  - (1+shape) / temp100) / (inequ * Scale)
        dl.dinequ <- ((log(zedd) * (shape - (1+shape)/temp100)) /
                         inequ - 1) / inequ
        dl.dshape <- -log(temp100) + 1/shape
        dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
        dinequ.deta <- dtheta.deta(inequ, .linequ, earg = .einequ)
        dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
        c(w) * cbind(dl.dscale * dscale.deta,
                     dl.dinequ * dinequ.deta, 
                     dl.dshape * dshape.deta)
    }), list( .lscale = lscale, .linequ = linequ,
              .lshape = lshape,
              .escale = escale, .einequ = einequ,
              .eshape = eshape))),
    weight = eval(substitute(expression({
        temp200 <- digamma(shape) - digamma(1) - 1
        d2scale.deta2 <- shape / ((inequ*Scale)^2 * (shape+2))
        d2inequ.deta2 <- (shape * (temp200^2 + trigamma(shape) + trigamma(1)
                             ) + 2*(temp200+1)) / (inequ^2 * (shape+2))
        d2shape.deta2 <- 1 / shape^2
        d2si.deta2 <- (shape*(-temp200) -1) / (inequ^2 * Scale * (shape+2))
        d2ss.deta2 <- -1 / ((inequ*Scale) * (shape+1))
        d2is.deta2 <- temp200 / (inequ*(shape+1))
        wz <- matrix(0, n, dimm(M))
        wz[, iam(1, 1, M)] <- dscale.deta^2 * d2scale.deta2
        wz[, iam(2, 2, M)] <- dinequ.deta^2 * d2inequ.deta2
        wz[, iam(3, 3, M)] <- dshape.deta^2 * d2shape.deta2
        wz[, iam(1, 2, M)] <- dscale.deta * dinequ.deta * d2si.deta2
        wz[, iam(1, 3, M)] <- dscale.deta * dshape.deta * d2ss.deta2
        wz[, iam(2, 3, M)] <- dinequ.deta * dshape.deta * d2is.deta2
        c(w) * wz
    }), list( .lscale = lscale, .linequ = linequ, .lshape = lshape,
              .escale = escale, .einequ = einequ, .eshape = eshape))))
}




 paretoIII <- function(location = 0,
                       lscale = "loge",
                       linequality = "loge",
                       iscale = NULL, iinequality = NULL) {

  if (!is.Numeric(location))
      stop("argument 'location' must be numeric")
  if (is.Numeric(iscale) && any(iscale <= 0))
      stop("argument 'iscale' must be positive")
  if (is.Numeric(iinequality) && any(iinequality <= 0))
      stop("argument 'iinequality' must be positive")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  linequ <- as.list(substitute(linequality))
  einequ <- link2list(linequ)
  linequ <- attr(einequ, "function.name")


  iinequ <- iinequality



  new("vglmff",
  blurb = c("Pareto(III) distribution F(y)=1-[1+((y - ", location,
            ")/scale)^(1/inequality)]^(-1),",
            "\n", "         y > ",
            location, ", scale > 0, inequality > 0, \n",
            "Links:    ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("inequality", linequ, earg = einequ ), "\n",
            "Mean:    location + scale * NA"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
    c(namesof("scale", .lscale , earg = .escale , tag = FALSE),
      namesof("inequ", .linequ, earg = .einequ, tag = FALSE))
    extra$location = location = .location

    if (any(y <= location))
      stop("the response must have values > than the 'location' argument")


    if (!length(etastart)) {
            inequ.init <- if (length( .iinequ)) .iinequ else  NULL
            scale.init <- if (length( .iscale )) .iscale else NULL
            if (!length(inequ.init) || !length(scale.init)) {
                probs <- (1:4)/5
                ytemp <- quantile(x = log(y-location), probs = probs)
                fittemp <- lsfit(x = logit(probs), y = ytemp, intercept = TRUE)
                if (!length(inequ.init))
                    inequ.init <- max(fittemp$coef["X"], 0.01)
                if (!length(scale.init))
                    scale.init <- exp(fittemp$coef["Intercept"])
            }
            etastart=cbind(
            theta2eta(rep(scale.init, length.out = n),
                     .lscale , earg = .escale ),
            theta2eta(rep(inequ.init, length.out = n),
                      .linequ,
                      earg = .einequ))
        }
  }), list( .location = location, .lscale = lscale,
            .linequ = linequ,
            .escale = escale, .einequ = einequ,
            .iscale = iscale, .iinequ = iinequ ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    location <- extra$location
    Scale      <- eta2theta(eta[, 1], .lscale     , earg = .escale )
    inequ <- eta2theta(eta[, 2], .linequ, earg = .einequ)
    location + Scale * NA
  }, list( .lscale = lscale, .linequ = linequ,
           .escale = escale, .einequ = einequ ))),
  last = eval(substitute(expression({
    misc$link <-    c("scale" = .lscale , "inequality" = .linequ)
    misc$earg <- list("scale" = .escale , "inequality" = .einequ)

    misc$location <- extra$location # Use this for prediction
  }), list( .lscale = lscale, .linequ = linequ,
            .escale = escale, .einequ = einequ ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        location <- extra$location
        Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
        inequ <- eta2theta(eta[, 2], .linequ, earg = .einequ)
        zedd <- (y - location) / Scale
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(c(w) * dparetoIII(x = y, location = location, scale=Scale,
                               inequ=inequ, log = TRUE))
        }
    }, list( .lscale = lscale, .linequ = linequ,
             .escale = escale, .einequ = einequ ))),
    vfamily = c("paretoIII"),
    deriv = eval(substitute(expression({
        location <- extra$location
        Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
        inequ <- eta2theta(eta[, 2], .linequ, earg = .einequ)
        shape <- 1
        zedd <- (y - location) / Scale
        temp100 <- 1 + zedd^(1/inequ)
        dl.dscale <- (shape  - (1+shape) / temp100) / (inequ * Scale)
        dl.dinequ <- ((log(zedd) * (shape - (1+shape)/temp100)) /
                         inequ - 1) / inequ
        dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
        dinequ.deta <- dtheta.deta(inequ, .linequ, earg = .einequ)
        c(w) * cbind(dl.dscale * dscale.deta,
                     dl.dinequ * dinequ.deta)
    }), list( .lscale = lscale, .linequ = linequ,
              .escale = escale, .einequ = einequ ))),
    weight = eval(substitute(expression({
        d2scale.deta2 <- 1 / ((inequ*Scale)^2 * 3)
        d2inequ.deta2 <- (1 + 2* trigamma(1)) / (inequ^2 * 3)
        wz <- matrix(0, n, M) # It is diagonal
        wz[, iam(1, 1, M)] <- dscale.deta^2 * d2scale.deta2
        wz[, iam(2, 2, M)] <- dinequ.deta^2 * d2inequ.deta2
        c(w) * wz
    }), list( .lscale = lscale, .linequ = linequ,
              .escale = escale, .einequ = einequ ))))
}





 paretoII <- function(location = 0,
                      lscale = "loge",
                      lshape = "loge",
                      iscale = NULL, ishape = NULL) {

  if (!is.Numeric(location))
    stop("argument 'location' must be numeric")
  if (is.Numeric(iscale) && any(iscale <= 0))
    stop("argument 'iscale' must be positive")
  if (is.Numeric(ishape) && any(ishape <= 0))
    stop("argument 'ishape' must be positive")


  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")




  new("vglmff",
  blurb = c("Pareto(II) distribution F(y)=1-[1+(y - ", location,
            ")/scale]^(-shape),",
            "\n", "         y > ",
            location, ", scale > 0,  shape > 0,\n",
            "Links:    ", namesof("scale", lscale, earg = escale), ", ",
                          namesof("shape", lshape, earg = eshape), "\n",
            "Mean:    location + scale * NA"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
      c(namesof("scale", .lscale , earg = .escale , tag = FALSE),
        namesof("shape", .lshape , earg = .eshape , tag = FALSE))

    extra$location <- location <- .location

    if (any(y <= location))
      stop("the response must have values > than the 'location' argument")

    if (!length(etastart)) {
            scale.init <- if (length( .iscale )) .iscale else NULL
            shape.init <- if (length( .ishape )) .ishape else  NULL
            if (!length(shape.init) || !length(scale.init)) {
                probs <- (1:4)/5
                scale.init.0 <- 1
                ytemp <- quantile(x = log(y-location+scale.init.0),
                                 probs = probs)
                fittemp <- lsfit(x = log1p(-probs), y = ytemp,
                                intercept = TRUE)
                if (!length(shape.init))
                    shape.init <- max(-1/fittemp$coef["X"], 0.01)
                if (!length(scale.init))
                    scale.init <- exp(fittemp$coef["Intercept"])
            }
            etastart=cbind(
            theta2eta(rep(scale.init, length.out = n),
                      .lscale , earg = .escale ),
            theta2eta(rep(shape.init, length.out = n),
                      .lshape , earg = .eshape ))
        }
    }), list( .location = location, .lscale = lscale,
              .escale = escale, .eshape = eshape, 
              .lshape = lshape, .iscale = iscale, .ishape = ishape ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        location <- extra$location
        Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
        shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
        location + Scale * NA
    }, list( .lscale = lscale, .lshape = lshape,
             .escale = escale, .eshape = eshape))),
    last = eval(substitute(expression({
        misc$link <-    c("scale" = .lscale , "shape" = .lshape)

        misc$earg <- list("scale" = .escale , "shape" = .eshape )

        misc$location <- extra$location # Use this for prediction
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        location <- extra$location
        Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
        shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
        zedd <- (y - location) / Scale
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(c(w) * dparetoII(x = y, location = location, scale=Scale,
                              shape = shape, log = TRUE))
        }
    }, list( .lscale = lscale, .lshape = lshape,
             .escale = escale, .eshape = eshape))),
    vfamily = c("paretoII"),
    deriv = eval(substitute(expression({
        location <- extra$location
        Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
        shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
        zedd <- (y - location) / Scale
        temp100 <- 1 + zedd
        dl.dscale <- (shape  - (1+shape) / temp100) / (1 * Scale)
        dl.dshape <- -log(temp100) + 1/shape
        dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
        dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
        c(w) * cbind(dl.dscale * dscale.deta,
                     dl.dshape * dshape.deta)
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape))),
    weight = eval(substitute(expression({
        d2scale.deta2 <- shape / (Scale^2 * (shape+2))
        d2shape.deta2 <- 1 / shape^2
        d2ss.deta2 <- -1 / (Scale * (shape+1))
        wz <- matrix(0, n, dimm(M))
        wz[, iam(1, 1, M)] <- dscale.deta^2 * d2scale.deta2
        wz[, iam(2, 2, M)] <- dshape.deta^2 * d2shape.deta2
        wz[, iam(1, 2, M)] <- dscale.deta * dshape.deta * d2ss.deta2
        c(w) * wz
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape))))
}





dpareto <- function(x, location, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  L <- max(length(x), length(location), length(shape)) 
  x <- rep(x, length.out = L);
  location <- rep(location, length.out = L);
  shape <- rep(shape, length.out = L)

  logdensity <- rep(log(0), length.out = L)
  xok <- (x > location)
  logdensity[xok] <- log(shape[xok]) + shape[xok] * log(location[xok]) -
                      (shape[xok]+1) * log(x[xok])
  if (log.arg) logdensity else exp(logdensity)
}


ppareto <- function(q, location, shape) {

  L <- max(length(q), length(location), length(shape))
  q <- rep(q, length.out = L);
  location <- rep(location, length.out = L);
  shape <- rep(shape, length.out = L)

  ans <- ifelse(q > location, 1 - (location/q)^shape, 0)
  ans[location <= 0] <- NaN
  ans[shape    <= 0] <- NaN
  ans
}


qpareto <- function(p, location, shape) {
  if (any(p <= 0) || any(p >= 1))
    stop("argument 'p' must be between 0 and 1")

  ans <- location / (1 - p)^(1/shape)
  ans[location <= 0] <- NaN
  ans[shape    <= 0] <- NaN
  ans
}


rpareto <- function(n, location, shape) {
  ans <- location / runif(n)^(1/shape)
  ans[location <= 0] <- NaN
  ans[shape    <= 0] <- NaN
  ans
}



 pareto1 <- function(lshape = "loge", location = NULL) {
  if (is.Numeric(location) && location <= 0)
      stop("argument 'location' must be positive")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  earg <- eshape


  new("vglmff",
  blurb = c("Pareto distribution ",
            "f(y) = shape * location^shape / y^(shape+1),",
            " 0<location<y, shape>0\n",
            "Link:    ", namesof("shape", lshape, earg = earg), "\n", "\n",
            "Mean:    location*shape/(shape-1) for shape>1"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)


    predictors.names <-
      namesof("shape", .lshape , earg = .earg , tag = FALSE)


    locationhat <- if (!length( .location)) {
      locationEstimated <- TRUE
      min(y) # - .smallno
    } else {
      locationEstimated <- FALSE
      .location
    }
    if (any(y < locationhat))
      stop("the value of location is too high ",
           "(requires 0 < location < min(y))")
    extra$location <- locationhat
    extra$locationEstimated <- locationEstimated

    if (!length(etastart)) {
        k.init <- (y + 1/8) / (y - locationhat + 1/8)
        etastart <- theta2eta(k.init, .lshape , earg = .earg )
    }
  }), list( .lshape = lshape, .earg = earg,
            .location = location ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    k <- eta2theta(eta, .lshape , earg = .earg )
    location <- extra$location
    ifelse(k > 1, k * location / (k-1), NA)
  }, list( .lshape = lshape, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c(k = .lshape)

    misc$earg <- list(k = .earg )

    misc$location <- extra$location # Use this for prediction
  }), list( .lshape = lshape, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    k <- eta2theta(eta, .lshape , earg = .earg )
    location <- extra$location
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {


      sum(c(w) * (log(k) + k * log(location) - (k+1) * log(y )))
    }
  }, list( .lshape = lshape, .earg = earg ))),
  vfamily = c("pareto1"),
  deriv = eval(substitute(expression({
    location <- extra$location
    k <- eta2theta(eta, .lshape , earg = .earg )
    dl.dk <- 1/k + log(location/y)
    dk.deta <- dtheta.deta(k, .lshape , earg = .earg )
    c(w) * dl.dk * dk.deta
  }), list( .lshape = lshape, .earg = earg ))),
  weight = eval(substitute(expression({
    ed2l.dk2 <- 1 / k^2
    wz <- c(w) * dk.deta^2 * ed2l.dk2
    wz
  }), list( .lshape = lshape, .earg = earg ))))
}





dtpareto <- function(x, lower, upper, shape, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(x))
    stop("bad input for argument 'x'")
  if (!is.Numeric(lower, positive = TRUE))
    stop("argument 'lower' must be positive")
  if (!is.Numeric(upper, positive = TRUE))
    stop("argument 'upper' must be positive")
  if (!is.Numeric(shape, positive = TRUE))
    stop("argument 'shape' must be positive")

  L <- max(length(x), length(lower), length(upper), length(shape))
  x <- rep(x, length.out = L);
  shape <- rep(shape, length.out = L)
  lower <- rep(lower, length.out = L);
  upper <- rep(upper, length.out = L);


  logdensity <- rep(log(0), length.out = L)
  xok <- (0 < lower) & (lower < x) & (x < upper) & (shape > 0)

  logdensity[xok] <- log(shape[xok]) + shape[xok] * log(lower[xok]) -
                     (shape[xok] + 1) * log(x[xok]) -
                     log1p(-(lower[xok] / upper[xok])^(shape[xok]))

  logdensity[shape <= 0] <- NaN
  logdensity[upper < lower] <- NaN
  logdensity[0 > lower] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}


ptpareto <- function(q, lower, upper, shape) {
  if (!is.Numeric(q))
    stop("bad input for argument 'q'")

  L <- max(length(q), length(lower), length(upper), length(shape)) 
  q <- rep(q, length.out = L);
  lower <- rep(lower, length.out = L);
  upper <- rep(upper, length.out = L);
  shape <- rep(shape, length.out = L)

  ans <- q * 0
  xok <- (0 < lower) & (lower < q) & (q < upper) & (shape > 0)
  ans[xok] <- (1 - (lower[xok]/q[xok])^shape[xok]) / (1 -
                  (lower[xok]/upper[xok])^shape[xok])
  ans[q >= upper] <- 1

  ans[upper < lower] <- NaN
  ans[lower <= 0] <- NaN
  ans[upper <= 0] <- NaN
  ans[shape <= 0] <- NaN

  ans
}


qtpareto <- function(p, lower, upper, shape) {
  if (!is.Numeric(p, positive = TRUE))
    stop("bad input for argument 'p'")
  if (max(p) >= 1)
    stop("argument 'p' must be in (0, 1)")

  ans <- lower / (1 - p*(1-(lower/upper)^shape))^(1/shape)
  ans[lower <= 0] <- NaN
  ans[upper <= 0] <- NaN
  ans[shape <= 0] <- NaN
  ans[upper <  lower] <- NaN
  ans
}


rtpareto <- function(n, lower, upper, shape) {

  ans <- qtpareto(p = runif(n), lower = lower, upper = upper, shape = shape)
  ans[lower <= 0] <- NaN
  ans[upper <= 0] <- NaN
  ans[shape <= 0] <- NaN
  ans
}




 tpareto1 <- function(lower, upper, lshape = "loge",
                      ishape = NULL, imethod = 1) {

  if (!is.Numeric(lower, positive = TRUE, allowable.length = 1))
    stop("bad input for argument 'lower'")
  if (!is.Numeric(upper, positive = TRUE, allowable.length = 1))
    stop("bad input for argument 'upper'")
  if (lower >= upper)
    stop("lower < upper is required")

  if (length(ishape) && !is.Numeric(ishape, positive = TRUE))
    stop("bad input for argument 'ishape'")



  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")
  earg <- eshape


  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  new("vglmff",
  blurb = c("Truncated Pareto distribution f(y) = shape * lower^shape /",
            "(y^(shape+1) * (1-(lower/upper)^shape)),",
            " 0 < lower < y < upper < Inf, shape>0\n",
            "Link:    ", namesof("shape", lshape, earg = earg), "\n", "\n",
            "Mean:    shape*lower^shape*(upper^(1-shape)-lower^(1-shape)) /",
                      " ((1-shape) * (1-(lower/upper)^shape))"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)




    predictors.names <- namesof("shape", .lshape , earg = .earg ,
                                tag = FALSE)
    if (any(y <= .lower))
      stop("the value of argument 'lower' is too high ",
           "(requires '0 < lower < min(y)')")

    extra$lower <- .lower
    if (any(y >= .upper))
        stop("the value of argument 'upper' is too low ",
             "(requires 'max(y) < upper')")
    extra$upper <- .upper

    if (!length(etastart)) {
      shape.init <- if (is.Numeric( .ishape )) 0 * y + .ishape else
      if ( .imethod == 2) {
        0 * y + median(rep((y + 1/8) / (y - .lower + 1/8), times = w))
      } else {
        tpareto1.Loglikfun <- function(shape, y, x, w, extraargs) {
           myratio <- .lower / .upper
           sum(c(w) * (log(shape) + shape * log( .lower) -
                    (shape+1) * log(y) - log1p(-myratio^shape)))
        }
        shape.grid <- 2^((-4):4)
        try.this <- getMaxMin(shape.grid, objfun = tpareto1.Loglikfun,
                             y = y,  x = x, w = w)
        try.this = rep(try.this, length.out = n)
        try.this
      }
      etastart <- theta2eta(shape.init, .lshape , earg = .earg )
    }
  }), list( .lshape = lshape, .earg = earg,
            .ishape = ishape,
            .imethod = imethod,
            .lower = lower, .upper = upper ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .earg )
    myratio <- .lower / .upper
    constprop <- shape * .lower^shape / (1 - myratio^shape)
    constprop * ( .upper^(1-shape) - .lower^(1-shape)) / (1-shape)
  }, list( .lshape = lshape, .earg = earg,
             .lower = lower, .upper = upper ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape = .lshape)

    misc$earg <- list(shape = .earg )

    misc$lower <- extra$lower
    misc$upper <- extra$upper
    misc$expected <- TRUE
  }), list( .lshape = lshape, .earg = earg,
            .lower = lower, .upper = upper ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .earg )
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {
      ans <- sum(c(w) * dtpareto(x = y, lower = .lower ,
                                upper = .upper ,
                                shape = shape, log = TRUE))
      ans
    }
  }, list( .lshape = lshape, .earg = earg,
           .lower = lower, .upper = upper ))),
  vfamily = c("tpareto1"),
  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .earg )
    myratio <- .lower / .upper
    myratio2 <-  myratio^shape
    tmp330 <- myratio2 * log(myratio) / (1 - myratio2)

    dl.dshape <- 1 / shape + log( .lower) - log(y) + tmp330 

    dshape.deta <- dtheta.deta(shape, .lshape , earg = .earg )

    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .earg = earg,
            .lower = lower, .upper = upper ))),
  weight = eval(substitute(expression({
    ned2l.dshape2 <- 1 / shape^2 - tmp330^2 / myratio2
    wz <- c(w) * dshape.deta^2 * ned2l.dshape2
    wz
  }), list( .lshape = lshape, .earg = earg,
            .lower = lower, .upper = upper ))))
}




erf <- function(x)
    2 * pnorm(x * sqrt(2)) - 1

erfc <- function(x)
    2 * pnorm(x * sqrt(2), lower.tail = FALSE)



 wald <- function(link.lambda = "loge", init.lambda = NULL) {

  link.lambda <- as.list(substitute(link.lambda))
  earg <- link2list(link.lambda)
  link.lambda <- attr(earg, "function.name")



  new("vglmff",
  blurb = c("Standard Wald distribution\n\n",
           "f(y) = sqrt(lambda/(2*pi*y^3)) * ",
           "exp(-lambda*(y-1)^2/(2*y)), y&lambda>0",
           "\n", 
           "Link:     ",
                         namesof("lambda", link.lambda, earg = earg), "\n",
           "Mean:     ", "1\n",
           "Variance: 1 / lambda"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1)


    predictors.names <-
      namesof("lambda", .link.lambda, earg = .earg , short = TRUE)


    if (!length(etastart)) {
      initlambda <- if (length( .init.lambda)) .init.lambda else
                   1 / (0.01 + (y-1)^2)
      initlambda <- rep(initlambda, length.out = n)
      etastart <-
        cbind(theta2eta(initlambda,
                        link = .link.lambda , earg = .earg ))
      }
  }), list( .link.lambda = link.lambda, .earg = earg,
           .init.lambda=init.lambda ))),
  linkinv = function(eta, extra = NULL) {
      0*eta + 1
  },
    last = eval(substitute(expression({
        misc$link <-    c(lambda = .link.lambda )

        misc$earg <- list(lambda = .earg )

    }), list( .link.lambda = link.lambda, .earg = earg ))),
    loglikelihood = eval(substitute(
             function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        lambda <- eta2theta(eta, link=.link.lambda, earg = .earg )
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else
        sum(c(w) * (0.5 * log(lambda/(2*pi*y^3)) - lambda * (y-1)^2 / (2*y)))
    }, list( .link.lambda = link.lambda, .earg = earg ))),
    vfamily = "wald",
    deriv = eval(substitute(expression({
        lambda <- eta2theta(eta, link=.link.lambda, earg = .earg )
        dl.dlambda <- 0.5 / lambda + 1 - 0.5 * (y + 1/y)
        dlambda.deta <- dtheta.deta(theta=lambda, link=.link.lambda, earg = .earg )
        c(w) * cbind(dl.dlambda * dlambda.deta)
    }), list( .link.lambda = link.lambda, .earg = earg ))),
    weight = eval(substitute(expression({
        d2l.dlambda2 <- 0.5 / (lambda^2)
        c(w) * cbind(dlambda.deta^2 * d2l.dlambda2)
    }), list( .link.lambda = link.lambda, .earg = earg ))))
}


 expexp <- function(lshape = "loge", lscale = "loge",
                    ishape = 1.1, iscale = NULL, # ishape cannot be 1
                    tolerance = 1.0e-6,
                    zero = NULL) {


  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")

  if (!is.Numeric(tolerance, positive = TRUE, allowable.length = 1) ||
      tolerance > 1.0e-2)
    stop("bad input for argument 'tolerance'")
  if (!is.Numeric(ishape, positive = TRUE))
      stop("bad input for argument 'ishape'")

  if (length(iscale) && !is.Numeric(iscale, positive = TRUE))
      stop("bad input for argument 'iscale'")

  ishape[ishape == 1] <- 1.1 # Fails in @deriv


  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")



  new("vglmff",
  blurb = c("Exponentiated Exponential Distribution\n",
             "Links:    ",
             namesof("shape", lshape, earg = eshape), ", ",
             namesof("scale", lscale, earg = escale),"\n",
             "Mean:     (digamma(shape+1)-digamma(1))/scale"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)



      predictors.names <-
      c(namesof("shape", .lshape , earg = .eshape, short = TRUE), 
        namesof("scale", .lscale , earg = .escale , short = TRUE))


      if (!length(etastart)) {
            shape.init <- if (!is.Numeric( .ishape, positive = TRUE))
                   stop("argument 'ishape' must be positive") else
                   rep( .ishape, length.out = n)
            scale.init <- if (length( .iscale ))
                        rep( .iscale, length.out = n) else
                        (digamma(shape.init+1) - digamma(1)) / (y+1/8)
            scale.init <- rep(weighted.mean(scale.init, w = w),
                             length.out = n)
            etastart <-
              cbind(theta2eta(shape.init, .lshape , earg = .eshape ),
                    theta2eta(scale.init, .lscale , earg = .escale ))
        }
    }), list( .lshape = lshape, .lscale = lscale,
              .iscale = iscale, .ishape = ishape,
              .eshape = eshape, .escale = escale))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        shape <- eta2theta(eta[, 1], .lshape , earg = .eshape )
        scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
        (digamma(shape+1)-digamma(1)) / scale
    }, list( .lshape = lshape, .lscale = lscale,
             .eshape = eshape, .escale = escale))),
  last = eval(substitute(expression({
    misc$link <-    c("shape" = .lshape , "scale" = .lscale )
    misc$earg <- list("shape" = .eshape , "scale" = .escale )

    misc$expected <- TRUE
  }), list( .lshape = lshape, .lscale = lscale,
            .eshape = eshape, .escale = escale))),
  loglikelihood= eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    shape <- eta2theta(eta[, 1], .lshape , earg = .eshape )
    scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else
    sum(c(w) * (log(shape) + log(scale) + 
             (shape-1)*log1p(-exp(-scale*y)) - scale*y))
  }, list( .lscale = lscale, .lshape = lshape,
           .eshape = eshape, .escale = escale))),
  vfamily = c("expexp"),
  deriv = eval(substitute(expression({
    shape <- eta2theta(eta[, 1], .lshape , earg = .eshape )
    scale <- eta2theta(eta[, 2], .lscale , earg = .escale )

    dl.dscale <- 1/scale + (shape-1)*y*exp(-scale*y) / (-expm1(-scale*y)) - y
    dl.dshape <- 1/shape + log1p(-exp(-scale*y))

    dscale.deta <- dtheta.deta(scale, .lscale , earg = .escale )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    c(w) * cbind(dl.dshape * dshape.deta,
                 dl.dscale * dscale.deta)
  }), list( .lshape = lshape, .lscale = lscale,
            .eshape = eshape, .escale = escale))),
  weight = eval(substitute(expression({
        d11 <- 1 / shape^2  # True for all shape
        d22 <- d12 <- rep(as.numeric(NA), length.out = n)
        index2 <- abs(shape - 2) > .tolerance  # index2 = shape != 1
        largeno <- 10000
        if (any(index2)) {
            Shape <- shape[index2]
            Shape[abs(Shape-1) < .tolerance] <- 1.001 # digamma(0) is undefined
            Scale <- scale[index2]
            tmp200 <- trigamma(1)-trigamma(Shape-1) +
                  (digamma(Shape-1)-digamma(1))^2    # Fails when Shape == 1
            tmp300 <- trigamma(1)-digamma(Shape)+(digamma(Shape)-digamma(1))^2
            d22[index2] <- (1 + Shape*(Shape-1)*tmp200/(Shape-2)) / Scale^2 +
                          Shape*tmp300 / Scale^2
        }
        if (any(!index2)) {
            Scale <- scale[!index2]
            d22[!index2] <- (1 + 4 * sum(1/(2 + (0:largeno))^3)) / Scale^2
        }

        index1 <- abs(shape - 1) > .tolerance  # index1 <- shape != 1
        if (any(index1)) {
            Shape <- shape[index1]
            Scale <- scale[index1]
            d12[index1] <- -(Shape*(digamma(Shape)-digamma(1))/(Shape-1) -
                          digamma(Shape+1) + digamma(1)) / Scale
        }
        if (any(!index1)) {
            Scale <- scale[!index1]
            d12[!index1] <- -sum(1/(2 + (0:largeno))^2) / Scale
        }
        wz <- matrix(0, n, dimm(M))
        wz[, iam(1, 1, M)] <- dshape.deta^2 * d11
        wz[, iam(2, 2, M)] <- dscale.deta^2 * d22
        wz[, iam(1, 2, M)] <- dscale.deta * dshape.deta * d12
        c(w) * wz
    }), list( .tolerance = tolerance ))))
}





 expexp1 <- function(lscale = "loge",
                     iscale = NULL,
                     ishape = 1) {

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")


  if (length(iscale) && !is.Numeric(iscale, positive = TRUE))
      stop("bad input for argument 'iscale'")



  new("vglmff",
  blurb = c("Exponentiated Exponential Distribution",
            " (profile likelihood estimation)\n",
            "Links:    ",
            namesof("scale", lscale, earg = escale), "\n",
            "Mean:     (digamma(shape+1)-digamma(1))/scale"),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)




    predictors.names <-
      namesof("scale", .lscale , earg = .escale , short = TRUE)

    if (length(w) != n ||
        !is.Numeric(w, integer.valued = TRUE, positive = TRUE))
      stop("argument 'weights' must be a vector of positive integers")

    if (!intercept.only)
      stop("this family function only works for an ",
           "intercept-only, i.e., y ~ 1")
    extra$yvector <- y
    extra$sumw <- sum(w)
    extra$w <- w

    if (!length(etastart)) {
      shape.init <- if (!is.Numeric( .ishape, positive = TRUE))
             stop("argument 'ishape' must be positive") else
             rep( .ishape, length.out = n)
      scaleinit <- if (length( .iscale ))
                  rep( .iscale, length.out = n) else
                  (digamma(shape.init+1) - digamma(1)) / (y+1/8)  
      etastart <-
        cbind(theta2eta(scaleinit, .lscale , earg = .escale ))
    }
  }), list( .lscale = lscale, .iscale = iscale, .ishape = ishape,
            .escale = escale))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    scale <- eta2theta(eta, .lscale , earg = .escale )
    temp7 <-  -expm1(-scale*extra$yvector)
    shape <- -extra$sumw / sum(extra$w*log(temp7)) # \gamma(\theta)
    (digamma(shape+1)-digamma(1)) / scale
  }, list( .lscale = lscale,
           .escale = escale))),
  last = eval(substitute(expression({
    misc$link <-    c("scale" = .lscale)
    misc$earg <- list("scale" = .escale )

    temp7 <-  -expm1(-scale*y)
    shape <- -extra$sumw / sum(w*log(temp7)) # \gamma(\theta)
    misc$shape <- shape   # Store the ML estimate here
    misc$pooled.weight <- pooled.weight
  }), list( .lscale = lscale, .escale = escale))),
  loglikelihood= eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    scale <- eta2theta(eta, .lscale , earg = .escale )
    temp7 <-  -expm1(-scale*y)
    shape <- -extra$sumw / sum(w*log(temp7)) # \gamma(\theta)
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else
    sum(c(w) * (log(shape) + log(scale) + 
             (shape-1)*log1p(-exp(-scale*y)) - scale*y))
  }, list( .lscale = lscale, .escale = escale))),
  vfamily = c("expexp1"),
  deriv = eval(substitute(expression({
    scale <- eta2theta(eta, .lscale , earg = .escale )

    temp6 <- exp(-scale*y)
    temp7 <- 1-temp6
    shape <- -extra$sumw / sum(w*log(temp7)) # \gamma(\theta)
    d1 <- 1/scale + (shape-1)*y*temp6/temp7 - y

    c(w) * cbind(d1 * dtheta.deta(scale, .lscale , earg = .escale ))
  }), list( .lscale = lscale, .escale = escale))),
  weight = eval(substitute(expression({
    d11 <- 1/scale^2  + y*(temp6/temp7^2) * ((shape-1) *
          (y*temp7+temp6) - y*temp6 / (log(temp7))^2)

    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M)] =
      dtheta.deta(scale, .lscale , earg = .escale )^2 * d11 -
      d2theta.deta2(scale, .lscale , earg = .escale ) * d1

    if (FALSE && intercept.only) {
      sumw <- sum(w)
      for(ii in 1:ncol(wz))
          wz[, ii] <- sum(wz[, ii]) / sumw
      pooled.weight <- TRUE
      wz <- c(w) * wz   # Put back the weights
    } else
      pooled.weight <- FALSE
    c(w) * wz
  }), list( .lscale = lscale, .escale = escale))))
}










 logistic2 <- function(llocation = "identity",
                       lscale = "loge",
                       ilocation = NULL, iscale = NULL,
                       imethod = 1, zero = -2) {

  ilocat <- ilocation


  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE))
    stop("bad input for argument 'zero'")
  if (length(iscale) && !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")



  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  new("vglmff",
  blurb = c("Two-parameter logistic distribution\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale",    lscale, earg = escale),
            "\n", "\n",
            "Mean:     location", "\n",
            "Variance: (pi * scale)^2 / 3"),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
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



    mynames1 <- paste("location", if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("scale",    if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .llocat , earg = .elocat , tag = FALSE),
          namesof(mynames2, .lscale , earg = .escale , tag = FALSE))[
          interleave.VGAM(M, M = Musual)]


    if (!length(etastart)) {
      if ( .imethod == 1) {
        locat.init <- y
        scale.init <- sqrt(3) * apply(y, 2, sd) / pi
      } else {
        locat.init <- scale.init <- NULL
        for(ii in 1:ncoly) {
          locat.init <- c(locat.init, median(rep(y[, ii], w[, ii])))
          scale.init <- c(scale.init, sqrt(3) * sum(w[, ii] *
                        (y[, ii] - locat.init[ii])^2) / (sum(w[, ii]) * pi))
        }
      }
      locat.init <- matrix(if (length( .ilocat )) .ilocat else
                          locat.init, n, ncoly, byrow = TRUE)
      if ( .llocat == "loge")
        locat.init <- abs(locat.init) + 0.001


      scale.init <- matrix(if (length( .iscale )) .iscale else
                          scale.init, n, ncoly, byrow = TRUE)

      etastart <- cbind(
        theta2eta(locat.init, .llocat , earg = .elocat ),
        theta2eta(scale.init, .lscale , earg = .escale ))[,
                        interleave.VGAM(M, M = Musual)]
    }
  }), list( .imethod = imethod,
            .elocat = elocat, .escale = escale,
            .llocat = llocat, .lscale = lscale,
            .ilocat = ilocat, .iscale = iscale ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    M <- ncol(eta)
    Musual <- 2
    ncoly <- M / Musual 
    eta2theta(eta[, (1:ncoly) * Musual - 1], .llocat , earg = .elocat )
  }, list( .llocat = llocat,
           .elocat = elocat ))),

  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <-
      c(rep( .llocat , length = ncoly),
        rep( .lscale , length = ncoly))[interleave.VGAM(M, M = Musual)]
    temp.names <- c(mynames1, mynames2)[
                    interleave.VGAM(M, M = Musual)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for(ii in 1:ncoly) {
      misc$earg[[Musual*ii-1]] <- .elocat
      misc$earg[[Musual*ii  ]] <- .escale
    }

    misc$Musual <- Musual
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .imethod = imethod,
             .llocat = llocat, .lscale = lscale,
             .elocat = elocat, .escale = escale))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    M <- ncol(eta)
    Musual <- 2
    ncoly <- M / Musual 

    locat <- eta2theta(eta[, (1:ncoly)*Musual-1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, (1:ncoly)*Musual  ], .lscale , earg = .escale )
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {
        sum(c(w) * dlogis(x = y, location = locat,
                          scale = Scale, log = TRUE))
    }
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale))),
  vfamily = c("logistic2"),
  deriv = eval(substitute(expression({
    Musual <- 2
    ncoly <- M / Musual 

    locat <- eta2theta(eta[, (1:ncoly)*Musual-1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, (1:ncoly)*Musual  ], .lscale , earg = .escale )

    zedd <- (y - locat) / Scale
    ezedd <- exp(-zedd)
    dl.dlocat <- (-expm1(-zedd)) / ((1 + ezedd) * Scale)
    dl.dscale <-  zedd * (-expm1(-zedd)) / ((1 + ezedd) * Scale) -
                 1 / Scale

    dlocat.deta <- dtheta.deta(locat, .llocat , earg = .elocat )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )

    c(w) * cbind(dl.dlocat * dlocat.deta,
                 dl.dscale * dscale.deta)[, interleave.VGAM(M, M = Musual)]
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale))),
  weight = eval(substitute(expression({
    ned2l.dlocat2 <- 1 / (3 * Scale^2)
    ned2l.dscale2 <- (3 + pi^2) / (9 * Scale^2)

    wz <- matrix(as.numeric(NA), nrow = n, ncol = M) # diagonal
    wz[, (1:ncoly) * Musual - 1] <- ned2l.dlocat2 * dlocat.deta^2
    wz[, (1:ncoly) * Musual    ] <- ned2l.dscale2 * dscale.deta^2

    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = ncoly)
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale))))
}







 negbinomial.size <- function(size = Inf,
                              lmu = "loge",
                              imu = NULL,
                              probs.y = 0.75,
                              imethod = 1,
                              shrinkage.init = 0.95, zero = NULL) {






  if (any(size <= 0))
    stop("bad input for argument 'size'")
  if (any(is.na(size)))
    stop("bad input for argument 'size'")



  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")




  if (length(imu) && !is.Numeric(imu, positive = TRUE))
    stop("bad input for argument 'imu'")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")
  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
     shrinkage.init < 0 ||
     shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")




  ans <- 
  new("vglmff",

  blurb = c("Negative-binomial distribution with size known\n\n",
            "Links:    ",
            namesof("mu",   lmu, earg = emu), "\n",
            "Mean:     mu\n",
            "Variance: mu * (1 + mu / size) for NB-2"),

  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)

  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 1,
         zero = .zero)
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({
    Musual <- 1

    if (any(y < 0))
      stop("negative values not allowed for the 'negbinomial.size' family")

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



    M <- Musual * ncol(y) 
    NOS <- ncoly <- ncol(y) # Number of species
    mynames1 <- paste("mu", if (NOS > 1) 1:NOS else "", sep = "")
    predictors.names <-
      namesof(mynames1, .lmu , earg = .emu , tag = FALSE)


    if (is.numeric( .mu.init ))
      MU.INIT <- matrix( .mu.init, nrow(y), ncol(y), byrow = TRUE)


    if (!length(etastart)) {
      mu.init <- y
      for(iii in 1:ncol(y)) {
        use.this <- if ( .imethod == 1) {
          weighted.mean(y[, iii], w[, iii]) + 1/16
        } else if ( .imethod == 3) {
          c(quantile(y[, iii], probs = .probs.y) + 1/16)
        } else {
          median(y[, iii]) + 1/16
        }

        if (is.numeric( .mu.init )) {
          mu.init[, iii] <- MU.INIT[, iii]
        } else {
          medabsres <- median(abs(y[, iii] - use.this)) + 1/32
          allowfun <- function(z, maxtol = 1)
            sign(z)*pmin(abs(z), maxtol)
          mu.init[, iii] <- use.this + (1 - .sinit) *
                           allowfun(y[, iii] - use.this,
                                    maxtol = medabsres)

          mu.init[, iii] <- abs(mu.init[, iii]) + 1 / 1024
        }
      } # of for(iii)


    kmat <- matrix( .size , n, NOS, byrow = TRUE)




    newemu <- .emu
    if ( .lmu == "nbcanlink") {
      newemu$size <- kmat
    }



    etastart <-
      cbind(theta2eta(mu.init , link = .lmu , earg = newemu ))
      }
  }), list( .lmu = lmu,
            .emu = emu,
            .mu.init = imu,
            .size = size, .probs.y = probs.y,
            .sinit = shrinkage.init,
            .zero = zero, .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Musual <- 1
    eta <- cbind(eta)
    NOS <- ncol(eta) / Musual
    n <- nrow(eta)
    kmat <- matrix( .size , n, NOS, byrow = TRUE)






    newemu <- .emu
    if ( .lmu == "nbcanlink") {
      newemu$size <- kmat
    }


    eta2theta(eta, .lmu , earg = newemu)
  }, list( .lmu = lmu,
           .emu = emu,
           .size = size ))),

  last = eval(substitute(expression({
    misc$link <- rep( .lmu , length = NOS)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for(ii in 1:NOS) {
      misc$earg[[ii]] <- newemu
    }


    misc$imethod <- .imethod 
    misc$expected <- TRUE
    misc$shrinkage.init <- .sinit
    misc$size <- kmat
  }), list( .lmu = lmu,
            .emu = emu,
            .sinit = shrinkage.init,
            .imethod = imethod ))),


  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    mu <- cbind(mu)
    y <- cbind(y)
    w <- cbind(w)
    eta <- cbind(eta)
    NOS <- ncol(eta)
    n   <- nrow(eta)
    kmat <- matrix( .size , n, NOS, byrow = TRUE)

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      ind1 <- is.finite(kmat)
      ans1 <- ans2 <- 0
      for (kk in 1:NOS) {
        ind1 <- is.finite(kmat[, kk])
        ans1 <- ans1 +
               sum(w[ind1] * dnbinom(x = y[ind1, kk], mu = mu[ind1, kk],
                            size = kmat[ind1, kk], log = TRUE))
        ans2 <- ans2 +
               sum(w[!ind1] * dpois(x = y[!ind1, kk], lambda  = mu[!ind1, kk],
                                    log = TRUE))
      }

      ans <- ans1 + ans2
      ans
    }
  }, list( .size = size ))),

  vfamily = c("negbinomial.size"),

  deriv = eval(substitute(expression({
    eta <- cbind(eta)
    NOS <- M <- ncol(eta)
    kmat <- matrix( .size , n, M, byrow = TRUE)



    newemu <- .emu
    if ( .lmu == "nbcanlink") {
      newemu$size <- kmat
    }


    dl.dmu <- y/mu - (y+kmat)/(kmat+mu)
    dl.dmu[!is.finite(dl.dmu)] <-  (y/mu)[!is.finite(dl.dmu)] - 1

    if ( .lmu == "nbcanlink")
      newemu$wrt.eta <- 1
    dmu.deta <- dtheta.deta(mu, .lmu , earg = newemu) # eta1

    myderiv <- c(w) * dl.dmu * dmu.deta
    myderiv
  }), list( .lmu = lmu, 
            .emu = emu,
           .size = size ))),

  weight = eval(substitute(expression({
    wz <- matrix(as.numeric(NA), n, M)  # wz is 'diagonal' 

    ned2l.dmu2 <- 1 / mu - 1 / (mu + kmat)
    wz <- dmu.deta^2 * ned2l.dmu2



    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }), list( .lmu = lmu ))))

  ans
}








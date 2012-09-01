# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.








VGAM.weights.function <- function(w, M, n) {


  ncolw = ncol(as.matrix(w))
  if (ncolw == 1) {
    wz = matrix(w, nrow = n, ncol = M) # w_i * diag(M)
  } else if (ncolw == M) {
    wz = as.matrix(w)
  } else if (ncolw < M && M > 1) {
    stop("ambiguous input for 'weights'")
  } else if (ncolw > M*(M+1)/2) {
    stop("too many columns")
  } else {
    wz = as.matrix(w)
  }
  wz
}












 gaussianff <- function(dispersion = 0, parallel = FALSE, zero = NULL)
{

  if (!is.Numeric(dispersion, allowable.length = 1) ||
      dispersion < 0)
    stop("bad input for argument 'dispersion'")
  estimated.dispersion <- dispersion == 0


  new("vglmff",
  blurb = c("Vector linear/additive model\n",
            "Links:    identity for Y1,...,YM"),
  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(1, M, 1), x, .parallel , constraints)
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .parallel = parallel, .zero = zero ))),
  deviance = function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    M <- if (is.matrix(y)) ncol(y) else 1
    n <- if (is.matrix(y)) nrow(y) else length(y)
    wz <- VGAM.weights.function(w = w, M = M, n = n)
    if (residuals) {
      if (M > 1) {
        U <- vchol(wz, M = M, n = n) 
        temp <- mux22(U, y-mu, M = M, upper = TRUE, as.matrix = TRUE)
        dimnames(temp) <- dimnames(y)
        temp
      } else (y-mu) * sqrt(wz)
    } else {
      rss.vgam(y-mu, wz = wz, M = M)
    }
  },

  infos = eval(substitute(function(...) {
    list(Musual = 1,
         zero = .zero)
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({
    if (is.R())
      assign("CQO.FastAlgorithm", TRUE, envir = VGAM::VGAMenv) else
      CQO.FastAlgorithm <<- TRUE
    if (any(function.name == c("cqo", "cao")) &&
       (length( .zero ) ||
       (is.logical( .parallel ) && .parallel )))
        stop("cannot handle non-default arguments for cqo() and cao()")

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    M = if (is.matrix(y)) ncol(y) else 1
    dy = dimnames(y)

    predictors.names <- if (!is.null(dy[[2]])) dy[[2]] else
                       paste("Y", 1:M, sep = "")

    if (!length(etastart)) 
      etastart = 0 * y
  }), list( .parallel = parallel, .zero = zero ))),
  linkinv = function(eta, extra = NULL) eta, 
  last = eval(substitute(expression({
    dy = dimnames(y)
    if (!is.null(dy[[2]]))
        dimnames(fit$fitted.values) = dy
    dpar = .dispersion
    if (!dpar) {
      wz = VGAM.weights.function(w = w, M = M, n = n)
      temp5 = rss.vgam(y-mu, wz = wz, M = M)
        dpar = temp5 / (length(y) -
        (if(is.numeric(ncol(X_vlm_save))) ncol(X_vlm_save) else 0))
    }
    misc$dispersion = dpar
    misc$default.dispersion = 0
    misc$estimated.dispersion = .estimated.dispersion

    misc$link = rep("identity", length = M)
    names(misc$link) = predictors.names
    misc$earg = vector("list", M)
    for (ilocal in 1:M)
      misc$earg[[ilocal]] <- list()
    names(misc$link) = predictors.names


    if (is.R()) {
      if (exists("CQO.FastAlgorithm", envir = VGAM::VGAMenv))
        rm("CQO.FastAlgorithm", envir = VGAM::VGAMenv)
    } else {
      while (exists("CQO.FastAlgorithm"))
        remove("CQO.FastAlgorithm")
    }

    misc$expected = TRUE
    misc$multipleResponses <- TRUE
  }), list( .dispersion = dispersion,
            .estimated.dispersion = estimated.dispersion ))),
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    M = if (is.matrix(y)) ncol(y) else 1
    n = if (is.matrix(y)) nrow(y) else length(y)
    wz = VGAM.weights.function(w = w, M = M, n = n)
    temp1 = rss.vgam(y-mu, wz = wz, M = M)



    if (M == 1 || ncol(wz) == M) {
      -0.5 * temp1 + 0.5 * sum(log(wz)) - n * (M / 2) * log(2*pi)
    } else {
      if (all(wz[1, ] == apply(wz, 2, min)) &&
          all(wz[1, ] == apply(wz, 2, max))) {
        onewz = m2adefault(wz[1, , drop = FALSE], M = M)
        onewz = onewz[, ,1]  # M x M

        logdet <- sum(log(eigen(onewz, symmetric = TRUE,
                                only.values = TRUE)$values))
        logretval <- -0.5 * temp1 + 0.5 * n * logdet -
                     n * (M / 2) * log(2*pi)
      logretval
    } else {
      logretval = -0.5 * temp1 - n * (M / 2) * log(2*pi)
      for (ii in 1:n) {
        onewz = m2adefault(wz[ii, , drop = FALSE], M = M)
        onewz = onewz[, ,1]  # M x M
          logdet <- sum(log(eigen(onewz, symmetric = TRUE,
                                  only.values = TRUE)$values))
          logretval = logretval + 0.5 * logdet
        }
        logretval
      }
    }
  },
  linkfun = function(mu, extra = NULL) mu,
  vfamily = "gaussianff",
  deriv = expression({
    wz = VGAM.weights.function(w = w, M = M, n = n)
    mux22(cc = t(wz), xmat = y-mu, M = M, as.matrix = TRUE)
  }),
  weight = expression({
    wz
  }))
}










dposnorm <- function(x, mean = 0, sd = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  L = max(length(x), length(mean), length(sd))
  x = rep(x, len = L);
  mean = rep(mean, len = L);
  sd = rep(sd, len = L);

  if (log.arg) {
    ifelse(x < 0, log(0), dnorm(x, mean = mean, sd = sd, log = TRUE) -
           pnorm(mean / sd, log.p = TRUE))
  } else {
    ifelse(x < 0, 0, dnorm(x = x, mean = mean, sd = sd) / pnorm(mean / sd))
  }
}


pposnorm <- function(q, mean = 0, sd = 1) {
  L = max(length(q), length(mean), length(sd))
  q = rep(q, len = L);
  mean = rep(mean, len = L);
  sd = rep(sd, len = L);
  ifelse(q < 0, 0, (pnorm(q, mean = mean, sd = sd) -
                    pnorm(0, mean = mean, sd = sd)) / pnorm(q = mean/sd))
}


qposnorm <- function(p, mean = 0, sd = 1) {
  if (!is.Numeric(p, positive = TRUE) || max(p) >= 1)
    stop("bad input for argument 'p'")
  qnorm(p = p + (1-p) * pnorm(0, mean = mean, sd = sd),
        mean = mean, sd = sd)
}


rposnorm <- function(n, mean = 0, sd = 1) {
  if (!is.Numeric(n, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'n'")
  mean = rep(mean, length = n)
  sd = rep(sd, length = n)
  qnorm(p = runif(n, min = pnorm(0, mean = mean, sd = sd)),
        mean = mean, sd = sd)
}



 posnormal1.control <- function(save.weight = TRUE, ...) {
    list(save.weight=save.weight)
}




 posnormal1 <- function(lmean = "identity", lsd = "loge",
                        imean = NULL, isd = NULL,
                        nsimEIM = 100, zero = NULL)
{
 warning("this VGAM family function is not working properly yet")


  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")



  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")
  if (length(isd) &&
      !is.Numeric(isd, positive = TRUE))
    stop("bad input for argument 'isd'")


  if (length(nsimEIM))
    if (!is.Numeric(nsimEIM, allowable.length = 1,
                    integer.valued = TRUE) ||
        nsimEIM <= 10)
      stop("argument 'nsimEIM' should be an integer greater than 10")


  new("vglmff",
  blurb = c("Positive (univariate) normal distribution\n\n",
          "Links:    ",
          namesof("mean", lmean, earg = emean, tag = TRUE), "; ",
          namesof("sd",   lsd,   earg = esd,   tag = TRUE)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(Musual = 2,
         par.names = c("mean", "sd"),
         zero = .zero )
  }, list( .zero = zero
         ))),




  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      c(namesof("mean", .lmean, earg = .emean, tag = FALSE),
        namesof("sd",   .lsd,   earg = .esd,   tag = FALSE))

    if (!length(etastart)) {
        init.me = if (length( .i.mean)) rep( .i.mean, len = n) else NULL
        init.sd = if (length( .i.sd  )) rep( .i.sd  , len = n) else NULL
            if (!length(init.me))
              init.me = rep(quantile(y, probs=0.40), len = n)
            if (!length(init.sd))
          init.sd = rep(sd(c(y)) * 1.2, len = n)
        etastart = cbind(theta2eta(init.me, .lmean, earg = .emean),
                         theta2eta(init.sd, .lsd,   earg = .esd ))
    }
  }), list( .lmean = lmean, .lsd = lsd,
            .i.mean = imean, .i.sd = isd,
            .emean = emean, .esd = esd
           ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mymu = eta2theta(eta[, 1], .lmean, earg = .emean )
    mysd = eta2theta(eta[, 2], .lsd,   earg = .esd )
    mymu + mysd * dnorm(-mymu/mysd) / pnorm(mymu/mysd)
  }, list( .lmean = lmean, .lsd = lsd,
           .emean = emean, .esd = esd
         ))),
  last = eval(substitute(expression({
    misc$link =    c("mean" = .lmean , "sd" = .lsd )
    misc$earg = list("mean" = .emean , "sd" = .esd )
    misc$expected = TRUE
    misc$nsimEIM = .nsimEIM
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    mymu = eta2theta(eta[, 1], .lmean, earg = .emean)
    mysd = eta2theta(eta[, 2], .lsd,   earg = .esd )
    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {

      sum(c(w) * dposnorm(x=y, m=mymu, sd = mysd, log = TRUE))
    }
  }, list( .lmean = lmean, .lsd = lsd,
           .emean = emean, .esd = esd ))),
  vfamily = c("posnormal1"),
  deriv = eval(substitute(expression({
    mymu = eta2theta(eta[, 1], .lmean, earg = .emean)
    mysd = eta2theta(eta[, 2], .lsd,  earg = .esd )

    zedd = (y-mymu) / mysd
    temp7 = dnorm(-mymu/mysd)
    temp8 = pnorm(mymu/mysd) * mysd

    dl.dmu = zedd / mysd^2 - temp7 / temp8
    dl.dsd = (mymu*temp7/temp8 + zedd^3 / mysd - 1) / mysd

    dmu.deta = dtheta.deta(mymu, .lmean, earg = .emean)
    dsd.deta = dtheta.deta(mysd, .lsd, earg = .esd )
    dthetas.detas = cbind(dmu.deta, dsd.deta)
    c(w) * dthetas.detas * cbind(dl.dmu, dl.dsd)
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd ))),
  weight = eval(substitute(expression({
    run.varcov = 0
    ind1 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    if (length( .nsimEIM )) {
        for(ii in 1:( .nsimEIM )) {
          ysim <- rposnorm(n, m=mymu, sd = mysd)
          zedd = (ysim-mymu) / mysd
          temp7 = dnorm(-mymu/mysd)
          temp8 = pnorm(mymu/mysd) * mysd
          dl.dmu = zedd / mysd^2 - temp7 / temp8
          dl.dsd = (mymu*temp7/temp8 + zedd^3 / mysd - 1) / mysd

          rm(ysim)
          temp3 = matrix(c(dl.dmu, dl.dsd), n, 2)
          run.varcov = ((ii-1) * run.varcov +
               temp3[, ind1$row.index]*temp3[, ind1$col.index]) / ii
      }
        wz = if (intercept.only)
            matrix(colMeans(run.varcov),
                   n, ncol(run.varcov), byrow = TRUE) else run.varcov

      wz = wz * dthetas.detas[, ind1$row] * dthetas.detas[, ind1$col]
      wz = c(w) * matrix(wz, n, dimm(M))
    } else {
      wz = matrix(as.numeric(NA), n, dimm(M))
      ed2l.dmu2 = (1 - temp7*mymu/temp8) / mysd^2  - (temp7/temp8)^2
      ed2l.dmusd = (temp7 /(mysd * temp8)) * (1 + (mymu/mysd)^2 +
                   mymu*temp7 / temp8)
      ed2l.dsd2 = 2 / mysd^2  - (temp7 * mymu /(mysd^2 * temp8)) *
                  (1 + (mymu/mysd)^2 + mymu*temp7/temp8)
      wz[, iam(1, 1, M)] = ed2l.dmu2  * dmu.deta^2
      wz[, iam(2, 2, M)] = ed2l.dsd2  * dsd.deta^2
      wz[, iam(1, 2, M)] = ed2l.dmusd * dsd.deta * dmu.deta
      wz = c(w) * wz
    }
    wz
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd,
            .nsimEIM = nsimEIM ))))
}




dbetanorm <- function(x, shape1, shape2, mean = 0, sd = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  ans =
  if (log.arg) {
    dnorm(x = x, mean = mean, sd = sd, log = TRUE) +
    (shape1-1) * pnorm(q = x, mean = mean, sd = sd, log.p = TRUE) +
    (shape2-1) * pnorm(q = x, mean = mean, sd = sd, log.p = TRUE,
                       lower.tail = FALSE) -
    lbeta(shape1, shape2)
  } else {
    dnorm(x = x, mean = mean, sd = sd) *
    pnorm(q = x, mean = mean, sd = sd)^(shape1-1) *
    pnorm(q = x, mean = mean, sd = sd,
          lower.tail = FALSE)^(shape2-1) / beta(shape1, shape2)
  }

  ans
}




pbetanorm <- function(q, shape1, shape2, mean = 0, sd = 1,
    lower.tail = TRUE, log.p = FALSE) {
    pbeta(q=pnorm(q = q, mean = mean, sd = sd),
                  shape1=shape1, shape2=shape2,
                  lower.tail = lower.tail, log.p = log.p)
}


qbetanorm <- function(p, shape1, shape2, mean = 0, sd = 1) {
  if (!is.Numeric(p, positive = TRUE) || max(p) >= 1)
    stop("bad input for argument 'p'")
  qnorm(p = qbeta(p = p, shape1 = shape1, shape2 = shape2),
        mean = mean, sd = sd)
}


rbetanorm <- function(n, shape1, shape2, mean = 0, sd = 1) {
  if (!is.Numeric(n, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'n'")
  qnorm(p = qbeta(p = runif(n), shape1 = shape1, shape2 = shape2),
        mean = mean, sd = sd)
}




dtikuv <- function(x, d, mean = 0, sigma = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if (!is.Numeric(d, allowable.length = 1) ||
      max(d) >= 2)
    stop("bad input for argument 'd'")

  L = max(length(x), length(mean), length(sigma))
  x = rep(x, len = L); mean = rep(mean, len = L);
  sigma = rep(sigma, len = L);
  hh = 2 - d
  KK = 1 / (1 + 1/hh + 0.75/hh^2)
  if (log.arg) {
    dnorm(x = x, mean = mean, sd = sigma, log = TRUE) + log(KK) +
    2 * log1p(((x-mean)/sigma)^2 / (2*hh))
  } else {
    dnorm(x = x, mean = mean, sd = sigma) * KK *
    (1 + ((x-mean)/sigma)^2 / (2*hh))^2
  }
}


ptikuv <- function(q, d, mean = 0, sigma = 1) {
  if (!is.Numeric(d, allowable.length = 1) ||
      max(d) >= 2)
    stop("bad input for argument 'd'")

  L = max(length(q), length(mean), length(sigma))
  q = rep(q, len = L); mean = rep(mean, len = L);
  sigma = rep(sigma, len = L);
  zedd1 = 0.5 * ((q - mean) / sigma)^2
  ans = q*0 + 0.5
  hh = 2 - d
  KK = 1 / (1 + 1/hh + 0.75/hh^2)
  if (any(lhs <- q < mean)) {
    ans[lhs] = ( KK/(2*sqrt(pi))) * (
    gamma(0.5) * (1 - pgamma(zedd1[lhs], 0.5)) +
    2 * gamma(1.5) * (1 - pgamma(zedd1[lhs], 1.5)) / hh +
    gamma(2.5) * (1 - pgamma(zedd1[lhs], 2.5)) / hh^2)
  }
  if (any(rhs <- q > mean)) {
    ans[rhs] = 1.0 - Recall(q = (2*mean[rhs] - q[rhs]), d = d,
               mean = mean[rhs], sigma = sigma[rhs])
  }
  ans
}


qtikuv <- function(p, d, mean = 0, sigma = 1, ...) {
  if (!is.Numeric(p, positive = TRUE) || max(p) >= 1)
    stop("bad input for argument 'p'")
  if (!is.Numeric(d, allowable.length = 1) || max(d) >= 2)
    stop("bad input for argument 'd'")
  if (!is.Numeric(mean))
    stop("bad input for argument 'mean'")
  if (!is.Numeric(sigma))
    stop("bad input for argument 'sigma'")
  L = max(length(p), length(mean), length(sigma))
  p = rep(p, len = L);
  mean = rep(mean, len = L);
  sigma = rep(sigma, len = L);
  ans = rep(0.0, len = L)

  myfun <- function(x, d, mean = 0, sigma = 1, p)
    ptikuv(q = x, d = d, mean = mean, sigma = sigma) - p
  for(i in 1:L) {
    Lower = ifelse(p[i] <= 0.5, mean[i] - 3 * sigma[i], mean[i])
    while (ptikuv(q = Lower, d = d, mean = mean[i],
                  sigma = sigma[i]) > p[i])
      Lower = Lower - sigma[i]
    Upper = ifelse(p[i] >= 0.5, mean[i] + 3 * sigma[i], mean[i])
    while (ptikuv(q = Upper, d = d, mean = mean[i],
                  sigma = sigma[i]) < p[i])
      Upper = Upper + sigma[i]
      ans[i] = uniroot(f = myfun, lower = Lower, upper = Upper,
                       d = d, p = p[i],
                       mean = mean[i], sigma = sigma[i], ...)$root
  }
  ans
}


rtikuv <- function(n, d, mean = 0, sigma = 1, Smallno = 1.0e-6) {
  if (!is.Numeric(n, positive = TRUE, integer.valued = TRUE))
    stop("bad input for argument 'n'")
  if (!is.Numeric(d, allowable.length = 1) || max(d) >= 2)
    stop("bad input for argument 'd'")
  if (!is.Numeric(mean, allowable.length = 1))
    stop("bad input for argument 'mean'")
  if (!is.Numeric(sigma, allowable.length = 1))
    stop("bad input for argument 'sigma'")
  if (!is.Numeric(Smallno, positive = TRUE, allowable.length = 1) ||
      Smallno > 0.01 ||
      Smallno < 2 * .Machine$double.eps)
      stop("bad input for argument 'Smallno'")
  ans = rep(0.0, len = n)

  ptr1 = 1; ptr2 = 0
  hh = 2 - d
  KK = 1 / (1 + 1/hh + 0.75/hh^2)
  ymax = ifelse(hh < 2,
                dtikuv(x = mean + sigma*sqrt(4 - 2*hh),
                       d = d, mean = mean, sigma = sigma),
                KK / (sqrt(2 * pi) * sigma))
  while (ptr2 < n) {
    Lower = mean - 5 * sigma
    while (ptikuv(q = Lower, d = d, mean = mean, sigma = sigma) > Smallno)
      Lower = Lower - sigma
    Upper = mean + 5 * sigma
    while (ptikuv(q = Upper, d = d, mean = mean, sigma = sigma) < 1-Smallno)
      Upper = Upper + sigma
    x = runif(2*n, min = Lower, max = Upper)
    index = runif(2*n, max = ymax) <
            dtikuv(x, d = d, mean = mean, sigma = sigma)
    sindex = sum(index)
    if (sindex) {
      ptr2 = min(n, ptr1 + sindex - 1)
      ans[ptr1:ptr2] = (x[index])[1:(1+ptr2-ptr1)]
      ptr1 = ptr2 + 1
    }
  }
  ans
}




 tikuv <- function(d, lmean = "identity", lsigma = "loge",
                  isigma = NULL, zero = 2)
{


  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsigma <- as.list(substitute(lsigma))
  e.sigma <- link2list(lsigma)
  l.sigma <- attr(e.sigma, "function.name")



  if (length(zero) &&
     (!is.Numeric(zero, integer.valued = TRUE, positive = TRUE) ||
     max(zero) > 2))
    stop("bad input for argument 'zero'")
  if (!is.Numeric(d, allowable.length = 1) || max(d) >= 2)
      stop("bad input for argument 'd'")



  new("vglmff",
  blurb = c("Short-tailed symmetric [Tiku and Vaughan (1999)] ",
            "distribution\n",
          "Link:     ",
          namesof("mean",  lmean,  earg = emean), ", ",
          namesof("sigma", l.sigma, earg = e.sigma),
          "\n", "\n",
          "Mean:     mean"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero)
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y)


    predictors.names <- 
      c(namesof("mean",  .lmean,  earg = .emean,  tag = FALSE),
        namesof("sigma", .l.sigma, earg = .e.sigma, tag = FALSE))


    if (!length(etastart)) {
      sigma.init = if (length(.i.sigma)) rep(.i.sigma, length = n) else {
        hh = 2 - .d
        KK = 1 / (1 + 1/hh + 0.75/hh^2)
        K2 = 1 + 3/hh + 15/(4*hh^2)
        rep(sqrt(var(y) / (KK*K2)), len = n)
      }
      mean.init = rep(weighted.mean(y, w), len = n) 
      etastart = cbind(theta2eta(mean.init,  .lmean,  earg = .emean),
                       theta2eta(sigma.init, .l.sigma, earg = .e.sigma))
    }
  }),list( .lmean = lmean, .l.sigma = l.sigma,
                             .i.sigma = isigma, .d = d,
           .emean = emean, .e.sigma = e.sigma ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .lmean, earg = .emean)
  }, list( .lmean = lmean,
           .emean = emean, .e.sigma = e.sigma ))),
  last = eval(substitute(expression({
      misc$link =    c("mean"= .lmean , "sigma"= .l.sigma )
      misc$earg = list("mean"= .emean , "sigma"= .e.sigma )
      misc$expected = TRUE
      misc$d = .d 
  }), list( .lmean = lmean, .l.sigma = l.sigma, .d = d,
            .emean = emean, .e.sigma = e.sigma ))),
  loglikelihood = eval(substitute(
      function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      mymu  = eta2theta(eta[, 1], .lmean,  earg = .emean)
      sigma = eta2theta(eta[, 2], .l.sigma, earg = .e.sigma)
      if (residuals) stop("loglikelihood residuals ",
                          "not implemented yet") else {
          sum(c(w) * dtikuv(x=y, d = .d , mean = mymu,
                         sigma = sigma, log = TRUE))
      }
  }, list( .lmean = lmean, .l.sigma = l.sigma, .d = d,
           .emean = emean, .e.sigma = e.sigma ))),
  vfamily = c("tikuv"),
  deriv = eval(substitute(expression({
    mymu  = eta2theta(eta[, 1], .lmean,  earg = .emean)
    sigma = eta2theta(eta[, 2], .l.sigma, earg = .e.sigma)

    dmu.deta = dtheta.deta(mymu, .lmean, earg = .emean)
    dsigma.deta = dtheta.deta(sigma, .l.sigma, earg = .e.sigma)

    zedd = (y - mymu) / sigma
    hh = 2 - .d 
    gzedd = zedd / (1 + 0.5*zedd^2 / hh)

    dl.dmu = zedd / sigma - 2 * gzedd / (hh*sigma)
    dl.dsigma = (zedd^2 - 1 - 2 * zedd * gzedd / hh) / sigma

    c(w) * cbind(dl.dmu    * dmu.deta,
                 dl.dsigma * dsigma.deta)
  }), list( .lmean = lmean, .l.sigma = l.sigma, .d = d,
            .emean = emean, .e.sigma = e.sigma ))),
  weight = eval(substitute(expression({
    ayy = 1 / (2*hh)
    Dnos = 1 - (2/hh) * (1 - ayy) / (1 + 2*ayy + 3*ayy^2)
    Dstar = -1 + 3 * (1 + 2*ayy + 11*ayy^2) / (1 + 2*ayy + 3*ayy^2)

    ned2l.dmymu2 = Dnos / sigma^2
    ned2l.dnu2   = Dstar / sigma^2

    wz = matrix(as.numeric(NA), n, M) # diagonal matrix
    wz[, iam(1, 1, M)] = ned2l.dmymu2 * dmu.deta^2
    wz[, iam(2, 2, M)] = ned2l.dnu2 * dsigma.deta^2
    c(w) * wz
  }), list( .lmean = lmean, .l.sigma = l.sigma,
            .emean = emean, .e.sigma = e.sigma ))))
}



dfnorm <- function(x, mean = 0, sd = 1, a1 = 1, a2=1) {
  if (!is.Numeric(a1, positive = TRUE) ||
      !is.Numeric(a2, positive = TRUE))
    stop("bad input for arguments 'a1' and 'a2'")
  if (any(a1 <= 0 | a2 <= 0))
    stop("arguments 'a1' and 'a2' must have positive values only")
  ans = dnorm(x = x/(a1*sd) - mean/sd)/(a1*sd) +
        dnorm(x = x/(a2*sd) + mean/sd)/(a2*sd)
  ans[x < 0] = 0
  ans[a1 <= 0 | a2 <= 0 | is.na(a1) | is.na(a2)] = NA
  ans
}


pfnorm <- function(q, mean = 0, sd = 1, a1 = 1, a2=1) {
  if (!is.Numeric(a1, positive = TRUE) ||
      !is.Numeric(a2, positive = TRUE))
    stop("bad input for arguments 'a1' and 'a2'")
  if (any(a1 <= 0 | a2 <= 0))
    stop("arguments 'a1' and 'a2' must have positive values only")
  L = max(length(q), length(mean), length(sd))
  q = rep(q, len = L);
  mean = rep(mean, len = L);
  sd = rep(sd, len = L);

  ifelse(q < 0, 0,
         pnorm(q =  q/(a1*sd) - mean/sd) -
         pnorm(q = -q/(a2*sd) - mean/sd))
}


qfnorm <- function(p, mean = 0, sd = 1, a1 = 1, a2 = 1, ...) {
  if (!is.Numeric(p, positive = TRUE) || max(p) >= 1)
    stop("bad input for argument 'p'")
  if (!is.Numeric(a1, positive = TRUE) ||
     !is.Numeric(a2, positive = TRUE))
    stop("bad input for arguments 'a1' and 'a2'")
  if (any(a1 <= 0 | a2 <= 0))
    stop("arguments 'a1' and 'a2' must have positive values only")

  L = max(length(p), length(mean), length(sd), length(a1), length(a2))
  p = rep(p, len = L);
  mean = rep(mean, len = L);
  sd = rep(sd, len = L);
  a1 = rep(a1, len = L);
  a2 = rep(a2, len = L);
  ans = rep(0.0, len = L)

  myfun <- function(x, mean = 0, sd = 1, a1 = 1, a2=2, p)
    pfnorm(q = x, mean = mean, sd = sd, a1 = a1, a2 = a2) - p
  for(i in 1:L) {
    mytheta = mean[i]/sd[i]
    EY = sd[i] * ((a1[i]+a2[i]) *
         (mytheta * pnorm(mytheta) + dnorm(mytheta)) -
         a2[i] * mytheta)
    Upper = 2 * EY
    while (pfnorm(q = Upper, mean = mean[i], sd = sd[i],
                  a1 = a1[i], a2 = a2[i]) < p[i])
        Upper = Upper + sd[i]
    ans[i] = uniroot(f = myfun, lower = 0, upper = Upper,
                     mean = mean[i],
                     sd = sd[i], a1 = a1[i], a2 = a2[i],
                     p = p[i], ...)$root
  }
  ans
}


rfnorm <- function(n, mean = 0, sd = 1, a1 = 1, a2=1) {
  if (!is.Numeric(n, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'n'")
  if (!is.Numeric(a1, positive = TRUE) ||
      !is.Numeric(a2, positive = TRUE))
    stop("bad input for arguments 'a1' and 'a2'")
  if (any(a1 <= 0 | a2 <= 0))
    stop("arguments 'a1' and 'a2' must have positive values only")
  X = rnorm(n, mean = mean, sd = sd)
  pmax(a1 * X, -a2*X)
}




 fnormal1 <- function(lmean = "identity", lsd = "loge",
                      imean = NULL,       isd = NULL,
                      a1 = 1, a2 = 1,
                      nsimEIM = 500, imethod = 1, zero = NULL)
{
  if (!is.Numeric(a1, positive = TRUE, allowable.length = 1) ||
      !is.Numeric(a2, positive = TRUE, allowable.length = 1))
    stop("bad input for arguments 'a1' and 'a2'")
  if (any(a1 <= 0 | a2 <= 0))
    stop("arguments 'a1' and 'a2' must each be a positive value")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")




  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")





  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")

  if (!is.Numeric(nsimEIM, allowable.length = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 10)
    stop("argument 'nsimEIM' should be an integer greater than 10")
  if (length(imean) && !is.Numeric(imean))
    stop("bad input for 'imean'")

  if (length(isd) && !is.Numeric(isd, positive = TRUE))
    stop("bad input for 'isd'")


  new("vglmff",
  blurb = c("(Generalized) folded univariate normal distribution\n\n",
          "Link:     ",
          namesof("mean", lmean, earg = emean, tag = TRUE), "; ",
          namesof("sd",   lsd,   earg = esd,   tag = TRUE)),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
        c(namesof("mean", .lmean, earg = .emean, tag = FALSE),
          namesof("sd",   .lsd,   earg = .esd,   tag = FALSE))

    if (!length(etastart)) {
        junk = lm.wfit(x = x, y=y, w = w)


 if (FALSE) {
        if ((ncol(cbind(w)) != 1) || any(w != round(w)))
            stop("'weights' must be a vector or a one-column matrix ",
                 "with integer values")
            m1d = meany = weighted.mean(y, w)
            m2d = weighted.mean(y^2, w)
            stddev = sqrt( sum(c(w) * junk$resid^2) / junk$df.residual )
            Ahat = m1d^2 / m2d
            thetahat = sqrt(max(1/Ahat -1, 0.1))
            mean.init = rep(if(length( .i.mean)) .i.mean else
                thetahat * sqrt((stddev^2 + meany^2) * Ahat), len = n)
            sd.init = rep(if(length( .i.sd)) .i.sd else
                sqrt((stddev^2 + meany^2) * Ahat), len = n)
}


        stddev = sqrt( sum(c(w) * junk$resid^2) / junk$df.residual )
        meany = weighted.mean(y, w)
        mean.init = rep(if(length( .i.mean)) .i.mean else
            {if( .imethod == 1) median(y) else meany}, len = n)
        sd.init = rep(if(length( .i.sd)) .i.sd else
            {if( .imethod == 1)  stddev else 1.2*sd(c(y))}, len = n)
        etastart = cbind(theta2eta(mean.init, .lmean, earg = .emean),
                         theta2eta(sd.init,   .lsd,   earg = .esd ))
    }
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd,
            .i.mean = imean, .i.sd = isd,
            .a1 = a1, .a2 = a2, .imethod = imethod ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        mymu = eta2theta(eta[, 1], .lmean, earg = .emean)
        mysd = eta2theta(eta[, 2], .lsd, earg = .esd )
        mytheta = mymu/mysd
        mysd * (( .a1+ .a2) * (mytheta * pnorm(mytheta) +
                dnorm(mytheta)) - .a2 * mytheta)
  }, list( .lmean = lmean, .lsd = lsd,
           .emean = emean, .esd = esd,
           .a1 = a1, .a2 = a2 ))),
  last = eval(substitute(expression({
    misc$link =    c("mu" = .lmean , "sd" = .lsd )

    misc$earg = list("mu" = .emean , "sd" = .esd )

    misc$multipleResponses <- FALSE
    misc$expected = TRUE
    misc$nsimEIM = .nsimEIM
    misc$simEIM = TRUE
    misc$imethod = .imethod
    misc$a1 = .a1
    misc$a2 = .a2
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd,
            .imethod = imethod, .nsimEIM = nsimEIM,
            .a1 = a1, .a2 = a2 ))),
  loglikelihood = eval(substitute(
      function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      mymu = eta2theta(eta[, 1], .lmean, earg = .emean)
      mysd = eta2theta(eta[, 2], .lsd, earg = .esd )
      a1vec = .a1
      a2vec = .a2
      if (residuals) stop("loglikelihood residuals ",
                          "not implemented yet") else {
            sum(c(w)*log(dnorm(x=y/(a1vec*mysd) - mymu/mysd)/(a1vec*mysd) +
                      dnorm(x=y/(a2vec*mysd) + mymu/mysd)/(a2vec*mysd)))
        }
    }, list( .lmean = lmean, .lsd = lsd,
             .emean = emean, .esd = esd, .a1 = a1, .a2 = a2 ))),
    vfamily = c("fnormal1"),
    deriv = eval(substitute(expression({
        mymu = eta2theta(eta[, 1], .lmean, earg = .emean)
        mysd = eta2theta(eta[, 2], .lsd, earg = .esd )
        dmu.deta = dtheta.deta(mymu, .lmean, earg = .emean)
        dsd.deta = dtheta.deta(mysd, .lsd, earg = .esd )
        a1vec = .a1
        a2vec = .a2
        d3 = deriv3(~ log((exp(-0.5*(y/(a1vec*mysd) - mymu/mysd)^2)/a1vec +
                           exp(-0.5*(y/(a2vec*mysd) +
                               mymu/mysd)^2)/a2vec)/(mysd*sqrt(2*pi))),
                    name=c("mymu","mysd"), hessian= FALSE)
        eval.d3 = eval(d3)
        dl.dthetas =  attr(eval.d3, "gradient")  # == cbind(dl.dmu, dl.dsd)
        DTHETA.detas = cbind(dmu.deta, dsd.deta)
        c(w) * DTHETA.detas * dl.dthetas
    }), list( .lmean = lmean, .lsd = lsd, .emean = emean, .esd = esd,
              .a1 = a1, .a2 = a2 ))),
    weight = eval(substitute(expression({
        de3 = deriv3(~ log((exp(-0.5*(ysim/(a1vec*mysd) -
                                 mymu/mysd)^2)/a1vec +
                            exp(-0.5*(ysim/(a2vec*mysd) +
                                 mymu/mysd)^2)/a2vec)/(mysd*sqrt(2*pi))),
                     name=c("mymu","mysd"), hessian= TRUE)
        run.mean = 0
        for(ii in 1:( .nsimEIM )) {
            ysim = rfnorm(n = n, mean = mymu, sd = mysd,
                          a1 = a1vec, a2 = a2vec)
            eval.de3 = eval(de3)
            d2l.dthetas2 =  attr(eval.de3, "hessian")
            rm(ysim)

            temp3 = matrix(0, n, dimm(M))
            for(ss in 1:M)
                for(tt in ss:M)
                    temp3[, iam(ss,tt, M)] =  -d2l.dthetas2[, ss,tt]

            run.mean = ((ii-1) * run.mean + temp3) / ii
        }

        wz = if (intercept.only)
            matrix(colMeans(run.mean), n, dimm(M), byrow = TRUE) else
            run.mean

        index0 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)
        wz = wz * DTHETA.detas[, index0$row] * DTHETA.detas[, index0$col]
        c(w) * wz
    }), list( .nsimEIM = nsimEIM, .a1 = a1, .a2 = a2 ))))
}





lqnorm.control <- function(trace = TRUE, ...) {
    list(trace = trace)
}





lqnorm <- function(qpower = 2,
                   link = "identity",
                   imethod = 1, imu = NULL, shrinkage.init = 0.95)
{


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")



  if (!is.Numeric(qpower, allowable.length = 1) || qpower <= 1)
    stop("bad input for argument 'qpower'")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")

  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
      shrinkage.init < 0 ||
      shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")



    new("vglmff",
    blurb = c("Minimizing the q-norm of residuals\n",
            "Links:    ",
            namesof("Y1", link, earg = earg, tag = TRUE)),
    initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    M = if (is.matrix(y)) ncol(y) else 1
    dy = dimnames(y)


    predictors.names <- if (!is.null(dy[[2]])) dy[[2]] else
                       paste("mu", 1:M, sep = "")
    predictors.names <- namesof(predictors.names, link = .link,
                               earg = .earg, short = TRUE)


    if (!length(etastart))  {
        meany = weighted.mean(y, w)
        mean.init = rep(if(length( .i.mu)) .i.mu else
            {if( .imethod == 2) median(y) else 
             if ( .imethod == 1) meany else
             .sinit * meany + (1 - .sinit) * y
            }, len = n)
        etastart = theta2eta(mean.init, link = .link, earg = .earg)
    }
  }), list( .imethod = imethod, .i.mu = imu,
            .sinit = shrinkage.init,
            .link = link, .earg = earg ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
      mu = eta2theta(eta, link = .link, earg = .earg)
      mu
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    dy = dimnames(y)
    if (!is.null(dy[[2]]))
        dimnames(fit$fitted.values) = dy
    misc$link = rep( .link, length = M)
    names(misc$link) = predictors.names

    misc$earg = list(mu = .earg)

    misc$qpower = .qpower
    misc$imethod = .imethod
    misc$objectiveFunction = sum( c(w) * (abs(y - mu))^(.qpower) )
  }), list( .qpower = qpower,
            .link = link, .earg = earg,
            .imethod = imethod ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, link = .link, earg = .earg)
  }, list( .link = link, .earg = earg ))),
  vfamily = "lqnorm",
  deriv = eval(substitute(expression({
    dmu.deta = dtheta.deta(theta=mu, link = .link, earg = .earg )
    myresid = y - mu
    signresid = sign(myresid)
    temp2 = (abs(myresid))^(.qpower-1)
    .qpower * c(w) * temp2 * signresid * dmu.deta
  }), list( .qpower = qpower, .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    temp3 = (abs(myresid))^(.qpower-2)
    wz = .qpower * (.qpower - 1) * c(w) * temp3 * dmu.deta^2
    wz
  }), list( .qpower = qpower, .link = link, .earg = earg ))))
}







dtobit <- function(x, mean = 0, sd = 1,
                   Lower = 0, Upper = Inf, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  L = max(length(x), length(mean), length(sd), length(Lower),
          length(Upper))
  x = rep(x, len = L);
  mean = rep(mean, len = L);
  sd = rep(sd, len = L);
  Lower = rep(Lower, len = L);
  Upper = rep(Upper, len = L);

  ans = dnorm(x = x, mean = mean, sd = sd, log = log.arg)
  ans[x <  Lower] = if (log.arg) log(0.0) else 0.0
  ans[x >  Upper] = if (log.arg) log(0.0) else 0.0


  ind3 <- x == Lower
  ans[ind3] = if (log.arg) {
                log(exp(ans[ind3]) +
                    pnorm(q = Lower[ind3], mean = mean[ind3], sd = sd[ind3]))
              } else {
                ans[ind3] +
                pnorm(q = Lower[ind3], mean = mean[ind3], sd = sd[ind3])
              }

  ind4 <- x == Upper
  ans[ind4] = if (log.arg) {
                log(exp(ans[ind4]) +
                    pnorm(q = Upper[ind4], mean = mean[ind4], sd = sd[ind4],
                          lower.tail = FALSE))
              } else {
                ans[ind4] +
                pnorm(q = Upper[ind4], mean = mean[ind4], sd = sd[ind4],
                      lower.tail = FALSE)
              }
  ans
}



ptobit <- function(q, mean = 0, sd = 1,
                  Lower = 0, Upper = Inf,
                  lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail) != 1)
    stop("argument 'lower.tail' must be a single logical")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("argument 'log.p' must be a single logical")

  L = max(length(q), length(mean), length(sd), length(Lower),
          length(Upper))
  q = rep(q, len = L);
  mean = rep(mean, len = L);
  sd = rep(sd, len = L);
  Lower = rep(Lower, len = L);
  Upper = rep(Upper, len = L);

  ans = pnorm(q = q, mean = mean, sd = sd, lower.tail = lower.tail)
  ind1 <- (q <  Lower)
  ans[ind1] = if (lower.tail) ifelse(log.p, log(0.0), 0.0) else
                              ifelse(log.p, log(1.0), 1.0)
  ind2 <- (Upper <= q)
  ans[ind2] = if (lower.tail) ifelse(log.p, log(1.0), 1.0) else
                              ifelse(log.p, log(0.0), 0.0)

  ans
}




qtobit <- function(p, mean = 0, sd = 1,
                  Lower = 0, Upper = Inf) {

  L = max(length(p), length(mean), length(sd), length(Lower),
          length(Upper))
  p = rep(p, len = L);
  mean = rep(mean, len = L);
  sd = rep(sd, len = L);
  Lower = rep(Lower, len = L);
  Upper = rep(Upper, len = L);

  ans = qnorm(p = p, mean = mean, sd = sd)
  pnorm.Lower = ptobit(q = Lower, mean = mean, sd = sd)
  pnorm.Upper = ptobit(q = Upper, mean = mean, sd = sd)

  ind1 <- (p <= pnorm.Lower)
  ans[ind1] = Lower[ind1]

  ind2 <- (pnorm.Upper <= p)
  ans[ind2] = Upper[ind2]

  ans
}






rtobit <- function(n, mean = 0, sd = 1,
                  Lower = 0, Upper = Inf) {

  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE,
                          allowable.length = 1, positive = TRUE))
            stop("bad input for argument 'n'") else n
  L = max(use.n, length(mean), length(sd), length(Lower),
          length(Upper))
  mean = rep(mean, len = L);
  sd = rep(sd, len = L);
  Lower = rep(Lower, len = L);
  Upper = rep(Upper, len = L);

  ans = rnorm(n = use.n, mean = mean, sd = sd)
  cenL <- (ans < Lower)
  ans[cenL] = Lower[cenL]
  cenU <- (ans > Upper)
  ans[cenU] = Upper[cenU]

  attr(ans, "Lower") <- Lower
  attr(ans, "Upper") <- Upper
  attr(ans, "cenL") <- cenL
  attr(ans, "cenU") <- cenU
  ans
}




tobit.control <- function(save.weight = TRUE, ...)
{
  list(save.weight = save.weight)
}


 tobit <- function(Lower = 0, Upper = Inf,
                   lmu = "identity",  lsd = "loge",
                   nsimEIM = 250,
                   imu = NULL,        isd = NULL,
                   type.fitted = c("uncensored", "censored", "mean.obs"),
                   imethod = 1, zero = -2) {








  lmu <- as.list(substitute(lmu))
  e.mu <- link2list(lmu)
  l.mu <- attr(e.mu, "function.name")

  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")



  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
    imethod > 2)
    stop("argument 'imethod' must be 1 or 2")
  if ( # length(Lower) != 1 || length(Upper) != 1 ||
    !is.numeric(Lower) ||
    !is.numeric(Upper) ||
    any(Lower >= Upper))
    stop("Lower and Upper must ",
         "be numeric with Lower < Upper")

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE))
    stop("bad input for argument 'zero'")
  if (!is.Numeric(nsimEIM, allowable.length = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 10)
    stop("argument 'nsimEIM' should be an integer greater than 10")

  if(mode(type.fitted) != "character" && mode(type.fitted) != "name")
        type.fitted <- as.character(substitute(type.fitted))
  type.fitted <- match.arg(type.fitted,
                           c("uncensored", "censored", "mean.obs"))[1]


  stdTobit = all(Lower == 0.0) &&
             all(!is.finite(Upper)) &&
             all(lmu == "identity")


  new("vglmff",
  blurb = c("Tobit model\n\n",
          "Links:    ",
          namesof("mu", l.mu, earg = e.mu, tag = TRUE), "; ",
          namesof("sd", lsd, earg = esd, tag = TRUE), "\n",
          "Mean:                 mu", "\n",
          "Conditional variance: sd^2"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)

  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero ,
         nsimEIM = .nsimEIM )
  }, list( .zero = zero, .nsimEIM = nsimEIM ))),

  initialize = eval(substitute(expression({
    Musual = 2


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    ncoly = ncol(y)
    M = Musual * ncoly

    Lowmat = matrix( .Lower , nrow = n, ncol = ncoly, byrow = TRUE)
    Uppmat = matrix( .Upper , nrow = n, ncol = ncoly, byrow = TRUE)

    extra$censoredL = (y <= Lowmat)
    extra$censoredU = (y >= Uppmat)
    if (any(y < Lowmat)) {
      warning("replacing response values less than the value ",
              .Lower , " by ", .Lower )
      y[y < Lowmat] = Lowmat[y < Lowmat]
    }
    if (any(y > Uppmat)) {
      warning("replacing response values greater than the value ",
              .Upper, " by ", .Upper)
      y[y > Uppmat] = Uppmat[y > Uppmat]
    }

    temp1.names <-
      if (ncoly == 1) "mu" else paste("mu", 1:ncoly, sep = "")
    temp2.names <-
      if (ncoly == 1) "sd" else paste("sd", 1:ncoly, sep = "")
    predictors.names <-
        c(namesof(temp1.names, .l.mu, earg = .e.mu, tag = FALSE),
          namesof(temp2.names, .lsd, earg = .esd, tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M = Musual)]

    if (!length(etastart)) {
      anyc <- cbind(extra$censoredL | extra$censoredU)
      i11 <- if ( .imethod == 1) anyc else FALSE # can be all data

      mu.init <-
      sd.init <- matrix(0.0, n, ncoly)
      for(ii in 1:ncol(y)) {
        use.i11 <- i11[, ii]
        mylm <- lm.wfit(x = cbind(x[!use.i11, ]),
                       y = y[!use.i11, ii], w = w[!use.i11, ii])
        sd.init[, ii] <- sqrt( sum(w[!use.i11, ii] * mylm$resid^2)
                              / mylm$df.residual ) * 1.5
        mu.init[!use.i11, ii] <- mylm$fitted.values
        if (any(anyc[, ii]))
          mu.init[anyc[, ii], ii] <- x[anyc[, ii],, drop = FALSE] %*%
                                    mylm$coeff
      }

      if (length( .i.mu ))
        mu.init <- matrix( .i.mu , n, ncoly, byrow = TRUE)
      if (length( .i.sd ))
        sd.init <- matrix( .i.sd , n, ncoly, byrow = TRUE)

      etastart <- cbind(theta2eta(mu.init, .l.mu, earg = .e.mu ),
                       theta2eta(sd.init, .lsd, earg = .esd ))

      etastart <- etastart[, interleave.VGAM(M, M = Musual), drop = FALSE]
    }
 }), list( .Lower = Lower, .Upper = Upper,
           .l.mu = l.mu, .lsd = lsd,
           .e.mu = e.mu, .esd = esd,
           .i.mu = imu, .i.sd = isd,
           .imethod = imethod ))),
  linkinv = eval(substitute( function(eta, extra = NULL) {
    Musual = 2
    ncoly = ncol(eta) / Musual
    mum = eta2theta(eta[, Musual*(1:ncoly)-1, drop=FALSE], .l.mu, earg = .e.mu )
    if ( .type.fitted == "uncensored")
      return(mum)

    Lowmat = matrix( .Lower, nrow = nrow(eta), ncol = ncoly, byrow = TRUE)
    Uppmat = matrix( .Upper, nrow = nrow(eta), ncol = ncoly, byrow = TRUE)
    if ( .type.fitted == "censored") {
      mum[mum < Lowmat] <- Lowmat[mum < Lowmat]
      mum[mum > Uppmat] <- Uppmat[mum > Uppmat]
      mum
    } else {

      sdm = eta2theta(eta[, Musual*(1:ncoly)-0, drop = FALSE],
                      .lsd , earg = .esd )
      zeddL = (Lowmat - mum) / sdm
      zeddU = (Uppmat - mum) / sdm
      Phi.L = pnorm(zeddL)
      phi.L = dnorm(zeddL)
      Phi.U = pnorm(zeddU)
      phi.U = dnorm(zeddU)
      mum * (Phi.U - Phi.L) +
      sdm * (phi.L - phi.U) +
      Lowmat *      Phi.L +
      Uppmat * (1 - Phi.U)
    }
  }, list( .l.mu = l.mu, .lsd = lsd,
           .e.mu = e.mu, .esd = esd,
           .Lower = Lower, .Upper = Upper,
           .type.fitted = type.fitted ))),
  last = eval(substitute(expression({

    temp0303 = c(rep( .l.mu, length = ncoly),
                 rep( .lsd, length = ncoly))
    names(temp0303) =
      c(if (ncoly == 1) "mu" else paste("mu", 1:ncoly, sep = ""),
        if (ncoly == 1) "sd" else paste("sd", 1:ncoly, sep = ""))
    temp0303 = temp0303[interleave.VGAM(M, M = Musual)]
    misc$link = temp0303 # Already named

    misc$earg = vector("list", M)
    names(misc$earg) = names(misc$link)
    for(ii in 1:ncoly) {
      misc$earg[[Musual*ii-1]] = .e.mu
      misc$earg[[Musual*ii  ]] = .esd
    }

    misc$multipleResponses <- TRUE
    misc$expected = TRUE
    misc$imethod = .imethod
    misc$nsimEIM = .nsimEIM
    misc$Musual = Musual
    misc$stdTobit = .stdTobit
    misc$Lower = Lowmat
    misc$Upper = Uppmat


    if ( .stdTobit ) {
      save.weight <- control$save.weight <- FALSE
      fit$weights <- NULL
    }


  }), list( .l.mu = l.mu, .lsd = lsd,
            .e.mu = e.mu, .esd = esd,
            .nsimEIM = nsimEIM, .imethod = imethod,
            .stdTobit = stdTobit,
            .Lower = Lower,
            .Upper = Upper ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Musual = 2
    y = cbind(y)
    ncoly = ncol(y)

    cenL = extra$censoredL
    cenU = extra$censoredU
    cen0 = !cenL & !cenU   # uncensored obsns
    Lowmat = matrix( .Lower , nrow = nrow(eta), ncol = ncoly, byrow = TRUE)
    Uppmat = matrix( .Upper , nrow = nrow(eta), ncol = ncoly, byrow = TRUE)


    mum = eta2theta(eta[, Musual*(1:ncoly)-1, drop = FALSE],
                    .l.mu , earg = .e.mu )
    sdm = eta2theta(eta[, Musual*(1:ncoly)-0, drop = FALSE],
                    .lsd , earg = .esd )

    ell0 = dnorm(  y[cen0], mean = mum[cen0], sd = sdm[cen0],
                 log = TRUE)
    ellL = pnorm(Lowmat[cenL], mean = mum[cenL], sd = sdm[cenL],
                 log.p = TRUE, lower.tail = TRUE)
    ellU = pnorm(Uppmat[cenU], mean = mum[cenU], sd = sdm[cenU],
                 log.p = TRUE, lower.tail = FALSE)

    wmat = matrix(w, nrow = nrow(eta), ncol = ncoly)
    if (residuals) {
      stop("loglikelihood residuals not ",
           "implemented yet") 
    } else {
      sum(wmat[cen0] * ell0) +
      sum(wmat[cenL] * ellL) +
      sum(wmat[cenU] * ellU)
    }
  }, list( .l.mu = l.mu, .lsd = lsd,
           .e.mu = e.mu, .esd = esd,
           .Lower = Lower, .Upper = Upper ))),
  vfamily = c("tobit"),
  deriv = eval(substitute(expression({
    Musual = 2
    y = cbind(y)
    ncoly = ncol(y)

    Lowmat = matrix( .Lower , nrow = n, ncol = ncoly, byrow = TRUE)
    Uppmat = matrix( .Upper , nrow = n, ncol = ncoly, byrow = TRUE)

    cenL = extra$censoredL
    cenU = extra$censoredU
    cen0 = !cenL & !cenU   # uncensored obsns

    mum = eta2theta(eta[, Musual*(1:ncoly)-1, drop = FALSE], .l.mu, earg = .e.mu )
    sdm = eta2theta(eta[, Musual*(1:ncoly)-0, drop = FALSE], .lsd, earg = .esd )

    zedd = (y - mum) / sdm
    dl.dmu = zedd / sdm
    dl.dsd = (zedd^2 - 1) / sdm

    dmu.deta = dtheta.deta(mum, .l.mu, earg = .e.mu )
    dsd.deta = dtheta.deta(sdm, .lsd, earg = .esd )

    if (any(cenL)) {
      mumL = Lowmat - mum
      temp21L = mumL[cenL] / sdm[cenL]
      PhiL = pnorm(temp21L)
      phiL = dnorm(temp21L)
      fred21 = phiL / PhiL
      dl.dmu[cenL] = -fred21 / sdm[cenL]
      dl.dsd[cenL] =  fred21 * (-mumL[cenL] / sdm[cenL]^2)
    }
    if (any(cenU)) {
      mumU = Uppmat - mum
      temp21U = mumU[cenU] / sdm[cenU]
      PhiU = pnorm(temp21U, lower.tail = FALSE)
      phiU = dnorm(temp21U)
      fred21 = -phiU / PhiU
      dl.dmu[cenU] = -fred21 / sdm[cenU]   # Negated
      dl.dsd[cenU] =  fred21 * (-mumU[cenU] / sdm[cenU]^2)
    }

    dthetas.detas = cbind(dmu.deta, dsd.deta)
    dThetas.detas = dthetas.detas[, interleave.VGAM(M, M = Musual)]

    myderiv = cbind(c(w) * dl.dmu,
                    c(w) * dl.dsd) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list( .l.mu = l.mu, .lsd = lsd,
            .e.mu = e.mu, .esd = esd,
            .Lower = Lower, .Upper = Upper ))),
  weight = eval(substitute(expression({

    wz = matrix(0.0, n, M + M - 1) # wz is 'tridiagonal'
    ind1 = iam(NA, NA, M = Musual, both = TRUE, diag = TRUE)


    if (is.numeric( .nsimEIM ) &&
        ! .stdTobit ) {


    run.varcov = 0

    for(spp. in 1:ncoly) {
      run.varcov = 0
      muvec = mum[, spp.]
      sdvec = sdm[, spp.]

      for(ii in 1:( .nsimEIM )) {
        ysim = rtobit(n = n, mean = muvec, sd = sdvec,
                      Lower = Lowmat[, spp.], Upper = Uppmat[, spp.])
        cenL = attr(ysim, "cenL")
        cenU = attr(ysim, "cenU")
        cen0 = !cenL & !cenU   # uncensored obsns

        zedd = (ysim - muvec) / sdvec
        dl.dmu =   zedd / sdvec
        dl.dsd = (zedd^2 - 1) / sdvec
      if (any(cenL)) {
        mumL = Lowmat[, spp.] - muvec
        temp21L = mumL[cenL] / sdvec[cenL]
        PhiL = pnorm(temp21L)
        phiL = dnorm(temp21L)
        fred21 = phiL / PhiL
        dl.dmu[cenL] = -fred21 / sdvec[cenL]
        dl.dsd[cenL] =  fred21 * (-mumL[cenL] / sdvec[cenL]^2)
      }
      if (any(cenU)) {
        mumU = Uppmat[, spp.] - muvec
        temp21U = mumU[cenU] / sdvec[cenU]
        PhiU = pnorm(temp21U, lower.tail = FALSE)
        phiU = dnorm(temp21U)
        fred21 = -phiU / PhiU
        dl.dmu[cenU] = -fred21 / sdvec[cenU]   # Negated
        dl.dsd[cenU] =  fred21 * (-mumU[cenU] / sdvec[cenU]^2)
      }

      rm(ysim)
      temp3 = cbind(dl.dmu, dl.dsd)
      run.varcov = run.varcov +
                   temp3[, ind1$row.index] *
                   temp3[, ind1$col.index]
    }
    run.varcov = run.varcov / .nsimEIM

    wz1 = if (intercept.only && FALSE)
        matrix(colMeans(run.varcov),
               n, ncol(run.varcov), byrow = TRUE) else
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

    } else {

      wz1 = matrix(0.0, n, dimm(Musual))
      for(spp. in 1:ncoly) {
        zedd  = (y[, spp.] - mum[, spp.]) / sdm[, spp.]
        zedd0 = (            mum[, spp.]) / sdm[, spp.]
        phivec = dnorm(zedd0)
        Phivec = pnorm(zedd0)

        wz1[, iam(1, 1, M = Musual)] =   -(phivec * zedd0 -
                                        phivec^2 / (1 - Phivec) -
                                        Phivec)
        wz1[, iam(2, 2, M = Musual)] =   -(phivec   * zedd0^3 +
                                        phivec   * zedd0 -
                                        phivec^2 * zedd0^2 / (1 - Phivec) -
                                        2 * Phivec)
        wz1[, iam(1, 2, M = Musual)] = +(phivec   * zedd0^2 +
                                        phivec   -
                                        phivec^2 * zedd0 / (1 - Phivec))

        wz1 = wz1 / sdm[, spp.]^2
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

    } # End of EIM


    temp = w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = ncoly)


    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = ncoly)
  }), list( .l.mu = lmu, .Lower = Lower, .Upper = Upper,
            .lsd = lsd,
            .stdTobit = stdTobit,
            .nsimEIM = nsimEIM ))))
} # End of tobit()





 normal1 <- function(lmean = "identity", lsd = "loge", lvar = "loge",
                     var.arg = FALSE,
                     imethod = 1,
                     isd = NULL,
                     parallel = FALSE,
                     intercept.apply = FALSE,
                     zero = -2)
{





  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")

  lvar <- as.list(substitute(lvar))
  e.var <- link2list(lvar)
  l.var <- attr(e.var, "function.name")





  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE))
      stop("bad input for argument 'zero'")


  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
      stop("argument 'imethod' must be 1 or 2 or 3 or 4")

  if (!is.logical(var.arg) || length(var.arg) != 1)
    stop("argument 'var.arg' must be a single logical")
  if (!is.logical(intercept.apply) ||
      length(intercept.apply) != 1)
    stop("argument 'intercept.apply' must be a single logical")


  if (is.logical(parallel) && parallel && length(zero))
    stop("set 'zero = NULL' if 'parallel = TRUE'")


  new("vglmff",
  blurb = c("Univariate normal distribution\n\n",
            "Links:    ",
            namesof("mean", lmean, earg = emean, tag = TRUE), "; ",
            if (var.arg)
            namesof("var",  l.var, earg = e.var, tag = TRUE) else
            namesof("sd" ,  lsd,  earg = esd,  tag = TRUE),
            "\n",
            if (var.arg) "Variance: var" else "Variance: sd^2"),



  constraints = eval(substitute(expression({

    constraints <- cm.vgam(matrix(1, M, 1), x, .parallel , constraints,
                          intercept.apply = .intercept.apply )

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero,
            .parallel = parallel, .intercept.apply = intercept.apply ))),


  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero)
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({
    orig.y <- y








    if (length(attr(orig.y, "Prior.Weights"))) {
      if (any(c(w) != 1))
        warning("replacing the 'weights' argument by the 'Prior.Weights'",
                "attribute of the response (probably due to Qvar()")


      w = attr(orig.y, "Prior.Weights")


      extra$attributes.y = attributes(orig.y)

    } else {
    }






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



    mynames1 <- paste("mean",
                      if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste(if ( .var.arg ) "var" else "sd",
                      if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lmean , earg = .emean , tag = FALSE),
          if ( .var.arg ) 
          namesof(mynames2, .l.var  , earg = .e.var  , tag = FALSE) else
          namesof(mynames2, .lsd   , earg = .esd   , tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M = Musual)]
    extra$predictors.names <- predictors.names


    if (!length(etastart)) {
      sdev.init <- mean.init <- matrix(0, n, ncoly)
      for (jay in 1:ncoly) {
        jfit <- lm.wfit(x = x,  y = y[, jay], w = w[, jay])
        mean.init[, jay] <- if ( .lmean == "loge")
                            pmax(1/1024, y[, jay]) else
          if( .imethod == 1) median(y[, jay]) else
          if( .imethod == 2) weighted.mean(y[, jay], w = w[, jay]) else
          if( .imethod == 3) weighted.mean(y[, jay], w = w[, jay]) *
                             0.5 + y[, jay] * 0.5 else
                                 mean(jfit$fitted)

        sdev.init[, jay] <-
          if( .imethod == 1) {
            sqrt( sum(w[, jay] *
                (y[, jay] - mean.init[, jay])^2) / sum(w[, jay]) )
          } else if( .imethod == 2) {
            if (jfit$df.resid > 0)
              sqrt( sum(w[, jay] * jfit$resid^2) / jfit$df.resid ) else
              sqrt( sum(w[, jay] * jfit$resid^2) / sum(w[, jay]) )
          } else if( .imethod == 3) {
            sqrt( sum(w[, jay] * 
                  (y[, jay] - mean.init[, jay])^2) / sum(w[, jay]) )
          } else {
            sqrt( sum(w[, jay] * abs(y[, jay] -
                                     mean.init[, jay])) / sum(w[, jay]) )
          }

        if (any(sdev.init[, jay] <= sqrt( .Machine$double.eps ) ))
          sdev.init[, jay] <- 1.01

      }


      if (length( .i.sd )) {
        sdev.init <- matrix( .i.sd , n, ncoly, byrow = TRUE)
      }


      etastart <-
        cbind(theta2eta(mean.init, .lmean , earg = .emean ),
              if ( .var.arg )
              theta2eta(sdev.init^2, .l.var , earg = .e.var ) else
              theta2eta(sdev.init  , .lsd  , earg = .esd  ))
      etastart <-
        etastart[, interleave.VGAM(ncol(etastart), M = Musual)]

      colnames(etastart) <- predictors.names
    }
  }), list( .lmean = lmean, .lsd = lsd, .l.var = l.var,
            .emean = emean, .esd = esd, .e.var = e.var,
                              .i.sd = isd,
            .var.arg = var.arg, .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Musual <- extra$Musual
    ncoly <- extra$ncoly
    eta2theta(eta[, Musual*(1:ncoly) - 1], .lmean , earg = .emean )
  }, list( .lmean = lmean,
           .emean = emean, .esd = esd , .e.var = e.var ))),

  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <- c(rep( .lmean , length = ncoly),
                   rep( .lsd   , length = ncoly))
    misc$link <- misc$link[interleave.VGAM(Musual * ncoly, M = Musual)]
    temp.names <- c(mynames1, mynames2)
    temp.names <- temp.names[interleave.VGAM(Musual * ncoly, M = Musual)]
    names(misc$link) <- temp.names


    misc$earg <- vector("list", Musual * ncoly)
    names(misc$earg) <- temp.names
    for(ii in 1:ncoly) {
      misc$earg[[Musual*ii-1]] <- .emean
      misc$earg[[Musual*ii  ]] <- if ( .var.arg) .e.var else .esd
    }
    names(misc$earg) <- temp.names

    misc$var.arg <- .var.arg
    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$imethod <- .imethod
    misc$multipleResponses <- TRUE
    misc$parallel <- .parallel
    misc$intercept.apply <- .intercept.apply
  }), list( .lmean = lmean, .lsd = lsd, .l.var = l.var,
            .emean = emean, .esd = esd, .e.var = e.var,
            .parallel = parallel, .intercept.apply = intercept.apply,
            .var.arg = var.arg, .imethod = imethod ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    ncoly <- extra$ncoly
    Musual <- extra$Musual
    if ( .var.arg ) {
      Varm <- eta2theta(eta[, Musual*(1:ncoly)], .l.var , earg = .e.var )
      sdev <- sqrt(Varm)
    } else {
      sdev <- eta2theta(eta[, Musual*(1:ncoly)], .lsd  , earg = .esd  )
    }
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dnorm(y, m = mu, sd = sdev, log = TRUE))
    }
  }, list( .lsd = lsd, .l.var = l.var,
           .esd = esd, .e.var = e.var,
           .var.arg = var.arg ))),
  vfamily = c("normal1"),
  deriv = eval(substitute(expression({
    ncoly <- extra$ncoly
    Musual <- extra$Musual

    mymu <- eta2theta(eta[, Musual*(1:ncoly) - 1], .lmean , earg = .emean )
    if ( .var.arg ) {
      Varm <- eta2theta(eta[, Musual*(1:ncoly)    ], .l.var , earg = .e.var )
      sdev <- sqrt(Varm)
    } else {
      sdev <- eta2theta(eta[, Musual*(1:ncoly)    ], .lsd  , earg = .esd  )
    }

    dl.dmu <- (y - mymu) / sdev^2
    if ( .var.arg ) {
      dl.dva <- -0.5 / Varm + 0.5 * (y - mymu)^2 / sdev^4
    } else {
      dl.dsd <- -1.0 / sdev +       (y - mymu)^2 / sdev^3
    }

    dmu.deta <- dtheta.deta(mymu, .lmean , earg = .emean )
    if ( .var.arg ) {
      dva.deta <- dtheta.deta(Varm, .l.var , earg = .e.var )
    } else {
      dsd.deta <- dtheta.deta(sdev, .lsd  , earg = .esd )
    }

    ans <- c(w) *
           cbind(dl.dmu * dmu.deta,
                 if ( .var.arg ) dl.dva * dva.deta else
                                 dl.dsd * dsd.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M = Musual)]
    ans
  }), list( .lmean = lmean, .lsd = lsd, .l.var = l.var,
            .emean = emean, .esd = esd, .e.var = e.var,
            .var.arg = var.arg ))),
  weight = eval(substitute(expression({
    wz <- matrix(as.numeric(NA), n, M) # diag matrix; y is 1-column too

    ned2l.dmu2 <- 1 / sdev^2
    if ( .var.arg ) {
      ned2l.dva2 <- 0.5 / Varm^2
    } else {
      ned2l.dsd2 <- 2 / sdev^2
    }

    wz[, Musual*(1:ncoly) - 1] <- ned2l.dmu2 * dmu.deta^2
    wz[, Musual*(1:ncoly)    ] <- if ( .var.arg ) {
      ned2l.dva2 * dva.deta^2
    } else {
      ned2l.dsd2 * dsd.deta^2
    }

    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = ncoly)
  }), list( .var.arg = var.arg ))))
}









 normal1.term <-
  function(linklist = NULL, # list(),
           earglist = NULL, # list(),
           lsd = "loge", lvar = "loge",
           esd = list(), evar = list(),
           var.arg = FALSE,
           imethod = 1,
           isd = NULL,
           ieta.coeffs = NULL,
           zero = "M")
{




 print("20120730; in normal1.term()")




  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")

  lvar <- as.list(substitute(lvar))
  e.var <- link2list(lvar)
  l.var <- attr(e.var, "function.name")




  if (is.character(zero) && zero != "M")
    stop("bad input for argument 'zero'")



  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
      stop("argument 'imethod' must be 1 or 2 or 3 or 4")


  if (!is.logical(var.arg) || length(var.arg) != 1)
    stop("argument 'var.arg' must be a single logical")



  new("vglmff",
  blurb = c("Univariate normal distribution with ",
            "varying coefficients links/constraints\n\n",
            "Links:    ",
            if (var.arg)
            namesof("var",  l.var, earg = e.var, tag = TRUE) else
            namesof("sd" ,  lsd,  earg = esd,  tag = TRUE), "; ",
            "\n",
            if (var.arg) "Variance: var" else "Variance: sd^2"),

  constraints = eval(substitute(expression({


    dotzero <- .zero
    if (is.character(dotzero) && dotzero == "M")
      dotzero <- M

    Musual <- M
    eval(negzero.expression)
  }), list( .zero = zero 
          ))),

  infos = eval(substitute(function(...) {
    list(Musual = NA,
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    asgn <- attr(x, "assign")
    nasgn <- names(asgn)
    asgn2 <- attr(Xm2, "assign")
    nasgn2 <- names(asgn2)


 print("head(x)")
 print( head(x) )
 print("head(Xm2)")
 print( head(Xm2) )


 print("attributes(x)")
 print( attributes(x) )
 print("attributes(Xm2)")
 print( attributes(Xm2) )





 print("names(constraints)")
 print( names(constraints) )
 print('nasgn')
 print( nasgn )
 print('nasgn2')
 print( nasgn2 )


    linklist <- .linklist
    Linklist <- vector("list", length(nasgn2))
    names(Linklist) <- nasgn2
    for (ilocal in 1:length(nasgn2))
      Linklist[[ilocal]] <- "identity"
    if (length( linklist ) > 0) {
      for (ilocal in 1:length(nasgn2))
        if (any(names(linklist) == nasgn2[ilocal]))
          Linklist[[ilocal]] <- linklist[[(nasgn2[ilocal])]]
    }
 print('linklist')
 print( linklist )
 print('Linklist')
 print( Linklist )
 print('unlist(Linklist)')
 print( unlist(Linklist) )







    orig.y <- y




    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    number.mlogit <- sum(unlist(Linklist) == "mlogit")
    print("number.mlogit")
    print( number.mlogit )
    if (number.mlogit == 1)
      stop('cannot have only one "mlogit"')


    ncoly <- ncol(y)
    Musual <- NA
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- ncol(Xm2) - (number.mlogit > 0) + 1
 print("M ,,,,,,,,,")
 print( M )
    extra$Xm2 <- Xm2



    cn.Xm2 <- colnames(Xm2)
    mynames1 <- NULL
    for (ilocal in 1:length(cn.Xm2))
      mynames1 <- c(mynames1,
                namesof(cn.Xm2[ilocal], Linklist[[ilocal]],
                        list(), tag = FALSE))

 print("mynames1")
 print( mynames1 )

    mynames2 <- paste(if ( .var.arg ) "var" else "sd",
                      if (ncoly > 1) 1:ncoly else "", sep = "")

    predictors.names <-
        c(mynames1,
          if ( .var.arg ) 
          namesof(mynames2, .l.var  , earg = .e.var  , tag = FALSE) else
          namesof(mynames2, .lsd   , earg = .esd   , tag = FALSE))
 print("predictors.names ,,,,,,,,,")
 print( predictors.names )
    extra$predictors.names <- predictors.names


    if (!length(etastart)) {
      sdev.init <- mean.init <- matrix(0, n, ncoly)
      for (jay in 1:ncoly) {
        jfit <- lm.wfit(x = Xm2,  y = y[, jay], w = w[, jay])
        mean.init[, jay] <- if ( mynames2 == "loge")
                            pmax(1/1024, y[, jay]) else
          if( .imethod == 1) median(y[, jay]) else
          if( .imethod == 2) weighted.mean(y[, jay], w = w[, jay]) else
          if( .imethod == 3) weighted.mean(y[, jay], w = w[, jay]) *
                             0.5 + y[, jay] * 0.5 else
                                 mean(jfit$fitted)

        sdev.init[, jay] <-
          if( .imethod == 1) {
            sqrt( sum(w[, jay] *
                (y[, jay] - mean.init[, jay])^2) / sum(w[, jay]) )
          } else if( .imethod == 2) {
            if (jfit$df.resid > 0)
              sqrt( sum(w[, jay] * jfit$resid^2) / jfit$df.resid ) else
              sqrt( sum(w[, jay] * jfit$resid^2) / sum(w[, jay]) )
          } else if( .imethod == 3) {
            sqrt( sum(w[, jay] * 
                  (y[, jay] - mean.init[, jay])^2) / sum(w[, jay]) )
          } else {
            sqrt( sum(w[, jay] * abs(y[, jay] -
                                     mean.init[, jay])) / sum(w[, jay]) )
          }

        if (any(sdev.init[, jay] <= sqrt( .Machine$double.eps ) ))
          sdev.init[, jay] <- 1.01

 print("head(sdev.init[, jay])9")
 print( head(sdev.init[, jay])  )
      }


      if (length( .i.sd )) {
        sdev.init <- matrix( .i.sd , n, ncoly, byrow = TRUE)
      }


      etastart <-
        cbind(eta.equi.probs,
              if ( .var.arg )
              theta2eta(sdev.init^2, .l.var , earg = .e.var ) else
              theta2eta(sdev.init  , .lsd  , earg = .esd  ))

      colnames(etastart) <- predictors.names
 print("head(etastart)9")
 print( head(etastart) )

      new.coeffs <- weighted.mean(y, w)
      extra$new.coeffs <- new.coeffs 

    }
  }), list( .linklist = linklist,
            .earglist = earglist,
            .lsd = lsd, .l.var = lvar,
            .esd = esd, .e.var = evar,
            .i.sd = isd,
            .ieta.coeffs = ieta.coeffs,
            .var.arg = var.arg, .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
 print("hi9")

    M <- ncol(eta)
    betas.matrix <- 1 / (1 + exp(-eta[, -M, drop = FALSE]))
    betas.matrix <- cbind(betas.matrix,
                          1 / (1 + rowSums(exp(eta[, -M, drop = FALSE]))))
 print("head(betas.matrix)1")
 print( head(betas.matrix) )

      betas.matrix <- cbind(extra$new.coeffs[1], betas.matrix)

 print("head(betas.matrix)2")
 print( head(betas.matrix) )
 print("head(extra$Xm2)")
 print( head(extra$Xm2) )



    rowSums(extra$Xm2 * betas.matrix)
  }, list( .linklist = linklist,
           .earglist = earglist,
           .esd = esd , .e.var = evar ))),

  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <- c(rep( "mlogit", length = M - 1),
                   rep( .lsd   , length = ncoly))
    temp.names <- c(mynames1, mynames2)
    names(misc$link) <- temp.names




    misc$var.arg <- .var.arg
    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$imethod <- .imethod
    misc$multipleResponses <- FALSE
  }), list( .linklist = linklist,
            .earglist = earglist,
            .lsd = lsd, .l.var = lvar,
            .esd = esd, .e.var = evar,
            .var.arg = var.arg, .imethod = imethod ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    ncoly <- extra$ncoly
    Musual <- 1 # extra$Musual
    if ( .var.arg ) {
      Varm <- eta2theta(eta[, Musual*(1:ncoly)], .l.var , earg = .e.var )
      sdev <- sqrt(Varm)
    } else {
      sdev <- eta2theta(eta[, Musual*(1:ncoly)], .lsd  , earg = .esd  )
    }
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dnorm(y, m = mu, sd = sdev, log = TRUE))
    }
  }, list( .lsd = lsd, .l.var = lvar,
           .esd = esd, .e.var = evar,
           .var.arg = var.arg ))),
  vfamily = c("normal1.term"),
  deriv = eval(substitute(expression({
 print("------ in @ deriv -------------")
    extra$new.coeffs <- new.coeffs

    ncoly <- extra$ncoly
    Musual <- 1 # extra$Musual

    if ( .var.arg ) {
      Varm <- eta2theta(eta[, Musual*(1:ncoly)    ], .l.var , earg = .e.var )
      sdev <- sqrt(Varm)
    } else {
      sdev <- eta2theta(eta[, Musual*(1:ncoly)    ], .lsd  , earg = .esd  )
    }



    betas.matrix <- 1 / (1 + exp(-eta[, -M, drop = FALSE]))
    betas.matrix <- cbind(betas.matrix,
                          1 / (1 + rowSums(exp(eta[, -M, drop = FALSE]))))
 print("head(betas.matrix)5")
 print( head(betas.matrix) )

    if ( !extra$sum1.intercept &&
          any(colnames(extra$X_LM) == "(Intercept)"))
      betas.matrix <- cbind(extra$new.coeffs[1], betas.matrix)

 print("head(betas.matrix)6")
 print( head(betas.matrix) )
 print("head(extra$Xm2)")
 print( head(extra$Xm2) )

    use.x <- if ( sum1.intercept )
             Xm2[, -ncol(Xm2), drop = FALSE] else
             Xm2[, -c(1, ncol(Xm2)), drop = FALSE]
    mymu <- rowSums(Xm2 * betas.matrix)
    dMu.deta <- mymu * (1 - mymu) * use.x


 print("head(mymu)9")
 print( head(mymu) )
 print("head(dMu.deta)9")
 print( head(dMu.deta) )
    if ( .var.arg ) {
      dl.dva <- -0.5 / Varm + 0.5 * (y - mymu)^2 / sdev^4
    } else {
      dl.dsd <- -1.0 / sdev +       (y - mymu)^2 / sdev^3
    }
    dl.dmu <- (y - mymu) / sdev^2


    if ( .var.arg ) {
      dva.deta <- dtheta.deta(Varm, .l.var , earg = .e.var )
    } else {
      dsd.deta <- dtheta.deta(sdev, .lsd  , earg = .esd )
    }

    ans <- c(w) *
           cbind(dl.dmu * dMu.deta,
                 if ( .var.arg ) dl.dva * dva.deta else
                                 dl.dsd * dsd.deta)
 print("head(deriv.ans)9")
 print( head(ans) )
    ans
  }), list( .linklist = linklist, .lsd = lsd, .l.var = lvar,
            .earglist = earglist, .esd = esd, .e.var = evar,
            .var.arg = var.arg ))),
  weight = eval(substitute(expression({
 print("------ in @ weight -------------")
    wz <- matrix(0, n, dimm(M)) # diag matrix; y is 1-column too
 print("head(wz)")
 print( head(wz) )

    if ( .var.arg ) {
      ned2l.dva2 <- 0.5 / Varm^2
    } else {
      ned2l.dsd2 <- 2 / sdev^2
    }





    wz[, iam(M, M, M = M)] <- if ( .var.arg ) {
      ned2l.dva2 * dva.deta^2
    } else {
      ned2l.dsd2 * dsd.deta^2
    }


    index = iam(NA, NA, M  , both = TRUE, diag = TRUE)
    indtw = iam(NA, NA, M-1, both = TRUE, diag = TRUE)
 print("index")
 print( index )
 print("indtw")
 print( indtw )

    
    twz = dMu.deta[, indtw$row.index, drop = FALSE] *
          dMu.deta[, indtw$col.index, drop = FALSE]
 print("head(twz)9------------------------------------------------")
 print( head(twz) )


    for (ilocal in 1:ncol(twz))
      wz[, iam(index$row.index[ilocal],
               index$col.index[ilocal], M = M)] <-
     twz[, iam(indtw$row.index[ilocal],
               indtw$col.index[ilocal], M = M-1)]


 print("head(wz)9------------------------------------------------")
 print( head(wz) )

    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = ncoly)
  }), list( .var.arg = var.arg ))))
} # End of normal1.term()





 lognormal <- function(lmeanlog = "identity", lsdlog = "loge",
                       zero = 2)
{




  lmulog <- as.list(substitute(lmeanlog))
  emulog <- link2list(lmulog)
  lmulog <- attr(emulog, "function.name")

  lsdlog <- as.list(substitute(lsdlog))
  esdlog <- link2list(lsdlog)
  lsdlog <- attr(esdlog, "function.name")




  if (length(zero) &&
     (!is.Numeric(zero, integer.valued = TRUE, positive = TRUE) ||
     zero > 2))
    stop("bad input for argument argument 'zero'")


  new("vglmff",
  blurb = c("Two-parameter (univariate) lognormal distribution\n\n",
          "Links:    ",
          namesof("meanlog", lmulog, earg = emulog, tag = TRUE), ", ",
          namesof("sdlog",   lsdlog, earg = esdlog, tag = TRUE)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              Is.positive.y = TRUE)



    predictors.names <-
        c(namesof("meanlog", .lmulog, earg = .emulog, tag = FALSE),
          namesof("sdlog",   .lsdlog, earg = .esdlog, tag = FALSE))

    if (!length(etastart)) {
      mylm = lm.wfit(x = x, y = log(y), w = w)
      sdlog.y.est = sqrt( sum(c(w) * mylm$resid^2) / mylm$df.residual )
      etastart = cbind(
        meanlog = rep(theta2eta(log(median(y)), .lmulog,
                                earg = .emulog), length = n),
        sdlog   = rep(theta2eta(sdlog.y.est, .lsdlog,
                                earg = .esdlog), length = n))
    }
  }), list( .lmulog = lmulog, .lsdlog = lsdlog,
            .emulog = emulog, .esdlog = esdlog ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mulog = eta2theta(eta[, 1], .lmulog , earg = .emulog )
    sdlog = eta2theta(eta[, 2], .lsdlog , earg = .esdlog )
    exp(mulog + 0.5 * sdlog^2)
  }, list( .lmulog = lmulog, .lsdlog = lsdlog,
           .emulog = emulog, .esdlog = esdlog ))),
  last = eval(substitute(expression({
    misc$link =    c("meanlog" = .lmulog , "sdlog" = .lsdlog )
    misc$earg = list("meanlog" = .emulog , "sdlog" = .esdlog )

    misc$expected = TRUE
  }), list( .lmulog = lmulog, .lsdlog = lsdlog,
            .emulog = emulog, .esdlog = esdlog ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    mulog = eta2theta(eta[, 1], .lmulog, earg = .emulog)
    sdlog = eta2theta(eta[, 2], .lsdlog,   earg = .esdlog)
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dlnorm(y, meanlog = mulog, sdlog = sdlog, log = TRUE))
    }
  }, list( .lmulog = lmulog, .lsdlog = lsdlog,
           .emulog = emulog, .esdlog = esdlog ))),
  vfamily = c("lognormal"),
  deriv = eval(substitute(expression({
    mulog = eta2theta(eta[, 1], .lmulog, earg = .emulog)
    sdlog = eta2theta(eta[, 2], .lsdlog, earg = .esdlog)

    dmulog.deta = dtheta.deta(mulog, .lmulog, earg = .emulog)
    dsdlog.deta = dtheta.deta(sdlog, .lsdlog,   earg = .esdlog)

    dl.dmulog = (log(y) - mulog) / sdlog^2
    dl.dsdlog = -1 / sdlog + (log(y) - mulog)^2 / sdlog^3

    c(w) * cbind(dl.dmulog * dmulog.deta, 
                 dl.dsdlog * dsdlog.deta)
  }), list( .lmulog = lmulog, .lsdlog = lsdlog,
            .emulog = emulog, .esdlog = esdlog ))),
  weight = expression({
    wz = matrix(as.numeric(NA), n, 2) # Diagonal!
    ned2l.dmulog2 = 1 / sdlog^2
    ned2l.dsdlog2 = 2 * ned2l.dmulog2

    wz[, iam(1, 1, M)] = ned2l.dmulog2 * dmulog.deta^2
    wz[, iam(2, 2, M)] = ned2l.dsdlog2 * dsdlog.deta^2

    wz = c(w) * wz
    wz
  }))
}






 lognormal3 <- function(lmeanlog = "identity", lsdlog = "loge",
                        powers.try = (-3):3,
                        delta = NULL, zero = 2)
{


  if (length(delta) &&
      !is.Numeric(delta, positive = TRUE))
    stop("bad input for argument argument 'delta'")



  lmulog <- as.list(substitute(lmeanlog))
  emulog <- link2list(lmulog)
  lmulog <- attr(emulog, "function.name")

  lsdlog <- as.list(substitute(lsdlog))
  esdlog <- link2list(lsdlog)
  lsdlog <- attr(esdlog, "function.name")




  if (length(zero) &&
     (!is.Numeric(zero, integer.valued = TRUE, positive = TRUE) ||
     zero > 3))
    stop("bad input for argument argument 'zero'")





  new("vglmff",
  blurb = c("Three-parameter (univariate) lognormal distribution\n\n",
          "Links:    ",
          namesof("meanlog", lmulog, earg = emulog, tag = TRUE), "; ",
          namesof("sdlog",   lsdlog, earg = esdlog, tag = TRUE), "; ",
          namesof("lambda", "identity", earg = list(), tag = TRUE)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y)



    predictors.names <- 
       c(namesof("meanlog", .lmulog, earg = .emulog, tag = FALSE), 
         namesof("sdlog",   .lsdlog, earg = .esdlog, tag = FALSE),
         "lambda")

    if (!length(etastart)) {
      miny = min(y)
      if (length( .delta)) {
        lambda.init = rep(miny- .delta, length = n)
      } else {
        pvalue.vec = NULL
        powers.try = .powers.try
        for(delta in 10^powers.try) {
          pvalue.vec = c(pvalue.vec,
                         shapiro.test(sample(log(y-miny+delta),
                         size=min(5000, length(y ))))$p.value) 
        }
        index.lambda = (1:length(powers.try))[pvalue.vec ==
                                              max(pvalue.vec)]
        lambda.init = miny - 10^powers.try[index.lambda]
      }
      mylm = lm.wfit(x = x, y=log(y-lambda.init), w = w)
      sdlog.y.est = sqrt( sum(c(w) * mylm$resid^2) / mylm$df.residual )
      etastart = cbind(mu = log(median(y - lambda.init)),
           sdlog = rep(theta2eta(sdlog.y.est, .lsdlog, earg = .esdlog),
                       length = n),
           lambda = lambda.init)
    }
  }), list( .lmulog = lmulog, .lsdlog = lsdlog,
            .emulog = emulog, .esdlog = esdlog,
            .delta = delta, .powers.try = powers.try ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mymu    = eta2theta(eta[, 1], .lmulog, earg = .emulog)
    sdlog   = eta2theta(eta[, 2], .lsdlog, earg = .esdlog)
    lambda  = eta2theta(eta[, 3], "identity", earg = list(theta = NULL))
    lambda + exp(mymu + 0.5 * sdlog^2)
  }, list( .lmulog = lmulog, .lsdlog = lsdlog,
           .emulog = emulog, .esdlog = esdlog ))),
  last = eval(substitute(expression({
    misc$link =    c("meanlog" = .lmulog,
                     "sdlog"   = .lsdlog,
                     "lambda"  = "identity")

    misc$earg = list("meanlog" = .emulog,
                     "sdlog"   = .esdlog,
                     "lambda"  = list())

    misc$expected = TRUE
  }), list( .lmulog = lmulog, .lsdlog = lsdlog,
            .emulog = emulog, .esdlog = esdlog ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    mymu   = eta2theta(eta[, 1], .lmulog , earg = .emulog)
    sdlog  = eta2theta(eta[, 2], .lsdlog , earg = .esdlog)
    lambda = eta2theta(eta[, 3], "identity", earg = list(theta = NULL))
    if (any(y < lambda))
      warning("bad 'y'")
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
          sum(c(w) * dlnorm(y - lambda, meanlog = mymu,
                            sdlog = sdlog, log = TRUE))
      }
  }, list( .lmulog = lmulog, .lsdlog = lsdlog,
           .emulog = emulog, .esdlog = esdlog ))),
  vfamily = c("lognormal3"),
  deriv = eval(substitute(expression({
    mymu   = eta2theta(eta[, 1], .lmulog,  earg = .emulog)
    sdlog  = eta2theta(eta[, 2], .lsdlog,    earg = .esdlog)
    lambda = eta2theta(eta[, 3], "identity", earg = list(theta = NULL))

    if (any(y < lambda))
        warning("bad 'y'")

    dl.dmymu <- (log(y-lambda)-mymu) / sdlog^2
    dl.dsdlog <- -1/sdlog + (log(y-lambda)-mymu)^2 / sdlog^3
    dl.dlambda <- (1 + (log(y-lambda)-mymu) / sdlog^2) / (y-lambda)

    dmymu.deta <- dtheta.deta(mymu, .lmulog, earg = .emulog)
    dsdlog.deta <- dtheta.deta(sdlog, .lsdlog, earg = .esdlog)
    dlambda.deta <- dtheta.deta(lambda, "identity", earg = list())

    c(w) * cbind(dl.dmymu   *   dmymu.deta, 
                 dl.dsdlog  *  dsdlog.deta, 
                 dl.dlambda * dlambda.deta)
  }), list( .lmulog = lmulog, .lsdlog = lsdlog,
            .emulog = emulog, .esdlog = esdlog ))),
  weight = expression({
    wz <- matrix(0, n, dimm(M))

    ned2l.dmymu2 <- 1 / sdlog^2
    ned2l.dsdlog <- 2 / sdlog^2
    temp9 <- exp(-mymu + sdlog^2 / 2)
    ned2l.dlambda2 <- exp(2*(-mymu+sdlog^2)) * (1+sdlog^2) / sdlog^2

    wz[, iam(1, 1, M)] <- ned2l.dmymu2 * dmymu.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dsdlog * dsdlog.deta^2
    wz[, iam(3, 3, M)] <- ned2l.dlambda2 * dlambda.deta^2
    wz[, iam(1, 3, M)] <- temp9 * dmymu.deta * dlambda.deta / sdlog^2
    wz[, iam(2, 3, M)] <- -2 * temp9 / sdlog * dsdlog.deta * dlambda.deta
    wz <- c(w) * wz
    wz
  }))
}





dsnorm <- function(x, location = 0, scale = 1, shape = 0, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if (!is.Numeric(scale, positive = TRUE))
    stop("bad input for argument 'scale'")
  zedd <- (x - location) / scale
  loglik <- log(2) + dnorm(zedd, log = TRUE) +
           pnorm(shape * zedd, log.p = TRUE) -
           log(scale)
  if (log.arg) {
    loglik
  } else {
    exp(loglik)
  }
}



rsnorm <- function(n, location = 0, scale = 1, shape=0) {
  if (!is.Numeric(n, positive = TRUE,
                  integer.valued = TRUE, allowable.length = 1))
    stop("bad input for argument 'n'")
  if (!is.Numeric(scale, positive = TRUE))
    stop("bad input for argument 'scale'")
  if (!is.Numeric(shape))
    stop("bad input for argument 'shape'")

  rho <- shape / sqrt(1 + shape^2)
  u0 <- rnorm(n)
  v <- rnorm(n)
  u1 <- rho*u0 + sqrt(1 - rho^2) * v
  location + scale * ifelse(u0 >= 0, u1, -u1)
}




 skewnormal1 <- function(lshape = "identity",
                         ishape = NULL,
                         nsimEIM = NULL)
{


  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  if (length(nsimEIM) &&
     (!is.Numeric(nsimEIM, allowable.length = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 10))
    stop("argument 'nsimEIM' should be an integer greater than 10")


  new("vglmff",
  blurb = c("1-parameter skew-normal distribution\n\n",
          "Link:     ",
          namesof("shape", lshape , earg = eshape ), "\n",
          "Mean:     shape * sqrt(2 / (pi * (1 + shape^2 )))\n",
          "Variance: 1-mu^2"),
  infos = eval(substitute(function(...) {
    list(Musual = 1,
         nsimEIM = .nsimEIM)
  }, list( .nsimEIM = nsimEIM ))),
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
      namesof("shape", .lshape , earg = .eshape , tag = FALSE)

    if (!length(etastart)) {
      init.shape <- if (length( .ishape))
        rep( .ishape, len = n) else {
        temp <- y
        index <- abs(y) < sqrt(2/pi)-0.01
        temp[!index] <- y[!index]
        temp[index] <- sign(y[index]) / sqrt(2/(pi*y[index]*y[index])-1)
        temp
      }
      etastart <- matrix(init.shape, n, ncol(y))
    }
  }), list( .lshape = lshape, .eshape = eshape,
            .ishape = ishape ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha <- eta2theta(eta, .lshape, earg = .eshape)
    alpha * sqrt(2/(pi * (1+alpha^2 )))
  }, list( .eshape = eshape, .lshape = lshape ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape = .lshape) 

    misc$earg <- list(shape = .eshape )

    misc$nsimEIM = .nsimEIM
      misc$expected <- (length( .nsimEIM ) > 0)
  }), list( .eshape = eshape, .lshape = lshape,
            .nsimEIM = nsimEIM ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    alpha <- mu / sqrt(2/pi - mu^2)
    theta2eta(alpha, .lshape, earg = .eshape)
  }, list( .eshape = eshape, .lshape = lshape ))),
  loglikelihood = eval(substitute(
     function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        alpha <- eta2theta(eta, .lshape, earg = .eshape)
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dsnorm(x = y, location = 0, scale = 1,
                        shape = alpha, log = TRUE))
    }
  }, list( .eshape = eshape, .lshape = lshape ))), 
  vfamily = c("skewnormal1"),
  deriv = eval(substitute(expression({
    alpha <- eta2theta(eta, .lshape, earg = .eshape)

    zedd <- y*alpha
    tmp76 <- pnorm(zedd)
    tmp86 <- dnorm(zedd)
    dl.dshape <- tmp86 * y / tmp76

    dshape.deta <- dtheta.deta(alpha, .lshape, earg = .eshape)

    c(w) * dl.dshape * dshape.deta
  }), list( .eshape = eshape, .lshape = lshape ))),
  weight = eval(substitute(expression({
    if ( length( .nsimEIM )) {
      run.mean = 0
      for(ii in 1:( .nsimEIM)) {
          ysim = rsnorm(n, location = 0, scale = 1, shape = alpha)
          zedd = ysim*alpha
          tmp76 = pnorm(zedd)
          tmp86 = dnorm(zedd)
          d2l.dshape2 = -ysim*ysim*tmp86*(tmp76*zedd+tmp86)/tmp76^2
          rm(ysim)
          run.mean = ((ii-1) * run.mean + d2l.dshape2) / ii
      }
      if (intercept.only)
        run.mean = mean(run.mean)
      wz =  -c(w) * (dshape.deta^2) * run.mean
    } else {
      d2shape.deta2 = d2theta.deta2(alpha, .lshape, earg = .eshape)
      d2l.dshape2 = -y*y * tmp86 * (tmp76 * zedd + tmp86) / tmp76^2
      wz = -(dshape.deta^2) * d2l.dshape2 - d2shape.deta2 * dl.dshape
      wz = c(w) * wz
    }
    wz
  }), list( .eshape = eshape, .lshape = lshape, .nsimEIM = nsimEIM ))))
}










# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.












N.hat.posbernoulli <-
  function(eta, link, earg = list(),
           R = NULL, w = NULL,
           X_vlm = NULL, Hlist = NULL,
           extra = list(),
           model.type = c("b", "t", "tb")
          ) {



  if (!is.null(w) && !all(1 == w))
    warning("estimate of N may be wrong when prior weights ",
            "are not all unity")

  model.type <- match.arg(model.type, c("b", "t", "tb"))[1]

  tau <-
    switch(model.type,
           "b"  = extra$tau,
           "t"  = ncol(eta),
           "tb" = (ncol(eta) + 1) / 2)
  if (length(extra$tau) && extra$tau != tau)
    warning("variable 'tau' is mistaken")  # Checking only

  jay.index <-
    switch(model.type,
           "b"  = rep(1, length = tau),  # Subset: 1 out of 1:2
           "t"  = 1:tau,  # All of them
           "tb" = 1:tau)  # Subset: first tau of them out of M = 2*tau-1
  prc <- eta2theta(eta[, jay.index], link, earg = earg)  # cap.probs
  QQQ <- exp(rowSums(log1p(-prc)))
  pibbeta <- exp(log1p(-QQQ))  # One.minus.QQQ
  N.hat <- sum(1 / pibbeta)  # Point estimate
  ss2 <- sum(QQQ / pibbeta^2)  # Assumes bbeta is known


  if (length(R)) {

    dvect <- matrix(0, length(pibbeta), ncol = ncol(X_vlm))
    M <- nrow(Hlist[[1]])
    n_lm <- nrow(X_vlm) / M  # Number of rows of the LM matrix
    dprc.deta <- dtheta.deta(prc, link, earg = earg)
    Hmatrices <- matrix(c(unlist(Hlist)), nrow = M)
    for (jay in 1:tau) {
      lapred.index <- jay.index[jay]
      Index0 <- Hmatrices[lapred.index, ] != 0
      X_lm_jay <- X_vlm[(0:(n_lm - 1)) * M + lapred.index, Index0,
                        drop = FALSE]

      dvect[, Index0] <-
      dvect[, Index0] + (QQQ / (1-prc[, jay])) * dprc.deta[, jay] * X_lm_jay
    }


   dvect <- dvect * (-1 / pibbeta^2)
   dvect <- colSums(dvect)  # Now a vector

    ncol_X_vlm <- nrow(R)
    rinv <- diag(ncol_X_vlm)
    rinv <- backsolve(R, rinv)
    rowlen <- drop(((rinv^2) %*% rep(1, ncol_X_vlm))^0.5)
    covun <- rinv %*% t(rinv)
    vecTF <- FALSE
    for (jay in 1:tau) {
      lapred.index <- jay.index[jay]
      vecTF <- vecTF | (Hmatrices[lapred.index, ] != 0)
    }
    vecTF.index <- (1:length(vecTF))[vecTF]
    covun <- covun[vecTF.index, vecTF.index, drop = FALSE]
    dvect <- dvect[vecTF.index, drop = FALSE]
  }
 
  list(N.hat    = N.hat,
       SE.N.hat = if (length(R)) sqrt(ss2 + t(dvect) %*% covun %*% dvect) else
                                 sqrt(ss2)
      )
}




aux.posbernoulli <- function(y, check.y = FALSE) {








  y <- as.matrix(y)
  if ((tau <- ncol(y)) == 1)
    stop("argument 'y' needs to be a matrix with at least two columns")
  if (check.y) {
    if (!all(y == 0 | y == 1 | y == 1/tau | is.na(y)))
      stop("response 'y' must contain 0s and 1s only")
  }

  zeddij <- cbind(0, t(apply(y, 1, cumsum))) # tau + 1 columns
  zij <- (0 + (zeddij > 0))[, 1:tau] # 0 or 1.
  if (length(colnames(y)))
    colnames(zij) <- colnames(y)

  cp1 <- numeric(nrow(y))
  for (jay in tau:1)
    cp1[y[, jay] > 0] <- jay
  if (any(cp1 == 0))
    warning("some individuals were never captured!")

  yr1i <- zeddij[, tau + 1] - 1
  list(cap.hist1 = zij,
       cap1      = cp1, # aka ti1
       y0i       = cp1 - 1,
       yr0i      = tau - cp1 - yr1i,
       yr1i      = yr1i)
}










rposbern <-
  function(n, nTimePts = 5, pvars = length(xcoeff),
           xcoeff = c(-2, 1, 2),
           cap.effect = -1,
           link = "logit",
           is.popn = FALSE,
           earg.link = FALSE) {








  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
               stop("bad input for argument 'n'") else n
  orig.n <- use.n
  if (!is.popn)
    use.n <- 1.50 * use.n + 100 # Bigger due to rejections

  if (pvars == 0)
    stop("argument 'pvars' must be at least one")
  if (pvars > length(xcoeff))
    stop("argument 'pvars' is too high")
  

  if (earg.link) {
    earg <- link
  } else {
    link <- as.list(substitute(link))
    earg <- link2list(link)
  }
  link <- attr(earg, "function.name")




  Ymatrix <- matrix(0, use.n, nTimePts,
                    dimnames = list(as.character(1:use.n),
                                    paste("y", 1:nTimePts, sep = "")))

  CHmatrix <- matrix(0, use.n, nTimePts,
                     dimnames = list(as.character(1:use.n),
                                     paste("ch", 0:(nTimePts-1), sep = "")))

  Xmatrix <- cbind(x1 = rep(1.0, len = use.n))
  if (pvars > 1)
    Xmatrix <- cbind(Xmatrix,
                     matrix(runif(n = use.n * (pvars-1)),
                            use.n, pvars - 1,
                            dimnames = list(as.character(1:use.n),
                                            paste("x", 2:pvars, sep = ""))))


  lin.pred.baseline <- xcoeff[1]
  if (pvars > 1)
    lin.pred.baseline <- lin.pred.baseline +
                         Xmatrix[, 2:pvars, drop = FALSE] %*%
                         xcoeff[2:pvars]
  sumrowy <- rep(0, length = use.n)

  for (jlocal in 1:nTimePts) {
    CHmatrix[, jlocal] <- as.numeric(sumrowy > 0)

    lin.pred <- lin.pred.baseline + (CHmatrix[, jlocal] >  0) * cap.effect

    Ymatrix[, jlocal] <-
      rbinom(use.n, size = 1,
             prob = eta2theta(lin.pred, link = link, earg = earg))

    sumrowy <- sumrowy + Ymatrix[, jlocal]
  }


  index0 <- (sumrowy == 0)
  if (all(!index0))
    stop("bug in this code: cannot handle no animals being caught")
   Ymatrix <-  Ymatrix[!index0, , drop = FALSE]
   Xmatrix <-  Xmatrix[!index0, , drop = FALSE]
  CHmatrix <- CHmatrix[!index0, , drop = FALSE]

  zCHmatrix <- matrix(0, nrow(CHmatrix), ncol(CHmatrix),
                      dimnames = list(as.character(1:nrow(CHmatrix)),
                      paste("zch", 0:(ncol(CHmatrix)-1), sep = "")))


  ans <- data.frame(Ymatrix, Xmatrix, CHmatrix, zCHmatrix,
                    Chistory = rep(0, length = nrow(Ymatrix)))


  if (!is.popn) {
    ans <- if (nrow(ans) >= orig.n) {
      ans[1:orig.n, ]
    } else {
      rbind(ans,
            Recall(n = orig.n - nrow(ans),
                   nTimePts = nTimePts, pvars = pvars,
                   xcoeff = xcoeff,
                   cap.effect = cap.effect,
                   link = earg, earg.link = TRUE))
    }
  }

  rownames(ans) <- as.character(1:nrow(ans))

  attr(ans, "pvars")      <- pvars
  attr(ans, "nTimePts")   <- nTimePts
  attr(ans, "cap.effect") <- cap.effect
  attr(ans, "is.popn")    <- is.popn
  attr(ans, "n")          <- n

  ans
}





  

dposbern <- function(x, prob, prob0 = prob, log = FALSE) {


  x     <- as.matrix(x)
  prob  <- as.matrix(prob)
  prob0 <- as.matrix(prob0)

  if (!is.logical(log.arg <- log) ||
      length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)
  if (ncol(x) < 2)
    stop("columns of argument 'x' should be 2 or more")


  logAA0 <- rowSums(log1p(-prob0))
  AA0 <- exp(logAA0)
  ell1 <- x * log(prob) + (1 - x) * log1p(-prob) - log1p(-AA0) / ncol(x)
  if (log.arg) ell1 else exp(ell1)
}








dposnegbin <- function(x, size, prob = NULL, munb = NULL, log = FALSE) {
  if (length(munb)) {
    if (length(prob))
      stop("'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  LLL <- max(length(x), length(prob), length(size))
  x    <- rep(x,    len = LLL);
  prob <- rep(prob, len = LLL);
  size <- rep(size, len = LLL);

  ans <- dnbinom(x = x, size = size, prob = prob, log = log.arg)
  index0 <- (x == 0)

  if (log.arg) {
    ans[ index0] <- log(0.0)
    ans[!index0] <- ans[!index0] - log1p(-dnbinom(x = 0 * x[!index0],
                    size = size[!index0], prob = prob[!index0]))
  } else {
    ans[ index0] <- 0.0
    ans[!index0] <- ans[!index0] / pnbinom(q = 0 * x[!index0],
                    size = size[!index0], prob = prob[!index0],
                    lower.tail = FALSE)
  }
  ans
}


pposnegbin <- function(q, size, prob = NULL, munb = NULL) {

  if (length(munb)) {
    if (length(prob))
      stop("'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }
  L <- max(length(q), length(prob), length(size))
  if (length(q)    != L)
    q    <- rep(q,    length.out = L);
  if (length(prob) != L)
    prob <- rep(prob, length.out = L);
  if (length(size) != L)
    size <- rep(size, length.out = L)

  ifelse(q < 1, 0,
        (pnbinom(q, size = size, prob = prob) -
         dnbinom(0, size = size, prob = prob))
       / pnbinom(0, size = size, prob = prob, lower.tail = FALSE))
}


qposnegbin <- function(p, size, prob = NULL, munb = NULL) {


  if (length(munb)) {
    if (length(prob))
      stop("'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  ans <- qnbinom(pnbinom(q = 0, size = size, prob = prob,
                         lower.tail = FALSE) * p +
                 dnbinom(x = 0, size = size, prob = prob),
                 size = size, prob = prob)
  ans[p >  1] <- NaN
  ans[p <  0] <- NaN
  ans[p == 1] <- Inf
  ans
}









posnegbinomial.control <- function(save.weight = TRUE, ...) {
  list(save.weight = save.weight)
}



 posnegbinomial <- function(lmunb = "loge", lsize = "loge",
                            isize = NULL, zero = -2,
                            nsimEIM = 250,
                            shrinkage.init = 0.95, imethod = 1) {

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
      stop("bad input for argument 'isize'")
  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
     shrinkage.init < 0 ||
     shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")


  lmunb <- as.list(substitute(lmunb))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")


  if (!is.Numeric(nsimEIM, allowable.length = 1,
                  positive = TRUE, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 30)
    warning("argument 'nsimEIM' should be greater than 30, say")


  new("vglmff",
  blurb = c("Positive-negative binomial distribution\n\n",
            "Links:    ",
            namesof("munb", lmunb, earg = emunb ), ", ",
            namesof("size", lsize, earg = esize ), "\n",
            "Mean:     munb / (1 - (size / (size + munb))^size)"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(Musual = 2,
         lmunb = .lmunb ,
         emunb = .emunb ,
         lsize = .lsize ,
         esize = .esize )
  }, list( .lmunb = lmunb, .lsize = lsize, .isize = isize,
            .emunb = emunb, .esize = esize,
            .sinit = shrinkage.init,
            .imethod = imethod ))),

  initialize = eval(substitute(expression({
    Musual <- 2

    if (any(y == 0))
      stop("there are zero values in the response")
    y <- as.matrix(y) 


    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y





    M <- Musual * ncol(y) 
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species

    predictors.names <- c(
      namesof(if (NOS == 1) "munb" else
              paste("munb", 1:NOS, sep = ""),
              .lmunb, earg = .emunb, tag = FALSE),
      namesof(if (NOS == 1) "size" else
              paste("size", 1:NOS, sep = ""),
              .lsize, earg = .esize, tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M = Musual)]

    if (!length(etastart)) {
      mu.init <- y
      for(iii in 1:ncol(y)) {
        use.this <- if ( .imethod == 1) {
          weighted.mean(y[, iii], w[, iii])
        } else {
          median(y[,iii])
        }
        mu.init[, iii] <- (1 - .sinit) * y[, iii] + .sinit * use.this
      }

      if ( is.Numeric( .isize )) {
        kmat0 <- matrix( .isize , nrow = n, ncol = NOS, byrow = TRUE)
      } else {
        posnegbinomial.Loglikfun =
            function(kmat, y, x, w, extraargs) {
            munb <- extraargs
              sum(w * dposnegbin(x = y, size = kmat, munb = munb,
                                 log = TRUE))
              }
            k.grid <- 2^((-6):6)
            kmat0 <- matrix(0, nrow = n, ncol = NOS)
            for(spp. in 1:NOS) {
              kmat0[, spp.] <- getMaxMin(k.grid,
                                objfun = posnegbinomial.Loglikfun,
                                y = y[, spp.], x = x, w = w[, spp.],
                                extraargs = mu.init[, spp.])
            }
      }
      p00 <- (kmat0 / (kmat0 + mu.init))^kmat0
      etastart <-
        cbind(
              theta2eta(mu.init * (1 - p00), .lmunb, earg = .emunb ),
              theta2eta(kmat0,               .lsize, earg = .esize ))
      etastart <- etastart[,interleave.VGAM(M, M = Musual), drop = FALSE]
    }
  }), list( .lmunb = lmunb, .lsize = lsize, .isize = isize,
            .emunb = emunb, .esize = esize,
            .sinit = shrinkage.init,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Musual <- 2
    NOS <- ncol(eta) / Musual
    munb <- eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                     .lmunb, earg = .emunb )
    kmat <- eta2theta(eta[, Musual*(1:NOS),   drop = FALSE],
                     .lsize, earg = .esize )
    po0 <- (kmat / (kmat + munb))^kmat
    munb / (1 - po0)
  }, list( .lsize = lsize, .lmunb = lmunb,
           .esize = esize, .emunb = emunb ))),
  last = eval(substitute(expression({
    temp0303 <- c(rep( .lmunb , length = NOS),
                  rep( .lsize , length = NOS))
    names(temp0303) =
       c(if (NOS == 1) "munb" else paste("munb", 1:NOS, sep = ""),
         if (NOS == 1) "size" else paste("size", 1:NOS, sep = ""))
    temp0303  <- temp0303[interleave.VGAM(M, M = Musual)]
    misc$link <- temp0303  # Already named

    misc$earg <- vector("list", Musual*NOS)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:NOS) {
      misc$earg[[Musual*ii-1]] <- .emunb
      misc$earg[[Musual*ii  ]] <- .esize
    }

    misc$nsimEIM <- .nsimEIM
    misc$imethod <- .imethod
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize,
            .nsimEIM = nsimEIM, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Musual <- 2
    NOS <- ncol(eta) / Musual
    munb <- eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                      .lmunb, earg = .emunb )
    kmat <- eta2theta(eta[, Musual*(1:NOS)  , drop = FALSE],
                      .lsize, earg = .esize )
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(w * dposnegbin(x = y, size = kmat, munb = munb, log = TRUE))
    }
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize ))),

  vfamily = c("posnegbinomial"),
  deriv = eval(substitute(expression({
    Musual <- 2
    NOS <- extra$NOS

    munb <- eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                      .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, Musual*(1:NOS)  , drop = FALSE],
                      .lsize , earg = .esize )

    dmunb.deta <- dtheta.deta(munb, .lmunb, earg = .emunb )
    dsize.deta <- dtheta.deta(kmat, .lsize, earg = .esize )
    NOS <- ncol(eta) / Musual


    tempk <- kmat / (kmat + munb)
    tempm <- munb / (kmat + munb)
    prob0  <- tempk^kmat
    oneminusf0  <- 1 - prob0
    df0.dmunb   <- -tempk * prob0
    df0.dkmat   <- prob0 * (tempm + log(tempk))
    df02.dmunb2 <- prob0 * tempk / (kmat + munb) - tempk * df0.dmunb
    df02.dkmat2 <- (prob0 / kmat) * tempm^2
    df02.dkmat.dmunb <- prob0 * (-tempk) * (tempm + log(tempk)) -
                        tempm * prob0 / (kmat + munb)


    dl.dmunb <- y / munb - (y + kmat) / (munb + kmat) +
                df0.dmunb / oneminusf0
    dl.dsize <- digamma(y + kmat) - digamma(kmat) -
                (y + kmat)/(munb + kmat) + 1 + log(tempk) +
                df0.dkmat / oneminusf0

    myderiv <- c(w) * cbind(dl.dmunb * dmunb.deta,
                            dl.dsize * dsize.deta)
    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize ))),
  weight = eval(substitute(expression({
    run.varcov =
    wz <- matrix(0.0, n, 2 * Musual * NOS - 1)




    if (FALSE) {
    usualmeanY <-  munb
    meanY <- usualmeanY / oneminusf0
    ed2l.dmu2 <- meanY / munb^2 -
                (meanY + kmat) / (munb + kmat)^2 -
                df02.dmunb2 / oneminusf0 -
                (df0.dmunb / oneminusf0)^2
    }





    {
      ind2 <- iam(NA, NA, M = Musual, both = TRUE, diag = TRUE)
      for(ii in 1:( .nsimEIM )) {
        ysim <- rposnegbin(n = n*NOS, mu = c(munb), size = c(kmat))
        dim(ysim) <- c(n, NOS)

        dl.dmunb <- ysim / munb - (ysim + kmat) / (munb + kmat) +
                    df0.dmunb / oneminusf0
        dl.dsize <- digamma(ysim + kmat) - digamma(kmat) -
                    (ysim + kmat) / (munb + kmat) + 1 + log(tempk) +
                    df0.dkmat / oneminusf0

        for(kk in 1:NOS) {
          temp2 <- cbind(dl.dmunb[, kk],
                         dl.dsize[, kk]) *
                   cbind(dmunb.deta[, kk],
                         dsize.deta[, kk])
          small.varcov <- temp2[, ind2$row.index] *
                          temp2[, ind2$col.index]

          run.varcov[, ((kk-1)*Musual+1):(kk*Musual)] =
          run.varcov[, ((kk-1)*Musual+1):(kk*Musual)] +
            c(small.varcov[, 1:Musual])
          run.varcov[, M + (kk-1)*Musual + 1] =
          run.varcov[, M + (kk-1)*Musual + 1] +
            c(small.varcov[, Musual + 1])
        }
      } # ii

      run.varcov <- cbind(run.varcov / .nsimEIM )
      wz <- if (intercept.only)
          matrix(colMeans(run.varcov),
                 n, ncol(run.varcov), byrow = TRUE) else run.varcov

    }

    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / Musual)
  }), list( .nsimEIM = nsimEIM ))))
}





dposgeom <- function(x, prob, log = FALSE) {
  dgeom(x - 1, prob = prob, log = log)
}


pposgeom <- function(q, prob) {
  if (!is.Numeric(prob, positive = TRUE))
    stop("bad input for argument 'prob'")
  L <- max(length(q), length(prob))
  if (length(q)    != L) q    = rep(q,    length.out = L);
  if (length(prob) != L) prob = rep(prob, length.out = L);
  ifelse(q < 1, 0,
        (pgeom(q, prob) -
         dgeom(0, prob))
       / pgeom(0, prob, lower.tail = FALSE))
}


qposgeom <- function(p, prob) {




  ans <- qgeom(pgeom(0, prob, lower.tail = FALSE) * p +
               dgeom(0, prob),
               prob = prob)
  ans[p >  1] <- NaN
  ans[p <  0] <- NaN
  ans[p == 1] <- Inf
  ans
}




rposgeom <- function(n, prob) {
  qgeom(p = runif(n, min = dgeom(0, prob)), prob)
}









dpospois <- function(x, lambda, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if (!is.Numeric(lambda, positive = TRUE))
    stop("bad input for argument 'lambda'")
  L <- max(length(x), length(lambda))
  x <- rep(x, len = L); lambda <- rep(lambda, len = L); 

  ans <- if (log.arg) {
    ifelse(x == 0, log(0.0), dpois(x, lambda, log = TRUE) -
           log1p(-exp(-lambda)))
  } else {
    ifelse(x == 0, 0, -dpois(x, lambda) / expm1(-lambda))
  }
  ans
}


ppospois <- function(q, lambda) {
  if (!is.Numeric(lambda, positive = TRUE))
    stop("bad input for argument 'lambda'")
  L <- max(length(q), length(lambda))
  if (length(q)      != L) q      <- rep(q,      length.out = L);
  if (length(lambda) != L) lambda <- rep(lambda, length.out = L);

  ifelse(q < 1, 0,
        (ppois(q, lambda) -
         dpois(0, lambda))
       / ppois(0, lambda, lower.tail = FALSE))
}


qpospois <- function(p, lambda) {


  ans <- qpois(ppois(0, lambda, lower.tail = FALSE) * p +
               dpois(0, lambda),
               lambda = lambda)

  ans[p >  1] <- NaN
  ans[p <  0] <- NaN
  ans[p == 1] <- Inf
  ans
}




rpospois <- function(n, lambda) {
  qpois(p = runif(n, min = dpois(0, lambda)), lambda)
}



rposnegbin <- function(n, size, prob = NULL, munb = NULL) {
  if (!is.null(munb)) {
    if (!is.null(prob))
        stop("'prob' and 'mu' both specified")
    qnbinom(p = runif(n,
                      min = dnbinom(0, size,              mu = munb)),
            size,              mu = munb)
  } else {
    qnbinom(p = runif(n,
                      min = dnbinom(0, size, prob = prob           )),
            size, prob = prob           )
  }
}




 pospoisson <- function(link = "loge", expected = TRUE,
                        ilambda = NULL, imethod = 1, zero = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.logical(expected) || length(expected) != 1)
    stop("bad input for argument 'expected'")
  if (length( ilambda) && !is.Numeric(ilambda, positive = TRUE))
    stop("bad input for argument 'ilambda'")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
    imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")

  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE))
    stop("bad input for argument 'zero'")



  new("vglmff",
  blurb = c("Positive-Poisson distribution\n\n",
            "Links:    ",
            namesof("lambda", link, earg = earg, tag = FALSE)),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 1
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 1,
         link = .link ,
         earg = .earg)
  }, list( .link = link, .earg = earg ))),

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
    Musual <- 1
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly



    mynames1 <- paste("lambda",
                      if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
      namesof(mynames1, .link , earg = .earg, tag = FALSE)

    if ( .imethod == 1) {
      lambda.init <- apply(y, 2, median) + 1/8
      lambda.init <- matrix(lambda.init, n, ncoly, byrow = TRUE)
    } else if ( .imethod == 2) {
      lambda.init <- apply(y, 2, weighted.mean, w = w) + 1/8
      lambda.init <- matrix(lambda.init, n, ncoly, byrow = TRUE)
    } else {
      lambda.init <- -y / expm1(-y)
    }
    if (length( .ilambda))
      lambda.init <- lambda.init * 0 + .ilambda

    if (!length(etastart))
      etastart <- theta2eta(lambda.init, .link , earg = .earg)
  }), list( .link = link, .earg = earg,
            .ilambda = ilambda, .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    lambda <- eta2theta(eta, .link , earg = .earg )
    -lambda / expm1(-lambda)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <- rep( .link , len = M)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for(ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .link = link, .earg = earg, .expected = expected ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    lambda <- eta2theta(eta, .link , earg = .earg ) 
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      sum(w * dpospois(x = y, lambda = lambda, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("pospoisson"),
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta, .link , earg = .earg ) 

    temp6 <- expm1(lambda)
    dl.dlambda <- y / lambda - 1 - 1 / temp6

    dlambda.deta <- dtheta.deta(lambda, .link , earg = .earg )

    c(w) * dl.dlambda * dlambda.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    if ( .expected ) {
      ned2l.dlambda2 <- (temp6 + 1) * (1/lambda - 1/temp6) / temp6
      wz <-  ned2l.dlambda2 * dlambda.deta^2
    } else {
      d2l.dlambda2 <- y / lambda^2 - (temp6 + 1) / temp6^2
      d2lambda.deta2 <- d2theta.deta2(lambda, .link , earg = .earg)
      wz <- (dlambda.deta^2) * d2l.dlambda2 - dl.dlambda * d2lambda.deta2
    }
    c(w) * wz
  }), list( .link = link, .earg = earg, .expected = expected ))))
}








pposbinom <- function(q, size, prob 
                     ) {


  if (!is.Numeric(prob, positive = TRUE)) 
    stop("no zero or non-numeric values allowed for argument 'prob'")
  L <- max(length(q), length(size), length(prob))
  if (length(q)      != L) q      <- rep(q,      length.out = L);
  if (length(size)   != L) size   <- rep(size,   length.out = L);
  if (length(prob)   != L) prob   <- rep(prob,   length.out = L);

  ifelse(q < 1, 0,
        (pbinom(q = q, size = size, prob = prob) -
         dbinom(x = 0, size = size, prob = prob))
       / pbinom(q = 0, size = size, prob = prob, lower.tail = FALSE))
}


qposbinom <- function(p, size, prob
                     ) {




  ans <- qbinom(pbinom(0, size, prob, lower.tail = FALSE) * p +
                dbinom(0, size, prob),
                size = size, prob = prob)

  ans[p >  1] <- NaN
  ans[p <  0] <- NaN
  ans[p == 1] <- size[p == 1]
  ans
}



rposbinom <- function(n, size, prob) {
  qbinom(p = runif(n, min = dbinom(0, size, prob)), size, prob)
}



dposbinom <- function(x, size, prob, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  L <- max(length(x), length(size), length(prob))
  x    <- rep(x,    len = L);
  size <- rep(size, len = L);
  prob <- rep(prob, len = L);

  answer <- NaN * x
  is0 <- (x == 0)
  ok2 <- (prob > 0) & (prob <= 1) &
         (size == round(size)) & (size > 0)

  answer <-        dbinom(x = x, size = size, prob = prob, log = TRUE) -
            log1p(-dbinom(x = 0, size = size, prob = prob))
  answer[!ok2] <- NaN
  if (log.arg) {
    answer[is0 & ok2] <- log(0.0)
  } else {
    answer <- exp(answer)
    answer[is0 & ok2] <- 0.0
  }
  answer
}







 posbinomial <-
  function(link = "logit",
           mv = FALSE, parallel = FALSE, zero = NULL) {


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")



  if (!is.logical(mv) || length(mv) != 1)
    stop("bad input for argument 'mv'")

  if (mv && length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE))
    stop("bad input for argument 'zero'")


  new("vglmff",
  blurb = c("Positive-binomial distribution\n\n",
            "Links:    ",
            if (mv)
            c(namesof("prob1", link, earg = earg, tag = FALSE),
              ",...,",
              namesof("probM", link, earg = earg, tag = FALSE)) else
            namesof("prob", link, earg = earg, tag = FALSE),
            "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(1, M, 1), x, .parallel , constraints)

    dotzero <- .zero
    Musual <- 1
    eval(negzero.expression)
  }), list( .parallel = parallel, .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(Musual = 1,
         zero = .zero)
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    mustart.orig <- mustart
    if ( .mv ) {
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

      extra$orig.w <- w
      mustart <- matrix(colSums(y) / colSums(w), # Not colSums(y * w)...
                        n, ncoly, byrow = TRUE)

    } else {
      eval(binomialff(link = .earg , # earg = .earg ,
                      earg.link = TRUE)@initialize)
    }


    if ( .mv ) {

      dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
      dn2 <- if (length(dn2)) {
        paste("E[", dn2, "]", sep = "")
      } else {
        paste("prob", 1:M, sep = "")
      }
      predictors.names <- namesof(if (M > 1) dn2 else
        "prob", .link , earg = .earg, short = TRUE)

      w <- matrix(w, n, ncoly)
      y <- y / w # Now sample proportion
    } else {
      predictors.names <-
        namesof("prob", .link , earg = .earg , tag = FALSE)
    }

    if (length(extra)) extra$w <- w else extra <- list(w = w)

    if (!length(etastart)) {
      mustart.use <- if (length(mustart.orig)) mustart.orig else mustart
      etastart <- cbind(theta2eta(mustart.use, .link , earg = .earg ))
    }
    mustart <- NULL
  }), list( .link = link,
            .earg = earg, .mv = mv ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    w <- extra$w
    binprob <- eta2theta(eta, .link , earg = .earg )
    nvec <- if ( .mv ) {
             w
           } else {
             if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
               round(w)
           }
    binprob / (1.0 - (1.0 - binprob)^nvec)
  },

  list( .link = link, .earg = earg, .mv = mv ))),
  last = eval(substitute(expression({
    extra$w   <- NULL # Kill it off 


    misc$link <- rep( .link , length = M)
    names(misc$link) <- if (M > 1) dn2 else "prob"

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$expected <- TRUE

    misc$mv   <- .mv
    w <- as.numeric(w)
  }), list( .link = link, .earg = earg, .mv = mv ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {

      ycounts <- if ( .mv ) {
                  round(y * extra$orig.w)
                } else {
                  if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                  y * w # Convert proportions to counts
                }
      nvec <- if ( .mv ) {
               w
             } else {
               if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                 round(w)
             }
      use.orig.w <- if (is.numeric(extra$orig.w)) extra$orig.w else 1
    binprob <- eta2theta(eta, .link , earg = .earg )

    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {
      sum(use.orig.w * dposbinom(x = ycounts, size = nvec,
                                 prob = binprob, log = TRUE))
    }
  }, list( .link = link, .earg = earg, .mv = mv ))),

  vfamily = c("posbinomial"),
  deriv = eval(substitute(expression({
    use.orig.w <- if (is.numeric(extra$orig.w)) extra$orig.w else
                  rep(1, n)

    nvec <- if ( .mv ) {
              w
            } else {
              if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)
            }
    binprob <- eta2theta(eta, .link , earg = .earg )
    dmu.deta <- dtheta.deta(binprob, .link , earg = .earg )

    temp1 <- 1 - (1 - binprob)^nvec
    temp2 <-     (1 - binprob)^2
    temp3 <-     (1 - binprob)^(nvec-2)

    dl.dmu <- y / binprob - (1 - y) / (1 - binprob) -
             (1 - binprob) * temp3 / temp1

    c(w) * dl.dmu * dmu.deta
  }), list( .link = link, .earg = earg, .mv = mv ))),
  weight = eval(substitute(expression({

    ned2l.dmu2 <- 1 / (binprob * temp1) +
                  (1 - mu) / temp2 -
                  (nvec-1) * temp3 / temp1 -
                  nvec * (temp2^(nvec-1)) / temp1^2



    wz <- c(w) * ned2l.dmu2 * dmu.deta^2
    wz
  }), list( .link = link, .earg = earg, .mv = mv ))))
}







 posbernoulli.t <-
  function(link = "logit",
           parallel.t = FALSE,
           apply.parint = TRUE,
           iprob = NULL) {






  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(iprob))
  if (!is.Numeric(iprob, positive = TRUE) ||
        max(iprob) >= 1)
    stop("argument 'iprob' must have values in (0, 1)")

  if (!is.logical(apply.parint) ||
      length(apply.parint) != 1)
    stop("argument 'apply.parint' must be a single logical")


  new("vglmff",
  blurb = c("(Multiple) positive-Bernoulli (capture-recapture) model ",
            "with temporal effects (M_t)\n\n",
            "Links:    ",
            namesof("prob1", link, earg = earg, tag = FALSE), ", ",
            namesof("prob2", link, earg = earg, tag = FALSE), ", ..., ",
            namesof("probM", link, earg = earg, tag = FALSE),
            "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(1, M, 1), x, .parallel.t , constraints,
                           apply.int = .apply.parint , #  TRUE,
                           cm.default = diag(M),
                           cm.intercept.default = diag(M))
  }), list( .parallel.t = parallel.t,
            .apply.parint = apply.parint ))),
  infos = eval(substitute(function(...) {
    list(Musual = 1,
         multipleResponses = TRUE,
         apply.parint = .apply.parint ,
         parallel.t = .parallel.t )
  }, list( .parallel.t = parallel.t,
           .apply.parint = apply.parint ))),

  initialize = eval(substitute(expression({
    Musual <- 1

    mustart.orig <- mustart
    y <- as.matrix(y)
    M <- ncoly <- ncol(y)
    extra$tau <- tau <- ncol(y)
    extra$orig.w <- w

    w <- matrix(w, n, ncoly)
    mustart <- matrix(colSums(y) / colSums(w),
                    n, ncol(y), byrow = TRUE)
    mustart[mustart == 0] <- 0.05
    mustart[mustart == 1] <- 0.95

    if (ncoly == 1)
      stop("the response is univariate, therefore use posbinomial()")






    if (!all(y == 0 | y == 1))
      stop("response must contain 0s and 1s only")
    if (!all(w == 1))
      stop("argument 'weight' must contain 1s only")



    dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
    dn2 <- if (length(dn2)) {
      paste("E[", dn2, "]", sep = "")
    } else {
      paste("prob", 1:M, sep = "")
    }


    predictors.names <- namesof(dn2, .link , earg = .earg, short = TRUE)


    if (length(extra)) extra$w <- w else extra <- list(w = w)

    if (!length(etastart)) {
      mustart.use <- if (length(mustart.orig)) {
        mustart.orig
      } else {
        mustart
      }
      etastart <- cbind(theta2eta(mustart.use, .link , earg = .earg ))
    }
    mustart <- NULL
  }), list( .link = link, .earg = earg ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    probs <- eta2theta(eta, .link , earg = .earg )
    logAA0 <- rowSums(log1p(-probs))
    AA0 <- exp(logAA0)
    AAA <- exp(log1p(-AA0))  # 1 - AA0
    probs / AAA
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    extra$w   <- NULL   # Kill it off 


    misc$link <- rep( .link , length = M)
    names(misc$link) <- if (M > 1) dn2 else "prob"

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:M) misc$earg[[ii]] <- .earg


    misc$mv           <- TRUE
    misc$iprob        <- .iprob



    R <- tfit$qr$qr[1:ncol_X_vlm, 1:ncol_X_vlm, drop = FALSE]
    R[lower.tri(R)] <- 0
    tmp6 <- N.hat.posbernoulli(eta = eta, link = .link , earg = .earg ,
                               R = R, w = w,
                               X_vlm = X_vlm_save, Hlist = constraints,
                               extra = extra, model.type = "t")
    extra$N.hat    <- tmp6$N.hat
    extra$SE.N.hat <- tmp6$SE.N.hat




    misc$parallel.t   <- .parallel.t
    misc$apply.parint <- .apply.parint
  }), list( .link = link, .earg = earg,
            .parallel.t = parallel.t,
            .apply.parint = apply.parint,
            .iprob = iprob ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {

    ycounts <- y
    use.orig.w <- if (length(extra$orig.w)) extra$orig.w else 1

    probs <- eta2theta(eta, .link , earg = .earg )

    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {

      sum(dposbern(x = ycounts, # size = 1, # Bernoulli trials
                   prob = probs, prob0 = probs, log = TRUE))


      sum(use.orig.w *
          dposbern(x = ycounts, # size = 1, # Bernoulli trials
                   prob = probs, prob0 = probs, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("posbernoulli.t"),
  deriv = eval(substitute(expression({

 

    probs <- eta2theta(eta, .link , earg = .earg )
    dprobs.deta <- dtheta.deta(probs, .link , earg = .earg )

    logAA0 <- rowSums(log1p(-probs))
    AA0 <- exp(logAA0)
    AAA <- exp(log1p(-AA0))  # 1 - AA0

    B_s <- AA0 / (1 - probs)
    B_st <- array(AA0, c(n, M, M))
    for(slocal in 1:(M-1))
      for(tlocal in (slocal+1):M)
        B_st[, slocal, tlocal] =
        B_st[, tlocal, slocal] <- B_s[, slocal] / (1 - probs[, tlocal])




    temp2 <-     (1 - probs)^2

    dl.dprobs <- y / probs - (1 - y) / (1 - probs) - B_s / AAA

    deriv.ans <- w * dl.dprobs * dprobs.deta
    deriv.ans
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({

    ed2l.dprobs2 <- 1 / (probs * AAA) + 1 / temp2 -
                probs / (AAA * temp2) - (B_s / AAA)^2

    wz <- matrix(as.numeric(NA), n, dimm(M))
    wz[, 1:M] <- ed2l.dprobs2 * (dprobs.deta^2)

    for(slocal in 1:(M-1))
      for(tlocal in (slocal+1):M)
        wz[, iam(slocal, tlocal, M = M)] <- dprobs.deta[, slocal] *
                                            dprobs.deta[, tlocal] *
                                            (B_st[,slocal,tlocal] +
                                             B_s [,slocal] *
                                             B_s [,tlocal] / AAA) / (-AAA)



    wz
  }), list( .link = link, .earg = earg ))))
}





 posbernoulli.b <-
  function(link = "logit",
           parallel.b = FALSE,  # TRUE,
           apply.parint = TRUE,
           icap.prob = NULL,
           irecap.prob = NULL
          ) {




  fit.type <- 1  # Currently only this is implemented

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(icap.prob))
  if (!is.Numeric(icap.prob, positive = TRUE) ||
        max(icap.prob) >= 1)
    stop("argument 'icap.prob' must have values in (0, 1)")
  if (length(irecap.prob))
  if (!is.Numeric(irecap.prob, positive = TRUE) ||
        max(irecap.prob) >= 1)
    stop("argument 'irecap.prob' must have values in (0, 1)")

  if (!is.logical(parallel.b) ||
      length(parallel.b) != 1)
    stop("argument 'parallel.b' must be a single logical")


  new("vglmff",
  blurb = c("(Multiple) positive-Bernoulli (capture-recapture) model ",
            "with behavioural effects (M_b)\n\n",
            "Links:    ",
            namesof("cap.prob",   link, earg = earg, tag = FALSE), ", ",
            namesof("recap.prob", link, earg = earg, tag = FALSE),
            "\n"),

  constraints = eval(substitute(expression({

    constraints <- cm.vgam(matrix(1, 2, 1), x = x,
                           bool = .parallel.b ,
                           constraints = constraints,
                           apply.int = .apply.parint ,  # TRUE, 
                           cm.default = matrix(1, 2, 1),
                           cm.intercept.default = cbind(1, 0:1))
  }), list( .parallel.b = parallel.b,
            .apply.parint = apply.parint ))),

  infos = eval(substitute(function(...) {
    list( Musual = 2,
         apply.parint = .apply.parint ,
         multipleResponses = FALSE)
  }, list(
           .apply.parint = apply.parint
         ))),

  initialize = eval(substitute(expression({
    Musual <- 2
    if (!is.matrix(y) || ncol(y) == 1)
      stop("the response appears to be univariate")

    if (!all(y == 0 | y == 1))
      stop("response must contain 0s and 1s only")

    orig.y <- y
    extra$orig.w <- w
    extra$tau <- tau <- ncol(y)
    mustart.orig <- mustart
    M <- 2


    tmp3 <- aux.posbernoulli(y)
    y0i        <- extra$y0i  <-       tmp3$y0i
    yr0i       <- extra$yr0i <-       tmp3$yr0i
    yr1i       <- extra$yr1i <-       tmp3$yr1i
    cap1       <- extra$cap1 <-       tmp3$cap1
    cap.hist1  <- extra$cap.hist1  <- tmp3$cap.hist1


    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.min = 2,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = ncol(y),
              maximize = TRUE)
    w <- temp5$w  # Retain the 0-1 response
    y <- temp5$y  # Retain the 0-1 response

    mustart <- matrix(colMeans(y), n, tau, byrow = TRUE)
    mustart <- (mustart + orig.y) / 2




    predictors.names <-
      c(namesof(  "cap.prob",  .link , earg = .earg, short = TRUE),
        namesof("recap.prob",  .link , earg = .earg, short = TRUE))

    if (tau >= 4) {
      pbd <- posbern.aux(tau = tau)
    }

    if (!length(etastart)) {
      mustart.use <- if (length(mustart.orig)) {
        mustart.orig
      } else {
        mustart
      }

      etastart <-
        cbind(theta2eta(rowMeans(mustart.use), .link , earg = .earg ),
              theta2eta(rowMeans(mustart.use), .link , earg = .earg ))

      if (length(   .icap.prob ))
        etastart[, 1] <- theta2eta(   .icap.prob , .link , earg = .earg )
      if (length( .irecap.prob ))
        etastart[, 2] <- theta2eta( .irecap.prob , .link , earg = .earg )
    }
    mustart <- NULL
  }), list( .link = link, .earg = earg,
              .icap.prob =   icap.prob,
            .irecap.prob = irecap.prob
          ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    cap.probs <- eta2theta(eta[, 1], .link , earg = .earg )
    rec.probs <- eta2theta(eta[, 2], .link , earg = .earg )
    cap.probs <- matrix(cap.probs, nrow(eta), extra$tau)
    rec.probs <- matrix(rec.probs, nrow(eta), extra$tau)
    tau <- extra$tau

    if ( .fit.type == 1) {
      fv <- rec.probs
      mat.index <- cbind(1:nrow(fv), extra$cap1)
      fv[mat.index] <- cap.probs[mat.index]
      fv[extra$cap.hist1 == 0] <- cap.probs[extra$cap.hist1 == 0]
    } else if ( .fit.type == 2) {
      fv <- cap.probs
    } else if ( .fit.type == 3) {
      fv <- rec.probs
    } else if ( .fit.type == 4) {
      stop("argument 'fit.type' unmatched")
    } else {
      stop("argument 'fit.type' unmatched")
    }
    fv
  }, list( .link = link,
           .fit.type = fit.type,
           .earg = earg ))),
  last = eval(substitute(expression({

    misc$link <- c( .link , .link )
    names(misc$link) <- predictors.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- .earg
    misc$earg[[2]] <- .earg

    misc$expected    <- TRUE
    misc$mv          <- TRUE
    misc$icap.prob   <- .icap.prob
    misc$irecap.prob <- .irecap.prob
    misc$parallel.b  <- .parallel.b
    misc$fit.type    <- .fit.type
    misc$multipleResponses <- FALSE
    if (tau >= 4) {
      misc$pbd       <- pbd  # Needed for vcov() post-analysis.
    }
    misc$apply.parint <- .apply.parint



    R <- tfit$qr$qr[1:ncol_X_vlm, 1:ncol_X_vlm, drop = FALSE]
    R[lower.tri(R)] <- 0
    tmp6 <- N.hat.posbernoulli(eta = eta, link = .link , earg = .earg ,
                               R = R, w = w,
                               X_vlm = X_vlm_save, Hlist = constraints,
                               extra = extra, model.type = "b")
    extra$N.hat    <- tmp6$N.hat
    extra$SE.N.hat <- tmp6$SE.N.hat


  }), list( .link = link, .earg = earg,
            .fit.type = fit.type,
            .parallel.b = parallel.b,
            .icap.prob =   icap.prob,
            .irecap.prob = irecap.prob,
            .apply.parint = apply.parint
          ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {

    ycounts <- y
    use.orig.w <- if (length(extra$orig.w)) extra$orig.w else 1

    cap.probs <- eta2theta(eta[, 1], .link , earg = .earg )
    rec.probs <- eta2theta(eta[, 2], .link , earg = .earg )
    cap.probs <- matrix(cap.probs, nrow(eta), extra$tau)
    rec.probs <- matrix(rec.probs, nrow(eta), extra$tau)

    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {
      sum(use.orig.w *
          dposbern(x = ycounts,  # Bernoulli trials
                   prob = mu, prob0 = cap.probs, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("posbernoulli.b"),
  deriv = eval(substitute(expression({
    cap.probs <- eta2theta(eta[, 1], .link , earg = .earg )
    rec.probs <- eta2theta(eta[, 2], .link , earg = .earg )
    y0i  <- extra$y0i
    yr0i <- extra$yr0i
    yr1i <- extra$yr1i
    cap1 <- extra$cap1
    tau  <- extra$tau

    dcapprobs.deta <- dtheta.deta(cap.probs, .link , earg = .earg )
    drecprobs.deta <- dtheta.deta(rec.probs, .link , earg = .earg )

    QQQ <- (1 - cap.probs)^tau
    dl.dcap <-   1  /      cap.probs -
               y0i  / (1 - cap.probs) -
               tau * ((1 - cap.probs)^(tau - 1)) / (1 - QQQ)

    dl.drec <- yr1i /      rec.probs -
               yr0i / (1 - rec.probs)


    deriv.ans <- c(w) * cbind(dl.dcap * dcapprobs.deta,
                              dl.drec * drecprobs.deta)
    deriv.ans
  }), list( .link = link, .earg = earg ))),

  weight = eval(substitute(expression({

    wz <- matrix(0, n, M) # Diagonal EIM


    if (tau == 2)
      wz[, iam(2, 2, M = M)] <- (cap.probs / (rec.probs * (1 - rec.probs) *
                                 (1 - QQQ))) * drecprobs.deta^2
    if (tau == 3)
      wz[, iam(2, 2, M = M)] <- (cap.probs * (3 - cap.probs) / (
                                 rec.probs * (1 - rec.probs) *
                                 (1 - QQQ))) * drecprobs.deta^2


    if (tau >= 4) {
                                   # rec.probs = rec.probs)
      eim.rec.tot <- 0
      for (ii in 1:nrow(pbd$part1.rec)) {
        if (pbd$ml..konst.rec[ii, 1] != 0)
          eim.rec.tot <- eim.rec.tot +
          pbd$ml..konst.rec[ii, 1] * ((  cap.probs)^pbd$part1.rec[ii, 1] *
                                      (1-cap.probs)^pbd$part1.rec[ii, 2] *
                                      (  rec.probs)^pbd$part1.rec[ii, 3] *
                                      (1-rec.probs)^pbd$part1.rec[ii, 4])
        if (pbd$ml..konst.rec[ii, 2] != 0)
          eim.rec.tot <- eim.rec.tot +
          pbd$ml..konst.rec[ii, 2] * ((  cap.probs)^pbd$part2.rec[ii, 1] *
                                      (1-cap.probs)^pbd$part2.rec[ii, 2] *
                                      (  rec.probs)^pbd$part2.rec[ii, 3] *
                                      (1-rec.probs)^pbd$part2.rec[ii, 4])
        if (pbd$ml..konst.rec[ii, 3] != 0)
          eim.rec.tot <- eim.rec.tot +
          pbd$ml..konst.rec[ii, 3] * ((  cap.probs)^pbd$part3.rec[ii, 1] *
                                      (1-cap.probs)^pbd$part3.rec[ii, 2] *
                                      (  rec.probs)^pbd$part3.rec[ii, 3] *
                                      (1-rec.probs)^pbd$part3.rec[ii, 4])
        if (pbd$ml..konst.rec[ii, 4] != 0)
          eim.rec.tot <- eim.rec.tot +
          pbd$ml..konst.rec[ii, 4] * ((  cap.probs)^pbd$part4.rec[ii, 1] *
                                      (1-cap.probs)^pbd$part4.rec[ii, 2] *
                                      (  rec.probs)^pbd$part4.rec[ii, 3] *
                                      (1-rec.probs)^pbd$part4.rec[ii, 4])
      }
      eim.rec.tot <- (eim.rec.tot / (1 - QQQ)) * drecprobs.deta^2
      wz[, iam(2, 2, M = M)] <- eim.rec.tot
    }





    dA.dcapprobs <- -tau * ((1 - QQQ) * (tau-1) * (1 - cap.probs)^(tau-2) +
                            tau * (1 - cap.probs)^(2*tau -2)) / (1 - QQQ)^2

    if (tau == 2)
      wz[, iam(1, 1, M = M)] <-
        ((2 - 3 * cap.probs + 2 * cap.probs^2) / ((1 - QQQ) *
        cap.probs * (1 - cap.probs)) + dA.dcapprobs) *
        dcapprobs.deta^2
    if (tau == 3)
      wz[, iam(1, 1, M = M)] <-
        ((3 + cap.probs * (-6 + cap.probs * (7 + cap.probs * (-3)))) / (
         (1 - QQQ) * cap.probs * (1 - cap.probs)) + dA.dcapprobs) *
        dcapprobs.deta^2


    if (tau >= 4) {

      eim.cap.tot <- 0
      for (ii in 1:nrow(pbd$part1.cap)) {
        if (pbd$ml..konst.cap[ii, 1] != 0)
          eim.cap.tot <- eim.cap.tot +
          pbd$ml..konst.cap[ii, 1] * ((  cap.probs)^pbd$part1.cap[ii, 1] *
                                      (1-cap.probs)^pbd$part1.cap[ii, 2] *
                                      (  rec.probs)^pbd$part1.cap[ii, 3] *
                                      (1-rec.probs)^pbd$part1.cap[ii, 4])
        if (pbd$ml..konst.cap[ii, 2] != 0)
          eim.cap.tot <- eim.cap.tot +
          pbd$ml..konst.cap[ii, 2] * ((  cap.probs)^pbd$part2.cap[ii, 1] *
                                      (1-cap.probs)^pbd$part2.cap[ii, 2] *
                                      (  rec.probs)^pbd$part2.cap[ii, 3] *
                                      (1-rec.probs)^pbd$part2.cap[ii, 4])
        if (pbd$ml..konst.cap[ii, 3] != 0)
          eim.cap.tot <- eim.cap.tot +
          pbd$ml..konst.cap[ii, 3] * ((  cap.probs)^pbd$part3.cap[ii, 1] *
                                      (1-cap.probs)^pbd$part3.cap[ii, 2] *
                                      (  rec.probs)^pbd$part3.cap[ii, 3] *
                                      (1-rec.probs)^pbd$part3.cap[ii, 4])
        if (pbd$ml..konst.cap[ii, 4] != 0)
          eim.cap.tot <- eim.cap.tot +
          pbd$ml..konst.cap[ii, 4] * ((  cap.probs)^pbd$part4.cap[ii, 1] *
                                      (1-cap.probs)^pbd$part4.cap[ii, 2] *
                                      (  rec.probs)^pbd$part4.cap[ii, 3] *
                                      (1-rec.probs)^pbd$part4.cap[ii, 4])
      }
      eim.cap.tot <- (eim.cap.tot / (1 - QQQ) + dA.dcapprobs) *
                     dcapprobs.deta^2
      wz[, iam(1, 1, M = M)] <- eim.cap.tot
    }


    wz <- c(w) * wz
    wz
  }), list( .link = link, .earg = earg ))))
}



posbern.aux <- function(tau) {

  y.all <- matrix(0, 2^tau - 0, tau)
  for (jlocal in 1:tau)
    y.all[, jlocal] <- c(rep(0, len = 2^(tau-jlocal)),
                         rep(1, len = 2^(tau-jlocal)))
  y.all <- y.all[-1, ]

  aux <- aux.posbernoulli(y.all, check.y = FALSE)


  nstar <- nrow(y.all)
    l.power.cap <- matrix(0, nstar, 4)
    l.konst.cap <- matrix(0, nstar, 4)
  ml..power.cap <- matrix(0, nstar, 4)
  ml..konst.cap <- matrix(0, nstar, 4)
    l.power.rec <- matrix(0, nstar, 4)
    l.konst.rec <- matrix(0, nstar, 4)
  ml..power.rec <- matrix(0, nstar, 4)
  ml..konst.rec <- matrix(0, nstar, 4)



  l.power.rec[, 3] <- -1
  l.power.rec[, 4] <- -1
  for (jlocal in 1:tau) {
    l.konst.rec[, 3] <-
    l.konst.rec[, 3] + ifelse(y.all[, jlocal] >  0 & jlocal > aux$cap1, 1, 0)
    l.konst.rec[, 4] <-
    l.konst.rec[, 4] - ifelse(y.all[, jlocal] == 0 & jlocal > aux$cap1, 1, 0)
  }



  ml..power.rec[, 3] <- -2
  ml..power.rec[, 4] <- -2
  ml..konst.rec[, 3] <-  l.konst.rec[, 3]
  ml..konst.rec[, 4] <- -l.konst.rec[, 4]



  mux.mat <- cbind(1, aux$y0i, aux$yr1i, aux$yr0i)
  part1.rec <- mux.mat + cbind(ml..power.rec[, 1], 0, 0, 0)
  part2.rec <- mux.mat + cbind(0, ml..power.rec[, 2], 0, 0)
  part3.rec <- mux.mat + cbind(0, 0, ml..power.rec[, 3], 0)
  part4.rec <- mux.mat + cbind(0, 0, 0, ml..power.rec[, 4])







  l.power.cap[, 1] <-  1
  l.power.cap[, 2] <- -1
  l.konst.cap[, 1] <-  1
  l.konst.cap[, 2] <- -aux$y0i



  ml..power.cap[, 1] <- -2
  ml..power.cap[, 2] <- -2
  ml..konst.cap[, 1] <-  1
  ml..konst.cap[, 2] <-  aux$y0i



  mux.mat <- cbind(1, aux$y0i, aux$yr1i, aux$yr0i)
  part1.cap <- mux.mat + cbind(ml..power.cap[, 1], 0, 0, 0)
  part2.cap <- mux.mat + cbind(0, ml..power.cap[, 2], 0, 0)
  part3.cap <- mux.mat + cbind(0, 0, ml..power.cap[, 3], 0)
  part4.cap <- mux.mat + cbind(0, 0, 0, ml..power.cap[, 4])




  list(   y.all       =  y.all,
          part1.cap   =  part1.cap,
          part2.cap   =  part2.cap,
          part3.cap   =  part3.cap,
          part4.cap   =  part4.cap,

          part1.rec   =  part1.rec,
          part2.rec   =  part2.rec,
          part3.rec   =  part3.rec,
          part4.rec   =  part4.rec,
          l.konst.cap =    l.konst.cap,
          l.power.cap =    l.power.cap,
        ml..konst.cap =  ml..konst.cap,
        ml..power.cap =  ml..power.cap,
          l.konst.rec =    l.konst.rec,
          l.power.rec =    l.power.rec,
        ml..konst.rec =  ml..konst.rec,
        ml..power.rec =  ml..power.rec)
}




 posbernoulli.tb <-
  function(link = "logit",
           parallel.t = FALSE,
           parallel.b = FALSE,
           apply.parint = FALSE,
           imethod = 1,
           iprob = NULL,
           dconst = 0.1,
           dpower = -2) {





  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
          max(iprob) >= 1)
      stop("argument 'iprob' must have values in (0, 1)")

  if (!is.logical(parallel.t) ||
      length(parallel.t) != 1)
    stop("argument 'parallel.t' must be a single logical")

  if (!is.logical(parallel.b) ||
      length(parallel.b) != 1)
    stop("argument 'parallel.b' must be a single logical")

  if (!is.logical(apply.parint) ||
      length(apply.parint) != 1)
    stop("argument 'apply.parint' must be a single logical")


  new("vglmff",
  blurb = c("(Multiple) positive-Bernoulli (capture-recapture) model\n",
            "with temporal and behavioural effects (M_{tb})\n\n",
            "Links:    ",
            namesof("cap.prob.1",     link, earg = earg, tag = FALSE), ", ",
            namesof("cap.prob.2",     link, earg = earg, tag = FALSE), ", ",
            ", ...,\n",
            namesof("cap.prob.tau",   link, earg = earg, tag = FALSE), ", ",
            namesof("recap.prob.2",   link, earg = earg, tag = FALSE),
            ", ...,\n",
            namesof("recap.prob.tau", link, earg = earg, tag = FALSE),
            "\n"),
  constraints = eval(substitute(expression({

    tmp8.mat <- cbind(c(1, rep(0, len = 2*(tau-1))),
                      rbind(rep(0, len = tau-1), diag(tau-1), diag(tau-1)))
    tmp9.mat <- cbind(c(rep(0, len = tau), rep(1, len = tau-1)))

    cmk_tb <- if ( .parallel.t ) matrix(1, M, 1) else tmp8.mat

    cm1_tb <-
      if ( ( .parallel.t ) &&  ( .parallel.b )) matrix(1, M, 1) else
      if ( ( .parallel.t ) && !( .parallel.b )) cbind(1, tmp9.mat) else
      if (!( .parallel.t ) &&  ( .parallel.b )) tmp8.mat else
      if (!( .parallel.t ) && !( .parallel.b )) cbind(tmp8.mat, tmp9.mat)


    constraints <- cm.vgam(cmk_tb, x = x,
                           bool = .parallel.t ,  # Same as .parallel.b
                           constraints = constraints,
                           apply.int = .apply.parint ,  # FALSE,  
                           cm.default = cmk_tb,
                           cm.intercept.default = cm1_tb)

  }), list( .parallel.t = parallel.t,
            .parallel.b = parallel.b,
            .apply.parint = apply.parint ))),
  infos = eval(substitute(function(...) {
    list(Musual = 2,
         multipleResponses = TRUE,
         imethod = .imethod ,
         dconst  = .dconst ,
         dpower  = .dpower ,
         apply.parint = .apply.parint ,
         parallel.t = .parallel.t ,
         parallel.b = .parallel.b )
  }, list( .parallel.t = parallel.t,
           .parallel.b = parallel.b,
           .imethod = imethod,
           .dconst = dconst,
           .dpower = dpower,
           .apply.parint = apply.parint ))),

  initialize = eval(substitute(expression({
    Musual <- 2  # Not quite true



    if (ncol(cbind(w)) > 1)
      stop("variable 'w' should be a vector or one-column matrix")
    w <- c(w)  # Make it a vector

    mustart.orig <- mustart
    y <- as.matrix(y)
    extra$tau     <- tau   <- ncol(y)
    extra$ncoly   <- ncoly <- ncol(y)
    extra$orig.w  <- w
    extra$ycounts <- y
    M <- Musual * tau - 1  # recap.prob.1 is unused


    if (!(ncoly %in% 2:3))
      stop("the response currently must be a two- or three-column matrix")



    mustart <- matrix(c(weighted.mean(y[, 1], w),
                        weighted.mean(y[, 2], w),
                        if (tau == 3) weighted.mean(y[, 3], w) else NULL),
                      n, tau, byrow = TRUE)
    mustart[mustart == 0] <- 0.05
    mustart[mustart == 1] <- 0.95





    if (!all(y == 0 | y == 1))
      stop("response must contain 0s and 1s only")


    tmp3 <- aux.posbernoulli(y)
    cap.hist1  <- extra$cap.hist1  <- tmp3$cap.hist1
    if (tau > 2) {
      yindex <- 4 * y[, 1] + 2 * y[, 2] + 1 * y[, 3]
      if (length(table(yindex)) != 2^tau - 1)
        warning("there should be ", 2^tau - 1, " patterns of 0s and 1s ",
                "in the response matrix. May crash.")

    }


    dn2.cap   <- paste("cap.prob.",   1:ncoly, sep = "")
    dn2.recap <- paste("recap.prob.", 2:ncoly, sep = "")

    predictors.names <- c(
      namesof(dn2.cap,   .link , earg = .earg, short = TRUE),
      namesof(dn2.recap, .link , earg = .earg, short = TRUE))

    if (length(extra)) extra$w <- w else extra <- list(w = w)

    if (!length(etastart)) {
      if ( .imethod == 1) {


        mu.init <- if (length( .iprob ))
                     matrix( .iprob , n, M, byrow = TRUE) else
                   if (length(mustart.orig))
                     matrix(rep(mustart.orig, length = n * M), n, M) else
                     matrix(rep(mustart, length = n * M), n, M)
        etastart <- theta2eta(mu.init, .link , earg = .earg ) # n x M
      } else {
        mu.init <- matrix(runif(n * M), n, M)
        etastart <- theta2eta(mu.init, .link , earg = .earg ) # n x M
      }
    }
    mustart <- NULL
  }), list( .link = link, .earg = earg,
            .iprob = iprob,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    tau <- extra$ncoly
    probs <- eta2theta(eta, .link , earg = .earg )
    prc <- probs[, 1:tau]
    prr <- cbind(0, probs[, (1+tau):ncol(probs)])  # 1st coln ignored

    probs.numer <- cbind(probs[, 1],
                         ifelse(extra$cap.hist1[, 2] == 1, prr[, 2], prc[, 2]))

    if (tau == 3)
      probs.numer <- cbind(probs.numer,
                           ifelse(extra$cap.hist1[, 3] == 1, prr[, 3], prc[, 3]))

    logQQQ <- rowSums(log1p(-prc))
    QQQ <- exp(logQQQ)
    AAA <- exp(log1p(-QQQ))  # 1 - QQQ
    probs.numer / AAA
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    extra$w   <- NULL   # Kill it off 


    misc$link <- rep( .link , length = M)
    names(misc$link) <- c(dn2.cap, dn2.recap)

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:M)
      misc$earg[[ii]] <- .earg


    misc$mv       <- TRUE
    misc$iprob    <- .iprob


    R <- tfit$qr$qr[1:ncol_X_vlm, 1:ncol_X_vlm, drop = FALSE]
    R[lower.tri(R)] <- 0
    tmp6 <- N.hat.posbernoulli(eta = eta, link = .link , earg = .earg ,
                               R = R, w = w,
                               X_vlm = X_vlm_save, Hlist = constraints,
                               extra = extra, model.type = "tb")
    extra$N.hat    <- tmp6$N.hat
    extra$SE.N.hat <- tmp6$SE.N.hat


    misc$parallel.t   <- .parallel.t
    misc$parallel.b   <- .parallel.b


    misc$dconst <- .dconst
    misc$dpower <- .dpower
    misc$working.ridge  <- c(rep(adjustment.posbern_tb, length = tau),
                             rep(0,                     length = tau-1))

    misc$apply.parint <- .apply.parint

  }), list( .link = link, .earg = earg,
            .apply.parint = apply.parint,
            .parallel.t = parallel.t,
            .parallel.b = parallel.b,
            .dconst = dconst,
            .dpower = dpower,
            .iprob = iprob ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {

    tau <- extra$ncoly
    ycounts <- y
    use.orig.w <- if (length(extra$orig.w)) extra$orig.w else 1

    probs <- eta2theta(eta, .link , earg = .earg )
    prc <- probs[, 1:tau]
    prr <- cbind(0, probs[, (1+tau):ncol(probs)])  # 1st coln ignored

    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {

    probs.numer <- cbind(probs[, 1],
                         ifelse(extra$cap.hist1[, 2] == 1, prr[, 2], prc[, 2]))
    if (tau == 3)
      probs.numer <- cbind(probs.numer,
                           ifelse(extra$cap.hist1[, 3] == 1, prr[, 3], prc[, 3]))

      sum(use.orig.w *
          dposbern(x = ycounts, # size = 1, # Bernoulli trials
                   prob = probs.numer, prob0 = prc, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("posbernoulli.tb"),
  deriv = eval(substitute(expression({
    tau <- extra$ncoly
    probs <- eta2theta(eta, .link , earg = .earg )
    prc <- probs[, 1:tau]
    prr <- cbind(0, probs[, (1+tau):ncol(probs)])  # 1st coln ignored

    logQQQ <- rowSums(log1p(-prc))
    QQQ <- exp(logQQQ)

    dprobs.deta <- dtheta.deta(probs, .link , earg = .earg )
    dprc.deta <- dprobs.deta[, 1:tau]
    dprr.deta <- cbind(0, dprobs.deta[, (1+tau):ncol(probs)])  # 1st coln ignored

    dQ.dprc   <- -QQQ / (1 - prc)


    d2Q.dprc <- array(0, c(n, tau, tau))
    for (jay in 1:(tau-1))
      for (kay in (jay+1):tau)
        d2Q.dprc[, jay, kay] <-
        d2Q.dprc[, kay, jay] <-  QQQ / ((1 - prc[, jay]) *
                                        (1 - prc[, kay]))

    if (tau == 2)
    dl.dpr <-  cbind(y[, 1] / prc[, 1] - (1 - y[, 1]) / (1 - prc[, 1]) +
                     dQ.dprc[, 1] / (1 - QQQ),
                     (1 - y[, 1]) *
                    (y[, 2] / prc[, 2] - (1 - y[, 2]) / (1 - prc[, 2])) +
                     dQ.dprc[, 2] / (1 - QQQ),
                          y[, 1]  *
                    (y[, 2] / prr[, 2] - (1 - y[, 2]) / (1 - prr[, 2])))

    if (tau == 3)
    dl.dpr <-  cbind(y[, 1] / prc[, 1] - (1 - y[, 1]) / (1 - prc[, 1]) +
                     dQ.dprc[, 1] / (1 - QQQ),

                     (1 - extra$cap.hist1[, 2]) *  # (1 - y[, 1]) *
                    (y[, 2] / prc[, 2] - (1 - y[, 2]) / (1 - prc[, 2])) +
                     dQ.dprc[, 2] / (1 - QQQ),

                     (1 - extra$cap.hist1[, 3]) *  # (1 - y[, 1]) * (1 - y[, 2]) *
                     y[, 3] / prc[, 3] +
                     dQ.dprc[, 3] / (1 - QQQ),

                     extra$cap.hist1[, 2] *  # y[, 1]  *
                    (y[, 2] / prr[, 2] - (1 - y[, 2]) / (1 - prr[, 2])),

                     extra$cap.hist1[, 3] *
                    (y[, 3] / prr[, 3] - (1 - y[, 3]) / (1 - prr[, 3]))
                    )

    deriv.ans <- c(w) * dl.dpr * dprobs.deta

    deriv.ans
  }), list( .link = link, .earg = earg ))),

  weight = eval(substitute(expression({
    wz <- matrix(0, n, sum(M:(M - (tau - 1))))

    cindex <- iam(NA, NA, M = M, both = TRUE)
    cindex$row.index <- rep(cindex$row.index, length = ncol(wz))
    cindex$col.index <- rep(cindex$col.index, length = ncol(wz))


    if (tau == 2) {
      wz[, iam(1, 1, M = M)] <-
               (1 - prc[, 1] * (1 - prc[, 2])) / (prc[, 1] * (1 - prc[, 1]) *
               (1 - QQQ)) -
              ((1 - prc[, 2]) / (1 - QQQ))^2
      wz[, iam(1, 1, M = M)] <- wz[, iam(1, 1, M = M)] * dprc.deta[, 1]^2

      wz[, iam(2, 2, M = M)] <- 
              (prc[, 1] * (1 - prc[, 1]) / (prc[, 2] * (1 - QQQ)^2)) *
               dprc.deta[, 2]^2

      wz[, iam(3, 3, M = M)] <-
              (prc[, 1] / (prr[, 2] * (1 - prr[, 2]) * (1 - QQQ))) *
               dprr.deta[, 2]^2
  
      wz[, iam(1, 2, M = M)] <- -dprc.deta[, 1] * dprc.deta[, 2] / (1 - QQQ)^2
    } else if (tau == 3) {

      wz[, iam(1, 1, M = M)] <-
        ((1 - prc[, 2]) * prc[, 3] + prc[, 2]) / ((1 - prc[, 1]) * (1 - QQQ)) +
         1 / (prc[, 1] * (1 - QQQ)) -
        (dQ.dprc[, 1] / (1 - QQQ))^2


      wz[, iam(2, 2, M = M)] <- 
        (1 - prc[, 1]) * (1 - prc[, 2] * (1 - prc[, 3])) / (
         prc[, 2] * (1 - prc[, 2]) * (1 - QQQ)) -
        (dQ.dprc[, 2] / (1 - QQQ))^2


      wz[, iam(3, 3, M = M)] <-
        (1 - prc[, 1]) * (1 - prc[, 2]) / (prc[, 3] * (1 - QQQ)) -
        (dQ.dprc[, 3] / (1 - QQQ))^2


      wz[, iam(4, 4, M = M)] <-
        prc[, 1] / (prr[, 2] * (1 - prr[, 2]) * (1 - QQQ))
  

      wz[, iam(5, 5, M = M)] <-
        (prc[, 1] + prc[, 2] * (1 - prc[, 1])) / (
         prr[, 3] * (1 - prr[, 3]) * (1 - QQQ))


      for (jay in 1:(tau-1))
        for (kay in (jay+1):tau)
          wz[, iam(jay, kay, M = M)] <-
            -(d2Q.dprc[, jay, kay] +
               dQ.dprc[, jay] *
               dQ.dprc[, kay] / (1 - QQQ)) / (1 - QQQ)


      wz <- wz * dprobs.deta[, cindex$row.index] *
                 dprobs.deta[, cindex$col.index]


    } else {
      stop("tau must equal 2 or 3")
    }


    adjustment.posbern_tb <- .dconst * iter^( .dpower )


     for (jay in 1:tau)
      wz[, iam(jay, jay, M = M)] <- wz[, iam(jay, jay, M = M)] +
                                    adjustment.posbern_tb



    c(w) * wz
  }), list( .link = link, .earg = earg,
            .dconst = dconst,
            .dpower = dpower
          ))))
}










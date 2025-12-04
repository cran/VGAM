# These functions are
# Copyright (C) 1998-2025 T.W. Yee, University of Auckland.
# All rights reserved.











dcard2 <- function(x, mu, rho2 = 0, log = FALSE) {
  if (!isFALSE(log.arg <- log) && !isTRUE(log))
    stop("bad input for argument 'log'")
  rm(log)

  L <- max(length(x), length(mu), length(rho2))
  if (length(x)    < L) x    <- rep_len(x,    L)
  if (length(mu)   < L) mu   <- rep_len(mu,   L)
  if (length(rho2) < L) rho2 <- rep_len(rho2, L)

  ox <- !is.finite(x)
  zero <- ox | x < 0 | x >= 2 * pi
  ans <- rep_len(log(0), L)
  if (any(!zero))
    ans[!zero] <- -log(2 * pi) + log1p(
         rho2[!zero] * cos(x[!zero] - mu[!zero]))
  if (any(ox)) ans[ox] <- log(0)
  ans[mu < 0 | mu >= 2 * pi] <- NaN # Outside PS
  ans[rho2 < 0 | rho2 > 1] <- NaN
  if (log.arg) ans else exp(ans)
}  # dcard2



pcard2 <-
    function(q, mu, rho2 = 0, lower.tail = TRUE,
             log.p = FALSE) {
  if (!isFALSE(lower.tail) && !isTRUE(lower.tail))
    stop("bad input for argument 'lower.tail'")
  if (!isFALSE(log.p) && !isTRUE(log.p))
    stop("bad input for argument 'log.p'")

  ans <- (q + rho2 * (sin(q - mu) + sin(mu))) / (2 * pi)
  ans[q <= 0] <- 0
  ans[q >= 2 * pi] <- 1
        
  ans[mu < 0 | mu >= 2 * pi] <- NaN # Outside PS
  ans[rho2 < 0 | rho2 > 1] <- NaN

  if (lower.tail) {
    if (log.p) log(ans) else ans
  } else { 
    if (log.p) log1p(-ans) else 1 - ans
  }
}  # pcard2





qcard2 <-
  function(p, mu, rho2 = 0,
           lower.tail = TRUE, log.p = FALSE,
           tol0 = 1e-6) {

  if (!isFALSE(log.p) && !isTRUE(log.p))
    stop("bad input for argument 'log.p'")
  if (lower.tail) {
    if (log.p) p <- exp(p)
  } else {  # Unperfect
    p <- if (log.p) -expm1(p) else 1 - p
  }
  if (!is.Numeric(tol0, length.arg = 1,
                  positive = TRUE) | tol0 > 0.001)
    stop("bad argument 'tol0'... too large?")

  L <- max(length(p), length(mu), length(rho2))
  if (length(p)    < L) p    <- rep_len(p,    L)
  if (length(mu)   < L) mu   <- rep_len(mu,   L)
  if (length(rho2) < L) rho2 <- rep_len(rho2, L)

  ans <- p + mu + rho2

  bad0 <- !is.finite(mu) | !is.finite(rho2) | 
          mu < 0 | mu >= 2 * pi | rho2 < 0 | rho2 > 1
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- numeric(L)
  hi <- rep(2 * pi, L)
  dont.iterate <- bad  # | is.easy

  foo <- function(q, mu, rho2, p)
    pcard2(q, mu, rho2) - p

  lhs <- dont.iterate  # | (p <= dfoldnorm(0, mean))

  if (any(!lhs)) {
    ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs],
                      mu = mu[!lhs], rho2 = rho2[!lhs],
                      p = p[!lhs], tol = tol0)
  }


  ans[!bad0 & !is.na(p) & p == 0] <- 0
  ans[!bad0 & !is.na(p) & p == 1] <- 2 * pi
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN

  ans
}  # qcard2




rcard2 <-
    function(n, mu, rho2, tol0 = 1e-4) {
  use.n <- if ((lenn <- length(n)) > 1) lenn else
           if (!is.Numeric(n, integer.valued = TRUE,
                length.arg = 1, positive = TRUE))
              stop("bad input for 'n'") else n
  if (!is.Numeric(tol0, length.arg = 1,
                  positive = TRUE) | tol0 > 0.05)
    stop("bad argument 'tol0'")

  L <- use.n
  if (length(mu)   < L) mu   <- rep_len(mu,   L)
  if (length(rho2) < L) rho2 <- rep_len(rho2, L)

  U <- runif(use.n)
  psi <- runif(use.n, 0, 2 * pi)
  pvec <- (1 + rho2 * cos(psi)) / 2
  omega <- psi  # For U <= pvec really
  ind1 <- U > pvec & psi <= pi
  omega[ind1] <-     pi - psi[ind1]
  ind2 <- U > pvec & psi >  pi
  omega[ind2] <- 3 * pi - psi[ind2]
  ans <- (omega + mu) %% (2 * pi)

  if (any((ind5 <- rho2 < tol0)))
    ans[ind5] <- runif(sum(ind5), 0, 2 * pi)

  ans[mu < 0 | mu >= 2 * pi] <- NaN  # Outside PS
  ans[rho2 < 0 | rho2 > 1] <- NaN

  ans
}  # rcard2






 cardioid2 <- function(
     lmu   = "extlogitlink(min = 0, max = 2*pi)",
     lrho2 = "logitlink",
     imu  = NULL, irho2 = NULL,
     gmu   = ppoints(8) * 2 * pi,
     grho2 = ppoints(8),
     zero = "rho2") {

  if (is.character(lmu))
    lmu <- substitute(y9, list(y9 = lmu))
  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  if (is.character(lrho2))
    lrho2 <- substitute(y9, list(y9 = lrho2))
  lrho2 <- as.list(substitute(lrho2))
  erho2 <- link2list(lrho2)
  lrho2 <- attr(erho2, "function.name")

 if (length(imu) &&
    (!is.Numeric(imu, positive = TRUE) ||
      any(imu > 2 * pi)))
    stop("bad input for argument 'imu'")
  if (length(irho2) &&
     (!is.Numeric(irho2, positive = TRUE) ||
        any(irho2 > 1)))
      stop("bad input for argument 'irho2'")

  F <- FALSE
  new("vglmff",
  blurb = c("Cardioid distribution\n\n",
            "Links:    ",
            namesof("mu",   lmu,     emu, tag = F), ", ",
            namesof("rho2", lrho2, erho2, tag = F), "\n",
            "Mean:     ",
            "pi+(rho2/(2*pi))*((2*pi-mu)*sin(2*pi-mu)",
            "+cos(2*pi-mu)-mu*sin(mu)-cos(mu))"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x,
             .zero , M = M, M1 = 2,
             predictors.names = predictors.names)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2, M = 2,
         Q1 = 1,
         dpqrfun = "card2",
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("mu", "rho2"),
         lmu   = .lmu   ,
         lrho2 = .lrho2 ,
         zero = .zero )
  },
  list( .zero = zero, .lmu = lmu, .lrho2 = lrho2))),

  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    if (n <= 10)
      warning("n may be too small for Fisher scoring")
    if (any(y < 0 | y >= 2 * pi))
      stop("the response must be in [0, 2*pi)")

    gmuuu <- ( .gmu   )
    imuuu <- ( .imu   )
    grho2 <- ( .grho2 )
    irho2 <- ( .irho2 )
    if (length(imuuu)) gmuuu <- imuuu
    if (length(irho2)) grho2 <- irho2
      
    predictors.names <- c(
      namesof("mu",   .lmu   , .emu   , tag = FALSE),
      namesof("rho2", .lrho2 , .erho2 , tag = FALSE))

    if (!length(etastart)) {
      c.ell <- function(mu, rho2, y, x, w, extraargs)
        sum(c(w) * dcard2(y, mu, rho2, log = TRUE))
      try.this <-
        grid.search2(gmuuu, grho2, objfun = c.ell,
                     y = y, w = w,
                     ret.objfun = TRUE)
      muuu.init <- try.this["Value1"]
      rho2.init <- try.this["Value2"]
      etastart <- matrix(byrow = TRUE, nr = n, nc = 2,
        c(theta2eta(muuu.init, .lmu   , .emu   ),
          theta2eta(rho2.init, .lrho2 , .erho2 )))
    }
  }), list( .lmu = lmu, .lrho2 = lrho2,
            .imu = imu, .irho2 = irho2,
            .gmu = gmu, .grho2 = grho2,
            .emu = emu, .erho2 = erho2 ))),
  linkinv = eval(substitute(function(eta, extra = NULL){
    mu   <- eta2theta(eta[, 1], .lmu   , .emu   )
    rho2 <- eta2theta(eta[, 2], .lrho2 , .erho2 )
    pi + (rho2 / (2 * pi)) *
    ((2 * pi - mu) * sin(2 * pi - mu) +
    cos(2 * pi - mu) - mu * sin(mu) - cos(mu))
  }, list( .lmu = lmu, .lrho2 = lrho2,
           .emu = emu, .erho2 = erho2 ))),
  last = eval(substitute(expression({
    misc$link <-    c("mu" = .lmu , "rho2" = .lrho2 )
    misc$earg <- list("mu" = .emu , "rho2" = .erho2 )
  }),
  list( .lmu = lmu, .lrho2 = lrho2,
        .emu = emu, .erho2 = erho2 ))),
  loglikelihood = eval(substitute(
  function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    mu   <- eta2theta(eta[, 1], .lmu   , .emu   )
    rho2 <- eta2theta(eta[, 2], .lrho2 , .erho2 )
    if (residuals) {
      stop("loglikelihood residuals not implemented")
    } else {
      ll.elts <- c(w) * dcard2(y, mu, rho2, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmu = lmu, .lrho2 = lrho2,
           .emu = emu, .erho2 = erho2 ))),
  vfamily = c("cardioid2"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
  mu   <- eta2theta(eta[, 1], .lmu   , .emu   )
  rho2 <- eta2theta(eta[, 2], .lrho2 , .erho2 )
  okay1 <-
  all(is.finite(mu  )) && all(0 <=   mu & mu   < 2*pi) &&
  all(is.finite(rho2)) && all(0 <  rho2 & rho2 < 1)
  okay1
  }, list( .lmu = lmu, .lrho2 = lrho2,
           .emu = emu, .erho2 = erho2 ))),

  simslot = eval(substitute(function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
       pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object); T <- TRUE; F <- FALSE
    mu   <- eta2theta(eta[, c(T, F)], .lmu   , .emu   )
    rho2 <- eta2theta(eta[, c(F, T)], .lrho2 , .erho2 )
    rcard2(nsim * length(mu), mu, rho2)
  }, list( .lmu = lmu, .lrho2 = lrho2,
           .emu = emu, .erho2 = erho2 ))),

  deriv = eval(substitute(expression({
    mu   <- eta2theta(eta[, 1], .lmu   , .emu   )
    rho2 <- eta2theta(eta[, 2], .lrho2 , .erho2 )
    dmu.deta   <- dtheta.deta(mu,   .lmu   , .emu   )
    drho2.deta <- dtheta.deta(rho2, .lrho2 , .erho2 )
    Dnm <-      1 + rho2 * cos(y - mu)
    dl.dmu   <-     rho2 * sin(y - mu) / Dnm
    dl.drho2 <-            cos(y - mu) / Dnm
    c(w) * cbind(dl.dmu  *  dmu.deta,
                 dl.drho2 * drho2.deta)
  }),
  list( .lmu = lmu, .lrho2 = lrho2,
        .emu = emu, .erho2 = erho2 ))),
  weight = eval(substitute(expression({
    tmp5 <- sqrt(abs(1 - rho2^2))
    ned2l.dmu2   <-  1 - tmp5
    ned2l.drho22 <- (1 - tmp5) / (tmp5 * rho2^2)
    wz <- cbind(ned2l.dmu2   *   dmu.deta^2,
                ned2l.drho22 * drho2.deta^2)
    c(w) * wz
  }),
  list( .lmu = lmu, .lrho2 = lrho2))))
}  # cardioid2







dvMF3 <-
  function(x, colatitude, longitude, concentration,
           byrow.arg = FALSE, log = FALSE) {
  if (!isFALSE(log.arg <- log) && !isTRUE(log))
    stop("bad input for argument 'log'")
  rm(log)

  if (is.vector(x) && length(x) == 2)
    x <- rbind(x)
  if (!is.matrix(x))
    x <- as.matrix(x)
  LLL <- max(nrow(x), length(colatitude),
             length(longitude), length(concentration))
  if (nrow(x) < LLL)
    x <- matrix(as.vector(x), LLL, 2, byrow = byrow.arg)
  if (length(longitude) < LLL)
    longitude <- rep_len(longitude, LLL)
  if (length(colatitude) < LLL)
    colatitude <- rep_len(colatitude, LLL)
  if (length(concentration) < LLL)
    concentration <- rep_len(concentration, LLL)
  
  bad0 <- !is.finite(colatitude) | !is.finite(longitude) |
          !is.finite(concentration)  # | concentration <= 0
  bad <- bad0 | !is.finite(rowSums(x))

  logpdf <- rowSums(x) + colatitude + longitude + concentration

  if (any(!bad)) {
    ind4 <- (1:LLL)[!bad]
    xsub <- x[ind4, 1:2, drop = FALSE]
    logpdf[!bad] <- log(concentration[!bad]) -
      log(4 * pi) - log(sinh(concentration[!bad])) +
      (concentration[!bad]) * (
       sin(xsub[, 1]) * sin(colatitude[!bad]) *
       cos(xsub[, 2] - longitude[!bad]) +
       cos(xsub[, 1]) * cos(colatitude[!bad]))
  }
  
  logpdf[!bad0 & is.infinite(rowSums(x))] <- log(0)
  logpdf[ bad0] <- NaN

  if (log.arg) logpdf else exp(logpdf)
}  # dvMF3()









 vMF3 <-
  function(lcolati = "extlogitlink(min = -pi, max = pi)",  #"identitylink",
           llongit = "extlogitlink(min = -pi, max = pi)",  #"identitylink",
           lconcen = "loglink",  # "logitlink",
           icolati = NULL, ilongit = NULL, iconcen = NULL,
           gcolati = exp(2*ppoints(5) - 1),
           glongit = exp(2*ppoints(5) - 1),
           gconcen = exp(2*ppoints(5) - 1),
           tol12 = 1.0e-4,
           zero = NULL) {
  if (is.character(lcolati))
    lcolati <- substitute(y9, list(y9 = lcolati))
  lcolati <- as.list(substitute(lcolati))
  ecolati <- link2list(lcolati)
  lcolati <- attr(ecolati, "function.name")

  if (is.character(llongit))
    llongit <- substitute(y9, list(y9 = llongit))
  llongit <- as.list(substitute(llongit))
  elongit <- link2list(llongit)
  llongit <- attr(elongit, "function.name")

  if (is.character(lconcen))
    lconcen <- substitute(y9, list(y9 = lconcen))
  lconcen <- as.list(substitute(lconcen))
  econcen <- link2list(lconcen)
  lconcen <- attr(econcen, "function.name")




  new("vglmff",
  blurb = c("von Mises-Fisher distribution on the sphere\n\n",
            "Links:    ",
            namesof("colati", lcolati, ecolati, tag = FALSE), ", ",
            namesof("longit", llongit, elongit, tag = FALSE), ", ",
            namesof("concen", lconcen, econcen, tag = FALSE), "\n",
            "Mean:     zz longit * beta(1 + 1 / colati, longit)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero ,
                     M = M, M1 = 3,
                     predictors.names = predictors.names)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 2,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("colati", "longit", "concen"),
         lcolati = .lcolati ,
         llongit = .llongit ,
         lconcen = .lconcen ,
         zero = .zero )
  },
  list( .zero = zero, .lcolati = lcolati,
        .llongit = llongit , .lconcen = lconcen ))),

  initialize = eval(substitute(expression({
    Q1 <- 2  # Bivariate response
    checklist <- w.y.check(w = w, y = y,  # Is.positive.y = TRUE,
                           ncol.w.max = Inf, ncol.y.max = Inf,
                           out.wy = TRUE, colsyperw = Q1,
                           maximize = TRUE)
    w <- checklist$w
    y <- checklist$y  # Now 'w' and 'y' have the correct dimension.

    extra$ncoly <- ncoly <- ncol(y)
    extra$M1 <- M1 <- 3
    M <- M1 * ncoly / Q1
    NOS <- M / M1
    mynames1 <- param.names("colati", ncoly, skip1 = TRUE)
    mynames2 <- param.names("longit", ncoly, skip1 = TRUE)
    mynames3 <- param.names("concen", ncoly, skip1 = TRUE)
    predictors.names <-
      c(namesof(mynames1, .lcolati , earg = .ecolati , tag = FALSE),
        namesof(mynames2, .llongit , earg = .elongit , tag = FALSE),
        namesof(mynames3, .lconcen , earg = .econcen , tag = FALSE))[
        interleave.VGAM(M, M1 = M1)]

    if (!length(etastart)) {
      colati.init <-
      longit.init <-
      concen.init <- matrix(NA_real_, n, NOS)

      vMF3.Loglikfun <- function(colati, y, x, w, extraargs) {
        mediany <- colSums(y * w) / colSums(w)
        longit <- log(0.5) / log1p(-(mediany^colati))
        concen <- log(0.5) / log1p(-(mediany^colati))
        sum(c(w) * dvMF3(y, colati = colati, longit = longit,
                         concen = concen, log = TRUE))
      }

      
      for (spp. in 1:NOS) {  # For each response 'y_spp.'... do:
        yvec <- y[, Q1 * spp. - (1:0)]  # A 2-coln matrix, actually
        wvec <- w[, spp.]
        gcolati <- ( .gcolati )
        glongit <- ( .glongit )
        gconcen <- ( .gconcen )


          if (length( .icolati )) gcolati <- rep_len( .icolati , NOS)
          if (length( .ilongit )) glongit <- rep_len( .ilongit , NOS)
          if (length( .iconcen )) gconcen <- rep_len( .iconcen , NOS)


          ll.vMF3 <- function(concenval, colati, longit,
                              x = x, y = y, w = w, extraargs) {
            ans <- sum(c(w) * dvMF3(x = y,
                                    concen = concenval,
                                    colati = colati,
                                    longit = longit, log = TRUE))
            ans
          }
        try.this <-
          grid.search3(gconcen, gcolati, glongit,
                       objfun = ll.vMF3,
                       y = yvec, w = wvec,
                       ret.objfun = TRUE)  # Last value is the loglik

          concen.init[, spp.] <- try.this["Value1" ]
          colati.init[, spp.] <- try.this["Value2" ]
          longit.init[, spp.] <- try.this["Value3" ]
      }  # End of for (spp. ...)




      etastart <-
        cbind(theta2eta(colati.init, .lcolati , earg = .ecolati ),
              theta2eta(longit.init, .llongit , earg = .elongit ),
              theta2eta(concen.init, .lconcen , earg = .econcen ))[,
                  interleave.VGAM(M, M1 = M1)]
    }
  }), list(
      .lcolati = lcolati, .llongit = llongit, .lconcen = lconcen,
      .icolati = icolati, .ilongit = ilongit, .iconcen = iconcen,
      .ecolati = ecolati, .elongit = elongit, .econcen = econcen,
      .gcolati = gcolati, .glongit = glongit, .gconcen = gconcen
      ))),
  linkinv =  eval(substitute(function(eta, extra = NULL) {
  colati=eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lcolati , .ecolati )
  longit=eta2theta(eta[, c(FALSE, TRUE, FALSE)], .llongit , .elongit )
  concen=eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lconcen , .econcen )
    longit
  }, list(
     .lcolati = lcolati, .llongit = llongit, .lconcen = lconcen,
     .ecolati = ecolati, .elongit = elongit, .econcen = econcen ))),
  last = eval(substitute(expression({
    misc$link <- c(rep_len( .lcolati , M / M1),
                   rep_len( .llongit , M / M1),
                   rep_len( .lconcen , M / M1))[
                   interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2, mynames3)[
                  interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:(M / M1)) {
      misc$earg[[M1*ii-2]] <- ( .ecolati )
      misc$earg[[M1*ii-1]] <- ( .elongit )
      misc$earg[[M1*ii  ]] <- ( .econcen )
    }
  }),
  list(
    .lcolati = lcolati, .llongit = llongit, .lconcen = lconcen,
    .ecolati = ecolati, .elongit = elongit, .econcen = econcen
  ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
  colati=eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lcolati , .ecolati )
  longit=eta2theta(eta[, c(FALSE, TRUE, FALSE)], .llongit , .elongit )
  concen=eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lconcen , .econcen )
    if (residuals) {
      stop("loglikelihood residuals not implemented")
    } else {
      ll.elts <- c(w) * dvMF3(y, colati, longitude = longit,
                              concentration = concen, log = TRUE)
      if (summation) sum(ll.elts) else ll.elts
    }
    },
  list( .lcolati = lcolati, .llongit = llongit,
        .lconcen = lconcen,
        .ecolati = ecolati, .elongit = elongit,
        .econcen = econcen ))),
  vfamily = c("vMF3"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    colati=eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lcolati , .ecolati )
    longit=eta2theta(eta[, c(FALSE, TRUE, FALSE)], .llongit , .elongit )
    concen=eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lconcen , .econcen )
    okay1 <-
      all(is.finite(colati)) &&  # all(0 < colati) &&
      all(is.finite(longit)) &&  # all(0 < longit) &&
      all(is.finite(concen))  # && all(0 < concen)
    if (!okay1)
      cat("not okay\n")
    okay1
  },
  list( .lcolati = lcolati, .llongit = llongit,
        .lconcen = lconcen,
        .ecolati = ecolati, .elongit = elongit,
        .econcen = econcen ))),
  simslot = eval(substitute(
  function(object, nsim) {
    eta <- predict(object)
    colati=eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lcolati , .ecolati )
    longit=eta2theta(eta[, c(FALSE, TRUE, FALSE)], .llongit , .elongit )
    concen=eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lconcen , .econcen )
    rvMF3(nsim * length(colati), colati = colati,
          longit = longit, concen = concen)
  },
  list( .lcolati = lcolati, .llongit = llongit,
       .lconcen = lconcen,
       .ecolati = ecolati, .elongit = elongit,
       .econcen = econcen ))),
  deriv = eval(substitute(expression({
    colati = eta2theta(eta[, c(TRUE, FALSE, FALSE)],
                       .lcolati , .ecolati )
    longit = eta2theta(eta[, c(FALSE, TRUE, FALSE)],
                     .llongit , .elongit )
    concen = eta2theta(eta[, c(FALSE, FALSE, TRUE)],
                       .lconcen , .econcen )
    dcolati.deta <- dtheta.deta(colati, .lcolati , .ecolati )
    dlongit.deta <- dtheta.deta(longit, .llongit , .elongit )
    dconcen.deta <- dtheta.deta(concen, .lconcen , .econcen )

    coth <- function(z) cosh(z) / sinh(z)
    Akappa <- coth(concen) - 1 / concen
    R.rem <- sin(y[, c(TRUE, FALSE)]) *
             sin(colati) *
             cos(y[, c(FALSE, TRUE)] - longit) +
             cos(y[, c(TRUE, FALSE)]) * cos(colati)


    dl.dcolati <- concen * (sin(y[, c(TRUE, FALSE)]) *
                            cos(colati) *
                            cos(y[, c(FALSE, TRUE)] - longit) -
                            cos(y[, c(TRUE, FALSE)]) * sin(colati))
    dl.dlongit <- concen * sin(y[, c(TRUE, FALSE)]) *
                  sin(colati) *
                  sin(y[, c(FALSE, TRUE)] - longit)
    dl.dconcen <- (-Akappa) + R.rem  # 1 / concen - coth(concen) + R.rem
    dl.deta <- c(w) * cbind(dl.dcolati * dcolati.deta,
                            dl.dlongit * dlongit.deta,
                            dl.dconcen * dconcen.deta)
    dl.deta[, interleave.VGAM(M, M1 = M1)]
  }),
  list( .lcolati = lcolati, .llongit = llongit,
       .lconcen = lconcen,
       .ecolati = ecolati, .elongit = elongit,
       .econcen = econcen ))),
  weight = eval(substitute(expression({
    ned2l.dcolati2 <- (-concen) * Akappa  # 1 - coth(concen) / concen
    ned2l.dlongit2 <- concen * Akappa * (sin(colati))^2
    ned2l.dconcen2 <- (1 / sinh(concen))^2 - (1 / concen)^2

    wz <-
      array(c(c(w) * ned2l.dcolati2 * dcolati.deta^2,
              c(w) * ned2l.dlongit2 * dlongit.deta^2,
              c(w) * ned2l.dconcen2 * dconcen.deta^2,
              numeric(n),
              numeric(n),
              numeric(n)),
              dim = c(n, M / M1, 6))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }),
  list( .lcolati = lcolati, .llongit = llongit,
       .lconcen = lconcen,
       .ecolati = ecolati, .elongit = elongit,
       .econcen = econcen, .tol12 = tol12 ))))
}  # vMF3







dcard <- function(x, mu, rho, log = FALSE) {
  if (!isFALSE(log.arg <- log) && !isTRUE(log))
    stop("bad input for argument 'log'")
  rm(log)

  L <- max(length(x), length(mu), length(rho))
  if (length(x)   < L) x   <- rep_len(x,   L)
  if (length(mu)  < L) mu  <- rep_len(mu,  L)
  if (length(rho) < L) rho <- rep_len(rho, L)

  ox <- !is.finite(x)
  zero <- ox | x < 0 | x >= 2 * pi
  ans <- rep_len(log(0), L)
  if (any(!zero))
    ans[!zero] <- -log(2 * pi) + log1p(2 *
         rho[!zero] * cos(x[!zero] - mu[!zero]))
  if (any(ox)) ans[ox] <- log(0)
  ans[mu  <     0] <- NaN
  ans[mu  >= 2*pi] <- NaN
  ans[abs(rho) > 0.5] <- NaN
  if (log.arg) ans else exp(ans)
}  # dcard



pcard <-
    function(q, mu, rho, lower.tail = TRUE,
             log.p = FALSE) {
  if (!isFALSE(lower.tail) && !isTRUE(lower.tail))
    stop("bad input for argument 'lower.tail'")
  if (!isFALSE(log.p) && !isTRUE(log.p))
    stop("bad input for argument 'log.p'")

  ans <- (q + 2 * rho * (sin(q - mu) + sin(mu))) / (2*pi)
  ans[q <= 0] <- 0
  ans[q >= (2*pi)] <- 1
        
  ans[mu < 0 | mu > 2*pi] <- NaN 
  ans[abs(rho) > 0.5] <- NaN

  if (lower.tail) {
    if (log.p) log(ans) else ans
  } else { 
    if (log.p) log1p(-ans) else 1 - ans
  }
}  # pcard





qcard <-
    function(p, mu, rho,
             lower.tail = TRUE, log.p = FALSE,
             tol0 = 1e-6) {

  if (!isFALSE(log.p) && !isTRUE(log.p))
    stop("bad input for argument 'log.p'")
  if (lower.tail) {
    if (log.p) p <- exp(p)
  } else {  # Unperfect
    p <- if (log.p) -expm1(p) else 1 - p
  }
  if (!is.Numeric(tol0, length.arg = 1,
                  positive = TRUE) | tol0 > 0.001)
    stop("bad argument 'tol0'... too large?")

  L <- max(length(p), length(mu), length(rho))
  if (length(p)   < L) p   <- rep_len(p,   L)
  if (length(mu)  < L) mu  <- rep_len(mu,  L)
  if (length(rho) < L) rho <- rep_len(rho, L)

  ans <- p + mu + rho

  bad0 <- !is.finite(mu) | !is.finite(rho) | 
          mu < 0 | mu >= 2 * pi | abs(rho) > 0.5
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- numeric(L)
  hi <- rep(2 * pi, L)
  dont.iterate <- bad  # | is.easy

  foo <- function(q, mu, rho, p)
    pcard(q, mu, rho) - p

  lhs <- dont.iterate  # | (p <= dfoldnorm(0, mean))

  if (any(!lhs)) {
    ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs],
                      mu = mu[!lhs], rho = rho[!lhs],
                      p = p[!lhs], tol = tol0)
  }


  ans[!bad0 & !is.na(p) & p == 0] <- 0
  ans[!bad0 & !is.na(p) & p == 1] <- 2 * pi
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN

  ans
}  # qcard



rcard <-
    function(n, mu, rho) {
  use.n <- if ((lenn <- length(n)) > 1) lenn else
           if (!is.Numeric(n, integer.valued = TRUE,
                length.arg = 1, positive = TRUE))
              stop("bad input for 'n'") else n


  qcard(runif(use.n), mu, rho)
}  # rcard








 cardioid <- function(
     lmu  = "extlogitlink(min = 0, max = 2*pi)",
     lrho = "extlogitlink(min = -0.5, max = 0.5)",
     imu = NULL, irho = NULL,
     gmu  = ppoints(16) * 2 * pi,
     grho = ppoints(16) - 0.5,
     zero = NULL) {

  if (is.character(lmu))
    lmu <- substitute(y9, list(y9 = lmu))
  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  if (is.character(lrho))
    lrho <- substitute(y9, list(y9 = lrho))
  lrho <- as.list(substitute(lrho))
  erho <- link2list(lrho)
  lrho <- attr(erho, "function.name")

 if (length(imu) &&
    (!is.Numeric(imu, positive = TRUE) ||
      any(imu > 2*pi)))
    stop("bad input for argument 'imu'")
  if (length(irho))
    if (!is.Numeric(irho) || max(abs(irho)) > 0.5)
      stop("bad input for argument 'irho'")

  F <- FALSE
  new("vglmff",
  blurb = c("Cardioid distribution\n\n",
            "Links:    ",
            namesof("mu",  lmu,  emu,  tag = F), ", ",
            namesof("rho", lrho, erho, tag = F), "\n",
            "Mean:     ",
            "pi + (rho/pi)*((2*pi-mu)*sin(2*pi-mu)+",
            "cos(2*pi-mu)-mu*sin(mu)-cos(mu))"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x,
             .zero , M = M, M1 = 2,
             predictors.names = predictors.names)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2, M = 2,
         Q1 = 1,
         dpqrfun = "card",
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("mu", "rho"),
         lmu  = .lmu  ,
         lrho = .lrho ,
         zero = .zero )
  },
  list( .zero = zero, .lmu = lmu, .lrho = lrho))),

  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    if (n <= 10)
      warning("n may be too small for Fisher scoring")
    if (any(y < 0 | y >= 2 * pi))
      stop("the response must be in [0, 2*pi)")

    gmuu <- ( .gmu  )
    imuu <- ( .imu  )
    grho <- ( .grho )
    irho <- ( .irho )
    if (length(imuu)) gmuu <- imuu
    if (length(irho)) grho <- irho
      
    predictors.names <- c(
      namesof("mu",  .lmu  , .emu  , tag = FALSE),
      namesof("rho", .lrho , .erho , tag = FALSE))

    if (!length(etastart)) {
      c.ell <- function(mu, rho, y, x, w, extraargs)
        sum(c(w) * dcard(y, mu, rho, log = TRUE))
      try.this <-
        grid.search2(gmuu, grho, objfun = c.ell,
                     y = y, w = w,
                     ret.objfun = TRUE)
      muu.init <-  try.this["Value1"]
      rho.init <-  try.this["Value2"]
      etastart <- matrix(byrow = TRUE, nr = n, nc = 2,
        c(theta2eta(muu.init, .lmu  , .emu  ),
          theta2eta(rho.init, .lrho , .erho )))
    }
  }), list( .lmu = lmu, .lrho = lrho,
            .imu = imu, .irho = irho,
            .gmu = gmu, .grho = grho,
            .emu = emu, .erho = erho ))),
  linkinv = eval(substitute(function(eta, extra = NULL){
    mu  <- eta2theta(eta[, 1], .lmu  , .emu  )
    rho <- eta2theta(eta[, 2], .lrho , .erho )
    pi + (rho / pi) *
    ((2 * pi - mu) * sin(2 * pi - mu) +
    cos(2 * pi - mu) - mu * sin(mu) - cos(mu))
  }, list( .lmu = lmu, .lrho = lrho,
           .emu = emu, .erho = erho ))),
  last = eval(substitute(expression({
    misc$link <-    c("mu" = .lmu , "rho" = .lrho )
    misc$earg <- list("mu" = .emu , "rho" = .erho )
  }),
  list( .lmu = lmu, .lrho = lrho,
        .emu = emu, .erho = erho ))),
  loglikelihood = eval(substitute(
  function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    mu  <- eta2theta(eta[, 1], .lmu  , .emu  )
    rho <- eta2theta(eta[, 2], .lrho , .erho )
    if (residuals) {
      stop("loglikelihood residuals not implemented")
    } else {
      ll.elts <- c(w) * dcard(y, mu, rho, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmu = lmu, .lrho = lrho,
           .emu = emu, .erho = erho ))),
  vfamily = c("cardioid"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    mu  <- eta2theta(eta[, 1], .lmu  , .emu  )
    rho <- eta2theta(eta[, 2], .lrho , .erho )
    okay1 <-
      all(is.finite(mu )) && all(0 <= mu & mu < 2*pi) &&
      all(is.finite(rho)) && all(abs(rho) < 0.5)
    okay1
  }, list( .lmu = lmu, .lrho = lrho,
           .emu = emu, .erho = erho ))),


  simslot = eval(substitute(function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
       pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object); T <- TRUE; F <- FALSE
    mu  <- eta2theta(eta[, c(T, F)], .lmu  , .emu  )
    rho <- eta2theta(eta[, c(F, T)], .lrho , .erho )
    rcard(nsim * length(mu), mu, rho)
  }, list( .lmu = lmu, .lrho = lrho,
           .emu = emu, .erho = erho ))),

  deriv = eval(substitute(expression({
    mu  <- eta2theta(eta[, 1], .lmu  , .emu  )
    rho <- eta2theta(eta[, 2], .lrho , .erho )
    dmu.deta  <- dtheta.deta(mu,  .lmu  , .emu  )
    drho.deta <- dtheta.deta(rho, .lrho , .erho )
    Dnm <- 1 + 2 * rho * cos(y - mu)
    dl.dmu  <- 2 * rho * sin(y - mu) / Dnm
    dl.drho <- 2 *       cos(y - mu) / Dnm
    c(w) * cbind(dl.dmu  *  dmu.deta,
                 dl.drho * drho.deta)
  }),
  list( .lmu = lmu, .lrho = lrho,
        .emu = emu, .erho = erho ))),
  weight = eval(substitute(expression({
    tmp5 <- sqrt(abs(1 - 4 * rho^2))
    ned2l.dmu2  <-  1 - tmp5
    ned2l.drho2 <- (1 - tmp5) / (tmp5 * rho^2)
    wz <- cbind(ned2l.dmu2  * dmu.deta^2,
                ned2l.drho2 * drho.deta^2)
    c(w) * wz
  }),
  list( .lmu = lmu, .lrho = lrho))))
}  # cardioid





 vonmises <-
  function(llocation = "extlogitlink(min = 0, max = 2*pi)",
           lscale  = "loglink",
           ilocation = NULL, iscale  = NULL,
           imethod = 1, zero = NULL) {
  if (is.character(llocation))
    llocation <- substitute(y9, list(y9 = llocation))
  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  if (is.character(lscale))
    lscale <- substitute(y9, list(y9 = lscale))
  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  ilocat <- ilocation


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")




  new("vglmff",
  blurb = c("Von Mises distribution\n\n",
            "Links:    ",
            namesof("location", llocat, elocat), ", ",
            namesof("scale",    lscale, escale),
            "\n", "\n",
            "Mean:     location"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x,
                     .zero , M = M, M1 = 2,
                     predictors.names = predictors.names)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("location", "scale"),
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y)


      predictors.names <-
        c(namesof("location", .llocat , .elocat , tag = FALSE),
          namesof("scale",    .lscale , .escale , tag = FALSE))

      if (!length(etastart)) {
        if ( .imethod == 1) {
          locat.init <- mean(y)
          rat10 <- sqrt((sum(w*cos(y )))^2 + sum(w*sin(y))^2) / sum(w)
          scale.init <- sqrt(1 - rat10)
        } else {
          locat.init <- median(y)
          scale.init <- sqrt(sum(w*abs(y - locat.init)) / sum(w))
        }

        locat.init <- rep_len(if (length( .ilocat )) .ilocat else
                              locat.init,n)
        scale.init <- rep_len(if (length( .iscale )) .iscale else 1, n)
        etastart <- cbind(
            theta2eta(locat.init, .llocat , .elocat ),
            theta2eta(scale.init, .lscale , .escale ))
      }
      y <- y %% (2*pi)  # Coerce after initial values have been computed
  }),
  list( .imethod = imethod, .ilocat = ilocat,
        .escale = escale, .elocat = elocat,
        .lscale = lscale, .llocat = llocat,
        .iscale = iscale ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat , earg = .elocat ) %% (2*pi)
  }, list( .escale = escale, .lscale = lscale,
           .llocat = llocat, .elocat = elocat ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat , scale = .lscale )
    misc$earg <- list(location = .elocat , scale = .escale )



  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    locat <- eta2theta(eta[, 1], .llocat , .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , .escale )

    if (residuals) {
      stop("loglikelihood residuals not implemented")
    } else {
      ll.elts <- c(w) * (Scale * cos(y - locat) -
                         log(mbesselI0(x = Scale)))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
    },
    list( .escale = escale, .lscale = lscale,
          .llocat = llocat, .elocat = elocat ))),
  vfamily = c("vonmises"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    locat <- eta2theta(eta[, 1], .llocat , .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , .escale )
    okay1 <- all(is.finite(locat)) && all(0 < locat & locat < 2*pi) &&
             all(is.finite(Scale)) && all(0 < Scale)
    okay1
  },
  list( .escale = escale, .lscale = lscale,
        .llocat = llocat, .elocat = elocat ))),
  deriv = eval(substitute(expression({
    locat <- eta2theta(eta[, 1], .llocat , .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , .escale )

    tmp6 <- mbesselI0(x = Scale, deriv = 2)
    dl.dlocat <- Scale * sin(y - locat)
    dl.dscale <- cos(y-locat) - tmp6[, 2] / tmp6[, 1]

    dlocat.deta <- dtheta.deta(locat, .llocat ,
                                 earg = .elocat )
    dscale.deta <- dtheta.deta(Scale, .lscale , .escale )

    c(w) * cbind(dl.dlocat * dlocat.deta,
                 dl.dscale * dscale.deta)
  }),
  list( .escale = escale, .lscale = lscale,
        .llocat = llocat, .elocat = elocat ))),
  weight = eval(substitute(expression({
    ned2l.dlocat2 <- Scale * tmp6[, 2] / tmp6[, 1]
    ned2l.dscale2 <- tmp6[, 3] / tmp6[, 1] -
                    (tmp6[, 2] / tmp6[, 1])^2

    wz <- matrix(0, nrow = n, ncol = 2)  # diagonal
    wz[, iam(1, 1, M)] <- ned2l.dlocat2 * dlocat.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dscale2 * dscale.deta^2
    c(w) * wz
  }),
  list( .escale = escale, .elocat = elocat,
        .lscale = lscale, .llocat = llocat ))))
}  # vonmises
















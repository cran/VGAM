# These functions are
# Copyright (C) 1998-2024 T.W. Yee, University of Auckland.
# All rights reserved.

















 dgensh <-
  function(x, shape,
           location = 0, scale = 1,
           tol0 = 1e-4,
           log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  L <- max(length(x), length(location),
           length(scale), length(shape))
  if (length(x)        != L) x        <- rep_len(x,        L)
  if (length(location) != L) location <- rep_len(location, L)
  if (length(scale)    != L) scale    <- rep_len(scale,    L)
  if (length(shape)    != L) shape    <- rep_len(shape,    L)
  z <- (x - location) / scale

  bad0 <- !is.finite(location) |
          !is.finite(scale) |
          !is.finite(shape) |
          scale <= 0 | shape <= -pi
  bad <- bad0 | !is.finite(x)
  c2 <- c1 <- aa <-
  logpdf <- x + shape + scale + location
 
  if (any(!bad)) {
    ind1 <- !bad & shape <= 0  # ,,,,,,,,,,,
    sintt <- ifelse(abs(shape[ind1]) < tol0,
                    rep(1, sum(ind1)),  # limit
                    sin(shape[ind1]) / shape[ind1])
    aa[ind1] <- cos(shape[ind1])
    c2[ind1] <- sqrt((pi^2 - shape[ind1]^2) / 3)
    c1[ind1] <- c2[ind1] * sintt
    ind2 <- !bad & shape > 0  # ,,,,,,,,,,,
    sinhtt <- ifelse(abs(shape[ind2]) < tol0,
                 rep(1, sum(ind2)),  # limit
                 sinh(shape[ind2]) / shape[ind2])
    aa[ind2] <- cosh(shape[ind2])
    c2[ind2] <- sqrt((pi^2 + shape[ind2]^2) / 3)
    c1[ind2] <- c2[ind2] * sinhtt
    logpdf[!bad] <-
       -log(exp( c2[!bad] * z[!bad]) +
            2 * aa[!bad] +
            exp(-c2[!bad] * z[!bad])) +
        log(c1[!bad]) - log(scale[!bad])
  }  # any(!bad)
  
  logpdf[!bad0 & is.infinite(x)] <- log(0)
  logpdf[ bad0] <- NaN
  if (log.arg) logpdf else exp(logpdf)
}  # dgensh




 pgensh <-
  function(q, shape,
           location = 0, scale = 1,
           tol0 = 1e-4,
           lower.tail = TRUE) {
  if (!is.logical(lower.tail) || length(lower.tail) != 1)
    stop("bad input for argument 'lower.tail'")

  L <- max(length(q), length(location),
           length(scale), length(shape))
  if (length(q)        != L) q        <- rep_len(q,        L)
  if (length(location) != L) location <- rep_len(location, L)
  if (length(scale)    != L) scale    <- rep_len(scale,    L)
  if (length(shape)    != L) shape    <- rep_len(shape,    L)
  z <- (q - location) / scale

  bad0 <- !is.finite(location) |
          !is.finite(scale) |
          !is.finite(shape) |
          scale <= 0 | shape <= -pi
  bad <- bad0 | !is.finite(q)
  c2 <-   # c1 <- aa <-
  ans <- z + shape + scale + location

  if (any(!bad)) {
    ind1 <- !bad & shape <= -tol0 &
                   -pi < shape  # ,,,,,,,,,,,
    if (any(ind1)) {
      c2[ind1] <- sqrt((pi^2 - shape[ind1]^2) / 3)
      tmp3a <- -(exp(c2[ind1] * z[ind1]) +
                 cos(shape[ind1])) / (
                 sin(shape[ind1]))
      tmp3a <- atan(1 / tmp3a)
      ans[ind1] <- 1 + tmp3a / shape[ind1]
      if (any(shape[ind1] < -pi / 2)) {
        tmp3b <- -(exp(c2[ind1] * abs(z[ind1])) +
                   cos(shape[ind1])) / (
                   sin(shape[ind1]))
        tmp3b <- atan(1 / tmp3b)
        ans3b <- 1 + tmp3b / shape[ind1]
      }  # any(shape[ind1] < -pi / 2)
      ans[ind1] <-
        ifelse(shape[ind1] < -pi / 2,
               ifelse(z[ind1] < 0, 
                      1 - ans3b, # Fix
                      ans[ind1]),
               ans[ind1])
    }  # any(ind1)
    ind2 <- !bad & tol0 <= shape  # ,,,,,,,,,,,
    acoth <- function(x)
        0.5 * log((-1 - x) / (1 - x))
    c2[ind2] <- sqrt((pi^2 + shape[ind2]^2) / 3)
    ans[ind2] <- 1 - acoth( (exp(c2[ind2] *
                                  z[ind2]) +
          cosh(shape[ind2])) / sinh(shape[ind2])
          ) / shape[ind2]
    ind3 <- !bad & abs(shape) <= tol0  # ,,,,,,,
    ans[ind3] <- logitlink(pi * z[ind3] / sqrt(3),
                           inverse = TRUE)
  }  # any(!bad)
  ans[!bad0 & is.infinite(q) & q > 0] <- 1
  ans[!bad0 & is.infinite(q) & q < 0] <- 0
  ans[ bad0] <- NaN
  if (lower.tail) ans else 1 - ans
}  # pgensh
  
  



 qgensh <-
  function(p, shape,
           location = 0, scale = 1,
           tol0 = 1e-4) {
  L <- max(length(p), length(location),
           length(scale), length(shape))
  if (length(p)        != L) p        <- rep_len(p,        L)
  if (length(location) != L) location <- rep_len(location, L)
  if (length(scale)    != L) scale    <- rep_len(scale,    L)
  if (length(shape)    != L) shape    <- rep_len(shape,    L)

  bad0 <- !is.finite(location) |
          !is.finite(scale) |
          !is.finite(shape) |
          scale <= 0 | shape <= -pi
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p
  c2 <-    # c1 <- aa <-
  ans <- p + shape + scale + location


  if (any(!bad)) {
    ind1 <- !bad & shape <= -tol0 &
                   -pi < shape  # ,,,,,,,,,,,
    c2[ind1] <- sqrt((pi^2 - shape[ind1]^2) / 3)
    ans[ind1] <- log(sin(shape[ind1] * p[ind1])
    / sin(shape[ind1] * (1 - p[ind1]))) / c2[ind1]
    ind2 <- !bad & tol0 <= shape  # ,,,,,,,,,,,
    c2[ind2] <- sqrt((pi^2 + shape[ind2]^2) / 3)
    ans[ind2] <- log(sinh(shape[ind2] * p[ind2])
    / sinh(shape[ind2] * (1 - p[ind2]))) / c2[ind2]
    ind3 <- !bad & abs(shape) <= tol0  # ,,,,,,,
    ans[ind3] <- logitlink(p[ind3]) * sqrt(3) / pi

    ans[!bad] <- location[!bad] +
                 scale[!bad] * ans[!bad]
  }  # any(!bad)

  ans[!bad0 & !is.na(p) & p == 0] <- -Inf
  ans[!bad0 & !is.na(p) & p == 1] <- Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgensh
  
  

rgensh <-
  function(n, shape,
           location = 0, scale = 1,
           tol0 = 1e-04) {

  use.n <- if ((length.n <- length(n)) > 1)
             length.n else n
  qgensh(runif(use.n), shape = shape,
         location = location, scale = scale,
         tol0 = tol0)         
}  # rgensh











 gensh <-
  function(shape,
           llocation = "identitylink",
           lscale = "loglink",
           zero = "scale",
           ilocation = NULL, iscale = NULL,
           imethod = 1,
           glocation.mux = exp((-4:4)/2),
           gscale.mux = exp((-4:4)/2),
           probs.y = 0.3,
           tol0 = 1e-4) {

  if (!is.Numeric(shape, length.arg = 1) ||
      shape < -pi)
    stop("bad input for argument 'shape'")
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



  new("vglmff",
  blurb = c("Generalized secant hyperbolic ",
        "distribution (Vaughan, 2002)\n",
        "f(y; aa, b, shape) = (c1/b) ",
        "/ (exp(c2*z) + 2*aa + exp(-c2*",
        "z)),", "\n",
        "z = (y - a) / b,\n",
   "aa = cos(shape) or cosh(shape),\n",
   "c1 = c2 * sin(shape) / shape or ",
        "c2 * sinh(shape) / shape,\n",
   "c2 = sqrt((pi^2 - shape^2) / 3) or ",
        "sqrt((pi^2 + shape^2) / 3),\n",
   "location = a, scale = b > 0, shape > -pi\n\n",
            "Links:    ",
       namesof("location", llocat, elocat), ", ",
       namesof("scale",    lscale, escale),
       "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints,
             x = x, .zero , M = M, M1 = 2,
             predictors.names = predictors.names)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         dpqrfun = "gensh",
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("location", "scale"),
         imethod = .imethod ,
         llocation  = .llocat ,
         lscale     = .lscale ,
         shape      = .shape ,
         zero = .zero )
  }, list( .zero = zero,
           .imethod = imethod ,
           .llocat  = llocat ,
           .lscale  = lscale ,
           .shape   = shape ))),

  initialize = eval(substitute(expression({
    M1 <- 2; T <- TRUE
    Q1 <- 1

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = FALSE,
              Is.integer.y = FALSE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    NOS <- ncoly <- ncol(y)  # Number of species
    M <- M1 * ncoly


    tm1 <- param.names("location", NOS, skip1 = T)
    tm2 <- param.names("scale",    NOS, skip1 = T)
    predictors.names <-
      c(namesof(tm1, .llocat , .elocat ,
                tag = FALSE),
        namesof(tm2, .lscale , .escale ,
                tag = FALSE))
    predictors.names <-
    predictors.names[interleave.VGAM(M, M1 = M1)]



    if (!length(etastart)) {
      lo.init <-
      sc.init <- matrix(NA_real_, n, NOS)
      if (length( .ilocat ))
      lo.init <- matrix( .ilocat , n, NOS, by = T)
      if (length( .iscale ))
      sc.init <- matrix( .iscale , n, NOS, by = T)

      for (spp. in 1:NOS) {  # For each 'y_spp.':
        yvec <- y[, spp.]
        wvec <- w[, spp.]
        mu.init <- switch( .imethod ,
           weighted.mean(yvec, w = wvec),
           median(yvec),  # More reliable?
           quantile(yvec, prob = .probs.y ),
           stop("argument 'imethod' unmatched"))

          glocat  <- .glocat.mux * mu.init
          gscale  <- .gscale.mux * abs(mu.init)
   if (length( .ilocat ))
     glocat  <- rep_len( .ilocat , NOS)
   if (length( .iscale ))
     gscale  <- rep_len( .iscale , NOS)


    ll.gensh <- function(scaleval, locn,
      x = x, y = y, w = w, extraargs) {
      ans <- sum(c(w) * dgensh(y,
                        scale = scaleval,
                        locat = locn,
                        shape = extraargs$shape,
                        log = TRUE))
      ans
    }
    try.this <-
      grid.search2(gscale, glocat,
          objfun = ll.gensh,
          y = yvec, w = wvec,
          extraargs = list(shape = .shape ),
          ret.objfun = TRUE)  # Last val is \ell

      sc.init[, spp.] <- try.this["Value1" ]
      lo.init[, spp.] <- try.this["Value2" ]


      }  # End of for (spp. ...)


  etastart <-
    cbind(theta2eta(lo.init, .llocat , .elocat ),
          theta2eta(sc.init, .lscale , .escale ))
  etastart <- etastart[, interleave.VGAM(M,
                                         M1 = M1)]
    }
  }),
  list( .llocat = llocat, .lscale = lscale,
        .elocat = elocat, .escale = escale,
        .ilocat = ilocat, .iscale = iscale,
        .imethod = imethod ,
        .glocat.mux = glocation.mux,
        .gscale.mux = gscale.mux,
        .shape      = shape,
        .probs.y = probs.y
        ))),
  linkinv = eval(substitute(
    function(eta, extra = NULL) {
      eta2theta(eta[, c(TRUE, FALSE)],
                .llocat , .elocat )
    },
  list( .llocat = llocat, .lscale = lscale,
        .elocat = elocat, .escale = escale ))),
  last = eval(substitute(expression({
    tmp34 <- c(rep_len( .llocat , NOS),
               rep_len( .lscale , NOS))
    names(tmp34) <- c(tm1, tm2)
    tmp34 <- tmp34[interleave.VGAM(M, M1 = M1)]
    misc$link <- tmp34  # Already named

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- ( .elocat )
      misc$earg[[M1*ii  ]] <- ( .escale )
    }
  }),
  list( .llocat = llocat, .lscale = lscale,
        .elocat = elocat, .escale = escale,
        .shape  = shape))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    TF1 <- c(TRUE, FALSE)
    TF2 <- c(FALSE, TRUE)
    a <- eta2theta(eta[, TF1], .llocat , .elocat )
    b <- eta2theta(eta[, TF2], .lscale , .escale )
    if (residuals) {
      stop("loglikelihood resids not implemented")
    } else {
      ll.elts <-
        c(w) * dgensh(y, loc = a, scale = b,
                      shape = .shape , log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  list( .llocat = llocat, .lscale = lscale,
        .shape = shape,
        .elocat = elocat, .escale = escale))),
  vfamily = c("gensh"),
  validparams = eval(substitute(
      function(eta, y, extra = NULL) {
    TF1 <- c(TRUE, FALSE)
    TF2 <- c(FALSE, TRUE)
    aa <- eta2theta(eta[, TF1], .llocat , .elocat )
    bb <- eta2theta(eta[, TF2], .lscale , .escale )
    kk <- c( .shape )
    okay1 <- 
        all(is.finite(bb)) && all(0 < bb) &&
        all(is.finite(aa))  # &&
    okay1
  },
  list( .llocat = llocat, .lscale = lscale,
        .shape  = shape,
        .elocat = elocat, .escale = escale))),




  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- M / M1
    TF1 <- c(TRUE, FALSE)
    TF2 <- c(FALSE, TRUE)
    locat <- eta2theta(eta[, TF1, drop = FALSE],
                       .llocat , .elocat )
    Scale <- eta2theta(eta[, TF2, drop = FALSE],
                       .lscale , .escale )
    shape <- c( .shape )
    zedd <- as.matrix((y - locat) / Scale)

    ind1 <- shape <= 0  # ,,,,,,,,,,,
    if (ind1) {
      sintt <- ifelse(abs(shape) < .tol0 ,
                      1,  # limit
                      sin(shape) / shape)
      aa <- cos(shape)
      c2 <- sqrt((pi^2 - shape^2) / 3)
      c1 <- c2 * sintt
    } else {
      sinhtt <- ifelse(abs(shape) < .tol0 ,
                       1,  # limit
                       sinh(shape) / shape)
      aa <- cosh(shape)
      c2 <- sqrt((pi^2 + shape^2) / 3)
      c1 <- c2 * sinhtt
    }


    dl.dzedd <-
        -c2 * (exp( c2 * zedd)  -
               exp(-c2 * zedd)) / (
       exp( c2 * zedd) + 2 * aa +
       exp(-c2 * zedd))
    dz.dlocat <-    -1 / Scale
    dz.dscale <- -zedd / Scale
    dl.dlocat <- dl.dzedd * dz.dlocat
    dl.dscale <- -1 / Scale + dl.dzedd * dz.dscale

    dlocat.deta <- dtheta.deta(locat, .llocat ,
                               .elocat )
    dscale.deta <- dtheta.deta(Scale, .lscale ,
                               .escale )

    myderiv <-
      c(w) * cbind(dl.dlocat * dlocat.deta,
                   dl.dscale * dscale.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }),
  list( .llocat = llocat, .lscale = lscale,
        .elocat = elocat, .escale = escale,
        .shape  = shape,  .tol0   = tol0))),
  weight = eval(substitute(expression({

    if (ind1) {
     ned2l.dlocat2 <-  c2^2 * (
        shape - sin(shape) * cos(shape)) / (
        2 * Scale^2 * shape * (sin(shape))^2)
      ned2l.dscale2 <-
        ((pi^2 - shape^2) / (sin(shape))^2 -
        ((pi^2 - 3 * shape^2) * cos(shape)) / (
         shape * sin(shape))) / ( 6 * Scale^2)
    } else {
      ned2l.dlocat2 <- -c2^2 * (
        shape - sinh(shape) * cosh(shape)) / (
        2 * Scale^2 * shape * (sinh(shape))^2)
      ned2l.dscale2 <-
        ((pi^2 + shape^2) / (sinh(shape))^2 -
        ((pi^2 + 3 * shape^2) * cosh(shape)) / (
         shape * sinh(shape))) / (-6 * Scale^2)
    }

    if (abs(shape) < .tol0 ) {
      ned2l.dlocat2 <- (pi^2 / 3) / (3 * Scale^2)
      ned2l.dscale2 <- (3 + pi^2) / (9 * Scale^2)
    }
    

    wz <- 
    array(c(c(w) * ned2l.dlocat2 * dlocat.deta^2,
            c(w) * ned2l.dscale2 * dscale.deta^2),
          dim = c(n, M / M1, 2))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }),
  list( .tol0 = tol0))))
}  # gensh
































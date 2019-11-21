# These functions are
# Copyright (C) 1998-2019 T.W. Yee, University of Auckland.
# All rights reserved.













moments.nbin.gait <-
  function(size.p, munb.p,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           byrow.ai = FALSE,  # For pobs.a and pstr.i
           size.a = size.p, size.i = size.p,
           munb.a = munb.p, munb.i = munb.p,
           mlm = TRUE,  # mlm = FALSE iff mixture = TRUE
           type.fitted = "All") {  # or "mean"

  rmlife <- if (is.finite(max.support)) NA else 0
  NOS <- 1
  nnn <- length(size.p)
  lhs.prob <- pnbinom(max.support, size.p, mu = munb.p)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)

  sumt.p <- matrix(0, nnn, NOS)
  SumT.p <- matrix(rmlife, nnn, NOS)
  if (ltrunc)
    for (tval in truncate) {
      pmf.p <- dnbinom(tval, size.p, mu = munb.p)
      sumt.p <- sumt.p + pmf.p  # Need tval<=max.support
      SumT.p <- SumT.p + pmf.p * tval
    }

  use.pobs.a <- use.pstr.i <- matrix(0, nnn, 1)  # So rowSums() works.
  aprd <- 0  # aprd is an innerprod
  SumA.x <-  # For innerprod (aprd) only
  SumA.a <- suma.a <-
  SumA.p <- suma.p <- matrix(0, nnn, NOS)
  if (mlm) {
    if (lalter)
      use.pobs.a <- matrix(pobs.a, nnn, lalter, byrow = byrow.ai)
    if (linfla)
      use.pstr.i <- matrix(pstr.i, nnn, linfla, byrow = byrow.ai)
  } else {
    if (lalter)
      use.pobs.a <- matrix(pobs.a, nnn, 1)
    if (linfla)
      use.pstr.i <- matrix(pstr.i, nnn, 1)
  }
  if (lalter) {
    for (jay in seq_len(lalter)) {
      aval <- alter[jay]
      pmf.x <- if (mlm) use.pobs.a[, jay] else rep(0, nnn)
      pmf.p <- dnbinom(aval, size.p, mu = munb.p)
      pmf.a <- dnbinom(aval, size.a, mu = munb.a)
      suma.p <- suma.p + pmf.p
      SumA.p <- SumA.p + pmf.p * aval
      suma.a <- suma.a + pmf.a
      SumA.a <- SumA.a + pmf.a * aval
      SumA.x <- SumA.x + pmf.x * aval
    }  # for jay
    aprd <- if (mlm) SumA.x else use.pobs.a * SumA.a / suma.a
  }  # lalter

  iprd <- 0  # iprd is an innerprod
  sumi.i <- SumI.i <-
  sumi.p <- SumI.p <- matrix(0, nnn, NOS)
  if (linfla) {
    for (jay in seq_len(linfla)) {
      ival <- inflate[jay]
      pmf.p <- if (mlm) use.pstr.i[, jay] else
               dnbinom(ival, size.p, mu = munb.p)
      pmf.i <- if (mlm) use.pstr.i[, jay] else
               dnbinom(ival, size.i, mu = munb.i)
      sumi.p <- sumi.p + pmf.p
      SumI.p <- SumI.p + pmf.p * ival
      sumi.i <- sumi.i + pmf.i
      SumI.i <- SumI.i + pmf.i * ival
    }  # for jay
    iprd <- if (mlm) SumI.p else use.pstr.i * SumI.i / sumi.i
  }  # linfla

  use.this <- if (mlm)
    (1 - rowSums(use.pobs.a) - rowSums(use.pstr.i)) else
    (1 - use.pobs.a - use.pstr.i)

  themean <- aprd + iprd + use.this *
    (munb.p - SumA.p - SumT.p) / (lhs.prob - suma.p - sumt.p)
  if (type.fitted == "mean") {
    return(themean)
  }
      
  ans <- list('lhs.prob' = lhs.prob,
              'rmlife'   = rmlife,
              'sumt.p'   = sumt.p,
              'SumT.p'   = SumT.p,
              'suma.p'   = suma.p,
              'SumA.p'   = SumA.p,
              'sumi.p'   = sumi.p,
              'SumI.p'   = SumI.p,
              'suma.a'   = suma.a,
              'SumA.a'   = SumA.a,
              'sumi.i'   = sumi.i,
              'SumI.i'   = SumI.i,
              'aprd'     = aprd,
              'iprd'     = iprd,
              'mean'     = themean)
  ans
}  # moments.nbin.gait







 gatnbinomial.mix <-
  function(alter = NULL,  # May need a (>=3)zz-vector to run
           truncate = NULL,
           zero = c("pobs.a", "size"),
           parallel = FALSE,  # TRUE applies to the intercept
           lpobs.a = "logitlink",
           lmunb.p = "loglink",
           lsize.p = "loglink",
           lmunb.a = "loglink",
           lsize.a = "loglink",
           type.fitted = c("mean", "pobs.a",
                           "Pobs.a", "prob.a", "prob.t"),
           imethod = 1,
           imunb.p = NULL,    isize.p = NULL,  # exp(1),
           imunb.a = imunb.p, isize.a = isize.p,
           ishrinkage = 0.95,
           probs.y = 0.35,

           cutoff.prob = 0.999,  # higher is better for large 'size'
           eps.trig = 1e-7,
           max.chunk.MB = 30  # max.memory = Inf is allowed
           ) {
  max.support <- Inf  # Fixed for now, as its intractable
  lpobs.a <- as.list(substitute(lpobs.a))
  epobs.a <- link2list(lpobs.a)
  lpobs.a <- attr(epobs.a, "function.name")

  lmunb.p <- as.list(substitute(lmunb.p))
  emunb.p <- link2list(lmunb.p)
  lmunb.p <- attr(emunb.p, "function.name")
  lsize.p <- as.list(substitute(lsize.p))
  esize.p <- link2list(lsize.p)
  lsize.p <- attr(esize.p, "function.name")

  lmunb.a <- as.list(substitute(lmunb.a))
  emunb.a <- link2list(lmunb.a)
  lmunb.a <- attr(emunb.a, "function.name")
  lsize.a <- as.list(substitute(lsize.a))
  esize.a <- link2list(lsize.a)
  lsize.a <- attr(esize.a, "function.name")

  gait.errorcheck(alter, inflate = NULL, truncate, max.support)

  lalter <- length(alter)
  ltrunc <- length(truncate)
  type.fitted <- match.arg(type.fitted, c("mean", "pobs.a",
                   "Pobs.a", "prob.a", "prob.t"))[1]
  tmp3 <- c(pobs.a = lpobs.a,
            munb.p = lmunb.p,
            size.p = lsize.p,
            munb.a = lmunb.a,
            size.a = lsize.a)
  blurb1 <- "N"
  if (lalter) blurb1 <- "Generally-altered n"
  if (ltrunc) blurb1 <- "Generally-truncated n"
  if (lalter && ltrunc) blurb1 <- "Generally-altered and -truncated n"
                
  new("vglmff",
  blurb = c(blurb1, "egative binomial distribution ",
            "(GAT-NB-NB mixture)\n\n",
            "Links:   ",
            namesof("pobs.a", lpobs.a, epobs.a, tag = FALSE), ", ",
            namesof("munb.p", lmunb.p, emunb.p, tag = FALSE), ", ",
            namesof("size.p", lsize.p, esize.p, tag = FALSE), ", ",
            "\n         ",
            namesof("munb.a", lmunb.a, emunb.a, tag = FALSE), ", ",
            namesof("size.a", lsize.a, esize.a, tag = FALSE), "\n"),

  constraints = eval(substitute(expression({
    tmp <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1), 5, 3)
    constraints <- cm.VGAM(tmp, x = x,
                           bool = .parallel ,
                           constraints = constraints,
                           apply.int = TRUE)  # FALSE
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 5)
  }), list( .zero = zero,
            .parallel = parallel ))),

  infos = eval(substitute(function(...) {
    list(M1 = 5,
         Q1 = 1,
         link = .tmp3 ,
         link1parameter = TRUE,
         mixture.links = FALSE,
         alter = .alter ,
         truncate = .truncate ,
         max.support = .max.support , 
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE
         parameters.names = c("pobs.a", "munb.p", "size.p",
                                        "munb.a", "size.a"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted,
           .alter = alter,
           .truncate = truncate,
           .max.support = max.support, 
           .tmp3 = tmp3 ))),

  initialize = eval(substitute(expression({
    truncate <- as.vector( .truncate )
    alter <- as.vector( .alter )
    ltrunc <- length(truncate)
    lalter <- length(alter)
    if (lalter < 3)
      stop("argument 'alter' must contain at least 3 values")
    M1 <- 5
    NOS <- NCOL(y)
    M <- NOS * M1

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (ltrunc && any(y %in% truncate))
      stop("some response values == values in argument 'truncate'")
    if ( .max.support < max(y))
      stop("some response values > than argument 'max.support'")

    y0 <- matrix(0, n, lalter)
    for (jay in seq(lalter))
      y0[, jay] <- as.numeric(y == alter[jay])
    extra$skip.these <- matrix(as.logical(y0), n, lalter)  # dim lost
      if (any((css <- colSums(extra$skip.these)) == 0))
        stop("some 'alter' argument values have no response values: ",
             paste(alter[css == 0], collapse = ", "))          

    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)  # May be NULL
    extra$M1 <- M1

    mynames1 <- param.names("pobs.a", ncoly, skip1 = TRUE)
    mynames2 <- param.names("munb.p", ncoly, skip1 = TRUE)
    mynames3 <- param.names("size.p", ncoly, skip1 = TRUE)
    mynames4 <- param.names("munb.a", ncoly, skip1 = TRUE)
    mynames5 <- param.names("size.a", ncoly, skip1 = TRUE)
    predictors.names <-
      c(namesof(mynames1, .lpobs.a , earg = .epobs.a , tag = FALSE),
        namesof(mynames2, .lmunb.p , earg = .emunb.p , tag = FALSE),
        namesof(mynames3, .lsize.p , earg = .esize.p , tag = FALSE),
        namesof(mynames4, .lmunb.a , earg = .emunb.a , tag = FALSE),
        namesof(mynames5, .lsize.a , earg = .esize.a , tag = FALSE))

    if (!length(etastart)) {
      munb.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                           imu = .imunb.p ,  # x = x,
                           ishrinkage = .ishrinkage ,
                           probs.y = .probs.y )
      try.this <- if (length( .isize.p )) .isize.p else
        0.25 * (0.01 + weighted.mean(y, w)^2) / (0.01 + var(y))
      size.init <- matrix(try.this, n, NCOL(munb.init))
      po.a.init <- rep(sum(css) / n, n)  # MLE for pobs.a (unwted)

      etastart <-
        cbind(theta2eta(po.a.init, .lpobs.a , earg = .epobs.a ),
              theta2eta(munb.init, .lmunb.p , earg = .emunb.p ),
              theta2eta(size.init, .lsize.p , earg = .esize.p ),
              theta2eta(munb.init, .lmunb.a , earg = .emunb.a ),
              theta2eta(size.init, .lsize.a , earg = .esize.a ))
    }
  }), list( .lmunb.p = lmunb.p, .lsize.p = lsize.p, .lpobs.a = lpobs.a,
            .lmunb.a = lmunb.a, .lsize.a = lsize.a,
            .emunb.p = emunb.p, .esize.p = esize.p, .epobs.a = epobs.a,
            .emunb.a = emunb.a, .esize.a = esize.a,
            .imunb.p = imunb.p, .isize.p = isize.p,
            .imunb.a = imunb.a, .isize.a = isize.a,
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .imethod = imethod,
            .alter = alter,
            .truncate = truncate, .max.support = max.support, 
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <-
      if (length(extra$type.fitted)) extra$type.fitted else {
        warning("cannot find 'type.fitted'. Returning the 'mean'.")
        "mean"  # Unconditional mean
      }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "pobs.a",
                       "Pobs.a", "prob.a", "prob.t"))[1]
    alter <- as.vector( .alter )
    lalter <- length(alter)
    M1 <- 5
    NOS <- NCOL(eta) / M1
    if ((NOS <- NCOL(eta) / M1) != 1)
      stop("Currently NOS must be 1")
    pobs.a <- eta2theta(eta[, 1], .lpobs.a , earg = .epobs.a )
    munb.p <- eta2theta(eta[, 2], .lmunb.p , earg = .emunb.p )
    size.p <- eta2theta(eta[, 3], .lsize.p , earg = .esize.p )
    munb.a <- eta2theta(eta[, 4], .lmunb.a , earg = .emunb.a )
    size.a <- eta2theta(eta[, 5], .lsize.a , earg = .esize.a )
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    max.support <- as.vector( .max.support )

    bits <- moments.nbin.gait(size.p = size.p, munb.p = munb.p,
                              alter = alter,
                              truncate = truncate, mlm = FALSE,
                              max.support = max.support,
                              pobs.a = pobs.a, munb.a = munb.a,
                              size.a = size.a)

    if (type.fitted == "Pobs.a") {
      mymat <-
        dnbinom(matrix(alter, NROW(eta), lalter, byrow = TRUE),
                size = matrix(size.a, NROW(eta), lalter),
                mu   = matrix(munb.a, NROW(eta), lalter)) / (
        c(bits[["suma.a"]]))
    }
  
    ans <- switch(type.fitted,
                  "mean"   = bits[["mean"]],
                  "pobs.a" = pobs.a,
                  "Pobs.a" = pobs.a * mymat,  # matrix
                  "prob.a" = bits[["suma.p"]],
                  "prob.t" = bits[["sumt.p"]])
   ans <-
   label.cols.y(ans,
                 colnames.y = if (type.fitted  == "Pobs.a")
                              as.character(alter) else extra$colnames.y,
                 NOS = NOS)
   ans
  } , list( .lmunb.p = lmunb.p, .lsize.p = lsize.p, .lpobs.a = lpobs.a,
            .lmunb.a = lmunb.a, .lsize.a = lsize.a,
            .emunb.p = emunb.p, .esize.p = esize.p, .epobs.a = epobs.a,
            .emunb.a = emunb.a, .esize.a = esize.a,
            .alter = alter,
            .truncate = truncate, .max.support = max.support ))),

  last = eval(substitute(expression({
    misc$link  <- c(pobs.a = .lpobs.a ,
                    munb.p = .lmunb.p ,
                    size.p = .lsize.p ,
                    munb.a = .lmunb.a ,
                    size.a = .lsize.a )

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- .epobs.a  #
    misc$earg[[2]] <- .emunb.p  #
    misc$earg[[3]] <- .esize.p  #
    misc$earg[[4]] <- .emunb.a  #
    misc$earg[[5]] <- .esize.a  #
  }), list( .lmunb.p = lmunb.p, .lsize.p = lsize.p, .lpobs.a = lpobs.a,
            .lmunb.a = lmunb.a, .lsize.a = lsize.a,
            .emunb.p = emunb.p, .esize.p = esize.p, .epobs.a = epobs.a,
            .emunb.a = emunb.a, .esize.a = esize.a,
            .alter = alter,
            .truncate = truncate, .max.support = max.support ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    pobs.a <- eta2theta(eta[, 1], .lpobs.a , earg = .epobs.a )
    munb.p <- eta2theta(eta[, 2], .lmunb.p , earg = .emunb.p )
    size.p <- eta2theta(eta[, 3], .lsize.p , earg = .esize.p )
    munb.a <- eta2theta(eta[, 4], .lmunb.a , earg = .emunb.a )
    size.a <- eta2theta(eta[, 5], .lsize.a , earg = .esize.a )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitnbinom.mix(y,
                          size.p = size.p, munb.p = munb.p,
                          size.a = size.a, munb.a = munb.a,
                          pobs.a = pobs.a, truncate = .truncate ,
                          max.support = .max.support ,
                          alter = .alter , log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  } , list( .lmunb.p = lmunb.p, .lsize.p = lsize.p, .lpobs.a = lpobs.a,
            .lmunb.a = lmunb.a, .lsize.a = lsize.a,
            .emunb.p = emunb.p, .esize.p = esize.p, .epobs.a = epobs.a,
            .emunb.a = emunb.a, .esize.a = esize.a,
            .alter = alter,
            .truncate = truncate, .max.support = max.support ))),
  vfamily = c("gatnbinomial.mix"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pobs.a <- eta2theta(eta[, 1], .lpobs.a , earg = .epobs.a )
    munb.p <- eta2theta(eta[, 2], .lmunb.p , earg = .emunb.p )
    size.p <- eta2theta(eta[, 3], .lsize.p , earg = .esize.p )
    munb.a <- eta2theta(eta[, 4], .lmunb.a , earg = .emunb.a )
    size.a <- eta2theta(eta[, 5], .lsize.a , earg = .esize.a )

    okay1 <- all(is.finite(munb.p)) && all(0 <  munb.p) &&
             all(is.finite(munb.a)) && all(0 <  munb.a) &&
             all(is.finite(size.p)) && all(0 <  size.p) &&
             all(is.finite(size.a)) && all(0 <  size.a) &&
             all(is.finite(pobs.a)) && all(0 <= pobs.a & pobs.a < 1)
    okay1
  } , list( .lmunb.p = lmunb.p, .lsize.p = lsize.p, .lpobs.a = lpobs.a,
            .lmunb.a = lmunb.a, .lsize.a = lsize.a,
            .emunb.p = emunb.p, .esize.p = esize.p, .epobs.a = epobs.a,
            .emunb.a = emunb.a, .esize.a = esize.a,
            .alter = alter,
            .truncate = truncate, .max.support = max.support ))),
  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs.a <- eta2theta(eta[, 1], .lpobs.a , earg = .epobs.a )
    munb.p <- eta2theta(eta[, 2], .lmunb.p , earg = .emunb.p )
    size.p <- eta2theta(eta[, 3], .lsize.p , earg = .esize.p )
    munb.a <- eta2theta(eta[, 4], .lmunb.a , earg = .emunb.a )
    size.a <- eta2theta(eta[, 5], .lsize.a , earg = .esize.a )
    rgaitnbinom.mix(nsim * length(size.p),
                      size.p = size.p, munb.p = munb.p,
                      size.a = size.a, munb.a = munb.a,
                      pobs.a = pobs.a, truncate = .truncate ,
                      max.support = .max.support ,
                      alter = .alter )
  } , list( .lmunb.p = lmunb.p, .lsize.p = lsize.p, .lpobs.a = lpobs.a,
            .lmunb.a = lmunb.a, .lsize.a = lsize.a,
            .emunb.p = emunb.p, .esize.p = esize.p, .epobs.a = epobs.a,
            .emunb.a = emunb.a, .esize.a = esize.a,
            .alter = alter,
            .truncate = truncate, .max.support = max.support ))),
  deriv = eval(substitute(expression({
    alter <- as.vector( .alter )
    lalter <- length(alter)
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    M1 <- 5
    NOS <- NCOL(eta) / M1  # extra$NOS
    if (NOS != 1) stop("Multiple responses not handled")
    max.support <- as.vector( .max.support )
    if (lalter) {
      is.altered <- rowSums(extra$skip.these) > 0
    }
    pobs.a <- eta2theta(eta[, 1], .lpobs.a , earg = .epobs.a )
    munb.p <- eta2theta(eta[, 2], .lmunb.p , earg = .emunb.p )
    size.p <- eta2theta(eta[, 3], .lsize.p , earg = .esize.p )
    munb.a <- eta2theta(eta[, 4], .lmunb.a , earg = .emunb.a )
    size.a <- eta2theta(eta[, 5], .lsize.a , earg = .esize.a )

    bits <- moments.nbin.gait(munb.p, size.p = size.p, alter = alter,
                              truncate = truncate, mlm = FALSE,
                              max.support = max.support,
                              pobs.a = pobs.a, munb.a = munb.a,
                              size.a = size.a)

    pmf.deriv1.m <- function(y, munb, size)
      (y / munb - 1) * dnbinom(y, size, mu = munb) / (
      1 + munb / size)
    pmf.deriv2.m <- function(y, munb, size)
      size * (munb * (1 + size) * (munb - 2 * y) +
              y * (y - 1) * size) *  # zz bad if size \approx Inf
      dnbinom(y, size, mu = munb) / (munb * (munb + size))^2
    pmf.deriv1.s <- function(y, munb, size)
      (digamma(y + size) - digamma(size) -
      log1p(munb / size) - (y - munb) / (munb + size)) *
      dnbinom(y, size, mu = munb)
    pmf.deriv2.s <- function(y, munb, size)
      (trigamma(y + size) - trigamma(size) +
      munb / (size * (munb + size)) + (y - munb) / (munb + size)^2) *
      dnbinom(y, size, mu = munb) +
      (digamma(y + size) - digamma(size) -
      log1p(munb / size) - (y - munb) / (munb + size)) *
      pmf.deriv1.s(y, munb, size)
    pmf.deriv.sm <- function(y, munb, size)
      ((y - munb) / (munb + size)^2) * dnbinom(y, size, mu = munb) +
      ((y / munb - 1) / (1 + munb / size)) *
      pmf.deriv1.s(y, munb, size)


    sumderiv1a.a <- sumderiv2a.a <-
    sumderivxa.p <- sumderivxa.a <-
    sumderiv1a.p <- sumderiv2a.p <- matrix(0, NROW(eta), NOS)
    deriv0matrix <-
    deriv1matrix <- deriv2matrix <- matrix(0, NROW(eta), lalter)
    SumDeriv1a.a <- SumDeriv2a.a <-
    SumDeriv1a.p <- SumDeriv2a.p <- matrix(0, NROW(eta), NOS)
    DerivxMatrix <-  # Deriv0Matrix <-
    Deriv1Matrix <- Deriv2Matrix <- matrix(0, NROW(eta), lalter)
    for (jay in seq(lalter)) {
      aval <- alter[jay]
      sumderiv1a.p <- sumderiv1a.p + pmf.deriv1.m(aval, munb.p, size.p)
      sumderiv2a.p <- sumderiv2a.p + pmf.deriv2.m(aval, munb.p, size.p)
      sumderiv1a.a <- sumderiv1a.a + pmf.deriv1.m(aval, munb.a, size.a)
      sumderiv2a.a <- sumderiv2a.a + pmf.deriv2.m(aval, munb.a, size.a)
      sumderivxa.p <- sumderivxa.p + pmf.deriv.sm(aval, munb.p, size.p)
      sumderivxa.a <- sumderivxa.a + pmf.deriv.sm(aval, munb.a, size.a)
      SumDeriv1a.p <- SumDeriv1a.p + pmf.deriv1.s(aval, munb.p, size.p)
      SumDeriv2a.p <- SumDeriv2a.p + pmf.deriv2.s(aval, munb.p, size.p)
      SumDeriv1a.a <- SumDeriv1a.a + pmf.deriv1.s(aval, munb.a, size.a)
      SumDeriv2a.a <- SumDeriv2a.a + pmf.deriv2.s(aval, munb.a, size.a)
      pmf.a <- dnbinom(aval, size.a, mu = munb.a)
      deriv0matrix[, jay] <- pmf.a
      deriv1matrix[, jay] <- pmf.deriv1.m(aval, munb.a, size.a) / pmf.a
      deriv2matrix[, jay] <- pmf.deriv2.m(aval, munb.a, size.a) / pmf.a
      Deriv1Matrix[, jay] <- pmf.deriv1.s(aval, munb.a, size.a) / pmf.a
      Deriv2Matrix[, jay] <- pmf.deriv2.s(aval, munb.a, size.a) / pmf.a
      DerivxMatrix[, jay] <- pmf.deriv.sm(aval, munb.a, size.a) / pmf.a


      
    }
    deriv0matrix <-  deriv0matrix / rowSums(deriv0matrix)  # Normalized
      
    sumderivxt.p <- sumderivxt.a <-
    sumderiv1t.a <- sumderiv2t.a <-
    sumderiv1t.p <- sumderiv2t.p <- matrix(0, NROW(eta), NOS)
    SumDeriv1t.a <- SumDeriv2t.a <-
    SumDeriv1t.p <- SumDeriv2t.p <- matrix(0, NROW(eta), NOS)
    if (ltrunc)
      for (tval in truncate) {
        sumderiv1t.p <- sumderiv1t.p + pmf.deriv1.m(tval, munb.p, size.p)
        sumderiv2t.p <- sumderiv2t.p + pmf.deriv2.m(tval, munb.p, size.p)
        sumderiv1t.a <- sumderiv1t.a + pmf.deriv1.m(tval, munb.a, size.a)
        sumderiv2t.a <- sumderiv2t.a + pmf.deriv2.m(tval, munb.a, size.a)
        sumderivxt.p <- sumderivxt.p + pmf.deriv.sm(tval, munb.p, size.p)
        SumDeriv1t.p <- SumDeriv1t.p + pmf.deriv1.s(tval, munb.p, size.p)
        SumDeriv2t.p <- SumDeriv2t.p + pmf.deriv2.s(tval, munb.p, size.p)
        SumDeriv1t.a <- SumDeriv1t.a + pmf.deriv1.s(tval, munb.a, size.a)
        SumDeriv2t.a <- SumDeriv2t.a + pmf.deriv2.s(tval, munb.a, size.a)
      }


    onempobs.a <- 1 - pobs.a
    dl.dmunb.p <- dl.dsize.p <-
    dl.dmunb.a <- dl.dsize.a <- zero0n <- rep(0, n)
    dl.dpobs.a <- ifelse(is.altered, 1 / pobs.a, -1 / onempobs.a) 
    for (jay in seq(lalter)) {
      aval <- alter[jay]
      dl.dmunb.a <- dl.dmunb.a +
        ifelse(extra$skip[, jay],
               deriv1matrix[, jay] - sumderiv1a.a / bits[["suma.a"]],
               zero0n)
      dl.dsize.a <- dl.dsize.a +
        ifelse(extra$skip[, jay],
               Deriv1Matrix[, jay] - SumDeriv1a.a / bits[["suma.a"]],
               zero0n)
    }  # jay

    Denom <- bits[["lhs.prob"]] - bits[["suma.p"]] - bits[["sumt.p"]]
    dl.dmunb.p <-
      ifelse(is.altered,
             zero0n,
             (y / munb.p - 1) / (1 + munb.p / size.p) +
             (sumderiv1a.p + sumderiv1t.p) / Denom)
    dl.dsize.p <-
      ifelse(is.altered,
             zero0n,
             digamma(y + size.p) - digamma(size.p) -
             log1p(munb.p / size.p) -
             (y - munb.p) / (munb.p + size.p) +
             (SumDeriv1a.p + SumDeriv1t.p) / Denom)

    dpobs.a.deta <- dtheta.deta(pobs.a, .lpobs.a , .epobs.a )
    dmunb.p.deta <- dtheta.deta(munb.p, .lmunb.p , .emunb.p )
    dsize.p.deta <- dtheta.deta(size.p, .lsize.p , .esize.p )
    dmunb.a.deta <- dtheta.deta(munb.a, .lmunb.a , .emunb.a )
    dsize.a.deta <- dtheta.deta(size.a, .lsize.a , .esize.a )
    ans <- cbind(dl.dpobs.a * dpobs.a.deta,
                 dl.dmunb.p * dmunb.p.deta,
                 dl.dsize.p * dsize.p.deta,
                 dl.dmunb.a * dmunb.a.deta,
                 dl.dsize.a * dsize.a.deta)
    c(w) * ans
  }), list( .lmunb.p = lmunb.p, .lsize.p = lsize.p, .lpobs.a = lpobs.a,
            .lmunb.a = lmunb.a, .lsize.a = lsize.a,
            .emunb.p = emunb.p, .esize.p = esize.p, .epobs.a = epobs.a,
            .emunb.a = emunb.a, .esize.a = esize.a,
            .alter = alter,
            .truncate = truncate, .max.support = max.support ))),
  weight = eval(substitute(expression({
    cond.EY.p <-
      (munb.p - bits[["SumA.p"]] - bits[["SumT.p"]]) / Denom


    ned2l.dmunb.p2 <- onempobs.a * (
      cond.EY.p * (1 / munb.p^2 - 1 / (munb.p + size.p)^2) -
      1 / (munb.p / sqrt(size.p) + sqrt(size.p))^2 -
      (sumderiv2a.p + sumderiv2t.p) / Denom -
     ((sumderiv1a.p + sumderiv1t.p) / Denom)^2)


    ned2l.dmunb.p.size.p2 <- onempobs.a * (
     -(cond.EY.p - munb.p) / (munb.p + size.p)^2 -
      (sumderivxa.p + sumderivxt.p) / Denom -
      (sumderiv1a.p + sumderiv1t.p) *
      (SumDeriv1a.p + SumDeriv1t.p) / Denom^2)


    ned2l.dmunb.a2 <- pobs.a *
      rowSums(deriv0matrix * (
             -deriv2matrix +
              deriv1matrix^2 +
             (c(sumderiv2a.a) / c(bits[["suma.a"]])) -
             (c(sumderiv1a.a) / c(bits[["suma.a"]]))^2))


    ned2l.dsize.a2 <- pobs.a *
      rowSums(deriv0matrix * (
             -Deriv2Matrix +
              Deriv1Matrix^2 +
             (c(SumDeriv2a.a) / c(bits[["suma.a"]])) -
             (c(SumDeriv1a.a) / c(bits[["suma.a"]]))^2))


    ned2l.dmunb.a.size.a2 <- pobs.a *
      rowSums(deriv0matrix * (
             -DerivxMatrix +
              deriv1matrix * Deriv1Matrix +
            c(sumderivxa.a) / c(bits[["suma.a"]]) -
            c(sumderiv1a.a) * c(SumDeriv1a.a) / c(bits[["suma.a"]])^2))



    wz11 <- if ( .lpobs.a == "logitlink") {
      pobs.a * (1 - pobs.a)
    } else {
      dpobs.a.deta <- dtheta.deta(pobs.a, .lpobs.a , .epobs.a )
      ned2l.dpobs.a2 <- 1 / (pobs.a * (1 - pobs.a))
      ned2l.dpobs.a2 * dpobs.a.deta^2
    }






    min.support <- 0  # Usual case
    min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else
    min.support
    max.chunk.MB <- ( .max.chunk.MB )
    MM12 <- M1 * (M1 + 1) / 2
    wz <- matrix(0, n, MM12)  # A full matrix

    ned2l.dsize.p2 <- matrix(0, n, NOS)
    
    ind2 <- matrix(FALSE, n, NOS)  # Used for SFS
    for (jay in 1:NOS) {
      eff.p <- sort(c( .cutoff.prob , 1 - .cutoff.prob ))
      Q.mins  <- min.support.use  # 1
      Q.mins2 <- pmax(qgaitnbinom.mix(eff.p[1],
                                    size.p = size.p,
                                    munb.p = munb.p,
                                    alter = alter,
                                    truncate = truncate,
                                    max.support = max.support) - 2,
                      Q.mins)
      Q.maxs <-  qgaitnbinom.mix(p       = eff.p[2] ,
                               size.p      = size.p,
                               munb.p      = munb.p,
                               alter       = .alter ,
                               truncate    = .truncate ,
                               max.support = .max.support ) + 10
      eps.trig <- as.vector( .eps.trig )
      Q.MAXS <- pmax(10, ceiling(1 / sqrt(eps.trig)), na.rm = TRUE)
      Q.maxs <- pmin(Q.maxs, Q.MAXS, na.rm = TRUE)


      ind1 <- if (max.chunk.MB > 0)
                (Q.maxs - Q.mins < max.support) else FALSE
      if ((NN <- sum(ind1)) > 0) {
        Object.Size <- NN * 8 * max(Q.maxs - Q.mins) / (2^20)
        n.chunks <- if (intercept.only) 1 else
                    max(1, ceiling( Object.Size / max.chunk.MB))
        chunk.rows <- ceiling(NN / n.chunks)
        ind2  <- ind1  # Save this
        wind2 <- which(ind1)


        upr.ptr <- 0
        lwr.ptr <- upr.ptr + 1
        while (lwr.ptr <= NN) {
          upr.ptr <- min(upr.ptr + chunk.rows, NN)
          sind2 <- wind2[lwr.ptr:upr.ptr]
          eim.gatnb.p <-
          EIM.GATNB.speciald(munb = munb.p[sind2],
            size        = size.p[sind2],
            munb.a = munb.a[sind2], size.a = size.a[sind2],
            alter = alter, truncate = truncate,
            max.support = max.support,
            pobs.a      = pobs.a[sind2],  # pobs.a, zz
            EY.cond     = cond.EY.p[sind2],
            y.min       = min.support,  # min(Q.mins   ),
            y.max       = max(Q.maxs[sind2]),
            cutoff.prob = .cutoff.prob ,
            intercept.only = intercept.only,
            mlm         = FALSE,
            extra.bit   = TRUE)
          ned2l.dsize.p2[sind2] <- eim.gatnb.p

          lwr.ptr <- upr.ptr + 1
        }  # while (lwr.ptr <= NN)
      }  # if ((NN <- sum(ind1)) > 0)
    }  # end of for (jay in 1:NOS)




    ned2l.dsize.p2 <- onempobs.a * (
      ned2l.dsize.p2 -
      (SumDeriv2a.p + SumDeriv2t.p) / Denom -
     ((SumDeriv1a.p + SumDeriv1t.p) / Denom)^2)






    wz <- wz.merge(wz11,
            cbind(ned2l.dmunb.p2 * dmunb.p.deta^2,
                  ned2l.dsize.p2 * dsize.p.deta^2,
                  ned2l.dmunb.a2 * dmunb.a.deta^2,
                  ned2l.dsize.a2 * dsize.a.deta^2,
                  ned2l.dmunb.p.size.p2 * dmunb.p.deta * dsize.p.deta,
                  zero0n,
                  ned2l.dmunb.a.size.a2 * dmunb.a.deta * dsize.a.deta),
                  M1 = 1, M2 = 4)


    
   


    wz
  }), list( .alter = alter,
            .truncate = truncate, .max.support = max.support,
            .cutoff.prob = cutoff.prob, .eps.trig = eps.trig,
            .max.chunk.MB = max.chunk.MB,
            .lpobs.a = lpobs.a
           ))))
}  # gatnbinomial.mix
















moments.pois.gait <-
  function(lambda.p,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           byrow.ai = FALSE,  # For pobs.a and pstr.i
           lambda.a = lambda.p, lambda.i = lambda.p,
           mlm = TRUE,  # mlm = FALSE iff mixture = TRUE
           type.fitted = "All") {  # or "mean"

  rmlife <- lambda.p * ppois(max.support - 1, lambda.p,
                             lower.tail = FALSE)
  NOS <- 1
  nnn <- length(lambda.p)
  lhs.prob <- ppois(max.support, lambda.p)  # 1 - P(upper tail)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)

  sumt.p <- matrix(0, nnn, NOS)
  SumT.p <- matrix(rmlife, nnn, NOS)
  if (ltrunc)
    for (tval in truncate) {
      pmf.p <- dpois(tval, lambda.p)
      sumt.p <- sumt.p + pmf.p  # Need tval<=max.support
      SumT.p <- SumT.p + pmf.p * tval
    }

  use.pobs.a <- use.pstr.i <- matrix(0, nnn, 1)  # So rowSums() works.
  aprd <- 0  # aprd is an innerprod
  SumA.x <-  # For innerprod (aprd) only
  SumA.a <- suma.a <-
  SumA.p <- suma.p <- matrix(0, nnn, NOS)
  if (mlm) {
    if (lalter)
      use.pobs.a <- matrix(pobs.a, nnn, lalter, byrow = byrow.ai)
    if (linfla)
      use.pstr.i <- matrix(pstr.i, nnn, linfla, byrow = byrow.ai)
  } else {
    if (lalter)
      use.pobs.a <- matrix(pobs.a, nnn, 1)
    if (linfla)
      use.pstr.i <- matrix(pstr.i, nnn, 1)
  }
  if (lalter) {
    for (jay in seq_len(lalter)) {
      aval <- alter[jay]
      pmf.x <- if (mlm) use.pobs.a[, jay] else rep(0, nnn)
      pmf.p <- dpois(aval, lambda.p)
      pmf.a <- dpois(aval, lambda.a)
      suma.p <- suma.p + pmf.p
      SumA.p <- SumA.p + pmf.p * aval
      suma.a <- suma.a + pmf.a
      SumA.a <- SumA.a + pmf.a * aval
      SumA.x <- SumA.x + pmf.x * aval
    }  # for jay
    aprd <- if (mlm) SumA.x else use.pobs.a * SumA.a / suma.a
  }  # lalter

  iprd <- 0  # iprd is an innerprod
  sumi.i <- SumI.i <-
  sumi.p <- SumI.p <- matrix(0, nnn, NOS)
  if (linfla) {
    for (jay in seq_len(linfla)) {
      ival <- inflate[jay]
      pmf.p <- if (mlm) use.pstr.i[, jay] else dpois(ival, lambda.p)
      pmf.i <- if (mlm) use.pstr.i[, jay] else dpois(ival, lambda.i)
      sumi.p <- sumi.p + pmf.p
      SumI.p <- SumI.p + pmf.p * ival
      sumi.i <- sumi.i + pmf.i
      SumI.i <- SumI.i + pmf.i * ival
    }  # for jay
    iprd <- if (mlm) SumI.p else use.pstr.i * SumI.i / sumi.i
  }  # linfla

  use.this <- if (mlm)
    (1 - rowSums(use.pobs.a) - rowSums(use.pstr.i)) else
    (1 - use.pobs.a - use.pstr.i)

  themean <- aprd + iprd + use.this *
    (lambda.p - SumA.p - SumT.p) / (lhs.prob - suma.p - sumt.p)
  if (type.fitted == "mean") {
    return(themean)
  }
      
  ans <- list('lhs.prob' = lhs.prob,
              'rmlife'   = rmlife,
              'sumt.p'   = sumt.p,
              'SumT.p'   = SumT.p,
              'suma.p'   = suma.p,
              'SumA.p'   = SumA.p,
              'sumi.p'   = sumi.p,
              'SumI.p'   = SumI.p,
              'suma.a'   = suma.a,
              'SumA.a'   = SumA.a,
              'sumi.i'   = sumi.i,
              'SumI.i'   = SumI.i,
              'aprd'     = aprd,
              'iprd'     = iprd,
              'mean'     = themean)
  ans
}  # moments.pois.gait






dgaitnbinom.mix <-
  function(x, size.p, prob.p = NULL, munb.p = NULL,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           munb.a = munb.p, munb.i = munb.p,
           log.arg = FALSE) {
  gait.errorcheck(alter, inflate, truncate, max.support)
  if ((is.prob <- as.logical(length(prob.p))) &&
      length(munb.p))
    stop("cannot specify both 'prob.p' and 'munb.p' arguments")
  if (!length(c(alter, inflate, truncate)) &&
      is.infinite(max.support)) {
    return(if (is.prob)
      dnbinom(x, size = size.p, prob = prob.p, log = log.arg) else
      dnbinom(x, size = size.p, mu   = munb.p, log = log.arg))
  }
  if (min(pobs.a, na.rm = TRUE) < 0 || 1 < max(pobs.a, na.rm = TRUE))
    stop("argument 'pobs.a' out of range")
  if (min(pstr.i, na.rm = TRUE) < 0 || 1 < max(pstr.i, na.rm = TRUE))
    stop("argument 'pstr.i' out of range")
  allx.a <- if (length(alter))   0:max(alter  ) else NULL
  allx.i <- if (length(inflate)) 0:max(inflate) else NULL
  use.ms <- if (is.finite(max.support)) {
    whatsleft <- setdiff(0:max.support, alter)
    if (length(whatsleft) > 0)
      max(whatsleft) else max.support
  } else {
    max.support
  }
  pmf1 <- dgaitnbinom.mlm(x,  # Inner distribution
            size = size.p,
            munb = munb.p,
            prob = prob.p,
            truncate = c(alter[alter < use.ms], truncate),  # No ivec
            max.support = use.ms)
  pmf2.a <- pmf2.i <- 0
  if (length(alter))
    pmf2.a <- dgaitnbinom.mlm(x,  # An outer distribution
                size = size.a,
                munb = munb.a,
                prob = prob.a,
                truncate = setdiff(allx.a, alter),
                max.support = max(alter))
  if (length(inflate))
    pmf2.i <- dgaitnbinom.mlm(x,  # An outer distribution
                size = size.i,
                munb = munb.i,
                prob = prob.i,
                truncate = setdiff(allx.i, inflate),
                max.support = max(inflate))
  ans <- pobs.a * pmf2.a + pstr.i * pmf2.i +
         (1 - pobs.a - pstr.i) * pmf1
  if (log.arg) log(ans) else ans
}  # dgaitnbinom.mix





pgaitnbinom.mix <-
  function(q, size.p, prob.p = NULL, munb.p = NULL,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           munb.a = munb.p, munb.i = munb.p) {
  gait.errorcheck(alter, inflate, truncate, max.support)
  if ((is.prob <- as.logical(length(prob.p)) &&
      length(munb.p)))
    stop("cannot specify both 'prob.p' and 'munb.p' arguments")
  if (!length(c(alter, inflate, truncate)) &&
      is.infinite(max.support)) {
    return(if (is.prob)
      pnbinom(q, size = size.p, prob = prob.p) else
      pnbinom(q, size = size.p, mu   = munb.p))
  }
  if (min(pobs.a, na.rm = TRUE) < 0 || 1 < max(pobs.a, na.rm = TRUE))
    stop("argument 'pobs.a' out of range")
  if (min(pstr.i, na.rm = TRUE) < 0 || 1 < max(pstr.i, na.rm = TRUE))
    stop("argument 'pstr.i' out of range")
  allx.a <- if (length(alter))   0:max(alter  ) else NULL
  allx.i <- if (length(inflate)) 0:max(inflate) else NULL
  use.ms <- if (is.finite(max.support)) {
    whatsleft <- setdiff(0:max.support, alter)
    if (length(whatsleft) > 0)
      max(whatsleft) else max.support
  } else {
    max.support
  }
  cdf1 <- pgaitnbinom.mlm(q,  # Inner distribution
            size = size.p,
            munb = munb.p,
            prob = prob.p,
            truncate = c(alter[alter < use.ms], truncate),  # No ivec
            max.support = use.ms)

  cdf2.a <- cdf2.i <- 0
  if (length(alter))
    cdf2.a <- pgaitnbinom.mlm(q,  # An outer distribution
                size = size.a,
                munb = munb.a,
                prob = prob.a,
                truncate = setdiff(allx.a, alter),
                max.support = max(alter))
  if (length(inflate))
    cdf2.i <- pgaitnbinom.mlm(q,  # An outer distribution
                size = size.i,
                munb = munb.i,
                prob = prob.i,
                truncate = setdiff(allx.i, inflate),
                max.support = max(inflate))
  ans <- pobs.a * cdf2.a + pstr.i * cdf2.i +
         (1 - pobs.a - pstr.i) * cdf1
  ans
}  # pgaitnbinom.mix






qgaitnbinom.mix <-
  function(p, size.p, prob.p = NULL, munb.p = NULL,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           munb.a = munb.p, munb.i = munb.p) {
  gait.errorcheck(alter, inflate, truncate, max.support)
  if ((is.prob <- as.logical(length(prob.p)) &&
      length(munb.p)))
    stop("cannot specify both 'prob.p' and 'munb.p' arguments")
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (!length(c(alter, inflate, truncate)) &&
      is.infinite(max.support)) {
    return(if (is.prob)
      qnbinom(p, size = size.p, prob = prob.p) else
      qnbinom(p, size = size.p, mu   = munb.p))
  }

  if (min(pobs.a, na.rm = TRUE) < 0 || 1 < max(pobs.a, na.rm = TRUE))
    stop("argument 'pobs.a' out of range")
  if (min(pstr.i, na.rm = TRUE) < 0 || 1 < max(pstr.i, na.rm = TRUE))
    stop("argument 'pstr.i' out of range")

  LLL <- max(length(p),
             length(size.p), length(prob.p), length(munb.p),
             length(size.a), length(prob.a), length(munb.a),
             length(size.i), length(prob.i), length(munb.i),
             length(pobs.a), length(pstr.i))
  if (length(p)        != LLL) p        <- rep_len(p,        LLL)
  if (length(size.p)   != LLL) size.p   <- rep_len(size.p,   LLL)
  if (length(size.a)   != LLL) size.a   <- rep_len(size.a,   LLL)
  if (length(size.i)   != LLL) size.i   <- rep_len(size.i,   LLL)
  if (is.prob) {
  if (length(prob.p)   != LLL) prob.p   <- rep_len(prob.p,   LLL)
  if (length(prob.a)   != LLL) prob.a   <- rep_len(prob.a,   LLL)
  if (length(prob.i)   != LLL) prob.i   <- rep_len(prob.i,   LLL)
  } else {
  if (length(munb.p)   != LLL) munb.p   <- rep_len(munb.p,   LLL)
  if (length(munb.a)   != LLL) munb.a   <- rep_len(munb.a,   LLL)
  if (length(munb.i)   != LLL) munb.i   <- rep_len(munb.i,   LLL)
  }
  if (length(pobs.a)   != LLL) pobs.a   <- rep_len(pobs.a,   LLL)
  if (length(pstr.i)   != LLL) pstr.i   <- rep_len(pstr.i,   LLL)

  min.support <- 0  # Usual case
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else
    min.support
  ans <- p + size.p + size.a + size.i + pobs.a + pstr.i
  ans <- ans + (if (is.prob) prob.p + prob.a + prob.i else
         munb.p + munb.a + munb.i)

  bad0 <- !is.finite(size.p) | size.p <= 0
  if ( is.prob)
    bad0 <- bad0 |
         (!is.finite(prob.p) | prob.p <= 0 | 1 <= prob.p)
  if (!is.prob)
    bad0 <- bad0 |
         (!is.finite(munb.p) | munb.p <= 0)
  if ( is.prob)
    bad0 <- bad0 | (lalter &
         (!is.finite(prob.a) | prob.a <= 0 | 1 <= prob.a))
  if (!is.prob)
    bad0 <- bad0 | (lalter &
         (!is.finite(munb.a) | munb.a <= 0))
  if ( is.prob)
    bad0 <- bad0 | (linfla &
         (!is.finite(prob.i) | prob.i <= 0 | 1 <= prob.i))
  if (!is.prob)
    bad0 <- bad0 | (linfla &
         (!is.finite(munb.i) | munb.i <= 0))
  bad0 <- bad0 | (lalter &
         (!is.finite(size.a) | size.a <= 0))
  bad0 <- bad0 | (linfla &
         (!is.finite(size.i) | size.i <= 0))
  bad0 <- bad0 | (lalter &
         (!is.finite(pobs.a) | pobs.a <= 0 | 1 <= pobs.a))
  bad0 <- bad0 | (linfla &
         (!is.finite(pstr.i) | pstr.i <= 0 | 1 <= pstr.i))
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <- if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else 2 * lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
          p <= pgaitnbinom.mix(hi,
                 size.p = size.p, prob.p = prob.p, munb.p = munb.p,
                 size.a = size.a, prob.a = prob.a, munb.a = munb.a,
                 size.i = size.i, prob.i = prob.i, munb.i = munb.i,
                               alter = alter, inflate = inflate,
                               truncate = truncate,
                               max.support = max.support)
  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {



    hi[!done] <- 2 * hi[!done] + 10.5  # Bug fixed
    hi <- pmin(max.support + 0.5, hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitnbinom.mix(hi[!done],
size.p = size.p[!done], prob.p = prob.p[!done], munb.p = munb.p[!done],
size.a = size.a[!done], prob.a = prob.a[!done], munb.a = munb.a[!done],
size.i = size.i[!done], prob.i = prob.i[!done], munb.i = munb.i[!done],
                                   pobs.a  = pobs.a[!done],
                                   pstr.i  = pstr.i[!done],
                                   alter   = alter,
                                   inflate = inflate,
                                   truncate = truncate,
                                   max.support = max.support))
    iter <- iter + 1
  }  # while

  hi <- pmin(max.support + 0.5, 2 * hi)  # 20191108



  
      foo <- function(q, size.p, prob.p = NULL, munb.p = NULL,
                      pobs.a = 0, pstr.i = 0,
                      size.a = size.p,
                      prob.a = prob.p,
                      munb.a = munb.p,
                      size.i = size.p,
                      prob.i = prob.p,
                      munb.i = munb.p,
                      alter = NULL,
                      inflate = NULL, truncate = NULL,
                      max.support = Inf, p)
    pgaitnbinom.mix(q,
                    size.p = size.p, prob.p = prob.p, munb.p = munb.p,
                    size.a = size.a, prob.a = prob.a, munb.a = munb.a,
                    size.i = size.i, prob.i = prob.i, munb.i = munb.i,
                    pobs.a = pobs.a, pstr.i = pstr.i,
                    alter = alter, inflate = inflate,
                    truncate = truncate,
                    max.support = max.support) - p

  lhs <- dont.iterate |
         p <= dgaitnbinom.mix(min.support.use,
                    size.p = size.p, prob.p = prob.p, munb.p = munb.p,
                    size.a = size.a, prob.a = prob.a, munb.a = munb.a,
                    size.i = size.i, prob.i = prob.i, munb.i = munb.i,
                    pobs.a = pobs.a, pstr.i = pstr.i,
                    alter  = alter, inflate  = inflate,
                    truncate = truncate,
                    max.support = max.support)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
size.p = size.p[!lhs], prob.p = prob.p[!lhs], munb.p = munb.p[!lhs],
size.a = size.a[!lhs], prob.a = prob.a[!lhs], munb.a = munb.a[!lhs],
size.i = size.i[!lhs], prob.i = prob.i[!lhs], munb.i = munb.i[!lhs],
                      pobs.a   = pobs.a[!lhs],
                      pstr.i   = pstr.i[!lhs],
                      alter    = alter,
                      inflate  = inflate, truncate = truncate,
                      max.support = max.support,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])


    tmp <-
      ifelse(pgaitnbinom.mix(faa,
size.p = size.p[!lhs], prob.p = prob.p[!lhs], munb.p = munb.p[!lhs],
size.a = size.a[!lhs], prob.a = prob.a[!lhs], munb.a = munb.a[!lhs],
size.i = size.i[!lhs], prob.i = prob.i[!lhs], munb.i = munb.i[!lhs],
                             pobs.a = pobs.a[!lhs],
                             pstr.i = pstr.i[!lhs],
                             alter = alter,
                             inflate = inflate,
                             truncate = truncate,
                             max.support = max.support) < p[!lhs] &
             p[!lhs] <= pgaitnbinom.mix(faa+1,
size.p = size.p[!lhs], prob.p = prob.p[!lhs], munb.p = munb.p[!lhs],
size.a = size.a[!lhs], prob.a = prob.a[!lhs], munb.a = munb.a[!lhs],
size.i = size.i[!lhs], prob.i = prob.i[!lhs], munb.i = munb.i[!lhs],
                                        pobs.a = pobs.a[!lhs],
                                        pstr.i = pstr.i[!lhs],
                                        alter = alter,
                                        inflate = inflate,
                                        truncate = truncate,
                                        max.support = max.support),
             faa+1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitnbinom.mix(min.support.use, size.p = size.p,
                    prob.p = prob.p, munb.p = munb.p,
                    size.a = size.a,
                    prob.a = prob.a,
                    munb.a = munb.a,
                    size.i = size.i,
                    prob.i = prob.i,
                    munb.i = munb.i,
                                pobs.a = pobs.a,
                                pstr.i = pstr.i,
                                alter = alter,
                                inflate = inflate,
                                truncate = truncate,
                                max.support = max.support)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitnbinom.mix





rgaitnbinom.mix <-
  function(n, size.p, prob.p = NULL, munb.p = NULL,
           alter = NULL, inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           munb.a = munb.p, munb.i = munb.p) {
    qgaitnbinom.mix(runif(n), size.p = size.p,
                    prob.p = prob.p, munb.p = munb.p,
                    alter = alter, inflate = inflate,
                    truncate = truncate,
                    max.support = max.support,
                    pobs.a = pobs.a, pstr.i = pstr.i,
                    size.a = size.a,
                    prob.a = prob.a,
                    munb.a = munb.a,
                    size.i = size.i,
                    prob.i = prob.i,
                    munb.i = munb.i)
}  # rgaitnbinom.mix








dgaitpois.mix <-
  function(x, lambda.p,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           lambda.a = lambda.p, lambda.i = lambda.p,
           log.arg = FALSE) {
  gait.errorcheck(alter, inflate, truncate, max.support)
  if (!length(c(alter, inflate, truncate)) &&
      is.infinite(max.support)) {
    return(dpois(x, lambda.p, log = log.arg))
  }
  if (min(pobs.a, na.rm = TRUE) < 0 || 1 < max(pobs.a, na.rm = TRUE))
    stop("argument 'pobs.a' out of range")
  if (min(pstr.i, na.rm = TRUE) < 0 || 1 < max(pstr.i, na.rm = TRUE))
    stop("argument 'pstr.i' out of range")
  allx.a <- if (length(alter))   0:max(alter  ) else NULL
  allx.i <- if (length(inflate)) 0:max(inflate) else NULL
  use.ms <- if (is.finite(max.support)) {
    whatsleft <- setdiff(0:max.support, alter)
    if (length(whatsleft) > 0)
      max(whatsleft) else max.support
  } else {
    max.support
  }
  pmf1 <- dgaitpois.mlm(x, lambda.p,  # Inner distribution
            truncate = c(alter[alter < use.ms], truncate),  # No ivec
            max.support = use.ms)

  pmf2.a <- pmf2.i <- 0
  if (length(alter))
    pmf2.a <- dgaitpois.mlm(x, lambda.a,  # Outer distribution
                          truncate = setdiff(allx.a, alter),
                          max.support = max(alter))
  if (length(inflate))
    pmf2.i <- dgaitpois.mlm(x, lambda.i,  # Outer distribution
                          truncate = setdiff(allx.i, inflate),
                          max.support = max(inflate))
  ans <- pobs.a * pmf2.a + pstr.i * pmf2.i +  # mixprob * pmf2ai +
         (1 - pobs.a - pstr.i) * pmf1
  if (log.arg) log(ans) else ans
}  # dgaitpois.mix








pgaitpois.mix <-
  function(q, lambda.p,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           lambda.a = lambda.p, lambda.i = lambda.p) {
  gait.errorcheck(alter, inflate, truncate, max.support)
  if (!length(c(alter, inflate, truncate)) &&
      is.infinite(max.support)) {
    return(ppois(q, lambda.p))
  }
  if (min(pobs.a, na.rm = TRUE) < 0 || 1 < max(pobs.a, na.rm = TRUE))
    stop("argument 'pobs.a' out of range")
  if (min(pstr.i, na.rm = TRUE) < 0 || 1 < max(pstr.i, na.rm = TRUE))
    stop("argument 'pstr.i' out of range")
  allx.a <- if (length(alter))   0:max(alter  ) else NULL
  allx.i <- if (length(inflate)) 0:max(inflate) else NULL
  use.ms <- if (is.finite(max.support)) {
    whatsleft <- setdiff(0:max.support, alter)
    if (length(whatsleft) > 0)
      max(whatsleft) else max.support
  } else {
    max.support
  }
  cdf1 <- pgaitpois.mlm(q, lambda.p,  # Inner distribution
            truncate = c(alter[alter < use.ms], truncate),  # No ivec
            max.support = use.ms)

  cdf2.a <- cdf2.i <- 0
  if (length(alter))
    cdf2.a <- pgaitpois.mlm(q, lambda.a,  # Outer distribution
                          truncate = setdiff(allx.a, alter),
                          max.support = max(alter))
  if (length(inflate))
    cdf2.i <- pgaitpois.mlm(q, lambda.i,  # Outer distribution
                          truncate = setdiff(allx.i, inflate),
                          max.support = max(inflate))
  ans <- pobs.a * cdf2.a + pstr.i * cdf2.i +  # mixprob * pmf2ai +
         (1 - pobs.a - pstr.i) * cdf1
  ans
}  # pgaitpois.mix









qgaitpois.mix <-
  function(p, lambda.p,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           lambda.a = lambda.p, lambda.i = lambda.p) {
  gait.errorcheck(alter, inflate, truncate, max.support)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (!length(c(alter, inflate, truncate)) &&
      is.infinite(max.support)) {
    return(qpois(p, lambda.p))
  }

  if (min(pobs.a, na.rm = TRUE) < 0 || 1 < max(pobs.a, na.rm = TRUE))
    stop("argument 'pobs.a' out of range")
  if (min(pstr.i, na.rm = TRUE) < 0 || 1 < max(pstr.i, na.rm = TRUE))
    stop("argument 'pstr.i' out of range")

  LLL <- max(length(p), length(lambda.p), length(lambda.a),
             length(lambda.i), length(pobs.a), length(pstr.i))
  if (length(p)        != LLL) p        <- rep_len(p,        LLL)
  if (length(lambda.p) != LLL) lambda.p <- rep_len(lambda.p, LLL)
  if (length(lambda.a) != LLL) lambda.a <- rep_len(lambda.a, LLL)
  if (length(lambda.i) != LLL) lambda.i <- rep_len(lambda.i, LLL)
  if (length(pobs.a)   != LLL) pobs.a   <- rep_len(pobs.a,   LLL)
  if (length(pstr.i)   != LLL) pstr.i   <- rep_len(pstr.i,   LLL)

  min.support <- 0  # Usual case
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else
    min.support
  ans <- p + lambda.p + lambda.a + lambda.i + pobs.a + pstr.i

  bad0 <- !is.finite(lambda.p) | lambda.p <= 0 |
          !is.finite(lambda.a) | lambda.a <= 0 |
          !is.finite(lambda.i) | lambda.i <= 0 |
          (lalter &
         (!is.finite(pobs.a)  | pobs.a  <= 0 | 1 <= pobs.a)) |
          (linfla &
         (!is.finite(pstr.i)  | pstr.i  <= 0 | 1 <= pstr.i))
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <- if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else 2 * lo + 10.5
  dont.iterate <- bad   # | is.finite(max.support)
  done <- dont.iterate |
          p <= pgaitpois.mix(hi, lambda.p = lambda.p,
                               pobs.a  = pobs.a, pstr.i  = pstr.i,
                               lambda.a = lambda.a,
                               lambda.i = lambda.i,
                               alter = alter, inflate = inflate,
                               truncate = truncate,
                               max.support = max.support)
  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    lo[!done] <- hi[!done]
    hi[!done] <- 2 * hi[!done] + 10.5  # Bug fixed
    hi <- pmin(max.support + 0.5, hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitpois.mix(hi[!done],
                                   lambda.p = lambda.p[!done],
                                   lambda.a = lambda.a[!done],
                                   lambda.i = lambda.i[!done],
                                   pobs.a  = pobs.a[!done],
                                   pstr.i  = pstr.i[!done],
                                   alter   = alter,
                                   inflate = inflate,
                                   truncate = truncate,
                                   max.support = max.support))
    iter <- iter + 1
  }  # while

      foo <- function(q, lambda.p,
                      pobs.a = 0, pstr.i = 0,
                      lambda.a = lambda.p,
                      lambda.i = lambda.p,
                      alter = NULL,
                      inflate = NULL, truncate = NULL,
                      max.support = Inf, p)
    pgaitpois.mix(q, lambda.p = lambda.p,
                    pobs.a  = pobs.a, pstr.i  = pstr.i,
                    lambda.a = lambda.a, lambda.i = lambda.i,
                    alter = alter,
                    inflate = inflate, truncate = truncate,
                    max.support = max.support) - p
  lhs <- dont.iterate |
         p <= dgaitpois.mix(min.support.use,
                              lambda.p = lambda.p,
                              lambda.a = lambda.a,
                              lambda.i = lambda.i,
                              pobs.a   = pobs.a,
                              pstr.i   = pstr.i,
                              alter    = alter,
                              inflate  = inflate,
                              truncate = truncate,
                              max.support = max.support)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                      pobs.a   = pobs.a[!lhs],
                      pstr.i   = pstr.i[!lhs],
                      lambda.p = lambda.p[!lhs],
                      lambda.a = lambda.a[!lhs],
                      lambda.i = lambda.i[!lhs],
                      alter    = alter,
                      inflate  = inflate, truncate = truncate,
                      max.support = max.support,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitpois.mix(faa, lambda.p = lambda.p[!lhs],
                             pobs.a = pobs.a[!lhs],
                             pstr.i = pstr.i[!lhs],
                             lambda.a = lambda.a[!lhs],
                             lambda.i = lambda.i[!lhs],
                             alter = alter,
                             inflate = inflate,
                             truncate = truncate,
                             max.support = max.support) < p[!lhs] &
             p[!lhs] <= pgaitpois.mix(faa+1,
                                        lambda.p = lambda.p[!lhs],
                                        pobs.a = pobs.a[!lhs],
                                        pstr.i = pstr.i[!lhs],
                                        lambda.a = lambda.a[!lhs],
                                        lambda.i = lambda.i[!lhs],
                                        alter = alter,
                                        inflate = inflate,
                                        truncate = truncate,
                                        max.support = max.support),
             faa+1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]


  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitpois.mix(min.support.use,
                                lambda.p = lambda.p,
                                pobs.a = pobs.a,
                                pstr.i = pstr.i,
                                lambda.a = lambda.a,
                                lambda.i = lambda.i,
                                alter = alter,
                                inflate = inflate,
                                truncate = truncate,
                                max.support = max.support)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitpois.mix






rgaitpois.mix <-
  function(n, lambda.p,
           alter = NULL, inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           lambda.a = lambda.p, lambda.i = lambda.p) {
    qgaitpois.mix(runif(n), lambda.p = lambda.p,
                    pobs.a = pobs.a, pstr.i = pstr.i,
                    lambda.a = lambda.a,
                    lambda.i = lambda.i,
                    alter = alter, inflate = inflate,
                    truncate = truncate,
                    max.support = max.support)
}  # rgaitpois.mix





EIM.gatpoisson <-
  function(lambda, sumderiv1, sumderiv2,
           suma, sumt, Sume, SumT,
           lhs.prob, onempobs.a) {
  ned2l.dlambda2 <- onempobs.a *
    ((lambda - Sume - SumT) / ((lhs.prob - suma - sumt) *
                               lambda^2) -
     sumderiv2 / (lhs.prob - suma - sumt) -
    (sumderiv1 / (lhs.prob - suma - sumt))^2)
  ned2l.dlambda2
}  # EIM.gatpoisson






 gatpoisson.mix <-
  function(alter = NULL,  # Must be assigned a (>=2)-vector to run
           truncate = NULL, max.support = Inf,  # Optional
           zero = "pobs.a",
           parallel = FALSE,  # TRUE applies to the intercept
           lpobs.a = "logitlink",
           llambda.p = "loglink",
           llambda.a = "loglink",
           type.fitted = c("mean", "pobs.a", "Pobs.a",
                           "prob.a", "prob.t", "lhs.prob"),
           imethod = 1,
           ilambda.p = NULL,
           ilambda.a = NULL, ishrinkage = 0.95,
           probs.y = 0.35) {
  gait.errorcheck(alter, inflate = NULL, truncate, max.support)
  
  lpobs.a <- as.list(substitute(lpobs.a))
  epobs.a <- link2list(lpobs.a)
  lpobs.a <- attr(epobs.a, "function.name")
  llambda.p <- as.list(substitute(llambda.p))
  elambda.p <- link2list(llambda.p)
  llambda.p <- attr(elambda.p, "function.name")
  llambda.a <- as.list(substitute(llambda.a))
  elambda.a <- link2list(llambda.a)
  llambda.a <- attr(elambda.a, "function.name")
  lalter <- length(alter)

  type.fitted <- match.arg(type.fitted,
    c("mean", "pobs.a", "Pobs.a", "prob.a", "prob.t", "lhs.prob"))[1]
  tmp3 <- c(pobs.a   = lpobs.a,
            lambda.p = llambda.p,
            lambda.a = llambda.a)

  new("vglmff",
  blurb = c("Generally-altered and -truncated Poisson distribution\n",
            "(GAT-Pois-Pois mixture)\n\n",
        "Links:    ",
        namesof("pobs.a",    lpobs.a,   epobs.a,   tag = FALSE), ", ",
        namesof("lambda.p",  llambda.p, elambda.p, tag = FALSE), ", ",
        namesof("lambda.a",  llambda.a, elambda.a, tag = FALSE), "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(c(1, 0, 0, 0, 1, 1), 3, 2), x = x,
                           bool = .parallel ,
                           constraints = constraints,
                           apply.int = TRUE)  # FALSE
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .zero = zero,
            .parallel = parallel ))),

  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         link = .tmp3 ,
         link1parameter = TRUE,
         mixture.links = FALSE,
         alter = .alter ,
         truncate = .truncate ,
         max.support = .max.support , 
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE,
         parameters.names = c("pobs.a", "lambda.p", "lambda.a"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted,
           .alter = alter,
           .truncate = truncate,
           .max.support = max.support, 
           .tmp3 = tmp3
         ))),

  initialize = eval(substitute(expression({
    truncate <- as.vector( .truncate )
    alter <- as.vector( .alter )
    lalter <- length(alter)
    if (lalter == 0)
      stop("argument 'alter' must not be NULL or empty")
    M1 <- 3
    NOS <- NCOL(y)  # Only 1 currently
    M <- NOS * M1

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (length(truncate) && any(y %in% .truncate ))
      stop("some response values == values in argument 'truncate'")
    if ( .max.support < max(y))
      stop("some response values are greater than the ",
           "'max.support' argument")

    y0 <- matrix(0, n, lalter)
    for (jay in seq(lalter))
      y0[, jay] <- as.numeric(y == alter[jay])
    extra$skip.these <- matrix(as.logical(y0), n, lalter)  # dim lost
      if (any((css <- colSums(extra$skip.these)) == 0))
        stop("some 'alter' argument values have no response values: ",
             paste(alter[css == 0], collapse = ", "))          

    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- ( .type.fitted )
    extra$colnames.y  <- colnames(y)
    extra$M1 <- M1


    mynames1 <- param.names("pobs.a",   ncoly, skip1 = TRUE)
    mynames2 <- param.names("lambda.p", ncoly, skip1 = TRUE)
    mynames3 <- param.names("lambda.a", ncoly, skip1 = TRUE)
    predictors.names <-
    c(namesof(mynames1, .lpobs.a   , earg = .epobs.a   , tag = FALSE),
      namesof(mynames2, .llambda.p , earg = .elambda.p , tag = FALSE),
      namesof(mynames3, .llambda.a , earg = .elambda.a , tag = FALSE))

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                             imu = .ilambda.p ,  # x = x,
                             ishrinkage = .ishrinkage ,
                             pos.only = TRUE,
                             probs.y = .probs.y )
      pobs.a.init <- rep(sum(css) / n, n)  # MLE for pobs.a (unwted)
      lamb.a.init <- if (length( .ilambda.a ))
        rep( .ilambda.a , n) else lambda.init
      
      etastart <-
        cbind(theta2eta(pobs.a.init, .lpobs.a   , earg = .epobs.a   ),
              theta2eta(lambda.init, .llambda.p , earg = .elambda.p ),
              theta2eta(lamb.a.init, .llambda.a , earg = .elambda.a ))
    }
  }), list(
  .llambda.p = llambda.p, .llambda.a = llambda.a, .lpobs.a = lpobs.a,
  .elambda.p = elambda.p, .elambda.a = elambda.a, .epobs.a = epobs.a,
  .ilambda.p = ilambda.p, .ilambda.a = ilambda.a,
  .ishrinkage = ishrinkage, .probs.y = probs.y,
  .imethod = imethod,
  .alter = alter,
  .truncate = truncate,
  .max.support = max.support, 
  .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <-
     if (length(extra$type.fitted)) extra$type.fitted else {
       warning("cannot find 'type.fitted'. Returning the 'mean'.")
                     "mean"
     }

    type.fitted <- match.arg(type.fitted,
      c("mean", "pobs.a", "Pobs.a", "prob.a", "prob.t",
        "lhs.prob"))[1]
    alter <- as.vector( .alter )
    lalter <- length(alter)
    M1 <- 3
    if ((NOS <- NCOL(eta) / M1) != 1)
      stop("Currently NOS must be 1")
    pobs.a   <- eta2theta(eta[, 1], .lpobs.a   , earg = .epobs.a   )
    lambda.p <- eta2theta(eta[, 2], .llambda.p , earg = .elambda.p )
    lambda.a <- eta2theta(eta[, 3], .llambda.a , earg = .elambda.a )
    truncate <- as.vector( .truncate )
    max.support <- as.vector( .max.support )

    bits <- moments.pois.gait(lambda.p = lambda.p, alter = alter,
                              truncate = truncate, mlm = FALSE,
                              max.support = max.support,
                              pobs.a = pobs.a, lambda.a = lambda.a)
    if (type.fitted == "Pobs.a") {
      mymat <-
        dpois(matrix(alter, NROW(eta), lalter, byrow = TRUE),
              matrix(lambda.a, NROW(eta), lalter)) / c(bits[["suma.a"]])
    }
  
    ans <- switch(type.fitted,
                  "mean"     = bits[["mean"]],
                  "pobs.a"   = pobs.a,
                  "Pobs.a"   = pobs.a * mymat,  # matrix
                  "prob.a"   = bits[["suma.p"]],
                  "prob.t"   = bits[["sumt.p"]],
                  "lhs.prob" = bits[["lhs.prob"]])
   ans <-
   label.cols.y(ans,
                 colnames.y = if (type.fitted  == "Pobs.a")
                              as.character(alter) else extra$colnames.y,
                 NOS = NOS)
   ans
  }, list(
  .llambda.p = llambda.p, .llambda.a = llambda.a, .lpobs.a = lpobs.a,
  .elambda.p = elambda.p, .elambda.a = elambda.a, .epobs.a = epobs.a,
  .truncate = truncate,
  .max.support = max.support, 
  .alter = alter ))),

  last = eval(substitute(expression({
    misc$link  <- c(pobs.a   = .lpobs.a   ,
                    lambda.p = .llambda.p ,
                    lambda.a = .llambda.a )

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- .epobs.a  #
    misc$earg[[2]] <- .elambda.p  #
    misc$earg[[3]] <- .elambda.a  #
  }), list(
  .llambda.p = llambda.p, .llambda.a = llambda.a, .lpobs.a = lpobs.a,
  .elambda.p = elambda.p, .elambda.a = elambda.a, .epobs.a = epobs.a,
  .alter = alter ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    pobs.a   <- eta2theta(eta[, 1], .lpobs.a   , earg = .epobs.a   )
    lambda.p <- eta2theta(eta[, 2], .llambda.p , earg = .elambda.p )
    lambda.a <- eta2theta(eta[, 3], .llambda.a , earg = .elambda.a )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitpois.mix(y, lambda.p = lambda.p, lambda.a = lambda.a,
                        pobs.a = pobs.a, truncate = .truncate ,
                        max.support = .max.support ,
                        alter = .alter , log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(
  .llambda.p = llambda.p, .llambda.a = llambda.a, .lpobs.a = lpobs.a,
  .elambda.p = elambda.p, .elambda.a = elambda.a, .epobs.a = epobs.a,
  .truncate = truncate,
  .max.support = max.support, 
  .alter = alter ))),
  vfamily = c("gatpoisson.mix"),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pobs.a   <- eta2theta(eta[, 1], .lpobs.a   , earg = .epobs.a   )
    lambda.p <- eta2theta(eta[, 2], .llambda.p , earg = .elambda.p )
    lambda.a <- eta2theta(eta[, 3], .llambda.a , earg = .elambda.a )
    okay1 <- all(is.finite(lambda.p)) && all(0 < lambda.p) &&
             all(is.finite(lambda.a)) && all(0 < lambda.a) &&
             all(is.finite(pobs.a)) && all(0 < pobs.a & pobs.a < 1)
    okay1
  }, list(
  .llambda.p = llambda.p, .llambda.a = llambda.a, .lpobs.a = lpobs.a,
  .elambda.p = elambda.p, .elambda.a = elambda.a, .epobs.a = epobs.a,
  .truncate = truncate, .max.support = max.support, 
  .alter = alter ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs.a   <- eta2theta(eta[, 1], .lpobs.a   , earg = .epobs.a   )
    lambda.p <- eta2theta(eta[, 2], .llambda.p , earg = .elambda.p )
    lambda.a <- eta2theta(eta[, 3], .llambda.a , earg = .elambda.a )
    rgaitpois.mix(nsim * length(lambda.p),
                    lambda.a = lambda.a, lambda.p = lambda.p,
                    pobs.a = pobs.a,
                    alter = .alter , truncate = .truncate ,
                    max.support = .max.support )
  }, list(
  .llambda.p = llambda.p, .llambda.a = llambda.a, .lpobs.a = lpobs.a,
  .elambda.p = elambda.p, .elambda.a = elambda.a, .epobs.a = epobs.a,
  .truncate = truncate, .max.support = max.support,
  .alter = alter ))),
  deriv = eval(substitute(expression({
    alter <- as.vector( .alter )
    lalter <- length(alter)
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    M1 <- 3
    NOS <- NCOL(eta) / M1  # extra$NOS
    if (lalter) {
      is.altered <- rowSums(extra$skip.these) > 0
    }
    
    pobs.a   <- eta2theta(eta[, 1], .lpobs.a   , earg = .epobs.a   )
    lambda.p <- eta2theta(eta[, 2], .llambda.p , earg = .elambda.p )
    lambda.a <- eta2theta(eta[, 3], .llambda.a , earg = .elambda.a )
    max.support <- ( .max.support )
    bits <- moments.pois.gait(lambda.p,
                              alter = alter, truncate = truncate,
                              max.support = max.support,
                              pobs.a = pobs.a,
                              lambda.a = lambda.a, mlm = FALSE)

    pmf.deriv1 <- function(y, lambda)
      dpois(y-1, lambda) - dpois(y, lambda)
    pmf.deriv2 <- function(y, lambda)
      dpois(y-2, lambda) - 2 * dpois(y-1, lambda) + dpois(y, lambda)
    sumderiv1a.a <- sumderiv2a.a <-
    sumderiv1a.p <- sumderiv2a.p <- matrix(0, NROW(eta), NOS)
    deriv0matrix <-
    deriv1matrix <- deriv2matrix <- matrix(0, NROW(eta), lalter)
    if (lalter) {
      for (jay in seq(lalter)) {
        aval <- alter[jay]
        sumderiv1a.p <- sumderiv1a.p + pmf.deriv1(aval, lambda.p)
        sumderiv2a.p <- sumderiv2a.p + pmf.deriv2(aval, lambda.p)
        sumderiv1a.a <- sumderiv1a.a + pmf.deriv1(aval, lambda.a)
        sumderiv2a.a <- sumderiv2a.a + pmf.deriv2(aval, lambda.a)
        pmf.a <- dpois(aval, lambda.a)
        deriv0matrix[, jay] <- pmf.a
        deriv1matrix[, jay] <- pmf.deriv1(aval, lambda.a) / pmf.a
        deriv2matrix[, jay] <- pmf.deriv2(aval, lambda.a) / pmf.a
      }
      deriv0matrix <-  deriv0matrix / rowSums(deriv0matrix)
    }  # lalter
      
    sumderiv1t.a <- sumderiv2t.a <-
    sumderiv1t.p <- sumderiv2t.p <- matrix(0, NROW(eta), NOS)
    if (ltrunc)
      for (tval in truncate) {
        sumderiv1t.p <- sumderiv1t.p + pmf.deriv1(tval, lambda.p)
        sumderiv2t.p <- sumderiv2t.p + pmf.deriv2(tval, lambda.p)
        sumderiv1t.a <- sumderiv1t.a + pmf.deriv1(tval, lambda.a)
        sumderiv2t.a <- sumderiv2t.a + pmf.deriv2(tval, lambda.a)
      }

    sumderiv1t.p <- sumderiv1t.p + dpois( .max.support    , lambda.p)
    sumderiv2t.p <- sumderiv2t.p + dpois( .max.support - 1, lambda.p) -
                                   dpois( .max.support    , lambda.p)
    sumderiv1t.a <- sumderiv1t.a + dpois( .max.support    , lambda.a)
    sumderiv2t.a <- sumderiv2t.a + dpois( .max.support - 1, lambda.a) -
                                   dpois( .max.support    , lambda.a)

    onempobs.a <- 1 - pobs.a  # 1 - Pr(altered)
    dl.dlambda.p <- dl.dlambda.a <- zero0n <- rep(0, n)
    dl.dpobs.a <- ifelse(is.altered, 1 / pobs.a, -1 / onempobs.a) 
    for (jay in seq(lalter)) {
      aval <- alter[jay]
      dl.dlambda.a <- dl.dlambda.a +
        ifelse(extra$skip.these[, jay],
               deriv1matrix[, jay] - sumderiv1a.a / bits[["suma.a"]],
               zero0n)
    }  # jay

    Denom <- bits[["lhs.prob"]] - bits[["suma.p"]] - bits[["sumt.p"]]
    dl.dlambda.p <-
      ifelse(is.altered,
             zero0n,
             y / lambda.p - 1 +
             (sumderiv1a.p + sumderiv1t.p) / Denom)

    dpobs.a.deta   <- dtheta.deta(pobs.a  , .lpobs.a   , .epobs.a   )
    dlambda.p.deta <- dtheta.deta(lambda.p, .llambda.p , .elambda.p )
    dlambda.a.deta <- dtheta.deta(lambda.a, .llambda.a , .elambda.a )
    ans <- cbind(dl.dpobs.a   * dpobs.a.deta,
                 dl.dlambda.p * dlambda.p.deta,
                 dl.dlambda.a * dlambda.a.deta)
    c(w) * ans
  }), list(
  .llambda.p = llambda.p, .llambda.a = llambda.a, .lpobs.a = lpobs.a,
  .elambda.p = elambda.p, .elambda.a = elambda.a, .epobs.a = epobs.a,
  .truncate = truncate, .max.support = max.support, 
  .alter = alter ))),
  weight = eval(substitute(expression({

    cond.EY.p <-
      (lambda.p - bits[["SumA.p"]] - bits[["SumT.p"]]) / Denom
    
    ned2l.dlambda.p2 <- onempobs.a * (cond.EY.p / lambda.p^2 -
      (sumderiv2a.p + sumderiv2t.p) / Denom -  # '-' is correct
     ((sumderiv1a.p + sumderiv1t.p) / Denom)^2)

    ned2l.dlambda.a2 <- pobs.a *
      rowSums(deriv0matrix * (
             -deriv2matrix +
              deriv1matrix^2 +
             (c(sumderiv2a.a) / c(bits[["suma.a"]])) -
             (c(sumderiv1a.a) / c(bits[["suma.a"]]))^2))

    wz4 <- if ( .lpobs.a == "logitlink") {
      pobs.a * (1 - pobs.a)
    } else {
      dpobs.a.deta <- dtheta.deta(pobs.a, .lpobs.a , .epobs.a )
      ned2l.dpobs.a2 <- 1 / (pobs.a * (1 - pobs.a))
      ned2l.dpobs.a2 * dpobs.a.deta^2
    }
    wz <- wz.merge(wz4,
                   cbind(ned2l.dlambda.p2 * dlambda.p.deta^2,
                         ned2l.dlambda.a2 * dlambda.a.deta^2),
                   M1 = 1, M2 = 2)
    c(w) * wz
  }), list( .lpobs.a = lpobs.a,
            .epobs.a = epobs.a))))
}  # gatpoisson.mix








 gait.errorcheck <-
  function(alter = NULL,
           inflate = NULL,
           truncate = NULL,
           max.support = Inf,
           min.support = 0) {
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)

  if (!is.numeric(max.support) || is.na(max.support) ||
      length(max.support) != 1 || max.support < min.support ||
      round(max.support) != max.support ||
      (length(truncate) && (
          min(truncate, na.rm = TRUE) < min.support ||
          max.support <= max(truncate, na.rm = TRUE))))
    stop("bad input for argument 'max.support' and/or ",
         "'truncate'")

  bothargs <- c(alter, inflate)
  allargs <- c(bothargs, truncate)
  if (lalter + linfla)
    if (!is.Numeric(bothargs, integer.valued = TRUE) ||
        any(bothargs < min.support) ||
        any(max.support < bothargs))
      stop("bad input for arguments 'alter' and/or 'inflate'")
  if (length(unique(allargs)) < lalter + linfla + ltrunc)
      stop("duplicate values found in arguments 'alter', ",
           "'inflate' and 'truncate'")
}  # gait.errorcheck






dgaitnbinom.mlm <-
  function(x, size, prob = NULL, munb = NULL,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE,
           log.arg = FALSE) {
  gait.errorcheck(alter, inflate, truncate, max.support)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (lalter + linfla + ltrunc == 0 && is.infinite(max.support))
    return(if (length(prob))
           dnbinom(x, size, prob = prob, log = log.arg) else
           dnbinom(x, size, mu   = munb, log = log.arg))


  if (any(pobs.a < 0 | 1 <= pobs.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.a'")
  if (any(pstr.i < 0 | 1 <= pstr.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(x), length(size), length(prob), length(munb))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(size)   != LLL) size   <- rep_len(size,   LLL)
  if (length(prob) &&
      length(prob)   != LLL) prob   <- rep_len(prob,   LLL)
  if (length(munb) &&
      length(munb)   != LLL) munb   <- rep_len(munb,   LLL)

  sumt <- 0  # Initialization to 0 important
  if (ltrunc)
    for (tval in truncate)
      sumt <- sumt + (if (length(prob))
              dnbinom(tval, size, prob = prob) else
              dnbinom(tval, size, mu   = munb))
  vecTF.t <- is.finite(x) & ((x %in% truncate) | (max.support < x))
  lhs.prob <- if (length(prob))
              pnbinom(max.support, size, prob = prob) else
              pnbinom(max.support, size, mu   = munb)
  denom.t <- lhs.prob - sumt  # No sumt on RHS

  if (log.arg) {
    logpmf <- ifelse(vecTF.t, log(0),
      if (length(prob))
        dnbinom(x, size, prob = prob, log = TRUE) - log(denom.t) else
        dnbinom(x, size, mu   = munb, log = TRUE) - log(denom.t))
  } else {  # dgtbinom
    pmf <- ifelse(vecTF.t, 0,
                  if (length(prob))
                    dnbinom(x, size, prob = prob) / denom.t else
                    dnbinom(x, size, mu   = munb) / denom.t)
  }

  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (lalter) {
    pobs.a <-  matrix(pobs.a, LLL, lalter, byrow = byrow.arg)
    sum.a <- .rowSums(pobs.a, LLL, lalter)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.a'")  # zz

    for (aval in alter)
      suma <- suma + (if (length(prob))
                      dnbinom(aval, size, prob = prob) else
                      dnbinom(aval, size, mu   = munb))

    for (jay in seq(lalter)) {
      aval <- alter[jay]
      if (any(vecTF <- is.finite(x) & aval == x)) {
        if (log.arg) {
          logpmf[vecTF] <- log(pobs.a[vecTF, jay])
        } else {
          pmf[vecTF] <- pobs.a[vecTF, jay]
        }
      }
      vecTF.a <- vecTF.a | vecTF
    }  # jay
  }  # lalter


  sum.i <- 0
  if (linfla) {
    pstr.i <-  matrix(pstr.i, LLL, linfla, byrow = byrow.arg)
    sum.i <- .rowSums(pstr.i, LLL, linfla)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.i'")
  }

  skip <- vecTF.t | vecTF.a  # Leave these alone
  if (log.arg) {
    logpmf[!skip] <- (log1p(-sum.a - sum.i) + (if (length(prob))
      dnbinom(x, size, prob = prob, log = TRUE) -
      log(lhs.prob - suma - sumt) else
      dnbinom(x, size, mu   = munb, log = TRUE) -
      log(lhs.prob - suma - sumt)))[!skip]
  } else {
      pmf[!skip] <- ((1 - sum.a - sum.i) * (if (length(prob))
        dnbinom(x, size, prob = prob) else
        dnbinom(x, size, mu   = munb)
        ) / (lhs.prob - suma - sumt))[!skip]
  }

  if (linfla) {
    for (jay in seq(linfla)) {
      ival <- inflate[jay]
      if (any(vecTF <- is.finite(x) & ival == x)) {
        if (log.arg) {
          logpmf[vecTF] <- log(exp(logpmf[vecTF]) + pstr.i[vecTF, jay])
        } else {
          pmf[vecTF] <- pmf[vecTF] + pstr.i[vecTF, jay]
        }
      }
    }  # jay
  }  # linfla

  if (log.arg) logpmf else pmf
}  # dgaitnbinom.mlm






pgaitnbinom.mlm <-
  function(q, size, prob = NULL, munb = NULL,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE) {
  gait.errorcheck(alter, inflate, truncate, max.support)

  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (lalter + linfla + ltrunc == 0 &&
      is.infinite(max.support) && 0 < max.support)
    return(if (length(prob))
           pnbinom(q, size, prob = prob) else
           pnbinom(q, size, mu   = munb))  # lower.tail, log.p


  if (any(pobs.a < 0 | 1 <= pobs.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.a'")
  if (any(pstr.i < 0 | 1 <= pstr.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(q), length(size), length(prob), length(munb))
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(size)   != LLL) size   <- rep_len(size,   LLL)
  if (length(prob) &&
      length(prob)   != LLL) prob   <- rep_len(prob,   LLL)
  if (length(munb) &&
      length(munb)   != LLL) munb   <- rep_len(munb,   LLL)


  sumt <- 0
  fudge.t <- numeric(LLL)
  lhs.prob <- if (length(prob))
              pnbinom(max.support, size, prob = prob) else
              pnbinom(max.support, size, mu   = munb)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      local.pmf <- if (length(prob))
              dnbinom(tval, size, prob = prob) else
              dnbinom(tval, size, mu   = munb)
      sumt <- sumt + local.pmf
      if (any(vecTF <- is.finite(q) & tval <= q))
        fudge.t[vecTF] <- fudge.t[vecTF] + local.pmf[vecTF]
    }
  }  # ltrunc

  sum.a <- suma <- 0  # numeric(LLL)
  fudge.a <- numeric(LLL)
  offset.a <- numeric(LLL)
  if (lalter) {
    pobs.a <-  matrix(pobs.a, LLL, lalter, byrow = byrow.arg)
    sum.a <- .rowSums(pobs.a, LLL, lalter)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.a'")

    for (jay in seq(lalter)) {
      aval <- alter[jay]
      local.pmf <- if (length(prob))
                     dnbinom(aval, size, prob = prob) else
                     dnbinom(aval, size, mu   = munb)
      suma <- suma + local.pmf
      if (any(vecTF <- is.finite(q) & aval <= q)) {
        offset.a[vecTF] <- offset.a[vecTF] + pobs.a[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + local.pmf[vecTF]
      }
    }  # jay
  }  # lalter

  sum.i <- 0
  offset.i <- numeric(LLL)
  if (linfla) {
    pstr.i <-  matrix(pstr.i, LLL, linfla, byrow = byrow.arg)
    sum.i <- .rowSums(pstr.i, LLL, linfla)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.i'")

    for (jay in seq(linfla)) {
      ival <- inflate[jay]
      local.pmf <- if (length(prob))
                     dnbinom(ival, size, prob = prob) else
                     dnbinom(ival, size, mu   = munb)
      if (any(vecTF <- is.finite(q) & ival <= q)) {
        offset.i[vecTF] <- offset.i[vecTF] + pstr.i[vecTF, jay]
      }
    }  # jay
  }  # linfla

  numer1 <- 1 - sum.i - sum.a
  denom1 <- lhs.prob - sumt - suma
  ans <- numer1 * ((if (length(prob))
         pnbinom(q, size, prob = prob) else
         pnbinom(q, size, mu   = munb)) - fudge.t -
         fudge.a) / denom1 + offset.i + offset.a

  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  ans
}  # pgaitnbinom.mlm






qgaitnbinom.mlm <-
  function(p, size, prob = NULL, munb = NULL,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE) {
  gait.errorcheck(alter, inflate, truncate, max.support)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (lalter + linfla + ltrunc == 0 && is.infinite(max.support))
    return(if (length(prob))
           qnbinom(p, size, prob = prob) else
           qnbinom(p, size, mu   = munb))  # lower.tail, log.p


  if (any(pobs.a < 0 | 1 <= pobs.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.a'")
  if (any(pstr.i < 0 | 1 <= pstr.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(p), length(size), length(prob), length(munb))
  if (length(p)      != LLL) p      <- rep_len(p,      LLL)
  if (length(size)   != LLL) size   <- rep_len(size,   LLL)
  if (length(prob) &&
      length(prob)   != LLL) prob   <- rep_len(prob,   LLL)
  if (length(munb) &&
      length(munb)   != LLL) munb   <- rep_len(munb,   LLL)

  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)
  pstr.i <- matrix(pstr.i, LLL, linfla, byrow = byrow.arg)

  min.support <- 0  # Usual case
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else min.support
  ans <- p + size + (if (length(prob)) prob else munb)

  bad0 <- !is.finite(size) | size <= 0
  if (length(prob))
    bad0 <- bad0 | !is.finite(prob) | prob <= 0 | 1 <= prob
  if (length(munb))
    bad0 <- bad0 | !is.finite(munb) | munb <= 0
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <- if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else 2 * lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgaitnbinom.mlm(hi, size = size, prob = prob,
                           munb = munb, alter = alter,
                           inflate = inflate, truncate = truncate,
                           pstr.i = pstr.i,
                           pobs.a = pobs.a, byrow.arg = FALSE,
                           max.support = max.support)
  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    lo[!done] <- hi[!done]
    hi[!done] <- 2 * hi[!done] + 10.5  # Bug fixed
    hi <- pmin(max.support, hi) + 0.5  # 20190921
    done[!done] <-
      (p[!done] <= pgaitnbinom.mlm(hi[!done], size[!done],
               prob = if (length(prob)) prob[!done] else NULL,
               munb = if (length(munb)) munb[!done] else NULL,
                               alter = alter,
                               inflate = inflate, truncate = truncate,
                               pstr.i = pstr.i[!done, , drop = FALSE],
                               pobs.a = pobs.a[!done, , drop = FALSE],
                               byrow.arg = FALSE,
                               max.support = max.support))
    iter <- iter + 1
  }  # while

      foo <- function(q, size, prob = NULL, munb = NULL,
                      alter = NULL,
                      inflate = NULL, truncate = NULL,
                      pstr.i = 0,
                      pobs.a = 0, byrow.arg = FALSE,
                      max.support = Inf, p)
    pgaitnbinom.mlm(q, size = size, prob = prob, munb = munb,
                  alter = alter,
                  inflate = inflate, truncate = truncate,
                  pstr.i = pstr.i,
                  pobs.a = pobs.a, byrow.arg = FALSE,
                  max.support = max.support) - p
  lhs <- dont.iterate |
      p <= dgaitnbinom.mlm(min.support.use, size = size,
                         prob = prob, munb = munb,
                         alter = alter,
                         inflate = inflate, truncate = truncate,
                         pstr.i = pstr.i,
                         pobs.a = pobs.a, byrow.arg = FALSE,
                         max.support = max.support)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                      size = size[!lhs],
               prob = if (length(prob)) prob[!lhs] else NULL,
               munb = if (length(munb)) munb[!lhs] else NULL,
                      alter = alter,
                      inflate = inflate, truncate = truncate,
                      pstr.i = pstr.i[!lhs, , drop = FALSE],
                      pobs.a = pobs.a[!lhs, , drop = FALSE],
                      byrow.arg = FALSE,
                      max.support = max.support,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitnbinom.mlm(faa, size[!lhs],
               prob = if (length(prob)) prob[!lhs] else NULL,
               munb = if (length(munb)) munb[!lhs] else NULL,
                         alter = alter,
                         inflate = inflate, truncate = truncate,
                         pstr.i = pstr.i[!lhs, , drop = FALSE],
                         pobs.a = pobs.a[!lhs, , drop = FALSE],
                         byrow.arg = FALSE,
                         max.support = max.support) < p[!lhs] &
             p[!lhs] <= pgaitnbinom.mlm(faa+1, size[!lhs],
               prob = if (length(prob)) prob[!lhs] else NULL,
               munb = if (length(munb)) munb[!lhs] else NULL,
                      alter = alter,
                      inflate = inflate, truncate = truncate,
                      pstr.i = pstr.i[!lhs, , drop = FALSE],
                      pobs.a = pobs.a[!lhs, , drop = FALSE],
                      byrow.arg = FALSE,
                      max.support = max.support),
             faa+1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)



  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]


  vecTF <- !bad0 & !is.na(p) &
      p <= dgaitnbinom.mlm(min.support.use, size = size,
                         prob = prob,
                         munb = munb,
                         alter = alter,
                         inflate = inflate, truncate = truncate,
                         pstr.i = pstr.i, pobs.a = pobs.a,
                         byrow.arg = FALSE,
                         max.support = max.support)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitnbinom.mlm





rgaitnbinom.mlm <-
  function(n, size = NULL, prob = NULL, munb = NULL,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE) {
  qgaitnbinom.mlm(runif(n), size = size, prob = prob, munb = munb,
                alter = alter, inflate = inflate,
                truncate = truncate,
                pobs.a = pobs.a, pstr.i = pstr.i,
                byrow.arg = byrow.arg,
                max.support = max.support)
}  # rgaitnbinom.mlm






GATNB.deriv012 <-
  function(munb, size,
           alter = NULL, truncate = NULL, max.support = Inf) {
  if (is.finite(max.support))
    stop("can only handle finite 'max.support'")
  lalter <- length(alter)
  ltrunc <- length(truncate)
  sumderiv1.munb <- sumderiv1.size <-
  sumderiv2.munb <- sumderiv2.size <- sumderiv2.both <-
  suma <- SumA <- matrix(0, NROW(munb), NCOL(munb))
  if (lalter + ltrunc)
    for (aval in c(alter, truncate)) {
      pmf <- cbind(dnbinom(aval, size, mu = munb))
      suma <- suma + pmf
      SumA <- SumA + pmf * aval
      dl.dmunb <- aval / munb - (1 + aval / size) / (1 + munb / size)
      dl.dsize <- digamma(aval + size) - digamma(size) -
        (aval - munb) / (munb + size) +
        log1p(-munb / (munb + size))
      d2l.dmunb2 <- (1 + aval / size) / (munb / sqrt(size) +
        sqrt(size))^2 - aval / munb^2
      d2l.dsize2 <- trigamma(aval + size) - trigamma(size) +
        (aval - munb) / (munb + size)^2 +
        munb / (size * (munb + size))
      d2l.dmunbsize <- (aval - munb) / (munb + size)^2
      sumderiv1.munb <- sumderiv1.munb + dl.dmunb * pmf
      sumderiv1.size <- sumderiv1.size + dl.dsize * pmf
      sumderiv2.munb <- sumderiv2.munb + (d2l.dmunb2 + dl.dmunb^2) * pmf
      sumderiv2.size <- sumderiv2.size + (d2l.dsize2 + dl.dsize^2) * pmf
      sumderiv2.both <- sumderiv2.both +
        (d2l.dmunbsize + dl.dmunb * dl.dsize) * pmf
    }  # for aval
  list('sumderiv0'      = suma,  # + sumt,
       'sumEat'         = SumA,  # + SumT,
       'sumderiv1.munb' = sumderiv1.munb,
       'sumderiv2.munb' = sumderiv2.munb, 
       'sumderiv1.size' = sumderiv1.size,
       'sumderiv2.size' = sumderiv2.size,
       'sumderiv2.both' = sumderiv2.both
      )
}  # GATNB.deriv012






EIM.GATNB.speciald <-
  function(munb, size,  # == munb.p, size.p, respectively
           munb.a = munb, size.a = size,
           alter = NULL, truncate = NULL, max.support = Inf,
           pobs.a = 0,
           EY.cond = NULL,  # = munb would be a good default zz
           y.min = 0,  # 20160201; must be an integer
           y.max = NULL,  # Must be an integer
           cutoff.prob = 0.999,
           intercept.only = FALSE,
           mlm = TRUE,  # 20191106; Added, now handles 2 variants
           extra.bit = TRUE) {


  n.orig <- length(munb)
  if (intercept.only) {
    munb <- munb[1]
    size <- size[1]
    pobs.a <- if (mlm) {
      if (is.matrix(pobs.a)) pobs.a[1, ] else  # 20191002
        stop("confused on what should be a matrix: pobs.a")
    } else {
      pobs.a[1]  # 20191106
    }
  }

  if (!is.numeric(y.max)) {
    eff.p <- sort(c(cutoff.prob, 1 - cutoff.prob))
    y.max <- if (mlm)
      max(round(qgaitnbinom.mlm(p = eff.p[2],
        size = size, alter = alter,
        truncate = truncate, max.support = max.support,
        pobs.a = pobs.a, byrow.arg = intercept.only,
        munb = munb) * 1.1)) + 30 else
      max(round(qgaitnbinom.mix(p = eff.p[2],
        munb.p = munb, munb.a = munb.a, 
        size.p = size, size.a = size.a,
        truncate = truncate, max.support = max.support,
        alter = alter, pobs.a = pobs.a) * 1.1)) + 30
  }

  Y.mat <- if (intercept.only) y.min:y.max else
    matrix(y.min:y.max, length(munb), y.max-y.min+1, byrow = TRUE)
 
  trigg.term <- if (intercept.only) {  # 1 x 1 so c() needed below
    if (mlm)
      dgaitnbinom.mlm(Y.mat, size = size, munb = munb,
                        alter = NULL, truncate = c(alter, truncate),
                        max.support = max.support,
                        byrow.arg = intercept.only,
                        pobs.a = pobs.a) %*%
      trigamma(Y.mat + size) else
      dgaitnbinom.mix(Y.mat,
                        munb.p = munb, munb.a = munb.a, 
                        size.p = size, size.a = size.a,
                        alter = NULL, truncate = c(alter, truncate),
                        max.support = max.support,
                        pobs.a = pobs.a) %*%
      trigamma(Y.mat + size)
  } else {
    if (mlm)
      .rowSums(dgaitnbinom.mlm(Y.mat, size = size, munb = munb,
                                 alter = NULL,
                                 truncate = c(alter, truncate),
                                 max.support = max.support,
                                 byrow.arg = intercept.only,
                                 pobs.a = pobs.a) *
                                 trigamma(Y.mat + c(size)),
             NROW(Y.mat), NCOL(Y.mat)) else
      .rowSums(dgaitnbinom.mix(Y.mat,
                                 munb.p = munb, munb.a = munb.a, 
                                 size.p = size, size.a = size.a,
                                 alter = NULL,
                                 truncate = c(alter, truncate),
                                 max.support = max.support,
                                 pobs.a = pobs.a) *
                                 trigamma(Y.mat + c(size)),
             NROW(Y.mat), NCOL(Y.mat))
  }
  ned2l.dk2 <- trigamma(size) - c(trigg.term)  # 1 x 1

  if (extra.bit)
    ned2l.dk2 <- ned2l.dk2 - munb / (size * (size + munb)) -
                 (EY.cond - munb) / (munb + size)^2
                 
  if (intercept.only)
    matrix(ned2l.dk2, n.orig, 1) else as.matrix(ned2l.dk2)
}  # EIM.GATNB.speciald













gatnbinomial.mlm.control <-
  function(save.weights = TRUE,
           summary.HDEtest = FALSE,  # Overwrites summary() default.
           ...) {
  list(save.weights = save.weights,
       summary.HDEtest = summary.HDEtest)
}


 gatnbinomial.mlm <-
  function(alter = NULL,
           truncate = NULL,  # max.support = Inf,
           zero = "size",
           lmunb = "loglink", lsize = "loglink",
           type.fitted = c("mean", "pobs.a",
               "Pobs.a", "prob.a", "prob.t"),
           imethod = 1,
           imunb = NULL, isize = exp(1), ishrinkage = 0.95,
           probs.y = 0.35,

           cutoff.prob = 0.999,  # higher is better for large 'size'
           eps.trig = 1e-7,
           max.chunk.MB = 30  # max.memory = Inf is allowed
           ) {
  max.support <- Inf  # Fixed for now
  lmunb <- as.list(substitute(lmunb))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")

  gait.errorcheck(alter, inflate = NULL, truncate, max.support)

  lalter <- length(alter)
  ltrunc <- length(truncate)
  type.fitted <- match.arg(type.fitted, c("mean", "pobs.a",
                   "Pobs.a", "prob.a", "prob.t"))[1]
  temp7 <- if (lalter) paste0("pobs", alter) else NULL
  tmp3 <- c(if (lalter) rep("multilogitlink", lalter) else NULL,
            munb = lmunb, size = lsize)
  names(tmp3) <- c(temp7, "munb", "size")
  blurb1 <- "N"
  if (lalter) blurb1 <- "Generally-altered n"
  if (ltrunc) blurb1 <- "Generally-truncated n"
  if (lalter && ltrunc) blurb1 <- "Generally-altered and -truncated n"
                
  new("vglmff",
  blurb = c(blurb1, "egative binomial distribution ",
            "(GAT-NB-MLM in general)\n\n",
            "Links:   ", if (lalter)
            paste("multilogitlink(cbind(",
            paste(temp7, collapse = ", "), ", ",
            "\n                              ",
            "1-", paste(temp7, collapse = "-"),
            "))", ", ", sep = ""),
            if (length(alter)) "\n         ",
            namesof("munb", lmunb, earg = emunb, tag = FALSE), ", ",
            namesof("size", lsize, earg = esize, tag = FALSE)),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = length( .alter ) + 2)
  }), list( .zero = zero,
            .alter = alter ))),

  infos = eval(substitute(function(...) {
    alter <- ( .alter )
    temp7 <- paste0("pobs", alter)
    list(M1 = length( .alter ) + 2,
         Q1 = 1,
         link = .tmp3 ,  # multilogitlink is multiparameter:
         link1parameter = !as.logical(length( .alter )), 
         mixture.links = TRUE,
         alter = .alter ,
         truncate = .truncate ,
         max.support = .max.support , 
         expected = TRUE,
         multipleResponses = !as.logical(length( .alter ) +
                                         length( .truncate )),
         parameters.names = c(if (length( .alter )) temp7 else NULL,
                              "munb", "size"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted,
           .lmunb = lmunb, .emunb = emunb,
           .alter = alter,
           .truncate = truncate,
           .max.support = max.support, 
           .tmp3 = tmp3
         ))),

  initialize = eval(substitute(expression({
    truncate <- ( .truncate )
    ltrunc <- length(truncate)
    alter <- ( .alter )
    lalter <- length(alter)
    M1 <- lalter + 2
    NOS <- NCOL(y)
    M <- NOS * M1

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = if (lalter + ltrunc) 1 else Inf,
              ncol.y.max = if (lalter + ltrunc) 1 else Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (ltrunc && any(y %in% truncate))
      stop("some response values == values in argument 'truncate'")
    if ( .max.support < max(y))
      stop("some response values > than argument 'max.support'")

    if (lalter) {  # Memory hungry
      extra$y0 <- y0 <- matrix(0, n, lalter)
      for (jay in seq(lalter))
        extra$y0[, jay] <- y0[, jay] <- as.numeric(y == alter[jay])
      extra$skip.these <-
            skip.these <- matrix(as.logical(y0), n, lalter)  # dim lost
      if (any((css <- colSums(skip.these)) == 0))
        stop("some 'alter' argument values have no response values: ",
             paste(alter[css == 0], collapse = ", "))          
    }
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)  # May be NULL
    extra$M1 <- M1


    mynames1 <- if (lalter) {
      temp7 <- paste0("pobs", alter)
      denom.char <- paste0("1-", paste0(temp7, collapse = "-"))
      paste("log(", temp7, "/(", denom.char, "))", sep = "")
    } else {
      NULL
    }
    mynames2 <- param.names("munb", NOS, skip1 = TRUE)
    mynames3 <- param.names("size", NOS, skip1 = TRUE)
    predictors.names <-
      c(        mynames1,
        namesof(mynames2, .lmunb , earg = .emunb , tag = FALSE),
        namesof(mynames3, .lsize , earg = .esize , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {
      munb.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                           imu = .imunb ,  # x = x,
                           ishrinkage = .ishrinkage ,
                           probs.y = .probs.y )
      size.init <- matrix( .isize , n, NCOL(munb.init))
      if (lalter) {
        phimat <- colMeans(skip.these)  # yettodo: weight this by 'w'
        phimat <- matrix(phimat, n, lalter, byrow = TRUE)
        etastart <- multilogitlink(cbind(phimat,
                                         abs(1 - rowSums(phimat))))
      }  # lalter
      etastart <-
        cbind(etastart,
              theta2eta(munb.init, .lmunb , earg = .emunb ),
              theta2eta(size.init, .lsize , earg = .esize ))[,
        interleave.VGAM(M, M1 = M1)]
    }
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize,
            .imunb = imunb, .isize = isize,
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .imethod = imethod,
            .alter = alter,
            .truncate = truncate, .max.support = max.support, 
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <-
      if (length(extra$type.fitted)) extra$type.fitted else {
        warning("cannot find 'type.fitted'. Returning the 'mean'.")
        "mean"  # Unconditional mean
      }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "pobs.a",
                       "Pobs.a", "prob.a", "prob.t"))[1]
    alter <- ( .alter )
    lalter <- length(alter)
    M1 <- length(alter) + 2
    NOS <- NCOL(eta) / M1
    munb <- eta2theta(eta[, (1:extra$NOS)*extra$M1 - 1, drop = FALSE],
                      .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, (1:extra$NOS)*extra$M1    , drop = FALSE],
                      .lsize , earg = .esize )
    truncate <- ( .truncate )
    ltrunc <- length(truncate)
    max.support <- ( .max.support )
    lhs.prob <- pnbinom(max.support, kmat, mu = munb)  # Usually 1
    sumt <- matrix(0, NROW(eta), NOS)
    SumT <- matrix(0, NROW(eta), NOS)  # zz rmlife needed here
    if (ltrunc)
      for (tval in truncate) {
        pmf <- dnbinom(tval, kmat, mu = munb)
        sumt <- sumt + pmf  # Need tval<=max.support
        SumT <- SumT + pmf * tval
      }  # for tval
    iprd <- numeric(NROW(eta))  # iprd is an innerprod
    suma <- SumA <- matrix(0, NROW(eta), NOS)
    pobs.a <- 0
    phimat <- NULL  # Nonsensical otherwise
    if (lalter) {
      index3 <- 1:(NCOL(eta) - 2)
      phimat <- multilogitlink(eta[, index3, drop = FALSE],
                refLevel = NCOL(eta) - 1,  # Assumes one response
                inverse = TRUE)
      ynames.Pobs.a <- c(as.character(alter), "(Others)")
      dimnames(phimat) <- list(rownames(eta), ynames.Pobs.a)
      pobs.a <- rowSums(phimat[, -NCOL(phimat), drop = FALSE])
      for (jay in seq_len(lalter)) {
        aval <- alter[jay]
        pmf <- dnbinom(aval, kmat, mu = munb)
        suma <- suma + pmf
        SumA <- SumA + pmf * aval
        iprd <- iprd + phimat[, jay] * aval
      }  # for jay
    }  # lalter

    ans <- switch(type.fitted,
      "mean"       = iprd + (1 - pobs.a) *
                     (munb - SumA - SumT) / (lhs.prob - suma - sumt),
      "pobs.a"     =     pobs.a,  # Pr(Y is altered)
      "Pobs.a"     = phimat,  # matrix
      "prob.a"     = suma,
      "prob.t"     = sumt)
    ans <-
    label.cols.y(ans,
                 colnames.y = if (type.fitted  == "Pobs.a")
                              ynames.Pobs.a else extra$colnames.y,
                 NOS = NOS)
    ans
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize,
           .alter = alter,
           .truncate = truncate, .max.support = max.support ))),

  last = eval(substitute(expression({
    tmp9 <- if (lalter) rep_len("multilogitlink", lalter) else NULL
    temp.names <- c(tmp9, rep_len( .lmunb , NOS), rep_len( .lsize , NOS))
    misc$link  <- temp.names
    names(misc$link) <- c(mynames1, mynames2, mynames3)[
      interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)

    if (lalter) {
      for (ii in seq(M1*NOS - 2)) {
        misc$earg[[ii]] <- list(M = M - 2,  # M * NOS,
                                refLevel = M - 1)  # M * NOS
      }
      misc$earg[[M1*NOS - 1]] <- .emunb  #
      misc$earg[[M1*NOS    ]] <- .esize  # Last one
    } else {
      for (ii in seq(NOS)) {
        misc$earg[[M1*ii - 1]] <- .emunb  #
        misc$earg[[M1*ii    ]] <- .esize  # Last one
      }
    }
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize,
            .alter = alter ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    lalter <- length( .alter )
    index3 <- 1:(NCOL(eta) - 2)
    if (length( .alter ))
      pobs.a <- multilogitlink(eta[, index3, drop = FALSE],
                  refLevel = NCOL(eta) - 1,  # Assumes 1 response
                  inverse = TRUE)
    munb <- eta2theta(eta[, (1:extra$NOS)*extra$M1 - 1, drop = FALSE],
                      .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, (1:extra$NOS)*extra$M1    , drop = FALSE],
                      .lsize , earg = .esize )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitnbinom.mlm(y, kmat, munb = munb, log = TRUE,
                      truncate = .truncate ,
                      max.support = .max.support ,
                      alter = .alter ,  # byrow.arg = FALSE,
                      pobs.a = if (length( .alter ))
                      pobs.a[, -NCOL(pobs.a), drop = FALSE] else 0)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize,
           .alter = alter,
           .truncate = truncate, .max.support = max.support ))),
  vfamily = c("gatnbinomial.mlm"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    index3 <- 1:(NCOL(eta) - 2)
    pobs.a <- if (length( .alter ))
        multilogitlink(eta[, index3, drop = FALSE],
                refLevel = NCOL(eta) - 1,  # Assumes one response
                inverse = TRUE) else 0  # An okay value
    munb <- eta2theta(eta[, (1:extra$NOS)*extra$M1 - 1, drop = FALSE],
                      .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, (1:extra$NOS)*extra$M1    , drop = FALSE],
                      .lsize , earg = .esize )
    okay1 <- all(is.finite(munb))   && all(0 <  munb) &&
             all(is.finite(kmat))   && all(0 <  kmat) &&
             all(is.finite(pobs.a)) && all(0 <= pobs.a & pobs.a < 1)
    okay1
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize,
           .alter = alter,
           .truncate = truncate, .max.support = max.support ))),
  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    extra <- object@extra

    index3 <- 1:(NCOL(eta) - 2)
    if (length( .alter ))
      pobs.a <- multilogitlink(eta[, index3, drop = FALSE],
                  refLevel = NCOL(eta) - 1,  # Assumes 1 response
                  inverse = TRUE) else 0  # An okay value
    munb <- eta2theta(eta[, (1:extra$NOS)*extra$M1 - 1, drop = FALSE],
                      .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, (1:extra$NOS)*extra$M1    , drop = FALSE],
                      .lsize , earg = .esize )
    rgaitnbinom.mlm(nsim * length(munb), kmat, munb = munb,
                  alter = .alter ,
                  pobs.a = if (length( .alter ))
                    pobs.a[, -NCOL(pobs.a), drop = FALSE] else 0,
                  truncate = .truncate ,
                  max.support = .max.support )
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize,
           .alter = alter,
           .truncate = truncate, .max.support = max.support ))),
  deriv = eval(substitute(expression({
    alter <- ( .alter )
    lalter <- length(alter)
    truncate <- ( .truncate )
    ltrunc <- length(truncate)
    M1 <- lalter + 2
    NOS <- NCOL(eta) / M1  # extra$NOS
    min.support <- 0  # Usual case
    min.support.use <- if (ltrunc)
      min(setdiff(min.support:(ltrunc+5), truncate)) else min.support
    max.support <- ( .max.support )
    if (lalter) {
      y0 <- extra$y0
      skip <- extra$skip.these
      is.altered <- rowSums(skip) > 0  # TRUE if any(y %in% avec)
    }

    index3 <- 1:(NCOL(eta) - 2)
    phimat <-
    pobs.a <- if (length( .alter ))
      multilogitlink(eta[, index3, drop = FALSE],
                refLevel = NCOL(eta) - 1,  # Assumes 1 response
                inverse = TRUE) else 0  # An okay value
    munb <- eta2theta(eta[, (1:extra$NOS)*extra$M1 - 1, drop = FALSE],
                      .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, (1:extra$NOS)*extra$M1    , drop = FALSE],
                      .lsize , earg = .esize )
    lhs.prob <- pgaitnbinom.mlm(max.support, kmat, munb = munb)

    sdlist <- GATNB.deriv012(munb, size = kmat,
                             alter, truncate, max.support)

    onempobs.a <- if (lalter)
      phimat[, ncol(phimat), drop = FALSE] else 1
    pobs.a <- 1 - onempobs.a

    
    dl.deta <- if (lalter) {
      skip - phimat[, -NCOL(phimat), drop = FALSE]
    } else NULL


    dl.dmunb <- (if (lalter) (1 - is.altered) else 1) *
      (y / munb - (1 + y / kmat) / (1 + munb / kmat) +
      sdlist$sumderiv1.munb / (lhs.prob - sdlist$sumderiv0))
    dl.dsize <- (if (lalter) (1 - is.altered) else 1) *
      (digamma(y + kmat) - digamma(kmat) -
      (y - munb) / (munb + kmat) +
      log1p(-munb / (munb + kmat)) +
      sdlist$sumderiv1.size / (lhs.prob - sdlist$sumderiv0))

    dmunb.deta <- dtheta.deta(munb, .lmunb , earg = .emunb )
    dsize.deta <- dtheta.deta(kmat, .lsize , earg = .esize )
    ans <- (c(w) * cbind(dl.deta,
                         dl.dmunb * dmunb.deta,
                         dl.dsize * dsize.deta))[,
                   interleave.VGAM(M1*NOS, M1 = M1)]
    ans
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize,
            .alter = alter,
            .truncate = truncate, .max.support = max.support ))),
  weight = eval(substitute(expression({
    max.chunk.MB <- ( .max.chunk.MB )
    MM12 <- M1 * (M1 + 1) / 2
    wz <- matrix(0, n, MM12)  # A full matrix
    EY.cond <- (munb - sdlist$sumEat) / (lhs.prob - sdlist$sumderiv0)

    ned2l.dmunb2 <-
      (EY.cond * (1 / munb^2 - 1 / (munb + kmat)^2) -
       1 / (munb / sqrt(kmat) + sqrt(kmat))^2) -
         sdlist$sumderiv2.munb / (lhs.prob - sdlist$sumderiv0) -
        (sdlist$sumderiv1.munb / (lhs.prob - sdlist$sumderiv0))^2


    ned2l.dmunbsize <- -(EY.cond - munb) / (munb + kmat)^2 -
       sdlist$sumderiv2.both / (lhs.prob - sdlist$sumderiv0) -
       sdlist$sumderiv1.munb *
       sdlist$sumderiv1.size / (lhs.prob - sdlist$sumderiv0)^2

    ned2l.dsize2 <- matrix(0, n, NOS)
    
    ind2 <- matrix(FALSE, n, NOS)  # Used for SFS
    for (jay in 1:NOS) {
      eff.p <- sort(c( .cutoff.prob , 1 - .cutoff.prob ))
      Q.mins  <- min.support.use  # 1
      Q.mins2 <- pmax(qgaitnbinom.mlm(eff.p[1],
                                    size = kmat[, jay],
                                    munb = munb[, jay],
                                    alter = alter,
                                    truncate = truncate,
                                    max.support = max.support) - 2,
                      Q.mins)
      Q.maxs <-  qgaitnbinom.mlm(p           = eff.p[2] ,
                               size        = kmat[, jay],
                               munb        = munb[, jay],
                               alter       = .alter ,
                               truncate    = .truncate ,
                               max.support = .max.support ) + 10
      eps.trig <- .eps.trig
      Q.MAXS <- pmax(10, ceiling(1 / sqrt(eps.trig)))
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
          eim.gatnb <-
          EIM.GATNB.speciald(munb = munb[sind2, jay],
            size        = kmat[sind2, jay],
            alter = alter, truncate = truncate,
            max.support = max.support,
            pobs.a      = if (lalter)
            phimat[sind2, -NCOL(phimat), drop = FALSE] else pobs.a,
            EY.cond     = EY.cond[sind2, jay],
            y.min       = min.support,  # min(Q.mins   ),
            y.max       = max(Q.maxs[sind2]),
            cutoff.prob = .cutoff.prob ,
            intercept.only = intercept.only,
            extra.bit   = TRUE)
          ned2l.dsize2[sind2, jay] <- eim.gatnb

          if (FALSE)
          if (any(eim.kk.TF <-       wz[sind2, M1*jay] <= 0 |
                               is.na(wz[sind2, M1*jay]))) {
            ind2[sind2[eim.kk.TF], jay] <- FALSE
          }
          lwr.ptr <- upr.ptr + 1
        }  # while (lwr.ptr <= NN)
      }  # if ((NN <- sum(ind1)) > 0)
    }  # end of for (jay in 1:NOS)




    ned2l.dsize2 <- ned2l.dsize2 -
 ( 1)*    sdlist$sumderiv2.size / (lhs.prob - sdlist$sumderiv0) -
 ( 1)*   (sdlist$sumderiv1.size / (lhs.prob - sdlist$sumderiv0))^2





    if (lalter) {
      wz4 <- matrix(0, n, MM12)  # A full matrix
      wz4 <- matrix(0, n, lalter*(lalter+1)/2)
      wz <- if (lalter > 0) {
        if (lalter == 1) {
          wz4[, 1] <- phimat[, 1] * (1 - phimat[, 1])
        } else {
          index <- iam(NA, NA, lalter, both = TRUE, diag = TRUE)
          wz4 <- -phimat[, index$row] * phimat[, index$col]
          wz4[, 1:lalter] <- wz4[, 1:lalter] + phimat[, -NCOL(phimat)]
        }
 try1 = c(w) * cbind(ned2l.dmunb2 * dmunb.deta^2,
                     ned2l.dsize2 * dsize.deta^2,
                     ned2l.dmunbsize * dmunb.deta *
                       dsize.deta) * c(onempobs.a)
        wz.merge(c(w) * wz4,
                 c(w) * cbind(ned2l.dmunb2 * dmunb.deta^2,
                              ned2l.dsize2 * dsize.deta^2,
                              ned2l.dmunbsize * dmunb.deta *
                              dsize.deta) * c(onempobs.a),
                 M1 = lalter, M2 = 2)  # rm.trailing.cols = F
      }  # lalter > 0
    } else {
      wz <- array(c(c(w) * onempobs.a * ned2l.dmunb2 * dmunb.deta^2,
                    c(w) * onempobs.a * ned2l.dsize2 * dsize.deta^2,
                    c(w) * onempobs.a * ned2l.dmunbsize *
                           dmunb.deta * dsize.deta),
                   dim = c(n, NOS, 3))
      wz <- arwz2wz(wz, M = M, M1 = M1, full.arg = TRUE)
    }
    wz
  }), list( .alter = alter,
            .truncate = truncate, .max.support = max.support,
            .cutoff.prob = cutoff.prob, .eps.trig = eps.trig,
            .max.chunk.MB = max.chunk.MB
           ))))
}  # gatnbinomial.mlm








 gatpoisson.mlm <-
  function(alter = NULL,
           truncate = NULL, max.support = Inf,
           zero = NULL,
           llambda = "loglink",
           type.fitted = c("mean", "lambda", "pobs.a",
               "Pobs.a", "prob.a", "prob.t"),
           imethod = 1,
           ilambda = NULL, ishrinkage = 0.95,
           probs.y = 0.35) {
  gait.errorcheck(alter, inflate = NULL, truncate, max.support)
  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")
  lalter <- length(alter)

  type.fitted <- match.arg(type.fitted,
    c("mean", "lambda", "pobs.a", "Pobs.a", "prob.a", "prob.t"))[1]
  temp7 <- if (lalter) paste0("pobs", alter) else NULL
  tmp3 <- c(if (lalter)
            rep("multilogitlink",  lalter) else NULL,
            lambda = llambda )
  names(tmp3) <- c(temp7, "lambda")


  new("vglmff",
  blurb = c("Generally-altered and -truncated Poisson distribution\n",
            "(GAT-Pois-MLM)",
            "\n\n",
            "Links:    ", if (lalter)
      paste("multilogitlink(cbind(",
            paste(temp7, collapse = ", "), ", ",
            "\n                           ",
            "1 - ", paste(temp7, collapse = " - "),
            "))", ", ", sep = ""),
            if (length(alter)) "\n          ",
            namesof("lambda", llambda, earg = elambda, tag = FALSE)),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = length( .alter ) + 1)
  }), list( .zero = zero,
            .alter = alter ))),

  infos = eval(substitute(function(...) {
    alter <- ( .alter )
    temp7 <- paste("pobs", alter, sep = "")
    list(M1 = length( .alter ) + 1,
         Q1 = 1,
         link = .tmp3 ,
         link1parameter = if (length( .alter ))
               FALSE else TRUE,  # multilogitlink is multiparameter
         mixture.links = TRUE,
         alter = .alter ,
         truncate = .truncate ,
         max.support = .max.support , 
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE,
         parameters.names = c(if (length( .alter )) temp7 else NULL,
                              "lambda"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted,
           .llambda = llambda, .elambda = elambda,
           .alter = alter,
           .truncate = truncate,
           .max.support = max.support, 
           .tmp3 = tmp3
         ))),

  initialize = eval(substitute(expression({
    truncate <- ( .truncate )
    ltrunc <- length(truncate)
    alter <- ( .alter )
    lalter <- length(alter)
    M1 <- lalter + 1
    NOS <- NCOL(y)
    M <- NOS * M1

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = if (lalter + ltrunc) 1 else Inf,
              ncol.y.max = if (lalter + ltrunc) 1 else Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (ltrunc && any(y %in% .truncate ))
      stop("some response values equal values of the 'truncate' ",
           "argument")
    if ( .max.support < max(y))
      stop("some response values are greater than the ",
           "'max.support' argument")

    if (lalter) {
      extra$y0 <- y0 <- matrix(0, n, lalter)
      for (jay in seq(lalter))
        extra$y0[, jay] <- y0[, jay] <- as.numeric(y %in% alter[jay])
      extra$skip.these <- skip.these <- matrix(as.logical(y0), n, lalter)
      if (any((css <- colSums(skip.these)) == 0))
        stop("some 'alter' argument values have no response values: ",
             paste(alter[css == 0], collapse = ", "))          
    }
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)
    extra$M1 <- M1


    if (lalter) {
      temp7 <- paste0("pobs", alter)
      denom.char <- paste0("1-", paste0(temp7, collapse = "-"))
      mynames1 <- paste("log(", temp7, "/(",
                        denom.char, "))", sep = "")
    } else {
      mynames1 <- NULL
    }
    mynames2 <- param.names("lambda", ncoly, skip1 = TRUE)
    predictors.names <-
        c(        mynames1,
          namesof(mynames2, .llambda , earg = .elambda ,
                  tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                             imu = .ilambda ,  # x = x,
                             ishrinkage = .ishrinkage ,
                             pos.only = TRUE,
                             probs.y = .probs.y )
      if (lalter) {
        phimat <- colMeans(skip.these)  # yettodo: weight this by 'w'
        phimat <- matrix(phimat, n, lalter, byrow = TRUE)
        etastart <-  multilogitlink(cbind(phimat, 1 - rowSums(phimat)))
      }
      etastart <-
        cbind(etastart,
              theta2eta(lambda.init, .llambda , earg = .elambda ))
    }
  }), list( .llambda = llambda,
            .elambda = elambda,
            .ilambda = ilambda,
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .imethod = imethod,
            .alter = alter,
            .truncate = truncate,
            .max.support = max.support, 
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <- if (length(extra$type.fitted))
                   extra$type.fitted else {
                      warning("cannot find 'type.fitted'. ",
                              "Returning the 'mean'.")
                      "mean"
                    }

    type.fitted <- match.arg(type.fitted,
       c("mean", "lambda", "pobs.a", "Pobs.a", "prob.a", "prob.t"))[1]
    eta <- as.matrix(eta)
    alter <- ( .alter )
    lalter <- length(alter)
    M1 <- length(alter) + 1
    NOS <- NCOL(eta) / M1
    lambda <- eta2theta(eta[, (1:NOS)*M1, drop = FALSE],
                        .llambda , earg = .elambda )
    truncate <- ( .truncate )
    ltrunc <- length(truncate)
    max.support <- ( .max.support )
    if (lalter) {
      phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                  refLevel = NCOL(eta),  # Assumes one response
                  inverse = TRUE)
      ynames.Pobs.a <- c(as.character(alter), "(Others)")
      dimnames(phimat) <- list(rownames(eta), ynames.Pobs.a)
      pobs.a <- rowSums(phimat[, -NCOL(eta), drop = FALSE])
    }  # lalter

    bits <- moments.pois.gait(lambda,
                              alter = alter, truncate = truncate,
                              max.support = max.support,
                              pobs.a = phimat[, -ncol(phimat)],
                              lambda.a = lambda, mlm = TRUE)
    ans <- switch(type.fitted,
      "mean"       = bits[["mean"]],
      "lambda"     = lambda,
      "pobs.a"     =     pobs.a,  # Pr(Y is altered)
      "Pobs.a"     =     phimat,  # matrix
      "prob.a"     = bits[["suma.p"]],
      "prob.t"     = bits[["sumt.p"]])

   ans <-
   label.cols.y(ans,
                 colnames.y = if (type.fitted  == "Pobs.a")
                              ynames.Pobs.a else extra$colnames.y,
                 NOS = NOS)
   ans
  }, list( .llambda = llambda,
           .elambda = elambda,
           .truncate = truncate,
           .max.support = max.support, 
           .alter = alter ))),

  last = eval(substitute(expression({
    tmp9 <- if (lalter) rep_len( "multilogitlink" , lalter) else NULL
    temp.names <- c(tmp9, rep_len( .llambda , NOS))
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)

    if (lalter)
      for (ii in seq(M1*NOS - 1)) {
          misc$earg[[ii]] <- list(M = M - 1,  # M * NOS,
                                  refLevel = M)  # M * NOS
      }
    misc$earg[[M1*NOS]] <- .elambda  # Last one
  }), list( .llambda = llambda, .elambda = elambda,
            .alter = alter ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    lalter <- length( .alter )
    eta <- as.matrix(eta)

    if (lalter)
      pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                  refLevel = NCOL(eta),  # Assumes one response
                  inverse = TRUE)
    lambda <- eta2theta(eta[, (1:extra$NOS)*extra$M1, drop = FALSE],
                        .llambda , earg = .elambda )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitpois.mlm(y, lambda = lambda, log = TRUE,
                    truncate = .truncate ,
                    max.support = .max.support ,
                    alter = .alter ,  # byrow.arg = FALSE,
                    pobs.a = if (length( .alter ))
                      pobs.a[, -NCOL(pobs.a), drop = FALSE] else 0)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llambda = llambda, .elambda = elambda,
           .truncate = truncate,
           .max.support = max.support, 
           .alter = alter ))),
  vfamily = c("gatpoisson.mlm"),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pobs.a <- if (length( .alter ))
        multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE) else 0.5  # An okay value
    eta <- as.matrix(eta)
    lambda <- eta2theta(eta[, (1:extra$NOS)*extra$M1, drop = FALSE],
                        .llambda , earg = .elambda )
    okay1 <- all(is.finite(lambda)) && all(0 < lambda) &&
             all(is.finite(pobs.a)) && all(0 < pobs.a & pobs.a < 1)
    okay1
  }, list( .llambda = llambda, .elambda = elambda,
           .truncate = truncate,
           .max.support = max.support, 
           .alter = alter ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs.a <- if (length( .alter ))
        multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE) else 0.5  # Any okay value
    NOS <- npred(object) / npred(object, M1 = TRUE)
    lambda <- eta2theta(eta[, (1:NOS)*M1, drop = FALSE],
                        .llambda , earg = .elambda )
    rgaitpois.mlm(nsim * length(lambda), lambda = lambda,
                pobs.a = pobs.a, alter = .alter ,
                truncate = .truncate ,
                max.support = .max.support )
  }, list( .llambda = llambda, .elambda = elambda,
           .truncate = truncate,
           .max.support = max.support, 
           .alter = alter ))),

  deriv = eval(substitute(expression({
    alter <- ( .alter )
    lalter <- length(alter)
    truncate <- ( .truncate )
    ltrunc <- length(truncate)
    eta <- as.matrix(eta)
    M1 <- lalter + 1
    NOS <- NCOL(eta) / M1  # extra$NOS
    if (lalter) {
      y0 <- extra$y0
      skip <- extra$skip.these
      is.altered <- rowSums(skip) > 0  # TRUE if any(y %in% avec)
    }
    
    phimat <- if (lalter)
        multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                       refLevel = NCOL(eta),  # Assumes 1 response
                       inverse = TRUE) else 0.5
    lambda <- eta2theta(eta[, (1:NOS)*M1, drop = FALSE],
                        .llambda , earg = .elambda )
    lhs.prob <- ppois( .max.support , lambda)  # Usually 1
    sumt <- matrix(0, NROW(eta), NOS)
    SumT <- matrix(0, NROW(eta), NOS)  # rmlife not needed here
    if (ltrunc)
      for (tval in truncate) {
        pmf <- dpois(tval, lambda)
        sumt <- sumt + pmf
        SumT <- SumT + pmf * tval
      }

    SumT <- SumT + lambda * ppois( .max.support - 1,
                                   lambda, lower.tail = FALSE)


    suma <- Sume <- matrix(0, NROW(eta), NOS)
    if (lalter)
      for (aval in alter) {
        pmf <- dpois(aval, lambda)
        suma <- suma + pmf
        Sume <- Sume + pmf * aval
      }
    pmf.deriv1 <- function(y, lambda)
      dpois(y-1, lambda) - dpois(y, lambda)
    pmf.deriv2 <- function(y, lambda)
      dpois(y-2, lambda) - 2 * dpois(y-1, lambda) + dpois(y, lambda)
    sumderiv1 <- sumderiv2 <- matrix(0, NROW(eta), NOS)
    if (lalter)
      for (aval in alter) {
        sumderiv1 <- sumderiv1 + pmf.deriv1(aval, lambda)
        sumderiv2 <- sumderiv2 + pmf.deriv2(aval, lambda)
      }
    if (ltrunc)
      for (tval in truncate) {
        sumderiv1 <- sumderiv1 + pmf.deriv1(tval, lambda)
        sumderiv2 <- sumderiv2 + pmf.deriv2(tval, lambda)
      }

    sumderiv1 <- sumderiv1 + dpois( .max.support    , lambda)
    sumderiv2 <- sumderiv2 + dpois( .max.support - 1, lambda) -
                             dpois( .max.support    , lambda)


    pobs.a <- if (lalter)
      rowSums(phimat[, -NCOL(phimat), drop = FALSE]) else 0

    if (lalter)
      dl.deta <- skip - phimat[, -M, drop = FALSE]
    dl.dlambda <-
      (if (lalter) (1 - is.altered) else 1) *
      (y / lambda - 1 + sumderiv1 / (lhs.prob - suma - sumt))
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    ans <- cbind(if (lalter) c(w) * dl.deta else NULL,
                 c(w) * dl.dlambda * dlambda.deta)
    ans
  }), list( .llambda = llambda, .elambda = elambda,
            .truncate = truncate,
            .max.support = max.support, 
            .alter = alter ))),

  weight = eval(substitute(expression({
    MM12 <- M1 * (M1 + 1) / 2

    ned2l.dlambda2 <- (1 - pobs.a) *
      ((lambda - Sume - SumT) / ((lhs.prob - suma - sumt) *
                                 lambda^2) -
       sumderiv2 / (lhs.prob - suma - sumt) -
      (sumderiv1 / (lhs.prob - suma - sumt))^2)

    if (lalter) {
      wz4 <- matrix(0.0, n, MM12)  # A full matrix
      if (lalter > 0) {
        if (lalter == 1) {
          wz4[, 1] <- phimat[, 1] * (1 - phimat[, 1])
        } else {
          index <- iam(NA, NA, M - 1, both = TRUE, diag = TRUE)
          wz4 <- -phimat[, index$row] * phimat[, index$col]
          wz4[, 1:lalter] <- wz4[, 1:lalter] + phimat[, -M]
        }
      wz <- wz.merge(c(w) * wz4,
                     c(w) * ned2l.dlambda2 * dlambda.deta^2,
                     M1 = M1 - 1, M2 = 1)  # rm.trailing.cols = F
      }  # lalter > 0
    } else {
     wz <- c(w) * ned2l.dlambda2 * dlambda.deta^2
    }
    wz
  }), list( .alter = alter,
            .truncate = truncate,
            .max.support = max.support
           ))))
}  # gatpoisson.mlm






dgaitbinom.mlm <-
  function(x, size, prob,
           alter = NULL,
           inflate = NULL,
           truncate = NULL,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE,
           log.arg = FALSE,
           .errorcheck = TRUE) {
  if ( .errorcheck )
    gait.errorcheck(alter, inflate, truncate,
                    max.support = min(size, na.rm = TRUE))
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (lalter + linfla + ltrunc == 0)
    return(dbinom(x, size, prob, log = log.arg))

  if (any(pobs.a < 0 | 1 <= pobs.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.a'")
  if (any(pstr.i < 0 | 1 <= pstr.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(x), length(size), length(prob))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(size)   != LLL) size   <- rep_len(size,   LLL)
  if (length(prob)   != LLL) prob   <- rep_len(prob,   LLL)

  max.support <- size
  sumt <- 0  # Initialization to 0 important
  if (ltrunc)
    for (tval in truncate)
      sumt <- sumt + dbinom(tval, size, prob)  # Need tval<=max.support
  vecTF.t <- is.finite(x) & (x %in% truncate)
  lhs.prob <- pbinom(max.support, size, prob)  # Usually 1
  denom.t <- lhs.prob - sumt  # No sumt on RHS

  if (log.arg) {
    logpmf <- ifelse(vecTF.t, log(0),
                     dbinom(x, size, prob, log = TRUE) - log(denom.t))
  } else {  # dgtbinom
    pmf <- ifelse(vecTF.t, 0, dbinom(x, size, prob) / denom.t)
  }

  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (lalter) {
    pobs.a <-  matrix(pobs.a, LLL, lalter, byrow = byrow.arg)
    sum.a <- .rowSums(pobs.a, LLL, lalter)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.a'")  # zz

    for (aval in alter)
      suma <- suma + dbinom(aval, size, prob)

    for (jay in seq(lalter)) {
      aval <- alter[jay]
      if (any(vecTF <- is.finite(x) & aval == x)) {
        if (log.arg) {
          logpmf[vecTF] <- log(pobs.a[vecTF, jay])
        } else {
          pmf[vecTF] <- pobs.a[vecTF, jay]
        }
      }
      vecTF.a <- vecTF.a | vecTF
    }  # jay
  }  # lalter


  sum.i <- 0
  if (linfla) {
    pstr.i <-  matrix(pstr.i, LLL, linfla, byrow = byrow.arg)
    sum.i <- .rowSums(pstr.i, LLL, linfla)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.i'")
  }

  skip <- vecTF.t | vecTF.a  # Leave these alone
  if (log.arg) {
    logpmf[!skip] <- (log1p(-sum.a - sum.i) +
                      dbinom(x, size, prob, log = TRUE) -
                      log(lhs.prob - suma - sumt))[!skip]
  } else {
    pmf[!skip] <- ((1 - sum.a - sum.i) * dbinom(x, size, prob) /
                   (lhs.prob - suma - sumt))[!skip]
  }

  if (linfla) {
    for (jay in seq(linfla)) {
      ival <- inflate[jay]
      if (any(vecTF <- is.finite(x) & ival == x)) {
        if (log.arg) {
          logpmf[vecTF] <- log(exp(logpmf[vecTF]) + pstr.i[vecTF, jay])
        } else {
          pmf[vecTF] <- pmf[vecTF] + pstr.i[vecTF, jay]
        }
      }
    }  # jay
  }  # linfla

  if (log.arg) logpmf else pmf
}  # dgaitbinom.mlm






pgaitbinom.mlm <-
  function(q, size, prob,
           alter = NULL,
           inflate = NULL,
           truncate = NULL,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE,
           .errorcheck = TRUE) {
  if ( .errorcheck )
    gait.errorcheck(alter, inflate, truncate,
                    max.support = min(size, na.rm = TRUE))
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (lalter + linfla + ltrunc == 0)
    return(pbinom(q, size, prob))  # lower.tail, log.p


  if (any(pobs.a < 0 | 1 <= pobs.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.a'")
  if (any(pstr.i < 0 | 1 <= pstr.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(q), length(size), length(prob))
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(size)   != LLL) size   <- rep_len(size,   LLL)
  if (length(prob)   != LLL) prob   <- rep_len(prob,   LLL)

  max.support <- size
  sumt <- 0
  fudge.t <- numeric(LLL)
  lhs.prob <- pbinom(max.support, size, prob)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      local.pmf <- dbinom(tval, size, prob)
      sumt <- sumt + local.pmf
      if (any(vecTF <- is.finite(q) & tval <= q))
        fudge.t[vecTF] <- fudge.t[vecTF] + local.pmf[vecTF]
    }
  }  # ltrunc

  sum.a <- suma <- 0  # numeric(LLL)
  fudge.a <- numeric(LLL)
  offset.a <- numeric(LLL)
  if (lalter) {
    pobs.a <-  matrix(pobs.a, LLL, lalter, byrow = byrow.arg)
    sum.a <- .rowSums(pobs.a, LLL, lalter)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.a'")

    for (jay in seq(lalter)) {
      aval <- alter[jay]
      local.pmf <- dbinom(aval, size, prob)
      suma <- suma + local.pmf
      if (any(vecTF <- is.finite(q) & aval <= q)) {
        offset.a[vecTF] <- offset.a[vecTF] + pobs.a[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + local.pmf[vecTF]
      }
    }  # jay
  }  # lalter

  sum.i <- 0
  offset.i <- numeric(LLL)
  if (linfla) {
    pstr.i <-  matrix(pstr.i, LLL, linfla, byrow = byrow.arg)
    sum.i <- .rowSums(pstr.i, LLL, linfla)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.i'")

    for (jay in seq(linfla)) {
      ival <- inflate[jay]
      local.pmf <- dbinom(ival, size, prob)
      if (any(vecTF <- is.finite(q) & ival <= q)) {
        offset.i[vecTF] <- offset.i[vecTF] + pstr.i[vecTF, jay]
      }
    }  # jay
  }  # linfla

  numer1 <- 1 - sum.i - sum.a
  denom1 <- lhs.prob - sumt - suma
  ans <- numer1 * (pbinom(q, size, prob) - fudge.t -
         fudge.a) / denom1 + offset.i + offset.a

  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  ans
}  # pgaitbinom.mlm







qgaitbinom.mlm <-
  function(p, size, prob,
           alter = NULL,
           inflate = NULL,
           truncate = NULL,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE) {
  gait.errorcheck(alter, inflate, truncate,
                  max.support = min(size, na.rm = TRUE))
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (lalter + linfla + ltrunc == 0)
    return(qbinom(p, size, prob))  # lower.tail = TRUE, log.p = FALSE


  if (any(pobs.a < 0 | 1 <= pobs.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.a'")
  if (any(pstr.i < 0 | 1 <= pstr.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(p), length(size), length(prob))
  if (length(p)      != LLL) p      <- rep_len(p,    LLL)
  if (length(size)   != LLL) size   <- rep_len(size, LLL)
  if (length(prob)   != LLL) prob   <- rep_len(prob, LLL)

  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)
  pstr.i <- matrix(pstr.i, LLL, linfla, byrow = byrow.arg)

  min.support <- 0
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else
    min.support
  ans <- p + size + prob

  bad0 <- !is.finite(size) | size <= 0 |
          !is.finite(prob) | prob <= 0 | 1 <= prob
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <- size + 0.5  # 2 * lo + 10.5
  dont.iterate <- bad

  foo <- function(q, size, prob, alter = NULL,
                  inflate = NULL, truncate = NULL,
                  pstr.i = 0,
                  pobs.a = 0, byrow.arg = FALSE,
                  p)
    pgaitbinom.mlm(q, size = size, prob = prob, alter = alter,
                 inflate = inflate, truncate = truncate,
                 pstr.i = pstr.i,
                 pobs.a = pobs.a, byrow.arg = FALSE,
                 .errorcheck = FALSE) - p
  lhs <- dont.iterate |
         p <= dgaitbinom.mlm(min.support.use,
                           size = size, prob = prob,
                           alter = alter,
                           inflate = inflate, truncate = truncate,
                           pstr.i = pstr.i,
                           pobs.a = pobs.a, byrow.arg = FALSE,
                           .errorcheck = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                      size = size[!lhs],
                      prob = prob[!lhs],
                      alter = alter,
                      inflate = inflate, truncate = truncate,
                      pstr.i = pstr.i[!lhs, , drop = FALSE],
                      pobs.a = pobs.a[!lhs, , drop = FALSE],
                      byrow.arg = FALSE,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitbinom.mlm(faa,
                          size = size[!lhs],
                          prob = prob[!lhs],
                          alter = alter,
                          inflate = inflate, truncate = truncate,
                          pstr.i = pstr.i[!lhs, , drop = FALSE],
                          pobs.a = pobs.a[!lhs, , drop = FALSE],
                          byrow.arg = FALSE,
                          .errorcheck = FALSE) < p[!lhs] &
             p[!lhs] <= pgaitbinom.mlm(faa+1,
                                     size = size[!lhs],
                                     prob = prob[!lhs],
                                     alter = alter,
                      inflate = inflate, truncate = truncate,
                      pstr.i = pstr.i[!lhs, , drop = FALSE],
                      pobs.a = pobs.a[!lhs, , drop = FALSE],
                      byrow.arg = FALSE,
                      .errorcheck = FALSE),
             faa+1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitbinom.mlm(min.support.use, size, prob,
                             alter = alter,
                             inflate = inflate, truncate = truncate,
                             pstr.i = pstr.i, pobs.a = pobs.a,
                             byrow.arg = FALSE,
                             .errorcheck = FALSE)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- size[!bad0 & !is.na(p) & p == 1]
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitbinom.mlm






rgaitbinom.mlm <-
  function(n, size, prob,
           alter = NULL,
           inflate = NULL,
           truncate = NULL,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE) {
    qgaitbinom.mlm(runif(n), size = size, prob = prob,
                alter = alter, inflate = inflate,
                truncate = truncate,
                pobs.a = pobs.a, pstr.i = pstr.i,
                byrow.arg = byrow.arg)
}  # rgaitbinom.mlm
















dgaitpois.mlm <-
  function(x, lambda,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE,
           log.arg = FALSE) {
  gait.errorcheck(alter, inflate, truncate, max.support)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (lalter + linfla + ltrunc == 0 && is.infinite(max.support))
    return(dpois(x, lambda, log = log.arg))


  if (any(pobs.a < 0 | 1 <= pobs.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.a'")
  if (any(pstr.i < 0 | 1 <= pstr.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(x), length(lambda))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)

  sumt <- 0  # Initialization to 0 important
  if (ltrunc)
    for (tval in truncate)
      sumt <- sumt + dpois(tval, lambda)  # Need tval <= max.support
  vecTF.t <- is.finite(x) & ((x %in% truncate) | (max.support < x))
  lhs.prob <- ppois(max.support, lambda)  # Usually 1
  denom.t <- lhs.prob - sumt  # No sumt on RHS

  if (log.arg) {
    logpmf <- ifelse(vecTF.t, log(0),
                     dpois(x, lambda, log = TRUE) - log(denom.t))
  } else {
    pmf <- ifelse(vecTF.t, 0, dpois(x, lambda) / denom.t)  # dgtpois
  }

  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (lalter) {
    pobs.a <-  matrix(pobs.a, LLL, lalter, byrow = byrow.arg)
    sum.a <- .rowSums(pobs.a, LLL, lalter)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.a'")  # zz

    for (aval in alter)
      suma <- suma + dpois(aval, lambda)

    for (jay in seq(lalter)) {
      aval <- alter[jay]
      if (any(vecTF <- is.finite(x) & aval == x)) {
        if (log.arg) {
          logpmf[vecTF] <- log(pobs.a[vecTF, jay])
        } else {
          pmf[vecTF] <- pobs.a[vecTF, jay]
        }
      }
      vecTF.a <- vecTF.a | vecTF
    }  # jay
  }  # lalter


  sum.i <- 0
  if (linfla) {
    pstr.i <-  matrix(pstr.i, LLL, linfla, byrow = byrow.arg)
    sum.i <- .rowSums(pstr.i, LLL, linfla)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.i'")
  }

  skip <- vecTF.t | vecTF.a  # Leave these alone
  if (log.arg) {
    logpmf[!skip] <- (log1p(-sum.a - sum.i) +
                      dpois(x, lambda, log = TRUE) -
                      log(lhs.prob - suma - sumt))[!skip]
  } else {
    pmf[!skip] <- ((1 - sum.a - sum.i) * dpois(x, lambda) /
                   (lhs.prob - suma - sumt))[!skip]
  }

  if (linfla) {
    for (jay in seq(linfla)) {
      ival <- inflate[jay]
      if (any(vecTF <- is.finite(x) & ival == x)) {
        if (log.arg) {
          logpmf[vecTF] <- log(exp(logpmf[vecTF]) + pstr.i[vecTF, jay])
        } else {
          pmf[vecTF] <- pmf[vecTF] + pstr.i[vecTF, jay]
        }
      }
    }  # jay
  }  # linfla

  if (log.arg) logpmf else pmf
}  # dgaitpois.mlm






pgaitpois.mlm <-
  function(q, lambda,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE) {
  gait.errorcheck(alter, inflate, truncate, max.support)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (lalter + linfla + ltrunc == 0 && is.infinite(max.support))
    return(ppois(q, lambda))  # lower.tail, log.p

  if (any(pobs.a < 0 | 1 <= pobs.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.a'")
  if (any(pstr.i < 0 | 1 <= pstr.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(q), length(lambda))
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)


  sumt <- 0
  fudge.t <- numeric(LLL)
  lhs.prob <- ppois(max.support, lambda)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      local.pmf <- dpois(tval, lambda)
      sumt <- sumt + local.pmf
      if (any(vecTF <- is.finite(q) & tval <= q))
        fudge.t[vecTF] <- fudge.t[vecTF] + local.pmf[vecTF]
    }
  }  # ltrunc

  sum.a <- suma <- 0  # numeric(LLL)
  fudge.a <- numeric(LLL)
  offset.a <- numeric(LLL)
  if (lalter) {
    pobs.a <-  matrix(pobs.a, LLL, lalter, byrow = byrow.arg)
    sum.a <- .rowSums(pobs.a, LLL, lalter)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.a'")

    for (jay in seq(lalter)) {
      aval <- alter[jay]
      local.pmf <- dpois(aval, lambda)
      suma <- suma + local.pmf
      if (any(vecTF <- is.finite(q) & aval <= q)) {
        offset.a[vecTF] <- offset.a[vecTF] + pobs.a[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + local.pmf[vecTF]
      }
    }  # jay
  }  # lalter

  sum.i <- 0
  offset.i <- numeric(LLL)
  if (linfla) {
    pstr.i <-  matrix(pstr.i, LLL, linfla, byrow = byrow.arg)
    sum.i <- .rowSums(pstr.i, LLL, linfla)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.i'")

    for (jay in seq(linfla)) {
      ival <- inflate[jay]
      local.pmf <- dpois(ival, lambda)
      if (any(vecTF <- is.finite(q) & ival <= q)) {
        offset.i[vecTF] <- offset.i[vecTF] + pstr.i[vecTF, jay]
      }
    }  # jay
  }  # linfla

  numer1 <- 1 - sum.i - sum.a
  denom1 <- lhs.prob - sumt - suma
  ans <- numer1 * (ppois(q, lambda) - fudge.t - fudge.a) / denom1 +
         offset.i + offset.a

  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  ans
}  # pgaitpois.mlm








qgaitpois.mlm <-
  function(p, lambda,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE) {
  gait.errorcheck(alter, inflate, truncate, max.support)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (lalter + linfla + ltrunc == 0 && is.infinite(max.support))
    return(qpois(p, lambda))  # lower.tail = TRUE, log.p = FALSE

  if (any(pobs.a < 0 | 1 <= pobs.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.a'")
  if (any(pstr.i < 0 | 1 <= pstr.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(p), length(lambda))
  if (length(p)      != LLL) p      <- rep_len(p,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)

  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)
  pstr.i <- matrix(pstr.i, LLL, linfla, byrow = byrow.arg)

  min.support <- 0  # Usual case
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else
    min.support
  ans <- p + lambda

  bad0 <- !is.finite(lambda) | lambda <= 0
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <- if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else 2 * lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
          p <= pgaitpois.mlm(hi, lambda, alter = alter,
                           inflate = inflate, truncate = truncate,
                           pstr.i = pstr.i,
                           pobs.a = pobs.a, byrow.arg = FALSE,
                           max.support = max.support)
  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    lo[!done] <- hi[!done]
    hi[!done] <- 2 * hi[!done] + 10.5  # Bug fixed
    hi <- pmin(max.support + 0.5, hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitpois.mlm(hi[!done], lambda[!done],
                               alter = alter,
                               inflate = inflate, truncate = truncate,
                               pstr.i = pstr.i[!done, , drop = FALSE],
                               pobs.a = pobs.a[!done, , drop = FALSE],
                               byrow.arg = FALSE,
                               max.support = max.support))
    iter <- iter + 1
  }

      foo <- function(q, lambda, alter = NULL,
                      inflate = NULL, truncate = NULL,
                      pstr.i = 0,
                      pobs.a = 0, byrow.arg = FALSE,
                      max.support = Inf, p)
    pgaitpois.mlm(q, lambda = lambda, alter = alter,
                inflate = inflate, truncate = truncate,
                pstr.i = pstr.i,
                pobs.a = pobs.a, byrow.arg = FALSE,
                max.support = max.support) - p
  lhs <- dont.iterate |
         p <= dgaitpois.mlm(min.support.use, lambda, alter = alter,
                          inflate = inflate, truncate = truncate,
                          pstr.i = pstr.i,
                          pobs.a = pobs.a, byrow.arg = FALSE,
                          max.support = max.support)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                      lambda = lambda[!lhs], alter = alter,
                      inflate = inflate, truncate = truncate,
                      pstr.i = pstr.i[!lhs, , drop = FALSE],
                      pobs.a = pobs.a[!lhs, , drop = FALSE],
                      byrow.arg = FALSE,
                      max.support = max.support,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitpois.mlm(faa, lambda[!lhs],  alter = alter,
                         inflate = inflate, truncate = truncate,
                         pstr.i = pstr.i[!lhs, , drop = FALSE],
                         pobs.a = pobs.a[!lhs, , drop = FALSE],
                         byrow.arg = FALSE,
                         max.support = max.support) < p[!lhs] &
             p[!lhs] <= pgaitpois.mlm(faa+1, lambda[!lhs],
                                    alter = alter,
                      inflate = inflate, truncate = truncate,
                      pstr.i = pstr.i[!lhs, , drop = FALSE],
                      pobs.a = pobs.a[!lhs, , drop = FALSE],
                      byrow.arg = FALSE,
                      max.support = max.support),
             faa+1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitpois.mlm(min.support.use, lambda, alter = alter,
                            inflate = inflate, truncate = truncate,
                            pstr.i = pstr.i, pobs.a = pobs.a,
                            byrow.arg = FALSE,
                            max.support = max.support)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitpois.mlm






rgaitpois.mlm <-
  function(n, lambda,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE) {
    qgaitpois.mlm(runif(n), lambda = lambda,
                alter = alter, inflate = inflate,
                truncate = truncate,
                pobs.a = pobs.a, pstr.i = pstr.i,
                byrow.arg = byrow.arg,
                max.support = max.support)
}  # rgaitpois.mlm







dgtpois <-
  function(x, lambda, truncate = 0,  # NULL,
           max.support = Inf,  # 20181018
           log = FALSE) {
  if (!is.numeric(max.support) ||
      length(max.support) != 1 ||
      round(max.support) != max.support ||
      (!is.null(truncate) && max.support <= max(truncate)))
    stop("bad input for argument 'max.support'")

  if (is.null(truncate) && is.infinite(max.support) &&
      0 < max.support)
    return(dpois(x, lambda, log = log))
    
  if (is.list(truncate))
    truncate <- unlist(truncate)
    
  if (is.infinite(max.support) &&
      (!is.Numeric(truncate, integer.valued = TRUE) ||
       any(truncate < 0)))
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")

  LLL <- max(length(x), length(lambda))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  sump <- numeric(LLL)
  for (tval in truncate)
    sump <- sump + dpois(tval, lambda)
  vecTF <- (x %in% truncate) | (max.support < x)

  denom <- if (is.finite(max.support)) {
    if (log) log(ppois(max.support, lambda) - sump) else
                 ppois(max.support, lambda) - sump
  } else {
    if (log) log1p(-sump) else 1 - sump
  }

  ans <- if (log) {
    ifelse(vecTF, -Inf,  # log(0),
           dpois(x, lambda, log = TRUE) - denom)
  } else {
    ifelse(vecTF, 0, dpois(x, lambda) / denom)
  }
    
  ans
}  # dgtpois



pgtpois <-
  function(q, lambda, truncate = 0,  # NULL,
           max.support = Inf) {
  if (is.null(truncate) && is.infinite(max.support) &&
      max.support > 0)
    return(ppois(q, lambda))
    
  if (is.list(truncate))
    truncate <- unlist(truncate)
    
  if (is.infinite(max.support) &&
      (!is.Numeric(truncate, integer.valued = TRUE) ||
       any(truncate < 0)))
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")
  if (!is.numeric(max.support) ||
      length(max.support) != 1 ||
      round(max.support) != max.support ||
      (!is.null(truncate) && max.support <= max(truncate)))
    stop("bad input for argument 'max.support'")

  LLL <- max(length(q), length(lambda))
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  sump <- numeric(LLL)
  for (tval in truncate)
    sump <- sump + dpois(tval, lambda)
  denom <- ppois(max.support, lambda) - sump

  numer <- ppois(q, lambda)
  for (tval in truncate) {
    if (any(vecTF <- tval <= q))
      numer[vecTF] <- numer[vecTF] - dpois(tval, lambda[vecTF])
  }
  ans <- numer / denom
  if (is.finite(max.support)) {
    ans[max.support <= q] <- 1
  }
  ans
}  # pgtpois



qgtpois <-
  function(p, lambda, truncate = 0,  # NULL,
           max.support = Inf) {
  if (is.null(truncate) && is.infinite(max.support) &&
      max.support > 0)
    return(qpois(p, lambda))

  if (is.list(truncate))
    truncate <- unlist(truncate)
    
  if (is.infinite(max.support) &&
      (!is.Numeric(truncate, integer.valued = TRUE) ||
       any(truncate < 0)))
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")
  if (!is.numeric(max.support) ||
      length(max.support) != 1 ||
      round(max.support) != max.support ||
      (!is.null(truncate) && max.support <= max(truncate)))
    stop("bad input for argument 'max.support'")

  LLL <- max(length(p), length(lambda))
  if (length(p)      != LLL) p      <- rep_len(p,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)

  min.support <- ifelse(any(truncate == 0), min(truncate+1), 0)

  ans <- p + lambda

  bad0 <- !is.finite(lambda) | lambda <= 0
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(min.support - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <-  if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else 2 * lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgtpois(hi, lambda, truncate = truncate,
                   max.support = max.support)
  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    lo[!done] <- hi[!done]
    hi[!done] <- 2 * hi[!done] + 10.5  # Bug fixed
    hi[!done] <- pmin(max.support + 1.5, hi[!done])
    done[!done] <-
      (p[!done] <= pgtpois(hi[!done], lambda[!done],
                             truncate = truncate,
                             max.support = max.support))
    iter <- iter + 1
  }

  foo <- function(q, lambda, truncate = 0, max.support = Inf, p)
    pgtpois(q, lambda = lambda, truncate = truncate,
              max.support = max.support) - p

  lhs <- dont.iterate |
         p <= dgtpois(min.support, lambda, truncate = truncate,
                        max.support = max.support)
  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                      lambda = lambda[!lhs],
                      truncate = truncate, max.support = max.support,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgtpois(faa, lambda[!lhs], truncate = truncate,
                       max.support = max.support) < p[!lhs] &
             p[!lhs] <= pgtpois(faa+1, lambda[!lhs],
                                  truncate = truncate,
                                  max.support = max.support),
             faa+1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (length(truncate))
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgtpois(min.support, lambda,
                          truncate = truncate,
                          max.support = max.support)
  ans[vecTF] <- min.support  # min.support[vecTF]

  ans[!bad0 & !is.na(p) & p == 0] <- min.support
  ans[!bad0 & !is.na(p) & p == 1] <- max.support
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgtpois



rgtpois <-
  function(n, lambda, truncate = 0,  # NULL
           max.support = Inf) {
    qgtpois(runif(n), lambda = lambda,
              truncate = truncate,
              max.support = max.support)
}  # rgtpois






 gtpoisson <-
  function(truncate = 0,  # NULL,  # 0  #
           zero = NULL,
           max.support = Inf,
           link = "loglink",
           type.fitted = c("mean", "lambda", "prob.t"),
           ilambda = NULL, imethod = 1) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")

  if (length(truncate)) {
    if (!is.Numeric(truncate, integer.valued = TRUE) ||
        any(truncate < 0))
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")
  }
  if (!is.numeric(max.support) ||
      length(max.support) != 1 ||
      round(max.support) != max.support ||
      (!is.null(truncate) && max.support <= max(truncate)))
    stop("bad input for argument 'max.support'")

  if (length( ilambda) && !is.Numeric(ilambda, positive = TRUE))
    stop("bad input for argument 'ilambda'")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "prob.t"))[1]




  new("vglmff",
  blurb = c("Generally-truncated Poisson distribution (GT-Pois)\n\n",
            "Links:    ",
            namesof("lambda", link, earg = earg, tag = FALSE),
            if (length( truncate )) 
            c("\nTruncated at:    ",
            paste( truncate , collapse = ", ")) else NULL,
            "\nMaximum support value:    ",
            max.support),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         max.support = .max.support , 
         multipleResponses = TRUE,
         parameters.names = c("lambda"),
         link = .link ,
         truncate = .truncate ,
         type.fitted  = .type.fitted ,
         earg = .earg )
  }, list( .link = link, .earg = earg,
           .max.support = max.support, 
           .truncate = truncate,
           .type.fitted = type.fitted ))),

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
    if (any(y %in% .truncate ))
      stop("some response values equal values of the 'truncate' ",
           "argument")
    if ( .max.support < max(y))
      stop("some response values are greater than the 'max.support' ",
           "argument")

    ncoly <- ncol(y)
    lengthtruncate <- length( .truncate )
    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("lambda", ncoly, skip1 = TRUE)
    predictors.names <- namesof(mynames1, .link , earg = .earg,
                                tag = FALSE)

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                             imu = .ilambda )

      etastart <- theta2eta(lambda.init, .link , earg = .earg)
    }
  }), list( .link = link, .earg = earg,
            .truncate = truncate,
            .ilambda = ilambda, .imethod = imethod,
            .max.support = max.support, 
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- NCOL(eta) / c(M1 = 1)
    type.fitted <- if (length(extra$type.fitted))
                       extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "prob.t"))[1]

    lambda <- eta2theta(eta, .link , earg = .earg )
    sump <- matrix(0, NROW(eta), NCOL(eta))
    lhs.prob <- ppois( .max.support , lambda)  # 1 - P(upper tail)
    rmlife <- lambda * ppois(max.support - 1, lambda,
                             lower.tail = FALSE)
    Sume <- matrix(rmlife, NROW(eta), NCOL(eta))  # Corrected
    for (tval in .truncate ) {
      pmf <- dpois(tval, lambda)
      sump <- sump + pmf
      Sume <- Sume + pmf * tval
    }
    ans <- switch(type.fitted,
       "mean"      = (lambda - Sume) / (lhs.prob - sump),
       "lambda"    = lambda,
       "prob.t"    = sump + (1 - lhs.prob))
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .link = link,
           .truncate = truncate,
           .max.support = max.support, 
           .earg = earg ))),

  last = eval(substitute(expression({
    misc$link <- rep_len( .link , M)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg
  }), list( .link = link, .earg = earg,
            .truncate = truncate ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    lambda <- eta2theta(eta, .link , earg = .earg )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dgtpois(y, lambda = lambda, log = TRUE,
                                  truncate = .truncate ,
                                  max.support = .max.support )
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link,
           .truncate = truncate,
           .max.support = max.support, 
           .earg = earg ))),
  vfamily = c("gtpoisson"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    lambda <- eta2theta(eta, .link , earg = .earg )
    okay1 <- all(is.finite(lambda)) && all(0 < lambda)
    okay1
  }, list( .link = link,
           .truncate = truncate,
           .max.support = max.support, 
           .earg = earg ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    lambda <- eta2theta(eta, .link , earg = .earg )
    rgtpois(nsim * length(lambda), lambda, truncate = .truncate ,
              max.support = .max.support )
  }, list( .link = link,
           .truncate = truncate,
           .max.support = max.support, 
           .earg = earg ))),
 



  deriv = eval(substitute(expression({


    lambda <- eta2theta(eta, .link , earg = .earg )

    sump <-
    Sume <- matrix(0, NROW(eta), NCOL(eta))
    for (tval in .truncate ) {
      pmf <- dpois(tval, lambda)
      sump <- sump + pmf
      Sume <- Sume + pmf * tval
    }
    sump <- sump + ppois( .max.support , lambda, lower.tail = FALSE)
    Sume <- Sume + lambda * ppois( .max.support - 1,
                                   lambda, lower.tail = FALSE)


    pmf.deriv1 <- function(y, lambda)
      dpois(y-1, lambda) - dpois(y, lambda)
    pmf.deriv2 <- function(y, lambda)
      dpois(y-2, lambda) - 2 * dpois(y-1, lambda) +
      dpois(y  , lambda)

    sumderiv1 <-
    sumderiv2 <- matrix(0, NROW(eta), NCOL(eta))
    for (tval in .truncate ) {
      sumderiv1 <- sumderiv1 + pmf.deriv1(tval, lambda)
      sumderiv2 <- sumderiv2 + pmf.deriv2(tval, lambda)
    }
    sumderiv1 <- sumderiv1 + dpois( .max.support    , lambda)
    sumderiv2 <- sumderiv2 + dpois( .max.support - 1, lambda) -
                             dpois( .max.support    , lambda)




    dl.dlambda <- y / lambda - 1 + sumderiv1 / (1 - sump)
    dlambda.deta <- dtheta.deta(lambda, .link , earg = .earg )
    c(w) * dl.dlambda * dlambda.deta
  }), list( .link = link, .earg = earg,
            .max.support = max.support, 
            .truncate = truncate ))),
  weight = eval(substitute(expression({


    ned2l.dlambda2 <- (lambda - Sume) / ((1 - sump) * lambda^2) -
       sumderiv2 / (1 - sump) -
      (sumderiv1 / (1 - sump))^2


    wz <- ned2l.dlambda2 * dlambda.deta^2
    c(w) * wz
  }), list( .link = link,
            .max.support = max.support ))))
}  # gtpoisson







dgtbinom <-
  function(x, size, prob, truncate = 0,  # NULL,
           log = FALSE) {
  if (is.null(truncate))
    return(dbinom(x, size, prob, log = log))
    
  LLL <- max(length(x), length(size), length(prob))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(size  ) != LLL) size   <- rep_len(size,   LLL)
  if (length(prob  ) != LLL) prob   <- rep_len(prob,   LLL)

  if (is.list(truncate))
    truncate <- unlist(truncate)
    
  if (!is.Numeric(truncate, integer.valued = TRUE) ||
      any(truncate < 0))  # || max(truncate) > size
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")
    
  sump <- numeric(LLL)
  for (tval in truncate)
    sump <- sump + dbinom(tval, size, prob)
  vecTF <- x %in% truncate
    
  ans <- if (log) {
    ifelse(vecTF, -Inf,  # log(0),
           dbinom(x, size, prob, log = TRUE) - log1p(-sump))
  } else {
    ifelse(vecTF, 0, dbinom(x, size, prob)/ (1 - sump))
  }
    
  ans
}  # dgtbinom



pgtbinom <-
  function(q, size, prob, truncate = 0) {
  if (is.null(truncate))
    return(pbinom(q, size, prob))

  if (is.list(truncate))
    truncate <- unlist(truncate)

  if (!is.Numeric(truncate, integer.valued = TRUE) ||
      any(truncate < 0) ||
      max(size) < max(truncate))
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")

  LLL <- max(length(q), length(size), length(prob))
  if (length(q)      != LLL) q      <- rep_len(q,    LLL)
  if (length(size  ) != LLL) size   <- rep_len(size, LLL)
  if (length(prob  ) != LLL) prob   <- rep_len(prob, LLL)

  sump <- numeric(LLL)
  for (tval in truncate)
    sump <- sump + dbinom(tval, size, prob)
  denom <- 1 - sump

  numer <- pbinom(q, size, prob)
  for (tval in truncate) {
    if (any(vecTF <- tval <= q))
      numer[vecTF] <- numer[vecTF] -
                      dbinom(tval, size[vecTF], prob[vecTF])
  }
  ans <- numer / denom
  ans
}  # pgtbinom



qgtbinom <-
  function(p, size, prob, truncate = 0) {
  if (is.null(truncate))
    return(qbinom(p, size, prob))

  if (is.list(truncate))
    truncate <- unlist(truncate)
    
  if (!is.Numeric(truncate, integer.valued = TRUE) ||
      any(truncate < 0) ||
      max(size) < max(truncate))
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")


  LLL <- max(length(p), length(size), length(prob))
  if (length(p)      != LLL) p      <- rep_len(p,    LLL)
  if (length(size  ) != LLL) size   <- rep_len(size, LLL)
  if (length(prob  ) != LLL) prob   <- rep_len(prob, LLL)

  min.support <- ifelse(any(truncate == 0), min(truncate+1), 0)

  ans <- p + size + prob

  bad0 <- !is.finite(prob) | prob <= 0 | 1 <= prob |
          !is.finite(size) | size <= 0 | round(size) != size
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(min.support - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <- size + 0.5
  dont.iterate <- bad

  foo <- function(q, size, prob, truncate = 0, p) {
    use.truncate <- truncate[truncate <= min(size)]
    use.truncate <- if (length(use.truncate) == 0) {
      NULL  # Otherwise numeric(0)
    } else {
      use.truncate[!is.na(use.truncate)]
    }
    pgtbinom(q, size, prob, truncate = use.truncate) - p
  }

  lhs <- dont.iterate |
         p <= dgtbinom(min.support, size, prob, truncate)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                      size = size[!lhs], prob = prob[!lhs],
                      truncate = truncate,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgtbinom(faa, size[!lhs], prob[!lhs],
                        truncate = truncate) < p[!lhs] &
             p[!lhs] <= pgtbinom(faa+1, size[!lhs], prob[!lhs],
                                   truncate = truncate),
             faa+1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (length(truncate))
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgtbinom(min.support, size, prob,
                           truncate = truncate)
  ans[vecTF] <- min.support  # min.support[vecTF]

  ans[!bad0 & !is.na(p) & p == 0] <- min.support
  ans[!bad0 & !is.na(p) & p == 1] <- size[!bad0 & !is.na(p) & p == 1]
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgtbinom




rgtbinom <-
  function(n, size, prob, truncate = 0) {
    qgtbinom(runif(n), size = size, prob = prob,
               truncate = truncate)
}  # rgtbinom






 gtbinomial <-
  function(truncate = 0,  # NULL,  # 0  #,
           zero = NULL,
           link = "logitlink",
           type.fitted = c("mean", "prob", "prob.t"),
           multiple.responses = FALSE, parallel = FALSE) {



  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.logical(multiple.responses) ||
      length(multiple.responses) != 1)
    stop("bad input for argument 'multiple.responses'")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "prob", "prob.t"))[1]


  new("vglmff",
  blurb = c("Generally-truncated binomial distribution (GT-Binom)\n\n",
            "Links:    ",
            if (multiple.responses)
            c(namesof("prob1", link, earg = earg, tag = FALSE),
              ",...,",
              namesof("probM", link, earg = earg, tag = FALSE)) else
              namesof("prob", link, earg = earg, tag = FALSE),
            if (length( truncate ))
            c("\nTruncated at:    ",
            paste( truncate , collapse = ", ")) else NULL),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints = constraints)

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .parallel = parallel, .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = .multiple.responses ,
         parameters.names = c("prob"),
         truncate = .truncate ,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .truncate = truncate,
           .type.fitted = type.fitted,
           .multiple.responses = multiple.responses ))),

  initialize = eval(substitute(expression({
    if ( .multiple.responses ) {
      if (is.factor(y))
        stop("response cannot be a factor ",
             "if 'multiple.responses = TRUE'")
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

      no.successes <- y * w
      if (any(abs(no.successes - round(no.successes)) > 0.001))
        warning("Number of successes is not integer-valued")
      if (any(abs(w - round(w)) > 0.0001))
        warning("non-integer 'weights' (number of trials)")
      if (any(y < 0 | 1 < y))
        stop("y values must be 0 <= y <= 1")
 
      if (!length(mustart) && !length(etastart))
        mustart <- matrix(colSums(no.successes) / colSums(w),
                          n, ncol(y), byrow = TRUE)
    } else {
      if (NCOL(y) == 1) {
        if (is.factor(y))
          y <- as.matrix(y != levels(y)[1])
        y <- as.matrix(y)
        w <- as.matrix(w)
        if (any(y < 0 | 1 < y))
          stop("response values 'y' must be 0 <= y <= 1")
        if (NCOL(w) != 1)
          stop("argument 'weights' has too many columns")
        no.successes <- y * w
        if (any(abs(no.successes - round(no.successes)) > 0.001))
          warning("Number of successes is not integer-valued")
        if (!length(mustart) && !length(etastart))
          mustart <- (0.5 + no.successes) / (1 + w)
      } else if (NCOL(y) == 2) {
        if (NCOL(w) != 1)
          stop("argument 'weights' has too many columns")
        if (!all(w == 1))
          stop("enter positive weights using argument 'form2'")
        if (min(y) < 0)
          stop("Negative data not allowed!")
        if (any(abs(y - round(y)) > 0.001))
          warning("Count data is not integer-valued")
        w <- cbind(y[, 1] + y[, 2])  # -> nvec
        y <- cbind(y[, 1] / w)
        no.successes <- y * w
        if (!length(mustart) && !length(etastart))
          mustart <- (0.5 + no.successes) / (1 + w)
        } else {
          stop("response 'y' must be a ",
               "vector of 0s and 1s or a vector of \n",
               "proportions and 'weight' specified,\nor a factor",
               " (first level = fail, other levels = success)",
               ",\n",
               "or a 2-column matrix where col 1 is the no. of ",
               "successes and col 2 is the no. of failures")
        }

    }  # Not multiple.responses


    NOS <- ncoly <- NCOL(y)
    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly



    extra$pwts2 <-
    pwts2 <- if (is.null(Ym2)) 1 else matrix(Ym2, n, NOS)
    extra$w <- w  # Use possibly in @linkinv
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)


    # Check that no truncated values exist in the response:
    if (any(round(no.successes) %in% .truncate ))
      stop("some response values equal values of the 'truncate' ",
           "argument")




    if ( .multiple.responses ) {
      dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
      dn2 <- if (length(dn2) == M) {
        paste("E[", dn2, "]", sep = "")
      } else {
        param.names("prob", M)
      }
      predictors.names <-
        namesof(if (M > 1) dn2 else "prob",
                .link , earg = .earg, short = TRUE)

    } else {
      predictors.names <-
        namesof("prob", .link , earg = .earg , tag = FALSE)
    }


    if (!length(etastart)) {
      etastart <- cbind(theta2eta(mustart, .link , earg = .earg ))
    }
    mustart <- NULL
  }), list( .link = link, .earg = earg,
            .truncate = truncate,
            .type.fitted = type.fitted,
            .multiple.responses = multiple.responses ))),


  linkinv = eval(substitute(function(eta, extra = NULL) {
    binprob <- eta2theta(eta, .link , earg = .earg )
    type.fitted <- if (length(extra$type.fitted))
                     extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }
    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "prob.t"))[1]
    Size <- round(extra$w)
    NOS <- NCOL(eta)
    sump <- Sume <- matrix(0, NROW(eta), NCOL(eta))
    for (tval in .truncate ) {
      pmf <- dbinom(tval, Size, binprob)
      sump <- sump + pmf
      Sume <- Sume + pmf * tval / Size  # correct
    }
    ans <- switch(type.fitted,
       "mean"   = (binprob - Sume) / (1 - sump),
       "prob"   = binprob,
       "prob.t" = sump)  # Pr(Y=truncatedvalue) as it were
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  },
  list( .link = link, .earg = earg,
        .truncate = truncate,
        .type.fitted = type.fitted,
        .multiple.responses = multiple.responses ))),
  last = eval(substitute(expression({

    misc$link <- rep_len( .link , M)
    names(misc$link) <- if (M > 1) dn2 else "prob"

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$multiple.responses   <- .multiple.responses
    w <- as.numeric(w)
  }), list( .link = link, .earg = earg,
            .truncate = truncate,
            .multiple.responses = multiple.responses ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    ycounts <- round(y * w)
    Size <- round(w)
    binprob <- eta2theta(eta, .link , earg = .earg )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
        answer <- c(extra$pwts2) *
          dgtbinom(ycounts, Size, truncate = .truncate ,
                     prob = binprob, log = TRUE)
      ll.elts <- answer
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg,
           .truncate = truncate,
           .multiple.responses = multiple.responses ))),

  vfamily = c("gtbinomial"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    binprob <- eta2theta(eta, .link , earg = .earg )
    okay1 <- all(is.finite(binprob)) &&
             all(0 < binprob & binprob < 1)
    okay1
  }, list( .link = link,
           .truncate = truncate,
           .earg = earg ))),





  simslot = eval(substitute(
  function(object, nsim) {
    pwts2 <- if (length(pwts2 <- object@extra$pwts2) > 0)
               pwts2 else
             if (length(pwts2 <- depvar(object, type = "lm2")) > 0)
               pwts2 else 1
    if (any(pwts2 != 1))
      warning("ignoring prior weights")
    if ( .multiple.responses )
      stop("cannot run simulate() when 'multiple.responses = TRUE'")
    eta <- predict(object)
    binprob <- c(eta2theta(eta, .link , earg = .earg ))
    Size <- c(round(weights(object, type = "prior")))
    rgtbinom(nsim * length(eta), size = Size, prob = binprob,
               truncate = .truncate )
  }, list( .link = link, .earg = earg,
           .truncate = truncate,
           .multiple.responses = multiple.responses ))),



  deriv = eval(substitute(expression({
    pwts2 <- if (is.numeric(extra$pwts2)) extra$pwts2 else 1
    ycounts <- round(y * w)
    Size <- round(w)
    binprob <- eta2theta(eta, .link , earg = .earg )
    dmu.deta <- dtheta.deta(binprob, .link , earg = .earg )

    sump <- Sume <- matrix(0, NROW(eta), NCOL(eta))
    for (tval in .truncate ) {
      pmf <- dbinom(tval, Size, binprob)
      sump <- sump + pmf
      Sume <- Sume + pmf * tval / Size  # correct
    }

    pmf.deriv1 <- function(ycount, size, prob)
      dbinom(ycount, size, prob) * size *
      (ycount / size - prob) / (prob * (1 - prob))
    pmf.deriv2 <- function(ycount, size, prob)
      dbinom(ycount, size, prob) * size * (
      size * ((    ycount / size  /      prob) -
              (1 - ycount / size) / (1 - prob))^2 -
           ycount / size  /      prob^2 -
      (1 - ycount / size) / (1 - prob)^2)

    sumderiv1 <- sumderiv2 <- matrix(0, NROW(eta), NCOL(eta))
    for (tval in .truncate ) {
      sumderiv1 <- sumderiv1 + pmf.deriv1(tval, Size, binprob)
      sumderiv2 <- sumderiv2 + pmf.deriv2(tval, Size, binprob)
    }

    dl.dmu <- Size *      y  /      binprob -
              Size * (1 - y) / (1 - binprob) +
              sumderiv1 / (1 - sump)  # - (1-binprob)*temp3/temp1
    c(pwts2) * dl.dmu * dmu.deta
  }), list( .link = link, .earg = earg,
            .truncate = truncate,
            .multiple.responses = multiple.responses ))),
  weight = eval(substitute(expression({
    EY <- (binprob - Sume) / (1 - sump)  # The expectation of Y
    ned2l.dmu2 <- Size * EY / binprob^2 +
                  Size * (1 - EY) / (1 - binprob)^2 -
       sumderiv2 / (1 - sump) -
      (sumderiv1 / (1 - sump))^2

    wz <- ned2l.dmu2 * dmu.deta^2  # Size already inside
    c(pwts2) * wz
  }), list( .link = link, .earg = earg,
            .truncate = truncate,
            .multiple.responses = multiple.responses ))))
}  # gtbinomial






dgapois <-
  function(x, lambda, alter = 0,
           pobs.a = 0, byrow.arg = FALSE,
           log = FALSE) {
  if (is.null(alter))
    return(dpois(x, lambda, log = log))
    
  if (is.list(alter))
    alter <- unlist(alter)

  log.arg <- log
  rm(log)

  if (!is.Numeric(alter, integer.valued = TRUE) ||
      any(alter < 0))
    stop("bad input for argument 'alter'")
  if (!identical(alter, (unique(alter))))
    stop("values of argument 'alter' must be unique")

  if (any(pobs.a < 0) || any(1 < pobs.a))
    stop("bad input for argument 'pobs.a'")

  LLL <- max(length(x), length(lambda))  #, length(pobs.a)
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)

  lalter <- length(alter)
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  suma <- rowSums(pobs.a)
  if (any(1 < suma))
    stop("bad input for argument 'pobs.a'")

  if (log.arg) {
    logpmf <- log1p(-suma) +
              dgtpois(x, lambda, truncate = alter, log = TRUE)
  } else {
    pmf <- exp(log1p(-suma)) * dgtpois(x, lambda, truncate = alter)
  }
 for (jay in seq(lalter)) {
    aval <- alter[jay]
    if (any(vecTF <- aval == x)) {
      if (log.arg) {
        logpmf[vecTF] <- log(pobs.a[vecTF, jay])
      } else {
        pmf[vecTF] <- pobs.a[vecTF, jay]
      }
    }
  }  # jay
  if (log.arg) logpmf else pmf
}  # dgapois




pgapois <-
  function(q, lambda, alter = 0,
           pobs.a = 0, byrow.arg = FALSE) {
  if (is.null(alter))
    return(ppois(q, lambda))
    
  if (is.list(alter))
    truncate <- unlist(alter)
    
  if (!is.Numeric(alter, integer.valued = TRUE) ||
      any(alter < 0))
    stop("bad input for argument 'alter'")
  if (!identical(alter, (unique(alter))))
    stop("values of argument 'alter' must be unique")

  if (any(pobs.a < 0) || any(1 < pobs.a))
    stop("bad input for argument 'pobs.a'")

  LLL <- max(length(q), length(lambda))  # , length(pobs.a)
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)


  lalter <- length(alter)
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  suma <- rowSums(pobs.a)
  if (any(1 < suma))
    stop("bad input for argument 'pobs.a'")
  numer1 <- 1 - suma

  offset <- numeric(LLL)
  for (jay in seq(lalter)) {
    aval <- alter[jay]
    if (any(vecTF <- aval <= q)) {
      offset[vecTF] <- offset[vecTF] + pobs.a[vecTF, jay]
    }
  }
  ans <- numer1 * pgtpois(q, lambda, truncate = alter) + offset
  ans
}  # pgapois




qgapois <-
  function(p, lambda, alter = 0,  # NULL,
           pobs.a = 0, byrow.arg = FALSE) {
  if (is.null(alter))
    return(qpois(p, lambda))

  if (is.list(alter))
    alter <- unlist(alter)
    
  if (!is.Numeric(alter, integer.valued = TRUE) ||
      any(alter < 0))
    stop("bad input for argument 'alter'")
  if (!identical(alter, (unique(alter))))
    stop("values of argument 'alter' must be unique")

  if (any(pobs.a < 0) || any(1 < pobs.a))
    stop("bad input for argument 'pobs.a'")

  LLL <- max(length(p), length(lambda))
  if (length(p)      != LLL) p      <- rep_len(p,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)

  lalter <- length(alter)
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  min.support <- 0
  ans <- p + lambda

  bad0 <- !is.finite(lambda) | lambda <= 0
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(min.support - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <- 2 * lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgapois(hi, lambda, alter = alter,
                   pobs.a = pobs.a, byrow.arg = FALSE)
  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    lo[!done] <- hi[!done]
    hi[!done] <- 2 * hi[!done] + 10.5  # Bug fixed
    done[!done] <-
      (p[!done] <= pgapois(hi[!done], lambda[!done],
                             alter = alter,
                             pobs.a = pobs.a[!done, ],
                             byrow.arg = FALSE))
    iter <- iter + 1
  }

  foo <- function(q, lambda, alter = 0,
                  pobs.a = 0, byrow.arg = FALSE, p)
    pgapois(q, lambda = lambda, alter = alter,
              pobs.a = pobs.a,
              byrow.arg = FALSE) - p

  lhs <- dont.iterate |
         p <= dgapois(min.support, lambda, alter = alter,
                        pobs.a = pobs.a, byrow.arg = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                      lambda = lambda[!lhs],
                      alter = alter, pobs.a = pobs.a[!lhs, ],
                      byrow.arg = FALSE, p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgapois(faa, lambda[!lhs],  alter = alter,
                       pobs.a = pobs.a[!lhs, ],
                       byrow.arg = FALSE) < p[!lhs] &
             p[!lhs] <= pgapois(faa+1, lambda[!lhs],
                                  alter = alter,
                                  pobs.a = pobs.a[!lhs, ],
                                  byrow.arg = FALSE),
             faa+1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)

  vecTF <- !bad0 & !is.na(p) &
           p <= dgapois(min.support, lambda, alter = alter,
                          pobs.a = pobs.a, byrow.arg = FALSE)
  ans[vecTF] <- min.support

  ans[!bad0 & !is.na(p) & p == 0] <- min.support
  ans[!bad0 & !is.na(p) & p == 1] <- Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgapois



rgapois <-
  function(n, lambda, alter = 0,  # NULL
           pobs.a = 0, byrow.arg = FALSE) {
  qgapois(runif(n), lambda = lambda, alter = alter,
            pobs.a = pobs.a, byrow.arg = byrow.arg)
}  # rgapois










if (FALSE)
gapoisson.control <-
  function(summary.HDEtest = FALSE,  # Overwrites the summary() default.
           ...) {
  list(summary.HDEtest = summary.HDEtest)
}



 gapoisson.mlm <-
  function(alter = 0,
           zero = NULL,
           llambda = "loglink",
           type.fitted = c("mean", "lambda", "pobs.a",
               "onempobs.a", "Pobs.a"),
           imethod = 1,
           ilambda = NULL, ishrinkage = 0.95,
           probs.y = 0.35) {


  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  if (is.list(alter))
    alter <- unlist(alter)
  if (is.character(alter) && alter == "")
    alter <- NULL
  if (length(alter)) {
    if (!is.Numeric(alter, integer.valued = TRUE))
      stop("bad input for argument 'alter'")
    if (any(alter < 0))
      stop("values for argument 'alter' must be nonnegative")
    if (!identical(alter, (unique(alter))))
      stop("values of argument 'alter' must be unique")
  } else {
    stop("argument 'alter' is effectively empty; ",
         "use the family function poissonff() instead")
  }
 
  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "pobs.a",
                             "onempobs.a", "Pobs.a"))[1]
  temp7 <- if (length(alter)) paste0("pobs", alter) else NULL
  tmp3 <- c(if (length(alter))
            rep("multilogitlink",  length(alter)) else NULL,
            lambda = llambda )
  names(tmp3) <- c(temp7, "lambda")


  new("vglmff",
  blurb = c("Generally-altered Poisson distribution\n",
            "(GA-Pois-MLM)",
            "\n\n",
            "Links:    ", if (length(alter))
      paste("multilogitlink(cbind(",
            paste(temp7, collapse = ", "),
            ", ",
            "\n                           ",
            "1 - ", paste(temp7, collapse = " - "),
            "))", ", ", sep = ""),
            if (length(alter)) "\n          ",
            namesof("lambda", llambda, earg = elambda, tag = FALSE)),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = length( .alter ) + 1)
  }), list( .zero = zero,
            .alter = alter ))),

  infos = eval(substitute(function(...) {
    alter <- ( .alter )
    temp7 <- paste("pobs", alter, sep = "")
    list(M1 = length( .alter ) + 1,
         Q1 = 1,
         link = .tmp3 ,
         link1parameter = if (length( .alter ))
               FALSE else TRUE,  # multilogitlink is multiparameter
         mixture.links = TRUE,
         alter = .alter ,
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE,
         parameters.names = c(if (length( .alter )) temp7 else NULL,
                              "lambda"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted,
           .llambda = llambda, .elambda = elambda,
           .alter = alter,
           .tmp3 = tmp3
         ))),

  initialize = eval(substitute(expression({
    alter <- ( .alter )
    lalter <- length(alter)
    M1 <- lalter + 1
    NOS <- NCOL(y)
    M <- NOS * M1

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,  # Inf,
              ncol.y.max = 1,  # Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    extra$y0 <- y0 <- matrix(0, n, lalter)
    for (jay in seq(lalter))
      extra$y0[, jay] <- y0[, jay] <- as.numeric(y %in% alter[jay])
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, lalter)
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)
    if (any((css <- colSums(skip.these)) == 0))
      stop("some 'alter' argument values have no response values: ",
           paste(alter[css == 0], collapse = ", "))          

    
    fillerChar <- " "
    use.refLevel <- M1+1  # Assumes only one response
    allbut.refLevel <- (1:(M+1))[-use.refLevel]
    predictors.names <-
      paste("log(promega[,", allbut.refLevel,
            "]", fillerChar, "/", fillerChar, "promega[,",
            use.refLevel, "])", sep = "")

    temp7 <- paste0("pobs", alter)
    denom.char <- paste0("1-", paste0(temp7, collapse = "-"))
    mynames1 <- paste("log(", temp7, "/(",
                      denom.char, "))", sep = "")
    mynames2 <- param.names("lambda", ncoly, skip1 = TRUE)
    predictors.names <-
        c(        mynames1,
          namesof(mynames2, .llambda , earg = .elambda ,
                  tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                             imu = .ilambda ,  # x = x,
                             ishrinkage = .ishrinkage ,
                             pos.only = TRUE,
                             probs.y = .probs.y )

      phimat <- colMeans(skip.these)  # yettodo: weight this by 'w'
      phimat <- matrix(phimat, n, lalter, byrow = TRUE)
      etastart <-  multilogitlink(cbind(phimat, 1 - rowSums(phimat)))
      etastart <-
        cbind(etastart,
              theta2eta(lambda.init, .llambda , earg = .elambda ))
 #   etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .llambda = llambda,
            .elambda = elambda,
            .ilambda = ilambda,
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .imethod = imethod,
            .alter = alter,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted))
                  extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "lambda", "pobs.a",
                               "onempobs.a", "Pobs.a"))[1]

    alter <- ( .alter )
    M1 <- length(alter) + 1
    NOS <- NCOL(eta) / M1

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    ynames.Pobs.a <- c(as.character(alter), "(Others)")
    dimnames(phimat) <- list(rownames(eta), ynames.Pobs.a)
    lambda <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                        .llambda , earg = .elambda )
    pobs.a <- rowSums(phimat[, -NCOL(eta), drop = FALSE])
    sump <- matrix(0, NROW(eta), NOS)
    Sume <- matrix(0, NROW(eta), NOS)
    for (aval in alter) {
      pmf <- dpois(aval, lambda)
      sump <- sump + pmf
      Sume <- Sume + pmf * aval
    }
   
    ans <- switch(type.fitted,
      "mean"       = (1 - pobs.a) * (lambda - Sume) / (1 - sump),
      "lambda"     = lambda,
      "pobs.a"     =     pobs.a,  # Pr(Y is altered)
      "onempobs.a" = 1 - pobs.a,  # Pr(Y is not altered)
      "Pobs.a"     =     phimat)  # matrix
    label.cols.y(ans,
                 colnames.y = if (type.fitted  == "Pobs.a")
                              ynames.Pobs.a else extra$colnames.y,
                 NOS = NOS)
  }, list( .llambda = llambda,
           .elambda = elambda,
           .alter = alter ))),
  last = eval(substitute(expression({
    temp.names <- c(rep_len( "multilogitlink" , lalter),
                    rep_len( .llambda , NOS))
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)

    for (ii in seq(M1*NOS - 1)) {
        misc$earg[[ii]] <- list(M = M - 1,  # M * NOS,
                                refLevel = M)  # M * NOS
    }
    misc$earg[[M1*NOS]] <- .elambda  # Last one
  }), list( .llambda = llambda, .elambda = elambda,
            .alter = alter ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    lambda <- cbind(eta2theta(eta[,  NCOL(eta), drop = FALSE],
                              .llambda, earg = .elambda ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgapois(y, lambda = lambda, log = TRUE,
                  alter = .alter ,  # byrow.arg = FALSE,
                  pobs.a = pobs.a[, -NCOL(pobs.a), drop = FALSE])

      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llambda = llambda, .elambda = elambda,
           .alter = alter ))),
  vfamily = c("gapoisson.mlm"),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    lambda <- eta2theta(eta[,  NCOL(eta), drop = FALSE],
                        .llambda , earg = .elambda )
    okay1 <- all(is.finite(lambda)) && all(0 < lambda) &&
             all(is.finite(pobs.a)) && all(0 < pobs.a & pobs.a < 1)
    okay1
  }, list( .llambda = llambda, .elambda = elambda,
           .alter = alter ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    lambda <- eta2theta(eta[, ncol(eta)],
                        .llambda , earg = .elambda )
    rgapois(nsim * length(lambda), lambda = lambda,
              pobs.a = pobs.a, alter = .alter )
  }, list( .llambda = llambda, .elambda = elambda,
           .alter = alter ))),


  deriv = eval(substitute(expression({
    alter <- ( .alter )
    lalter <- length(alter)
    M1 <- lalter + 1
    NOS <- ncol(eta) / M1  # extra$NOS
    y0 <- extra$y0
    skip <- extra$skip.these
    is.altered <- rowSums(skip) > 0  # TRUE if (any) y %in% avec

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    lambda <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                        .llambda , earg = .elambda )


    sump <- Sume <- matrix(0, NROW(eta), NOS)
    for (aval in alter) {
      pmf <- dpois(aval, lambda)
      sump <- sump + pmf
      Sume <- Sume + pmf * aval
    }
    pmf.deriv1 <- function(y, lambda)
      dpois(y-1, lambda) - dpois(y, lambda)
    pmf.deriv2 <- function(y, lambda)
      dpois(y-2, lambda) - 2 * dpois(y-1, lambda) +
      dpois(y  , lambda)
    sumderiv1 <- sumderiv2 <- matrix(0, NROW(eta), NOS)
    for (aval in alter) {
      sumderiv1 <- sumderiv1 + pmf.deriv1(aval, lambda)
      sumderiv2 <- sumderiv2 + pmf.deriv2(aval, lambda)
    }

    pobs.a <- rowSums(phimat[, -NCOL(phimat), drop = FALSE])
    onempobs.a <- 1 - pobs.a

    dl.deta <- skip - phimat[, -M, drop = FALSE]

    dl.dlambda <- (1 - is.altered) *
      (y / lambda - 1 + sumderiv1 / (1 - sump))
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    ans <- cbind(c(w) * dl.deta,
                 c(w) * dl.dlambda * dlambda.deta)
    ans
  }), list( .llambda = llambda, .elambda = elambda,
            .alter = alter ))),


  weight = eval(substitute(expression({
    MM12 <- M1 * (M1 + 1) / 2

      ned2l.dlambda2 <- onempobs.a *
        ((lambda - Sume) / ((1 - sump) * lambda^2) -
         sumderiv2 / (1 - sump) -
        (sumderiv1 / (1 - sump))^2)

    wz4 <- matrix(0.0, n, MM12)  # A full matrix
    if (lalter > 0) {
      if (lalter == 1) {
        wz4[, 1] <- phimat[, 1] * (1 - phimat[, 1])
      } else {
        index <- iam(NA, NA, M - 1, both = TRUE, diag = TRUE)
        wz4 <- -phimat[, index$row] * phimat[, index$col]
        wz4[, 1:lalter] <- wz4[, 1:lalter] + phimat[, -M]
      }
    }
    wz <- wz.merge(c(w) * wz4,
                   c(w) * ned2l.dlambda2 * dlambda.deta^2,
                   M1 = M1 - 1, M2 = 1)  # rm.trailing.cols = F
    wz
  }), list( .alter = alter ))))
}  # gapoisson.mlm







dgipois <-
  function(x, lambda, inflate = 0,
           pstr.i = 0, byrow.arg = FALSE,
           log.arg = FALSE) {

  if (is.null(inflate))
    return(dpois(x, lambda, log = log.arg))
    
  if (is.list(inflate))
    inflate <- unlist(inflate)
    
  if (!is.Numeric(inflate, integer.valued = TRUE) ||
      any(inflate < 0))
    stop("bad input for argument 'inflate'")
  if (!identical(inflate, (unique(inflate))))
    stop("values of argument 'inflate' must be unique")

  if (any(pstr.i < 0) || any(1 < pstr.i))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(x), length(lambda))  # , length(pstr.i)
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)

  linflate <- length(inflate)
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)
  sumi <- .rowSums(pstr.i, LLL, linflate)
  if (any(1 < sumi))
    stop("bad input for argument 'pstr.i'")


  numer1 <- 1 - sumi

  pmf <- numer1 * dpois(x, lambda) 
  for (jay in seq(linflate)) {
    ival <- inflate[jay]
    if (any(vecTF <- ival == x)) {
      pmf[vecTF] <- pmf[vecTF] + pstr.i[vecTF, jay]
    }
  }

  if (log.arg) log(pmf) else pmf
}  # dgipois





pgipois <-
  function(q, lambda, inflate = 0,
           pstr.i = 0, byrow.arg = FALSE) {
  if (is.null(inflate))
    return(ppois(q, lambda))
    
  if (is.list(inflate))
    inflate <- unlist(inflate)
    
  if (!is.Numeric(inflate, integer.valued = TRUE) ||
      any(inflate < 0))
    stop("bad input for argument 'inflate'")
  if (!identical(inflate, (unique(inflate))))
    stop("values of argument 'inflate' must be unique")

  if (any(pstr.i < 0) || any(1 < pstr.i))
    stop("bad input for argument 'pstr.i'")
  
  LLL <- max(length(q), length(lambda))  # , length(pstr.i)
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)

  linflate <- length(inflate)
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)

  sumi <- rowSums(pstr.i)
  if (any(1 < sumi))
    stop("bad input for argument 'pstr.i'")
  numer1 <- 1 - sumi

  offset <- numeric(LLL)
  for (jay in seq(linflate)) {
    ival <- inflate[jay]
    if (any(vecTF <- ival <= q)) {
      offset[vecTF] <- offset[vecTF] + pstr.i[vecTF, jay]
    }
  }
  ans <- numer1 * ppois(q, lambda) + offset
  ans
}  # pgipois



qgipois <-
  function(p, lambda, inflate = 0,  # NULL,
           pstr.i = 0, byrow.arg = FALSE) {
  if (is.null(inflate))
    return(qpois(p, lambda))

  if (is.list(inflate))
    inflate <- unlist(inflate)
    
  if (!is.Numeric(inflate, integer.valued = TRUE) ||
      any(inflate < 0))
    stop("bad input for argument 'inflate'")
  if (!identical(inflate, (unique(inflate))))
    stop("values of argument 'inflate' must be unique")

  if (any(pstr.i < 0) || any(1 < pstr.i))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(p), length(lambda))
  if (length(p)      != LLL) p      <- rep_len(p,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)

  linflate <- length(inflate)
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)

  min.support <- 0
  ans <- p + lambda

  bad0 <- !is.finite(lambda) | lambda <= 0
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(min.support - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <- 2 * lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgipois(hi, lambda, inflate = inflate,
                   pstr.i = pstr.i, byrow.arg = FALSE)
  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    lo[!done] <- hi[!done]
    hi[!done] <- 2 * hi[!done] + 10.5  # Bug fixed
    done[!done] <-
      (p[!done] <= pgipois(hi[!done], lambda[!done],
                             inflate = inflate,
                             pstr.i = pstr.i[!done, ],
                             byrow.arg = FALSE))
    iter <- iter + 1
  }

  foo <- function(q, lambda, inflate = 0,
                  pstr.i = 0, byrow.arg = FALSE, p)
    pgipois(q, lambda = lambda, inflate = inflate,
              pstr.i = pstr.i,
              byrow.arg = FALSE) - p

  lhs <- dont.iterate |
         p <= dgipois(min.support, lambda, inflate = inflate,
              pstr.i = pstr.i         ,
              byrow.arg = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                      lambda = lambda[!lhs],
                      inflate = inflate,
                      pstr.i = pstr.i[!lhs, ],
                      byrow.arg = FALSE,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgipois(faa, lambda[!lhs],  inflate = inflate,
                       pstr.i = pstr.i[!lhs, ],
                       byrow.arg = FALSE) < p[!lhs] &
             p[!lhs] <= pgipois(faa+1, lambda[!lhs],
                                  inflate = inflate,
                                  pstr.i = pstr.i[!lhs, ],
                                  byrow.arg = FALSE),
             faa+1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)

  vecTF <- !bad0 & !is.na(p) &
           p <= dgipois(min.support, lambda,
                          inflate = inflate,
                          pstr.i = pstr.i,
                          byrow.arg = FALSE)
  ans[vecTF] <- min.support

  ans[!bad0 & !is.na(p) & p == 0] <- min.support
  ans[!bad0 & !is.na(p) & p == 1] <- Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgipois



rgipois <-
  function(n, lambda, inflate = 0,  # NULL
           pstr.i = 0, byrow.arg = FALSE) {
    qgipois(runif(n), lambda = lambda, inflate = inflate,
              pstr.i = pstr.i, byrow.arg = byrow.arg)
}  # rgipois








if (FALSE)
gipoisson.control <-
  function(summary.HDEtest = FALSE,  # Overwrites summary() default.
           ...) {
  list(summary.HDEtest = summary.HDEtest)
}



 gipoisson.mlm <-
  function(inflate = 0,
           zero = NULL,
           llambda = "loglink",
           type.fitted = c("mean", "lambda", "pstr.i",
                           "onempstr.i", "Pstr.i"),
           imethod = 1,
           mux.inflate = 0.5,
           ipstr0 = NULL, ilambda = NULL, ishrinkage = 0.95,
           probs.y = 0.35) {

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  if (is.list(inflate))
    inflate <- unlist(inflate)
  if (is.character(inflate) && inflate == "")
    inflate <- NULL
  if (length(inflate)) {
    if (!is.Numeric(inflate, integer.valued = TRUE))
      stop("bad input for argument 'inflate'")
    if (any(inflate < 0))
      stop("values for argument 'inflate' must be nonnegative")
    if (!identical(inflate, (unique(inflate))))
      stop("values of argument 'inflate' must be unique")
  } else {
    stop("argument 'inflate' is effectively empty; ",
         "use the family function poissonff() instead")
  }
 
  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "pstr.i",
                             "onempstr.i", "Pstr.i"))[1]
  temp7 <- if (length(inflate)) paste0("pstr", inflate) else NULL
  tmp3 <- c(if (length(inflate))
            rep("multilogitlink",  length(inflate)) else NULL,
            prob = llambda )
  names(tmp3) <- c(temp7, "lambda")


  new("vglmff",
  blurb = c("Generally-inflated Poisson distribution\n",
            "(GI-Pois-MLM)\n\n",
            "Links:    ", if (length(inflate))
      paste("multilogitlink(cbind(",
            paste(temp7, collapse = ", "),
            ", ",
            "\n                           ",
            "1 - ", paste(temp7, collapse = " - "),
            "))", ", ", sep = ""),
            if (length(inflate)) "\n          ",
            namesof("lambda", llambda, earg = elambda, tag = FALSE)),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = length( .inflate ) + 1)
  }), list( .zero = zero,
            .inflate = inflate ))),

  infos = eval(substitute(function(...) {
    inflate <- ( .inflate )
    temp7 <- paste("pstr", inflate, sep = "")
    list(M1 = length( .inflate ) + 1,
         Q1 = 1,
         link = .tmp3 ,
         link1parameter = if (length( .inflate ))
              FALSE else TRUE,  # multilogitlink is multiparameter
         mixture.links = TRUE,
         inflate = .inflate ,
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE,
         parameters.names = c(if (length( .inflate )) temp7 else NULL,
                              "lambda"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted,
           .llambda = llambda, .elambda = elambda,
           .inflate = inflate,
           .tmp3 = tmp3
         ))),

  initialize = eval(substitute(expression({
    inflate <- ( .inflate )
    linflate <- length(inflate)
    M1 <- linflate + 1
    NOS <- NCOL(y)
    M <- NOS * M1

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,  # Inf,
              ncol.y.max = 1,  # Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    extra$y0 <- y0 <- matrix(0, n, linflate)
    for (jay in seq(linflate))
      extra$y0[, jay] <- y0[, jay] <- as.numeric(y == inflate[jay])
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)
    if (any((css <- colSums(y0)) == 0))
      stop("some 'inflate' argument values have no response ",
           "values: ",  paste(inflate[css == 0], collapse = ", "))          

    
    fillerChar <- " "
    use.refLevel <- M1+1  # Assumes only one response
    allbut.refLevel <- (1:(M+1))[-use.refLevel]
    predictors.names <-
      paste("log(prphi[,", allbut.refLevel,
            "]", fillerChar, "/", fillerChar, "prphi[,",
            use.refLevel, "])", sep = "")

    temp7 <- paste0("pstr", inflate)



    if (FALSE) {
    mynames1 <- paste("multinomiallink1(", temp7, ")", sep = "")
    mynames2 <- param.names("lambda", ncoly, skip1 = TRUE)
    predictors.names <-
        c(        mynames1,
          namesof(mynames2, .llambda , earg = .elambda , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]
    }  # FALSE


    denom.char <- paste0("1-", paste0(temp7, collapse = "-"))
    mynames1 <- paste("log(", temp7, "/(",
                      denom.char, "))", sep = "")
    mynames2 <- param.names("lambda", ncoly, skip1 = TRUE)
    predictors.names <-
        c(mynames1,
          namesof(mynames2, .llambda , earg = .elambda ,
                  tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                             imu = .ilambda ,  # x = x,
                             ishrinkage = .ishrinkage ,
                             pos.only = TRUE,
                             probs.y = .probs.y )

      phimat <- colSums(c(w) * y0) / c(colSums(w))  # Wted by 'w'.
      phimat <- phimat * ( .mux.inflate )
      phimat <- matrix(phimat, n, linflate, byrow = TRUE)
      etastart <-  multilogitlink(cbind(phimat, 1 - rowSums(phimat)))
      etastart <-
        cbind(etastart,
              theta2eta(lambda.init, .llambda , earg = .elambda ))
    }
  }), list( .llambda = llambda, .elambda = elambda,
            .ipstr0 = ipstr0, .ilambda = ilambda,
            .mux.inflate = mux.inflate,
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .imethod = imethod, .inflate = inflate,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <- if (length(extra$type.fitted))
                       extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "lambda", "pstr.i",
                               "onempstr.i", "Pstr.i"))[1]

    inflate <- ( .inflate )
    M1 <- length(inflate) + 1
    NOS <- NCOL(eta) / M1

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    ynames.Pstr.i <- c(as.character(inflate), "(Others)")
    dimnames(phimat) <- list(rownames(eta), ynames.Pstr.i)
    lambda <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                        .llambda , earg = .elambda )
    pstr.i <- rowSums(phimat[, -NCOL(eta), drop = FALSE])
    Sume <- matrix(0, NROW(eta), NOS)  # Different defn of Sume
    for (jay in seq(length(inflate)))
      Sume <- Sume + phimat[, jay] * inflate[jay]

    ans <- switch(type.fitted,
      "mean"       = (1 - pstr.i) * lambda + Sume,
      "lambda"     = lambda,
      "pstr.i"     =     pstr.i,  # P(Y is structurally inflated)
      "onempstr.i" = 1 - pstr.i,  # P(Y is not structurally inflated)
      "Pstr.i"     =     phimat)  # matrix
    label.cols.y(ans,
                 colnames.y = if (type.fitted  == "Pstr.i")
                              ynames.Pstr.i else extra$colnames.y,
                 NOS = NOS)
  }, list( .llambda = llambda,
           .elambda = elambda,
           .inflate = inflate ))),

  last = eval(substitute(expression({
    temp.names <- c(rep_len( "multilogitlink" , linflate),
                    rep_len( .llambda , NOS))
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)

    for (ii in seq(M1*NOS - 1)) {
        misc$earg[[ii]] <- list(M = M - 1,  # M * NOS,
                                refLevel = M)  # M * NOS
    }
    misc$earg[[M1*NOS]] <- .elambda  # Last one
  }), list( .llambda = llambda, .elambda = elambda,
            .inflate = inflate ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pstr.i <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    lambda <- cbind(eta2theta(eta[, NCOL(eta), drop = FALSE],
                              .llambda, earg = .elambda ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgipois(y, lambda = lambda, log.arg = TRUE,
                  inflate = .inflate ,  # byrow.arg = FALSE,
                  pstr.i = pstr.i[, -NCOL(pstr.i), drop = FALSE])

      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llambda = llambda, .elambda = elambda,
           .inflate = inflate ))),
  vfamily = c("gipoisson.mlm"),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pstr.i <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    lambda <- eta2theta(eta[,  NCOL(eta), drop = FALSE],
                        .llambda , earg = .elambda )
    okay1 <- all(is.finite(lambda)) && all(0 < lambda) &&
             all(is.finite(pstr.i)) && all(0 < pstr.i & pstr.i < 1)
    okay1
  }, list( .llambda = llambda, .elambda = elambda,
           .inflate = inflate ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pstr.i <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    lambda <- eta2theta(eta[, ncol(eta)],
                        .llambda , earg = .elambda )
    rgipois(nsim * length(lambda), lambda = lambda,
              pstr.i = pstr.i, inflate = .inflate )
  }, list( .llambda = llambda, .elambda = elambda,
           .inflate = inflate ))),

  deriv = eval(substitute(expression({
    inflate <- ( .inflate )
    linflate <- length(inflate)
    M1 <- linflate + 1
    NOS <- ncol(eta) / M1  # extra$NOS
    y0 <- extra$y0  # Restore the value
    is.inflated <- rowSums(y0) > 0  # TRUE if (any) y %in% ivec

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    lambda <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                        .llambda , earg = .elambda )
    pstr.i <- rowSums(phimat[, -NCOL(phimat), drop = FALSE])
    onempstr.i <- 1 - pstr.i

    pmf.deriv1 <- function(y, lambda)
      dpois(y-1, lambda) - dpois(y, lambda)
    pmf.deriv2 <- function(y, lambda)
      dpois(y-2, lambda) - 2 * dpois(y-1, lambda) +
      dpois(y  , lambda)

    dl.deta <- -phimat[, -M, drop = FALSE]
    for (jay in seq(linflate)) {
      dl.deta[, jay] <- dl.deta[, jay] + y0[, jay] *
      phimat[, jay] / (phimat[, jay] +
      onempstr.i * dpois(inflate[jay], lambda))
    }  # for jay

    dl.dlambda <- (1 - is.inflated) * (y / lambda - 1)
    for (jay in seq(linflate)) {
      dl.dlambda <- dl.dlambda + y0[, jay] *
        onempstr.i * pmf.deriv1(inflate[jay], lambda) / (
        phimat[, jay] + onempstr.i * dpois(inflate[jay], lambda))
    }  # for jay
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    ans <- cbind(dl.deta,
                 dl.dlambda * dlambda.deta) * c(w)
    ans
  }), list( .llambda = llambda, .elambda = elambda,
            .inflate = inflate ))),


  weight = eval(substitute(expression({
    MM12 <- M1 * (M1 + 1) / 2
    wz <- matrix(0.0, n, MM12)  # A full matrix

  if (FALSE) { 
      ned2l.dlambda2.old <- onempstr.i / lambda
      for (jay in seq(linflate))
        ned2l.dlambda2.old <- ned2l.dlambda2.old -
          (1 - inflate[jay] / lambda)^2 *
          phimat[, jay] * onempstr.i * dpois(inflate[jay], lambda) / (
          phimat[, jay] + onempstr.i * dpois(inflate[jay], lambda))
  }


  if (TRUE) { 
      ned2l.dlambda2 <- 1 / lambda
      for (jay in seq(linflate))
        ned2l.dlambda2 <- ned2l.dlambda2 -
          dpois(inflate[jay], lambda) * inflate[jay] / lambda^2
      for (jay in seq(linflate))
        ned2l.dlambda2 <- ned2l.dlambda2 -
          pmf.deriv2(inflate[jay], lambda) +
          onempstr.i  * (pmf.deriv1(inflate[jay], lambda)^2) / (
          phimat[, jay] + onempstr.i * dpois(inflate[jay], lambda))
      ned2l.dlambda2 <- ned2l.dlambda2 * onempstr.i
  }

    wz[, iam(M1, M1, M = M1)] <-
      c(w) * ned2l.dlambda2 * dlambda.deta^2

    for (jay in seq(linflate)) {
      wz[, iam(jay, M1, M = M1)] <- c(w) * dlambda.deta *
        phimat[, jay] * onempstr.i *
        pmf.deriv1(inflate[jay], lambda) / (
        phimat[, jay] + onempstr.i * dpois(inflate[jay], lambda))
    }  # for jay


    if (linflate > 0) {
      index <- iam(NA, NA, M = M-1, both = TRUE, diag = TRUE)
      wz4 <- -phimat[, index$row, drop = FALSE] *
              phimat[, index$col, drop = FALSE]
      wz4[, 1:linflate] <- wz4[, 1:linflate] + phimat[, -M]
      for (jay in seq(linflate)) {
        wz4[, jay] <- wz4[, jay] -
          phimat[, jay] * onempstr.i * dpois(inflate[jay], lambda) / (
          phimat[, jay] + onempstr.i * dpois(inflate[jay], lambda))
      }  # for jay
    }  # linflate > 0
    wz4 <- as.matrix(wz4)  #  Needed when linflate == 1.
    if (linflate > 0)
      for (jay in seq(linflate))
        for (kay in jay:linflate)
           wz[, iam(jay, kay, M = M1  )] <- c(w) *
          wz4[, iam(jay, kay, M = M1-1)]
    wz
  }), list( .inflate = inflate ))))
}  # gipoisson.mlm







dgibinom <-
  function(x, size, prob, inflate = 0,
           pstr.i = 0, byrow.arg = FALSE,
           log.arg = FALSE) {
  if (is.null(inflate))
    return(dbinom(x, size, prob, log = log.arg))
    
  if (is.list(inflate))
    inflate <- unlist(inflate)
    
  if (!is.Numeric(inflate, integer.valued = TRUE) ||
      any(inflate < 0) ||
      max(size) < max(inflate))
    stop("bad input for argument 'inflate'")
  if (!identical(inflate, (unique(inflate))))
    stop("values of argument 'inflate' must be unique")

  if (any(pstr.i < 0) || any(1 < pstr.i))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(x), length(size), length(prob))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(size  ) != LLL) size   <- rep_len(size,   LLL)
  if (length(prob  ) != LLL) prob   <- rep_len(prob,   LLL)

  linflate <- length(inflate)
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)
  sumi <- rowSums(pstr.i)
  if (any(1 < sumi))
    stop("bad input for argument 'pstr.i'")
  onempstr.i <- 1 - sumi
  pmf <- onempstr.i * dbinom(x, size, prob) 
  for (jay in seq(linflate)) {
    if (any(vecTF <- inflate[jay] == x))
      pmf[vecTF] <- pmf[vecTF] + pstr.i[vecTF, jay]
  }

  if (log.arg) log(pmf) else pmf
}  # dgibinom



pgibinom <-
  function(q, size, prob, inflate = 0,
           pstr.i = 0, byrow.arg = FALSE) {
  if (is.null(inflate))
    return(pbinom(q, size, prob))
    
  if (is.list(inflate))
    inflate <- unlist(inflate)
    
  if (!is.Numeric(inflate, integer.valued = TRUE) ||
      any(inflate < 0))
    stop("bad input for argument 'inflate'")
  if (!identical(inflate, (unique(inflate))))
    stop("values of argument 'inflate' must be unique")

  if (any(pstr.i < 0) || any(1 < pstr.i))
    stop("bad input for argument 'pstr.i'")
  
  LLL <- max(length(q), length(size), length(prob))
  if (length(q)      != LLL) q      <- rep_len(q,    LLL)
  if (length(size)   != LLL) size   <- rep_len(size, LLL)
  if (length(prob)   != LLL) prob   <- rep_len(prob, LLL)

  linflate <- length(inflate)
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)

  sumi <- rowSums(pstr.i)
  if (any(1 < sumi))
    stop("bad input for argument 'pstr.i'")
  numer1 <- 1 - sumi

  offset <- numeric(LLL)
  for (jay in seq(linflate)) {
    ival <- inflate[jay]
    if (any(vecTF <- ival <= q)) {
      offset[vecTF] <- offset[vecTF] + pstr.i[vecTF, jay]
    }
  }
  ans <- numer1 * pbinom(q, size, prob) + offset
  ans
}  # pgibinom



qgibinom <-
  function(p, size, prob, inflate = 0,  # NULL,
           pstr.i = 0, byrow.arg = FALSE) {
  if (is.null(inflate))
    return(qbinom(p, size, prob))

  if (is.list(inflate))
    inflate <- unlist(inflate)
    
  if (!is.Numeric(inflate, integer.valued = TRUE) ||
      any(inflate < 0))
    stop("bad input for argument 'inflate'")
  if (!identical(inflate, (unique(inflate))))
    stop("values of argument 'inflate' must be unique")

  if (any(pstr.i < 0) || any(1 < pstr.i))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(p), length(size), length(prob))
  if (length(p)      != LLL) p      <- rep_len(p,    LLL)
  if (length(size)   != LLL) size   <- rep_len(size, LLL)
  if (length(prob)   != LLL) prob   <- rep_len(prob, LLL)

  linflate <- length(inflate)
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)

  min.support <- 0
  ans <- p + size + prob

  bad0 <- !is.finite(prob) | prob <= 0 | 1 <= prob |
          !is.finite(size) | size <= 0 | round(size) != size
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(min.support - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <- size + 0.5
  dont.iterate <- bad

  foo <- function(q, size, prob, inflate = 0,
                  pstr.i = 0, byrow.arg = FALSE, p)
    pgibinom(q, size, prob = prob, inflate = inflate,
               pstr.i = pstr.i, byrow.arg = FALSE) - p

  lhs <- dont.iterate |
         p <= dgibinom(min.support, size, prob, inflate = inflate,
              pstr.i = pstr.i, byrow.arg = FALSE)
  
  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                      size = size[!lhs], prob = prob[!lhs],
                      inflate = inflate, pstr.i = pstr.i[!lhs, ],
                      byrow.arg = FALSE, p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgibinom(faa, size = size[!lhs], prob[!lhs],
                        inflate = inflate, pstr.i = pstr.i[!lhs, ],
                        byrow.arg = FALSE) < p[!lhs] &
             p[!lhs] <= pgibinom(faa+1, size = size[!lhs], prob[!lhs],
                                   inflate = inflate,
                                   pstr.i = pstr.i[!lhs, ],
                                   byrow.arg = FALSE),
             faa+1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)

  vecTF <- !bad0 & !is.na(p) &
           p <= dgibinom(min.support, size, prob,
                           inflate = inflate,
                           pstr.i = pstr.i, byrow.arg = FALSE)
  ans[vecTF] <- min.support

  ans[!bad0 & !is.na(p) & p == 0] <- min.support
  ans[!bad0 & !is.na(p) & p == 1] <- size[!bad0 & !is.na(p) & p == 1]
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgibinom





rgibinom <-
  function(n, size, prob, inflate = 0,  # NULL
           pstr.i = 0, byrow.arg = FALSE) {
    qgibinom(runif(n), size = size, prob, inflate = inflate,
              pstr.i = pstr.i, byrow.arg = byrow.arg)
}  # rgibinom





 if (FALSE)
 gibinomial <-
  function(inflate = 0,
           zero = NULL,
           lprob = "logitlink",
           type.fitted = c("mean", "prob", "pstr.i",
                           "onempstr.i", "Pstr.i"),
           imethod = 1,
           mux.inflate = 0.5,
           ipstr0 = NULL, iprob = NULL, ishrinkage = 0.95,
           probs.y = 0.35) {
  multiple.responses = FALSE  # Added 20190417, yettodo later

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  if (is.list(inflate))
    inflate <- unlist(inflate)
  if (is.character(inflate) && inflate == "")
    inflate <- NULL
  if (length(inflate)) {
    if (!is.Numeric(inflate, integer.valued = TRUE))
      stop("bad input for argument 'inflate'")
    if (any(inflate < 0))
      stop("values for argument 'inflate' must be nonnegative")
    if (!identical(inflate, (unique(inflate))))
      stop("values of argument 'inflate' must be unique")
  } else {
    stop("argument 'inflate' is effectively empty; ",
         "use the family function binomialff() instead")
  }
 
  type.fitted <- match.arg(type.fitted,
                           c("mean", "prob", "pstr.i",
                             "onempstr.i", "Pstr.i"))[1]
  temp7 <- if (length(inflate)) paste0("pstr", inflate) else NULL
  tmp3 <- c(if (length(inflate))
            rep("multilogitlink",  length(inflate)) else NULL,
            prob = lprob )
  names(tmp3) <- c(temp7, "prob")


  new("vglmff",
  blurb = c("Generally-inflated binomial distribution\n",
            "(GI-Binom-MLM)\n\n",
            "Links:    ", if (length(inflate))
      paste("multilogitlink(cbind(",
            paste(temp7, collapse = ", "),
            ", ",
            "\n                           ",
            "1 - ", paste(temp7, collapse = " - "),
            "))", ", ", sep = ""),
            if (length(inflate)) "\n          ",
            namesof("prob", lprob, earg = eprob, tag = FALSE)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = length( .inflate ) + 1)
  }), list( .zero = zero,
            .inflate = inflate ))),

  infos = eval(substitute(function(...) {
    inflate <- ( .inflate )
    temp7 <- paste("pstr", inflate, sep = "")
    list(M1 = length( .inflate ) + 1,
         Q1 = 1,
         link = .tmp3 ,
         link1parameter = if (length( .inflate ))
              FALSE else TRUE,  # multilogitlink is multiparameter
         mixture.links = TRUE,
         inflate = .inflate ,
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE, for later possibly
         parameters.names = c(if (length( .inflate )) temp7 else NULL,
                              "prob"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted,
           .lprob = lprob, .eprob = eprob,
           .inflate = inflate,
           .tmp3 = tmp3
         ))),

  initialize = eval(substitute(expression({


    if ( .multiple.responses ) {
      if (is.factor(y))
        stop("response cannot be a factor ",
             "if 'multiple.responses = TRUE'")
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

      no.successes <- w * y
      if (any(abs(w - round(w)) > 0.0001))
        warning("non-integer 'weights' (number of trials)")
      if (any(abs(no.successes - round(no.successes)) > 0.0001))
        warning("non-integer number of successes")
      if (any(y < 0 | 1 < y))
        stop("y values must be 0 <= y <= 1")
 
      if (!length(mustart) && !length(etastart))
        mustart <- matrix(colSums(no.successes) / colSums(w),
                          n, ncol(y), byrow = TRUE)
    } else {
      if (NCOL(y) == 1) {
        if (is.factor(y))
          y <- as.matrix(y != levels(y)[1])
        y <- as.matrix(y)
        w <- as.matrix(w)
        if (any(y < 0 | 1 < y))
          stop("response values 'y' must be 0 <= y <= 1")
        if (NCOL(w) != 1)
          stop("argument 'weights' has too many columns")
        no.successes <- y * w
        if (any(abs(no.successes - round(no.successes)) > 0.001))
          warning("Number of successes is not integer-valued")
        if (!length(mustart) && !length(etastart))
          mustart <- (0.5 + no.successes) / (1 + w)
      } else if (NCOL(y) == 2) {
        if (NCOL(w) != 1)
          stop("argument 'weights' has too many columns")
        if (!all(w == 1))
          stop("enter positive weights using argument 'form2'")
        if (min(y) < 0)
          stop("Negative data not allowed!")
        if (any(abs(y - round(y)) > 0.001))
          warning("Count data is not integer-valued")
        w <- cbind(y[, 1] + y[, 2])  # -> nvec
        y <- cbind(y[, 1] / w)
        no.successes <- y * w
        if (!length(mustart) && !length(etastart))
          mustart <- (0.5 + no.successes) / (1 + w)
        } else {
          stop("response 'y' must be a ",
               "vector of 0s and 1s or a vector of \n",
               "proportions and 'weight' specified,\nor a factor",
               " (first level = fail, other levels = success)",
               ",\n",
               "or a 2-column matrix where col 1 is the no. of ",
               "successes and col 2 is the no. of failures")
        }
    }  # Not multiple.responses





    inflate <- ( .inflate )
    linflate <- length(inflate)
    M1 <- linflate + 1
    NOS <- ncoly <- NCOL(y)
    M <- NOS * M1
    extra$NOS <- extra$ncoly <- ncoly  # Number of species
    extra$M1 <- M1


    extra$pwts2 <-
    pwts2 <- if (is.null(Ym2)) 1 else matrix(Ym2, n, NOS)
    extra$w <- w  # Use possibly in @linkinv
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)



    if (!all(inflate %in% c(unique(round(no.successes)))))
      stop("some values of the the 'inflate' argument ",
           "not found in the response.")


    extra$y0 <- y0 <- matrix(0, n, linflate)
    for (jay in seq(linflate))
      extra$y0[, jay] <- y0[, jay] <-
      as.numeric(round(no.successes) == inflate[jay])
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <-
      matrix(as.logical(y0), n, linflate)
    if (any((css <- colSums(y0)) == 0))
      stop("some 'inflate' argument values have no response ",
           "values: ", paste(inflate[css == 0], collapse = ", "))



    if ( .multiple.responses ) {
      warning("this code chunk is incomplete")
      dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
      dn2 <- if (length(dn2) == M) {
        paste("E[", dn2, "]", sep = "")
      } else {
        param.names("prob", M)
      }
      predictors.names <-
        namesof(if (M > 1) dn2 else "prob",
                .link , earg = .earg, short = TRUE)
    } else {
      temp7 <- paste0("pstr", inflate)
      denom.char <- paste0("1-", paste0(temp7, collapse = "-"))
      mynames1 <- paste("log(", temp7, "/(",
                        denom.char, "))", sep = "")
      mynames2 <- param.names("prob", ncoly, skip1 = TRUE)
      predictors.names <-
          c(        mynames1,
            namesof(mynames2, .lprob , earg = .eprob , tag = FALSE))[
            interleave.VGAM(M1*NOS, M1 = M1)]
  }


    if (!length(etastart)) {
      prob.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                           imu = .iprob , pos.only = TRUE)

      phimat <- colMeans(skip.these)  # yettodo: weight this by 'w'
      phimat <- matrix(phimat, n, linflate, byrow = TRUE)
      etastart <-  multilogitlink(cbind(phimat, 1 - rowSums(phimat)))
      etastart <-
        cbind(etastart,
              theta2eta(prob.init, .lprob , earg = .eprob ))
      mustart <- NULL  # 20181203
    }
  }), list( .lprob = lprob, .eprob = eprob,
            .ipstr0 = ipstr0, .iprob = iprob,
            .mux.inflate = mux.inflate,
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .imethod = imethod, .inflate = inflate,
            .multiple.responses = multiple.responses,
            .type.fitted = type.fitted ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Size <- round(extra$w)
    type.fitted <- if (length(extra$type.fitted))
                     extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "pstr.i",
                       "onempstr.i", "Pstr.i"))[1]

    inflate <- ( .inflate )
    M1 <- length(inflate) + 1
    NOS <- NCOL(eta) / M1

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes 1 response
                inverse = TRUE)
    ynames.Pstr.i <- c(as.character(inflate), "(Others)")
    dimnames(phimat) <- list(rownames(eta), ynames.Pstr.i)
    prob <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                      .lprob , earg = .eprob )
    colnames(prob) <- NULL  # Was, e.g., "logitlink(prob)".
    pstr.i <- rowSums(phimat[, -NCOL(eta), drop = FALSE])
    Sume <- matrix(0, NROW(eta), NOS)  # Different defn of Sume
    for (jay in seq(length(inflate)))
      Sume <- Sume + phimat[, jay] * inflate[jay] / Size

    ans <- switch(type.fitted,
      "mean"       = (1 - pstr.i) * prob + Sume,
      "prob"       = prob,
      "pstr.i"     =     pstr.i,  # P(Y is structurally inflated)
      "onempstr.i" = 1 - pstr.i,  # P(Y is not structurally inflated)
      "Pstr.i"     =     phimat)  # matrix
    label.cols.y(ans,
                 colnames.y = if (type.fitted  == "Pstr.i")
                              ynames.Pstr.i else extra$colnames.y,
                 NOS = NOS)
  }, list( .lprob = lprob, .eprob = eprob,
           .inflate = inflate ))),

  last = eval(substitute(expression({
    temp.names <- c(rep_len("multilogitlink" , linflate),
                    rep_len( .lprob , NOS))
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)

    for (ii in seq(M1*NOS - 1)) {
        misc$earg[[ii]] <- list(M = M - 1,  # M * NOS,
                                refLevel = M)  # M * NOS
    }
    misc$earg[[M1*NOS]] <- .eprob  # Last one
  }), list( .lprob = lprob, .eprob = eprob,
            .inflate = inflate ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    ycounts <- round(y * w)
    Size <- round(w)

    pstr.i <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    prob <- cbind(eta2theta(eta[, NCOL(eta), drop = FALSE],
                            .lprob, earg = .eprob ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(extra$pwts2) *
        dgibinom(ycounts, Size, prob = prob, log = TRUE,
                   inflate = .inflate , byrow.arg = FALSE,
                   pstr.i = pstr.i[, -NCOL(pstr.i), drop = FALSE])

      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .eprob = eprob,
           .inflate = inflate ))),
  vfamily = c("gibinomial"),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    Size <- round(extra$w)
    okay2 <- max( .inflate ) <= max(Size)  # yettodo: do this once
    pstr.i <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    prob <- eta2theta(eta[,  NCOL(eta), drop = FALSE],
                      .lprob , earg = .eprob )
    okay1 <- all(is.finite(prob)) &&
             all(0 < prob   & prob   < 1) &&
             all(is.finite(pstr.i)) &&
             all(0 < pstr.i & pstr.i < 1)
    okay1 && okay2
  }, list( .lprob = lprob, .eprob = eprob,
           .inflate = inflate ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts2 <- if (length(pwts2 <- object@extra$pwts2) > 0)
               pwts2 else
             if (length(pwts2 <- depvar(object, type = "lm2")) > 0)
               pwts2 else 1
    if (any(pwts2 != 1))
      warning("ignoring prior weights")
    if ( .multiple.responses )
      stop("cannot run simulate() when 'multiple.responses = TRUE'")
    Size <- round(weights(object, type = "prior"))

    eta <- predict(object)
    pstr.i <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    prob <- eta2theta(eta[, ncol(eta)], .lprob , earg = .eprob )
    rgibinom(nsim * length(prob), Size, prob = prob,
    pstr.i = kronecker(matrix(1, nsim, 1), pstr.i[, -ncol(pstr.i)]),
               inflate = .inflate )
  }, list( .lprob = lprob, .eprob = eprob,
           .multiple.responses = multiple.responses,
           .inflate = inflate ))),



  deriv = eval(substitute(expression({
    pwts2 <- if (is.numeric(extra$pwts2)) extra$pwts2 else 1
    Size <- round(w)

    inflate <- ( .inflate )
    linflate <- length(inflate)
    M1 <- linflate + 1
    NOS <- ncol(eta) / M1  # extra$NOS
    y0 <- extra$y0  # Restore the value
    skip <- extra$skip.these
    is.inflated <- rowSums(y0) > 0  # TRUE if (any) y %in% ivec

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    prob <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                      .lprob , earg = .eprob )
    pstr.i <- rowSums(phimat[, -NCOL(phimat), drop = FALSE])
    onempstr.i <- 1 - pstr.i

    pmf.deriv1 <- function(ycount, size, prob)
      dbinom(ycount, size, prob) * size *
      (ycount / size - prob) / (prob * (1 - prob))
    pmf.deriv2 <- function(ycount, size, prob)
      dbinom(ycount, size, prob) * size * (
      size * ((    ycount / size  /      prob) -
              (1 - ycount / size) / (1 - prob))^2 -
           ycount / size  /      prob^2 -
      (1 - ycount / size) / (1 - prob)^2)

    sumderiv1 <- sumderiv2 <- matrix(0, NROW(eta), NCOL(eta))
    for (ival in inflate) {
      sumderiv1 <- sumderiv1 + c(pmf.deriv1(ival, Size, prob))
      sumderiv2 <- sumderiv2 + c(pmf.deriv2(ival, Size, prob))
    }

    Delta <- matrix(0, n, linflate)
    for (jay in seq(linflate))
      Delta[, jay] <- phimat[, jay] + onempstr.i *
        dbinom(inflate[jay], Size, prob)
    
    dl.dphi <- matrix(-(1 - is.inflated) / onempstr.i, n, linflate)
    for (jay in seq(linflate)) {
      dl.dphi[, jay] <- dl.dphi[, jay] + y0[, jay] *
        (y0[, jay] - dbinom(inflate[jay], Size, prob)) / Delta[, jay]
    }  # for jay
    dl.deta <- matrix(0, n, linflate)
    for (jay in seq(linflate)) {
      for (ess in seq(linflate)) {
        dl.deta[, jay] <- dl.deta[, jay] + phimat[, ess] *
          ((ess == jay) - phimat[, jay]) * dl.dphi[, ess]
    }  # for ess
  }  # jay

    dl.dprob <- (1 - is.inflated) *
      Size * (y / prob - 1) / (1 - prob)
    for (jay in seq(linflate)) {
      dl.dprob <- dl.dprob + y0[, jay] * onempstr.i * 
        pmf.deriv1(inflate[jay], Size, prob = prob) / Delta[, jay]
    }  # for jay
    dprob.deta <- dtheta.deta(prob, .lprob , earg = .eprob )
    ans <- cbind(dl.deta,
                 dl.dprob * dprob.deta) * c(pwts2)
    ans
  }), list( .lprob = lprob, .eprob = eprob,
            .inflate = inflate ))),




  weight = eval(substitute(expression({



    if (linflate > 0) {
      index <- iam(NA, NA, M = M, both = TRUE)
      dphi.deta <- -phimat[, index$row, drop = FALSE] *
                    phimat[, index$col, drop = FALSE]
      dphi.deta[, 1:M] <- dphi.deta[, 1:M] + phimat[, 1:M]
    for (jay in seq(M))
      dphi.deta[, iam(jay, M, M)] <- 0
    }  # linflate > 0











    sump <- Sume <- 0  # For one response
    for (tval in inflate) {
      pmf <- dbinom(tval, Size, prob)
      sump <- sump + pmf
      Sume <- Sume + pmf * tval / Size  # correct
    }
    EY <- (prob - Sume) / (1 - sump)
    tmp7 <- Size * (EY / prob^2 + (1 - EY) / (1 - prob)^2)
    ned2l.dprob2 <- (1 - sump) * tmp7
    for (jay in seq(linflate))
      ned2l.dprob2 <- ned2l.dprob2 -
        pmf.deriv2(inflate[jay], Size, prob) + onempstr.i *
       (pmf.deriv1(inflate[jay], Size, prob))^2 / Delta[, jay]
    ned2l.dprob2 <- ned2l.dprob2 * onempstr.i
    wz <- matrix(0, n, dimm(M))
    wz[, iam(M1, M1, M1)] <-
      c(pwts2) * ned2l.dprob2 * dprob.deta^2  # As usual









    for (jay in seq(linflate)) {
      wz[, iam(jay, M1, M1)] <- 0  # Already initialized
      for (ess in seq(linflate)) {
        wz[, iam(jay, M1, M1)] <-
        wz[, iam(jay, M1, M1)] + c(pwts2) *
          dprob.deta *
          pmf.deriv1(inflate[ess], Size, prob) *
          (phimat[, ess] + onempstr.i * (ess == jay)) / Delta[, ess]
      }  # for ess
    }  # for jay





wz.ver1 <- matrix(0, n, linflate)
for (jay in seq(linflate))
  for (ess in seq(linflate))
  wz.ver1[, jay] <-
  wz.ver1[, jay] + wz[, iam(ess, M1, M1)] *
            dphi.deta[, iam(ess, jay, M1)]


 if (FALSE) {
}








      






    if (linflate >= 1 && TRUE) {
      wz7 <- matrix((1 - sump) / onempstr.i, n, dimm(M))
      for (bbb in 1:linflate) {
        for (vee in bbb:linflate) {
          for (ess in seq(linflate)) {
             wz7[, iam(bbb, vee, M = M)] <-
             wz7[, iam(bbb, vee, M = M)] +
       ((ess == bbb) - dbinom(inflate[ess], Size, prob)) *
       ((ess == vee) - dbinom(inflate[ess], Size, prob)) / Delta[, ess]
          }  # for ess
        }  # for vvv
      }  # for bbb
      wz7 <- wz7 * c(pwts2)
    }  # linflate > 1

    for (jay in seq(M1))
      wz7[, iam(jay, M1, M1)] <- 0  # Entire last column

    wz7 <- muxtXAX(wz7, dphi.deta, M)

    

    for (jay in seq(linflate))
      wz7[, iam(jay, M, M)] <- wz.ver1[, jay]  # See below too
    wz7[, iam(M, M, M = M)] <-
      c(pwts2) * ned2l.dprob2 * dprob.deta^2  # As usual
  
 




    wz7
  }), list( .inflate = inflate ))))
}  # gibinomial










dgabinom <-
  function(x, size, prob, alter = 0,
           pobs.a = 0, byrow.arg = FALSE,
           log = FALSE) {
  if (is.null(alter))
    return(dbinom(x, size, prob, log = log))
    
  if (is.list(alter))
    alter <- unlist(alter)

  log.arg <- log
  rm(log)

  if (!is.Numeric(alter, integer.valued = TRUE) ||
      any(alter < 0) ||
      max(size) < max(alter))
    stop("bad input for argument 'alter'")
  if (!identical(alter, (unique(alter))))
    stop("values of argument 'alter' must be unique")

  if (any(pobs.a < 0) || any(1 < pobs.a))
    stop("bad input for argument 'pobs.a'")

  LLL <- max(length(x), length(size), length(prob))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(size  ) != LLL) size   <- rep_len(size,   LLL)
  if (length(prob  ) != LLL) prob   <- rep_len(prob,   LLL)

  lalter <- length(alter)
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  suma <- rowSums(pobs.a)
  if (any(1 < suma))
    stop("bad input for argument 'pobs.a'")

  if (log.arg) {
    logpmf <- log1p(-suma) +
              dgtbinom(x, size, prob, truncate = alter, log = TRUE)
  } else {
    pmf <- exp(log1p(-suma)) *
           dgtbinom(x, size, prob, truncate = alter)
  }
 for (jay in seq(lalter)) {
    aval <- alter[jay]
    if (any(vecTF <- aval == x)) {
      if (log.arg) {
        logpmf[vecTF] <- log(pobs.a[vecTF, jay])
      } else {
        pmf[vecTF] <- pobs.a[vecTF, jay]
      }
    }
  }  # jay
  if (log.arg) logpmf else pmf
}  # dgabinom







pgabinom <-
  function(q, size, prob, alter = 0,
           pobs.a = 0, byrow.arg = FALSE) {
  if (is.null(alter))
    return(pbinom(q, size, prob))
    
  if (is.list(alter))
    truncate <- unlist(alter)
    
  if (!is.Numeric(alter, integer.valued = TRUE) ||
      any(alter < 0) ||
      max(size) < max(alter))
    stop("bad input for argument 'alter'")
  if (!identical(alter, (unique(alter))))
    stop("values of argument 'alter' must be unique")

  if (any(pobs.a < 0) || any(1 < pobs.a))
    stop("bad input for argument 'pobs.a'")

  LLL <- max(length(q), length(size), length(prob))
  if (length(q)      != LLL) q      <- rep_len(q,     LLL)
  if (length(size  ) != LLL) size   <- rep_len(size,  LLL)
  if (length(prob  ) != LLL) prob   <- rep_len(prob,  LLL)

  lalter <- length(alter)
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  suma <- rowSums(pobs.a)
  if (any(1 < suma))
    stop("bad input for argument 'pobs.a'")
  numer1 <- 1 - suma

  offset <- numeric(LLL)
  for (jay in seq(lalter)) {
    aval <- alter[jay]
    if (any(vecTF <- aval <= q)) {
      offset[vecTF] <- offset[vecTF] + pobs.a[vecTF, jay]
    }
  }
  ans <- numer1 * pgtbinom(q, size, prob, truncate = alter) + offset
  ans
}  # pgabinom



qgabinom <-
  function(p, size, prob, alter = 0,  # NULL,
           pobs.a = 0, byrow.arg = FALSE) {
  if (is.null(alter))
    return(qbinom(p, size, prob))

  if (is.list(alter))
    alter <- unlist(alter)

  if (!is.Numeric(alter, integer.valued = TRUE) ||
      any(alter < 0) ||
      max(size) < max(alter))
    stop("bad input for argument 'alter'")
  if (!identical(alter, (unique(alter))))
    stop("values of argument 'alter' must be unique")

  if (any(pobs.a < 0) || any(1 < pobs.a))
    stop("bad input for argument 'pobs.a'")

  LLL <- max(length(p), length(size), length(prob))
  if (length(p)      != LLL) p      <- rep_len(p,     LLL)
  if (length(size  ) != LLL) size   <- rep_len(size,  LLL)
  if (length(prob  ) != LLL) prob   <- rep_len(prob,  LLL)

  lalter <- length(alter)
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  min.support <- 0
  ans <- p + size + prob

  bad0 <- !is.finite(size) | size <  0 | round(size) != size |
          !is.finite(prob) | prob <= 0 | 1 <= prob
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  lo <- rep_len(min.support - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <- size + 10.5
  dont.iterate <- bad

  foo <- function(q, size, prob, alter = 0,
                  pobs.a = 0, byrow.arg = FALSE, p) {
    use.alter <- alter[alter <= min(size)]
    use.alter <- if (length(use.alter) == 0) {
      NULL  # Otherwise numeric(0)
    } else {
      use.alter[!is.na(use.alter)]
    }
    pgabinom(q, size = size, prob = prob, alter = alter,
               pobs.a = pobs.a, byrow.arg = FALSE) - p
  }

  lhs <- dont.iterate |
         p <= dgabinom(min.support, size = size, prob = prob,
                         alter = alter,
                         pobs.a = pobs.a, byrow.arg = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                      size = size[!lhs], prob = prob[!lhs],
                      alter = alter, pobs.a = pobs.a[!lhs, ],
                      byrow.arg = FALSE, p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgabinom(faa, size = size[!lhs], prob[!lhs],
                        alter = alter,
                        pobs.a = pobs.a[!lhs, ],
                        byrow.arg = FALSE) < p[!lhs] &
             p[!lhs] <= pgabinom(faa+1, size = size[!lhs],
                                   prob[!lhs], alter = alter,
                                   pobs.a = pobs.a[!lhs, ],
                                   byrow.arg = FALSE),
             faa+1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)

  vecTF <- !bad0 & !is.na(p) &
           p <= dgabinom(min.support, size = size, prob = prob,
                           alter = alter,
                           pobs.a = pobs.a, byrow.arg = FALSE)
  ans[vecTF] <- min.support

  ans[!bad0 & !is.na(p) & p == 0] <- min.support
  ans[!bad0 & !is.na(p) & p == 1] <- size[!bad0 & !is.na(p) & p == 1]
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgabinom





rgabinom <-
  function(n, size, prob, alter = 0,  # NULL
           pobs.a = 0, byrow.arg = FALSE) {
    qgabinom(runif(n), size, prob = prob, alter = alter,
               pobs.a = pobs.a, byrow.arg = byrow.arg)
}  # rgabinom






if (FALSE)
gabinomial.control <-
  function(summary.HDEtest = FALSE,
           ...) {  # Overwrites the summary() default.
  list(summary.HDEtest = summary.HDEtest)
}



 gabinomial.mlm <-
  function(alter = 0,  # NULL,  # 0  #,
           zero = NULL,  # Was zero = 2 prior to 20130917
           lprob  = "logitlink",
           type.fitted = c("mean", "prob", "pobs.a", "Pobs.a"),
           imethod = 1,
           iprob = NULL
          ) {


  multiple.responses = FALSE  # Added 20190417, yettodo later

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  if (is.list(alter))
    alter <- unlist(alter)
  if (is.character(alter) && alter == "")
    alter <- NULL
  if (length(alter)) {
    if (!is.Numeric(alter, integer.valued = TRUE))
      stop("bad input for argument 'alter'")
    if (any(alter < 0))
      stop("values for argument 'alter' must be nonnegative")
    if (!identical(alter, (unique(alter))))
      stop("values of argument 'alter' must be unique")
  } else {
    stop("argument 'alter' is effectively empty; ",
         "use the family function binomialff() instead")
  }

  type.fitted <- match.arg(type.fitted,
                           c("mean", "prob", "pobs.a", "Pobs.a"))[1]
  temp7 <- if (length(alter)) paste0("pobs", alter) else NULL
  tmp3 <- c(if (length(alter))
            rep("multilogitlink",  length(alter)) else NULL,
            prob = lprob )
  names(tmp3) <- c(temp7, "prob")
 

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      any(iprob >= 1))
    stop("argument 'iprob' is out of range")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Generally-altered binomial distribution\n",
            "(GA-Binom-MLM)\n\n",
            "Links:    ", if (length(alter))
      paste("multilogitlink(cbind(",
            paste(temp7, collapse = ", "),
            ", ",
            "\n                               ",
            "1 - ", paste(temp7, collapse = " - "),
            "))", ", ", sep = ""),
            if (length(alter)) "\n          ",
            namesof("prob", lprob, earg = eprob, tag = FALSE)),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = length( .alter ) + 1)
  }), list( .zero = zero,
            .alter = alter ))),

  infos = eval(substitute(function(...) {
    alter <- ( .alter )
    temp7 <- paste("pobs", alter, sep = "")
    list(M1 = length( .alter ) + 1,
         Q1 = 1,  # Proportion when simplified
         link = .tmp3 ,
         link1parameter = if (length( .alter ))
              FALSE else TRUE,  # multilogitlink is multiparameter
         mixture.links = TRUE,
         alter = .alter ,
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE, zz possibly for later
         parameters.names = c(if (length( .alter )) temp7 else NULL,
                              "prob"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted,
           .lprob = lprob, .eprob = eprob,
           .alter = alter,
           .tmp3 = tmp3
         ))),

  initialize = eval(substitute(expression({

    if ( .multiple.responses ) {
      if (is.factor(y))
        stop("response cannot be a factor ",
             "if 'multiple.responses = TRUE'")
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

      no.successes <- w * y
      if (any(abs(w - round(w)) > 0.0001))
        warning("non-integer 'weights' (number of trials)")
      if (any(abs(no.successes - round(no.successes)) > 0.0001))
        warning("non-integer number of successes")
      if (any(y < 0 | 1 < y))
        stop("y values must be 0 <= y <= 1")
 
      if (!length(mustart) && !length(etastart))
        mustart <- matrix(colSums(no.successes) / colSums(w),
                          n, ncol(y), byrow = TRUE)
    } else {
      if (NCOL(y) == 1) {
        if (is.factor(y))
          y <- as.matrix(y != levels(y)[1])
        y <- as.matrix(y)
        w <- as.matrix(w)
        if (any(y < 0 | 1 < y))
          stop("response values 'y' must be 0 <= y <= 1")
        if (NCOL(w) != 1)
          stop("argument 'weights' has too many columns")
        no.successes <- y * w
        if (any(abs(no.successes - round(no.successes)) > 0.001))
          warning("Number of successes is not integer-valued")
        if (!length(mustart) && !length(etastart))
          mustart <- (0.5 + no.successes) / (1 + w)
      } else if (NCOL(y) == 2) {
        if (NCOL(w) != 1)
          stop("argument 'weights' has too many columns")
        if (!all(w == 1))
          stop("enter positive weights using argument 'form2'")
        if (min(y) < 0)
          stop("Negative data not allowed!")
        if (any(abs(y - round(y)) > 0.001))
          warning("Count data is not integer-valued")
        w <- cbind(y[, 1] + y[, 2])  # -> nvec
        y <- cbind(y[, 1] / w)
        no.successes <- y * w
        if (!length(mustart) && !length(etastart))
          mustart <- (0.5 + no.successes) / (1 + w)
        } else {
          stop("response 'y' must be a ",
               "vector of 0s and 1s or a vector of \n",
               "proportions and 'weight' specified,\nor a factor",
               " (first level = fail, other levels = success)",
               ",\n",
               "or a 2-column matrix where col 1 is the no. of ",
               "successes and col 2 is the no. of failures")
        }
    }  # Not multiple.responses

    alter <- ( .alter )
    lalter <- length(alter)
    M1 <- lalter + 1
    NOS <- ncoly <- NCOL(y)
    M <- NOS * M1
    extra$NOS <- extra$ncoly <- ncoly  # Number of species
    extra$M1 <- M1

    extra$pwts2 <-
    pwts2 <- if (is.null(Ym2)) 1 else matrix(Ym2, n, NOS)
    extra$w <- w  # Use possibly in @linkinv
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)


    if (!all(alter %in% c(unique(round(no.successes)))))
      stop("some values of the the 'alter' argument ",
           "not found in the response.")



    extra$y0 <- y0 <- matrix(0, n, lalter)
    for (jay in seq(lalter))
    extra$y0[, jay] <- y0[, jay] <-
      as.numeric(round(no.successes) == alter[jay])
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, lalter)
    if (any((css <- colSums(skip.these)) == 0))
      stop("some 'alter' argument values have no response values: ",
           paste(alter[css == 0], collapse = ", "))          

    
    if ( .multiple.responses ) {
      warning("this code chunk is incomplete")
      dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
      dn2 <- if (length(dn2) == M) {
        paste("E[", dn2, "]", sep = "")
      } else {
        param.names("prob", M)
      }
      predictors.names <-
        namesof(if (M > 1) dn2 else "prob",
                .link , earg = .earg, short = TRUE)
    } else {
      temp7 <- paste("pobs", alter, sep = "")
      denom.char <- paste0("1-", paste0(temp7, collapse = "-"))
      mynames1 <- paste("log(", temp7, "/(",
                        denom.char, "))", sep = "")
      mynames2 <- param.names("prob", ncoly, skip1 = TRUE)
      predictors.names <-
          c(        mynames1,
            namesof(mynames2, .lprob , earg = .eprob , tag = FALSE))[
            interleave.VGAM(M1*NOS, M1 = M1)]
    }


     if (!length(etastart)) {
      prob.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                           imu = .iprob ,  # x = x,
                           pos.only = TRUE
                           )

      phimat <- colMeans(skip.these)  # yettodo: weight this by 'w'
      phimat <- matrix(phimat, n, lalter, byrow = TRUE)
      etastart <-  multilogitlink(cbind(phimat, 1 - rowSums(phimat)))
      etastart <-
        cbind(etastart,
              theta2eta(prob.init, .lprob , earg = .eprob ))
      mustart <- NULL  # 20181203
    }
  }), list( .lprob = lprob, .eprob = eprob, .iprob = iprob,
            .imethod = imethod, .alter = alter,
            .multiple.responses = multiple.responses,
            .type.fitted = type.fitted ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Size <- c(round(extra$w))

   type.fitted <- if (length(extra$type.fitted))
                  extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "pobs.a",
                       "onempobs.a", "Pobs.a"))[1]

    alter <- ( .alter )
    M1 <- length(alter) + 1
    NOS <- NCOL(eta) / M1  # 1 for gabinomial()

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes 1 response
                inverse = TRUE)
    ynames.Pobs.a <- c(as.character(alter), "(Others)")
    dimnames(phimat) <- list(rownames(eta), ynames.Pobs.a)
    prob <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                      .lprob , earg = .eprob )
    colnames(prob) <- NULL  # Was, e.g., "logitlink(prob)".
    pobs.a <- rowSums(phimat[, -NCOL(eta), drop = FALSE])
    sump <- Sume <- matrix(0, NROW(eta), NOS)
    for (aval in alter) {
      pmf <- dbinom(aval, Size, prob)
      sump <- sump + pmf
      Sume <- Sume + pmf * aval / Size  # On a proportions scale
    }
    ans <- switch(type.fitted,
      "mean"       = (1 - pobs.a) * (prob - Sume) / (1 - sump) +
  colSums(alter * t(phimat[, -ncol(phimat), drop = FALSE] / Size)),
      "prob"       = prob,
      "pobs.a"     =     pobs.a,  # Pr(Y is altered)
      "onempobs.a" = 1 - pobs.a,  # Pr(Y is not altered)
      "Pobs.a"     =     phimat)  # matrix
    label.cols.y(ans,
                 colnames.y = if (type.fitted  == "Pobs.a")
                              ynames.Pobs.a else extra$colnames.y,
                 NOS = NOS)
  }, list( .lprob = lprob, .eprob = eprob,
           .alter = alter ))),
  last = eval(substitute(expression({

    temp.names <- c(rep_len( "multilogitlink" , lalter),
                    rep_len( .lprob , NOS))
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)

    for (ii in seq(M1*NOS - 1)) {
        misc$earg[[ii]] <- list(M = M - 1,  # M * NOS,
                                refLevel = M)  # M * NOS
    }
    misc$earg[[M1*NOS]] <- .eprob  # Last one
  }), list( .lprob = lprob, .eprob = eprob,
            .alter = alter ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    ycounts <- round(y * w)
    Size <- round(w)
    pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes 1 y
                inverse = TRUE)
    prob <- cbind(eta2theta(eta[,  NCOL(eta), drop = FALSE],
                            .lprob , earg = .eprob ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(extra$pwts2) *
        dgabinom(ycounts, size = Size, prob, log = TRUE,
                   alter = .alter , byrow.arg = FALSE,
                   pobs.a = pobs.a[, -NCOL(pobs.a), drop = FALSE])

      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .eprob = eprob,
           .alter = alter ))),
  vfamily = c("gabinomial.mlm"),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    Size <- round(extra$w)
    okay2 <- max( .alter ) <= max(Size)  # yettodo: do this once
    pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    prob <- eta2theta(eta[,  NCOL(eta), drop = FALSE],
                      .lprob , earg = .eprob )
    okay1 <- all(is.finite(prob))   &&
             all(0 < prob   & prob   < 1) &&
             all(is.finite(pobs.a)) &&
             all(0 < pobs.a & pobs.a < 1)
    okay1 && okay2
  }, list( .lprob = lprob, .eprob = eprob,
           .alter = alter ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts2 <- if (length(pwts2 <- object@extra$pwts2) > 0)
               pwts2 else
             if (length(pwts2 <- depvar(object, type = "lm2")) > 0)
               pwts2 else 1
    if (any(pwts2 != 1))
      warning("ignoring prior weights")
    if ( .multiple.responses )
      stop("cannot run simulate() when 'multiple.responses = TRUE'")
    Size <- round(weights(object, type = "prior"))

    eta <- predict(object)
    pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    prob <- eta2theta(eta[, ncol(eta)], .lprob , earg = .eprob )
    rgabinom(nsim * length(prob), Size, prob = prob,
    pobs.a = kronecker(matrix(1, nsim, 1), pobs.a[, -ncol(pobs.a)]),
               alter = .alter )
  }, list( .lprob = lprob, .eprob = eprob,
            .multiple.responses = multiple.responses,
           .alter = alter ))),


  deriv = eval(substitute(expression({
    pwts2 <- if (is.numeric(extra$pwts2)) extra$pwts2 else 1
    ycounts <- round(y * w)
    Size <- round(w)

    alter <- ( .alter )
    lalter <- length(alter)
    M1 <- lalter + 1
    NOS <- ncol(eta) / M1  # extra$NOS
    y0 <- extra$y0
    skip <- extra$skip.these
    is.altered <- rowSums(skip) > 0  # TRUE if (any) y %in% avec

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE)
    prob <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                      .lprob , earg = .eprob )


    sump <- Sume <- matrix(0, NROW(eta), NOS)
    for (aval in alter) {
      pmf <- dbinom(aval, Size, prob)
      sump <- sump + pmf
      Sume <- Sume + pmf * aval / Size  # On a proportions scale
    }

    pmf.deriv1 <- function(ycount, size, prob)
      dbinom(ycount, size, prob) * size *
      (ycount / size - prob) / (prob * (1 - prob))
    pmf.deriv2 <- function(ycount, size, prob)
      dbinom(ycount, size, prob) * size * (
      size * (((   ycount / size) /      prob) -
              (1 - ycount / size) / (1 - prob))^2 -
      (    ycount / size) /      prob^2 -
      (1 - ycount / size) / (1 - prob)^2)

    sumderiv1 <- sumderiv2 <- matrix(0, NROW(eta), NOS)
    for (aval in alter) {
      sumderiv1 <- sumderiv1 + pmf.deriv1(aval, Size, prob)
      sumderiv2 <- sumderiv2 + pmf.deriv2(aval, Size, prob)
    }

    pobs.a <- rowSums(phimat[, -NCOL(phimat), drop = FALSE])
    onempobs.a <- 1 - pobs.a

    dl.deta <- skip - phimat[, -M, drop = FALSE]

    dl.dprob <- (1 - is.altered) * (
      c(w) *      y         /      prob  -
      c(w) * (1 - y       ) / (1 - prob) +
                sumderiv1 / (1 - sump))
    dprob.deta <- dtheta.deta(prob, .lprob , earg = .eprob )
    ans <- cbind(c(pwts2) * c(w) * dl.deta,
                 c(pwts2) *        dl.dprob * dprob.deta)
    ans
  }), list( .lprob = lprob, .eprob = eprob,
            .alter = alter ))),


  weight = eval(substitute(expression({
    MM12 <- M1 * (M1 + 1) / 2

      ned2l.dprob2 <- onempobs.a *
        ((Size / prob^2) *
         (prob - Sume) / (1 - sump) +
         (1 - (prob - Sume) / (1 - sump)) * Size / (1 - prob)^2 -
         sumderiv2 / (1 - sump) - (sumderiv1 / (1 - sump))^2)

    wz4 <- matrix(0.0, n, MM12)  # A full matrix
    use.refLevel <- M
    if (lalter > 0) {
      if (lalter == 1) {
        wz4[, 1] <-      phimat[, 1] * (1 - phimat[, 1])
      } else {
        index <- iam(NA, NA, M-1, both = TRUE, diag = TRUE)
        wz4 <- -phimat[, index$row] * phimat[, index$col]
        wz4[, 1:lalter] <- wz4[, 1:lalter] + phimat[, -M]
      }
    }
    wz4 <- as.matrix(wz4)  #  Needed when lalter == 1
    wz <- wz.merge(c(pwts2) * c(w) * wz4,
                   c(pwts2) * ned2l.dprob2 * dprob.deta^2,
                   M1 = M1 - 1, M2 = 1)  # rm.trailing.cols = FALSE
    wz
  }), list( .alter = alter ))))
}  # gabinomial.mlm






muxtXAX <- function(A, X, M) {

  if ((n <- nrow(X)) != nrow(A))
    stop("arguments 'X' and 'A' are not conformable")

  ans1 <- array(0, c(n, M, M))
  for (eye in 1:M)
    for (jay in 1:M)
      for (kay in 1:M)
        ans1[, eye, jay] <- ans1[, eye, jay] +
          A[, iam(eye, kay, M)] * X[, iam(kay, jay, M)]

  ans2 <- matrix(0, n, dimm(M))
  for (eye in 1:M)
    for (jay in eye:M)
      for (kay in 1:M)
        ans2[, iam(eye, jay, M)] <- ans2[, iam(eye, jay, M)] +
        X[, iam(kay, eye, M)] * ans1[, kay, jay]
  ans2
}  # muxtXAX















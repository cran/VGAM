# These functions are
# Copyright (C) 1998-2020 T.W. Yee, University of Auckland.
# All rights reserved.


















 moments.zeta.gait <-
  function(shape.p,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           byrow.ai = FALSE,  # For pobs.a and pstr.i
           shape.a = shape.p, shape.i = shape.p,
           mlm = TRUE,  # mlm = FALSE iff mixture = TRUE
           type.fitted = "All") {  # or "mean"
  NOS <- 1
  nnn <- length(shape.p)
  rmlife <- if (is.finite(max.support)) zeta(shape.p) * (1 -
    pzeta(max.support, shape.p - 1)) / zeta(shape.p + 1) else
    numeric(nnn)
  rmlife[shape.p <= 1] <- NA  # NA or Inf, not sure

  lhs.prob <- pzeta(max.support, shape.p)  # 1 - Pr(upper tail)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)

  sumt.p <- matrix(0, nnn, NOS)  # Does not include upper RHS tail
  SumT.p <- matrix(rmlife, nnn, NOS)  # Does include upper RHS tail
  if (ltrunc)
    for (tval in truncate) {
      pmf.p <- dzeta(tval, shape.p)
      sumt.p <- sumt.p + pmf.p  # Need tval <= max.support
      SumT.p <- SumT.p + pmf.p * tval
    }

  use.pobs.a <- use.pstr.i <- matrix(0, nnn, 1)  # So rowSums() works.
  aprd <- 0  # aprd is an innerprod
  SumA.x <-  # For innerprod (aprd) only
  SumA.a <- suma.a <-
  SumA.p <- suma.p <- matrix(0, nnn, NOS)
  if (mlm) {  # use.pobs.a and use.pstr.i have the right dimensions.
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
      pmf.p <- dzeta(aval, shape.p)
      pmf.a <- dzeta(aval, shape.a)
      suma.p <- suma.p + pmf.p
      SumA.p <- SumA.p + pmf.p * aval
      suma.a <- suma.a + pmf.a
      SumA.a <- SumA.a + pmf.a * aval
      SumA.x <- SumA.x + pmf.x * aval
    }  # for jay
    aprd <- if (mlm) SumA.x else use.pobs.a * SumA.a / suma.a
  }  # lalter

  iprd <- 0  # iprd is an innerprod
  SumI.x <-  # For innerprod (iprd) only
  sumi.i <- SumI.i <-
  sumi.p <- SumI.p <- matrix(0, nnn, NOS)
  if (linfla) {
    for (jay in seq_len(linfla)) {
      ival <- inflate[jay]
      pmf.x <- if (mlm) use.pstr.i[, jay] else rep(0, nnn)
      pmf.p <- dzeta(ival, shape.p)  # Changed
      pmf.i <- dzeta(ival, shape.i)  # Changed
      sumi.p <- sumi.p + pmf.p
      SumI.p <- SumI.p + pmf.p * ival
      sumi.i <- sumi.i + pmf.i
      SumI.i <- SumI.i + pmf.i * ival
      SumI.x <- SumI.x + pmf.x * ival
    }  # for jay
    iprd <- if (mlm) SumI.x else use.pstr.i * SumI.i / sumi.i
  }  # linfla

  use.this <- if (mlm)
    (1 - rowSums(use.pobs.a) - rowSums(use.pstr.i)) else
    (1 - use.pobs.a - use.pstr.i)

  mean.true.p <- zeta(shape.p) / zeta(shape.p + 1)
  mean.true.p[shape.p <= 1] <- NA  # Does not exist
  themean <- aprd + iprd + use.this *
    (mean.true.p - SumA.p - SumT.p) / (lhs.prob - suma.p - sumt.p)
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
}  # moments.zeta.gait







 gaitzeta.mix <-
  function(alter = NULL,
           inflate = NULL,  
           truncate = NULL,  #  max.support = Inf,
           zero = c("pobs.a", "pstr.i"),  # Pruned later if necessary
           eq.ap = FALSE,  # TRUE applies to the intercept
           eq.ip = FALSE,  # TRUE applies to the intercept
           lshape.p = "loglink",
           lpobs.a = "logitlink",
           lshape.a = "loglink",
           lpstr.i = "logitlink",
           lshape.i = "loglink",
         type.fitted = c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
                         "prob.a", "prob.i", "prob.t", "lhs.prob"),
           imethod = 1,
           gshape.p = exp(-7 * ppoints(12)),
           ishape.p = NULL,  # Higher is better than lower
           ishape.a = NULL, ishape.i = NULL,
           ipobs.a = NULL,  # 0.25, 
           ipstr.i = NULL,  # 0.25, 
           ishrinkage = 0.95,
           probs.y = 0.35) {
  max.support = Inf  # For now; zz
  lowsup <- 1
  semigait.errorcheck(alter, inflate, truncate, max.support,
                      min.support = lowsup)

  if (!is.logical(eq.ap) || length(eq.ap) != 1)
    stop("argument 'eq.ap' must be a single logical")
  if (!is.logical(eq.ip) || length(eq.ip) != 1)
    stop("argument 'eq.ip' must be a single logical")

  lshape.p <- as.list(substitute(lshape.p))
  eshape.p <- link2list(lshape.p)
  lshape.p <- attr(eshape.p, "function.name")
  lpobs.a <- as.list(substitute(lpobs.a))  # \omega
  epobs.a <- link2list(lpobs.a)
  lpobs.a <- attr(epobs.a, "function.name")
  lshape.a <- as.list(substitute(lshape.a))
  eshape.a <- link2list(lshape.a)
  lshape.a <- attr(eshape.a, "function.name")

  lpstr.i <- as.list(substitute(lpstr.i))  # \phi
  epstr.i <- link2list(lpstr.i)
  lpstr.i <- attr(epstr.i, "function.name")
  lshape.i <- as.list(substitute(lshape.i))
  eshape.i <- link2list(lshape.i)
  lshape.i <- attr(eshape.i, "function.name")

  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)

  if (is.vector(zero) && is.character(zero) && length(zero) == 2) {
    if (lalter && linfla == 0)
      zero <- setdiff(zero, "pstr.i")
    if (linfla && lalter == 0)
      zero <- setdiff(zero, "pobs.a")
  }

  if (lalter + linfla == 0) zero <- NULL
  lshape.p.save <- lshape.p
  if (lalter + linfla + ltrunc == 0 && is.infinite(max.support))
    return(eval(substitute(
           zetaff(lshape = .lshape.p.save , zero = .zero ),
           list( .lshape.p.save = lshape.p.save, .zero = zero))))


  if (lalter == 1 && eq.ap)
    warning("Less than 2 altered values and hence 1 shape, so ",
            "setting 'eq.ap = TRUE' is meaningless")
  if (linfla == 1 && eq.ip)
    warning("Less than 2 inflated values and hence 1 shape, so ",
            "setting 'eq.ip = TRUE' is meaningless")

  type.fitted <- match.arg(type.fitted,
    c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
      "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]

  tmp3.TF <- c(TRUE, lalter > 0, lalter > 1, linfla > 0, linfla > 1)
  tmp3 <- c(shape.p = lshape.p,  # Version 2
            pobs.a   = lpobs.a, shape.a = lshape.a,
            pstr.i   = lpstr.i, shape.i = lshape.i)  # [tmp3.TF]
      
  blurb1 <- "Z"
  if (lalter) blurb1 <- "Generally-altered Z"
  if (linfla) blurb1 <- "Generally-inflated Z"
  if (ltrunc) blurb1 <- "Generally-truncated Z"
  if ( lalter &&  linfla && !ltrunc)
    blurb1 <- "Generally-altered and -inflated Z"
  if ( lalter && !linfla &&  ltrunc)
    blurb1 <- "Generally-altered and -truncated Z"
  if (!lalter &&  linfla &&  ltrunc)
    blurb1 <- "Generally-inflated and -truncated Z"
  if ( lalter &&  linfla &&  ltrunc)
    blurb1 <- "Generally-altered, -inflated and -truncated Z"

  new("vglmff",
  blurb = c(blurb1, "eta regression\n",
            "(GAIT-Zeta(shape.p)-Zeta(shape.a)-Zeta(shape.i) ",
            "mixture generally)\n\n",
            "Links:    ",
            namesof("shape.p",  lshape.p, eshape.p, tag = FALSE),
            if (lalter > 0) c(", ",
            namesof("pobs.a",    lpobs.a,  epobs.a, tag = FALSE)),
            if (lalter > 1) c(", ",
            namesof("shape.a",  lshape.a, eshape.a, tag = FALSE)),
            if (lalter && linfla) ", \n        ",
            if (linfla > 0) c(  if (lalter) "  " else ", ",
            namesof("pstr.i",    lpstr.i,  epstr.i, tag = FALSE)),
            if (linfla > 1) c(", ",
            namesof("shape.i",  lshape.i, eshape.i, tag = FALSE)),
            "\n"),
  constraints = eval(substitute(expression({

    use.mat <- if ( ( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,1,0,1,  0,1,0,0,0,  0,0,0,1,0), 5, 3) else
    if (!( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,0,0,1,  0,1,0,0,0,  0,0,1,0,0,
               0,0,0,1,0), 5, 4) else
    if ( ( .eq.ap ) && !( .eq.ip ))
      matrix(c(1,0,1,0,0,  0,1,0,0,0,  0,0,0,1,0,
               0,0,0,0,1), 5, 4) else
      diag(5)

    use.mat <- use.mat[ .tmp3.TF , , drop = FALSE]
    M1 <- nrow(use.mat)
    for (jay in ncol(use.mat):1)
      if (all(use.mat[, jay] == 0))
        use.mat <- use.mat[, -jay, drop = FALSE]
    
    constraints <- cm.VGAM(use.mat, x = x,
                           bool = .eq.ap | .eq.ip ,
                           constraints = constraints,
                           apply.int = TRUE)  # FALSE

    if ( .lalter + .linfla > 0)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = M1)
  }), list( .zero = zero, .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
            .lalter = lalter, .linfla = linfla,
            .eq.ap = eq.ap, .eq.ip = eq.ip ))),

  infos = eval(substitute(function(...) {
    list(Q1 = 1,
         M1 = sum( .tmp3.TF ),
         link = c( .tmp3 )[ .tmp3.TF ] ,
         link1parameter = TRUE,
         mixture.links = FALSE,
         alter = as.vector( .alter ),
         inflate = as.vector( .inflate ),
         truncate = as.vector( .truncate ),
         max.support = as.vector( .max.support ), 
         Support  = c( .lowsup , Inf, 1),
         eq.ap = .eq.ap , eq.ip = .eq.ip ,
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE,
         parameters.names = names(c( .tmp3 )[ .tmp3.TF ]),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero, .tmp3 = tmp3, .tmp3.TF = tmp3.TF,
           .type.fitted = type.fitted, .lowsup = lowsup,
           .eq.ap = eq.ap, .eq.ip = eq.ip,
           .alter = alter, .inflate = inflate,
           .truncate = truncate, .max.support = max.support ))),

  initialize = eval(substitute(expression({
    truncate <- as.vector( .truncate )
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    M1 <- sum( .tmp3.TF )
    NOS <- NCOL(y)  # Only 1 currently
    M <- NOS * M1
    tmp3 <- ( .tmp3 )
    tmp3.TF <- ( .tmp3.TF )
    ntmp3 <- names(tmp3)

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    glist <- y.gait.check(alter, inflate, truncate, .max.support , y)
    css.a <- glist$css.a
    css.i <- glist$css.i
    extra$skip.a <- glist$skip.a
    extra$skip.i <- glist$skip.i

    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- ( .type.fitted )
    extra$colnames.y  <- colnames(y)
    extra$M1 <- M1

   predictors.names <-
    c(namesof(ntmp3[1], .lshape.p , earg = .eshape.p , tag = FALSE),
      namesof(ntmp3[2], .lpobs.a  , earg = .epobs.a  , tag = FALSE),
      namesof(ntmp3[3], .lshape.a , earg = .eshape.a , tag = FALSE),
      namesof(ntmp3[4], .lpstr.i  , earg = .epstr.i  , tag = FALSE),
      namesof(ntmp3[5], .lshape.i , earg = .eshape.i , tag = FALSE))[
      tmp3.TF]  # Relevant subset



    

    if (!length(etastart)) {
      shape.p.init <- if (length( .ishape.p )) .ishape.p else {
        zetaff.Loglikfun <- function(shapeval, y, x, w, extraargs) {
          sum(c(w) * dzeta(x = y, shape = shapeval, log = TRUE))
        }
        shape.p.grid <- ( .gshape.p )
        grid.search(shape.p.grid, objfun = zetaff.Loglikfun,
                    y = y, w = w)
      }
      shape.p.init <- rep(shape.p.init, length = n)

      pobs.a.init <-  # Unneeded? Low is best
      pstr.i.init <- 0.025 / (1 + lalter + linfla)
      if (lalter > 0) {
        pobs.a.init <-  rep(if (length( .ipobs.a )) .ipobs.a else
                            sum(css.a) / n, length = n)
      }
      if (lalter > 1)
        shape.a.init <- if (length( .ishape.a ))
          rep( .ishape.a , n) else shape.p.init

      if (linfla > 0)
        pstr.i.init <-  rep(if (length( .ipstr.i )) .ipstr.i else
                            0.5 * sum(css.i) / n, length = n)
      if (linfla > 1)
        shape.i.init <- if (length( .ishape.i ))
          rep( .ishape.i , n) else shape.p.init

      while (any((vecTF <- pobs.a.init + pstr.i.init > 0.95))) {
        pobs.a.init[vecTF] <- 0.875 * pobs.a.init[vecTF]
        pstr.i.init[vecTF] <- 0.875 * pstr.i.init[vecTF]
      }

      etastart <-
        cbind(theta2eta(shape.p.init, .lshape.p , earg = .eshape.p ),
              if (lalter <= 0) NULL else
              theta2eta(pobs.a.init,  .lpobs.a  , earg = .epobs.a  ),
              if (lalter <= 1) NULL else
              theta2eta(shape.a.init, .lshape.a , earg = .eshape.a ),
              if (linfla <= 0) NULL else
              theta2eta(pstr.i.init,  .lpstr.i  , earg = .epstr.i  ),
              if (linfla <= 1) NULL else
              theta2eta(shape.i.init, .lshape.i , earg = .eshape.i ))
    }
  }), list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
    .ishape.p = ishape.p, .ishape.i = ishape.i, .ipstr.i = ipstr.i,
                          .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                          .eshape.a = eshape.a, .epobs.a = epobs.a,
    .ishape.p = ishape.p, .ishape.a = ishape.a, .ipobs.a = ipobs.a,
    .gshape.p = gshape.p, 
    .ishrinkage = ishrinkage, .probs.y = probs.y,
    .alter = alter, .inflate = inflate,
    .truncate = truncate, .max.support = max.support,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .imethod = imethod, .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <-
     if (length(extra$type.fitted)) extra$type.fitted else {
       warning("cannot find 'type.fitted'. Returning the 'mean'.")
                     "mean"
     }
   
    type.fitted <-
      match.arg(type.fitted,
                c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
                  "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]
    eta <- as.matrix(eta)  # Needed when linfla == 0
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate) 
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    max.support <- ( .max.support )
    M1 <- sum( .tmp3.TF )
    if ((NOS <- NCOL(eta) / M1) != 1)
      stop("Currently NOS must be 1")

    iptr <- 1  # Unpack
    shape.p <- eta2theta(eta[, iptr], .lshape.p , earg = .eshape.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a  , earg = .epobs.a )
    shape.a <- if (lalter <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.a , earg = .eshape.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i  , earg = .epstr.i )
    shape.i <- if (linfla <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.i , earg = .eshape.i )
    Bits <- moments.zeta.gait(shape.p = shape.p, mlm = FALSE,
                              alter = alter, inflate = inflate,
                              truncate = truncate,
                              max.support = max.support,
                              pobs.a = pobs.a, shape.a = shape.a,
                              pstr.i = pstr.i, shape.i = shape.i)
    if (type.fitted == "Pobs.a") {
      if (lalter == 0) stop("no altered values!")
      proportion.mat <-
        dzeta(matrix(alter,   NROW(eta), lalter, byrow = TRUE),
              matrix(shape.a, NROW(eta), lalter)) / c(Bits[["suma.a"]])
    }
    if (type.fitted == "Pstr.i") {
      if (linfla == 0) stop("no inflated values!")
      proportion.mat <-
        dzeta(matrix(inflate, nrow(eta), linfla, byrow = TRUE),
              matrix(shape.i, nrow(eta), linfla)) / c(Bits[["sumi.i"]])
    }
  
    ans <- switch(type.fitted,
      "mean"     = Bits[["mean"]],
      "pobs.a"   = pobs.a,
      "pstr.i"   = pstr.i,
      "Pobs.a"   = pobs.a * proportion.mat,  # matrix
      "Pstr.i"   = c(pstr.i) * proportion.mat +
                   c(1 - pobs.a - pstr.i) *
           dzeta(matrix(inflate, nrow(eta), linfla, byrow = TRUE),
                 matrix(shape.p, nrow(eta), linfla)) / c(
           Bits[["lhs.prob"]] - Bits[["suma.p"]] - Bits[["sumt.p"]]),
      "prob.a"   = Bits[["suma.p"]],
      "prob.i"   = Bits[["sumi.p"]],
      "prob.t"   = Bits[["sumt.p"]],
      "lhs.prob" = Bits[["lhs.prob"]])
   ynames.Pobs.a <- as.character(alter)    # Works with NULLs
   ynames.Pstr.i <- as.character(inflate)  # Works with NULLs
   label.cols.y(ans,
        colnames.y = switch(type.fitted,
                            "Pobs.a" = ynames.Pobs.a,
                            "Pstr.i" = ynames.Pstr.i,
                            extra$colnames.y),
        NOS = NOS)
  }, list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                          .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                          .eshape.a = eshape.a, .epobs.a = epobs.a,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .truncate = truncate, .max.support = max.support, 
    .alter = alter, .inflate = inflate ))),

  last = eval(substitute(expression({
    misc$link  <- c( .tmp3 )[ .tmp3.TF ]
    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    misc$earg[[(iptr <- 1)]] <- .eshape.p  #
    if (lalter > 0)
      misc$earg[[(iptr <- iptr + 1)]] <- .epobs.a  #
    if (lalter > 1)
      misc$earg[[(iptr <- iptr + 1)]] <- .eshape.a  #
    if (linfla > 0)
      misc$earg[[(iptr <- iptr + 1)]] <- .epstr.i  #
    if (linfla > 1)
      misc$earg[[(iptr <- iptr + 1)]] <- .eshape.i  #
  }), list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                            .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                            .eshape.a = eshape.a, .epobs.a = epobs.a,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .alter = alter, .inflate = inflate ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    lalter <- length(as.vector( .alter ))
    linfla <- length(as.vector( .inflate ))
    if (lalter + linfla == 0)
      eta <- as.matrix(eta)
    iptr <- 1  # Unpack
    shape.p <- eta2theta(eta[, iptr], .lshape.p , earg = .eshape.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a  , earg = .epobs.a )
    shape.a <- if (lalter <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.a , earg = .eshape.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i  , earg = .epstr.i )
    shape.i <- if (linfla <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.i , earg = .eshape.i )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitzeta(y, shape.p, shape.a = shape.a, shape.i = shape.i,
                  pobs.mix.a = pobs.a, pstr.mix.i = pstr.i, log = TRUE,
                  truncate = .truncate , max.support = .max.support ,
                  alter.mix = .alter , inflate.mix = .inflate )
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                          .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                          .eshape.a = eshape.a, .epobs.a = epobs.a,
    .truncate = truncate, .max.support = max.support, 
    .alter = alter, .inflate = inflate ))),
  vfamily = c("gaitzeta.mix"),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    lalter <- length(as.vector( .alter ))
    linfla <- length(as.vector( .inflate ))
    if (lalter + linfla == 0)
      eta <- as.matrix(eta)
    iptr <- 1  # Unpack
    shape.p <- eta2theta(eta[, iptr], .lshape.p , earg = .eshape.p )
    pobs.a   <- if (lalter <= 0) rep(0.25, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a , earg = .epobs.a )
    shape.a <- if (lalter <= 1) shape.p else
    eta2theta(eta[, (iptr <- iptr + 1)], .lshape.a  , earg = .eshape.a )
    pstr.i   <- if (linfla <= 0) rep(0.25, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i , earg = .epstr.i )
    shape.i <- if (linfla <= 1) shape.p else
    eta2theta(eta[, (iptr <- iptr + 1)], .lshape.i  , earg = .eshape.i )

    okay1 <- all(is.finite(shape.p)) && all(0 < shape.p) &&
             all(is.finite(shape.a)) && all(0 < shape.a) &&
             all(is.finite(shape.i)) && all(0 < shape.i) &&
             all(is.finite(pobs.a)) && all(0 < pobs.a & pobs.a < 1) &&
             all(is.finite(pstr.i)) && all(0 < pstr.i & pstr.i < 1) &&
             all(pobs.a + pstr.i < 1)
    okay1
  }, list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                          .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                          .eshape.a = eshape.a, .epobs.a = epobs.a,
    .truncate = truncate, .max.support = max.support, 
    .alter = alter, .inflate = inflate ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    lalter <- length(as.vector( .alter ))
    linfla <- length(as.vector( .inflate ))
    iptr <- 1  # Unpack
    shape.p <- eta2theta(eta[, iptr], .lshape.p , earg = .eshape.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a  , earg = .epobs.a )
    shape.a <- if (lalter <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.a , earg = .eshape.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i  , earg = .epstr.i )
    shape.i <- if (linfla <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.i , earg = .eshape.i )
    rgaitzeta(nsim * length(shape.p), shape.p = shape.p,
             shape.a = shape.a, shape.i = shape.i,
             pobs.mix.a = pobs.a, pstr.mix.i = pstr.i,
             alter.mix = .alter , inflate.mix = .inflate ,
             truncate = .truncate , max.support = .max.support )
  }, list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                          .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                          .eshape.a = eshape.a, .epobs.a = epobs.a,
    .truncate = truncate, .max.support = max.support,
    .alter = alter, .inflate = inflate ))),
  deriv = eval(substitute(expression({
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    M1 <- sum( .tmp3.TF )
    NOS <- NCOL(eta) / M1  # extra$NOS
    max.support <- as.vector( .max.support )

    ans2 <- matrix(0, n, M1)
    Offset.i <- if (lalter == 1) 2 else if (lalter > 1) 3 else 1


    is.altered <- if (lalter)
      rowSums(extra$skip.a) > 0 else rep(FALSE, n)
    is.inflated <- if (linfla)
      rowSums(extra$skip.i) > 0 else rep(FALSE, n)


    if (lalter + linfla == 0)
      eta <- as.matrix(eta)
    iptr <- 1  # Unpack
    shape.p <- eta2theta(eta[, iptr], .lshape.p , earg = .eshape.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a  , earg = .epobs.a )
    shape.a <- if (lalter <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.a , earg = .eshape.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i  , earg = .epstr.i )
    shape.i <- if (linfla <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.i , earg = .eshape.i )

    bits <- moments.zeta.gait(shape.p, alter = alter,
                              inflate = inflate, truncate = truncate,
                              max.support = max.support,
                              pobs.a = pobs.a, pstr.i = pstr.i,
                              shape.a = shape.a, shape.i = shape.i,
                              mlm = FALSE)
    Denom <- c(bits[["lhs.prob"]] - bits[["suma.p"]] - bits[["sumt.p"]])

    pmf.deriv1 <- function(y, shape) {
      -dzeta(y, shape) * (log(y) +
        zeta(shape + 1, deriv = 1) / zeta(shape + 1))
    }
    pmf.deriv2 <- function(y, shape) {
      tmp2 <- zeta(shape + 1, deriv = 1) / zeta(shape + 1)
      dzeta(y, shape) * ((log(y) + tmp2)^2 -
        zeta(shape + 1, deriv = 2) / zeta(shape + 1) + tmp2^2)
    }
    
    sumderiv1a.p <- sumderiv2a.p <- matrix(0, n, NOS)
    if (lalter > 0) {
      deriv0matrix.a <-
      deriv1matrix.a <- deriv2matrix.a <- matrix(0, n, lalter)
      sumderiv1a.a <- sumderiv2a.a <- matrix(0, n, NOS)
      for (jay in seq(lalter)) {
        aval <- alter[jay]
        sumderiv1a.p <- sumderiv1a.p + pmf.deriv1(aval, shape.p)
        sumderiv2a.p <- sumderiv2a.p + pmf.deriv2(aval, shape.p)
        sumderiv1a.a <- sumderiv1a.a + pmf.deriv1(aval, shape.a)
        sumderiv2a.a <- sumderiv2a.a + pmf.deriv2(aval, shape.a)
        pmf.a <- dzeta(aval, shape.a)
        deriv0matrix.a[, jay] <- pmf.a
        deriv1matrix.a[, jay] <- pmf.deriv1(aval, shape.a)
        deriv2matrix.a[, jay] <- pmf.deriv2(aval, shape.a)
      }
    }  # lalter > 0





    if (linfla) {
      deriv0matrix.i <-  # wrt inflated distribution
      deriv1matrix.i <- deriv2matrix.i <- matrix(0, n, linfla)
      Deriv0Matrix.i <-  # wrt parent distribution
      Deriv1Matrix.i <- Deriv2Matrix.i <- matrix(0, n, linfla)
      if (linfla > 0)
        for (jay in seq(linfla)) {
          ival <- inflate[jay]
          pmf.i <- dzeta(ival, shape.i)
          deriv0matrix.i[, jay] <- pmf.i
          deriv1matrix.i[, jay] <- pmf.deriv1(ival, shape.i)
          deriv2matrix.i[, jay] <- pmf.deriv2(ival, shape.i)
          pmf.p <- dzeta(ival, shape.p)
          Deriv0Matrix.i[, jay] <- pmf.p
          Deriv1Matrix.i[, jay] <- pmf.deriv1(ival, shape.p)
          Deriv2Matrix.i[, jay] <- pmf.deriv2(ival, shape.p)  # /pmf.p
        }  # jay
    }  # linfla




    sumderiv1t.a <- sumderiv2t.a <-
    sumderiv1t.i <- sumderiv2t.i <-
    sumderiv1t.p <- sumderiv2t.p <- matrix(0, n, NOS)
    if (ltrunc)
      for (tval in truncate) {
        sumderiv1t.p <- sumderiv1t.p + pmf.deriv1(tval, shape.p)
        sumderiv2t.p <- sumderiv2t.p + pmf.deriv2(tval, shape.p)
        sumderiv1t.a <- sumderiv1t.a + pmf.deriv1(tval, shape.a)
        sumderiv2t.a <- sumderiv2t.a + pmf.deriv2(tval, shape.a)
        sumderiv1t.i <- sumderiv1t.i + pmf.deriv1(tval, shape.i)
        sumderiv2t.i <- sumderiv2t.i + pmf.deriv2(tval, shape.i)
      }

    onembothprobs <- 1 - pobs.a - pstr.i  # 1-Pr(altered)-Pr(inflated)
    if (is.finite(max.support)) {
      stop("20200323; intractable")
    }


    fred0.p <- zeta(shape.p + 1)
    fred1.p <- zeta(shape.p + 1, deriv = 1)
    dl.dshape.p <- -log(y) - fred1.p / fred0.p  # Usual formula

    zero0n <- rep(0, n)
    dl.dshape.a <-
    dl.dshape.i <- zero0n  # Replace some elts below
    dl.dshape.p <- dl.dshape.p +
      (sumderiv1a.p + sumderiv1t.p) / Denom  # \notin A, I, T
    dl.dshape.p[is.altered] <- 0
    dl.dpobs.a <-  # Replace some elts below
    dl.dpstr.i <- -1 / onembothprobs  # \notin A, I, T
    Delta <- matrix(0, n, NOS)  # If linfla == 0

    if (linfla > 0) {
      d0A.io <- deriv0matrix.i / c(bits[["sumi.i"]])
      d0B.pi <- Deriv0Matrix.i / Denom
      Delta <- pstr.i * d0A.io + onembothprobs * d0B.pi


      d1A.io <- deriv1matrix.i / c(bits[["sumi.i"]]) -
                deriv0matrix.i * rowSums(deriv1matrix.i) / c(
                    bits[["sumi.i"]])^2
      d2A.io <- deriv2matrix.i / c(bits[["sumi.i"]]) -
                2 * deriv1matrix.i * rowSums(deriv1matrix.i) / c(
                    bits[["sumi.i"]])^2 -
                deriv0matrix.i * rowSums(deriv2matrix.i) / c(
                    bits[["sumi.i"]])^2 +
                2 * deriv0matrix.i * rowSums(deriv1matrix.i)^2 / c(
                    bits[["sumi.i"]])^3
      d1B.pi <- Deriv1Matrix.i / Denom +
        Deriv0Matrix.i * c(sumderiv1t.p + sumderiv1a.p) / Denom^2
      d2B.pi <- Deriv2Matrix.i / Denom +
        2 * Deriv1Matrix.i * c(sumderiv1t.p + sumderiv1a.p) / Denom^2 +
        Deriv0Matrix.i * c(sumderiv2t.p + sumderiv2a.p) / Denom^2 +
        2 * Deriv0Matrix.i * (c(sumderiv1t.p + sumderiv1a.p)^2) / Denom^3
    }  # linfla > 0
    
    if (lalter) {
      dl.dpobs.a[is.altered] <- 1 / pobs.a[is.altered]
      dl.dpstr.i[is.altered] <- 0

      d0A.al <- deriv0matrix.a / c(bits[["suma.a"]])
      d1A.al <- deriv1matrix.a / c(bits[["suma.a"]]) -
                deriv0matrix.a * rowSums(deriv1matrix.a) / c(
                    bits[["suma.a"]])^2
      d2A.al <- deriv2matrix.a / c(bits[["suma.a"]]) -
                2 * deriv1matrix.a * rowSums(deriv1matrix.a) / c(
                    bits[["suma.a"]])^2 -
                deriv0matrix.a * rowSums(deriv2matrix.a) / c(
                    bits[["suma.a"]])^2 +
                2 * deriv0matrix.a * rowSums(deriv1matrix.a)^2 / c(
                    bits[["suma.a"]])^3

      for (jay in seq(lalter)) {
        aval <- alter[jay]
        is.alte.j <- extra$skip.a[, jay]  # Logical vector
        tmp2 <- d1A.al[, jay] / d0A.al[, jay]
        dl.dshape.a[is.alte.j] <- tmp2[is.alte.j]
      }  # jay
    }  # lalter

    if (linfla > 0)
      for (jay in seq(linfla)) {
        ival <- inflate[jay]
        is.infl.j <- extra$skip.i[, jay]  # Logical vector
        tmp7 <- onembothprobs * d1B.pi[, jay] / Delta[, jay]
        dl.dshape.p[is.infl.j] <- tmp7[is.infl.j]
        tmp2 <- -d0B.pi[, jay] / Delta[, jay]
        dl.dpobs.a[is.infl.j] <- tmp2[is.infl.j]
        tmp8 <- (d0A.io[, jay] - d0B.pi[, jay]) / Delta[, jay]
        dl.dpstr.i[is.infl.j] <- tmp8[is.infl.j]
        if (linfla > 1) {
          tmp2 <- pstr.i * d1A.io[, jay] / Delta[, jay]
          dl.dshape.i[is.infl.j] <- tmp2[is.infl.j]
        }
      }  # jay


    dshape.p.deta <- dtheta.deta(shape.p, .lshape.p , .eshape.p )
    ans2[, 1] <- dl.dshape.p * dshape.p.deta
    dpobs.a.deta <- dtheta.deta(pobs.a, .lpobs.a , .epobs.a )
    if (lalter > 0) {
      ans2[, 2] <- dl.dpobs.a * dpobs.a.deta
    }
    if (lalter > 1) {
      dshape.a.deta <- dtheta.deta(shape.a, .lshape.a , .eshape.a )
      ans2[, 3] <- dl.dshape.a * dshape.a.deta
    }
    if (linfla > 0) {
      dpstr.i.deta <- dtheta.deta(pstr.i, .lpstr.i , .epstr.i )
      ans2[, Offset.i + 1] <- dl.dpstr.i * dpstr.i.deta
    }
    if (linfla > 1) {
      dshape.i.deta <- dtheta.deta(shape.i, .lshape.i , .eshape.i )
      ans2[, Offset.i + 2] <- dl.dshape.i * dshape.i.deta
    }
    c(w) * ans2
  }), list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                          .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                          .eshape.a = eshape.a, .epobs.a = epobs.a,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .truncate = truncate, .max.support = max.support, 
    .alter = alter , .inflate = inflate ))),

  weight = eval(substitute(expression({
    onemrsDelta <- onembothprobs *  # (onempstr.i - pobs.a) *
                   (1 - c(bits[["sumi.p"]]) / Denom)

    ned2l.dshape.p.shape.a <- zero0n  # Final; nothing to do
    ned2l.dpobs.a.shape.p  <- zero0n  # mB overwritten below [4279]
    ned2l.dpobs.a.shape.a  <- zero0n  # Final; nothing to do
    ned2l.dpobs.a.pstr.i   <- zero0n  # mB overwritten below
    ned2l.dpobs.a.shape.i  <- zero0n  # mB overwritten below
    ned2l.dpstr.i.shape.p  <- zero0n  # mB overwritten below
    ned2l.dpstr.i.shape.a  <- zero0n  # Final; nothing to do
    ned2l.dpstr.i.shape.i  <- zero0n  # mB overwritten below
    ned2l.dshape.a.shape.i <- zero0n  # Final; nothing to do

    ned2l.dshape.p2 <-
      onemrsDelta * (zeta(shape.p + 1, deriv = 2) / fred0.p -
                    (fred1.p / fred0.p)^2 -
                     c(sumderiv2t.p + sumderiv2a.p) / Denom -
                    (c(sumderiv1t.p + sumderiv1a.p) / Denom)^2)
    if (linfla > 0)
      ned2l.dshape.p2 <- ned2l.dshape.p2 + onembothprobs *
        rowSums(onembothprobs * d1B.pi^2 / Delta - d2B.pi)

    if (lalter > 0) {
      wz22 <- if ( .lpobs.a == "logitlink" && linfla == 0) {
        pobs.a * (1 - pobs.a)
      } else {
        ned2l.dpobs.a2 <- 1 / pobs.a +
          onemrsDelta / onembothprobs^2 + (if (linfla > 0)
          rowSums(d0B.pi^2 / Delta) else 0)
        ned2l.dpobs.a2 * dpobs.a.deta^2
      }
    }
    if (lalter > 1)
      ned2l.dshape.a2 <- pobs.a * (
        rowSums(deriv1matrix.a^2  / deriv0matrix.a) / c(
        bits[["suma.a"]]) -
        (c(sumderiv1a.a) / c(bits[["suma.a"]]))^2)

    if (linfla > 0)
      ned2l.dpstr.i2 <-
        onemrsDelta / onembothprobs^2 +
        (if (linfla > 0)  # zz > 1 ??
        rowSums((deriv0matrix.i / c(bits[["sumi.i"]]) -
                 Deriv0Matrix.i / Denom)^2 / Delta) else 0)

    if (linfla > 0)
      ned2l.dpstr.i.shape.p <-
        rowSums(d1B.pi * (1 + onembothprobs *
                         (d0A.io - d0B.pi) / Delta))

    if (linfla > 1)
      ned2l.dshape.i2 <- pstr.i *
        rowSums(pstr.i * (d1A.io^2) / Delta - d2A.io)

    
    if (linfla > 1)
      ned2l.dpstr.i.shape.i <-
        rowSums(d1A.io * (pstr.i * (d0A.io - d0B.pi) / Delta - 1))


    if (linfla > 1)
      ned2l.dshape.p.shape.i <- pstr.i * onembothprobs *
                                  rowSums(d1A.io * d1B.pi / Delta)


    if (lalter > 0 && linfla > 0)
      ned2l.dpobs.a.pstr.i <- onemrsDelta / onembothprobs^2 -
        rowSums(d0B.pi * (d0A.io - d0B.pi) / Delta)


    if (linfla > 1)
      ned2l.dpobs.a.shape.i <-
        rowSums(-pstr.i * d1A.io * d0B.pi / Delta)

    
    if (linfla > 0)
      ned2l.dpobs.a.shape.p <-
        rowSums(d1B.pi * (1 - onembothprobs * d0B.pi / Delta))
   

    wz <- matrix(0, n, if (lalter && !linfla) M1 else M1*(M1+1)/2)

    wz[, iam(1, 1, M)] <- ned2l.dshape.p2 * dshape.p.deta^2

    if (linfla > 0)
      wz[, iam(1, 2, M)] <-
        ned2l.dpobs.a.shape.p * dpobs.a.deta * dshape.p.deta
    if (lalter > 0)
      wz[, iam(2, 2, M)] <- wz22
    if (lalter > 1)
      wz[, iam(3, 3, M)] <- ned2l.dshape.a2 * dshape.a.deta^2

    if (linfla > 0) {
      wz[, iam(Offset.i + 1, Offset.i + 1, M)] <-
        ned2l.dpstr.i2 * dpstr.i.deta^2
      wz[, iam(Offset.i + 1, 1, M)] <-
        ned2l.dpstr.i.shape.p * dpstr.i.deta * dshape.p.deta
      if (lalter > 0)
        wz[, iam(Offset.i + 1, 2, M)] <-
          ned2l.dpobs.a.pstr.i * dpobs.a.deta * dpstr.i.deta
    }

    if (linfla > 1) {
      wz[, iam(Offset.i + 2, Offset.i + 2, M)] <-
        ned2l.dshape.i2 * dshape.i.deta^2
      wz[, iam(Offset.i + 2, 1, M)] <-
        ned2l.dshape.p.shape.i * dshape.p.deta * dshape.i.deta
      wz[, iam(Offset.i + 2, 2, M)] <-
        ned2l.dpobs.a.shape.i * dpobs.a.deta * dshape.i.deta
      wz[, iam(Offset.i + 2, Offset.i + 1, M)] <-
        ned2l.dpstr.i.shape.i * dpstr.i.deta * dshape.i.deta
    }

    c(w) * wz
  }), list( .lpobs.a = lpobs.a, .lpstr.i = lpstr.i))))
}  # gaitzeta.mix

















 dgaitpois <-
  function(x, lambda.p,
           alter.mix = NULL,
           alter.mlm = NULL,
           inflate.mix = NULL,
           inflate.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix.a = 0,  # vector
           pobs.mlm.a = 0,  # matrix
           pstr.mix.i = 0,  # vector
           pstr.mlm.i = 0,  # matrix
           lambda.a = lambda.p, lambda.i = lambda.p,
           byrow.arg = FALSE,  # Applies to 'pobs.mlm.a' & 'pstr.mlm.i'
           deflation = FALSE,  # Single logical
           log.arg = FALSE) {
  lowsup <- 0
  gait.errorcheck(alter.mix, alter.mlm, inflate.mix, inflate.mlm,
                  truncate, max.support)
  lalter.mix <- length(alter.mix)
  lalter.mlm <- length(alter.mlm)
  linfla.mix <- length(inflate.mix)
  linfla.mlm <- length(inflate.mlm)
  ltrunc     <- length(truncate)
  if (lalter.mix + lalter.mlm + linfla.mix + linfla.mlm + ltrunc == 0 &&
      is.infinite(max.support))
    return(dpois(x, lambda.p, log = log.arg))


  if (lalter.mix == 0) pobs.mix.a <- 0
  if (lalter.mlm == 0) pobs.mlm.a <- 0
  if (linfla.mix == 0) pstr.mix.i <- 0
  if (linfla.mlm == 0) pstr.mlm.i <- 0
 
  if (any(pobs.mix.a < 0 | 1 <= pobs.mix.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix.a'")
  if (any(pobs.mlm.a < 0 | 1 <= pobs.mlm.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm.a'")
  if (any(pstr.mix.i < 0 | 1 <= pstr.mix.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix.i'")
  if (any(1 <= pstr.mlm.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm.i'")
  if (!deflation && any(pstr.mlm.i < 0, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm.i'")

  LLL <- max(length(lambda.p),   length(x),
             length(pobs.mix.a), length(lambda.a),
             length(pstr.mix.i), length(lambda.i))  # added
  if (length(x)          != LLL) x          <- rep_len(x,          LLL)
  if (length(lambda.p)   != LLL) lambda.p   <- rep_len(lambda.p,   LLL)
  if (length(lambda.a)   != LLL) lambda.a   <- rep_len(lambda.a,   LLL)
  if (length(lambda.i)   != LLL) lambda.i   <- rep_len(lambda.i,   LLL)
  if (length(pobs.mix.a) != LLL) pobs.mix.a <- rep_len(pobs.mix.a, LLL)
  if (length(pstr.mix.i) != LLL) pstr.mix.i <- rep_len(pstr.mix.i, LLL)



  sumt <- 0  # Initialization to 0 important
  if (ltrunc)
    for (tval in truncate)
      sumt <- sumt + dpois(tval, lambda.p)  # Need tval <= max.support
  vecTF.t <- is.finite(x) & ((x %in% truncate) | (max.support < x))
  lhs.prob <- ppois(max.support, lambda.p)  # Usually 1
  denom.t <- lhs.prob - sumt  # No sumt on RHS

    pmf0 <- ifelse(vecTF.t, 0, dpois(x, lambda.p) / denom.t)  # dgtpois


  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (lalter.mlm) {
    pobs.mlm.a <-  matrix(pobs.mlm.a, LLL, lalter.mlm,
                          byrow = byrow.arg)
    sum.a <- .rowSums(pobs.mlm.a, LLL, lalter.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm.a'")  # zz

    for (aval in alter.mlm)
      suma <- suma + dpois(aval, lambda.p)  # Part i

    for (jay in seq(lalter.mlm)) {
      aval <- alter.mlm[jay]
      if (any(vecTF <- is.finite(x) & aval == x)) {
          pmf0[vecTF] <- pobs.mlm.a[vecTF, jay]
      }
      vecTF.a <- vecTF.a | vecTF  # Cumulative
    }  # jay
  }  # lalter.mlm



  pmf2.a <- pmf2.i <- 0
  if (lalter.mix) {
    allx.a <- lowsup:max(alter.mix)
    pmf2.a <- dgaitpois(x, lambda.a,  # Outer distribution---mlm type
                        truncate = setdiff(allx.a, alter.mix),
                        max.support = max(alter.mix))
    for (aval in alter.mix) {
      suma <- suma + dpois(aval, lambda.p)  # Part ii added; cumulative
      vecTF <- is.finite(x) & aval == x
      pmf0[vecTF] <- 0  # added; the true values are assigned below
      vecTF.a <- vecTF.a | vecTF  # Cumulative; added
    }
  }

  if (linfla.mix) {
    allx.i <- if (length(inflate.mix)) lowsup:max(inflate.mix) else NULL
    pmf2.i <- dgaitpois(x, lambda.i,  # Outer distribution---mlm type
                        truncate = setdiff(allx.i, inflate.mix),
                        max.support = max(inflate.mix))
  }




  sum.i <- 0
  if (linfla.mlm) {
    pstr.mlm.i <-  matrix(pstr.mlm.i, LLL, linfla.mlm,
                          byrow = byrow.arg)
    sum.i <- .rowSums(pstr.mlm.i, LLL, linfla.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm.i'")
  }  # linfla.mlm

  skip <- vecTF.t | vecTF.a  # Leave these values alone
  tmp6 <- 1 - sum.a - sum.i - pobs.mix.a - pstr.mix.i
  if (linfla.mlm) {
    if (deflation) {
      tmp0 <- lhs.prob - suma - sumt
      for (jay in 1:linfla.mlm) {
        vecTF <- is.finite(x) & inflate.mlm[jay] == x
        pmf.i <- dpois(inflate.mlm[jay], lambda.p[vecTF])
        if (any(pstr.mlm.i[vecTF, jay] <
                -(tmp6[vecTF] + pstr.mlm.i[vecTF, jay]) * pmf.i / (
                  tmp0[vecTF] - pmf.i), na.rm = TRUE)) {
          warning("too much deflation in argument 'pstr.mlm.i'. ",
                  "Returning NA")
          tmp6[vecTF] <- NA
        }
      }  # for
    } else {
      if (any(tmp6[!skip] < 0, na.rm = TRUE)) {
        warning("the sum of variables 'sum.a', 'sum.i', 'pobs.mix.a' ",
                "and 'pstr.mix.i' exceeds unity. Returning NA")
        tmp6[!skip & tmp6 < 0] <- NA
      }
    }  # deflation
  }  # linfla.mlm


  pmf0[!skip] <- (tmp6 *
    dpois(x, lambda.p) / (lhs.prob - suma - sumt))[!skip]  # added


  if (linfla.mlm) {
    for (jay in seq(linfla.mlm)) {
      ival <- inflate.mlm[jay]
      if (any(vecTF <- is.finite(x) & ival == x)) {
          pmf0[vecTF] <- pmf0[vecTF] + pstr.mlm.i[vecTF, jay]
      }
    }  # jay
  }  # linfla.mlm


  pmf0 <- pmf0 + pobs.mix.a * pmf2.a + pstr.mix.i * pmf2.i




if (FALSE) {
  allx.a <- if (length(alter.mix))   lowsup:max(alter.mix  ) else NULL
  allx.i <- if (length(inflate.mix)) lowsup:max(inflate.mix) else NULL
  use.ms <- if (is.finite(max.support)) {
    whatsleft <- setdiff(lowsup:max.support, alter.mix)
    if (length(whatsleft) > 0)
      max(whatsleft) else max.support
  } else {
    max.support
  }
  pmf1 <- dgaitpois(x, lambda.p,  # Inner distribution
                    truncate = c(alter.mix[alter.mix < use.ms],
                                 alter.mlm[alter.mlm < use.ms],  # added
                                 truncate),  # No ivec
                    max.support = use.ms)
  if (linfla.mix)
    pmf2.i <- dgaitpois(x, lambda.i,  # Outer distribution---mlm type
                        truncate = setdiff(allx.i, inflate.mix),
                        max.support = max(inflate.mix))
  pmf0.check <- pobs.mix.a * pmf2.a +
                pstr.mix.i * pmf2.i +  # mixprob * pmf2ai +
                (1 - sum.a - sum.i - pobs.mix.a - pstr.mix.i) * pmf1
}  # TRUE or FALSE

  if (log.arg) log(pmf0) else pmf0
}  # dgaitpois







 gait.errorcheck <-
  function(alter.mix = NULL, alter.mlm = NULL,
           inflate.mix = NULL, inflate.mlm = NULL,
           truncate = NULL,
           max.support = Inf,
           min.support = 0) {
  lalter.mix <- length(alter.mix)
  lalter.mlm <- length(alter.mlm)
  linfla.mix <- length(inflate.mix)
  linfla.mlm <- length(inflate.mlm)
  ltrunc <- length(truncate)

  if (!is.numeric(max.support) || is.na(max.support) ||
      length(max.support) != 1 || max.support < min.support ||
      round(max.support) != max.support ||
      (length(truncate) && (
          min(truncate, na.rm = TRUE) < min.support ||
          max.support <= max(truncate, na.rm = TRUE))))
    stop("bad input for argument 'max.support' and/or ",
         "'truncate'")

  allargs <- c(alter.mix, alter.mlm, inflate.mix, inflate.mlm)
  allargs <- c(allargs, truncate)  # No NA, NaN, -Inf or Inf allowed
  if (lalter.mix + lalter.mlm + linfla.mix + linfla.mlm)
    if (!is.Numeric(allargs, integer.valued = TRUE) ||
        any(allargs < min.support) ||
        any(max.support < allargs))
      stop("bad input for arguments 'alter.mix', 'alter.mlm',",
           " 'inflate.mix' and/or 'inflate.mlm'")
  if (length(unique(allargs)) < lalter.mix + lalter.mlm +
                                linfla.mix + linfla.mlm + ltrunc)
      stop("duplicate values found in arguments 'alter.mix', ",
           "'alter.mlm', 'inflate.mix', 'inflate.mlm' and 'truncate'")
}  # gait.errorcheck






 pgaitpois <-
  function(q, lambda.p,
           alter.mix = NULL,
           alter.mlm = NULL,
           inflate.mix = NULL,
           inflate.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix.a = 0,
           pobs.mlm.a = 0,
           pstr.mix.i = 0,
           pstr.mlm.i = 0,
           lambda.a = lambda.p, lambda.i = lambda.p,
           byrow.arg = FALSE) {
  lowsup <- 0
  gait.errorcheck(alter.mix, alter.mlm, inflate.mix, inflate.mlm,
                  truncate, max.support)
  lalter.mix <- length(alter.mix)
  lalter.mlm <- length(alter.mlm)
  linfla.mix <- length(inflate.mix)
  linfla.mlm <- length(inflate.mlm)
  ltrunc     <- length(truncate)
  if (lalter.mix + lalter.mlm + linfla.mix + linfla.mlm + ltrunc == 0 &&
      is.infinite(max.support))
    return(ppois(q, lambda.p))  # lower.tail, log.p


  if (lalter.mix == 0) pobs.mix.a <- 0
  if (lalter.mlm == 0) pobs.mlm.a <- 0
  if (linfla.mix == 0) pstr.mix.i <- 0
  if (linfla.mlm == 0) pstr.mlm.i <- 0

  if (any(pobs.mix.a < 0 | 1 <= pobs.mix.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix.a'")
  if (any(pobs.mlm.a < 0 | 1 <= pobs.mlm.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm.a'")
  if (any(pstr.mix.i < 0 | 1 <= pstr.mix.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix.i'")
  if (any(pstr.mlm.i < 0 | 1 <= pstr.mlm.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm.i'")

  LLL <- max(length(lambda.p),   length(q),
             length(pobs.mix.a), length(lambda.a),
             length(pstr.mix.i), length(lambda.i))  # added
  offset.a <- offset.i <- Offset.a <- Offset.i <- numeric(LLL)
  if (length(q)          != LLL) q          <- rep_len(q,          LLL)
  if (length(lambda.p)   != LLL) lambda.p   <- rep_len(lambda.p,   LLL)
  if (length(lambda.a)   != LLL) lambda.a   <- rep_len(lambda.a,   LLL)
  if (length(lambda.i)   != LLL) lambda.i   <- rep_len(lambda.i,   LLL)
  if (length(pobs.mix.a) != LLL) pobs.mix.a <- rep_len(pobs.mix.a, LLL)
  if (length(pstr.mix.i) != LLL) pstr.mix.i <- rep_len(pstr.mix.i, LLL)


  sumt <- 0
  fudge.t <- numeric(LLL)
  lhs.prob <- ppois(max.support, lambda.p)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      pmf.p <- dpois(tval, lambda.p)
      sumt <- sumt + pmf.p
      if (any(vecTF <- is.finite(q) & tval <= q))
        fudge.t[vecTF] <- fudge.t[vecTF] + pmf.p[vecTF]
    }
  }  # ltrunc

  sum.a <- suma <- 0  # numeric(LLL)
  fudge.a <- numeric(LLL)
  if (lalter.mlm) {
    pobs.mlm.a <- matrix(pobs.mlm.a, LLL, lalter.mlm, byrow = byrow.arg)
    sum.a <- .rowSums(pobs.mlm.a, LLL, lalter.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm.a'")

    for (jay in seq(lalter.mlm)) {
      aval <- alter.mlm[jay]
      pmf.p <- dpois(aval, lambda.p)
      suma <- suma + pmf.p  # cumulative; part i
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        offset.a[vecTF] <- offset.a[vecTF] + pobs.mlm.a[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # lalter.mlm

  sum.i <- 0
  if (linfla.mlm) {
    pstr.mlm.i <- matrix(pstr.mlm.i, LLL, linfla.mlm, byrow = byrow.arg)
    sum.i <- .rowSums(pstr.mlm.i, LLL, linfla.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm.i'")

    for (jay in seq(linfla.mlm)) {
      ival <- inflate.mlm[jay]
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        offset.i[vecTF] <- offset.i[vecTF] + pstr.mlm.i[vecTF, jay]
      }
    }  # jay
  }  # linfla.mlm



  use.pobs.mix.a <- 0
  if (lalter.mix) {
    use.pobs.mix.a <- matrix(0, LLL, lalter.mix)
    for (jay in seq(lalter.mix)) {
      aval <- alter.mix[jay]
      pmf.a <- dpois(aval, lambda.a)
      pmf.p <- dpois(aval, lambda.p)
      use.pobs.mix.a[, jay] <- pmf.a
      suma <- suma + pmf.p  # cumulative; part ii
    }
    use.pobs.mix.a <- pobs.mix.a *
                      use.pobs.mix.a / rowSums(use.pobs.mix.a)

    for (jay in seq(lalter.mix)) {
      aval <- alter.mix[jay]
      pmf.p <- dpois(aval, lambda.p)
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        Offset.a[vecTF] <- Offset.a[vecTF] + use.pobs.mix.a[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # lalter.mix

  use.pstr.mix.i <- 0
  if (linfla.mix) {
    use.pstr.mix.i <- matrix(0, LLL, linfla.mix)
    for (jay in seq(linfla.mix)) {
      ival <- inflate.mix[jay]
      use.pstr.mix.i[, jay] <- dpois(ival, lambda.i)
    }
    use.pstr.mix.i <- pstr.mix.i *
                      use.pstr.mix.i / rowSums(use.pstr.mix.i)

    for (jay in seq(linfla.mix)) {
      ival <- inflate.mix[jay]
      pmf.p <- dpois(ival, lambda.p)
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        Offset.i[vecTF] <- Offset.i[vecTF] + use.pstr.mix.i[vecTF, jay]
      }
    }  # jay
  }  # linfla.mix

  numer1 <- 1 - sum.i - sum.a - pstr.mix.i - pobs.mix.a
  denom1 <- lhs.prob - sumt - suma
  ans <- numer1 * (ppois(q, lambda.p) - fudge.t - fudge.a) / denom1 +
         offset.i + offset.a + Offset.i + Offset.a
  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  ans
}  # pgaitpois






 qgaitpois <-
  function(p, lambda.p,
           alter.mix = NULL,
           alter.mlm = NULL,
           inflate.mix = NULL,
           inflate.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix.a = 0,
           pobs.mlm.a = 0,
           pstr.mix.i = 0,
           pstr.mlm.i = 0,
           lambda.a = lambda.p, lambda.i = lambda.p,
           byrow.arg = FALSE) {
  lowsup <- 0
  gait.errorcheck(alter.mix, alter.mlm, inflate.mix, inflate.mlm,
                  truncate, max.support)
  lalter.mix <- length(alter.mix)
  lalter.mlm <- length(alter.mlm)
  linfla.mix <- length(inflate.mix)
  linfla.mlm <- length(inflate.mlm)
  ltrunc     <- length(truncate)
  if (lalter.mix + lalter.mlm + linfla.mix + linfla.mlm + ltrunc == 0 &&
      is.infinite(max.support))
    return(qpois(p, lambda.p))  # lower.tail = TRUE, log.p = FALSE


  if (lalter.mix == 0) pobs.mix.a <- 0
  if (lalter.mlm == 0) pobs.mlm.a <- 0
  if (linfla.mix == 0) pstr.mix.i <- 0
  if (linfla.mlm == 0) pstr.mlm.i <- 0
 
  if (any(pobs.mix.a < 0 | 1 <= pobs.mix.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix.a'")
  if (any(pobs.mlm.a < 0 | 1 <= pobs.mlm.a, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm.a'")
  if (any(pstr.mix.i < 0 | 1 <= pstr.mix.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix.i'")
  if (any(pstr.mlm.i < 0 | 1 <= pstr.mlm.i, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm.i'")


  LLL <- max(length(lambda.p),   length(p),
             length(pobs.mix.a), length(lambda.a),
             length(pstr.mix.i), length(lambda.i))  # added
  if (length(p)          != LLL) p          <- rep_len(p,          LLL)
  if (length(lambda.p)   != LLL) lambda.p   <- rep_len(lambda.p,   LLL)
  if (length(lambda.a)   != LLL) lambda.a   <- rep_len(lambda.a,   LLL)
  if (length(lambda.i)   != LLL) lambda.i   <- rep_len(lambda.i,   LLL)
  if (length(pobs.mix.a) != LLL) pobs.mix.a <- rep_len(pobs.mix.a, LLL)
  if (length(pstr.mix.i) != LLL) pstr.mix.i <- rep_len(pstr.mix.i, LLL)

  pobs.mlm.a <- matrix(pobs.mlm.a, LLL, max(lalter.mlm, 1),
                       byrow = byrow.arg)
  pstr.mlm.i <- matrix(pstr.mlm.i, LLL, max(linfla.mlm, 1),
                       byrow = byrow.arg)

  min.support <- lowsup  # Usual case; same as lowsup
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else min.support
  ans <- p + lambda.p

  bad0 <- !is.finite(lambda.p) | lambda.p <= 0
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  Lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- Lo  # True at lhs
  Hi <- if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else 2 * Lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgaitpois(Hi, lambda.p,
                   alter.mix = alter.mix, alter.mlm = alter.mlm,
                   inflate.mix = inflate.mix, inflate.mlm = inflate.mlm,
                   truncate = truncate, max.support = max.support,
                   pstr.mix.i = pstr.mix.i, pobs.mix.a = pobs.mix.a,
                   pstr.mlm.i = pstr.mlm.i, pobs.mlm.a = pobs.mlm.a,
                   lambda.a = lambda.a, lambda.i = lambda.i,
                   byrow.arg = FALSE)

  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    Lo[!done] <- Hi[!done]
    Hi[!done] <- 2 * Hi[!done] + 10.5  # Bug fixed
    Hi <- pmin(max.support + 0.5, Hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitpois(Hi[!done], lambda.p[!done],
                       alter.mix = alter.mix, alter.mlm = alter.mlm,
                       inflate.mix = inflate.mix,
                       inflate.mlm = inflate.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix.a = pobs.mix.a[!done],
                       pstr.mix.i = pstr.mix.i[!done],
                       pobs.mlm.a = pobs.mlm.a[!done, , drop = FALSE],
                       pstr.mlm.i = pstr.mlm.i[!done, , drop = FALSE],
                       lambda.a = lambda.a[!done],
                       lambda.i = lambda.i[!done],
                       byrow.arg = FALSE))
    iter <- iter + 1
  }

      foo <- function(q, lambda.p,
                      alter.mix = NULL, alter.mlm = NULL,
                      inflate.mix = NULL, inflate.mlm = NULL,
                      truncate = NULL, max.support = Inf,
                      pobs.mix.a = 0, pstr.mix.i = 0,
                      pobs.mlm.a = 0, pstr.mlm.i = 0,
                      lambda.a = lambda.p, lambda.i = lambda.p,
                      byrow.arg = FALSE, p)
      pgaitpois(q, lambda.p = lambda.p,
                       alter.mix = alter.mix, alter.mlm = alter.mlm,
                       inflate.mix = inflate.mix,
                       inflate.mlm = inflate.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix.a = pobs.mix.a,
                       pstr.mix.i = pstr.mix.i,
                       pobs.mlm.a = pobs.mlm.a,
                       pstr.mlm.i = pstr.mlm.i,
                       lambda.a = lambda.a, lambda.i = lambda.i,
                       byrow.arg = FALSE) - p

      lhs <- dont.iterate |
        p <= dgaitpois(min.support.use, lambda.p = lambda.p,
                       alter.mix = alter.mix, alter.mlm = alter.mlm,
                       inflate.mix = inflate.mix,
                       inflate.mlm = inflate.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix.a = pobs.mix.a,
                       pstr.mix.i = pstr.mix.i,
                       pobs.mlm.a = pobs.mlm.a,
                       pstr.mlm.i = pstr.mlm.i,
                       lambda.a = lambda.a, lambda.i = lambda.i,
                       byrow.arg = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, Lo[!lhs], Hi[!lhs], tol = 1/16,
                      lambda.p = lambda.p[!lhs],
                      alter.mix = alter.mix, alter.mlm = alter.mlm,
                      inflate.mix = inflate.mix,
                      inflate.mlm = inflate.mlm,
                      truncate = truncate, max.support = max.support,
                      pstr.mix.i = pstr.mix.i[!lhs],
                      pstr.mlm.i = pstr.mlm.i[!lhs, , drop = FALSE],
                      pobs.mix.a = pobs.mix.a[!lhs],
                      pobs.mlm.a = pobs.mlm.a[!lhs, , drop = FALSE],
                      lambda.a = lambda.a[!lhs],
                      lambda.i = lambda.i[!lhs],
                      byrow.arg = FALSE,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitpois(faa, lambda.p[!lhs],
                       alter.mix = alter.mix, alter.mlm = alter.mlm,
                       inflate.mix = inflate.mix,
                       inflate.mlm = inflate.mlm,
                       truncate = truncate, max.support = max.support,
                       pstr.mix.i = pstr.mix.i[!lhs],
                       pstr.mlm.i = pstr.mlm.i[!lhs, , drop = FALSE],
                       pobs.mix.a = pobs.mix.a[!lhs],
                       pobs.mlm.a = pobs.mlm.a[!lhs, , drop = FALSE],
                       lambda.a = lambda.a[!lhs],
                       lambda.i = lambda.i[!lhs],
                       byrow.arg = FALSE) < p[!lhs] &
             p[!lhs] <= pgaitpois(faa + 1, lambda.p[!lhs],
                       alter.mix = alter.mix, alter.mlm = alter.mlm,
                       inflate.mix = inflate.mix,
                       inflate.mlm = inflate.mlm,
                       truncate = truncate, max.support = max.support,
                       pstr.mix.i = pstr.mix.i[!lhs],
                       pstr.mlm.i = pstr.mlm.i[!lhs, , drop = FALSE],
                       pobs.mix.a = pobs.mix.a[!lhs],
                       pobs.mlm.a = pobs.mlm.a[!lhs, , drop = FALSE],
                       lambda.a = lambda.a[!lhs],
                       lambda.i = lambda.i[!lhs],
                       byrow.arg = FALSE),
             faa + 1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitpois(min.support.use, lambda.p,
                       alter.mix = alter.mix, alter.mlm = alter.mlm,
                       inflate.mix = inflate.mix,
                       inflate.mlm = inflate.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix.a = pobs.mix.a,
                       pstr.mix.i = pstr.mix.i,
                       pobs.mlm.a = pobs.mlm.a,
                       pstr.mlm.i = pstr.mlm.i,
                       lambda.a = lambda.a, lambda.i = lambda.i,
                       byrow.arg = FALSE)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitpois





 rgaitpois <-
  function(n, lambda.p,
           alter.mix = NULL,
           alter.mlm = NULL,
           inflate.mix = NULL,
           inflate.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix.a = 0,  # vector
           pobs.mlm.a = 0,  # matrix
           pstr.mix.i = 0,  # vector
           pstr.mlm.i = 0,  # matrix
           lambda.a = lambda.p, lambda.i = lambda.p,
           byrow.arg = FALSE) {
    qgaitpois(runif(n), lambda.p,
              alter.mix = alter.mix,
              alter.mlm = alter.mlm,
              inflate.mix = inflate.mix,
              inflate.mlm = inflate.mlm,
              truncate = truncate, max.support = max.support,
              pobs.mix.a = pobs.mix.a,
              pobs.mlm.a = pobs.mlm.a,
              pstr.mix.i = pstr.mix.i,
              pstr.mlm.i = pstr.mlm.i,
              lambda.a = lambda.a, lambda.i = lambda.i,
              byrow.arg = byrow.arg)
}  # rgaitpois













specialsvglm <-
  function(object, ...) {
  infos <- object@family@infos()
  ans <- list(alter    = infos$alter,
              inflate  = infos$inflate,
              truncate = infos$truncate)
  if (is.numeric(tmp7 <- infos$max.support))
    ans <- c(ans, max.support = tmp7)
  ans
}  # specialsvglm



if (!isGeneric("altered"))
  setGeneric("altered", function(object, ...)
             standardGeneric("altered"),
             package = "VGAM")
setMethod("altered", "vglm",
          function(object, ...)
          specialsvglm(object, ...)$alter)


if (!isGeneric("inflated"))
  setGeneric("inflated", function(object, ...)
             standardGeneric("inflated"),
             package = "VGAM")
setMethod("inflated", "vglm",
          function(object, ...)
          specialsvglm(object, ...)$inflate)


if (!isGeneric("truncated"))
  setGeneric("truncated", function(object, ...)
             standardGeneric("truncated"),
             package = "VGAM")
setMethod("truncated", "vglm",
          function(object, ...) {
          ans <- specialsvglm(object, ...)
          if (any(names(ans) == "max.support"))
            ans[c("truncate", "max.support")] else
            ans[["truncate"]]
          })


if (!isGeneric("specials"))
  setGeneric("specials", function(object, ...)
             standardGeneric("specials"),
             package = "VGAM")
setMethod("specials", "vglm",
          function(object, ...)
          specialsvglm(object, ...))










 y.gait.check <-
  function(alter = NULL, inflate = NULL,
           truncate = NULL, max.support = Inf, y,
           min.support = 0) {
  lalter <- length(alter)
  linfla <- length(inflate)

  n <- length(y)
  css.a <- css.i <- skip.a <- skip.i <- NULL  # Default

  if (length(truncate) && any(y %in% truncate))
    stop("some response values == values in argument 'truncate'")
  if (max.support < max(y))
    stop("some response values are greater than the ",
         "'max.support' argument")


  if (lalter > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0 <- matrix(0, n, lalter)
    for (jay in seq(lalter))
      y0[, jay] <- as.numeric(y == alter[jay])
    skip.a <- matrix(as.logical(y0), n, lalter)  # dim lost
    if (any((css.a <- colSums(skip.a)) == 0))
      stop("some 'alter' argument values have no response values: ",
           paste(alter[css.a == 0], collapse = ", "))          
  }  # lalter
  if (linfla > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0 <- matrix(0, n, linfla)
    for (jay in seq(linfla))
      y0[, jay] <- as.numeric(y == inflate[jay])
    skip.i <- matrix(as.logical(y0), n, linfla)  # dim lost
      if (any((css.i <- colSums(skip.i)) == 0))
      stop("some 'inflate' argument values have no response values: ",
           paste(inflate[css.i == 0], collapse = ", "))          
  }  # linfla
  list(css.a = css.a, skip.a = skip.a,
       css.i = css.i, skip.i = skip.i)
}  # y.gait.check






 gaitlog.mix <-
  function(alter = NULL,
           inflate = NULL,  
           truncate = NULL,   max.support = Inf,
           zero = c("pobs.a", "pstr.i"),  # Pruned later if necessary
           eq.ap = FALSE,  # TRUE applies to the intercept
           eq.ip = FALSE,  # TRUE applies to the intercept
           lshape.p = "logitlink",
           lpobs.a = "logitlink",
           lshape.a = "logitlink",
           lpstr.i = "logitlink",
           lshape.i = "logitlink",
         type.fitted = c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
                         "prob.a", "prob.i", "prob.t", "lhs.prob"),
           imethod = 1,
           gshape.p = -expm1(-7 * ppoints(12)),  # ppoints(4) inadequate
           ishape.p = NULL,  # Higher is better than lower
           ishape.a = NULL, ishape.i = NULL,
           ipobs.a = NULL,  # 0.25, 
           ipstr.i = NULL,  # 0.25, 
           ishrinkage = 0.95,
           probs.y = 0.35) {
  lowsup <- 1
  semigait.errorcheck(alter, inflate, truncate, max.support,
                  min.support = lowsup)

  if (!is.logical(eq.ap) || length(eq.ap) != 1)
    stop("argument 'eq.ap' must be a single logical")
  if (!is.logical(eq.ip) || length(eq.ip) != 1)
    stop("argument 'eq.ip' must be a single logical")

  lshape.p <- as.list(substitute(lshape.p))
  eshape.p <- link2list(lshape.p)
  lshape.p <- attr(eshape.p, "function.name")
  lpobs.a <- as.list(substitute(lpobs.a))  # \omega
  epobs.a <- link2list(lpobs.a)
  lpobs.a <- attr(epobs.a, "function.name")
  lshape.a <- as.list(substitute(lshape.a))
  eshape.a <- link2list(lshape.a)
  lshape.a <- attr(eshape.a, "function.name")

  lpstr.i <- as.list(substitute(lpstr.i))  # \phi
  epstr.i <- link2list(lpstr.i)
  lpstr.i <- attr(epstr.i, "function.name")
  lshape.i <- as.list(substitute(lshape.i))
  eshape.i <- link2list(lshape.i)
  lshape.i <- attr(eshape.i, "function.name")

  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)

  if (is.vector(zero) && is.character(zero) && length(zero) == 2) {
    if (lalter && linfla == 0)
      zero <- setdiff(zero, "pstr.i")
    if (linfla && lalter == 0)
      zero <- setdiff(zero, "pobs.a")
  }

  if (lalter + linfla == 0) zero <- NULL
  lshape.p.save <- lshape.p
  if (lalter + linfla + ltrunc == 0 && is.infinite(max.support))
    return(eval(substitute(
           logff(lshape = .lshape.p.save , zero = .zero ),
           list( .lshape.p.save = lshape.p.save, .zero = zero))))


  if (lalter == 1 && eq.ap)
    warning("Less than 2 altered values and hence 1 shape, so ",
            "setting 'eq.ap = TRUE' is meaningless")
  if (linfla == 1 && eq.ip)
    warning("Less than 2 inflated values and hence 1 shape, so ",
            "setting 'eq.ip = TRUE' is meaningless")

  type.fitted <- match.arg(type.fitted,
    c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
      "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]

  tmp3.TF <- c(TRUE, lalter > 0, lalter > 1, linfla > 0, linfla > 1)
  tmp3 <- c(shape.p = lshape.p,  # Version 2
            pobs.a   = lpobs.a, shape.a = lshape.a,
            pstr.i   = lpstr.i, shape.i = lshape.i)  # [tmp3.TF]
      
  blurb1 <- "P"
  if (lalter) blurb1 <- "Generally-altered l"
  if (linfla) blurb1 <- "Generally-inflated l"
  if (ltrunc) blurb1 <- "Generally-truncated l"
  if ( lalter &&  linfla && !ltrunc)
    blurb1 <- "Generally-altered and -inflated l"
  if ( lalter && !linfla &&  ltrunc)
    blurb1 <- "Generally-altered and -truncated l"
  if (!lalter &&  linfla &&  ltrunc)
    blurb1 <- "Generally-inflated and -truncated l"
  if ( lalter &&  linfla &&  ltrunc)
    blurb1 <- "Generally-altered, -inflated and -truncated l"

  new("vglmff",
  blurb = c(blurb1, "ogarithmic regression\n",
            "(GAIT-Log(shape.p)-Log(shape.a)-Log(shape.i) ",
            "mixture generally)\n\n",
            "Links:    ",
            namesof("shape.p",  lshape.p, eshape.p, tag = FALSE),
            if (lalter > 0) c(", ",
            namesof("pobs.a",    lpobs.a,   epobs.a,   tag = FALSE)),
            if (lalter > 1) c(", ",
            namesof("shape.a",  lshape.a, eshape.a, tag = FALSE)),
            if (lalter && linfla) ", \n        ",
            if (linfla > 0) c(  if (lalter) "  " else ", ",
            namesof("pstr.i",    lpstr.i,   epstr.i,   tag = FALSE)),
            if (linfla > 1) c(", ",
            namesof("shape.i",  lshape.i, eshape.i, tag = FALSE)),
            "\n"),
  constraints = eval(substitute(expression({

    use.mat <- if ( ( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,1,0,1,  0,1,0,0,0,  0,0,0,1,0), 5, 3) else
    if (!( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,0,0,1,  0,1,0,0,0,  0,0,1,0,0,
               0,0,0,1,0), 5, 4) else
    if ( ( .eq.ap ) && !( .eq.ip ))
      matrix(c(1,0,1,0,0,  0,1,0,0,0,  0,0,0,1,0,
               0,0,0,0,1), 5, 4) else
      diag(5)

    use.mat <- use.mat[ .tmp3.TF , , drop = FALSE]
    M1 <- nrow(use.mat)
    for (jay in ncol(use.mat):1)
      if (all(use.mat[, jay] == 0))
        use.mat <- use.mat[, -jay, drop = FALSE]
    
    constraints <- cm.VGAM(use.mat, x = x,
                           bool = .eq.ap | .eq.ip ,
                           constraints = constraints,
                           apply.int = TRUE)  # FALSE

    if ( .lalter + .linfla > 0)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = M1)
  }), list( .zero = zero, .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
            .lalter = lalter, .linfla = linfla,
            .eq.ap = eq.ap, .eq.ip = eq.ip ))),

  infos = eval(substitute(function(...) {
    list(Q1 = 1,
         M1 = sum( .tmp3.TF ),
         link = c( .tmp3 )[ .tmp3.TF ] ,
         link1parameter = TRUE,
         mixture.links = FALSE,
         alter = as.vector( .alter ),
         inflate = as.vector( .inflate ),
         truncate = as.vector( .truncate ),
         max.support = as.vector( .max.support ), 
         Support  = c( .lowsup , Inf, 1),
         eq.ap = .eq.ap , eq.ip = .eq.ip ,
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE,
         parameters.names = names(c( .tmp3 )[ .tmp3.TF ]),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero, .tmp3 = tmp3, .tmp3.TF = tmp3.TF,
           .type.fitted = type.fitted, .lowsup = lowsup,
           .eq.ap = eq.ap, .eq.ip = eq.ip,
           .alter = alter, .inflate = inflate,
           .truncate = truncate, .max.support = max.support ))),

  initialize = eval(substitute(expression({
    truncate <- as.vector( .truncate )
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    M1 <- sum( .tmp3.TF )
    NOS <- NCOL(y)  # Only 1 currently
    M <- NOS * M1
    tmp3 <- ( .tmp3 )
    tmp3.TF <- ( .tmp3.TF )
    ntmp3 <- names(tmp3)

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    glist <- y.gait.check(alter, inflate, truncate, .max.support , y)
    css.a <- glist$css.a
    css.i <- glist$css.i
    extra$skip.a <- glist$skip.a
    extra$skip.i <- glist$skip.i

    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- ( .type.fitted )
    extra$colnames.y  <- colnames(y)
    extra$M1 <- M1

   predictors.names <-
    c(namesof(ntmp3[1], .lshape.p , earg = .eshape.p , tag = FALSE),
      namesof(ntmp3[2], .lpobs.a  , earg = .epobs.a  , tag = FALSE),
      namesof(ntmp3[3], .lshape.a , earg = .eshape.a , tag = FALSE),
      namesof(ntmp3[4], .lpstr.i  , earg = .epstr.i  , tag = FALSE),
      namesof(ntmp3[5], .lshape.i , earg = .eshape.i , tag = FALSE))[
      tmp3.TF]  # Relevant subset



    

    if (!length(etastart)) {
      shape.p.init <- if (length( .ishape.p )) .ishape.p else {
        logff.Loglikfun <- function(shapeval, y, x, w, extraargs) {
          sum(c(w) * dlog(x = y, shape = shapeval, log = TRUE))
        }
        shape.p.grid <- ( .gshape.p )
        grid.search(shape.p.grid, objfun = logff.Loglikfun, y = y, w = w)
      }
      shape.p.init <- rep(shape.p.init, length = n)

        pobs.a.init <-  # Unneeded? Low is best
        pstr.i.init <- 0.025 / (1 + lalter + linfla)
      if (lalter > 0) {
        pobs.a.init <-  rep(if (length( .ipobs.a )) .ipobs.a else
                            sum(css.a) / n, length = n)
      }
      if (lalter > 1)
        shape.a.init <- if (length( .ishape.a ))
          rep( .ishape.a , n) else shape.p.init

      if (linfla > 0)
        pstr.i.init <-  rep(if (length( .ipstr.i )) .ipstr.i else
                            0.5 * sum(css.i) / n, length = n)
      if (linfla > 1)
        shape.i.init <- if (length( .ishape.i ))
          rep( .ishape.i , n) else shape.p.init

      while (any((vecTF <- pobs.a.init + pstr.i.init > 0.95))) {
        pobs.a.init[vecTF] <- 0.875 * pobs.a.init[vecTF]
        pstr.i.init[vecTF] <- 0.875 * pstr.i.init[vecTF]
      }

      etastart <-
        cbind(theta2eta(shape.p.init, .lshape.p , earg = .eshape.p ),
              if (lalter <= 0) NULL else
              theta2eta(pobs.a.init,  .lpobs.a  , earg = .epobs.a  ),
              if (lalter <= 1) NULL else
              theta2eta(shape.a.init, .lshape.a , earg = .eshape.a ),
              if (linfla <= 0) NULL else
              theta2eta(pstr.i.init,  .lpstr.i  , earg = .epstr.i  ),
              if (linfla <= 1) NULL else
              theta2eta(shape.i.init, .lshape.i , earg = .eshape.i ))
    }
  }), list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
    .ishape.p = ishape.p, .ishape.i = ishape.i, .ipstr.i = ipstr.i,
                          .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                          .eshape.a = eshape.a, .epobs.a = epobs.a,
    .ishape.p = ishape.p, .ishape.a = ishape.a, .ipobs.a = ipobs.a,
    .gshape.p = gshape.p, 
    .ishrinkage = ishrinkage, .probs.y = probs.y,
    .alter = alter, .inflate = inflate,
    .truncate = truncate, .max.support = max.support,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .imethod = imethod, .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <-
     if (length(extra$type.fitted)) extra$type.fitted else {
       warning("cannot find 'type.fitted'. Returning the 'mean'.")
                     "mean"
     }
   
    type.fitted <-
      match.arg(type.fitted,
                c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
                  "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]
    eta <- as.matrix(eta)  # Needed when linfla == 0
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate) 
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    max.support <- ( .max.support )
    M1 <- sum( .tmp3.TF )
    if ((NOS <- NCOL(eta) / M1) != 1)
      stop("Currently NOS must be 1")

    iptr <- 1  # Unpack
    shape.p <- eta2theta(eta[, iptr], .lshape.p , earg = .eshape.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a  , earg = .epobs.a )
    shape.a <- if (lalter <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.a , earg = .eshape.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i  , earg = .epstr.i )
    shape.i <- if (linfla <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.i , earg = .eshape.i )
    Bits <- moments.log.gait(shape.p = shape.p, mlm = FALSE,
                             alter = alter, inflate = inflate,
                             truncate = truncate,
                             max.support = max.support,
                             pobs.a = pobs.a, shape.a = shape.a,
                             pstr.i = pstr.i, shape.i = shape.i)
    if (type.fitted == "Pobs.a") {
      if (lalter == 0) stop("no altered values!")
      proportion.mat <-
        dlog(matrix(alter,   NROW(eta), lalter, byrow = TRUE),
             matrix(shape.a, NROW(eta), lalter)) / c(Bits[["suma.a"]])
    }
    if (type.fitted == "Pstr.i") {
      if (linfla == 0) stop("no inflated values!")
      proportion.mat <-
        dlog(matrix(inflate, nrow(eta), linfla, byrow = TRUE),
             matrix(shape.i, nrow(eta), linfla)) / c(Bits[["sumi.i"]])
    }
  
    ans <- switch(type.fitted,
      "mean"     = Bits[["mean"]],
      "pobs.a"   = pobs.a,
      "pstr.i"   = pstr.i,
      "Pobs.a"   = pobs.a * proportion.mat,  # matrix
      "Pstr.i"   = c(pstr.i) * proportion.mat +
                   c(1 - pobs.a - pstr.i) *
           dlog(matrix(inflate, nrow(eta), linfla, byrow = TRUE),
                matrix(shape.p, nrow(eta), linfla)) / c(
           Bits[["lhs.prob"]] - Bits[["suma.p"]] - Bits[["sumt.p"]]),
      "prob.a"   = Bits[["suma.p"]],
      "prob.i"   = Bits[["sumi.p"]],
      "prob.t"   = Bits[["sumt.p"]],
      "lhs.prob" = Bits[["lhs.prob"]])
   ynames.Pobs.a <- as.character(alter)    # Works with NULLs
   ynames.Pstr.i <- as.character(inflate)  # Works with NULLs
   label.cols.y(ans,
        colnames.y = switch(type.fitted,
                            "Pobs.a" = ynames.Pobs.a,
                            "Pstr.i" = ynames.Pstr.i,
                            extra$colnames.y),
        NOS = NOS)
  }, list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                          .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                          .eshape.a = eshape.a, .epobs.a = epobs.a,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .truncate = truncate, .max.support = max.support, 
    .alter = alter, .inflate = inflate ))),

  last = eval(substitute(expression({
    misc$link  <- c( .tmp3 )[ .tmp3.TF ]
    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    misc$earg[[(iptr <- 1)]] <- .eshape.p  #
    if (lalter > 0)
      misc$earg[[(iptr <- iptr + 1)]] <- .epobs.a  #
    if (lalter > 1)
      misc$earg[[(iptr <- iptr + 1)]] <- .eshape.a  #
    if (linfla > 0)
      misc$earg[[(iptr <- iptr + 1)]] <- .epstr.i  #
    if (linfla > 1)
      misc$earg[[(iptr <- iptr + 1)]] <- .eshape.i  #
  }), list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                            .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                            .eshape.a = eshape.a, .epobs.a = epobs.a,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .alter = alter, .inflate = inflate ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    lalter <- length(as.vector( .alter ))
    linfla <- length(as.vector( .inflate ))
    if (lalter + linfla == 0)
      eta <- as.matrix(eta)
    iptr <- 1  # Unpack
    shape.p <- eta2theta(eta[, iptr], .lshape.p , earg = .eshape.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a  , earg = .epobs.a )
    shape.a <- if (lalter <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.a , earg = .eshape.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i  , earg = .epstr.i )
    shape.i <- if (linfla <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.i , earg = .eshape.i )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitlog(y, shape.p, shape.a = shape.a, shape.i = shape.i,
                 pobs.mix.a = pobs.a, pstr.mix.i = pstr.i,
                 truncate = .truncate , max.support = .max.support ,
                 alter.mix = .alter , inflate.mix = .inflate , log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                          .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                          .eshape.a = eshape.a, .epobs.a = epobs.a,
    .truncate = truncate, .max.support = max.support, 
    .alter = alter, .inflate = inflate ))),
  vfamily = c("gaitlog.mix"),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    lalter <- length(as.vector( .alter ))
    linfla <- length(as.vector( .inflate ))
    if (lalter + linfla == 0)
      eta <- as.matrix(eta)
    iptr <- 1  # Unpack
    shape.p <- eta2theta(eta[, iptr], .lshape.p , earg = .eshape.p )
    pobs.a   <- if (lalter <= 0) rep(0.25, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a , earg = .epobs.a )
    shape.a <- if (lalter <= 1) shape.p else
    eta2theta(eta[, (iptr <- iptr + 1)], .lshape.a  , earg = .eshape.a )
    pstr.i   <- if (linfla <= 0) rep(0.25, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i , earg = .epstr.i )
    shape.i <- if (linfla <= 1) shape.p else
    eta2theta(eta[, (iptr <- iptr + 1)], .lshape.i  , earg = .eshape.i )

    okay1 <- all(is.finite(shape.p)) && all(0 < shape.p) &&
                                        all(shape.p < 1) &&
             all(is.finite(shape.a)) && all(0 < shape.a) &&
                                        all(shape.a < 1) &&
             all(is.finite(shape.i)) && all(0 < shape.i) &&
                                        all(shape.i < 1) &&
             all(is.finite(pobs.a)) && all(0 < pobs.a & pobs.a < 1) &&
             all(is.finite(pstr.i)) && all(0 < pstr.i & pstr.i < 1) &&
             all(pobs.a + pstr.i < 1)
    okay1
  }, list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                          .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                          .eshape.a = eshape.a, .epobs.a = epobs.a,
    .truncate = truncate, .max.support = max.support, 
    .alter = alter, .inflate = inflate ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    lalter <- length(as.vector( .alter ))
    linfla <- length(as.vector( .inflate ))
    iptr <- 1  # Unpack
    shape.p <- eta2theta(eta[, iptr], .lshape.p , earg = .eshape.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a  , earg = .epobs.a )
    shape.a <- if (lalter <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.a , earg = .eshape.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i  , earg = .epstr.i )
    shape.i <- if (linfla <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.i , earg = .eshape.i )
    rgaitlog(nsim * length(shape.p), shape.p = shape.p,
             shape.a = shape.a, shape.i = shape.i,
             pobs.mix.a = pobs.a, pstr.mix.i = pstr.i,
             alter.mix = .alter , inflate.mix = .inflate ,
             truncate = .truncate , max.support = .max.support )
  }, list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                          .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                          .eshape.a = eshape.a, .epobs.a = epobs.a,
    .truncate = truncate, .max.support = max.support,
    .alter = alter, .inflate = inflate ))),
  deriv = eval(substitute(expression({
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    M1 <- sum( .tmp3.TF )
    NOS <- NCOL(eta) / M1  # extra$NOS
    max.support <- as.vector( .max.support )

    ans2 <- matrix(0, n, M1)
    Offset.i <- if (lalter == 1) 2 else if (lalter > 1) 3 else 1


    is.altered <- if (lalter)
      rowSums(extra$skip.a) > 0 else rep(FALSE, n)
    is.inflated <- if (linfla)
      rowSums(extra$skip.i) > 0 else rep(FALSE, n)


    if (lalter + linfla == 0)
      eta <- as.matrix(eta)
    iptr <- 1  # Unpack
    shape.p <- eta2theta(eta[, iptr], .lshape.p , earg = .eshape.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a  , earg = .epobs.a )
    shape.a <- if (lalter <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.a , earg = .eshape.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i  , earg = .epstr.i )
    shape.i <- if (linfla <= 1) shape.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .lshape.i , earg = .eshape.i )

    bits <- moments.log.gait(shape.p, alter = alter,
                             inflate = inflate, truncate = truncate,
                             max.support = max.support,
                             pobs.a = pobs.a, pstr.i = pstr.i,
                             shape.a = shape.a, shape.i = shape.i,
                             mlm = FALSE)
    Denom <- c(bits[["lhs.prob"]] - bits[["suma.p"]] - bits[["sumt.p"]])

    pmf.deriv1 <- function(y, shape) {
      A8 <- -1 / log1p(-shape)
      deriv0 <- A8 * (shape^y) / y
      A8 * (shape^(y-1) - deriv0 / (1 - shape))
    }
    pmf.deriv2 <- function(y, shape) {
      A8 <- -1 / log1p(-shape)
      A8prime <- -(A8^2) / (1 - shape)
      deriv0 <- A8 * (shape^y) / y
      deriv1 <- A8 * (shape^(y-1) - deriv0 / (1 - shape))
      A8prime * (shape^(y-1) - deriv0 / (1 - shape)) +
      A8 * ((y - 1) * shape^(y - 2) - deriv0 / (1 - shape)^2 -
            deriv1 / (1 - shape)) 
    }
    
    sumderiv1a.p <- sumderiv2a.p <- matrix(0, n, NOS)
    if (lalter > 0) {
      deriv0matrix.a <-
      deriv1matrix.a <- deriv2matrix.a <- matrix(0, n, lalter)
      sumderiv1a.a <- sumderiv2a.a <- matrix(0, n, NOS)
      for (jay in seq(lalter)) {
        aval <- alter[jay]
        sumderiv1a.p <- sumderiv1a.p + pmf.deriv1(aval, shape.p)
        sumderiv2a.p <- sumderiv2a.p + pmf.deriv2(aval, shape.p)
        sumderiv1a.a <- sumderiv1a.a + pmf.deriv1(aval, shape.a)
        sumderiv2a.a <- sumderiv2a.a + pmf.deriv2(aval, shape.a)
        pmf.a <- dlog(aval, shape.a)
        deriv0matrix.a[, jay] <- pmf.a
        deriv1matrix.a[, jay] <- pmf.deriv1(aval, shape.a)
        deriv2matrix.a[, jay] <- pmf.deriv2(aval, shape.a)
      }
    }  # lalter > 0





    if (linfla) {
      deriv0matrix.i <-  # wrt inflated distribution
      deriv1matrix.i <- deriv2matrix.i <- matrix(0, n, linfla)
      Deriv0Matrix.i <-  # wrt parent distribution
      Deriv1Matrix.i <- Deriv2Matrix.i <- matrix(0, n, linfla)
      if (linfla > 0)
        for (jay in seq(linfla)) {
          ival <- inflate[jay]
          pmf.i <- dlog(ival, shape.i)
          deriv0matrix.i[, jay] <- pmf.i
          deriv1matrix.i[, jay] <- pmf.deriv1(ival, shape.i)
          deriv2matrix.i[, jay] <- pmf.deriv2(ival, shape.i)
          pmf.p <- dlog(ival, shape.p)
          Deriv0Matrix.i[, jay] <- pmf.p
          Deriv1Matrix.i[, jay] <- pmf.deriv1(ival, shape.p)
          Deriv2Matrix.i[, jay] <- pmf.deriv2(ival, shape.p)  # /pmf.p
        }  # jay
    }  # linfla




    sumderiv1t.a <- sumderiv2t.a <-
    sumderiv1t.i <- sumderiv2t.i <-
    sumderiv1t.p <- sumderiv2t.p <- matrix(0, n, NOS)
    if (ltrunc)
      for (tval in truncate) {
        sumderiv1t.p <- sumderiv1t.p + pmf.deriv1(tval, shape.p)
        sumderiv2t.p <- sumderiv2t.p + pmf.deriv2(tval, shape.p)
        sumderiv1t.a <- sumderiv1t.a + pmf.deriv1(tval, shape.a)
        sumderiv2t.a <- sumderiv2t.a + pmf.deriv2(tval, shape.a)
        sumderiv1t.i <- sumderiv1t.i + pmf.deriv1(tval, shape.i)
        sumderiv2t.i <- sumderiv2t.i + pmf.deriv2(tval, shape.i)
      }

    onembothprobs <- 1 - pobs.a - pstr.i  # 1-Pr(altered)-Pr(inflated)
    A8.p <- -1 / log1p(-shape.p)
    A8.a <- -1 / log1p(-shape.a)
    A8.i <- -1 / log1p(-shape.i)
    if (is.finite(max.support)) {
    tmp1.p <- A8.p * (shape.p^max.support -
      (1 - plog(max.support, shape.p))) / (1 - shape.p)
    sumderiv1t.p <- sumderiv1t.p + tmp1.p
    sumderiv2t.p <- sumderiv2t.p + (A8.p / (1 - shape.p)) * (
      (shape.p^max.support) / (1 - shape.p) +
      max.support * shape.p^(max.support - 1) -
      (1 - plog(max.support, shape.p)) / (1 - shape.p) - 2 * tmp1.p)


    tmp1.a <- A8.a * (shape.a^max.support -
      (1 - plog(max.support, shape.a))) / (1 - shape.a)
    sumderiv1t.a <- sumderiv1t.a + tmp1.a
    sumderiv2t.a <- sumderiv2t.a + (A8.a / (1 - shape.a)) * (
      (shape.a^max.support) / (1 - shape.a) +
      max.support * shape.a^(max.support - 1) -
      (1 - plog(max.support, shape.a)) / (1 - shape.a) - 2 * tmp1.a)
        
    tmp1.i <- A8.i * (shape.i^max.support -
      (1 - plog(max.support, shape.i))) / (1 - shape.i)
    sumderiv1t.i <- sumderiv1t.i + tmp1.i
    sumderiv2t.i <- sumderiv2t.i + (A8.i / (1 - shape.i)) * (
      (shape.i^max.support) / (1 - shape.i) +
      max.support * shape.i^(max.support - 1) -
      (1 - plog(max.support, shape.i)) / (1 - shape.i) - 2 * tmp1.i)
    }  # is.finite(max.support)


    zero0n <- rep(0, n)
    dl.dshape.a <-
    dl.dshape.i <- zero0n  # Replace some elts below
    dl.dshape.p <- -A8.p / (1 - shape.p) + y / shape.p +
      (sumderiv1a.p + sumderiv1t.p) / Denom  # \notin A, I, T
    dl.dshape.p[is.altered] <- 0
    dl.dpobs.a <-  # Replace some elts below
    dl.dpstr.i <- -1 / onembothprobs  # \notin A, I, T
    Delta <- matrix(0, n, NOS)  # If linfla == 0

    if (linfla > 0) {
      d0A.io <- deriv0matrix.i / c(bits[["sumi.i"]])
      d0B.pi <- Deriv0Matrix.i / Denom
      Delta <- pstr.i * d0A.io + onembothprobs * d0B.pi


      d1A.io <- deriv1matrix.i / c(bits[["sumi.i"]]) -
                deriv0matrix.i * rowSums(deriv1matrix.i) / c(
                    bits[["sumi.i"]])^2
      d2A.io <- deriv2matrix.i / c(bits[["sumi.i"]]) -
                2 * deriv1matrix.i * rowSums(deriv1matrix.i) / c(
                    bits[["sumi.i"]])^2 -
                deriv0matrix.i * rowSums(deriv2matrix.i) / c(
                    bits[["sumi.i"]])^2 +
                2 * deriv0matrix.i * rowSums(deriv1matrix.i)^2 / c(
                    bits[["sumi.i"]])^3
      d1B.pi <- Deriv1Matrix.i / Denom +
        Deriv0Matrix.i * c(sumderiv1t.p + sumderiv1a.p) / Denom^2
      d2B.pi <- Deriv2Matrix.i / Denom +
        2 * Deriv1Matrix.i * c(sumderiv1t.p + sumderiv1a.p) / Denom^2 +
        Deriv0Matrix.i * c(sumderiv2t.p + sumderiv2a.p) / Denom^2 +
        2 * Deriv0Matrix.i * (c(sumderiv1t.p + sumderiv1a.p)^2) / Denom^3
    }  # linfla > 0
    
    if (lalter) {
      dl.dpobs.a[is.altered] <- 1 / pobs.a[is.altered]
      dl.dpstr.i[is.altered] <- 0

      d0A.al <- deriv0matrix.a / c(bits[["suma.a"]])
      d1A.al <- deriv1matrix.a / c(bits[["suma.a"]]) -
                deriv0matrix.a * rowSums(deriv1matrix.a) / c(
                    bits[["suma.a"]])^2
      d2A.al <- deriv2matrix.a / c(bits[["suma.a"]]) -
                2 * deriv1matrix.a * rowSums(deriv1matrix.a) / c(
                    bits[["suma.a"]])^2 -
                deriv0matrix.a * rowSums(deriv2matrix.a) / c(
                    bits[["suma.a"]])^2 +
                2 * deriv0matrix.a * rowSums(deriv1matrix.a)^2 / c(
                    bits[["suma.a"]])^3

      for (jay in seq(lalter)) {
        aval <- alter[jay]
        is.alte.j <- extra$skip.a[, jay]  # Logical vector
        tmp2 <- d1A.al[, jay] / d0A.al[, jay]
        dl.dshape.a[is.alte.j] <- tmp2[is.alte.j]
      }  # jay
    }  # lalter

    if (linfla > 0)
      for (jay in seq(linfla)) {
        ival <- inflate[jay]
        is.infl.j <- extra$skip.i[, jay]  # Logical vector
        tmp7 <- onembothprobs * d1B.pi[, jay] / Delta[, jay]
        dl.dshape.p[is.infl.j] <- tmp7[is.infl.j]
        tmp2 <- -d0B.pi[, jay] / Delta[, jay]
        dl.dpobs.a[is.infl.j] <- tmp2[is.infl.j]
        tmp8 <- (d0A.io[, jay] - d0B.pi[, jay]) / Delta[, jay]
        dl.dpstr.i[is.infl.j] <- tmp8[is.infl.j]
        if (linfla > 1) {
          tmp2 <- pstr.i * d1A.io[, jay] / Delta[, jay]
          dl.dshape.i[is.infl.j] <- tmp2[is.infl.j]
        }
      }  # jay


    dshape.p.deta <- dtheta.deta(shape.p, .lshape.p , .eshape.p )
    ans2[, 1] <- dl.dshape.p * dshape.p.deta
    dpobs.a.deta <- dtheta.deta(pobs.a, .lpobs.a , .epobs.a )
    if (lalter > 0) {
      ans2[, 2] <- dl.dpobs.a * dpobs.a.deta
    }
    if (lalter > 1) {
      dshape.a.deta <- dtheta.deta(shape.a, .lshape.a , .eshape.a )
      ans2[, 3] <- dl.dshape.a * dshape.a.deta
    }
    if (linfla > 0) {
      dpstr.i.deta <- dtheta.deta(pstr.i, .lpstr.i , .epstr.i )
      ans2[, Offset.i + 1] <- dl.dpstr.i * dpstr.i.deta
    }
    if (linfla > 1) {
      dshape.i.deta <- dtheta.deta(shape.i, .lshape.i , .eshape.i )
      ans2[, Offset.i + 2] <- dl.dshape.i * dshape.i.deta
    }
    c(w) * ans2
  }), list(
    .lshape.p = lshape.p, .lshape.i = lshape.i, .lpstr.i = lpstr.i,
    .eshape.p = eshape.p, .eshape.i = eshape.i, .epstr.i = epstr.i,
                            .lshape.a = lshape.a, .lpobs.a = lpobs.a,
                            .eshape.a = eshape.a, .epobs.a = epobs.a,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .truncate = truncate, .max.support = max.support, 
    .alter = alter , .inflate = inflate ))),

  weight = eval(substitute(expression({
    mean.true.p <- A8.p * shape.p / (1 - shape.p)
    cond.EY.p <- (mean.true.p - bits[["SumA.p"]] - bits[["SumI.p"]] -
                  bits[["SumT.p"]]) / c(
       bits[["lhs.prob"]] - bits[["suma.p"]] - bits[["sumi.p"]] -
       bits[["sumt.p"]])
    onemrsDelta <- onembothprobs *  # (onempstr.i - pobs.a) *
                   (1 - c(bits[["sumi.p"]]) / Denom)

    ned2l.dshape.p.shape.a <- zero0n  # Final; nothing to do
    ned2l.dpobs.a.shape.p  <- zero0n  # mB overwritten below [4279]
    ned2l.dpobs.a.shape.a  <- zero0n  # Final; nothing to do
    ned2l.dpobs.a.pstr.i   <- zero0n  # mB overwritten below
    ned2l.dpobs.a.shape.i  <- zero0n  # mB overwritten below
    ned2l.dpstr.i.shape.p  <- zero0n  # mB overwritten below
    ned2l.dpstr.i.shape.a  <- zero0n  # Final; nothing to do
    ned2l.dpstr.i.shape.i  <- zero0n  # mB overwritten below
    ned2l.dshape.a.shape.i <- zero0n  # Final; nothing to do

    ned2l.dshape.p2 <-
      onemrsDelta * (cond.EY.p / shape.p^2 +
                     A8.p * (1 - A8.p) / (1 - shape.p)^2 -
                     c(sumderiv2t.p + sumderiv2a.p) / Denom -
                    (c(sumderiv1t.p + sumderiv1a.p) / Denom)^2)
    if (linfla > 0)
      ned2l.dshape.p2 <- ned2l.dshape.p2 + onembothprobs *
        rowSums(onembothprobs * d1B.pi^2 / Delta - d2B.pi)

    if (lalter > 0) {
      wz22 <- if ( .lpobs.a == "logitlink" && linfla == 0) {
        pobs.a * (1 - pobs.a)
      } else {
        ned2l.dpobs.a2 <- 1 / pobs.a +
          onemrsDelta / onembothprobs^2 + (if (linfla > 0)
          rowSums(d0B.pi^2 / Delta) else 0)
        ned2l.dpobs.a2 * dpobs.a.deta^2
      }
    }
    if (lalter > 1)
      ned2l.dshape.a2 <- pobs.a * (
        rowSums(deriv1matrix.a^2  / deriv0matrix.a) / c(
        bits[["suma.a"]]) -
        (c(sumderiv1a.a) / c(bits[["suma.a"]]))^2)

    if (linfla > 0)
      ned2l.dpstr.i2 <-
        onemrsDelta / onembothprobs^2 +
        (if (linfla > 0)  # zz > 1 ??
        rowSums((deriv0matrix.i / c(bits[["sumi.i"]]) -
                 Deriv0Matrix.i / Denom)^2 / Delta) else 0)

    if (linfla > 0)
      ned2l.dpstr.i.shape.p <-
        rowSums(d1B.pi * (1 + onembothprobs *
                         (d0A.io - d0B.pi) / Delta))

    if (linfla > 1)
      ned2l.dshape.i2 <- pstr.i *
        rowSums(pstr.i * (d1A.io^2) / Delta - d2A.io)

    
    if (linfla > 1)
      ned2l.dpstr.i.shape.i <-
        rowSums(d1A.io * (pstr.i * (d0A.io - d0B.pi) / Delta - 1))


    if (linfla > 1)
      ned2l.dshape.p.shape.i <- pstr.i * onembothprobs *
                                  rowSums(d1A.io * d1B.pi / Delta)


    if (lalter > 0 && linfla > 0)
      ned2l.dpobs.a.pstr.i <- onemrsDelta / onembothprobs^2 -
        rowSums(d0B.pi * (d0A.io - d0B.pi) / Delta)


    if (linfla > 1)
      ned2l.dpobs.a.shape.i <-
        rowSums(-pstr.i * d1A.io * d0B.pi / Delta)

    
    if (linfla > 0)
      ned2l.dpobs.a.shape.p <-
        rowSums(d1B.pi * (1 - onembothprobs * d0B.pi / Delta))
   

    wz <- matrix(0, n, if (lalter && !linfla) M1 else M1*(M1+1)/2)

    wz[, iam(1, 1, M)] <- ned2l.dshape.p2 * dshape.p.deta^2

    if (linfla > 0)
      wz[, iam(1, 2, M)] <-
        ned2l.dpobs.a.shape.p * dpobs.a.deta * dshape.p.deta
    if (lalter > 0)
      wz[, iam(2, 2, M)] <- wz22
    if (lalter > 1)
      wz[, iam(3, 3, M)] <- ned2l.dshape.a2 * dshape.a.deta^2

    if (linfla > 0) {
      wz[, iam(Offset.i + 1, Offset.i + 1, M)] <-
        ned2l.dpstr.i2 * dpstr.i.deta^2
      wz[, iam(Offset.i + 1, 1, M)] <-
        ned2l.dpstr.i.shape.p * dpstr.i.deta * dshape.p.deta
      if (lalter > 0)
        wz[, iam(Offset.i + 1, 2, M)] <-
          ned2l.dpobs.a.pstr.i * dpobs.a.deta * dpstr.i.deta
    }

    if (linfla > 1) {
      wz[, iam(Offset.i + 2, Offset.i + 2, M)] <-
        ned2l.dshape.i2 * dshape.i.deta^2
      wz[, iam(Offset.i + 2, 1, M)] <-
        ned2l.dshape.p.shape.i * dshape.p.deta * dshape.i.deta
      wz[, iam(Offset.i + 2, 2, M)] <-
        ned2l.dpobs.a.shape.i * dpobs.a.deta * dshape.i.deta
      wz[, iam(Offset.i + 2, Offset.i + 1, M)] <-
        ned2l.dpstr.i.shape.i * dpstr.i.deta * dshape.i.deta
    }

    c(w) * wz
  }), list( .lpobs.a = lpobs.a, .lpstr.i = lpstr.i))))
}  # gaitlog.mix


















 gaitpoisson.mlm <-
  function(alter = NULL, inflate = NULL,
           truncate = NULL, max.support = Inf,  # Optional
           zero = c("pobs", "pstr"),  # Pruned of irrelevant values
           parallel.ap = FALSE,  # TRUE applies to the intercept
           parallel.ip = FALSE,  # TRUE applies to the intercept
           llambda = "loglink",
           type.fitted = c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
                           "prob.a", "prob.i", "prob.t", "lhs.prob"),
           imethod = 1,
           mux.inflate = 0.5,
           ilambda = NULL, ishrinkage = 0.95,
           probs.y = 0.35) {
  semigait.errorcheck(alter, inflate, truncate, max.support)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)
  if (lalter && linfla)
    stop("currently arguments 'alter' and 'inflate' cannot ",
         "both be assigned values")


  if (is.vector(zero) && is.character(zero) &&
      length(zero) == 2) {
    if (lalter)
      zero <- setdiff(zero, "pstr")
    if (linfla)
      zero <- setdiff(zero, "pobs")
  }

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")
  llambda.save <- llambda

  if (lalter + linfla == 0) zero <- NULL
  if (lalter + linfla + ltrunc == 0 && is.infinite(max.support))
    return(eval(substitute(
           poissonff(link = .llambda.save , zero = .zero ),
           list( .llambda.save = llambda.save, .zero = zero))))


  if (lalter && lalter <= 1 &&
    is.logical(parallel.ap) && parallel.ap)
    warning("Less than 2 altered values and hence 1 lambda, so ",
            "setting 'parallel.ap = TRUE' is meaningless")
  if (linfla && linfla <= 1 &&
      is.logical(parallel.ip) && parallel.ip)
    warning("Less than 2 inflated values and hence 1 lambda, so ",
            "setting 'parallel.ip = TRUE' is meaningless")


  type.fitted <- match.arg(type.fitted,
    c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
      "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]
  temp7 <- if (lalter) paste0("pobs", alter) else NULL
  if (linfla) temp7 <- paste0("pstr", inflate)
  tmp3 <- c(lambda = llambda,
            if (lalter) rep("multilogitlink",  lalter) else NULL,
            if (linfla) rep("multilogitlink",  linfla) else NULL)
  names(tmp3) <- c("lambda", temp7)





  
  blurb1 <- "P"
  if (lalter) blurb1 <- "Generally-altered P"
  if (linfla) blurb1 <- "Generally-inflated P"
  if (ltrunc) blurb1 <- "Generally-truncated P"
  if (lalter && ltrunc) blurb1 <- "Generally-altered and -truncated P"
  if (linfla && ltrunc) blurb1 <- "Generally-inflated and -truncated P"

  new("vglmff",
  blurb = c(blurb1, "oisson regression\n",
            "(GAIT-Pois(lambda)-MLM-MLM generally)\n\n",
            "Links:    ",
            namesof("lambda", llambda, earg = elambda, tag = FALSE),
            if (lalter)
            paste(",\n          multilogitlink(cbind(",
            paste(temp7, collapse = ", "), ", ",
            "\n                               ",
            "1 - ", paste(temp7, collapse = " - "),
            "))", sep = "") else NULL,
            if (length(inflate))
            paste(",\n          ",
            "multilogitlink(cbind(",
            paste(temp7, collapse = ", "),
            ", ",
            "\n                           ",
            "1 - ", paste(temp7, collapse = " - "),
            "))", ", ", sep = "") else NULL),

  constraints = eval(substitute(expression({


    if ( .lalter > 1) {
      use.mat <- matrix(0, 1 + .lalter , 2)
      use.mat[ 1, 1] <- 1  # lambda by itself
      use.mat[-1, 2] <- 1  # All MLM probabilities are equal
      constraints <- cm.VGAM(use.mat, x = x, bool = .parallel.ap ,
                             constraints = constraints,
                             apply.int = TRUE)  # FALSE
    }
  
    if ( .lalter )
      constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1 + .lalter )

    if ( .linfla > 1) {
      use.mat <- matrix(0, 1 + .linfla , 2)
      use.mat[ 1, 1] <- 1  # lambda by itself
      use.mat[-1, 2] <- 1  # All MLM probabilities are equal
      constraints <- cm.VGAM(use.mat, x = x, bool = .parallel.ip ,
                             constraints = constraints,
                             apply.int = TRUE)  # FALSE
    }
  
    if ( .linfla )
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1 + .linfla )
  }), list( .zero = zero,
            .parallel.ap = parallel.ap, .parallel.ip = parallel.ip,
            .lalter = lalter, .linfla = linfla ))),

  infos = eval(substitute(function(...) {
    alter <- as.vector( .alter )
    lalter <- length(alter)
    if (lalter)
      temp7 <- paste("pobs", alter, sep = "")
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    if (linfla)
      temp7 <- paste("pstr", inflate, sep = "")
    list(M1 = 1 + lalter + linfla,
         Q1 = 1,
         link = .tmp3 ,
         link1parameter = as.logical(lalter <= 1) &&
                          as.logical(linfla <= 1),
         mixture.links  = as.logical(lalter >  1) ||
                          as.logical(linfla >  1),
         alter = as.vector( .alter ),
         inflate = as.vector( .inflate ),
         truncate = as.vector( .truncate ),
         max.support = as.vector( .max.support ),
         expected = TRUE,
         multipleResponses = FALSE,  # poissonff called for TRUE case
         parameters.names = c("lambda",
                              if (lalter) temp7 else NULL,
                              if (linfla) temp7 else NULL),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted,
           .llambda = llambda, .elambda = elambda,
           .alter = alter, .inflate = inflate,
           .truncate = truncate, .max.support = max.support,
           .tmp3 = tmp3
         ))),

  initialize = eval(substitute(expression({
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    M1 <- 1 + lalter + linfla
    NOS <- NCOL(y)
    M <- NOS * M1

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,  # Since max.support = 9 is possible
              ncol.y.max = 1,
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





    if (linfla) {
      extra$y0 <- y0 <- matrix(0, n, linfla)
      for (jay in seq(linfla))
        extra$y0[, jay] <- y0[, jay] <- as.numeric(y == inflate[jay])
      extra$skip.these <- skip.these <- matrix(as.logical(y0), n, linfla)
      if (any((css <- colSums(skip.these)) == 0))
        stop("some 'inflate' argument values have no response values: ",
             paste(inflate[css == 0], collapse = ", "))          
    }


    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)
    extra$M1 <- M1

    
    if (TRUE) {
    fillerChar <- " "
    use.refLevel <- M1+1  # Assumes only one response
    allbut.refLevel <- (1:(M+1))[-use.refLevel]
    predictors.names <-
      paste("log(pstr.i[,", allbut.refLevel,
            "]", fillerChar, "/", fillerChar, "pstr.i[,",
            use.refLevel, "])", sep = "")

    temp7 <- paste0("pstr", inflate)  # "pstr" if is.null(inflate)
    }


    if (FALSE) {
    mynames1 <- paste("multinomiallink1(", temp7, ")", sep = "")
    mynames2 <- param.names("lambda", ncoly, skip1 = TRUE)
    predictors.names <-
        c(        mynames1,
          namesof(mynames2, .llambda , earg = .elambda , tag = FALSE))
    }  # FALSE


    mynames1 <- param.names("lambda", ncoly, skip1 = TRUE)
    mynames2 <- if (lalter) {
      temp7 <- paste0("pobs", alter)
      denom.char <- paste0("1-", paste0(temp7, collapse = "-"))
      paste("log(", temp7, "/(", denom.char, "))", sep = "")
    } else if (linfla) {
      denom.char <- paste0("1-", paste0(temp7, collapse = "-"))
      paste("log(", temp7, "/(", denom.char, "))", sep = "")
    } else NULL
    predictors.names <-
        c(namesof(mynames1, .llambda , earg = .elambda ,
                  tag = FALSE),
          mynames2)





    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                             imu = .ilambda ,  # x = x,
                             ishrinkage = .ishrinkage ,
                             pos.only = TRUE, probs.y = .probs.y )

      if (lalter) {
        phimat <- colMeans(skip.these)  # yettodo: weight this by 'w'
        phimat <- matrix(phimat, n, lalter, byrow = TRUE)
        etastart <-  multilogitlink(cbind(phimat, 1 - rowSums(phimat)))
      }

      if (linfla) {
        phimat <- colSums(c(w) * y0) / c(colSums(w))  # Wted by 'w'.
        onemphis <- 1 - sum(phimat)

        phimat <- phimat * (1 - onemphis)
        if (FALSE)
        phimat <- phimat * ( .mux.inflate )
        phimat <- matrix(phimat, n, linfla, byrow = TRUE)
        etastart <-  multilogitlink(cbind(phimat,
                                          1 - rowSums(phimat)))
      }

      
      etastart <-
        cbind(theta2eta(lambda.init, .llambda , earg = .elambda ),
              etastart)
    }
  }), list( .llambda = llambda, .elambda = elambda,
            .ilambda = ilambda,
            .mux.inflate = mux.inflate,  # .gpstr0 = gpstr0,
            .truncate = truncate, .max.support = max.support, 
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .imethod = imethod,
            .alter = alter, .inflate = inflate,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <-
      if (length(extra$type.fitted)) extra$type.fitted else {
        warning("cannot find 'type.fitted'. Returning the 'mean'.")
        "mean"
      }
    pobs.a <- omegamat <- 0     # In case lalter == 0
    pstr.i <- phimat   <- 0     # In case linfla == 0
    type.fitted <-
      match.arg(type.fitted,
                c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
                  "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]
    eta <- as.matrix(eta)
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    max.support <- ( .max.support )
    M1 <- 1 + lalter + linfla
    NOS <- NCOL(eta) / M1
    lambda <-
      eta2theta(if (lalter + linfla) eta[, 1, drop = FALSE] else eta,
                .llambda , earg = .elambda )

    if (lalter) {
      omegamat <- multilogitlink(eta[, -1, drop = FALSE],
                    refLevel = NCOL(eta),  # Assumes 1 response
                    inverse = TRUE)  # rowSums == 1
      ynames.Pobs.a <- c(as.character(alter), "(Others)")
      dimnames(omegamat) <- list(rownames(eta), ynames.Pobs.a)
    }  # lalter
    
    if (linfla) {
      phimat <- multilogitlink(eta[, -1, drop = FALSE],
                  refLevel = NCOL(eta),  # Assumes 1 response
                  inverse = TRUE)  # rowSums == 1
      ynames.Pstr.i <- c(as.character(inflate), "(Others)")
      dimnames(phimat) <- list(rownames(eta), ynames.Pstr.i)
    }  # linfla
    Bits <-
      moments.pois.gait(lambda, mlm = TRUE,
                        pobs.a = omegamat[, -ncol(omegamat)],
                        pstr.i =   phimat[, -ncol(phimat)],
                        alter = alter, inflate = inflate,
                        truncate = truncate, max.support = max.support)

    ans <- switch(type.fitted,
      "mean"       = Bits[["mean"]],  # Unconditional mean
      "pobs.a"  = omegamat,  #[, -ncol(omegamat)], # cbind(omega[1:L_A])
      "pstr.i"  =   phimat,  #[, -ncol(phimat  )], # cbind(  phi[1:L_I])
      "Pobs.a"  = omegamat,  #[, -ncol(omegamat)], # Same as "pobs.a"
      "Pstr.i"  =   phimat[, -ncol(phimat)] +
                     (1 - phimat[, "(Others)"]) *
                     dpois(matrix(inflate,  nrow(eta), linfla,
                                  byrow = TRUE),
                     matrix(lambda, nrow(eta), linfla)) / c(
                     Bits[["lhs.prob"]] - Bits[["suma.p"]] -
                     Bits[["sumt.p"]]),
      "prob.a"     = Bits[["suma.p"]],
      "prob.i"     = Bits[["sumi.p"]],
      "prob.t"     = Bits[["sumt.p"]],  # Truncated probs w/o RHS
      "lhs.prob"   = Bits[["lhs.prob"]])  # 1 - Pr(y <= max.support)
    label.cols.y(ans,
        colnames.y = switch(type.fitted,
                            "pobs.a" =,
                            "Pobs.a" = ynames.Pobs.a,
                            "pstr.i" = ynames.Pstr.i,
                     "Pstr.i" = ynames.Pstr.i[-length(ynames.Pstr.i)],
                            extra$colnames.y),
        NOS = NOS)
  }, list( .llambda = llambda, .elambda = elambda,
           .truncate = truncate, .max.support = max.support, 
           .alter = alter, .inflate = inflate ))),

  last = eval(substitute(expression({
    tmp9 <- if (lalter) rep_len("multilogitlink", lalter) else
            if (linfla) rep_len("multilogitlink", linfla) else NULL
    misc$link  <- c( .llambda , tmp9)
    names(misc$link) <-
      c(mynames1, mynames2)  # [interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1)
    names(misc$earg) <- names(misc$link)

    misc$earg[[1]] <- .elambda  # First one
    if (lalter) {
      for (ii in seq(M1 - 1)) {
          misc$earg[[1 + ii]] <- list(M = M - 1,  # M * NOS,
                                      refLevel = M)  # M * NOS
      }  # ii
    }  # lalter
    if (linfla) {
      for (ii in seq(M1 - 1)) {
          misc$earg[[1 + ii]] <- list(M = M - 1,  # M * NOS,
                                      refLevel = M)  # M * NOS
      }  # ii
    }  # linfla
  }), list( .llambda = llambda, .elambda = elambda,
            .alter = alter, .inflate = inflate ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    eta <- as.matrix(eta)
    alter <- as.vector ( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    lambda <-
      eta2theta(if (lalter + linfla) eta[, 1, drop = FALSE] else eta,
                .llambda , earg = .elambda )
    if (lalter)
      pobs.a <- multilogitlink(eta[, -1, drop = FALSE],
                  refLevel = NCOL(eta),  # Assumes one response
                  inverse = TRUE)
    if (linfla)
      pstr.i <- multilogitlink(eta[, -1, drop = FALSE],
                  refLevel = NCOL(eta),  # Assumes 1 response
                  inverse = TRUE)

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitpois(y, lambda.p = lambda, log = TRUE,
                  truncate = .truncate ,  # byrow.arg = FALSE,
                  max.support = .max.support ,
                  alter.mlm = alter , inflate.mlm = inflate,
                  pobs.mlm.a = if (lalter)
                    pobs.a[, -NCOL(pobs.a), drop = FALSE] else 0,
                  pstr.mlm.i = if (linfla)
                    pstr.i[, -NCOL(pstr.i), drop = FALSE] else 0)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llambda = llambda, .elambda = elambda,
           .truncate = truncate, .max.support = max.support, 
           .alter = alter, .inflate = inflate ))),
  vfamily = c("gaitpoisson.mlm"),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    alter <- as.vector ( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    pobs.a <- if (lalter)
        multilogitlink(eta[, -1, drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE) else 0.5  # An okay value
    pstr.i <- if (linfla)
        multilogitlink(eta[, -1, drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE) else 0.5  # An okay value
    eta <- as.matrix(eta)
    lambda <-
      eta2theta(if (lalter + linfla) eta[, 1, drop = FALSE] else eta,
                .llambda , earg = .elambda )
    okay1 <- all(is.finite(lambda)) && all(0 < lambda) &&
             all(lambda < .max.support ) &&
             all(is.finite(pobs.a)) && all(0 < pobs.a & pobs.a < 1) &&
             all(is.finite(pstr.i)) && all(0 < pstr.i & pstr.i < 1)
    okay1
  }, list( .llambda = llambda, .elambda = elambda,
           .truncate = truncate, .max.support = max.support, 
           .alter = alter, .inflate = inflate ))),

  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    alter <- as.vector ( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    lambda <- eta2theta(if (lalter + linfla)
                        eta[, 1, drop = FALSE] else eta,
                        .llambda , earg = .elambda )
    pobs.a <- if (lalter)
        multilogitlink(eta[, -1, drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE) else 0.5  # Any okay value
    pstr.i <- if (linfla)
        multilogitlink(eta[, -1, drop = FALSE],
                refLevel = NCOL(eta),  # Assumes one response
                inverse = TRUE) else 0  # Last coln == "(Other)"
    rgaitpois(nsim * length(lambda), lambda.p = lambda,
        pstr.mlm.i = if (linfla) pstr.i[, -NCOL(pstr.i)] else pstr.i,
        pobs.mlm.a = pobs.a, alter = alter, inflate = inflate,
        truncate = .truncate , max.support = .max.support )
  }, list( .llambda = llambda, .elambda = elambda,
           .truncate = truncate, .max.support = max.support, 
           .alter = alter, .inflate = inflate ))),

  deriv = eval(substitute(expression({
    eta <- as.matrix(eta)
    max.support <- ( .max.support )
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    M1 <- 1 + lalter + linfla
    NOS <- ncol(eta) / M1  # extra$NOS
    if (NOS != 1) stop("can only handle 1 response")

    is.altered <- if (lalter) {
      y0 <- extra$y0
      skip <- extra$skip.these
      rowSums(skip) > 0  # TRUE if any(y %in% avec)
    } else rep(FALSE, NROW(eta))


    is.inflated <- if (linfla) {  # Unpack
      y0 <- extra$y0
      skip <- extra$skip.these
      rowSums(skip) > 0  # TRUE if any(y %in% ivec)
    } else rep(FALSE, NROW(eta))

    lambda <-
      eta2theta(if (lalter + linfla) eta[, 1, drop = FALSE] else eta,
                .llambda , earg = .elambda )
    phimat <- if (lalter)
        multilogitlink(eta[, -1, drop = FALSE],
                       refLevel = NCOL(eta),  # Assumes 1 response
                       inverse = TRUE) else 0
    pstr.i <-  # phimat is the same as pstr.i, Pr(inflated)
      if (linfla)
      multilogitlink(eta[, -1, drop = FALSE],
                     refLevel = NCOL(eta),  # Assumes 1 response
                     inverse = TRUE) else 0
    onempstr.i <- if (linfla)
                    pstr.i[, ncol(pstr.i)] else 1  # Vector
    pstr.i <- if (linfla) pstr.i[, -ncol(pstr.i), drop = FALSE] else 0
    pmf.deriv1 <- function(y, lambda)
      dpois(y-1, lambda) - dpois(y, lambda)
    pmf.deriv2 <- function(y, lambda)
      dpois(y-2, lambda) - 2 * dpois(y-1, lambda) + dpois(y, lambda)

    bits <-
      moments.pois.gait(lambda, pstr.i = pstr.i,
                        pobs.a = if (lalter)
                        phimat[, -NCOL(phimat), drop = FALSE] else phimat,
                        alter = alter, inflate = inflate,
                        truncate = truncate, max.support = max.support,
                        lambda.a = lambda,
                        mlm = TRUE)  # lambda.i = lambda,


    sumderiv1a.p <- sumderiv2a.p <- matrix(0, NROW(eta), NOS)
    if (lalter)
      for (aval in alter) {
        sumderiv1a.p <- sumderiv1a.p + pmf.deriv1(aval, lambda)
        sumderiv2a.p <- sumderiv2a.p + pmf.deriv2(aval, lambda)
      }

    
    if (linfla)
      Denom <- c(bits[["lhs.prob"]] - bits[["sumt.p"]])  # As alter==NULL

    if (linfla > 0) {
      Deriv0Matrix <-  # wrt parent distribution
      Deriv1Matrix <- Deriv2Matrix <- matrix(0, NROW(eta), linfla)
      for (jay in seq(linfla)) {
        ival <- inflate[jay]
        pmf.p <- dpois(ival, lambda)
        Deriv0Matrix[, jay] <- pmf.p
        Deriv1Matrix[, jay] <- pmf.deriv1(ival, lambda)
        Deriv2Matrix[, jay] <- pmf.deriv2(ival, lambda)  # / pmf.p
      }  # jay
    }  # linfla
      
    sumderiv1t.p <- sumderiv2t.p <- matrix(0, NROW(eta), NOS)
    if (ltrunc)
      for (tval in truncate) {
        sumderiv1t.p <- sumderiv1t.p + pmf.deriv1(tval, lambda)
        sumderiv2t.p <- sumderiv2t.p + pmf.deriv2(tval, lambda)
      }

    sumderiv1t.p <- sumderiv1t.p + dpois( .max.support    , lambda)
    sumderiv2t.p <- sumderiv2t.p + dpois( .max.support - 1, lambda) -
                                   dpois( .max.support    , lambda)

    pobs.a <- if (lalter)
      1 - phimat[, NCOL(phimat), drop = FALSE] else 0
    if (lalter) {
      dl.deta <- skip - phimat[, -M, drop = FALSE]

      Denom <- bits[["lhs.prob"]] - bits[["suma.p"]] - bits[["sumt.p"]]
      dl.dlambda <- (if (lalter) (1 - is.altered) else 1) *
        (y / lambda - 1 + (sumderiv1a.p  + sumderiv1t.p) / Denom)

    }  # lalter






    if (linfla > 0) {
      d0B.pi <- Deriv0Matrix / Denom
      Delta <- pstr.i + onempstr.i * d0B.pi  # n x linfla
      d1A.io <- matrix(0, n, linfla)
      d2A.io <- matrix(0, n, linfla)
      d1B.pi <- Deriv1Matrix / Denom +
                Deriv0Matrix * c(sumderiv1t.p) / Denom^2
      d2B.pi <- Deriv2Matrix / Denom +
                2 * Deriv1Matrix * c(sumderiv1t.p) / Denom^2 +
                Deriv0Matrix * c(sumderiv2t.p) / Denom^2 +
                2 * Deriv0Matrix * c(sumderiv1t.p^2) / Denom^3

      tmp5b <- y / lambda - 1
      zero0n <- rep(0, n)
      dl.dlambda <- ifelse(is.inflated,
                           zero0n,  # Replace soon
                           tmp5b + c(sumderiv1t.p) / Denom)
      tmp0 <- -1 / onempstr.i
    }  # linfla > 0

    
    if (linfla > 0) {
      dl.dpstr.i <- matrix(tmp0, n, linfla)  # Some elts 2 b replaced

      for (jay in seq(linfla)) {
        ival <- inflate[jay]
        is.infl.j <- extra$skip.these[, jay]  # Logical vector
        tmp7 <- onempstr.i * d1B.pi[, jay] / Delta[, jay]
        dl.dlambda[is.infl.j] <- tmp7[is.infl.j]



        tmp9 <- (0 - d0B.pi[, jay]) / Delta[, jay]
          dl.dpstr.i[is.infl.j, ] <- tmp9[is.infl.j]
        tmp8 <- (1 - d0B.pi[, jay]) / Delta[, jay]
        dl.dpstr.i[is.infl.j, jay] <- tmp8[is.infl.j]
      }  # jay

      dl.deta <- matrix(0, n, linfla)
      for (jay in seq(linfla)) {
        for (sss in seq(linfla)) {
          dl.deta[, jay] <- dl.deta[, jay] +
            pstr.i[, sss] * ((sss == jay) - pstr.i[, jay]) *
            dl.dpstr.i[, sss]
        }  # sss
      }  # jay
    }  # linfla > 0


    dlambda.deta <- dtheta.deta(lambda, .llambda , .elambda )
    ans <- cbind(dl.dlambda * dlambda.deta,
                 if (lalter) c(w) * dl.deta else NULL,
                 if (linfla) c(w) * dl.deta else NULL)
    ans
  }), list( .llambda = llambda, .elambda = elambda,
            .truncate = truncate, .max.support = max.support, 
            .alter = alter, .inflate = inflate ))),

  weight = eval(substitute(expression({
    cond.EY.p <- 
      (lambda - c(bits[["SumI.p"]]) - c(bits[["SumT.p"]])) / c(
       bits[["lhs.prob"]] - bits[["sumi.p"]] - bits[["sumt.p"]])


    onemrsDelta.0 <- onempstr.i * (1 - c(bits[["sumi.p"]]) / Denom)
    if (linfla > 0) {
      onemrsDelta.1 <- 1 - rowSums(Delta)  # linfla > 0 needed
      if (max(abs(onemrsDelta.1 - onemrsDelta.0)) > 1e-7)
        warning("check: the two estimates of onemrsDelta differ")
    }



    onemrsDelta <- onemrsDelta.0  # .1 was good for ned2l.dee.i2



    ned2l.dlambda2 <- if (lalter) {
      (1 - pobs.a) *
      ((lambda - bits[["SumA.p"]] -  bits[["SumT.p"]]) / (
       Denom * lambda^2) -
        (sumderiv2t.p + sumderiv2a.p) / Denom -
       ((sumderiv1a.p + sumderiv1t.p) / Denom)^2)
     } else {
    onemrsDelta.0 * (cond.EY.p / lambda^2 -
      sumderiv2t.p / Denom - (sumderiv1t.p / Denom)^2) +
      (if (linfla > 0) onempstr.i *
         rowSums(onempstr.i * d1B.pi^2 / Delta - d2B.pi) else 0)
    }


    if (lalter > 0) {
      MM12 <- M1 * (M1 + 1) / 2  # Assumes NOS == 1.
      wz4 <- matrix(0, n, MM12)  # A full matrix
      if (lalter == 1) {
        wz4[, 1] <- phimat[, 1] * (1 - phimat[, 1])
      } else {  # lalter > 1
        index <- iam(NA, NA, M - 1, both = TRUE, diag = TRUE)
        wz4 <- -phimat[, index$row] * phimat[, index$col]
        wz4[, 1:lalter] <- wz4[, 1:lalter] + phimat[, -M]
      }
        wz <- wz.merge(ned2l.dlambda2 * dlambda.deta^2,
                       wz4,
                       M1 = 1, M2 = M1 - 1)  # rm.trailing.cols = F
    } # lalter > 0



    

    if (linfla > 0) {
      ned2l.dpstr.i2 <-
        matrix(onemrsDelta.1 / onempstr.i^2, n, linfla*(linfla+1)/2)
      for (uuu in seq(linfla))
        for (sss in seq(linfla))
          ned2l.dpstr.i2[, iam(uuu, uuu, linfla)] <-
          ned2l.dpstr.i2[, iam(uuu, uuu, linfla)] +
            ((sss == uuu) - d0B.pi[, sss])^2 / Delta[, sss]
      if (linfla > 1) {
        for (uuu in 1:(linfla-1))
          for (vvv in (uuu+1):linfla)
            for (sss in seq(linfla))
              ned2l.dpstr.i2[, iam(uuu, vvv, linfla)] <-
              ned2l.dpstr.i2[, iam(uuu, vvv, linfla)] +
              ((sss == uuu) - d0B.pi[, sss]) *
              ((sss == vvv) - d0B.pi[, sss]) / Delta[, sss]
      }  # if (linfla > 1)
    }  # linfla > 0

    
    if (linfla > 0) {
      ned2l.dpd <- matrix(0, n, linfla)
      for (vvv in seq(linfla))
        for (sss in seq(linfla))
          ned2l.dpd[, vvv] <- ned2l.dpd[, vvv] +
          d1B.pi[, sss] * (1 + onempstr.i *
          (max(0, sss == vvv) - d0B.pi[, sss]) / Delta[, sss])
    }  # linfla > 0


    if (linfla > 0) {
      ned2l.dpeta <- matrix(0, n, linfla)
      for (jay in seq(linfla)) {
        for (sss in seq(linfla)) {
          ned2l.dpeta[, jay] <- ned2l.dpeta[, jay] +
            pstr.i[, sss] * (max(0, sss == jay) - pstr.i[, jay]) *
            ned2l.dpd[, sss]
        }  # sss
      }  # jay
    }  # linfla > 0


    if (linfla > 0) {
      ned2l.dee.i2 <- matrix(0, n, linfla*(linfla+1)/2)
      for (uuu in seq(linfla)) {
        for (vvv in uuu:linfla) {
          for (sss in seq(linfla)) {
            for (ttt in seq(linfla)) {
              ned2l.dee.i2[, iam(uuu, vvv, linfla)] <-
              ned2l.dee.i2[, iam(uuu, vvv, linfla)] +
                pstr.i[, sss] * (max(0, sss == uuu) - pstr.i[, uuu]) *
                ned2l.dpstr.i2[, iam(sss, ttt, linfla)] *
                pstr.i[, ttt] * (max(0, ttt == vvv) - pstr.i[, vvv])
            }  # ttt
          }  # sss
        }  # vvv
      }  # uuu
    }  # linfla > 0


    if (linfla > 0) {
      wz <- wz.merge(ned2l.dlambda2 * dlambda.deta^2,
                     ned2l.dee.i2,
                     M1 = 1, M2 = linfla, rm.trailing.cols = FALSE)
      for (jay in seq(linfla))
        wz[, iam(1, 1+jay, M)] <- dlambda.deta * ned2l.dpeta[, jay]
    }


    if (lalter + linfla == 0) {  # Ordinary Poisson if truncate == NULL
      wz <- ned2l.dlambda2 * dlambda.deta^2
    }

    c(w) * wz
  }), list( .alter = alter, .inflate = inflate ))))
}  # gaitpoisson.mlm













 gaitpoisson.mix <-
  function(alter = NULL,
           inflate = NULL,  
           truncate = NULL, max.support = Inf,
           zero = c("pobs.a", "pstr.i"),  # Pruned later if necessary
           eq.ap = FALSE,  # TRUE applies to the intercept
           eq.ip = FALSE,  # TRUE applies to the intercept
           llambda.p = "loglink",
           lpobs.a = "logitlink",
           llambda.a = "loglink",
           lpstr.i = "logitlink",
           llambda.i = "loglink",
           type.fitted = c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
                           "prob.a", "prob.i", "prob.t", "lhs.prob"),
           imethod = 1,
           ilambda.p = NULL,
           ilambda.a = NULL, ilambda.i = NULL,
           ipobs.a = NULL,  # 0.25, 
           ipstr.i = NULL,  # 0.25, 
           ishrinkage = 0.95,
           probs.y = 0.35) {
  semigait.errorcheck(alter, inflate, truncate, max.support)
  lowsup <- 0

  if (!is.logical(eq.ap) || length(eq.ap) != 1)
    stop("argument 'eq.ap' must be a single logical")
  if (!is.logical(eq.ip) || length(eq.ip) != 1)
    stop("argument 'eq.ip' must be a single logical")

  llambda.p <- as.list(substitute(llambda.p))
  elambda.p <- link2list(llambda.p)
  llambda.p <- attr(elambda.p, "function.name")
  lpobs.a <- as.list(substitute(lpobs.a))  # \omega
  epobs.a <- link2list(lpobs.a)
  lpobs.a <- attr(epobs.a, "function.name")
  llambda.a <- as.list(substitute(llambda.a))
  elambda.a <- link2list(llambda.a)
  llambda.a <- attr(elambda.a, "function.name")

  lpstr.i <- as.list(substitute(lpstr.i))  # \phi
  epstr.i <- link2list(lpstr.i)
  lpstr.i <- attr(epstr.i, "function.name")
  llambda.i <- as.list(substitute(llambda.i))
  elambda.i <- link2list(llambda.i)
  llambda.i <- attr(elambda.i, "function.name")

  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)

  if (is.vector(zero) && is.character(zero) && length(zero) == 2) {
    if (lalter && linfla == 0)
      zero <- setdiff(zero, "pstr.i")
    if (linfla && lalter == 0)
      zero <- setdiff(zero, "pobs.a")
  }

  if (lalter + linfla == 0) zero <- NULL
  llambda.p.save <- llambda.p
  if (lalter + linfla + ltrunc == 0 && is.infinite(max.support))
    return(eval(substitute(
           poissonff(link = .llambda.p.save , zero = .zero ),
           list( .llambda.p.save = llambda.p.save, .zero = zero))))


  if (lalter == 1 && eq.ap)
    warning("Less than 2 altered values and hence 1 lambda, so ",
            "setting 'eq.ap = TRUE' is meaningless")
  if (linfla == 1 && eq.ip)
    warning("Less than 2 inflated values and hence 1 lambda, so ",
            "setting 'eq.ip = TRUE' is meaningless")

  type.fitted <- match.arg(type.fitted,
    c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
      "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]

  tmp3.TF <- c(TRUE, lalter > 0, lalter > 1, linfla > 0, linfla > 1)
  tmp3 <- c(lambda.p = llambda.p,  # Version 2
            pobs.a   = lpobs.a, lambda.a = llambda.a,
            pstr.i   = lpstr.i, lambda.i = llambda.i)  # [tmp3.TF]
      
  blurb1 <- "P"
  if (lalter) blurb1 <- "Generally-altered P"
  if (linfla) blurb1 <- "Generally-inflated P"
  if (ltrunc) blurb1 <- "Generally-truncated P"
  if ( lalter &&  linfla && !ltrunc)
    blurb1 <- "Generally-altered and -inflated P"
  if ( lalter && !linfla &&  ltrunc)
    blurb1 <- "Generally-altered and -truncated P"
  if (!lalter &&  linfla &&  ltrunc)
    blurb1 <- "Generally-inflated and -truncated P"
  if ( lalter &&  linfla &&  ltrunc)
    blurb1 <- "Generally-altered, -inflated and -truncated P"

  new("vglmff",
  blurb = c(blurb1, "oisson regression\n",
            "(GAIT-Pois(lambda.p)-Pois(lambda.a)-Pois(lambda.i) ",
            "mixture generally)\n\n",
            "Links:    ",
            namesof("lambda.p",  llambda.p, elambda.p, tag = FALSE),
            if (lalter > 0) c(", ",
            namesof("pobs.a",    lpobs.a,   epobs.a,   tag = FALSE)),
            if (lalter > 1) c(", ",
            namesof("lambda.a",  llambda.a, elambda.a, tag = FALSE)),
            if (lalter && linfla) ", \n        ",
            if (linfla > 0) c(  if (lalter) "  " else ", ",
            namesof("pstr.i",    lpstr.i,   epstr.i,   tag = FALSE)),
            if (linfla > 1) c(", ",
            namesof("lambda.i",  llambda.i, elambda.i, tag = FALSE)),
            "\n"),
  constraints = eval(substitute(expression({

    use.mat <- if ( ( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,1,0,1,  0,1,0,0,0,  0,0,0,1,0), 5, 3) else
    if (!( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,0,0,1,  0,1,0,0,0,  0,0,1,0,0,
               0,0,0,1,0), 5, 4) else
    if ( ( .eq.ap ) && !( .eq.ip ))
      matrix(c(1,0,1,0,0,  0,1,0,0,0,  0,0,0,1,0,
               0,0,0,0,1), 5, 4) else
      diag(5)

    use.mat <- use.mat[ .tmp3.TF , , drop = FALSE]
    M1 <- nrow(use.mat)
    for (jay in ncol(use.mat):1)
      if (all(use.mat[, jay] == 0))
        use.mat <- use.mat[, -jay, drop = FALSE]
    
    constraints <- cm.VGAM(use.mat, x = x,
                           bool = .eq.ap | .eq.ip ,
                           constraints = constraints,
                           apply.int = TRUE)  # FALSE

    if ( .lalter + .linfla > 0)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = M1)
  }), list( .zero = zero, .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
            .lalter = lalter, .linfla = linfla,
            .eq.ap = eq.ap, .eq.ip = eq.ip ))),

  infos = eval(substitute(function(...) {
    list(Q1 = 1,
         M1 = sum( .tmp3.TF ),
         link = c( .tmp3 )[ .tmp3.TF ] ,
         link1parameter = TRUE,
         mixture.links = FALSE,
         alter = as.vector( .alter ),
         inflate = as.vector( .inflate ),
         truncate = as.vector( .truncate ),
         max.support = as.vector( .max.support ), 
         Support  = c( .lowsup , Inf, 1),
         eq.ap = .eq.ap , eq.ip = .eq.ip ,
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE,
         parameters.names = names(c( .tmp3 )[ .tmp3.TF ]),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero, .tmp3 = tmp3, .tmp3.TF = tmp3.TF,
           .type.fitted = type.fitted,
           .eq.ap = eq.ap, .eq.ip = eq.ip,
           .alter = alter, .inflate = inflate, .lowsup = lowsup,
           .truncate = truncate, .max.support = max.support ))),

  initialize = eval(substitute(expression({
    truncate <- as.vector( .truncate )
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    M1 <- sum( .tmp3.TF )
    NOS <- NCOL(y)  # Only 1 currently
    M <- NOS * M1
    tmp3 <- ( .tmp3 )
    tmp3.TF <- ( .tmp3.TF )
    ntmp3 <- names(tmp3)

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
    glist <- y.gait.check(alter, inflate, truncate, .max.support , y)
    css.a <- glist$css.a
    css.i <- glist$css.i
    extra$skip.a <- glist$skip.a
    extra$skip.i <- glist$skip.i

    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- ( .type.fitted )
    extra$colnames.y  <- colnames(y)
    extra$M1 <- M1

   predictors.names <-
    c(namesof(ntmp3[1], .llambda.p , earg = .elambda.p , tag = FALSE),
      namesof(ntmp3[2], .lpobs.a   , earg = .epobs.a   , tag = FALSE),
      namesof(ntmp3[3], .llambda.a , earg = .elambda.a , tag = FALSE),
      namesof(ntmp3[4], .lpstr.i   , earg = .epstr.i   , tag = FALSE),
      namesof(ntmp3[5], .llambda.i , earg = .elambda.i , tag = FALSE))[
      tmp3.TF]  # Relevant subset



    

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                             imu = .ilambda.p ,  # x = x, pos.only = T,
                             ishrinkage = .ishrinkage ,
                             probs.y = .probs.y )
      pobs.a.init <-  # Unneeded? Low is best
      pstr.i.init <- 0.025 / (1 + lalter + linfla)
                               
      if (lalter > 0) {
        pobs.a.init <-  rep(if (length( .ipobs.a )) .ipobs.a else
                            sum(css.a) / n, length = n)
      }
      if (lalter > 1)
        lamb.a.init <- if (length( .ilambda.a ))
          rep( .ilambda.a , n) else lambda.init

      if (linfla > 0)
        pstr.i.init <-  rep(if (length( .ipstr.i )) .ipstr.i else
                            0.5 * sum(css.i) / n, length = n)
      if (linfla > 1)
        lamb.i.init <- if (length( .ilambda.i ))
          rep( .ilambda.i , n) else lambda.init

      while (any((vecTF <- pobs.a.init + pstr.i.init > 0.95))) {
        pobs.a.init[vecTF] <- 0.875 * pobs.a.init[vecTF]
        pstr.i.init[vecTF] <- 0.875 * pstr.i.init[vecTF]
      }

      etastart <-
        cbind(theta2eta(lambda.init, .llambda.p , earg = .elambda.p ),
              if (lalter <= 0) NULL else
              theta2eta(pobs.a.init, .lpobs.a   , earg = .epobs.a   ),
              if (lalter <= 1) NULL else
              theta2eta(lamb.a.init, .llambda.a , earg = .elambda.a ),
              if (linfla <= 0) NULL else
              theta2eta(pstr.i.init, .lpstr.i   , earg = .epstr.i   ),
              if (linfla <= 1) NULL else
              theta2eta(lamb.i.init, .llambda.i , earg = .elambda.i ))
    }
  }), list(
    .llambda.p = llambda.p, .llambda.i = llambda.i, .lpstr.i = lpstr.i,
    .elambda.p = elambda.p, .elambda.i = elambda.i, .epstr.i = epstr.i,
    .ilambda.p = ilambda.p, .ilambda.i = ilambda.i, .ipstr.i = ipstr.i,
                            .llambda.a = llambda.a, .lpobs.a = lpobs.a,
                            .elambda.a = elambda.a, .epobs.a = epobs.a,
    .ilambda.p = ilambda.p, .ilambda.a = ilambda.a, .ipobs.a = ipobs.a,
    .ishrinkage = ishrinkage, .probs.y = probs.y,
    .alter = alter, .inflate = inflate,
    .truncate = truncate, .max.support = max.support,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .imethod = imethod, .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <-
     if (length(extra$type.fitted)) extra$type.fitted else {
       warning("cannot find 'type.fitted'. Returning the 'mean'.")
                     "mean"
     }
   
    type.fitted <-
      match.arg(type.fitted,
                c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
                  "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]
    eta <- as.matrix(eta)  # Needed when linfla == 0
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate) 
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    max.support <- ( .max.support )
    M1 <- sum( .tmp3.TF )
    if ((NOS <- NCOL(eta) / M1) != 1)
      stop("Currently NOS must be 1")

    iptr <- 1  # Unpack
    lambda.p <- eta2theta(eta[, iptr], .llambda.p , earg = .elambda.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a   , earg = .epobs.a )
    lambda.a <- if (lalter <= 1) lambda.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .llambda.a , earg = .elambda.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i   , earg = .epstr.i )
    lambda.i <- if (linfla <= 1) lambda.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .llambda.i , earg = .elambda.i )
    Bits <- moments.pois.gait(lambda.p = lambda.p, mlm = FALSE,
                              alter = alter, inflate = inflate,
                              truncate = truncate,
                              max.support = max.support,
                              pobs.a = pobs.a, lambda.a = lambda.a,
                              pstr.i = pstr.i, lambda.i = lambda.i)
    if (type.fitted == "Pobs.a") {
      if (lalter == 0) stop("no altered values!")
      proportion.mat <-
        dpois(matrix(alter,    NROW(eta), lalter, byrow = TRUE),
              matrix(lambda.a, NROW(eta), lalter)) / c(Bits[["suma.a"]])
    }
    if (type.fitted == "Pstr.i") {
      if (linfla == 0) stop("no inflated values!")
      proportion.mat <-
        dpois(matrix(inflate,  nrow(eta), linfla, byrow = TRUE),
              matrix(lambda.i, nrow(eta), linfla)) / c(Bits[["sumi.i"]])
    }
  
    ans <- switch(type.fitted,
      "mean"     = Bits[["mean"]],
      "pobs.a"   = pobs.a,
      "pstr.i"   = pstr.i,
      "Pobs.a"   = pobs.a * proportion.mat,  # matrix
      "Pstr.i"   = c(pstr.i) * proportion.mat +
                   c(1 - pobs.a - pstr.i) *
           dpois(matrix(inflate,  nrow(eta), linfla, byrow = TRUE),
                 matrix(lambda.p, nrow(eta), linfla)) / c(
           Bits[["lhs.prob"]] - Bits[["suma.p"]] - Bits[["sumt.p"]]),
      "prob.a"   = Bits[["suma.p"]],
      "prob.i"   = Bits[["sumi.p"]],
      "prob.t"   = Bits[["sumt.p"]],
      "lhs.prob" = Bits[["lhs.prob"]])
   ynames.Pobs.a <- as.character(alter)    # Works with NULLs
   ynames.Pstr.i <- as.character(inflate)  # Works with NULLs
   label.cols.y(ans,
        colnames.y = switch(type.fitted,
                            "Pobs.a" = ynames.Pobs.a,
                            "Pstr.i" = ynames.Pstr.i,
                            extra$colnames.y),
        NOS = NOS)
  }, list(
    .llambda.p = llambda.p, .llambda.i = llambda.i, .lpstr.i = lpstr.i,
    .elambda.p = elambda.p, .elambda.i = elambda.i, .epstr.i = epstr.i,
                            .llambda.a = llambda.a, .lpobs.a = lpobs.a,
                            .elambda.a = elambda.a, .epobs.a = epobs.a,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .truncate = truncate, .max.support = max.support, 
    .alter = alter, .inflate = inflate ))),

  last = eval(substitute(expression({
    misc$link  <- c( .tmp3 )[ .tmp3.TF ]
    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    misc$earg[[(iptr <- 1)]] <- .elambda.p  #
    if (lalter > 0)
      misc$earg[[(iptr <- iptr + 1)]] <- .epobs.a  #
    if (lalter > 1)
      misc$earg[[(iptr <- iptr + 1)]] <- .elambda.a  #
    if (linfla > 0)
      misc$earg[[(iptr <- iptr + 1)]] <- .epstr.i  #
    if (linfla > 1)
      misc$earg[[(iptr <- iptr + 1)]] <- .elambda.i  #
  }), list(
    .llambda.p = llambda.p, .llambda.i = llambda.i, .lpstr.i = lpstr.i,
    .elambda.p = elambda.p, .elambda.i = elambda.i, .epstr.i = epstr.i,
                            .llambda.a = llambda.a, .lpobs.a = lpobs.a,
                            .elambda.a = elambda.a, .epobs.a = epobs.a,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .alter = alter, .inflate = inflate ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    lalter <- length(as.vector( .alter ))
    linfla <- length(as.vector( .inflate ))
    if (lalter + linfla == 0)
      eta <- as.matrix(eta)
    iptr <- 1  # Unpack
    lambda.p <- eta2theta(eta[, iptr], .llambda.p , earg = .elambda.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a   , earg = .epobs.a )
    lambda.a <- if (lalter <= 1) lambda.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .llambda.a , earg = .elambda.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i   , earg = .epstr.i )
    lambda.i <- if (linfla <= 1) lambda.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .llambda.i , earg = .elambda.i )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitpois(y, lambda.p, lambda.a = lambda.a, lambda.i = lambda.i,
                  pobs.mix.a = pobs.a, pstr.mix.i = pstr.i,
                  truncate = .truncate , max.support = .max.support ,
                  alter.mix = .alter , inflate.mix = .inflate , log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(
    .llambda.p = llambda.p, .llambda.i = llambda.i, .lpstr.i = lpstr.i,
    .elambda.p = elambda.p, .elambda.i = elambda.i, .epstr.i = epstr.i,
                            .llambda.a = llambda.a, .lpobs.a = lpobs.a,
                            .elambda.a = elambda.a, .epobs.a = epobs.a,
    .truncate = truncate, .max.support = max.support, 
    .alter = alter, .inflate = inflate ))),
  vfamily = c("gaitpoisson.mix"),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    lalter <- length(as.vector( .alter ))
    linfla <- length(as.vector( .inflate ))
    if (lalter + linfla == 0)
      eta <- as.matrix(eta)
    iptr <- 1  # Unpack
    lambda.p <- eta2theta(eta[, iptr], .llambda.p , earg = .elambda.p )
    pobs.a   <- if (lalter <= 0) rep(0.25, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a   , earg = .epobs.a )
    lambda.a <- if (lalter <= 1) lambda.p else
    eta2theta(eta[, (iptr <- iptr + 1)], .llambda.a , earg = .elambda.a )
    pstr.i   <- if (linfla <= 0) rep(0.25, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i   , earg = .epstr.i )
    lambda.i <- if (linfla <= 1) lambda.p else
    eta2theta(eta[, (iptr <- iptr + 1)], .llambda.i , earg = .elambda.i )

    okay1 <- all(is.finite(lambda.p)) && all(0 < lambda.p) &&
             all(is.finite(lambda.a)) && all(0 < lambda.a) &&
             all(is.finite(lambda.i)) && all(0 < lambda.i) &&
             all(is.finite(pobs.a)) && all(0 < pobs.a & pobs.a < 1) &&
             all(is.finite(pstr.i)) && all(0 < pstr.i & pstr.i < 1) &&
             all(pobs.a + pstr.i < 1)
    okay1
  }, list(
    .llambda.p = llambda.p, .llambda.i = llambda.i, .lpstr.i = lpstr.i,
    .elambda.p = elambda.p, .elambda.i = elambda.i, .epstr.i = epstr.i,
                            .llambda.a = llambda.a, .lpobs.a = lpobs.a,
                            .elambda.a = elambda.a, .epobs.a = epobs.a,
    .truncate = truncate, .max.support = max.support, 
    .alter = alter, .inflate = inflate ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    lalter <- length(as.vector( .alter ))
    linfla <- length(as.vector( .inflate ))
    iptr <- 1  # Unpack
    lambda.p <- eta2theta(eta[, iptr], .llambda.p , earg = .elambda.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a   , earg = .epobs.a )
    lambda.a <- if (lalter <= 1) lambda.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .llambda.a , earg = .elambda.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i   , earg = .epstr.i )
    lambda.i <- if (linfla <= 1) lambda.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .llambda.i , earg = .elambda.i )
    rgaitpois(nsim * length(lambda.p), lambda.p = lambda.p,
              lambda.a = lambda.a, lambda.i = lambda.i,
              pobs.mix.a = pobs.a, pstr.mix.i = pstr.i,
              alter.mix = .alter , inflate.mix = .inflate ,
              truncate = .truncate , max.support = .max.support )
  }, list(
    .llambda.p = llambda.p, .llambda.i = llambda.i, .lpstr.i = lpstr.i,
    .elambda.p = elambda.p, .elambda.i = elambda.i, .epstr.i = epstr.i,
                            .llambda.a = llambda.a, .lpobs.a = lpobs.a,
                            .elambda.a = elambda.a, .epobs.a = epobs.a,
    .truncate = truncate, .max.support = max.support,
    .alter = alter, .inflate = inflate ))),
  deriv = eval(substitute(expression({
    alter <- as.vector( .alter )
    lalter <- length(alter)
    inflate <- as.vector( .inflate )
    linfla <- length(inflate)
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    M1 <- sum( .tmp3.TF )
    NOS <- NCOL(eta) / M1  # extra$NOS
    max.support <- as.vector( .max.support )

    ans2 <- matrix(0, n, M1)
    Offset.i <- if (lalter == 1) 2 else if (lalter > 1) 3 else 1


    is.altered <- if (lalter)
      rowSums(extra$skip.a) > 0 else rep(FALSE, n)
    is.inflated <- if (linfla)
      rowSums(extra$skip.i) > 0 else rep(FALSE, n)


    if (lalter + linfla == 0)
      eta <- as.matrix(eta)
    iptr <- 1  # Unpack
    lambda.p <- eta2theta(eta[, iptr], .llambda.p , earg = .elambda.p )
    pobs.a   <- if (lalter <= 0) rep(0, NROW(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpobs.a   , earg = .epobs.a )
    lambda.a <- if (lalter <= 1) lambda.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .llambda.a , earg = .elambda.a )
    pstr.i   <- if (linfla <= 0) numeric(nrow(eta)) else
      eta2theta(eta[, (iptr <- iptr + 1)], .lpstr.i   , earg = .epstr.i )
    lambda.i <- if (linfla <= 1) lambda.p else
      eta2theta(eta[, (iptr <- iptr + 1)], .llambda.i , earg = .elambda.i )

    bits <- moments.pois.gait(lambda.p, alter = alter,
                              inflate = inflate, truncate = truncate,
                              max.support = max.support,
                              pobs.a = pobs.a, pstr.i = pstr.i,
                              lambda.a = lambda.a, lambda.i = lambda.i,
                              mlm = FALSE)
    Denom <- c(bits[["lhs.prob"]] - bits[["suma.p"]] - bits[["sumt.p"]])

    pmf.deriv1 <- function(y, lambda)
      dpois(y-1, lambda) - dpois(y, lambda)
    pmf.deriv2 <- function(y, lambda)
      dpois(y-2, lambda) - 2 * dpois(y-1, lambda) + dpois(y, lambda)

    sumderiv1a.p <- sumderiv2a.p <- matrix(0, n, NOS)
    if (lalter > 0) {
      deriv0matrix.a <-
      deriv1matrix.a <- deriv2matrix.a <- matrix(0, n, lalter)
      sumderiv1a.a <- sumderiv2a.a <- matrix(0, n, NOS)
      for (jay in seq(lalter)) {
        aval <- alter[jay]
        sumderiv1a.p <- sumderiv1a.p + pmf.deriv1(aval, lambda.p)
        sumderiv2a.p <- sumderiv2a.p + pmf.deriv2(aval, lambda.p)
        sumderiv1a.a <- sumderiv1a.a + pmf.deriv1(aval, lambda.a)
        sumderiv2a.a <- sumderiv2a.a + pmf.deriv2(aval, lambda.a)
        pmf.a <- dpois(aval, lambda.a)
        deriv0matrix.a[, jay] <- pmf.a
        deriv1matrix.a[, jay] <- pmf.deriv1(aval, lambda.a)
        deriv2matrix.a[, jay] <- pmf.deriv2(aval, lambda.a)
      }
    }  # lalter > 0





    if (linfla) {
      deriv0matrix.i <-  # wrt inflated distribution
      deriv1matrix.i <- deriv2matrix.i <- matrix(0, n, linfla)
      Deriv0Matrix.i <-  # wrt parent distribution
      Deriv1Matrix.i <- Deriv2Matrix.i <- matrix(0, n, linfla)
      if (linfla > 0)
        for (jay in seq(linfla)) {
          ival <- inflate[jay]
          pmf.i <- dpois(ival, lambda.i)
          deriv0matrix.i[, jay] <- pmf.i
          deriv1matrix.i[, jay] <- pmf.deriv1(ival, lambda.i)
          deriv2matrix.i[, jay] <- pmf.deriv2(ival, lambda.i)
          pmf.p <- dpois(ival, lambda.p)
          Deriv0Matrix.i[, jay] <- pmf.p
          Deriv1Matrix.i[, jay] <- pmf.deriv1(ival, lambda.p)
          Deriv2Matrix.i[, jay] <- pmf.deriv2(ival, lambda.p)  # /pmf.p
        }  # jay
    }  # linfla




    sumderiv1t.a <- sumderiv2t.a <-
    sumderiv1t.i <- sumderiv2t.i <-
    sumderiv1t.p <- sumderiv2t.p <- matrix(0, n, NOS)
    if (ltrunc)
      for (tval in truncate) {
        sumderiv1t.p <- sumderiv1t.p + pmf.deriv1(tval, lambda.p)
        sumderiv2t.p <- sumderiv2t.p + pmf.deriv2(tval, lambda.p)
        sumderiv1t.a <- sumderiv1t.a + pmf.deriv1(tval, lambda.a)
        sumderiv2t.a <- sumderiv2t.a + pmf.deriv2(tval, lambda.a)
        sumderiv1t.i <- sumderiv1t.i + pmf.deriv1(tval, lambda.i)
        sumderiv2t.i <- sumderiv2t.i + pmf.deriv2(tval, lambda.i)
      }

    sumderiv1t.p <- sumderiv1t.p + dpois( .max.support    , lambda.p)
    sumderiv2t.p <- sumderiv2t.p + dpois( .max.support - 1, lambda.p) -
                                   dpois( .max.support    , lambda.p)
    sumderiv1t.a <- sumderiv1t.a + dpois( .max.support    , lambda.a)
    sumderiv2t.a <- sumderiv2t.a + dpois( .max.support - 1, lambda.a) -
                                   dpois( .max.support    , lambda.a)
    sumderiv1t.i <- sumderiv1t.i + dpois( .max.support    , lambda.i)
    sumderiv2t.i <- sumderiv2t.i + dpois( .max.support - 1, lambda.i) -
                                   dpois( .max.support    , lambda.i)

    onembothprobs <- 1 - pobs.a - pstr.i  # 1-Pr(altered)-Pr(inflated)

    zero0n <- rep(0, n)
    dl.dlambda.a <-
    dl.dlambda.i <- zero0n  # Replace some elts below
    dl.dlambda.p <- y / lambda.p - 1 +
      (sumderiv1a.p + sumderiv1t.p) / Denom  # \notin A, I, T
    dl.dlambda.p[is.altered] <- 0
    dl.dpobs.a <-  # Replace some elts below
    dl.dpstr.i <- -1 / onembothprobs  # \notin A, I, T
    Delta <- matrix(0, n, NOS)  # If linfla == 0

    if (linfla > 0) {
      d0A.io <- deriv0matrix.i / c(bits[["sumi.i"]])
      d0B.pi <- Deriv0Matrix.i / Denom
      Delta <- pstr.i * d0A.io + onembothprobs * d0B.pi


      d1A.io <- deriv1matrix.i / c(bits[["sumi.i"]]) -
                deriv0matrix.i * rowSums(deriv1matrix.i) / c(
                    bits[["sumi.i"]])^2
      d2A.io <- deriv2matrix.i / c(bits[["sumi.i"]]) -
                2 * deriv1matrix.i * rowSums(deriv1matrix.i) / c(
                    bits[["sumi.i"]])^2 -
                deriv0matrix.i * rowSums(deriv2matrix.i) / c(
                    bits[["sumi.i"]])^2 +
                2 * deriv0matrix.i * rowSums(deriv1matrix.i)^2 / c(
                    bits[["sumi.i"]])^3
      d1B.pi <- Deriv1Matrix.i / Denom +
        Deriv0Matrix.i * c(sumderiv1t.p + sumderiv1a.p) / Denom^2
      d2B.pi <- Deriv2Matrix.i / Denom +
        2 * Deriv1Matrix.i * c(sumderiv1t.p + sumderiv1a.p) / Denom^2 +
        Deriv0Matrix.i * c(sumderiv2t.p + sumderiv2a.p) / Denom^2 +
        2 * Deriv0Matrix.i * (c(sumderiv1t.p + sumderiv1a.p)^2) / Denom^3
    }  # linfla > 0
    
    if (lalter) {
      dl.dpobs.a[is.altered] <- 1 / pobs.a[is.altered]
      dl.dpstr.i[is.altered] <- 0

      d0A.al <- deriv0matrix.a / c(bits[["suma.a"]])
      d1A.al <- deriv1matrix.a / c(bits[["suma.a"]]) -
                deriv0matrix.a * rowSums(deriv1matrix.a) / c(
                    bits[["suma.a"]])^2
      d2A.al <- deriv2matrix.a / c(bits[["suma.a"]]) -
                2 * deriv1matrix.a * rowSums(deriv1matrix.a) / c(
                    bits[["suma.a"]])^2 -
                deriv0matrix.a * rowSums(deriv2matrix.a) / c(
                    bits[["suma.a"]])^2 +
                2 * deriv0matrix.a * rowSums(deriv1matrix.a)^2 / c(
                    bits[["suma.a"]])^3

      for (jay in seq(lalter)) {
        aval <- alter[jay]
        is.alte.j <- extra$skip.a[, jay]  # Logical vector
        tmp2 <- d1A.al[, jay] / d0A.al[, jay]
        dl.dlambda.a[is.alte.j] <- tmp2[is.alte.j]
      }  # jay
    }  # lalter

    if (linfla > 0)
      for (jay in seq(linfla)) {
        ival <- inflate[jay]
        is.infl.j <- extra$skip.i[, jay]  # Logical vector
        tmp7 <- onembothprobs * d1B.pi[, jay] / Delta[, jay]
        dl.dlambda.p[is.infl.j] <- tmp7[is.infl.j]
        tmp2 <- -d0B.pi[, jay] / Delta[, jay]
        dl.dpobs.a[is.infl.j] <- tmp2[is.infl.j]
        tmp8 <- (d0A.io[, jay] - d0B.pi[, jay]) / Delta[, jay]
        dl.dpstr.i[is.infl.j] <- tmp8[is.infl.j]
        if (linfla > 1) {
          tmp2 <- pstr.i * d1A.io[, jay] / Delta[, jay]
          dl.dlambda.i[is.infl.j] <- tmp2[is.infl.j]
        }
      }  # jay


    dlambda.p.deta <- dtheta.deta(lambda.p, .llambda.p , .elambda.p )
    ans2[, 1] <- dl.dlambda.p * dlambda.p.deta
    dpobs.a.deta <- dtheta.deta(pobs.a, .lpobs.a , .epobs.a )
    if (lalter > 0) {
      ans2[, 2] <- dl.dpobs.a * dpobs.a.deta
    }
    if (lalter > 1) {
      dlambda.a.deta <- dtheta.deta(lambda.a, .llambda.a , .elambda.a )
      ans2[, 3] <- dl.dlambda.a * dlambda.a.deta
    }
    if (linfla > 0) {
      dpstr.i.deta <- dtheta.deta(pstr.i, .lpstr.i , .epstr.i )
      ans2[, Offset.i + 1] <- dl.dpstr.i * dpstr.i.deta
    }
    if (linfla > 1) {
      dlambda.i.deta <- dtheta.deta(lambda.i, .llambda.i , .elambda.i )
      ans2[, Offset.i + 2] <- dl.dlambda.i * dlambda.i.deta
    }
    c(w) * ans2
  }), list(
    .llambda.p = llambda.p, .llambda.i = llambda.i, .lpstr.i = lpstr.i,
    .elambda.p = elambda.p, .elambda.i = elambda.i, .epstr.i = epstr.i,
                            .llambda.a = llambda.a, .lpobs.a = lpobs.a,
                            .elambda.a = elambda.a, .epobs.a = epobs.a,
    .tmp3 = tmp3, .tmp3.TF = tmp3.TF, 
    .truncate = truncate, .max.support = max.support, 
    .alter = alter , .inflate = inflate ))),

  weight = eval(substitute(expression({
    cond.EY.p <- (lambda.p - bits[["SumA.p"]] - bits[["SumI.p"]] -
                  bits[["SumT.p"]]) / c(
       bits[["lhs.prob"]] - bits[["suma.p"]] - bits[["sumi.p"]] -
       bits[["sumt.p"]])
    onemrsDelta <- onembothprobs *  # (onempstr.i - pobs.a) *
                   (1 - c(bits[["sumi.p"]]) / Denom)

    ned2l.dlambda.p.lambda.a <- zero0n  # Final; nothing to do
    ned2l.dpobs.a.lambda.p   <- zero0n  # mB overwritten below [4279]
    ned2l.dpobs.a.lambda.a   <- zero0n  # Final; nothing to do
    ned2l.dpobs.a.pstr.i     <- zero0n  # mB overwritten below
    ned2l.dpobs.a.lambda.i   <- zero0n  # mB overwritten below
    ned2l.dpstr.i.lambda.p   <- zero0n  # mB overwritten below
    ned2l.dpstr.i.lambda.a   <- zero0n  # Final; nothing to do
    ned2l.dpstr.i.lambda.i   <- zero0n  # mB overwritten below
    ned2l.dlambda.a.lambda.i <- zero0n  # Final; nothing to do

    ned2l.dlambda.p2 <- onemrsDelta * (cond.EY.p / lambda.p^2 -
                         c(sumderiv2t.p + sumderiv2a.p) / Denom -
                        (c(sumderiv1t.p + sumderiv1a.p) / Denom)^2)
    if (linfla > 0)
      ned2l.dlambda.p2 <- ned2l.dlambda.p2 + onembothprobs *
        rowSums(onembothprobs * d1B.pi^2 / Delta - d2B.pi)

    if (lalter > 0) {
      wz22 <- if ( .lpobs.a == "logitlink" && linfla == 0) {
        pobs.a * (1 - pobs.a)
      } else {
        ned2l.dpobs.a2 <- 1 / pobs.a +
          onemrsDelta / onembothprobs^2 + (if (linfla > 0)
          rowSums(d0B.pi^2 / Delta) else 0)
        ned2l.dpobs.a2 * dpobs.a.deta^2
      }
    }
    if (lalter > 1)
      ned2l.dlambda.a2 <- pobs.a * (
        rowSums(deriv1matrix.a^2  / deriv0matrix.a) / c(
        bits[["suma.a"]]) -
        (c(sumderiv1a.a) / c(bits[["suma.a"]]))^2)

    if (linfla > 0)
      ned2l.dpstr.i2 <-
        onemrsDelta / onembothprobs^2 +
        (if (linfla > 0)  # zz > 1 ??
        rowSums((deriv0matrix.i / c(bits[["sumi.i"]]) -
                 Deriv0Matrix.i / Denom)^2 / Delta) else 0)

    if (linfla > 0)
      ned2l.dpstr.i.lambda.p <-
        rowSums(d1B.pi * (1 + onembothprobs *
                         (d0A.io - d0B.pi) / Delta))

    if (linfla > 1)
      ned2l.dlambda.i2 <- pstr.i *
        rowSums(pstr.i * (d1A.io^2) / Delta - d2A.io)

    
    if (linfla > 1)
      ned2l.dpstr.i.lambda.i <-
        rowSums(d1A.io * (pstr.i * (d0A.io - d0B.pi) / Delta - 1))


    if (linfla > 1)
      ned2l.dlambda.p.lambda.i <- pstr.i * onembothprobs *
                                  rowSums(d1A.io * d1B.pi / Delta)


    if (lalter > 0 && linfla > 0)
      ned2l.dpobs.a.pstr.i <- onemrsDelta / onembothprobs^2 -
        rowSums(d0B.pi * (d0A.io - d0B.pi) / Delta)


    if (linfla > 1)
      ned2l.dpobs.a.lambda.i <-
        rowSums(-pstr.i * d1A.io * d0B.pi / Delta)

    
    if (linfla > 0)
      ned2l.dpobs.a.lambda.p <-
        rowSums(d1B.pi * (1 - onembothprobs * d0B.pi / Delta))
   

    wz <- matrix(0, n, if (lalter && !linfla) M1 else M1*(M1+1)/2)

    wz[, iam(1, 1, M)] <- ned2l.dlambda.p2 * dlambda.p.deta^2

    if (linfla > 0)
      wz[, iam(1, 2, M)] <-
        ned2l.dpobs.a.lambda.p * dpobs.a.deta * dlambda.p.deta
    if (lalter > 0)
      wz[, iam(2, 2, M)] <- wz22
    if (lalter > 1)
      wz[, iam(3, 3, M)] <- ned2l.dlambda.a2 * dlambda.a.deta^2

    if (linfla > 0) {
      wz[, iam(Offset.i + 1, Offset.i + 1, M)] <-
        ned2l.dpstr.i2 * dpstr.i.deta^2
      wz[, iam(Offset.i + 1, 1, M)] <-
        ned2l.dpstr.i.lambda.p * dpstr.i.deta * dlambda.p.deta
      if (lalter > 0)
        wz[, iam(Offset.i + 1, 2, M)] <-
          ned2l.dpobs.a.pstr.i * dpobs.a.deta * dpstr.i.deta
    }

    if (linfla > 1) {
      wz[, iam(Offset.i + 2, Offset.i + 2, M)] <-
        ned2l.dlambda.i2 * dlambda.i.deta^2
      wz[, iam(Offset.i + 2, 1, M)] <-
        ned2l.dlambda.p.lambda.i * dlambda.p.deta * dlambda.i.deta
      wz[, iam(Offset.i + 2, 2, M)] <-
        ned2l.dpobs.a.lambda.i * dpobs.a.deta * dlambda.i.deta
      wz[, iam(Offset.i + 2, Offset.i + 1, M)] <-
        ned2l.dpstr.i.lambda.i * dpstr.i.deta * dlambda.i.deta
    }

    c(w) * wz
  }), list( .lpobs.a = lpobs.a, .lpstr.i = lpstr.i))))
}  # gaitpoisson.mix







 moments.log.gait <-
  function(shape.p,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           byrow.ai = FALSE,  # For pobs.a and pstr.i
           shape.a = shape.p, shape.i = shape.p,
           mlm = TRUE,  # mlm = FALSE iff mixture = TRUE
           type.fitted = "All") {  # or "mean"
  A8.p <- -1 / log1p(-shape.p)
  rmlife <- A8.p * (shape.p^(max.support + 1)) / (1 - shape.p)
  NOS <- 1
  nnn <- length(shape.p)
  lhs.prob <- plog(max.support, shape.p)  # 1 - Pr(upper tail)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)

  sumt.p <- matrix(0, nnn, NOS)  # Does not include upper RHS tail
  SumT.p <- matrix(rmlife, nnn, NOS)  # Does include upper RHS tail
  if (ltrunc)
    for (tval in truncate) {
      pmf.p <- dlog(tval, shape.p)
      sumt.p <- sumt.p + pmf.p  # Need tval <= max.support
      SumT.p <- SumT.p + pmf.p * tval
    }

  use.pobs.a <- use.pstr.i <- matrix(0, nnn, 1)  # So rowSums() works.
  aprd <- 0  # aprd is an innerprod
  SumA.x <-  # For innerprod (aprd) only
  SumA.a <- suma.a <-
  SumA.p <- suma.p <- matrix(0, nnn, NOS)
  if (mlm) {  # use.pobs.a and use.pstr.i have the right dimensions.
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
      pmf.p <- dlog(aval, shape.p)
      pmf.a <- dlog(aval, shape.a)
      suma.p <- suma.p + pmf.p
      SumA.p <- SumA.p + pmf.p * aval
      suma.a <- suma.a + pmf.a
      SumA.a <- SumA.a + pmf.a * aval
      SumA.x <- SumA.x + pmf.x * aval
    }  # for jay
    aprd <- if (mlm) SumA.x else use.pobs.a * SumA.a / suma.a
  }  # lalter

  iprd <- 0  # iprd is an innerprod
  SumI.x <-  # For innerprod (iprd) only
  sumi.i <- SumI.i <-
  sumi.p <- SumI.p <- matrix(0, nnn, NOS)
  if (linfla) {
    for (jay in seq_len(linfla)) {
      ival <- inflate[jay]
      pmf.x <- if (mlm) use.pstr.i[, jay] else rep(0, nnn)
      pmf.p <- dlog(ival, shape.p)  # Changed
      pmf.i <- dlog(ival, shape.i)  # Changed
      sumi.p <- sumi.p + pmf.p
      SumI.p <- SumI.p + pmf.p * ival
      sumi.i <- sumi.i + pmf.i
      SumI.i <- SumI.i + pmf.i * ival
      SumI.x <- SumI.x + pmf.x * ival
    }  # for jay
    iprd <- if (mlm) SumI.x else use.pstr.i * SumI.i / sumi.i
  }  # linfla

  use.this <- if (mlm)
    (1 - rowSums(use.pobs.a) - rowSums(use.pstr.i)) else
    (1 - use.pobs.a - use.pstr.i)

  mean.true.p <- A8.p * shape.p / (1 - shape.p)
  themean <- aprd + iprd + use.this *
    (mean.true.p - SumA.p - SumT.p) / (lhs.prob - suma.p - sumt.p)
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
}  # moments.log.gait






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
  lhs.prob <- ppois(max.support, lambda.p)  # 1 - Pr(upper tail)
  lalter <- length(alter)
  linfla <- length(inflate)
  ltrunc <- length(truncate)

  sumt.p <- matrix(0, nnn, NOS)  # Does not include upper RHS tail
  SumT.p <- matrix(rmlife, nnn, NOS)  # Does include upper RHS tail
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
  if (mlm) {  # use.pobs.a and use.pstr.i have the right dimensions.
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
  SumI.x <-  # For innerprod (iprd) only
  sumi.i <- SumI.i <-
  sumi.p <- SumI.p <- matrix(0, nnn, NOS)
  if (linfla) {
    for (jay in seq_len(linfla)) {
      ival <- inflate[jay]
      pmf.x <- if (mlm) use.pstr.i[, jay] else rep(0, nnn)
      pmf.p <- dpois(ival, lambda.p)  # Changed
      pmf.i <- dpois(ival, lambda.i)  # Changed
      sumi.p <- sumi.p + pmf.p
      SumI.p <- SumI.p + pmf.p * ival
      sumi.i <- sumi.i + pmf.i
      SumI.i <- SumI.i + pmf.i * ival
      SumI.x <- SumI.x + pmf.x * ival
    }  # for jay
    iprd <- if (mlm) SumI.x else use.pstr.i * SumI.i / sumi.i
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
  semigait.errorcheck(alter, inflate, truncate, max.support)
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

 warning("this is based on an old algorithm")



      semigait.errorcheck(alter, inflate, truncate, max.support)
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







 pgaitnbinom.mix <-
  function(q, size.p, prob.p = NULL, munb.p = NULL,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           munb.a = munb.p, munb.i = munb.p) {
  semigait.errorcheck(alter, inflate, truncate, max.support)
  if ((is.prob <- as.logical(length(prob.p)) &&
      length(munb.p)))
    stop("cannot specify both 'prob.p' and 'munb.p' arguments")
  if (!length(c(alter, inflate, truncate)) &&
      is.infinite(max.support)) {
    return(if (is.prob)
      pnbinom(q, size = size.p, prob = prob.p) else
      pnbinom(q, size = size.p, mu   = munb.p))
  }

  LLL <- max(length(q),
             length(size.p), length(prob.p), length(munb.p),
             length(size.a), length(prob.a), length(munb.a),
             length(size.i), length(prob.i), length(munb.i),
             length(pobs.a), length(pstr.i))
  use.pobs.a <- 0
  if (length(alter)) {
    use.pobs.a <- matrix(0, LLL, length(alter))
    for (jay in seq(length(alter))) {
      aval <- alter[jay]
      use.pobs.a[, jay] <- if (length(prob.a))
        dnbinom(aval, size.a, prob = prob.a) else
        dnbinom(aval, size.a, mu   = munb.a)
    }
    use.pobs.a <- pobs.a * use.pobs.a / rowSums(use.pobs.a)
  }

  use.pstr.i <- 0
  if (length(inflate)) {
    use.pstr.i <- matrix(0, LLL, length(inflate))
    for (jay in seq(length(inflate))) {
      ival <- inflate[jay]
      use.pstr.i[, jay] <- if (length(prob.a))
        dnbinom(ival, size.i, prob = prob.i) else
        dnbinom(ival, size.i, mu   = munb.i)
    }
    use.pstr.i <- pstr.i * use.pstr.i / rowSums(use.pstr.i)
  }

  cdf1 <-
    pgaitnbinom.mlm(q, size.p,  # Inner distribution
                    prob = prob.p, munb = munb.p,
                    alter = alter,
                    inflate = inflate,
                    truncate = truncate,
                    max.support = max.support,
                    pobs.a = use.pobs.a,  # byrow.arg = FALSE,
                    pstr.i = use.pstr.i)
  cdf1
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
  semigait.errorcheck(alter, inflate, truncate, max.support)
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

  min.support <- 0  # Usual case; same as lowsup
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















 gatnbinomial.mlm <-
  function(alter = NULL,
           truncate = NULL,  # max.support = Inf,
           zero = "size",
           lmunb = "loglink", lsize = "loglink",
           type.fitted = c("mean", "pobs.a", "Pobs.a",
                       "prob.a", "prob.i", "prob.t", "lhs.prob"),
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

  semigait.errorcheck(alter, inflate = NULL, truncate, max.support)

  lalter <- length(alter)
  ltrunc <- length(truncate)
  type.fitted <- match.arg(type.fitted, c("mean", "pobs.a", "Pobs.a",
                       "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]
  temp7 <- if (lalter) paste0("pobs", alter) else NULL
  tmp3 <- c(munb = lmunb, size = lsize,
            if (lalter) rep("multilogitlink", lalter) else NULL)
  names(tmp3) <- c("munb", "size", temp7)
       
  blurb1 <- "N"
  if (lalter) blurb1 <- "Generally-altered n"
  if (ltrunc) blurb1 <- "Generally-truncated n"
  if (lalter && ltrunc) blurb1 <- "Generally-altered and -truncated n"
       
  new("vglmff",
  blurb = c(blurb1, "egative binomial regression\n",
            "(GAT-NB-MLM generally)\n\n",
            "Links:   ",
            namesof("munb", lmunb, earg = emunb, tag = FALSE), ", ",
            namesof("size", lsize, earg = esize, tag = FALSE),
            if (lalter)
            paste(",\n         multilogitlink(cbind(",
            paste(temp7, collapse = ", "), ", ",
            "\n                              ",
            "1-", paste(temp7, collapse = "-"),
            "))", sep = "")
      ),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2 + length( .alter ))
  }), list( .zero = zero,
            .alter = alter ))),

  infos = eval(substitute(function(...) {
    alter <- as.vector( .alter )
    temp7 <- paste0("pobs", alter)
    list(M1 = 2 + length( .alter ),
         Q1 = 1,
         link = .tmp3 ,  # multilogitlink is multiparameter:
         link1parameter = as.logical(length( .alter ) <= 2),
         mixture.links = as.logical(length( .alter ) > 2),
         alter = as.vector( .alter ),
         truncate = as.vector( .truncate ),
         max.support = as.vector( .max.support ), 
         expected = TRUE,
         multipleResponses = !as.logical(length( .alter ) +
                                         length( .truncate )),
         parameters.names = c("munb", "size",
                              if (length( .alter )) temp7 else NULL),
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
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    alter <- as.vector( .alter )
    lalter <- length(alter)
    M1 <- 2 + lalter
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


    mynames3 <- if (lalter) {
      temp7 <- paste0("pobs", alter)
      denom.char <- paste0("1-", paste0(temp7, collapse = "-"))
      paste("log(", temp7, "/(", denom.char, "))", sep = "")
    } else {
      NULL
    }
    mynames1 <- param.names("munb", NOS, skip1 = TRUE)
    mynames2 <- param.names("size", NOS, skip1 = TRUE)
    predictors.names <-
      c(namesof(mynames1, .lmunb , earg = .emunb , tag = FALSE),
        namesof(mynames2, .lsize , earg = .esize , tag = FALSE),
          mynames3)[
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
        cbind(theta2eta(munb.init, .lmunb , earg = .emunb ),
              theta2eta(size.init, .lsize , earg = .esize ),
              etastart)  # [interleave.VGAM(M, M1 = M1)]
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
                     c("mean", "pobs.a", "Pobs.a",
                       "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]
    alter <- as.vector( .alter )
    lalter <- length(alter)
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    max.support <- ( .max.support )
    M1 <- 2 + length(alter)
    if ((NOS <- NCOL(eta) / M1) != 1) stop("can only handle NOS == 1")
    munb.p <- eta2theta(eta[, 1, drop = FALSE], .lmunb , earg = .emunb )
    size.p <- eta2theta(eta[, 2, drop = FALSE], .lsize , earg = .esize )
    if (lalter) {
      index3 <- 3:NCOL(eta)
      omegamat <- multilogitlink(eta[, index3, drop = FALSE],
                  refLevel = NCOL(eta) - 1,  # Assumes 1 response
                  inverse = TRUE)
      ynames.Pobs.a <- c(as.character(alter), "(Others)")
      dimnames(omegamat) <- list(rownames(eta), ynames.Pobs.a)
    }

    Bits <- moments.nbin.gait(size.p = size.p, munb.p = munb.p,
                              alter = alter,
                              truncate = truncate, mlm = TRUE,
                              max.support = max.support,
                              pobs.a = omegamat[, -ncol(omegamat)])

    ans <- switch(type.fitted,
      "mean"       = Bits[["mean"]],
      "pobs.a"     = omegamat,
      "Pobs.a"     = omegamat,  #[, -ncol(omegamat)],  # Same as "pobs.a"
      "prob.a"     = Bits[["suma.p"]],  # suma
      "prob.i"     = Bits[["sumi.p"]],
      "prob.t"     = Bits[["sumt.p"]],  # Truncated probs w/o RHS
      "lhs.prob"   = Bits[["lhs.prob"]])  # 1 - Pr(y <= max.support)
    label.cols.y(ans,
        colnames.y = switch(type.fitted,
                            "pobs.a" =,
                            "Pobs.a" = ynames.Pobs.a,
                            extra$colnames.y),
        NOS = NOS)
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize,
           .alter = alter,
           .truncate = truncate, .max.support = max.support ))),

  last = eval(substitute(expression({
    tmp9 <- if (lalter) rep_len("multilogitlink", lalter) else NULL
    temp.names <- c( .lmunb , .lsize , tmp9)
    misc$link  <- temp.names
    names(misc$link) <- c(mynames1, mynames2, mynames3)
    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)

    misc$earg[[1]] <- .emunb  #
    misc$earg[[2]] <- .esize  #
    if (lalter) {
      for (ii in seq(M1 - 2)) {
        misc$earg[[2 + ii]] <- list(M = M - 2,  # M * NOS,
                                    refLevel = M - 1)  # M * NOS
      }  # ii
    }  # lalter
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize,
            .alter = alter ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    lalter <- length( .alter )
    index3 <- 3:NCOL(eta)
    if (length( .alter ))
      pobs.a <- multilogitlink(eta[, index3, drop = FALSE],
                  refLevel = NCOL(eta) - 1,  # Assumes 1 response
                  inverse = TRUE)
    munb <- eta2theta(eta[, 1, drop = FALSE], .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, 2, drop = FALSE], .lsize , earg = .esize )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitnbinom.mlm(y, kmat, munb = munb, log = TRUE,
                        truncate = .truncate ,
                        max.support = .max.support ,
                        alter = .alter ,  # byrow.arg = FALSE,
                        pobs.a = if (lalter)
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
    index3 <- 3:NCOL(eta)
    pobs.a <- if (length( .alter ))
        multilogitlink(eta[, index3, drop = FALSE],
                refLevel = NCOL(eta) - 1,  # Assumes one response
                inverse = TRUE) else 0  # An okay value
    munb <- eta2theta(eta[, 1, drop = FALSE], .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, 2, drop = FALSE], .lsize , earg = .esize )
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

    index3 <- 3:NCOL(eta)
    pobs.a <- if (length( .alter ))
        multilogitlink(eta[, index3, drop = FALSE],
                refLevel = NCOL(eta) - 1,  # Assumes one response
                inverse = TRUE) else 0  # An okay value
    munb <- eta2theta(eta[, 1, drop = FALSE], .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, 2, drop = FALSE], .lsize , earg = .esize )
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
    alter <- as.vector( .alter )
    lalter <- length(alter)
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    M1 <- 2 + lalter
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

    index3 <- 3:NCOL(eta)
    phimat <-
    pobs.a <- if (length( .alter ))
        multilogitlink(eta[, index3, drop = FALSE],
                refLevel = NCOL(eta) - 1,  # Assumes 1 response
                inverse = TRUE) else 0  # An okay value
    munb <- eta2theta(eta[, 1, drop = FALSE], .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, 2, drop = FALSE], .lsize , earg = .esize )
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
    ans <- c(w) * cbind(dl.dmunb * dmunb.deta,
                        dl.dsize * dsize.deta,
                        dl.deta)
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
              phimat[sind2, -NCOL(phimat), drop = FALSE] else
              rbind(pobs.a),  # rbind added for when pobs.a = 0
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
        wz.merge(c(w) * cbind(ned2l.dmunb2 * dmunb.deta^2,
                              ned2l.dsize2 * dsize.deta^2,
                              ned2l.dmunbsize * dmunb.deta *
                              dsize.deta) * c(onempobs.a),
                 c(w) * wz4,
                 M1 = 2, M2 = lalter)  # rm.trailing.cols = F
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






 dgaitbinom.mlm <-
  function(x, size, prob,
           alter = NULL,
           inflate = NULL,
           truncate = NULL,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE,
           log.arg = FALSE,
           .errorcheck = TRUE) {
  if ( .errorcheck )
    semigait.errorcheck(alter, inflate, truncate,
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
    semigait.errorcheck(alter, inflate, truncate,
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
  semigait.errorcheck(alter, inflate, truncate,
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
  blurb = c("Generally-truncated binomial regression (GT-Binom)\n\n",
            "Links:    ",
            if (multiple.responses)
            c(namesof("prob1", link, earg = earg, tag = FALSE),
            ",...,",
            namesof("probM", link, earg = earg, tag = FALSE)) else
            namesof("prob",  link, earg = earg, tag = FALSE),
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
  blurb = c("Generally-altered binomial regression\n",
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
    alter <- as.vector( .alter )
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

    alter <- as.vector( .alter )
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

    alter <- as.vector( .alter )
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

    alter <- as.vector( .alter )
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








 gatnbinomial.mix <-
  function(alter = NULL,  # May need a (>=3)zz-vector to run
           truncate = NULL,
           zero = c("pobs.a", "size"),
           parallel = FALSE,  # TRUE applies to the intercept
           lmunb.p = "loglink",
           lsize.p = "loglink",
           lpobs.a = "logitlink",
           lmunb.a = "loglink",
           lsize.a = "loglink",
           type.fitted = c("mean", "pobs.a", "Pobs.a",
             "prob.a", "prob.i", "prob.t", "lhs.prob"),
           imethod = 1,
           imunb.p = NULL,    isize.p = NULL,  # exp(1),
           ipobs.a = NULL,
           imunb.a = imunb.p, isize.a = isize.p,
           ishrinkage = 0.95,
           probs.y = 0.35,

           cutoff.prob = 0.999,  # higher is better for large 'size'
           eps.trig = 1e-7,
           max.chunk.MB = 30  # max.memory = Inf is allowed
           ) {
  max.support <- Inf  # Fixed for now, as its intractable
  lmunb.p <- as.list(substitute(lmunb.p))
  emunb.p <- link2list(lmunb.p)
  lmunb.p <- attr(emunb.p, "function.name")
  lsize.p <- as.list(substitute(lsize.p))
  esize.p <- link2list(lsize.p)
  lsize.p <- attr(esize.p, "function.name")

  lpobs.a <- as.list(substitute(lpobs.a))
  epobs.a <- link2list(lpobs.a)
  lpobs.a <- attr(epobs.a, "function.name")

  lmunb.a <- as.list(substitute(lmunb.a))
  emunb.a <- link2list(lmunb.a)
  lmunb.a <- attr(emunb.a, "function.name")
  lsize.a <- as.list(substitute(lsize.a))
  esize.a <- link2list(lsize.a)
  lsize.a <- attr(esize.a, "function.name")

  semigait.errorcheck(alter, inflate = NULL, truncate, max.support)

  lalter <- length(alter)
  ltrunc <- length(truncate)
  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs.a", "Pobs.a",
             "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]
  tmp3 <- c(munb.p = lmunb.p,
            size.p = lsize.p,
            pobs.a = lpobs.a,
            munb.a = lmunb.a,
            size.a = lsize.a)[seq(switch(as.character(lalter),
                                         "0" = 2, "1" = 3, 5))]
  blurb1 <- "N"
  if (lalter) blurb1 <- "Generally-altered n"
  if (ltrunc) blurb1 <- "Generally-truncated n"
  if (lalter && ltrunc) blurb1 <- "Generally-altered and -truncated n"
                
  new("vglmff",
  blurb = c(blurb1, "egative binomial regression\n",
            "(GAT-NB-NB mixture generally)\n\n",
            "Links:   ",
            namesof("munb.p", lmunb.p, emunb.p, tag = FALSE), ", ",
            namesof("size.p", lsize.p, esize.p, tag = FALSE),
            if (lalter > 0) c(", ",
            namesof("pobs.a", lpobs.a, epobs.a, tag = FALSE)),
            if (lalter > 1) c(", ",
            namesof("munb.a", lmunb.a, emunb.a, tag = FALSE), ", ",
            "\n         ",
            namesof("size.a", lsize.a, esize.a, tag = FALSE)),
            "\n"),

  constraints = eval(substitute(expression({
    use.mat <- if ( .lalter > 1)
      matrix(c(1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0), 5, 3) else
      if ( .lalter > 0) diag(3) else diag(2)
    if ( .lalter > 1)
    constraints <- cm.VGAM(use.mat, x = x,
                           bool = .parallel ,
                           constraints = constraints,
                           apply.int = TRUE)  # FALSE
    if (lalter > 0)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                   predictors.names = predictors.names,
                   M1 = switch(as.character(lalter), "0"=2, "1"=3, 5))
  }), list( .zero = zero, .lalter = lalter,
            .parallel = parallel ))),

  infos = eval(substitute(function(...) {
    list(M1 =  switch(as.character(length( .alter )),
                      "0" = 2, "1" = 3, 5),
         Q1 = 1,
         link = .tmp3 ,
         link1parameter = TRUE,
         mixture.links = FALSE,
         alter = as.vector( .alter ),
         truncate = as.vector( .truncate ),
         max.support = as.vector( .max.support ), 
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE
         parameters.names = names (.tmp3 ),
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
    if (lalter == 2)
      warning("argument 'alter' should not contain 2 values")
    M1 <- switch(as.character(lalter), "0" = 2, "1" = 3, 5)
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

    if (lalter > 0) {
      y0 <- matrix(0, n, lalter)
      for (jay in seq(lalter))
        y0[, jay] <- as.numeric(y == alter[jay])
      extra$skip.these <- matrix(as.logical(y0), n, lalter)  # dim lost
      if (any((css <- colSums(extra$skip.these)) == 0))
        stop("some 'alter' argument values have no response values: ",
             paste(alter[css == 0], collapse = ", "))          
    }

    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)  # May be NULL
    extra$M1 <- M1

    mynames1 <- param.names("munb.p", ncoly, skip1 = TRUE)
    mynames2 <- param.names("size.p", ncoly, skip1 = TRUE)
    mynames3 <- param.names("pobs.a", ncoly, skip1 = TRUE)
    mynames4 <- param.names("munb.a", ncoly, skip1 = TRUE)
    mynames5 <- param.names("size.a", ncoly, skip1 = TRUE)
    predictors.names <-
      c(namesof(mynames1, .lmunb.p , earg = .emunb.p , tag = FALSE),
        namesof(mynames2, .lsize.p , earg = .esize.p , tag = FALSE),
        if (lalter <= 0) NULL else
        namesof(mynames3, .lpobs.a , earg = .epobs.a , tag = FALSE),
        if (lalter <= 1) NULL else
        namesof(mynames4, .lmunb.a , earg = .emunb.a , tag = FALSE),
        if (lalter <= 1) NULL else
        namesof(mynames5, .lsize.a , earg = .esize.a , tag = FALSE))

    if (!length(etastart)) {
      munb.init <- if (length( .imunb.p )) rep( .imunb.p , n) else
        Init.mu(y = y, w = w, imethod = .imethod ,
                imu = .imunb.p ,  # x = x,
                ishrinkage = .ishrinkage ,
                probs.y = .probs.y )
      try.this <- if (length( .isize.p )) .isize.p else
        0.25 * (0.01 + weighted.mean(y, w)^2) / (0.01 + var(y))
      size.init <- matrix(try.this, n, NCOL(munb.init))
      if (lalter > 0)
        po.a.init <- if (length( .ipobs.a )) rep( .ipobs.a , n) else
                     rep(sum(css) / n, n)  # MLE for pobs.a (unwted)
      etastart <-
        cbind(theta2eta(munb.init, .lmunb.p , earg = .emunb.p ),
              theta2eta(size.init, .lsize.p , earg = .esize.p ),
              if (lalter <= 0) NULL else
              theta2eta(po.a.init, .lpobs.a , earg = .epobs.a ),
              if (lalter <= 1) NULL else
              theta2eta(munb.init, .lmunb.a , earg = .emunb.a ),
              if (lalter <= 1) NULL else
              theta2eta(size.init, .lsize.a , earg = .esize.a ))
    }
  }), list( .lmunb.p = lmunb.p, .lsize.p = lsize.p, .lpobs.a = lpobs.a,
            .lmunb.a = lmunb.a, .lsize.a = lsize.a,
            .emunb.p = emunb.p, .esize.p = esize.p, .epobs.a = epobs.a,
            .emunb.a = emunb.a, .esize.a = esize.a,
            .imunb.p = imunb.p, .isize.p = isize.p, .ipobs.a = ipobs.a,
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
                     c("mean", "pobs.a", "Pobs.a",
                       "prob.a", "prob.i", "prob.t", "lhs.prob"))[1]
    truncate <- as.vector( .truncate )
    ltrunc <- length(truncate)
    max.support <- as.vector( .max.support )
    alter <- as.vector( .alter )
    lalter <- length(alter)
    M1 <- switch(as.character(lalter), "0" = 2, "1" = 3, 5)
    NOS <- ncol(eta) / M1
    if ((NOS <- ncol(eta) / M1) != 1)
      stop("Currently NOS must be 1")
    munb.p <- eta2theta(eta[, 1], .lmunb.p , earg = .emunb.p )
    size.p <- eta2theta(eta[, 2], .lsize.p , earg = .esize.p )
    pobs.a <- if (lalter <= 0) rep(0, nrow(eta)) else
              eta2theta(eta[, 3], .lpobs.a , earg = .epobs.a )
    munb.a <- if (lalter <= 1) munb.p else
              eta2theta(eta[, 4], .lmunb.a , earg = .emunb.a )
    size.a <- if (lalter <= 1) size.p else
              eta2theta(eta[, 5], .lsize.a , earg = .esize.a )

    Bits <- moments.nbin.gait(size.p = size.p, munb.p = munb.p,
                              alter = alter,
                              truncate = truncate, mlm = FALSE,
                              max.support = max.support,
                              pobs.a = pobs.a, munb.a = munb.a,
                              size.a = size.a)

    if (type.fitted == "Pobs.a") {
      proportion.mat <-
        dnbinom(matrix(alter, NROW(eta), lalter, byrow = TRUE),
                size = matrix(size.a, NROW(eta), lalter),
                mu   = matrix(munb.a, NROW(eta), lalter)) / (
        c(Bits[["suma.a"]]))
    }
  
    ans <- switch(type.fitted,
       "mean"     = Bits[["mean"]],
       "pobs.a"   = pobs.a,
       "Pobs.a"   = pobs.a * proportion.mat,  # matrix
       "prob.a"   = Bits[["suma.p"]],
       "prob.i"   = Bits[["sumi.p"]],
       "prob.t"   = Bits[["sumt.p"]],
      "lhs.prob"  = Bits[["lhs.prob"]])
    ynames.Pobs.a <- as.character(alter)  # Works for NULLs
    label.cols.y(ans,
        colnames.y = switch(type.fitted,
                            "pobs.a" =,
                            "Pobs.a" = ynames.Pobs.a,
                            extra$colnames.y),
        NOS = NOS)
  } , list( .lmunb.p = lmunb.p, .lsize.p = lsize.p, .lpobs.a = lpobs.a,
            .lmunb.a = lmunb.a, .lsize.a = lsize.a,
            .emunb.p = emunb.p, .esize.p = esize.p, .epobs.a = epobs.a,
            .emunb.a = emunb.a, .esize.a = esize.a,
            .alter = alter,
            .truncate = truncate, .max.support = max.support ))),

  last = eval(substitute(expression({
    misc$link <- c(munb.p = .lmunb.p , size.p = .lsize.p )
    if (lalter > 0)
      misc$link <- c(misc$link, pobs.a   = .lpobs.a )
    if (lalter > 1)
      misc$link <- c(misc$link, munb.a = .lmunb.a , size.a = .lsize.a )

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- .emunb.p  #
    misc$earg[[2]] <- .esize.p  #
    if (lalter > 0)
      misc$earg[[3]] <- .epobs.a  #
    if (lalter > 1) {
      misc$earg[[4]] <- .emunb.a  #
      misc$earg[[5]] <- .esize.a  #
    }
  }), list( .lmunb.p = lmunb.p, .lsize.p = lsize.p, .lpobs.a = lpobs.a,
            .lmunb.a = lmunb.a, .lsize.a = lsize.a,
            .emunb.p = emunb.p, .esize.p = esize.p, .epobs.a = epobs.a,
            .emunb.a = emunb.a, .esize.a = esize.a,
            .alter = alter,
            .truncate = truncate, .max.support = max.support ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    alter <- as.vector( .alter )
    lalter <- length(alter)
    munb.p <- eta2theta(eta[, 1], .lmunb.p , earg = .emunb.p )
    size.p <- eta2theta(eta[, 2], .lsize.p , earg = .esize.p )
    pobs.a <- if (lalter <= 0) rep(0, nrow(eta)) else
              eta2theta(eta[, 3], .lpobs.a , earg = .epobs.a )
    munb.a <- if (lalter <= 1) munb.p else
              eta2theta(eta[, 4], .lmunb.a , earg = .emunb.a )
    size.a <- if (lalter <= 1) size.p else
              eta2theta(eta[, 5], .lsize.a , earg = .esize.a )

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
    alter <- as.vector( .alter )
    lalter <- length(alter)
    munb.p <- eta2theta(eta[, 1], .lmunb.p , earg = .emunb.p )
    size.p <- eta2theta(eta[, 2], .lsize.p , earg = .esize.p )
    pobs.a <- if (lalter <= 0) rep(0.5, nrow(eta)) else
              eta2theta(eta[, 3], .lpobs.a , earg = .epobs.a )
    munb.a <- if (lalter <= 1) munb.p else
              eta2theta(eta[, 4], .lmunb.a , earg = .emunb.a )
    size.a <- if (lalter <= 1) size.p else
              eta2theta(eta[, 5], .lsize.a , earg = .esize.a )
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
    munb.p <- eta2theta(eta[, 1], .lmunb.p , earg = .emunb.p )
    size.p <- eta2theta(eta[, 2], .lsize.p , earg = .esize.p )
    pobs.a <- if (lalter <= 0) rep(0, nrow(eta)) else
              eta2theta(eta[, 3], .lpobs.a , earg = .epobs.a )
    munb.a <- if (lalter <= 1) munb.p else
              eta2theta(eta[, 4], .lmunb.a , earg = .emunb.a )
    size.a <- if (lalter <= 1) size.p else
              eta2theta(eta[, 5], .lsize.a , earg = .esize.a )
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
    M1 <- switch(as.character(lalter), "0" = 2, "1" = 3, 5)
    NOS <- NCOL(eta) / M1  # extra$NOS
    if (NOS != 1) stop("Multiple responses not handled")
    max.support <- as.vector( .max.support )
    is.altered <- if (lalter)
      rowSums(extra$skip.these) > 0 else rep(FALSE, n)
    
    munb.p <- eta2theta(eta[, 1], .lmunb.p , earg = .emunb.p )
    size.p <- eta2theta(eta[, 2], .lsize.p , earg = .esize.p )
    pobs.a <- if (lalter <= 0) rep(0, nrow(eta)) else
              eta2theta(eta[, 3], .lpobs.a , earg = .epobs.a )
    munb.a <- if (lalter <= 1) munb.p else
              eta2theta(eta[, 4], .lmunb.a , earg = .emunb.a )
    size.a <- if (lalter <= 1) size.p else
              eta2theta(eta[, 5], .lsize.a , earg = .esize.a )

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
    sumderiv1a.p <- sumderiv2a.p <- matrix(0, nrow(eta), NOS)
    deriv0matrix <-
    deriv1matrix <- deriv2matrix <- matrix(0, nrow(eta), lalter)
    SumDeriv1a.a <- SumDeriv2a.a <-
    SumDeriv1a.p <- SumDeriv2a.p <- matrix(0, nrow(eta), NOS)
    DerivxMatrix <-  # Deriv0Matrix <-
    Deriv1Matrix <- Deriv2Matrix <- matrix(0, nrow(eta), lalter)
    if (lalter) {
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

    }  # jay
    deriv0matrix <-  deriv0matrix / rowSums(deriv0matrix)  # Normalized
    }  # lalter

    sumderivxt.p <- sumderivxt.a <-
    sumderiv1t.a <- sumderiv2t.a <-
    sumderiv1t.p <- sumderiv2t.p <- matrix(0, nrow(eta), NOS)
    SumDeriv1t.a <- SumDeriv2t.a <-
    SumDeriv1t.p <- SumDeriv2t.p <- matrix(0, nrow(eta), NOS)
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
    if (lalter > 1)
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
    ans <- cbind(dl.dmunb.p * dmunb.p.deta,
                 dl.dsize.p * dsize.p.deta,
                 if (lalter <= 0) NULL else dl.dpobs.a * dpobs.a.deta,
                 if (lalter <= 1) NULL else dl.dmunb.a * dmunb.a.deta,
                 if (lalter <= 1) NULL else dl.dsize.a * dsize.a.deta)
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


    if (lalter > 1)
    ned2l.dmunb.a2 <- pobs.a *
      rowSums(deriv0matrix * (
             -deriv2matrix +
              deriv1matrix^2 +
             (c(sumderiv2a.a) / c(bits[["suma.a"]])) -
             (c(sumderiv1a.a) / c(bits[["suma.a"]]))^2))


    if (lalter > 1)
    ned2l.dsize.a2 <- pobs.a *
      rowSums(deriv0matrix * (
             -Deriv2Matrix +
              Deriv1Matrix^2 +
             (c(SumDeriv2a.a) / c(bits[["suma.a"]])) -
             (c(SumDeriv1a.a) / c(bits[["suma.a"]]))^2))


    if (lalter > 1)
    ned2l.dmunb.a.size.a2 <- pobs.a *
      rowSums(deriv0matrix * (
             -DerivxMatrix +
              deriv1matrix * Deriv1Matrix +
            c(sumderivxa.a) / c(bits[["suma.a"]]) -
            c(sumderiv1a.a) * c(SumDeriv1a.a) / c(bits[["suma.a"]])^2))



    if (lalter > 0)
    wz11 <- if ( .lpobs.a == "logitlink") {
      pobs.a * (1 - pobs.a)
    } else {
      dpobs.a.deta <- dtheta.deta(pobs.a, .lpobs.a , .epobs.a )
      ned2l.dpobs.a2 <- 1 / (pobs.a * (1 - pobs.a))
      ned2l.dpobs.a2 * dpobs.a.deta^2
    }






    min.support <- 0  # Usual case; same as lowsup
    min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else
    min.support
    max.chunk.MB <- ( .max.chunk.MB )

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
      Q.maxs <-  qgaitnbinom.mix(p           = eff.p[2] ,
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
          eim.gatnb.p <- EIM.GATNB.speciald(
            munb        = munb.p[sind2],
            size        = size.p[sind2],
            munb.a = munb.a[sind2], size.a = size.a[sind2],
            alter = alter, truncate = truncate,
            max.support = max.support,
            pobs.a      = 0,  # Importantly not pobs.a[sind2],
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





    ned2l.dsize.p2 <- ned2l.dsize.p2 -
      (SumDeriv2a.p + SumDeriv2t.p) / Denom -
     ((SumDeriv1a.p + SumDeriv1t.p) / Denom)^2


    ned2l.dsize.p2 <- onempobs.a * ned2l.dsize.p2



    wz <- matrix(0, n, 2*M1-1)
    wz[, iam(1, 1, M1)] <- ned2l.dmunb.p2 * dmunb.p.deta^2
    wz[, iam(2, 2, M1)] <- ned2l.dsize.p2 * dsize.p.deta^2
    wz[, iam(1, 2, M1)] <- ned2l.dmunb.p.size.p2 * dmunb.p.deta *
                                                   dsize.p.deta
    if (lalter > 0)
      wz[, iam(3, 3, M1)] <- wz11
    if (lalter > 1) {
      wz[, iam(4, 4, M1)] <- ned2l.dmunb.a2 * dmunb.a.deta^2
      wz[, iam(5, 5, M1)] <- ned2l.dsize.a2 * dsize.a.deta^2
      wz[, iam(4, 5, M1)] <- ned2l.dmunb.a.size.a2 * dmunb.a.deta *
                                                     dsize.a.deta
    }

    wz
  }), list( .alter = alter,
            .truncate = truncate, .max.support = max.support,
            .cutoff.prob = cutoff.prob, .eps.trig = eps.trig,
            .max.chunk.MB = max.chunk.MB, .lpobs.a = lpobs.a
           ))))
}  # gatnbinomial.mix









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









 semigait.errorcheck <-
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
}  # semigait.errorcheck








 dgaitnbinom.mlm <-
  function(x, size, prob = NULL, munb = NULL,
           alter = NULL,
           inflate = NULL,
           truncate = NULL, max.support = Inf,
           pobs.a = 0, pstr.i = 0, byrow.arg = FALSE,
           log.arg = FALSE) {
  semigait.errorcheck(alter, inflate, truncate, max.support)
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
  semigait.errorcheck(alter, inflate, truncate, max.support)

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
  semigait.errorcheck(alter, inflate, truncate, max.support)
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












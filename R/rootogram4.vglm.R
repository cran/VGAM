# These functions are
# Copyright (C) 1998-2020 T.W. Yee, University of Auckland.
# All rights reserved.








 rootogram4vglm <-
  function(object, newdata = NULL, breaks = NULL, max = NULL, 
           xlab = NULL, main = NULL, width = NULL,
           ...) {


  vfamily <- object@family@vfamily[1]






  mf <- if (is.null(newdata)) {
    model.frame(object)
  } else {
   mt <- terms(object)
   model.frame(mt, newdata, na.action = na.omit)
  }
  if (is.null(newdata))
    mt <- attr(mf, "terms")
  y <- Y <- model.response(mf)
  if (!is.factor(Y))
    y <- Y <- as.matrix(Y)
  n.lm <- NROW(Y)
  w <- Wts <- model.weights(mf)
  if (length(Wts) == 0L)
    w <- Wts <- rep(1, n.lm)  # Safest (uses recycling and is a vector)



  Q1.infos <- M.infos <- M1.infos <- 1  # Default really
  infos.fun <- object@family@infos
  infos.list <- infos.fun()
  if (is.list(infos.list) && length(infos.list$M))
    M.infos <- infos.list$M
  if (is.list(infos.list) && length(infos.list$M1))
    M1.infos <- infos.list$M1
  if (is.list(infos.list) && length(infos.list$Q1))
    Q1.infos <- infos.list$Q1
  if ((NOS <- M.infos / M1.infos) != 1)
    stop("can only handle one response")
  if (Q1.infos != 1)
    stop("Q1 must be unity")
  M <- M.infos



  link1parameter <- infos.list$link1parameter
  if (is.null(link1parameter))
    link1parameter <- TRUE  # The default, for ordinary 1-par links
  eta.mat <- predict(object)
  n.LM <- NROW(eta.mat)
  Param.mat <- matrix(NA_real_, n.LM, M)
  mylinks <- linkfun(object)  # Of length 1 for GLMs, char only



  multipleResponses <-
    if (is.logical((tmp3 <- infos.list$multipleResponses))) tmp3 else
      FALSE
  mixture.links <-
    if (is.logical((tmp3 <- infos.list$mixture.links))) tmp3 else
      FALSE


  lowsup <- 0L  # Default (at least for count distributions)

      
  GAITffs <-
      c("gatpoisson.mix", "gatpoisson.mlm",
        "gitpoisson.mix", "gitpoisson.mlm",
        "gaitlog.mix",
        "gaitpoisson.mix", "gaitpoisson.mlm",
        "gatnbinomial.mix", "gatnbinomial.mlm",
        "gitnbinomial.mix", "gitnbinomial.mlm")
  if (vfamily %in% GAITffs) {
    alter <- infos.list$alter
    inflate <- infos.list$inflate
    truncate <- infos.list$truncate
    if (is.numeric(tmp9 <- infos.list$lowsup))
      lowsup <- tmp9
    max.support <- infos.list$max.support  # Often Inf but never NULL
    if ((M1 <- ncol(eta.mat)) != M1.infos)
      stop("confused about the variable 'M1'")
    pobs.a <- pstr.i <- 0  # Initialize; maybe NULL would be safer
    if (length(alter))
      pobs.a <- fitted(object, type.fitted = "pobs.a")
    if (length(inflate))
      pstr.i <- fitted(object, type.fitted = "pstr.i")
  }  # GAITffs



      if (!mixture.links &&  # !multipleResponses && 
          link1parameter) {
    for (jay in 1:M) {      # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      Param.mat[, jay] <-
      Param.vec <- eta2theta(eta.mat[, jay], mylinks[jay],
                             earg = object@misc$earg[[jay]])
    }  # for (jay)
  } else {
  }


  



  if (!(vfamily %in% c("uninormal", "binomialff",
                       "posbinomial",
                       "zabinomial", "zabinomialff",
                       "zibinomial", "zibinomialff"))) {
    max0 <- if (is.null(max)) max(1.5 * max(y[w > 0]), 20L) else max
    obsrvd <- as.vector(xtabs(w ~ factor(y, levels = lowsup:max0)))
   At <- lowsup:max0
  }



  pmat <- matrix(NA_real_, n.lm, length(At))  # Important initialization

  if (vfamily == "bellff") {
    for (i in At)
      pmat[, i + 1L] <- dbell(i, shape = Param.mat[, 1])
  }
  if (vfamily == "borel.tanner") {
    for (i in At)
      pmat[, i + 1L] <- dbort(i, a = Param.mat[, 1],
                              Qsize = object@misc$Qsize)
  }
  if (vfamily == "diffzeta") {
    for (i in At)
      pmat[, i + 1L] <- ddiffzeta(i, shape = Param.mat[, 1],
                                  start = object@misc$start)
  }
  if (vfamily == "genpoisson") {
    for (i in At)
      pmat[, i + 1L] <- dgenpois(i, theta = Param.mat[, 2],
                                 lambda = Param.mat[, 1])
  }
  if (vfamily == "gatnbinomial.mix") {
    for (i in At)
      pmat[, i + 1L] <- dgaitnbinom.mix(i, alter = alter,
                                        truncate = truncate,
                                        max.support = max.support,
                                        size.p = Param.mat[, 2],
                                        munb.p = Param.mat[, 1],
                                        size.a = Param.mat[, 5],
                                        munb.a = Param.mat[, 4],
                                        pobs.a = Param.mat[, 3])
  }
  if (vfamily == "gatnbinomial.mlm") {
    munb <- eta2theta(eta.mat[, 1], mylinks[1],
                      earg = object@misc$earg[[1]])
    size <- eta2theta(eta.mat[, 2], mylinks[2],
                      earg = object@misc$earg[[2]])
    for (i in At)
      pmat[, i + 1L] <- dgaitnbinom.mlm(i, alter = alter,
                                        truncate = truncate,
                                        max.support = max.support,
                                        munb = munb, size = size,
                                        pobs.a = pobs.a[, -ncol(pobs.a)],
                                        byrow.arg = FALSE)
  }
  if (vfamily == "gaitlog.mix") {
    is.altered <- as.logical(length(alter))
    is.inflated <- as.logical(length(inflate))
    index.a <- which(infos.list$parameters.names == "shape.a")
    index.i <- which(infos.list$parameters.names == "shape.i")
    shape.a <- Param.mat[, if (is.altered)  index.a else 1]
    shape.i <- Param.mat[, if (is.inflated) index.i else 1]
    for (i in At)
      pmat[, i + 0L] <-
        dgaitlog(i, alter.mix = alter, inflate.mix = inflate,
                 truncate = truncate,
                 shape.p = Param.mat[, 1],
                 pobs.mix.a = pobs.a,
                 pstr.mix.i = pstr.i,
                 shape.a = shape.a,
                 shape.i = shape.i)
  }
  if (vfamily == "gaitpoisson.mix") {
    is.altered <- as.logical(length(alter))
    is.inflated <- as.logical(length(inflate))
    index.a <- which(infos.list$parameters.names == "lambda.a")
    index.i <- which(infos.list$parameters.names == "lambda.i")
    lambda.a <- Param.mat[, if (is.altered)  index.a else 1]
    lambda.i <- Param.mat[, if (is.inflated) index.i else 1]
    for (i in At)
      pmat[, i + 1L] <-
        dgaitpois(i, alter.mix = alter, inflate.mix = inflate,
                  truncate = truncate,
                  max.support = max.support,
                  lambda.p = Param.mat[, 1],
                  pobs.mix.a = pobs.a,
                  pstr.mix.i = pstr.i,
                  lambda.a = lambda.a,
                  lambda.i = lambda.i)
  }
  if (vfamily == "gaitpoisson.mlm") {
    is.altered <- as.logical(length(alter))
    is.inflated <- as.logical(length(inflate))
    Pobs.a <- if (is.altered)
      fitted(object, type.fitted = "Pobs.a") else 0
    Pstr.i <- if (is.inflated)
      fitted(object, type.fitted = "Pstr.i") else 0
 # Param.mat is almost all NAs, so dont use it at all.
    lambda.p <- eta2theta(eta.mat[, 1], mylinks[1],
                          earg = object@misc$earg[[1]])
    for (i in At) aaa <-
      pmat[, i + 1L] <-
        dgaitpois(i, alter.mlm = alter, inflate.mlm = inflate,
                  truncate = truncate,  # byrow.arg = F,
                  max.support = max.support,
                  lambda.p = lambda.p,  # There are 2 twists:
                  pobs.mlm.a = if (is.altered)
                    Pobs.a[, -NCOL(Pobs.a)] else 0,
                  pstr.mlm.i = if (is.inflated)
                    Pstr.i[, -NCOL(Pstr.i)] else 0)
  }
  if (vfamily == "gatpoisson.mix") {
    for (i in At)
      pmat[, i + 1L] <- dgaitpois(i, alter.mix = alter,
                                  truncate = truncate,
                                  max.support = max.support,
                                  lambda.p   = Param.mat[, 1],
                                  pobs.mix.a = Param.mat[, 2],
                                  lambda.a   = Param.mat[, 3])
  }
  if (vfamily == "gatpoisson.mlm") {
    lambda <- eta2theta(eta.mat[, 1], mylinks[1],
                        earg = object@misc$earg[[1]])
    for (i in At)
      pmat[, i + 1L] <- dgaitpois(i, alter.mlm = alter,
                                  truncate = truncate,
                                  max.support = max.support,
                                  lambda.p = lambda,
                                  pobs.mlm.a = pobs.a[, -ncol(pobs.a)],
                                  byrow.arg = FALSE)
  }
  if (vfamily == "gitpoisson.mix") {
    for (i in At)
      pmat[, i + 1L] <- dgaitpois(i, inflate.mix = inflate,
                                  truncate = truncate,
                                  max.support = max.support,
                                  lambda.p   = Param.mat[, 1],
                                  pstr.mix.i = Param.mat[, 2],
                                  lambda.i   = Param.mat[, 3])
  }
  if (vfamily == "gitpoisson.mlm") {
    lambda <- eta2theta(eta.mat[, 1], mylinks[1],
                        earg = object@misc$earg[[1]])
    for (i in At)
      pmat[, i + 1L] <- dgaitpois(i, inflate.mlm = inflate,
                                  truncate = truncate,
                                  max.support = max.support,
                                  lambda.p = lambda,
                                  pstr.mlm.i = pstr.i[, -ncol(pstr.i)],
                                  byrow.arg = FALSE)
  }
  if (vfamily == "geometric") {
    for (i in At)
      pmat[, i + 1L] <- dgeom(i, prob = Param.mat[, 1])
  }
  if (vfamily == "hzeta") {
    for (i in At)
      pmat[, i + 1L] <- dhzeta(i, shape = Param.mat[, 1])
  }
  if (vfamily == "negbinomial") {
    for (i in At)
      pmat[, i + 1L] <- dnbinom(i, size = Param.mat[, 2],
                                mu = Param.mat[, 1])
  }
  if (vfamily == "negbinomial,size") {
    sizevec <- object@misc$size
    if (length(unique(sizevec)) > 1)
      stop("the size values must be all the same")
    for (i in At)
      pmat[, i + 1L] <- dnbinom(i, size = sizevec[1],
                                mu = Param.mat[, 1])
  }
  if (vfamily == "poissonff") {
    for (i in At)
      pmat[, i + 1L] <- dpois(i, lambda = Param.mat[, 1])
  }
  if (vfamily == "polya") {
    for (i in At)
      pmat[, i + 1L] <- dnbinom(i, size = Param.mat[, 2],
                                prob = Param.mat[, 1])
  }
  if (vfamily == "polyaR") {
    for (i in At)
      pmat[, i + 1L] <- dnbinom(i, size = Param.mat[, 1],
                                prob = Param.mat[, 2])
  }
  if (vfamily == "posnegbinomial") {
    for (i in At)
      pmat[, i + 1L] <- dposnegbin(i, size = Param.mat[, 2],
                                   munb = Param.mat[, 1])
  }
  if (vfamily == "truncgeometric") {
    upper.limit <- object@extra$upper.limit
    if (length(unique(upper.limit)) > 1)
      stop("the upper.limit values must be all the same")
    prob <- Param.mat[, 1]
    for (i in At)
      pmat[, i + 1L] <- (dgeom(i, prob = prob) /
                        (1 - (1.0 - prob)^(1 + upper.limit)))
  }
  if (vfamily == "yulesimon") {
    for (i in At)
      pmat[, i + 1L] <- dyules(i, shape = Param.mat[, 1])
  }     
  if (vfamily == "zanegbinomial") {
    for (i in At)
      pmat[, i + 1L] <- dzanegbin(i, pobs0 = Param.mat[, 1],
                                  munb = Param.mat[, 2],
                                  size = Param.mat[, 3])
  }
  if (vfamily == "zanegbinomialff") {
    for (i in At)
      pmat[, i + 1L] <- dzanegbin(i, pobs0 = 1 - Param.mat[, 3],
                                  munb = Param.mat[, 1],
                                  size = Param.mat[, 2])
  }
  if (vfamily == "zapoisson") {
    for (i in At)
      pmat[, i + 1L] <- dzapois(i, pobs0 = Param.mat[, 1],
                                lambda = Param.mat[, 2])
  }
  if (vfamily == "zapoissonff") {
    for (i in At)
      pmat[, i + 1L] <- dzapois(i, pobs0 = 1 - Param.mat[, 2],
                                lambda = Param.mat[, 1])
  }
  if (vfamily == "zigeometric") {
    for (i in At)
      pmat[, i + 1L] <- dzigeom(i, pstr0 = Param.mat[, 1],
                                prob = Param.mat[, 2])
  }
  if (vfamily == "zigeometricff") {
    for (i in At)
      pmat[, i + 1L] <- dzigeom(i, pstr0 = 1 - Param.mat[, 2],
                                prob = Param.mat[, 1])
  }
  if (vfamily == "zinegbinomial") {
    for (i in At)
      pmat[, i + 1L] <- dzinegbin(i, pstr0 = Param.mat[, 1],
                                  munb = Param.mat[, 2],
                                  size = Param.mat[, 3])
  }
  if (vfamily == "zinegbinomialff") {
    for (i in At)
      pmat[, i + 1L] <- dzinegbin(i, pstr0 = 1 - Param.mat[, 3],
                                  munb = Param.mat[, 1],
                                  size = Param.mat[, 2])
  }
  if (vfamily == "zipf") {
    for (i in At)
      pmat[, i + 1L] <- dzipf(i, N = object@misc$N,  # zz,
                              shape = Param.mat[, 1])
  }
  if (vfamily == "zipoisson") {
    for (i in At)
      pmat[, i + 1L] <- dzipois(i, pstr0 = Param.mat[, 1],
                                lambda = Param.mat[, 2])
  }
  if (vfamily == "zipoissonff") {
    for (i in At)
      pmat[, i + 1L] <- dzipois(i, pstr0 = 1 - Param.mat[, 2],
                                lambda = Param.mat[, 1])
  }
  if (vfamily == "zz") {
    for (i in At)
      pmat[, i + 1L] <- dpois(i, lambda = Param.mat[, 1])
  }




  if (!all(is.na(pmat))) {
    expctd <- colSums(pmat * w)
 ## try to guess a good maximum
    if (is.null(max)) {
      max <- if (all(expctd >= 1L)) max0 else
             max(ceiling(mean(y)), min(which(expctd < 1L)) - 1L)
      max <- min(max, length(expctd) - 1L)
    }  # is.null(max)
    breaks <- ((lowsup - 1L):max) + 0.5
    obsrvd <- obsrvd[1L:(length(breaks) - 1L)]
    expctd <- expctd[1L:(length(breaks) - 1L)]
  }  # Not uninormal or binomialff





  
  if (vfamily == "uninormal") {  # Continuous distn
    mu <- Param.mat[, 1]
    s <- if (infos.list$var.arg)
         sqrt(Param.mat[, 2]) else Param.mat[, 2]
    if (is.null(breaks)) 
      breaks <- "Sturges"
    breaks <- hist(y[w > 0], plot = FALSE, breaks = breaks)$breaks
    obsrvd <- as.vector(xtabs(w ~ cut(y, breaks, include.lowest = TRUE)))
    pmat <- matrix(NA, nrow = length(y), ncol = length(breaks) - 1L)
    for (i in 1L:ncol(pmat))
      pmat[, i] <- pnorm(breaks[i + 1L], mean = mu, sd = s) -
                   pnorm(breaks[i],      mean = mu, sd = s)
    expctd <- colSums(pmat * w)
  }

 if (vfamily == "binomialff") {
    if (NCOL(y) < 2L) 
      y <- cbind(y, 1L - y)
    size <- unique(rowSums(y))
    if (length(size) > 1L) 
      stop("rootogram4 only applicable to binomial ",
           "distributions with same size")
    At <- 0L:size
    breaks <- -1L:size + 0.5
    obsrvd <- as.vector(xtabs(w ~ factor(y[, 1L], levels = At)))
    pmat <- matrix(NA, length(mu), length(At))
    for (i in At)
      pmat[, i + 1L] <- dbinom(i, prob = mu, size = size)
    expctd <- colSums(pmat * w)
  }

 if (vfamily %in% c("posbinomial")) {
    if (NCOL(y) < 2L) 
      y <- cbind(y, 1L - y)
    size <- unique(rowSums(y))
    if (length(size) > 1L) 
        stop("rootogram4 only applicable to posbinomial",
             " distributions with same size")
    At <- 0L:size
    breaks <- -1L:size + 0.5
    obsrvd <- as.vector(xtabs(w ~ factor(y[, 1L], levels = At)))
    pmat <- matrix(NA, length(mu), length(At))
    for (i in At)
      pmat[, i + 1L] <- dposbinom(i, prob = Param.mat[, 1],
                                  size = size)
    expctd <- colSums(pmat * w)
  }

 if (vfamily %in% c("zabinomial", "zabinomialff")) {
    if (NCOL(y) < 2L) 
      y <- cbind(y, 1L - y)
    size <- unique(rowSums(y))
    if (length(size) > 1L) 
        stop("rootogram4 only applicable to zabinomial",
             ifelse(vfamily == "zabinomial", " ", "ff "),
             "distributions with same size")
    At <- 0L:size
    breaks <- -1L:size + 0.5
    obsrvd <- as.vector(xtabs(w ~ factor(y[, 1L], levels = At)))
    pmat <- matrix(NA, length(mu), length(At))
    for (i in At)
      pmat[, i + 1L] <-
        dzabinom(i,
             prob = Param.mat[, ifelse(vfamily == "zabinomial", 2, 1)],
             size = size,
             pobs0 = if (vfamily == "zabinomial") Param.mat[, 1] else
                     1 - Param.mat[, 2])
    expctd <- colSums(pmat * w)
  }

 if (vfamily %in% c("zibinomial", "zibinomialff")) {
    if (NCOL(y) < 2L) 
      y <- cbind(y, 1L - y)
    size <- unique(rowSums(y))
    if (length(size) > 1L) 
        stop("rootogram4 only applicable to zibinomial",
             ifelse(vfamily == "zibinomial", " ", "ff "),
             "distributions with same size")
    At <- (0L:size) + lowsup
    breaks <- -1L:size + 0.5
    obsrvd <- as.vector(xtabs(w ~ factor(y[, 1L], levels = At)))
    pmat <- matrix(NA, length(mu), length(At))
    for (i in At)
      pmat[, i + 1L] <-
        dzibinom(i,
             prob = Param.mat[, ifelse(vfamily == "zibinomial", 2, 1)],
             size = size,
             pstr0 = if (vfamily == "zibinomial") Param.mat[, 1] else
                     1 - Param.mat[, 2])
    expctd <- colSums(pmat * w)
  }


  if (all(is.na(pmat)))
    stop("family '", vfamily, "' currently not supported")

  
  if (is.null(xlab))
    xlab <- as.character(attr(mt, "variables"))[2L]
  if (is.null(main))
    main <- deparse(substitute(object))






  rootogram0.default(obsrvd, expctd, breaks = breaks, xlab = xlab, 
                     main = main, width = if (vfamily == "uninormal") 
                     1 else 0.9,
                     lowsup = lowsup, ...)
}  # rootogram4vglm












if (!isGeneric("rootogram4"))
  setGeneric("rootogram4", function(object, ...)
             standardGeneric("rootogram4"),
             package = "VGAM")


setMethod("rootogram4", "vglm",
          function(object, ...)
          rootogram4vglm(object, ...))



















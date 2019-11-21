# These functions are
# Copyright (C) 1998-2019 T.W. Yee, University of Auckland.
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
  Ans.infos <- infos.fun()
  if (is.list(Ans.infos) && length(Ans.infos$M))
    M.infos <- Ans.infos$M
  if (is.list(Ans.infos) && length(Ans.infos$M1))
    M1.infos <- Ans.infos$M1
  if (is.list(Ans.infos) && length(Ans.infos$Q1))
    Q1.infos <- Ans.infos$Q1
  if ((NOS <- M.infos / M1.infos) != 1)
    stop("can only handle one response")
  if (Q1.infos != 1)
    stop("Q1 must be unity")
  M <- M.infos



  link1parameter <- Ans.infos$link1parameter
  if (is.null(link1parameter))
    link1parameter <- TRUE  # The default, for ordinary 1-par links
  if (!link1parameter)
    stop("can only handle 1-parameter links")
  eta.mat <- predict(object)
  n.LM <- NROW(eta.mat)
  Param.mat <- matrix(NA_real_, n.LM, M)
  mylinks <- linkfun(object)  # Of length 1 for GLMs, char only

  for (jay in 1:M) {      # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    Param.mat[, jay] <-
    Param.vec <- eta2theta(eta.mat[, jay], mylinks[jay],
                           earg = object@misc$earg[[jay]])
  }  # for (jay)
      

  if (!(vfamily %in% c("uninormal", "binomialff",
                       "posbinomial",
                       "zabinomial", "zabinomialff",
                       "zibinomial", "zibinomialff"))) {
    max0 <- if (is.null(max)) max(1.5 * max(y[w > 0]), 20L) else max
    max0 <- if (is.null(max)) max(1.5 * max(y[w > 0]), 20L) else max
    obsrvd <- as.vector(xtabs(w ~ factor(y, levels = 0L:max0)))
    at <- 0L:max0
  }

  pmat <- matrix(NA_real_, n.lm, length(at))

  if (vfamily == "bellff") {
    for (i in at)
      pmat[, i + 1L] <- dbell(i, shape = Param.mat[, 1])
  }
  if (vfamily == "borel.tanner") {
    for (i in at)
      pmat[, i + 1L] <- dbort(i, a = Param.mat[, 1],
                              Qsize = object@misc$Qsize)
  }
  if (vfamily == "diffzeta") {
    for (i in at)
      pmat[, i + 1L] <- ddiffzeta(i, shape = Param.mat[, 1],
                                  start = object@misc$start)
  }
  if (vfamily == "genpoisson") {
    for (i in at)
      pmat[, i + 1L] <- dgenpois(i, theta = Param.mat[, 2],
                                 lambda = Param.mat[, 1])
  }
  if (vfamily == "geometric") {
    for (i in at)
      pmat[, i + 1L] <- dgeom(i, prob = Param.mat[, 1])
  }
  if (vfamily == "hzeta") {
    for (i in at)
      pmat[, i + 1L] <- dhzeta(i, shape = Param.mat[, 1])
  }
  if (vfamily == "negbinomial") {
    for (i in at)
      pmat[, i + 1L] <- dnbinom(i, size = Param.mat[, 2],
                                mu = Param.mat[, 1])
  }
  if (vfamily == "negbinomial,size") {
    sizevec <- object@misc$size
    if (length(unique(sizevec)) > 1)
      stop("the size values must be all the same")
    for (i in at)
      pmat[, i + 1L] <- dnbinom(i, size = sizevec[1],
                                mu = Param.mat[, 1])
  }
  if (vfamily == "poissonff") {
    for (i in at)
      pmat[, i + 1L] <- dpois(i, lambda = Param.mat[, 1])
  }
  if (vfamily == "polya") {
    for (i in at)
      pmat[, i + 1L] <- dnbinom(i, size = Param.mat[, 2],
                                prob = Param.mat[, 1])
  }
  if (vfamily == "polyaR") {
    for (i in at)
      pmat[, i + 1L] <- dnbinom(i, size = Param.mat[, 1],
                                prob = Param.mat[, 2])
  }
  if (vfamily == "posnegbinomial") {
    for (i in at)
      pmat[, i + 1L] <- dposnegbin(i, size = Param.mat[, 2],
                                   munb = Param.mat[, 1])
  }
  if (vfamily == "truncgeometric") {
    upper.limit <- object@extra$upper.limit
    if (length(unique(upper.limit)) > 1)
      stop("the upper.limit values must be all the same")
    prob <- Param.mat[, 1]
    for (i in at)
      pmat[, i + 1L] <- (dgeom(i, prob = prob) /
                        (1 - (1.0 - prob)^(1 + upper.limit)))
  }
  if (vfamily == "yulesimon") {
    for (i in at)
      pmat[, i + 1L] <- dyules(i, shape = Param.mat[, 1])
  }     
  if (vfamily == "zanegbinomial") {
    for (i in at)
      pmat[, i + 1L] <- dzanegbin(i, pobs0 = Param.mat[, 1],
                                  munb = Param.mat[, 2],
                                  size = Param.mat[, 3])
  }
  if (vfamily == "zanegbinomialff") {
    for (i in at)
      pmat[, i + 1L] <- dzanegbin(i, pobs0 = 1 - Param.mat[, 3],
                                  munb = Param.mat[, 1],
                                  size = Param.mat[, 2])
  }
  if (vfamily == "zapoisson") {
    for (i in at)
      pmat[, i + 1L] <- dzapois(i, pobs0 = Param.mat[, 1],
                                lambda = Param.mat[, 2])
  }
  if (vfamily == "zapoissonff") {
    for (i in at)
      pmat[, i + 1L] <- dzapois(i, pobs0 = 1 - Param.mat[, 2],
                                lambda = Param.mat[, 1])
  }
  if (vfamily == "zigeometric") {
    for (i in at)
      pmat[, i + 1L] <- dzigeom(i, pstr0 = Param.mat[, 1],
                                prob = Param.mat[, 2])
  }
  if (vfamily == "zigeometricff") {
    for (i in at)
      pmat[, i + 1L] <- dzigeom(i, pstr0 = 1 - Param.mat[, 2],
                                prob = Param.mat[, 1])
  }
  if (vfamily == "zinegbinomial") {
    for (i in at)
      pmat[, i + 1L] <- dzinegbin(i, pstr0 = Param.mat[, 1],
                                  munb = Param.mat[, 2],
                                  size = Param.mat[, 3])
  }
  if (vfamily == "zinegbinomialff") {
    for (i in at)
      pmat[, i + 1L] <- dzinegbin(i, pstr0 = 1 - Param.mat[, 3],
                                  munb = Param.mat[, 1],
                                  size = Param.mat[, 2])
  }
  if (vfamily == "zipf") {
    for (i in at)
      pmat[, i + 1L] <- dzipf(i, N = object@misc$N,  # zz,
                              shape = Param.mat[, 1])
  }
  if (vfamily == "zipoisson") {
    for (i in at)
      pmat[, i + 1L] <- dzipois(i, pstr0 = Param.mat[, 1],
                                lambda = Param.mat[, 2])
  }
  if (vfamily == "zipoissonff") {
    for (i in at)
      pmat[, i + 1L] <- dzipois(i, pstr0 = 1 - Param.mat[, 2],
                                lambda = Param.mat[, 1])
  }
  if (vfamily == "zz") {
    for (i in at)
      pmat[, i + 1L] <- dpois(i, lambda = Param.mat[, 1])
  }

  if (!all(is.na(pmat))) {
    expctd <- colSums(pmat * w)
 ## try to guess a good maximum
    if (is.null(max)) {
      max <- if (all(expctd >= 1L)) max0 else
             max(ceiling(mean(y)), min(which(expctd < 1L)) - 1L)
      max <- min(max, length(expctd) - 1L)
    }
    breaks <- -1L:max + 0.5
    obsrvd <- obsrvd[1L:(max + 1L)]
    expctd <- expctd[1L:(max + 1L)]
  }  # Not uninormal or binomialff



  
  if (vfamily == "uninormal") {  # Continuous distn
    mu <- Param.mat[, 1]
    s <- if (Ans.infos$var.arg)
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
    at <- 0L:size
    breaks <- -1L:size + 0.5
    obsrvd <- as.vector(xtabs(w ~ factor(y[, 1L], levels = at)))
    pmat <- matrix(NA, length(mu), length(at))
    for (i in at)
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
    at <- 0L:size
    breaks <- -1L:size + 0.5
    obsrvd <- as.vector(xtabs(w ~ factor(y[, 1L], levels = at)))
    pmat <- matrix(NA, length(mu), length(at))
    for (i in at)
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
    at <- 0L:size
    breaks <- -1L:size + 0.5
    obsrvd <- as.vector(xtabs(w ~ factor(y[, 1L], levels = at)))
    pmat <- matrix(NA, length(mu), length(at))
    for (i in at)
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
    at <- 0L:size
    breaks <- -1L:size + 0.5
    obsrvd <- as.vector(xtabs(w ~ factor(y[, 1L], levels = at)))
    pmat <- matrix(NA, length(mu), length(at))
    for (i in at)
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
                      1 else 0.9, ...)
}  # rootogram4vglm












if (!isGeneric("rootogram4"))
  setGeneric("rootogram4", function(object, ...)
             standardGeneric("rootogram4"),
             package = "VGAM")


setMethod("rootogram4", "vglm",
          function(object, ...)
          rootogram4vglm(object, ...))



















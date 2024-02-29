# These functions are
# Copyright (C) 1998-2024 T.W. Yee, University of Auckland.
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



  Q1.infos <- M.infos <- M1.infos <- npred(object)  # Default really
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

      


  GAITDffs <-
      c("gaitdpoisson", "gaitdlog", "gaitdzeta",
         "gaitdnbinomial")




      
  if (vfamily %in% GAITDffs) {
    spvals <- specials(object)
    a.mix  <- spvals$a.mix  # Might be NULL
    a.mlm  <- spvals$a.mlm  # Might be NULL
    i.mix  <- spvals$i.mix  # Might be NULL
    i.mlm  <- spvals$i.mlm  # Might be NULL
    d.mix  <- spvals$d.mix  # Might be NULL
    d.mlm  <- spvals$d.mlm  # Might be NULL
    truncate <- spvals$truncate  # Might be NULL
    max.support <- spvals$max.support  # Often Inf
    if (length(tmp9 <- infos.list$Support) >= 3)
      lowsup <- tmp9[1]  # Replace lowsup


    if (infos.list$MM1 > 2)
      stop("can only handle 1 & 2-parameter distributions currently")




    if ((M1 <- ncol(eta.mat)) != M1.infos)
      stop("confused about the variable 'M1'")
    pobs.mix <- pstr.mix <- pdip.mix <-
    pobs.mlm <- pstr.mlm <- pdip.mlm <- 0  # Initialize;
    if (la.mix <- length(a.mix))
      pobs.mix <- fitted(object, type.fitted = "pobs.mix")
    if (la.mlm <- length(a.mlm))
      pobs.mlm <- fitted(object, type.fitted = "pobs.mlm")
    if (li.mix <- length(i.mix))
      pstr.mix <- fitted(object, type.fitted = "pstr.mix")
    if (li.mlm <- length(i.mlm))
      pstr.mlm <- fitted(object, type.fitted = "pstr.mlm")
    if (ld.mix <- length(d.mix))
      pdip.mix <- fitted(object, type.fitted = "pdip.mix")
    if (ld.mlm <- length(d.mlm))
      pdip.mlm <- fitted(object, type.fitted = "pdip.mlm")




    thetanames <- infos.list$baseparams.argnames
    Thetas.z <- fitted(object, type.fitted = paste0(thetanames[1], "s"))
    Theta.p1 <- Thetas.z[, 1]  # Always
    Theta.a1 <- Theta.i1 <- Theta.d1 <- Theta.p1  # Needed
    tmp3.TF <- !is.na(rowSums(object@extra$indeta))
    if (infos.list$MM1 == 1) {
      if (tmp3.TF[ 3])
        Theta.a1 <- Thetas.z[, paste0(thetanames[1], ".a")]
      if (tmp3.TF[ 5])
        Theta.i1 <- Thetas.z[, paste0(thetanames[1], ".i")]
      if (tmp3.TF[ 7])
        Theta.d1 <- Thetas.z[, paste0(thetanames[1], ".d")]
    } else {  # infos.list$MM1 == 2
      if (tmp3.TF[ 4])
        Theta.a1 <- Thetas.z[, paste0(thetanames[1], ".a")]
      if (tmp3.TF[ 7])
        Theta.i1 <- Thetas.z[, paste0(thetanames[1], ".i")]
      if (tmp3.TF[10])
        Theta.d1 <- Thetas.z[, paste0(thetanames[1], ".d")]
    }




    if (infos.list$MM1 == 2) {
      ind.flip <- 2   # ifelse(flip.args, 1, 2)
      Thetas.z <- fitted(object,
                         type.fitted = paste0(thetanames[ind.flip], "s"))
      Theta.p2 <- Thetas.z[, paste0(thetanames[ind.flip], ".p")]
      Theta.a2 <- Theta.i2 <- Theta.d2 <- Theta.p2  # Needed
      if (tmp3.TF[ 4])
        Theta.a2 <- Thetas.z[, paste0(thetanames[ind.flip], ".a")]
      if (tmp3.TF[ 8])
        Theta.i2 <- Thetas.z[, paste0(thetanames[ind.flip], ".i")]
      if (tmp3.TF[11])
        Theta.d2 <- Thetas.z[, paste0(thetanames[ind.flip], ".d")]
    }  # infos.list$MM1 == 2
  }  # GAITDffs





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


      
  if (vfamily %in% GAITDffs) {



      

    baseparams.argnames <- infos.list$baseparams.argnames
    if (!length(baseparams.argnames))
      stop("cannot determine any base parameter argument name")
    if (infos.list$MM1 > 2)
      stop("cannot handle MM1 > 2")
    alist <- list(  # x = xx,  # theta.p,
                  a.mix = a.mix, a.mlm = a.mlm,
                  i.mix = i.mix, i.mlm = i.mlm,
                  d.mix = d.mix, d.mlm = d.mlm,
                  truncate = truncate, max.support = max.support,
                  byrow.aid = FALSE,  # Because of reconstruction
                  pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
                  pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
                  pdip.mix = pdip.mix, pdip.mlm = pdip.mlm)
    alist[[paste0(baseparams.argnames[1], ".p")]] <- Theta.p1
    alist[[paste0(baseparams.argnames[1], ".a")]] <- Theta.a1
    alist[[paste0(baseparams.argnames[1], ".i")]] <- Theta.i1
    alist[[paste0(baseparams.argnames[1], ".d")]] <- Theta.d1
    if (infos.list$MM1 == 2) {
      alist[[paste0(baseparams.argnames[2], ".p")]] <- Theta.p2
      alist[[paste0(baseparams.argnames[2], ".a")]] <- Theta.a2
      alist[[paste0(baseparams.argnames[2], ".i")]] <- Theta.i2
      alist[[paste0(baseparams.argnames[2], ".d")]] <- Theta.d2
    }
    dlist <- alist
    dlist$log <- FALSE
    dfun <- paste0("dgaitd", infos.list$parent.name[2])

    for (i in seq_along(At)) {
      dlist$x <- At[i]
      pmat[, i] <- do.call(dfun, dlist)  # i + lowsup - 1L
    }
  }  # vfamily %in% GAITDffs




      







  if (vfamily == "genpoisson0") {
    for (i in At)
      pmat[, i + 1L] <- dgenpois0(i, theta = Param.mat[, 1],
                                  lambda = Param.mat[, 2])
  }
  if (vfamily == "genpoisson1") {
    for (i in At)
      pmat[, i + 1L] <- dgenpois1(i, meanpar = Param.mat[, 1],
                                  dispind = Param.mat[, 2])
  }
  if (vfamily == "genpoisson2") {
    for (i in At)
      pmat[, i + 1L] <- dgenpois2(i, meanpar = Param.mat[, 1],
                                  disppar = Param.mat[, 2])
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
      pmat[, i + 1L] <- dgaitdnbinom(i, Param.mat[, 2],
                                     munb.p = Param.mat[, 1])
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
      pmat[, i + 1L] <- dgaitdbinom(i, size, Param.mat[, 1], truncate = 0)
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



















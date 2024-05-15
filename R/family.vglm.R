# These functions are
# Copyright (C) 1998-2024 T.W. Yee, University of Auckland.
# All rights reserved.




if (FALSE)
family.vglm <- function(object, ...)
    object$vfamily


if (FALSE)
print.vfamily <- function(x, ...) {
  f <- x$vfamily
  if (is.null(f))
    stop("not a VGAM family function")

  nn <- x$blurb
  if (is.null(nn))
    invisible(return(x))

  cat("Family: ", f[1], "\n")
  if (length(f)>1)
    cat("Classes:", paste(f, collapse = ", "), "\n")
  cat("\n")

  for (ii in seq_along(nn))
    cat(nn[ii])
  cat("\n")
  invisible(return(x))
}















GHfun <- function(n) {
  vals <- sqrt((1:(n - 1)) / 2)
  mat <- matrix(0, n, n)
  mat[col(mat) == row(mat) + 1] <- vals
  mat[col(mat) == row(mat) - 1] <- vals
  ans <- eigen(mat, symmetric = TRUE)
  list(nodes   = rev(ans$values),
       weights = sqrt(pi) * rev(ans$vectors[1, ]^2))
}  # GHfun








 N1binomial <-
  function(lmean = "identitylink",
           lsd = "loglink",
           lvar = "loglink",
           lprob = "logitlink",
           lapar = "rhobitlink",
           zero = c(if (var.arg) "var" else "sd",
                    "apar"),
           nnodes = 20,  # GH nodes
           copula = "gaussian",
           var.arg = FALSE,
           imethod = 1,
           isd = NULL,
           iprob = NULL, iapar = NULL) {


  stopifnot(is.numeric(nnodes),
            length(nnodes) == 1,
            round(nnodes) == nnodes,
            nnodes >= 5)

  copula <- match.arg(copula, c("gaussian"))[1]
  isdev <- isd

  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsdev <- as.list(substitute(lsd))
  esdev <- link2list(lsdev)
  lsdev <- attr(esdev, "function.name")

  lvare <- as.list(substitute(lvar))
  evare <- link2list(lvare)
  lvare <- attr(evare, "function.name")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lapar <- as.list(substitute(lapar))
  eapar <- link2list(lapar)
  lapar <- attr(eapar, "function.name")

  if (!is.Numeric(imethod, length.arg = 1,
      integer.valued = TRUE, positive = TRUE) ||
      imethod > 4)
      stop("arg 'imethod' is not 1, 2, 3 or 4")

  if (!is.logical(var.arg) ||
      length(var.arg) != 1)
    stop("arg 'var.arg' must be a single logical")
  if (var.arg)
    stop("currently 'var.arg' must be FALSE")

  new("vglmff",
  blurb = c("Univariate normal and binomial copula\n\n",
            "Links:    ",
    namesof("mean", lmean, emean, tag = TRUE), "; ",
    if (var.arg)
    namesof("var",  lvare, evare, tag = TRUE) else
    namesof("sd" ,  lsdev, esdev, tag = TRUE), "; ",
    namesof("prob", lprob, eprob, tag = TRUE), "; ",
    namesof("apar", lapar, eapar, tag = TRUE), "\n",
    if (var.arg) "Variance: var" else
                 "Variance: sd^2"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints,
       x = x, .zero , M = M, M1 = 4,
       predictors.names = predictors.names)
  }),
  list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 4,
         Q1 = 2,
         copula = copula,
         dpqrfun = "N1binom",  # No pq.
         expected = TRUE,
         hadof = FALSE,
         imethod = .imethod ,
         multipleResponses = FALSE,
         parameters.names = c("mean",
             if ( .var.arg ) "var" else "sd",
             "prob", "apar"),
         var.arg = .var.arg ,
         zero = .zero )
  },
  list( .zero = zero, .copula = copula,
        .imethod = imethod,
        .var.arg = var.arg))),

  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 2,
              ncol.y.max = 2,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (NCOL(y) != 2)
      stop("response does not have 2 columns")
    if (!all(w == 1))
      stop("all prior weights must be unity")
    if (!all(y[, 2] %in% 0:1))
      stop("2nd column of y must comprise 0s and 1s")
    if (abs(mean(y[, 2]) - 0.5) >= 0.5)
      stop("0 < mean(y[, 2]) < 1 is needed")

    ncoly <- ncol(y)
    M1 <- 4
    Q1 <- 2  # Number of responses is ncoly / Q1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    extra$Q1 <- Q1
    extra$nYs <- nYs <- ncoly / Q1  # Number of responses
    M <- M1 * nYs

    mynames1 <- param.names("mean", nYs, skip1 = TRUE)
    mynames2 <- param.names(if ( .var.arg ) "var" else "sd",
                            nYs, skip1 = TRUE)
    mynames3 <- param.names("prob", nYs, skip1 = TRUE)
    mynames4 <- param.names("apar", nYs, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .lmean , .emean , tag = FALSE),
          if ( .var.arg )
          namesof(mynames2, .lvare , .evare , tag = FALSE) else
          namesof(mynames2, .lsdev , .esdev , tag = FALSE),
          namesof(mynames3, .lprob , .eprob , tag = FALSE),
          namesof(mynames4, .lapar , .eapar , tag = FALSE))

    predictors.names <- predictors.names[interleave.VGAM(M,
                                         M1 = M1)]
    extra$predictors.names <- predictors.names

    GH.info <- GHfun( .nnodes )  # Wrt exp(-x^2)
    gh.nodes <- GH.info$nodes * sqrt(2)  # Wrt dnorm()
    gh.wghts <- GH.info$weights / sqrt(pi)
    extra$gh.nodes <- gh.nodes
    extra$gh.wghts <- gh.wghts


    if (!length(etastart)) {
      sdev.init <- mean.init <- matrix(0, n, nYs)
      for (jay in 1:nYs) {
        jfit <- lm.wfit(x, y[, jay], w[, jay])
        mean.init[, jay] <- if ( .lmean == "loglink")
                            pmax(1/1024, y[, jay]) else
          if ( .imethod == 1) median(y[, jay]) else
          if ( .imethod == 2)
            weighted.mean(y[, jay], w = w[, jay]) else
          if ( .imethod == 3)
            weighted.mean(y[, jay], w = w[, jay]) *
                          0.5 + y[, jay] * 0.5 else
                          mean(jfit$fitted)

        sdev.init[, jay] <-
          if ( .imethod == 1) {
            sqrt( sum(w[, jay] *
                (y[, jay] - mean.init[, jay])^2) / sum(w[, jay]) )
          } else if ( .imethod == 2) {
            if (jfit$df.resid > 0)
              sqrt( sum(w[, jay] * jfit$resid^2)
                    / jfit$df.resid ) else
              sqrt( sum(w[, jay] * jfit$resid^2)
                    / sum(w[, jay]) )
          } else if ( .imethod == 3) {
            sqrt( sum(w[, jay] *
                  (y[, jay] - mean.init[, jay])^2)
                 / sum(w[, jay]) )
          } else {
            sqrt( sum(w[, jay] * abs(y[, jay] -
                  mean.init[, jay]))
                 / sum(w[, jay]) )
          }

        if (any(sdev.init[, jay] <= sqrt( .Machine$double.eps ) ))
          sdev.init[, jay] <- 1.01

      }


      if (length( .isdev )) {
        sdev.init <- matrix( .isdev , n, nYs, byrow = TRUE)
      }
      prob.init <- if (length( .iprob )) {
        matrix( .iprob , n, nYs, byrow = TRUE)
      } else {
        rep(mean(y[, 2]), n)
      }
      apar.init <- if (length( .iapar )) {
        matrix( .iapar , n, nYs, byrow = TRUE)
      } else {
        rep(0.1, n)
      }

      etastart <-
        cbind(theta2eta(mean.init,   .lmean , .emean ),
              if ( .var.arg )
              theta2eta(sdev.init^2, .lvare , .evare ) else
              theta2eta(sdev.init  , .lsdev , .esdev ),
              theta2eta(prob.init  , .lprob , .eprob ),
              theta2eta(apar.init  , .lapar , .eapar ))

      colnames(etastart) <- predictors.names
    }
  }),
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .lprob = lprob, .lapar = lapar,
        .eprob = eprob, .eapar = eapar,
        .lmean = lmean,
        .emean = emean, .iprob = iprob,
        .isdev = isdev, .iapar = iapar,
        .nnodes = nnodes,
        .var.arg = var.arg, .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    if (extra$nYs > 1)
      stop("multiple ys unallowed")
    cbind(eta2theta(eta[, 1], .lmean , .emean ),
          eta2theta(eta[, 3], .lprob , .eprob ))
  },
  list( .lprob = lprob, .lmean = lmean,
        .eprob = eprob, .emean = emean ))),

  last = eval(substitute(expression({
    M1 <- extra$M1
    nYs <- extra$nYs
    temp.names <- c(mynames1, mynames2, mynames3, mynames4)
    temp.names <- temp.names[interleave.VGAM(M1 * nYs,
                                             M1 = M1)]
    misc$link <- rep_len( .lmean , M1 * nYs)
    misc$earg <- vector("list", M1 * nYs)
    names(misc$link) <- names(misc$earg) <- temp.names
    for (ii in 1:nYs) {
      misc$link[ M1*ii-3 ] <- ( .lmean )
      misc$link[ M1*ii-2 ] <- if ( .var.arg )
                   .lvare else ( .lsdev )
      misc$link[ M1*ii-1 ] <- ( .lprob )
      misc$link[ M1*ii   ] <- ( .lapar )
      misc$earg[[M1*ii-3]] <- ( .emean )
      misc$earg[[M1*ii-2]] <- if ( .var.arg )
                   .evare else ( .esdev )
      misc$earg[[M1*ii-1]] <- ( .eprob )
      misc$earg[[M1*ii  ]] <- ( .eapar )
    }
  }),
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .lprob = lprob, .lapar = lapar,
        .eprob = eprob, .eapar = eapar,
        .lmean = lmean,
        .emean = emean,
        .var.arg = var.arg))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    mymu <- eta2theta(eta[, 1], .lmean , .emean )
    sdev <- if ( .var.arg ) {
      sqrt(eta2theta(eta[, 2], .lvare , .evare ))
    } else {
      eta2theta(eta[, 2], .lsdev , .esdev )
    }
    prob <- eta2theta(eta[, 3], .lprob , .eprob )
    apar <- eta2theta(eta[, 4], .lapar , .eapar )
    if (residuals) {
      stop("loglik resids not implemented")
    } else {
        ll.elts <- c(w) *
            dN1binom(y[, 1], y[, 2], mymu, sdev,
                     prob, apar, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .lprob = lprob, .lapar = lapar,
        .eprob = eprob, .eapar = eapar,
        .lmean = lmean,
        .emean = emean,
        .var.arg = var.arg ))),
  vfamily = c("N1binomial"),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    mymu <- eta2theta(eta[, 1], .lmean , .emean )
    sdev <- Varm <- 111
    if ( .var.arg ) {
      Varm <- eta2theta(eta[, 2], .lvare , .evare )
    } else {
      sdev <- eta2theta(eta[, 2], .lsdev , .esdev )
    }
    prob <- eta2theta(eta[, 3], .lprob , .eprob )
    apar <- eta2theta(eta[, 4], .lapar , .eapar )
    okay1 <-
        all(is.finite(mymu)) &&
        all(is.finite(sdev)) && all(0 < sdev) &&
        all(is.finite(Varm)) && all(0 < Varm) &&
        all(is.finite(prob)) && all(0 < prob) &&
                                all(prob < 1) &&
        all(is.finite(apar)) && all(abs(apar) < 1)
    okay1
  },
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .lprob = lprob, .lapar = lapar,
        .eprob = eprob, .eapar = eapar,
        .lmean = lmean,
        .emean = emean,
        .var.arg = var.arg ))),




  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    if (ncol((fv <- fitted(object))) != 2)
      stop("fitted(object) has not got 2 cols")
    mean <- fv[, 1]
    prob <- fv[, 2]
    nn <- NROW(fv)
    eta <- predict(object)
    sdev <- if ( .var.arg ) {
      sqrt(eta2theta(eta[, 2], .lvare , .evare ))
    } else {
      eta2theta(eta[, 2], .lvare , .evare )
    }
    apar <- eta2theta(eta[, 4], .lapar , .eapar )
    aaa <- array(c(
       rN1binom(nn * nsim, mean, sdev, prob, apar)),
       dim = c(nn, nsim, 2))
    aaa <- aperm(aaa, c(1, 3, 2))
    attr(aaa, "Verbatim") <- TRUE  # Removed later
    aaa
  },
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .lprob = lprob, .lapar = lapar,
        .eprob = eprob, .eapar = eapar,
        .lmean = lmean,
        .emean = emean,
        .var.arg = var.arg ))),




  deriv = eval(substitute(expression({
    mymu <- eta2theta(eta[, 1], .lmean , .emean )
    sdev <- if ( .var.arg ) {
      sqrt(eta2theta(eta[, 2], .lvare , .evare ))
    } else {
      eta2theta(eta[, 2], .lsdev , .esdev )
    }
    muu2 <- eta2theta(eta[, 3], .lprob , .eprob )
    apar <- eta2theta(eta[, 4], .lapar , .eapar )

    zedd <- (y[, 1] - mymu) / sdev
    dl.dmean <- zedd / sdev
    dl.dvarr <- dl.dsdev <- NULL  # 4 cbind() below
    if ( .var.arg ) {
      dl.dvarr <- -0.5 / Varm + 0.5 * zedd^2
    } else {
      dl.dsdev <- (zedd^2 - 1) / sdev
    }

    dmean.deta <- dtheta.deta(mymu, .lmean , .emean )
    dvarr.deta <- dsdev.deta <- NULL  # 4 cbind() below
    if ( .var.arg ) {
      dvarr.deta <- dtheta.deta(Varm, .lvare , .evare )
    } else {
      dsdev.deta <- dtheta.deta(sdev, .lsdev , .esdev )
    }
    dmuu2.deta <- dtheta.deta(muu2, .lprob , .eprob )
    dapar.deta <- dtheta.deta(apar, .lapar , .eapar )

    Prob <- pfun.N1b(zedd, muu2, apar)  # Delta
    dProb.dmuu2 <- pfun.N1b(zedd, muu2, apar,
                            Integrand = FALSE,
                            1, "a")
    dProb.dzedd <- pfun.N1b(zedd, muu2, apar,
                            Integrand = FALSE,
                            1, "b")
    dProb.dapar <- pfun.N1b(zedd, muu2, apar,
                            Integrand = FALSE,
                            1, "apar")
    dzedd.dmean <- (-1) / sdev
    dzedd.dsdev <- (-zedd) / sdev
    dProb.dmean <- dProb.dzedd * dzedd.dmean
    dProb.dsdev <- dProb.dzedd * dzedd.dsdev
    if ( .var.arg )
      warning("yettodo: compute dl.dvare here")
    tmpID <- (y[, 2] == 1) /      Prob -
             (y[, 2] == 0) / (1 - Prob)
    dl.dmean <- dl.dmean + tmpID * dProb.dmean
    dl.dsdev <- dl.dsdev + tmpID * dProb.dsdev
    dl.dmuu2 <- tmpID * dProb.dmuu2
    dl.dapar <- tmpID * dProb.dapar

    dthetas.detas  <-  # Useful for wz too
        cbind(dmean.deta,
              dsdev.deta, dvarr.deta,  # Only 1
              dmuu2.deta, dapar.deta)
    ans <- c(w) * dthetas.detas *
      cbind(dl.dmean,
            dl.dvarr, dl.dsdev,  # Only 1
            dl.dmuu2, dl.dapar)


    ans
  }),
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .lprob = lprob, .lapar = lapar,
        .eprob = eprob, .eapar = eapar,
        .lmean = lmean,
        .emean = emean,
        .var.arg = var.arg ))),
  weight = eval(substitute(expression({
    dimmM <- dimm(M)
    indlist <- iam(NA, NA, M, both = TRUE)
    ind.rows <- indlist$row  # Length == dimmM
    ind.cols <- indlist$col
    gh.nodes <- extra$gh.nodes
    gh.wghts <- extra$gh.wghts  # Normalized
    mmm1 <- 2 * qnorm(muu2) * apar / (
            (1 + apar^2))
    sss1 <- sqrt((1 - apar^2) / (1 + apar^2))


    wz <- matrix(0, n, dimmM)
    for (ii in seq( .nnodes )) {
      node.use <- mmm1 + sss1 * gh.nodes[ii]
      node.use <- mymu + sdev * node.use
      wz <- wz + gh.wghts[ii] *
        eimjk.N1b(node.use,  # vector
                  Integrand = TRUE,
                  apar = apar, mu2 = muu2,
                  mean = mymu, sdev = sdev,
                  jay = NULL, kay = NULL)  # All
    }
    
      


    bigconst <-
      exp(-qnorm(muu2) / (1 + apar^2)) *
      sqrt((1 - apar^2) / (1 + apar^2)) / (
       (2 * pi))
    wz <- wz * bigconst


    ned2l.dmean2 <- 1 / sdev^2
    if ( .var.arg ) {
      ned2l.dvarr2 <- 0.5 / Varm^2
    } else {
      ned2l.dsdev2 <- 2 / sdev^2
    }
    wz[, iam(1, 1, M)] <- 
    wz[, iam(1, 1, M)] + ned2l.dmean2
    wz[, iam(2, 2, M)] <- 
    wz[, iam(2, 2, M)] + (if ( .var.arg )
      ned2l.dvarr2 else ned2l.dsdev2)


    wz <- wz * dthetas.detas[, ind.rows] *
               dthetas.detas[, ind.cols]
    colnames(wz) <- NULL  # Tidy up

    c(w) * wz
  }),
  list( .var.arg = var.arg,
        .nnodes = nnodes ))))
}  #  N1binomial()







 pfun.N1b <-
    function(b, a, apar = 0,
             deriv = 0,
             which = "b",
             Integrand = FALSE  # Include dnorm()
             ) {


  if (deriv == 0) return(pnorm(Qfun.N1b(b, a, apar)))
  if (deriv == 1) return(
    (if (Integrand) 1 else
      dnorm(Qfun.N1b(b, a, apar))) *
    Qfun.N1b(b, a, apar, deriv, which))
  stop("arg 'deriv' not in 0:1")
}  # pfun.N1b





 Qfun.N1b <-
    function(b, a, apar = 0,
             deriv = 0,
             which = "b") {
  if (deriv == 0) return(
    (qnorm(a) - apar * b) / sqrt(1 - apar^2))
  if (which == "a") {
    return(1 / (dnorm(qnorm(a)) *
                sqrt(1 - apar^2)))
  }
  if (which == "b") {
    return(-apar / sqrt(1 - apar^2))
  }
  if (which == "apar") {
    Q <- Qfun.N1b(b, a, apar, deriv = 0)
    return((apar * Q / sqrt(1 - apar^2) - b)
           / sqrt(1 - apar^2))
  }
  stop("args 'deriv' and 'which' unrecognized")
}  # Qfun.N1b







 dN1binom <-
    function(x1, x2,
             mean = 0, sd = 1,  # for x1
             prob,  # == E(x2)==mu2,  # size,
             apar = 0,
             copula = "gaussian",
             log = FALSE) {
  copula <- match.arg(copula, c("gaussian"))[1]
  if (!is.logical(log.arg <- log) ||
      length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  check1 <- all(x2 %in% 0:1, na.rm = TRUE)
  if (!check1)
    stop("arg 'x2' must have 1 or 0 values only")

  L <- max(length(x1),   length(x2),
           length(mean), length(sd),
           length(prob), length(apar))
  if (length(x1)   != L) x1   <- rep_len(x1,   L)
  if (length(x2)   != L) x2   <- rep_len(x2,   L)
  if (length(mean) != L) mean <- rep_len(mean, L)
  if (length(sd)   != L) sd   <- rep_len(sd,   L)
  if (length(prob) != L) prob <- rep_len(prob, L)
  if (length(apar) != L) apar <- rep_len(apar, L)

  logdensity <- dnorm(x1, mean, sd, log = TRUE)
  Prob <- pfun.N1b((x1 - mean) / sd, prob, apar)
  x20 <- x2 == 0
  logdensity[x20] <-
  logdensity[x20] + log1p(-Prob[x20])
  x21 <- x2 == 1
  logdensity[x21] <-
  logdensity[x21] + log(Prob[x21])
  if (log.arg) logdensity else exp(logdensity)
}  # dN1binom



 rN1binom <-
    function(n,
             mean = 0, sd = 1,
             prob,
             apar = 0,
             copula = "gaussian") {
        
  copula <- match.arg(copula, c("gaussian"))[1]
  use.n <-
    if ((length.n <- length(n)) > 1) length.n else
    if (!is.Numeric(n, integer.valued = TRUE,
                    length.arg = 1, positive = TRUE))
      stop("bad input for argument 'n'") else n
  L <- use.n
  if (length(mean) != L) mean <- rep_len(mean, L)
  if (length(sd)   != L) sd   <- rep_len(sd,   L)
  if (length(prob) != L) prob <- rep_len(prob, L)
  if (length(apar) != L) apar <- rep_len(apar, L)

  x1 <- rnorm(use.n, mean, sd)  # 1st marginal same
  pdf0 <- dN1binom(x1, 0, mean, sd, prob, apar)
  pdf1 <- dN1binom(x1, 1, mean, sd, prob, apar)

  Prob <- pdf1 / (pdf0 + pdf1)
  x2 <- rbinom(use.n, size = 1, Prob)
  cbind(X1 = x1, X2 = x2)
}  # rN1binom










 eimjk.N1b <-
  function(y1,  # abscissae
           mean = 0, sdev = 1,
           mu2,  # aka a or prob
           apar = 0,  # aka alpha
           Integrand = TRUE,  # Omit dnorm()
           jay = NULL,    # Used if which = "b"
           kay = NULL) {  # %in% 1:2 each

    zedd <- (y1 - mean) / sdev
    muu2 <- mu2
    dzedd.dmean <- (-1) / sdev
    dzedd.dsdev <- (-zedd) / sdev
    Prob <- pfun.N1b(zedd, a = muu2,
                     apar = apar)  # Delta
    dProb.dmuu2 <- pfun.N1b(zedd, muu2,
                            Integrand = TRUE,
                            apar, 1, "a")
    dProb.dzedd <- pfun.N1b(zedd, muu2, apar,
                            Integrand = TRUE,
                            1, "b")
    dProb.dapar <- pfun.N1b(zedd, muu2, apar,
                            Integrand = TRUE,
                            1, "apar")
    dProb.dmean <- dProb.dzedd * dzedd.dmean
    dProb.dsdev <- dProb.dzedd * dzedd.dsdev


    dmean.deta <- dsdev.deta <- 1  # Artificial here
    dmuu2.deta <- dapar.deta <- 1  # Artificial here

    ned2l.dmeanmuu2 <- dProb.dmean * dProb.dmuu2
    ned2l.dmeanapar <- dProb.dmean * dProb.dapar
    ned2l.dsdevmuu2 <- dProb.dsdev * dProb.dmuu2
    ned2l.dsdevapar <- dProb.dsdev * dProb.dapar
    ned2l.dmuu22 <- dProb.dmuu2^2
    ned2l.dapar2 <- dProb.dapar^2
    ned2l.dmuu2apar <- dProb.dmuu2 * dProb.dapar
    M <- 4
    n <- length(y1)  # Artificial
    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M)] <- (dProb.dmean * dmean.deta)^2
    wz[, iam(2, 2, M)] <- (dProb.dsdev * dsdev.deta)^2
    wz[, iam(3, 3, M)] <- ned2l.dmuu22 * dmuu2.deta^2
    wz[, iam(4, 4, M)] <- ned2l.dapar2 * dapar.deta^2
    wz[, iam(1, 2, M)] <- dProb.dmean * dmean.deta *
                          dProb.dsdev * dsdev.deta
    wz[, iam(1, 3, M)] <- ned2l.dmeanmuu2 *
                          dmean.deta * dmuu2.deta
    wz[, iam(1, 4, M)] <- ned2l.dmeanapar *
                          dmean.deta * dapar.deta
    wz[, iam(2, 3, M)] <- ned2l.dsdevmuu2 *
                          dsdev.deta * dmuu2.deta
    wz[, iam(2, 4, M)] <- ned2l.dsdevapar *
                          dsdev.deta * dapar.deta
    wz[, iam(3, 4, M)] <- ned2l.dmuu2apar *
                          dmuu2.deta * dapar.deta


  Delta <- Prob
  const3 <- (if (Integrand) 1 else
            dnorm(y1, mean, sdev)) / (
            (Delta * (1 - Delta)))
  const3 * (if (length(jay) && length(kay))
  wz[, iam(jay, kay, M)] else wz)
}  # eimjk.N1b










 N1poisson <-
  function(lmean = "identitylink",
           lsd = "loglink",
           lvar = "loglink",
           llambda = "loglink",
           lapar = "rhobitlink",
           zero = c(if (var.arg) "var" else "sd",
                    "apar"),
           doff = 5,  # -log1p(10),  # ok: <0
           nnodes = 20,  # GH nodes
           copula = "gaussian",
           var.arg = FALSE,
           imethod = 1,
           isd = NULL,
           ilambda = NULL, iapar = NULL) {


  stopifnot(is.numeric(nnodes),
            length(nnodes) == 1,
            round(nnodes) == nnodes,
            nnodes >= 5)

  copula <- match.arg(copula, c("gaussian"))[1]
  isdev <- isd

  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsdev <- as.list(substitute(lsd))
  esdev <- link2list(lsdev)
  lsdev <- attr(esdev, "function.name")

  lvare <- as.list(substitute(lvar))
  evare <- link2list(lvare)
  lvare <- attr(evare, "function.name")

  llamb <- as.list(substitute(llambda))
  elamb <- link2list(llamb)
  llamb <- attr(elamb, "function.name")

  lapar <- as.list(substitute(lapar))
  eapar <- link2list(lapar)
  lapar <- attr(eapar, "function.name")

  if (!is.Numeric(imethod, length.arg = 1,
      integer.valued = TRUE, positive = TRUE) ||
      imethod > 4)

    if (!is.Numeric(doff,  # positive = TRUE), 
                    length.arg = 1) ||
        doff == 0)
    stop("arg 'doff' is unsuitable")

  if (!is.logical(var.arg) ||
      length(var.arg) != 1)
    stop("arg 'var.arg' must be a single logical")
  if (var.arg)
    stop("currently 'var.arg' must be FALSE")

  new("vglmff",
  blurb = c("Univariate normal and Poisson copula\n\n",
            "Links:    ",
    namesof("mean",   lmean, emean, tag = TRUE), "; ",
    if (var.arg)
    namesof("var",    lvare, evare, tag = TRUE) else
    namesof("sd" ,    lsdev, esdev, tag = TRUE), "; ",
    namesof("lambda", llamb, elamb, tag = TRUE), "; ",
    namesof("apar",   lapar, eapar, tag = TRUE), "\n",
    if (var.arg) "Variance: var" else
                 "Variance: sd^2"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints,
       x = x, .zero , M = M, M1 = 4,
       predictors.names = predictors.names)
  }),
  list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 4,
         Q1 = 2,
         copula = copula,
         doff = .doff ,
         dpqrfun = "N1binom",  # No pq.
         expected = TRUE,
         hadof = FALSE,
         imethod = .imethod ,
         multipleResponses = FALSE,
         parameters.names = c("mean",
             if ( .var.arg ) "var" else "sd",
             "lambda", "apar"),
         var.arg = .var.arg ,
         zero = .zero )
  },
  list( .zero = zero, .copula = copula,
        .imethod = imethod, .doff = doff,
        .var.arg = var.arg))),

  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 2,
              ncol.y.max = 2,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (NCOL(y) != 2)
      stop("response does not have 2 columns")
    if (!all(w == 1))
      stop("all prior weights must be unity")
    if (!all(y[, 2] == round(y[, 2])))
      stop("2nd column of y must comprise 0s and 1s")
    if (min(y[, 2]) < 0)
      stop("some y[, 2] < 0 detected")

    ncoly <- ncol(y)
    M1 <- 4
    Q1 <- 2  # Number of responses is ncoly / Q1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    extra$Q1 <- Q1
    extra$nYs <- nYs <- ncoly / Q1  # Number of responses
    M <- M1 * nYs

    mynames1 <- param.names("mean", nYs, skip1 = TRUE)
    mynames2 <- param.names(if ( .var.arg ) "var" else "sd",
                            nYs, skip1 = TRUE)
    mynames3 <- param.names("lambda", nYs, skip1 = TRUE)
    mynames4 <- param.names("apar", nYs, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .lmean , .emean , tag = FALSE),
          if ( .var.arg )
          namesof(mynames2, .lvare , .evare , tag = FALSE) else
          namesof(mynames2, .lsdev , .esdev , tag = FALSE),
          namesof(mynames3, .llamb , .elamb , tag = FALSE),
          namesof(mynames4, .lapar , .eapar , tag = FALSE))

    predictors.names <- predictors.names[interleave.VGAM(M,
                                         M1 = M1)]
    extra$predictors.names <- predictors.names

    GH.info <- GHfun( .nnodes )  # Wrt exp(-x^2)
    gh.nodes <- GH.info$nodes * sqrt(2)  # Wrt dnorm()
    gh.wghts <- GH.info$weights / sqrt(pi)
    extra$gh.nodes <- gh.nodes
    extra$gh.wghts <- gh.wghts


    if (!length(etastart)) {
      sdev.init <- mean.init <- matrix(0, n, nYs)
      for (jay in 1:nYs) {
        jfit <- lm.wfit(x, y[, jay], w[, jay])
        mean.init[, jay] <- if ( .lmean == "loglink")
                            pmax(1/1024, y[, jay]) else
          if ( .imethod == 1) median(y[, jay]) else
          if ( .imethod == 2)
            weighted.mean(y[, jay], w = w[, jay]) else
          if ( .imethod == 3)
            weighted.mean(y[, jay], w = w[, jay]) *
                          0.5 + y[, jay] * 0.5 else
                          mean(jfit$fitted)

        sdev.init[, jay] <-
          if ( .imethod == 1) {
            sqrt( sum(w[, jay] *
                (y[, jay] - mean.init[, jay])^2) / sum(w[, jay]) )
          } else if ( .imethod == 2) {
            if (jfit$df.resid > 0)
              sqrt( sum(w[, jay] * jfit$resid^2)
                    / jfit$df.resid ) else
              sqrt( sum(w[, jay] * jfit$resid^2)
                    / sum(w[, jay]) )
          } else if ( .imethod == 3) {
            sqrt( sum(w[, jay] *
                  (y[, jay] - mean.init[, jay])^2)
                 / sum(w[, jay]) )
          } else {
            sqrt( sum(w[, jay] * abs(y[, jay] -
                  mean.init[, jay]))
                 / sum(w[, jay]) )
          }

        if (any(sdev.init[, jay] <= sqrt( .Machine$double.eps ) ))
          sdev.init[, jay] <- 1.01

      }


      if (length( .isdev )) {
        sdev.init <- matrix( .isdev , n, nYs, byrow = TRUE)
      }
      lamb.init <- if (length( .ilamb )) {
        matrix( .ilamb , n, nYs, byrow = TRUE)
      } else {
        rep(mean(y[, 2]), n)
      }
      apar.init <- if (length( .iapar )) {
        matrix( .iapar , n, nYs, byrow = TRUE)
      } else {
        rep(0.05, n)
      }

      etastart <-
        cbind(theta2eta(mean.init,   .lmean , .emean ),
              if ( .var.arg )
              theta2eta(sdev.init^2, .lvare , .evare ) else
              theta2eta(sdev.init  , .lsdev , .esdev ),
              theta2eta(lamb.init  , .llamb , .elamb ),
              theta2eta(apar.init  , .lapar , .eapar ))

      colnames(etastart) <- predictors.names
    }
  }),
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .llamb = llamb, .lapar = lapar,
        .elamb = elamb, .eapar = eapar,
        .lmean = lmean,
        .emean = emean, .ilamb = ilambda,
        .isdev = isdev, .iapar = iapar,
        .doff  = doff,  .nnodes = nnodes, 
        .var.arg = var.arg, .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    if (extra$nYs > 1)
      stop("multiple ys unallowed")
    cbind(eta2theta(eta[, 1], .lmean , .emean ),
          eta2theta(eta[, 3], .llamb , .elamb ))
  },
  list( .llamb = llamb, .lmean = lmean,
        .elamb = elamb, .emean = emean ))),

  last = eval(substitute(expression({
    M1 <- extra$M1
    nYs <- extra$nYs
    temp.names <- c(mynames1, mynames2, mynames3, mynames4)
    temp.names <- temp.names[interleave.VGAM(M1 * nYs,
                                             M1 = M1)]
    misc$link <- rep_len( .lmean , M1 * nYs)
    misc$earg <- vector("list", M1 * nYs)
    names(misc$link) <- names(misc$earg) <- temp.names
    for (ii in 1:nYs) {
      misc$link[ M1*ii-3 ] <- ( .lmean )
      misc$link[ M1*ii-2 ] <- if ( .var.arg )
                   .lvare else ( .lsdev )
      misc$link[ M1*ii-1 ] <- ( .llamb )
      misc$link[ M1*ii   ] <- ( .lapar )
      misc$earg[[M1*ii-3]] <- ( .emean )
      misc$earg[[M1*ii-2]] <- if ( .var.arg )
                   .evare else ( .esdev )
      misc$earg[[M1*ii-1]] <- ( .elamb )
      misc$earg[[M1*ii  ]] <- ( .eapar )
    }
  }),
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .llamb = llamb, .lapar = lapar,
        .elamb = elamb, .eapar = eapar,
        .lmean = lmean, .doff = doff,
        .emean = emean,
        .var.arg = var.arg))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    mymu <- eta2theta(eta[, 1], .lmean , .emean )
    sdev <- if ( .var.arg ) {
      sqrt(eta2theta(eta[, 2], .lvare , .evare ))
    } else {
      eta2theta(eta[, 2], .lsdev , .esdev )
    }
    Lamb <- eta2theta(eta[, 3], .llamb , .elamb )
    apar <- eta2theta(eta[, 4], .lapar , .eapar )
    if (residuals) {
      stop("loglik resids not implemented")
    } else {
      ll.elts <- c(w) *
        dN1pois(y[, 1], y[, 2], mymu, sdev,
                Lamb, apar = apar,
                doff =  .doff , log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .llamb = llamb, .lapar = lapar,
        .elamb = elamb, .eapar = eapar,
        .lmean = lmean, .doff = doff,
        .emean = emean,
        .var.arg = var.arg ))),
  vfamily = c("N1poisson"),

  validparams = eval(substitute(function(eta,
                     y, extra = NULL) {
    mymu <- eta2theta(eta[, 1], .lmean , .emean )
    sdev <- Varm <- 111
    if ( .var.arg ) {
      Varm <- eta2theta(eta[, 2], .lvare , .evare )
    } else {
      sdev <- eta2theta(eta[, 2], .lsdev , .esdev )
    }
    Lamb <- eta2theta(eta[, 3], .llamb , .elamb )
    apar <- eta2theta(eta[, 4], .lapar , .eapar )
    okay1 <-
        all(is.finite(mymu)) &&
        all(is.finite(sdev)) && all(0 < sdev) &&
        all(is.finite(Varm)) && all(0 < Varm) &&
        all(is.finite(Lamb)) && all(0 < Lamb) &&
        all(is.finite(apar)) && all(abs(apar) < 1)
    okay1
  },
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .llamb = llamb, .lapar = lapar,
        .elamb = elamb, .eapar = eapar,
        .lmean = lmean, .doff = doff,
        .emean = emean,
        .var.arg = var.arg ))),




  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    if (ncol((fv <- fitted(object))) != 2)
      stop("fitted(object) has not got 2 cols")
    mean <- fv[, 1]
    Lamb <- fv[, 2]
    nn <- NROW(fv)
    eta <- predict(object)
    sdev <- if ( .var.arg ) {
      sqrt(eta2theta(eta[, 2], .lvare , .evare ))
    } else {
      eta2theta(eta[, 2], .lvare , .evare )
    }
    apar <- eta2theta(eta[, 4], .lapar , .eapar )
    aaa <- array(c(
       rN1pois(nn * nsim, mean, sdev, Lamb,
               doff = .doff , apar)),
       dim = c(nn, nsim, 2))
    aaa <- aperm(aaa, c(1, 3, 2))
    attr(aaa, "Verbatim") <- TRUE  # Removed later
    aaa
  },
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .llamb = llamb, .lapar = lapar,
        .elamb = elamb, .eapar = eapar,
        .lmean = lmean, .doff = doff,
        .emean = emean,
        .var.arg = var.arg ))),




  deriv = eval(substitute(expression({
    mymu <- eta2theta(eta[, 1], .lmean , .emean )
    sdev <- if ( .var.arg ) {
      sqrt(eta2theta(eta[, 2], .lvare , .evare ))
    } else {
      eta2theta(eta[, 2], .lsdev , .esdev )
    }
    Lamb <-
    muu2 <- eta2theta(eta[, 3], .llamb , .elamb )
    apar <- eta2theta(eta[, 4], .lapar , .eapar )
    doff <- ( .doff )  # Nonzero

    zedd <- (y[, 1] - mymu) / sdev
    dl.dmean <- zedd / sdev
    dl.dvarr <- dl.dsdev <- NULL  # 4 cbind() below
    if ( .var.arg ) {
      dl.dvarr <- -0.5 / Varm + 0.5 * zedd^2
    } else {
      dl.dsdev <- (zedd^2 - 1) / sdev
    }

    dmean.deta <- dtheta.deta(mymu, .lmean , .emean )
    dvarr.deta <- dsdev.deta <- NULL  # 4 cbind() below
    if ( .var.arg ) {
      dvarr.deta <- dtheta.deta(Varm, .lvare , .evare )
    } else {
      dsdev.deta <- dtheta.deta(sdev, .lsdev , .esdev )
    }
    dmuu2.deta <- dtheta.deta(muu2, .llamb , .elamb )
    dapar.deta <- dtheta.deta(apar, .lapar , .eapar )

    Dstr <- if (doff > 0)
      muu2^(2/3) / (muu2^(2/3) + abs(doff)) else
      log1p(muu2) / (abs(doff) + log1p(muu2))
    dDstr.dmuu2 <- if (doff > 0)
      (2 / 3) * abs(doff) / (muu2^(1/3) *
      (muu2^(2/3) + abs(doff))^2) else
      abs(doff) / ((1 + muu2) * (
      abs(doff) + log1p(muu2))^2)

    Prob <- pfun.N1p(zedd,
                     a = Dstr,  # Was muu2,
                     apar)  # Delta
    dProb.dmuu2 <- pfun.N1p(zedd,
                            a = Dstr,  # Was muu2,
                            Integrand = FALSE,
                            apar, 1, "a") *
                   dDstr.dmuu2  # New
    dProb.dzedd <- pfun.N1p(zedd,
                            a = Dstr,  # Was muu2,
                            Integrand = FALSE,
                            apar, 1, "b")
    dProb.dapar <- pfun.N1p(zedd,
                            a = Dstr,  # Was muu2,
                            Integrand = FALSE,
                            apar, 1, "apar")
    dzedd.dmean <- (-1) / sdev
    dzedd.dsdev <- (-zedd) / sdev
    dProb.dmean <- dProb.dzedd * dzedd.dmean
    dProb.dsdev <- dProb.dzedd * dzedd.dsdev
    if ( .var.arg )
      warning("yettodo: compute dl.dvare here")
    tmpodds <- abs(doff) * Prob / (1 - Prob)
    dl2.dProb <- if (doff > 0)
      (1.5 / (1 - Prob)) *  # Change
      (y[, 2] / Prob - abs(doff)^1.5 *
       sqrt(Prob) / (1 - Prob)^1.5) else
      abs(doff) *
      (y[, 2] / expm1(tmpodds) - 1) *
       exp(tmpodds) / (1 - Prob)^2
    dl.dmean <- dl.dmean + dl2.dProb * dProb.dmean
    dl.dsdev <- dl.dsdev + dl2.dProb * dProb.dsdev
    dl.dmuu2 <- dl2.dProb * dProb.dmuu2
    dl.dapar <- dl2.dProb * dProb.dapar

    dthetas.detas  <-  # Useful for wz too
        cbind(dmean.deta,
              dsdev.deta, dvarr.deta,  # Only 1
              dmuu2.deta, dapar.deta)
    ans <- c(w) * dthetas.detas *
      cbind(dl.dmean,
            dl.dvarr, dl.dsdev,  # Only 1
            dl.dmuu2, dl.dapar)


    ans
  }),
  list( .lsdev = lsdev, .lvare = lvare,
        .esdev = esdev, .evare = evare,
        .llamb = llamb, .lapar = lapar,
        .elamb = elamb, .eapar = eapar,
        .lmean = lmean, .doff = doff,
        .emean = emean,
        .var.arg = var.arg ))),
  weight = eval(substitute(expression({
    dimmM <- dimm(M)
    indlist <- iam(NA, NA, M, both = TRUE)
    ind.rows <- indlist$row  # Length == dimmM
    ind.cols <- indlist$col
    gh.nodes <- extra$gh.nodes
    gh.wghts <- extra$gh.wghts  # Normalized
    mmm1 <- 2 * qnorm(Dstr) *  # Change (was muu2)
            apar / (1 + apar^2)
    sss1 <- sqrt((1 - apar^2) / (1 + apar^2))


    wz <- matrix(0, n, dimmM)
    for (ii in seq( .nnodes )) {
      node.use <- mmm1 + sss1 * gh.nodes[ii]
      node.use <- mymu + sdev * node.use
      wz <- wz + gh.wghts[ii] *
        eimjk.N1p(node.use,  # vector
                  Integrand = TRUE,
                  apar = apar, mu2 = muu2,
                  mean = mymu, sdev = sdev,
                  doff = .doff ,  # May be negative
                  jay = NULL, kay = NULL)  # All
    }
    
      


    bigconst <-
      exp(-qnorm(Dstr) / (  # Change
      1 + apar^2)) *
      sqrt((1 - apar^2) / (1 + apar^2)) / (
       (2 * pi))
    wz <- wz * bigconst


    ned2l.dmean2 <- 1 / sdev^2
    if ( .var.arg ) {
      ned2l.dvarr2 <- 0.5 / Varm^2
    } else {
      ned2l.dsdev2 <- 2 / sdev^2
    }
    wz[, iam(1, 1, M)] <- 
    wz[, iam(1, 1, M)] + ned2l.dmean2
    wz[, iam(2, 2, M)] <- 
    wz[, iam(2, 2, M)] + (if ( .var.arg )
      ned2l.dvarr2 else ned2l.dsdev2)


    wz <- wz * dthetas.detas[, ind.rows] *
               dthetas.detas[, ind.cols]
    colnames(wz) <- NULL  # Tidy up

    c(w) * wz
  }),
  list( .var.arg = var.arg, .doff = doff,
        .nnodes = nnodes ))))
}  #  N1poisson()







 pfun.N1p <-
    function(b, a, apar = 0,
             deriv = 0,
             which = "b",
             Integrand = FALSE  # Include dnorm()
             ) {


  if (deriv == 0) return(pnorm(Qfun.N1p(b, a, apar)))
  if (deriv == 1) return(
    (if (Integrand) 1 else
      dnorm(Qfun.N1p(b, a, apar))) *
    Qfun.N1p(b, a, apar, deriv, which))
  stop("arg 'deriv' not in 0:1")
}  # pfun.N1p





 Qfun.N1p <-
    function(b, a, apar = 0,
             deriv = 0,
             which = "b") {
  if (deriv == 0) return(
    (qnorm(a) - apar * b) / sqrt(1 - apar^2))
  if (which == "a") {
    return(1 / (dnorm(qnorm(a)) *
                sqrt(1 - apar^2)))
  }
  if (which == "b") {
    return(-apar / sqrt(1 - apar^2))
  }
  if (which == "apar") {
    Q <- Qfun.N1p(b, a, apar, deriv = 0)
    return((apar * Q / sqrt(1 - apar^2) - b)
           / sqrt(1 - apar^2))
  }
  stop("args 'deriv' and 'which' unrecognized")
}  # Qfun.N1p







 dN1pois <-
    function(x1, x2,
             mean = 0, sd = 1,  # for x1
             lambda,  # == E(x2)==mu2,  # size,
             apar = 0,
             doff = 5,  # -log1p(10),
             copula = "gaussian",
             log = FALSE) {
  Lamb <- lambda
  copula <- match.arg(copula, c("gaussian"))[1]
  if (!is.logical(log.arg <- log) ||
      length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)
  if (!is.numeric(doff) || length(doff) != 1 ||
      doff == 0)
    stop("bad input for argument 'doff'")


  L <- max(length(x1),   length(x2),
           length(mean), length(sd),
           length(Lamb), length(apar))
  if (length(x1)   != L) x1   <- rep_len(x1,   L)
  if (length(x2)   != L) x2   <- rep_len(x2,   L)
  if (length(mean) != L) mean <- rep_len(mean, L)
  if (length(sd)   != L) sd   <- rep_len(sd,   L)
  if (length(Lamb) != L) Lamb <- rep_len(Lamb, L)
  if (length(apar) != L) apar <- rep_len(apar, L)

  logdensity <- dnorm(x1, mean, sd, log = TRUE)
  txlamb <- if (doff > 0)
     Lamb^(2/3) / (abs(doff) + Lamb^(2/3)) else
     log1p(Lamb) / (abs(doff) + log1p(Lamb))
  Delt <- pfun.N1p((x1 - mean) / sd, txlamb, apar)
  use.l <- if (doff > 0)  # use.lambda
     (abs(doff) * Delt / (1 - Delt))^1.5 else
     expm1(abs(doff) * Delt / (1 - Delt))
  logdensity <- logdensity +
                dpois(x2, use.l, log = TRUE)
  if (log.arg) logdensity else exp(logdensity)
}  # dN1pois



 rN1pois <-
    function(n,
             mean = 0, sd = 1,
             lambda,
             apar = 0,
             doff = 5,  # -log1p(10),
             copula = "gaussian") {

  Lamb <- lambda
  copula <- match.arg(copula, c("gaussian"))[1]
  use.n <-
    if ((length.n <- length(n)) > 1) length.n else
    if (!is.Numeric(n, integer.valued = TRUE,
                    length.arg = 1, positive = TRUE))
      stop("bad input for argument 'n'") else n
  if (!is.numeric(doff) || length(doff) != 1 ||
      doff == 0)
    stop("bad input for argument 'doff'")
  L <- use.n
  if (length(mean) != L) mean <- rep_len(mean, L)
  if (length(sd)   != L) sd   <- rep_len(sd,   L)
  if (length(Lamb) != L) Lamb <- rep_len(Lamb, L)
  if (length(apar) != L) apar <- rep_len(apar, L)

  x1 <- rnorm(use.n, mean, sd)  # 1st marginal same
  txlamb <- if (doff > 0)  # Fed into qnorm()
       Lamb^(2/3) / (abs(doff) + Lamb^(2/3)) else
       log1p(Lamb) / (abs(doff) + log1p(Lamb))
  Delt <- pfun.N1p((x1 - mean) / sd, txlamb, apar)
  use.l <- if (doff > 0)  # use.lambda
     (abs(doff) * Delt / (1 - Delt))^1.5 else
     expm1(abs(doff) * Delt / (1 - Delt))
  x2 <- rpois(use.n, use.l)
  cbind(X1 = x1, X2 = x2)
}  # rN1pois









 eimjk.N1p <-
  function(y1,  # abscissae
           mean = 0, sdev = 1,
           mu2,  # aka Lamb
           apar = 0,  # aka alpha
           Integrand = TRUE,  # Omit dnorm()
           doff = 5,  # -log1p(10),  # May be <0
           jay = NULL,    # Used if which = "b"
           kay = NULL) {  # %in% 1:2 each

    zedd <- (y1 - mean) / sdev
    muu2 <- mu2  # lamb2
    dzedd.dmean <- (-1) / sdev
    dzedd.dsdev <- (-zedd) / sdev

    Dstr <- if (doff > 0)
       muu2^(2/3) / (muu2^(2/3) + abs(doff)) else
       log1p(muu2) / (abs(doff) + log1p(muu2))
    dDstr.dmuu2 <- if (doff > 0)
        (2 / 3) * abs(doff) / (muu2^(1/3) *
        (muu2^(2/3) + abs(doff))^2) else
        abs(doff) / ((1 + muu2) *
        (abs(doff) + log1p(muu2))^2)

    Prob <- pfun.N1p(zedd,
                     a = Dstr,  # Was muu2,
                     apar = apar)  # Delta
    dProb.dmuu2 <- pfun.N1p(zedd,
                            a = Dstr,  # Was muu2,
                            Integrand = TRUE,
                            apar, 1, "a") *
                   dDstr.dmuu2  # New
    dProb.dzedd <- pfun.N1p(zedd,
                            a = Dstr,  # Was muu2,
                            Integrand = TRUE,
                            apar, 1, "b")
    dProb.dapar <- pfun.N1p(zedd,
                            a = Dstr,  # was muu2,
                            Integrand = TRUE,
                            apar, 1, "apar")
    dProb.dmean <- dProb.dzedd * dzedd.dmean
    dProb.dsdev <- dProb.dzedd * dzedd.dsdev


    dmean.deta <- dsdev.deta <- 1  # Artificial here
    dmuu2.deta <- dapar.deta <- 1  # Artificial here

    ned2l.dmeanmuu2 <- dProb.dmean * dProb.dmuu2
    ned2l.dmeanapar <- dProb.dmean * dProb.dapar
    ned2l.dsdevmuu2 <- dProb.dsdev * dProb.dmuu2
    ned2l.dsdevapar <- dProb.dsdev * dProb.dapar
    ned2l.dmuu22 <- dProb.dmuu2^2  # orig. too
    ned2l.dapar2 <- dProb.dapar^2
    ned2l.dmuu2apar <- dProb.dmuu2 * dProb.dapar
    M <- 4
    n <- length(y1)  # Artificial
    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M)] <- (dProb.dmean * dmean.deta)^2
    wz[, iam(2, 2, M)] <- (dProb.dsdev * dsdev.deta)^2
    wz[, iam(3, 3, M)] <- ned2l.dmuu22 * dmuu2.deta^2
    wz[, iam(4, 4, M)] <- ned2l.dapar2 * dapar.deta^2
    wz[, iam(1, 2, M)] <- dProb.dmean * dmean.deta *
                          dProb.dsdev * dsdev.deta
    wz[, iam(1, 3, M)] <- ned2l.dmeanmuu2 *
                          dmean.deta * dmuu2.deta
    wz[, iam(1, 4, M)] <- ned2l.dmeanapar *
                          dmean.deta * dapar.deta
    wz[, iam(2, 3, M)] <- ned2l.dsdevmuu2 *
                          dsdev.deta * dmuu2.deta
    wz[, iam(2, 4, M)] <- ned2l.dsdevapar *
                          dsdev.deta * dapar.deta
    wz[, iam(3, 4, M)] <- ned2l.dmuu2apar *
                          dmuu2.deta * dapar.deta


  Delta <- Prob
  const3 <- if (doff > 0) {
    (if (Integrand) 1 else
     dnorm(y1, mean, sdev)) *
  (9/4) * abs(doff)^1.5 / (1 - Delta)^4
  } else {
    tmpodds <- abs(doff) * Delta / (1 - Delta)
    (if (Integrand) 1 else
      dnorm(y1, mean, sdev)) * doff^2 * 
      exp(2 * tmpodds) / (expm1(tmpodds) * 
      (1 - Delta)^4)
  }


  const3 * (if (length(jay) && length(kay))
  wz[, iam(jay, kay, M)] else wz)
}  # eimjk.N1p





































































# These functions are
# Copyright (C) 1998-2022 T.W. Yee, University of Auckland.
# All rights reserved.



















 multinomial <-
  function(zero = NULL, parallel = FALSE,
           nointercept = NULL, refLevel = "(Last)",
           imethod = 1,  # 20211008; DAMLMs need better init vals
           imu = NULL, byrow.arg = FALSE,  # 20211105
           whitespace = FALSE) {


      

  d.mlm <- NULL  # 20211105


  if (length(refLevel) != 1)
    stop("the length of 'refLevel' must be one")

  if ( is.numeric(refLevel) &&
      !is.Numeric(refLevel, integer.valued = TRUE, positive = TRUE))
    stop("argument 'refLevel' is not a positive integer")

  if (is.character(refLevel)) {
    if (refLevel == "(Last)")
      refLevel <- -1
  }
  if (is.factor(refLevel)) {
    if (is.ordered(refLevel))
      warning("argument 'refLevel' is from an ordered factor")

    refLevel <- as.character(refLevel) == levels(refLevel)
    refLevel <- (seq_along(refLevel))[refLevel]
    if (!is.Numeric(refLevel, length.arg = 1,
                    integer.valued = TRUE, positive = TRUE))
      stop("could not coerce 'refLevel' into a single ",
           "positive integer")
  }


  if ((LLL <- length(d.mlm)) > 0) {
    sd.mlm <- unique(sort(d.mlm))  # Values must be unique
    if (!is.Numeric(d.mlm, integer.valued = FALSE,
                    length.arg = length(sd.mlm), positive = FALSE))
      stop("bad input for 'd.mlm'")
    if (is.numeric(refLevel) && any(refLevel == sd.mlm))
      stop("cannot deflate the reference level")
    d.mlm <- as.vector(sd.mlm)  # Replacement
  }





  stopifnot(is.logical(whitespace) &&
            length(whitespace) == 1)
  fillerChar <- ifelse(whitespace, " ", "")


  new("vglmff",
      blurb = c(ifelse(length(d.mlm) > 0,
                       "Deflated-altered m", "M"),
            "ultinomial logit model\n\n",
            ifelse(length(d.mlm) > 0, "(Altered) ", ""),
            "Links:    ",
         if (is.numeric(refLevel)) {
         if (refLevel < 0) {
           ifelse(whitespace,
                  "log(mu[,j] / mu[,M+1]), j = 1:M,\n",
                  "log(mu[,j]/mu[,M+1]), j=1:M,\n")
         } else {
            if (refLevel == 1) {
              paste("log(mu[,", "j]", fillerChar, "/", fillerChar,
                    "mu[,", refLevel, "]), j",
                    fillerChar, "=", fillerChar, "2:(M+1),\n",
                    sep = "")
            } else {
              paste("log(mu[,", "j]", fillerChar, "/",
                   "mu[,", refLevel, "]), j",
                    fillerChar, "=", fillerChar, "c(1:", refLevel-1,
                    ",", fillerChar, refLevel+1, ":(M+1)),\n",
                    sep = "")
            }
         }
         } else {  # refLevel is character
           paste("log(mu[,", "j]", fillerChar, "/",
                "mu[,'", refLevel, "']), j",
                 fillerChar, " != '", fillerChar, refLevel,
                 "',\n",
                 sep = "")
         },
         "Variance: ",
           ifelse(whitespace,
                  "mu[,j] * (1 - mu[,j]); -mu[,j] * mu[,k]",
                  "mu[,j]*(1-mu[,j]); -mu[,j]*mu[,k]")),

  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           apply.int = TRUE,
                           constraints = constraints)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero ,
                     M = M, M1 = M,
                     predictors.names = predictors.names)
    constraints <- cm.nointercept.VGAM(constraints, x,
                                       .nointercept , M)
  }), list( .parallel = parallel, .zero = zero,
            .d.mlm = d.mlm, .nointercept = nointercept ))),

  deviance = Deviance.categorical.data.vgam,

  infos = eval(substitute(function(...) {
    list(parallel = .parallel ,
         refLevel = .refLevel ,  # original
         M1 = -1,
         link = "multilogitlink",
         link1parameter = FALSE,  # The link is multiparameter, but
         mixture.links = FALSE,  # the link is not a mixture.
         expected = TRUE,
         d.mlm = .d.mlm ,
         imethod = .imethod ,
         multipleResponses = FALSE,
         parameters.names = as.character(NA),
         zero = .zero )
  }, list( .zero = zero,
           .refLevel = refLevel, .imethod = imethod,
           .d.mlm = d.mlm,
           .parallel = parallel ))),
  initialize = eval(substitute(expression({

    if (is.factor(y) && is.ordered(y))
      warning("response should be nominal, not ordinal")


    delete.zero.colns <- TRUE
    eval(process.categorical.data.VGAM)

    if ( .imethod == 2) {  # MLE for intercept-only.
      mustart <- matrix(colMeans(y), n, NCOL(y), byrow = TRUE)
    }
    if ( .imethod > 2) stop("argument 'imethod' unmatched")
    

    M <- ncol(y) - 1
    use.refLevel <- if (is.numeric( .refLevel )) {
      if ( .refLevel < 0) M + 1 else ( .refLevel )
    } else {  # Is character. Match it with the response levels.
      tmp6 <- match( .refLevel , colnames(y))
      if (is.na(tmp6))
        stop("could not match argument 'refLevel' with any ",
             "columns of the response matrix")
      tmp6
    }


    if (use.refLevel > (M + 1))
      stop("argument 'refLevel' has a value that is too high")
    extra$use.refLevel <- use.refLevel  # Used in all other slots.
    sd.mlm <- extra$d.mlm <- as.vector( .d.mlm )  # May be NULL, okay.
    pmone <- rep_len(1, M + 1)
    if (length(sd.mlm)) {
      pmone[sd.mlm] <- -1
    }
    signvec <- extra$signvec.damlm <- pmone


    allbut.refLevel <- seq(M + 1)[-use.refLevel]
    predictors.names <-
      paste("log(mu[,", allbut.refLevel,
            "]", .fillerChar, "/", .fillerChar,
            "mu[,", use.refLevel, "])", sep = "")
    if (length(sd.mlm) > 0 && any(signvec[-use.refLevel] < 0)){
      ind.d <- sd.mlm
      ind.d[use.refLevel < ind.d] <-
      ind.d[use.refLevel < ind.d] - 1  # Scrunch RHS up
      predictors.names[ind.d] <- paste0("-", predictors.names[ind.d])
    }

    extra$colnames.y <- colnames(y)


    imu <- as.vector( .imu )
    if (length(imu)) {
      mustart <- matrix(imu, n, NCOL(y), byrow = .byrow.arg )
    }
  }), list( .refLevel = refLevel, .fillerChar = fillerChar,
            .d.mlm = d.mlm, .imu = imu, .byrow.arg = byrow.arg,
            .imethod = imethod, .whitespace = whitespace ))),
  linkinv = eval(substitute( function(eta, extra = NULL) {
    if (anyNA(eta))
      warning("there are NAs in eta in slot inverse")
    ans <-
      multilogitlink(eta,  # .refLevel , d.mlm = .d.mlm ,
                     refLevel = extra$use.refLevel, inverse = TRUE)
    if (anyNA(ans))
      warning("there are NAs here in slot linkinv")
    if (min(ans) == 0 || max(ans) == 1)
      warning("fitted probabilities numerically 0 or 1 occurred")

    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = 1)
  }), list( .refLevel = refLevel, .d.mlm = d.mlm )),
  last = eval(substitute(expression({
    misc$link <- "multilogitlink"

    misc$earg <- list(multilogitlink = list(
      M = M,
      refLevel = use.refLevel
    ))

    dy <- dimnames(y)
    if (!is.null(dy[[2]]))
      dimnames(fit$fitted.values) <- dy

    misc$nointercept <- ( .nointercept )
    misc$parallel <- ( .parallel )
    misc$refLevel <- use.refLevel
    misc$refLevel.orig <- ( .refLevel )
    misc$zero <- ( .zero )
  }), list( .refLevel = refLevel, .nointercept = nointercept,
            .parallel = parallel,
            .d.mlm = d.mlm, .zero = zero ))),

  linkfun = eval(substitute( function(mu, extra = NULL) {

    multilogitlink(mu,   # d.mlm = .d.mlm ,
                   refLevel = extra$use.refLevel)
  }), list( .refLevel = refLevel, .d.mlm = d.mlm )),

  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ycounts <- if (is.numeric(extra$orig.w))
                   y * w / extra$orig.w else
                   y * w  # Convert proportions to counts
      nvec <- if (is.numeric(extra$orig.w))
                round(w / extra$orig.w) else round(w)

      smallno <- 1.0e4 * .Machine$double.eps
      if (max(abs(ycounts - round(ycounts))) > smallno)
        warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)
      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dmultinomial(x = ycounts, size = nvec, prob = mu,
                      log = TRUE, dochecking = FALSE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("multinomial", "VGAMcategorical"),




  hadof = eval(substitute(
  function(eta, extra = list(),
           linpred.index = 1,
           w = 1, dim.wz = c(NROW(eta), NCOL(eta) * (NCOL(eta)+1)/2), 
           deriv = 1, ...) {


  multinomial.eim.deriv1 <- function(mu, jay, w) {
    M <- ncol(mu)  # Not ncol(mu) - 1 coz 1 coln has been deleted
    MM12 <- M * (M + 1) / 2  # Full matrix
    wz1 <- matrix(0, NROW(mu), MM12)
    ind5 <- iam(NA, NA, M = M, both = TRUE)
    for (i in 1:MM12) {
      i1 <- ind5$row.index[i]
      j1 <- ind5$col.index[i]
      if (i1 == jay && j1 == jay) {
        wz1[, iam(i1, j1, M = M)] <-  (1 - 2 * mu[, i1])  # *
      }
      if (i1 == jay && j1 != jay) {
        wz1[, iam(i1, j1, M = M)] <- -mu[, j1]  # -(1 - 2 * mu[, i1]) *
      }
      if (i1 != jay && j1 == jay) {
        wz1[, iam(i1, j1, M = M)] <- -mu[, i1]  # -(1 - 2 * mu[, j1]) *
      }
    }  # for (i)
    c(w) * wz1
  }  # multinomial.eim.deriv1



  multinomial.eim.deriv2 <- function(mu, jay, w) {
    M <- ncol(mu)  # Not ncol(mu) - 1 since 1 coln has been deleted
    zz.yettodo <- 0  # NA_real_
    MM12 <- M * (M + 1) / 2  # Full matrix
    wz1 <- matrix(NA_real_, NROW(mu), MM12)
    ind5 <- iam(NA, NA, M = M, both = TRUE)
    for (i in 1:MM12) {
      i1 <- ind5$row.index[i]
      j1 <- ind5$col.index[i]
      if (i1 == jay && j1 == jay) {
        wz1[, iam(i1, j1, M = M)] <-  zz.yettodo
      }
      if (i1 != jay && j1 == i1) {
        wz1[, iam(i1, j1, M = M)] <- -mu[, i1] * mu[, jay] *
                                 (1 - 2 * mu[, i1] - 2 * mu[, jay] +
                                      6 * mu[, i1] * mu[, jay])
                                            
      }
      if (i1 == jay && j1 != jay) {
        wz1[, iam(i1, j1, M = M)] <- zz.yettodo
      }
      if (i1 != jay && j1 == jay) {
        wz1[, iam(i1, j1, M = M)] <- zz.yettodo
      }
      if (any(is.na(wz1[, iam(i1, j1, M = M)]))) {
        wz1[, iam(i1, j1, M = M)] <-       2 * mu[, i1]  * 
                                               mu[, j1] * mu[, jay] *
                                      (1 - 3 * mu[, jay])
      }
    }  # for (i)
 cat("\n\n\n\n")
    c(w) * wz1
  }  # multinomial.eim.deriv2


    M <- NCOL(eta)
    use.refLevel <- extra$use.refLevel  # Restore its value
    if (!is.numeric(use.refLevel)) {
      warning("variable 'use.refLevel' cannot be found. ",
              "Trying the original value.")
      use.refLevel <- .refLevel  # Only if numeric...
      if (use.refLevel == "(Last)")
        use.refLevel <- M+1
    }
    mu.use <- multilogitlink(eta, refLevel = use.refLevel,
                             inverse = TRUE)
    mu.use <- pmax(mu.use, .Machine$double.eps * 1.0e-0)


    index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
    myinc <- (index$row.index >= use.refLevel)
    index$row.index[myinc] <- index$row.index[myinc] + 1
    myinc <- (index$col.index >= use.refLevel)
    index$col.index[myinc] <- index$col.index[myinc] + 1


    switch(as.character(deriv),
       "0" =  {
         wz <- -mu.use[, index$row] * mu.use[, index$col]
         wz[, 1:M] <- wz[, 1:M] + mu.use[, -use.refLevel ]
         c(w) * wz
       },
       "1" = {
         multinomial.eim.deriv1(mu.use[, -use.refLevel, drop = FALSE],
                                jay = linpred.index, w = w)
       },
       "2" = {
         multinomial.eim.deriv2(mu.use[, -use.refLevel, drop = FALSE],
                                jay = linpred.index, w = w)
       },
       stop("argument 'deriv' must be 0 or 1 or 2"))
  }, list( .refLevel = refLevel,  # End of @hadof
           .d.mlm = d.mlm ))),




  validparams = eval(substitute(function(eta, y, extra = NULL) {
    probs <-
      multilogitlink(eta, refLevel = extra$use.refLevel,
                     # d.mlm = extra$d.mlm,
                     inverse = TRUE)  # ( .refLevel )
    okay1 <- all(is.finite(probs)) && all(0 < probs)
    if (!length(extra$d.mlm))
      okay1 <- okay1 && all(probs < 1)
    okay1 <- all(is.finite(probs)) && all(0 < probs & probs < 1)
    okay1
  }, list( .refLevel = refLevel, .d.mlm = d.mlm ))),
  deriv = eval(substitute(expression({
    signvec <- extra$signvec.damlm
    use.refLevel <- extra$use.refLevel  # Restore its value
    d.mlm <- extra$d.mlm  # Deflated excluding baseline gp
    ansd <- ( y[, -use.refLevel, drop = FALSE] -
             mu[, -use.refLevel, drop = FALSE])  # AMLM


    if (length(d.mlm) > 0) {  # DAMLM3
      yy.use <-  y[, -use.refLevel, drop = FALSE]
      mu.use <- mu[, -use.refLevel, drop = FALSE]
      ind.d4 <- which(signvec[-use.refLevel] < 0)
      dl.dprob.damlm <-
        yy.use[, ind.d4, drop = FALSE] / (
        mu.use[, ind.d4, drop = FALSE]) -
        c(y[, use.refLevel] / mu[, use.refLevel])  # Recycling colns
      ansd[, ind.d4] <- ansd[, ind.d4] -  # DAMLM3
        2 * mu.use[, ind.d4, drop = FALSE] * dl.dprob.damlm
    }  # length(d.mlm) > 0  # DAMLM3
    c(w) * ansd
  }), list( .refLevel = refLevel, .d.mlm = d.mlm ))),
  weight = eval(substitute(expression({
    mytiny <- (mu <       sqrt(.Machine$double.eps)) |
              (mu > 1.0 - sqrt(.Machine$double.eps))

    if (M == 1) {
      wz <- mu[, 3 - use.refLevel] * (1 - mu[, 3 - use.refLevel])
    } else {  # M > 1
      index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
      myinc <- (index$row.index >= use.refLevel)
      index$row.index[myinc] <- index$row.index[myinc] + 1
      myinc <- (index$col.index >= use.refLevel)
      index$col.index[myinc] <- index$col.index[myinc] + 1
      wz <- -mu[, index$row] * mu[, index$col]
      wz[, 1:M] <- wz[, 1:M] + mu[, -use.refLevel ]  # AMLM
    }


    if (length(d.mlm) > 0) {  # DAMLM3
      for (rowval in ind.d4) {
        ind3 <- iam(rowval, 1:M, M = M)  # Both rows AND colns
        wz[, ind3] <- wz[, ind3] - c(
          2 * mu.use[, rowval] / mu[, use.refLevel])
      }  # rowval
      for (colval in ind.d4) {
        ind3 <- iam(1:M, colval, M = M)
        wz[, ind3] <- wz[, ind3] - c(
          2 * mu.use[, colval] / mu[, use.refLevel])
      }  # colval
      for (diagval in ind.d4) {
        ind3 <- iam(diagval, diagval, M = M)  # scalar
        wz[, ind3] <- wz[, ind3] +
          4 * (mu.use[, diagval, drop = FALSE]^2) * c(
          1 / mu.use[, diagval] + 0 / mu[, use.refLevel]) + c(
          2 * mu.use[, diagval] / mu[, use.refLevel])  # rm once
      }  # diagval
 





    }  # DAMLM3





    atiny <- (mytiny %*% rep(1, ncol(mu))) > 0
    if (any(atiny)) {
      if (M == 1)
        wz[atiny] <- wz[atiny] * (1 + .Machine$double.eps^0.5) +
                                      .Machine$double.eps else
        wz[atiny, 1:M] <- .Machine$double.eps +
        wz[atiny, 1:M] * (1 + .Machine$double.eps^0.5)
    }  # atiny
    c(w) * wz
  }), list( .refLevel = refLevel, .d.mlm = d.mlm ))))
}  # multinomial(), was damultinomial()







 process.categorical.data.VGAM <- expression({






    extra$y.integer <- TRUE


  if (!all(w == 1))
    extra$orig.w <- w

  if (!is.matrix(y)) {
    yf <- as.factor(y)
    lev <- levels(yf)
    llev <- length(lev)
    nn <- length(yf)
    y <- matrix(0, nn, llev)
    y[cbind(1:nn,as.vector(unclass(yf)))] <- 1
    dimnames(y) <- list(names(yf), lev)

    if (llev <= 1)
      stop("the response matrix does not have 2 or more columns")
  } else {
    nn <- nrow(y)
  }


  nvec <- rowSums(y)

  if (min(y) < 0 || any(round(y) != y))
    stop("the response must be non-negative counts (integers)")

  if (!exists("delete.zero.colns") ||
      (exists("delete.zero.colns") && delete.zero.colns)) {
    sumy2 <- colSums(y)
    if (any(index <- sumy2 == 0)) {
      y <- y[, !index, drop = FALSE]
      sumy2 <- sumy2[!index]
      if (all(index) || ncol(y) <= 1)
        stop("'y' matrix has 0 or 1 columns")
      warning("Deleted ", sum(!index),
              " columns of the response matrix due to zero counts")
    }
  }


  if (any(miss <- (nvec == 0))) {
    smiss <- sum(miss)
    warning("Deleted ", smiss,
            " rows of the response matrix due to zero counts")
    x <- x[!miss,, drop = FALSE]
    y <- y[!miss,, drop = FALSE]
    w <- cbind(w)
    w <- w[!miss,, drop = FALSE]

    nvec <- nvec[!miss]
    nn <- nn - smiss
  }

  w <- w * nvec

  nvec[nvec == 0] <- 1
  y <- prop.table(y, 1)   # Convert to proportions


  if (length(mustart) + length(etastart) == 0) {
      mustart <- y + (1 / ncol(y) - y) / nvec
  }
})





Deviance.categorical.data.vgam <-
  function(mu, y, w, residuals = FALSE, eta, extra = NULL,
           summation = TRUE) {




  if (ncol(y) == 1 || ncol(mu) == 1)
    stop("arguments 'y' and 'mu' must have at least 2 columns")



  double.eps <- sqrt( .Machine$double.xmin )
  devy <- y
  nonz <- (y != 0)
  devy[nonz] <- y[nonz] * log(y[nonz])

  devmu <- 0 * y  # filler; y*log(mu) gives a warning (fixed up anyway).
  if (any(smallmu <- (mu * (1 - mu) < double.eps))) {
    warning("fitted values close to 0 or 1")
    smu <- mu[smallmu]
    smy <-  y[smallmu]
    smu <- ifelse(smu < double.eps, double.eps, smu)
    devmu[smallmu] <- smy * log(smu)
  }
  devmu[!smallmu] <- y[!smallmu] * log(mu[!smallmu])

  devi <- 2 * (devy - devmu)

  if (residuals) {
    M <- if (is.matrix(eta)) ncol(eta) else 1
    if (M > 1)
      return(NULL)
    devi <- devi %*% rep_len(1, ncol(devi))  # Dev = \sum_i devi[i]
    return(c(sign(y[, 1] - mu[, 1]) * sqrt(abs(devi) * w)))
  } else {
    dev.elts <- c(w) * devi
    if (summation) {
      sum(dev.elts)
    } else {
      dev.elts
    }
  }
}







 dmultinomial <-
    function(x, size = NULL, prob, log = FALSE,
             dochecking = TRUE, smallno = 1.0e-7) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  x <- as.matrix(x)
  prob <- as.matrix(prob)
  if (((K <- ncol(x)) <= 1) ||
       ncol(prob) != K)
    stop("arguments 'x' and 'prob' must be matrices with ",
         "two or more columns")
  if (dochecking) {
    if (min(prob) < 0)
      stop("argument 'prob' contains some negative values")
    if (any(abs((rsprob <- rowSums(prob)) - 1) > smallno))
      stop("some rows of 'prob' do not add to unity")
    if (any(abs(x - round(x)) > smallno))
      stop("argument 'x' should be integer-valued")
    if (length(size)) {
      if (any(abs(size - rowSums(x)) > smallno))
        stop("rowSums(x) does not agree with argument 'size'")
    } else {
      size <- round(rowSums(x))
    }
  } else {
    if (!length(size))
      size <- round(rowSums(prob))
  }
  logdensity <- lfactorial(size) + rowSums(x*log(prob) - lfactorial(x))
  if (log.arg) logdensity else exp(logdensity)
}  # dmultinomial()










 sratio <- function(link = "logitlink",
                    parallel = FALSE, reverse = FALSE, zero = NULL,
                    whitespace = FALSE) {
  link <- as.list(substitute(link))
  earg  <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.logical(reverse) || length(reverse) != 1)
    stop("argument 'reverse' must be a single logical")

  stopifnot(is.logical(whitespace) &&
            length(whitespace) == 1)
  fillerChar <- ifelse(whitespace, " ", "")


  new("vglmff",
  blurb = c("Stopping ratio model\n\n",
            "Links:    ",
            namesof(if (reverse)
            ifelse(whitespace, "P[Y = j+1|Y <= j+1]",
                               "P[Y=j+1|Y<=j+1]") else
            ifelse(whitespace, "P[Y = j|Y >= j]",
                               "P[Y=j|Y>=j]"),
                   link, earg = earg), "\n",
            "Variance: ",
            ifelse(whitespace,
                   "mu[,j] * (1 - mu[,j]); -mu[,j] * mu[,k]",
                   "mu[,j]*(1-mu[,j]); -mu[,j]*mu[,k]")),
  infos = eval(substitute(function(...) {
    list(M1 = NA,  # zz -1?
         Q1 = NA,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = as.character(NA),
         parallel = .parallel ,
         reverse = .reverse ,
         whitespace = .whitespace ,
         zero = .zero ,
         link = .link )
  }, list( .link = link,
           .zero = zero,
           .parallel = parallel,
           .reverse = reverse,
           .whitespace = whitespace ))),

  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints = constraints)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero ,
                     M = M, M1 = M,
                     predictors.names = predictors.names)
  }), list( .parallel = parallel, .zero = zero ))),
  deviance = Deviance.categorical.data.vgam,

  initialize = eval(substitute(expression({

    if (is.factor(y) && !is.ordered(y))
      warning("response should be ordinal---see ordered()")



    delete.zero.colns <- TRUE
    eval(process.categorical.data.VGAM)
    extra$wy.prod <- TRUE
    M <- ncol(y) - 1

    mynames <- if ( .reverse )
      paste("P[Y", .fillerChar, "=", .fillerChar, 2:(M+1), "|Y",
             .fillerChar, "<=", .fillerChar, 2:(M+1), "]", sep = "") else
      paste("P[Y", .fillerChar, "=", .fillerChar, 1:M,     "|Y",
             .fillerChar, ">=", .fillerChar, 1:M,     "]", sep = "")
    predictors.names <-
      namesof(mynames, .link , short = TRUE, earg = .earg )
    y.names <- param.names("mu", M+1)

    extra$mymat <- if ( .reverse ) tapplymat1(y, "cumsum") else
                  tapplymat1(y[, ncol(y):1], "cumsum")[, ncol(y):1]

    extra$colnames.y  <- colnames(y)
  }), list( .earg = earg, .link = link, .reverse = reverse,
            .fillerChar = fillerChar,
            .whitespace = whitespace ))),

  linkinv = eval(substitute( function(eta, extra = NULL) {
    if (!is.matrix(eta))
      eta <- as.matrix(eta)
    fv.mat <- if ( .reverse ) {
      M <- NCOL(eta)
      djr <- eta2theta(eta, .link , earg = .earg )
      temp <- tapplymat1(1 - djr[, M:1], "cumprod")[, M:1]
      cbind(1, djr) * cbind(temp, 1)
    } else {
      dj <- eta2theta(eta, .link , earg = .earg )
      temp <- tapplymat1(1 - dj, "cumprod")
      cbind(dj, 1) * cbind(1, temp)
    }
    label.cols.y(fv.mat, colnames.y = extra$colnames.y, NOS = 1)
  }, list( .earg = earg, .link = link, .reverse = reverse) )),
  last = eval(substitute(expression({
    misc$link <- rep_len( .link , M)
    names(misc$link) <- mynames

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$parameters <- mynames
    misc$reverse <- .reverse
    misc$fillerChar <- .fillerChar
    misc$whitespace <- .whitespace

    extra <- list()  # kill what was used
  }), list( .earg = earg, .link = link, .reverse = reverse,
            .fillerChar = fillerChar,
            .whitespace = whitespace ))),
  linkfun = eval(substitute( function(mu, extra = NULL) {
    cump <- tapplymat1(mu, "cumsum")
    if ( .reverse ) {
      djr <- mu[, -1] / cump[, -1]
      theta2eta(djr, .link , earg = .earg )
    } else {
      M <- ncol(mu) - 1
      dj <- if (M == 1) mu[, 1] else
           mu[, 1:M] / (1 - cbind(0, cump[, 1:(M-1)]))
      theta2eta(dj, .link , earg = .earg )
    }
  }, list( .earg = earg, .link = link, .reverse = reverse) )),
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ycounts <- if (is.numeric(extra$orig.w))
                 y * w / extra$orig.w else
                 y * w  # Convert proportions to counts
      nvec <- if (is.numeric(extra$orig.w))
              round(w / extra$orig.w) else round(w)

      smallno <- 1.0e4 * .Machine$double.eps
      if (max(abs(ycounts - round(ycounts))) > smallno)
        warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dmultinomial(x = ycounts, size = nvec, prob = mu,
                      log = TRUE, dochecking = FALSE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    },
  vfamily = c("sratio", "VGAMordinal", "VGAMcategorical"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    djr <- eta2theta(eta, .link , earg = .earg )  # dj or djr
    okay1 <- all(is.finite(djr)) && all(0 < djr & djr < 1)
    okay1
  }, list( .earg = earg, .link = link, .reverse = reverse) )),
  deriv = eval(substitute(expression({
    if (!length(extra$mymat)) {
      extra$mymat <- if ( .reverse ) tapplymat1(y, "cumsum") else
                     tapplymat1(y[, ncol(y):1], "cumsum")[, ncol(y):1]
    }
    if ( .reverse ) {
      djr <- eta2theta(eta, .link , earg = .earg )
      Mp1 <- ncol(extra$mymat)
      c(w) * (y[, -1] / djr - extra$mymat[, -Mp1] / (1 - djr)) *
        dtheta.deta(djr, .link , earg = .earg )
    } else {
      dj <- eta2theta(eta, .link , earg = .earg )
      c(w) * (y[, -ncol(y)] / dj - extra$mymat[, -1] / (1 - dj)) *
        dtheta.deta(dj, .link , earg = .earg )
    }
  }), list( .earg = earg, .link = link, .reverse = reverse) )),
  weight = eval(substitute(expression({
    if ( .reverse ) {
      cump <- tapplymat1(mu, "cumsum")
      ddjr.deta <- dtheta.deta(djr, .link , earg = .earg )
      wz <- c(w) * ddjr.deta^2 *
           (mu[, -1] / djr^2 + cump[, 1:M] / (1 - djr)^2)
    } else {
      ccump <- tapplymat1(mu[, ncol(mu):1], "cumsum")[, ncol(mu):1]
      ddj.deta <- dtheta.deta(dj, .link , earg = .earg )
      wz <- c(w) * ddj.deta^2 *
           (mu[, 1:M] / dj^2 + ccump[, -1] / (1 - dj)^2)
    }

    wz
  }), list( .earg = earg, .link = link, .reverse = reverse ))))
}  # sratio()






 cratio <- function(link = "logitlink",
                    parallel = FALSE, reverse = FALSE, zero = NULL,
                    whitespace = FALSE) {
  link <- as.list(substitute(link))
  earg  <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.logical(reverse) || length(reverse) != 1)
    stop("argument 'reverse' must be a single logical")

  stopifnot(is.logical(whitespace) &&
            length(whitespace) == 1)
  fillerChar <- ifelse(whitespace, " ", "")


  new("vglmff",
  blurb = c("Continuation ratio model\n\n",
            "Links:    ",
            namesof(if (reverse)
            ifelse(whitespace, "P[Y < j+1|Y <= j+1]",
                               "P[Y<j+1|Y<=j+1]") else
            ifelse(whitespace, "P[Y > j|Y >= j]",
                               "P[Y>j|Y>=j]"),
                   link, earg = earg),
            "\n",
            "Variance: ",
            ifelse(whitespace,
                   "mu[,j] * (1 - mu[,j]); -mu[,j] * mu[,k]",
                   "mu[,j]*(1-mu[,j]); -mu[,j]*mu[,k]")),

  infos = eval(substitute(function(...) {
    list(M1 = NA,  # zz -1?
         Q1 = NA,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = as.character(NA),
         parallel = .parallel ,
         reverse = .reverse ,
         whitespace = .whitespace ,
         zero = .zero ,
         link = .link )
  }, list( .link = link,
           .zero = zero,
           .parallel = parallel,
           .reverse = reverse,
           .whitespace = whitespace ))),


  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints = constraints)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero ,
                     M = M, M1 = M,
                     predictors.names = predictors.names)
  }), list( .parallel = parallel, .zero = zero ))),


  deviance = Deviance.categorical.data.vgam,

  initialize = eval(substitute(expression({

    if (is.factor(y) && !is.ordered(y))
      warning("response should be ordinal---see ordered()")



    delete.zero.colns <- TRUE
    eval(process.categorical.data.VGAM)
    M <- ncol(y) - 1

    mynames <- if ( .reverse )
      paste("P[Y", .fillerChar, "<", .fillerChar, 2:(M+1), "|Y",
            .fillerChar, "<=", .fillerChar, 2:(M+1), "]", sep = "") else
      paste("P[Y", .fillerChar, ">", .fillerChar, 1:M,     "|Y",
            .fillerChar, ">=", .fillerChar, 1:M,     "]", sep = "")
    predictors.names <-
      namesof(mynames, .link , earg = .earg , short = TRUE)
    y.names <-  param.names("mu", M+1)

    extra$mymat <- if ( .reverse )
                   tapplymat1(y, "cumsum") else
                   tapplymat1(y[, ncol(y):1], "cumsum")[, ncol(y):1]

    extra$colnames.y  <- colnames(y)
  }), list( .earg = earg, .link = link, .reverse = reverse,
            .fillerChar = fillerChar,
            .whitespace = whitespace ))),

  linkinv = eval(substitute( function(eta, extra = NULL) {
    if (!is.matrix(eta))
      eta <- as.matrix(eta)
    fv.mat <- if ( .reverse ) {
      M <- ncol(eta)
      djrs <- eta2theta(eta, .link , earg = .earg )
      temp <- tapplymat1(djrs[, M:1], "cumprod")[, M:1]
      cbind(1, 1 - djrs) * cbind(temp, 1)
    } else {
      djs <- eta2theta(eta, .link , earg = .earg )
      temp <- tapplymat1(djs, "cumprod")
      cbind(1 - djs, 1) * cbind(1, temp)
    }
    label.cols.y(fv.mat, colnames.y = extra$colnames.y, NOS = 1)
  }, list( .earg = earg, .link = link, .reverse = reverse) )),
  last = eval(substitute(expression({

    misc$link <- rep_len( .link , M)
    names(misc$link) <- mynames

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$parameters <- mynames
    misc$reverse <- .reverse
    misc$fillerChar <- .fillerChar
    misc$whitespace <- .whitespace


    extra <- list()  # kill what was used
  }), list( .earg = earg, .link = link, .reverse = reverse,
            .fillerChar = fillerChar,
            .whitespace = whitespace ))),
  linkfun = eval(substitute( function(mu, extra = NULL) {
    cump <- tapplymat1(mu, "cumsum")
    if ( .reverse ) {
      djrs <- 1 - mu[, -1] / cump[, -1]
      theta2eta(djrs, .link , earg = .earg )
    } else {
      M <- ncol(mu) - 1
      djs <- if (M == 1) 1 - mu[, 1] else
             1 - mu[, 1:M] / (1 - cbind(0, cump[, 1:(M-1)]))
      theta2eta(djs, .link , earg = .earg )
    }
  }, list( .earg = earg, .link = link, .reverse = reverse) )),
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ycounts <- if (is.numeric(extra$orig.w))
                 y * w / extra$orig.w else
                 y * w  # Convert proportions to counts
      nvec <- if (is.numeric(extra$orig.w))
              round(w / extra$orig.w) else round(w)

      smallno <- 1.0e4 * .Machine$double.eps
      if (max(abs(ycounts - round(ycounts))) > smallno)
        warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dmultinomial(x = ycounts, size = nvec, prob = mu,
                      log = TRUE, dochecking = FALSE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("cratio", "VGAMordinal", "VGAMcategorical"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    djrs <- eta2theta(eta, .link , earg = .earg )  # djs or djrs
    okay1 <- all(is.finite(djrs)) && all(0 < djrs & djrs < 1)
    okay1
  }, list( .earg = earg, .link = link, .reverse = reverse) )),

  deriv = eval(substitute(expression({
    if (!length(extra$mymat)) {
      extra$mymat <- if ( .reverse ) tapplymat1(y, "cumsum") else
                     tapplymat1(y[, ncol(y):1], "cumsum")[, ncol(y):1]
    }
    if ( .reverse ) {
      djrs <- eta2theta(eta, .link , earg = .earg )
      Mp1 <- ncol(extra$mymat)
      -c(w) * (y[, -1]/(1 - djrs) - extra$mymat[, -Mp1]/djrs) *
        dtheta.deta(djrs, .link , earg = .earg )
    } else {
      djs <- eta2theta(eta, .link , earg = .earg )
      -c(w) * (y[, -ncol(y)]/(1 - djs) - extra$mymat[, -1]/djs) *
        dtheta.deta(djs, .link , earg = .earg )
    }
  }), list( .earg = earg, .link = link, .reverse = reverse) )),

  weight = eval(substitute(expression({
    if ( .reverse ) {
      cump <- tapplymat1(mu, "cumsum")
      ddjrs.deta <- dtheta.deta(djrs, .link , earg = .earg )
      wz <- c(w) * ddjrs.deta^2 *
            (mu[, -1] / (1 - djrs)^2 + cump[, 1:M] / djrs^2)
    } else {
      ccump <- tapplymat1(mu[, ncol(mu):1], "cumsum")[, ncol(mu):1]
      ddjs.deta <- dtheta.deta(djs, .link , earg = .earg )
      wz <- c(w) * ddjs.deta^2 *
            (mu[, 1:M] / (1 - djs)^2 + ccump[, -1] / djs^2)
    }
    wz
  }), list( .earg = earg, .link = link, .reverse = reverse ))))
}  # cratio()






 vglm.multinomial.deviance.control <-
  function(maxit = 31, panic = FALSE, ...) {
  if (maxit < 1) {
      warning("bad value of maxit; using 31 instead")
      maxit <- 31
  }
  list(maxit = maxit, panic = as.logical(panic)[1])
}


 vglm.multinomial.control <-
  function(maxit = 21, panic = FALSE,
           criterion = c("aic1", "aic2", names( .min.criterion.VGAM )),
           ...) {
  if (mode(criterion) != "character" && mode(criterion) != "name")
    criterion <- as.character(substitute(criterion))
  criterion <- match.arg(criterion,
      c("aic1", "aic2", names( .min.criterion.VGAM )))[1]

  if (maxit < 1) {
    warning("bad value of maxit; using 21 instead")
    maxit <- 21
  }
  list(maxit = maxit,
       panic = as.logical(panic)[1],
       criterion = criterion,
       min.criterion = c("aic1" = FALSE, "aic2" = TRUE,
                         .min.criterion.VGAM))
}


 vglm.VGAMcategorical.control <-
  function(maxit = 30,
           trace = FALSE,
           panic = TRUE, ...) {
  if (maxit < 1) {
    warning("bad value of maxit; using 200 instead")
    maxit <- 200
  }
  list(maxit = maxit,
       trace = as.logical(trace)[1],
       panic = as.logical(panic)[1])
}  # vglm.VGAMcategorical.control()








 cumulative <-
  function(link = "logitlink",
           parallel = FALSE,  # Does not apply to the intercept
           reverse = FALSE,
           multiple.responses = FALSE,
           whitespace = FALSE) {


  apply.parint <- FALSE


  link <- as.list(substitute(link))
  earg  <- link2list(link)
  link <- attr(earg, "function.name")



  stopifnot(is.logical(whitespace) &&
            length(whitespace) == 1)
  fillerChar <- ifelse(whitespace, " ", "")


  if (!is.logical(multiple.responses) || length(multiple.responses) != 1)
    stop("argument 'multiple.responses' must be a single logical")
  if (!is.logical(reverse) || length(reverse) != 1)
    stop("argument 'reverse' must be a single logical")


  new("vglmff",
  blurb = if ( multiple.responses )
          c(paste("Multivariate cumulative", link, "model\n\n"),
          "Links:   ",
          namesof(if (reverse)
                  ifelse(whitespace, "P[Y1 >= j+1]", "P[Y1>=j+1]") else
                  ifelse(whitespace, "P[Y1 <= j]",   "P[Y1<=j]"),
                  link, earg = earg),
          ", ...") else
          c(paste("Cumulative", link, "model\n\n"),
          "Links:   ",
          namesof(if (reverse)
                  ifelse(whitespace, "P[Y >= j+1]", "P[Y>=j+1]") else
                  ifelse(whitespace, "P[Y <= j]",   "P[Y<=j]"),
                  link, earg = earg)),
  infos = eval(substitute(function(...) {
    list(M1 = NA,  # zz -1?
         Q1 = NA,
         expected = TRUE,
         hadof = TRUE,
         multipleResponses = .multiple.responses ,
         parameters.names = as.character(NA),
         parallel = .parallel ,
         reverse = .reverse ,
         whitespace = .whitespace ,
         link =  .link )
  }, list( .link = link,
           .parallel = parallel,
           .multiple.responses = multiple.responses,
           .reverse = reverse,
           .whitespace = whitespace ))),

  constraints = eval(substitute(expression({
    if ( .multiple.responses ) {
      if ( !length(constraints) ) {
        Llevels <- extra$Llevels
        NOS <- extra$NOS
        Hk.matrix <- kronecker(diag(NOS), matrix(1,Llevels-1,1))
        constraints <- cm.VGAM(Hk.matrix, x = x,
                               bool = .parallel ,
                               apply.int = .apply.parint ,
                               constraints = constraints)
      }
    } else {
      constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                             bool = .parallel ,
                             apply.int = .apply.parint ,
                             constraints = constraints)
    }
  }), list( .parallel = parallel,
            .multiple.responses = multiple.responses,
            .apply.parint = apply.parint ))),
  deviance = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {

    answer <-
    if ( .multiple.responses ) {
      totdev <- 0
      NOS <- extra$NOS
      Llevels <- extra$Llevels
      for (iii in 1:NOS) {
        cindex <- (iii-1)*(Llevels-1) + 1:(Llevels-1)
        aindex <- (iii-1)*(Llevels  ) + 1:(Llevels)
        totdev <- totdev +
                  Deviance.categorical.data.vgam(
                    mu = mu[, aindex, drop = FALSE],
                    y = y[, aindex, drop = FALSE], w = w,
                    residuals = residuals,
                    eta = eta[, cindex, drop = FALSE],
                    extra = extra,
                    summation = TRUE)
      }
      totdev
    } else {
      Deviance.categorical.data.vgam(mu = mu, y = y, w = w,
                                     residuals = residuals,
                                     eta = eta, extra = extra,
                                     summation = TRUE)
    }
    answer
  }, list( .earg = earg, .link = link,
           .multiple.responses = multiple.responses ) )),

  initialize = eval(substitute(expression({

    if (colnames(x)[1] != "(Intercept)")
      warning("there seems to be no intercept term!")


    if (is.factor(y) && !is.ordered(y))
      warning("response should be ordinal---see ordered()")


    extra$multiple.responses <- .multiple.responses
    if ( .multiple.responses ) {
      checkCut(y)  # Check the input; stops if there is an error.
      if (any(w != 1) || NCOL(w) != 1)
          stop("the 'weights' argument must be a vector of all ones")
      Llevels <- max(y)
      delete.zero.colns <- FALSE
      orig.y <- cbind(y)  # Convert y into a matrix if necessary
      NOS <- NCOL(orig.y)
      use.y <- use.mustart <- NULL
      for (iii in 1:NOS) {
        y <- as.factor(orig.y[,iii])
        eval(process.categorical.data.VGAM)
        use.y <- cbind(use.y, y)
        use.mustart <- cbind(use.mustart, mustart)
      }
      mustart <- use.mustart
      y <- use.y  # n x (Llevels*NOS)
      M <- NOS * (Llevels-1)
      mynames <- y.names <- NULL
      for (iii in 1:NOS) {
        Y.names <- paste("Y", iii, sep = "")
        mu.names <- paste("mu", iii, ".", sep = "")
        mynames <- c(mynames, if ( .reverse )
          paste("P[", Y.names, ">=", 2:Llevels,     "]", sep = "") else
          paste("P[", Y.names, "<=", 1:(Llevels-1), "]", sep = ""))
        y.names <- c(y.names, param.names(mu.names, Llevels))
      }

      predictors.names <-
        namesof(mynames, .link , short = TRUE, earg = .earg )

      extra$NOS <- NOS
      extra$Llevels <- Llevels
  } else {


      delete.zero.colns <- TRUE

      eval(process.categorical.data.VGAM)
      M <- ncol(y) - 1
      mynames <- if ( .reverse )
        paste("P[Y", .fillerChar , ">=", .fillerChar,
              2:(1+M), "]", sep = "") else
        paste("P[Y", .fillerChar , "<=", .fillerChar,
              1:M,     "]", sep = "")

      predictors.names <-
        namesof(mynames, .link , short = TRUE, earg = .earg )
      y.names <- param.names("mu", M+1)

      if (NCOL(w) == 1) {
          if (length(mustart) && all(c(y) %in% c(0, 1)))
            for (iii in 1:ncol(y))
                mustart[,iii] <- weighted.mean(y[,iii], w)
      }

    extra$colnames.y  <- colnames(y)
  }
  }), list( .reverse = reverse, .multiple.responses = multiple.responses,
            .link = link, .earg = earg,
            .fillerChar = fillerChar,
            .whitespace = whitespace ))),


    linkinv = eval(substitute( function(eta, extra = NULL) {
        answer <-
        if ( .multiple.responses ) {
          NOS <- extra$NOS
          Llevels <- extra$Llevels
          fv.mat <- matrix(0, nrow(eta), NOS * Llevels)
          for (iii in 1:NOS) {
            cindex <- (iii-1)*(Llevels-1) + 1:(Llevels-1)
            aindex <- (iii-1)*(Llevels) + 1:(Llevels)
            if ( .reverse ) {
              ccump <- cbind(1,
                             eta2theta(eta[, cindex, drop = FALSE],
                                       .link , earg = .earg ))
              fv.mat[, aindex] <-
                  cbind(-tapplymat1(ccump, "diff"),
                        ccump[, ncol(ccump)])
            } else {
              cump <- cbind(eta2theta(eta[, cindex, drop = FALSE],
                                      .link ,
                                      earg = .earg ),
                            1)
              fv.mat[, aindex] <-
                  cbind(cump[, 1], tapplymat1(cump, "diff"))
            }
          }
          label.cols.y(fv.mat, NOS = NOS,
                       colnames.y = if (is.null(extra$colnames.y))
                                    NULL else
                         rep_len(extra$colnames.y, ncol(fv.mat)))
        } else {
          fv.mat <- if ( .reverse ) {
            ccump <- cbind(1, eta2theta(eta, .link , earg = .earg ))
            cbind(-tapplymat1(ccump, "diff"), ccump[, ncol(ccump)])
          } else {
            cump <- cbind(eta2theta(eta, .link , earg = .earg ), 1)
            cbind(cump[, 1], tapplymat1(cump, "diff"))
          }
          label.cols.y(fv.mat, colnames.y = extra$colnames.y, NOS = 1)
        }
        answer
    }, list( .reverse = reverse,
             .link = link, .earg = earg,
             .multiple.responses = multiple.responses ))),

  last = eval(substitute(expression({
    if ( .multiple.responses ) {
      misc$link <- .link
      misc$earg <- list( .earg )

    } else {
      misc$link <- rep_len( .link , M)
      names(misc$link) <- mynames

      misc$earg <- vector("list", M)
      names(misc$earg) <- names(misc$link)
      for (ii in 1:M)
        misc$earg[[ii]] <- .earg

    }

    misc$fillerChar <- .fillerChar
    misc$whitespace <- .whitespace

    misc$parameters <- mynames
    misc$reverse <- .reverse
    misc$parallel <- .parallel
    misc$multiple.responses <- .multiple.responses
  }), list(
            .reverse = reverse, .parallel = parallel,
            .link = link, .earg = earg,
            .fillerChar = fillerChar,
            .multiple.responses = multiple.responses,
            .whitespace = whitespace ))),

  linkfun = eval(substitute( function(mu, extra = NULL) {
    answer <-
    if ( .multiple.responses ) {
      NOS <- extra$NOS
      Llevels <- extra$Llevels
      eta.matrix <- matrix(0, nrow(mu), NOS*(Llevels-1))
      for (iii in 1:NOS) {
        cindex <- (iii-1)*(Llevels-1) + 1:(Llevels-1)
        aindex <- (iii-1)*(Llevels) + 1:(Llevels)
        cump <- tapplymat1(as.matrix(mu[, aindex]), "cumsum")
        eta.matrix[, cindex] <-
            theta2eta(if ( .reverse ) 1-cump[, 1:(Llevels-1)] else
                  cump[, 1:(Llevels-1)], .link , earg = .earg )
      }
      eta.matrix
    } else {
      cump <- tapplymat1(as.matrix(mu), "cumsum")
      M <- NCOL(mu) - 1
      theta2eta(if ( .reverse ) 1-cump[, 1:M] else cump[, 1:M],
                .link , earg = .earg )
    }
    answer
  }, list(
           .link = link, .earg = earg,
           .reverse = reverse,
           .multiple.responses = multiple.responses ))),
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                 y * w # Convert proportions to counts
      nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)

      smallno <- 1.0e4 * .Machine$double.eps
      if (max(abs(ycounts - round(ycounts))) > smallno)
        warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dmultinomial(x = ycounts, size = nvec, prob = mu,
                      log = TRUE, dochecking = FALSE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("cumulative", "VGAMordinal", "VGAMcategorical"),


  hadof = eval(substitute(
  function(eta, extra = list(),
           linpred.index = 1,
           w = 1, dim.wz = c(NROW(eta), NCOL(eta) * (NCOL(eta)+1)/2), 
           deriv = 1, ...) {
  if ( .multiple.responses )
    return(NA_real_, dim.wz[1], dim.wz[2])



  cumulative.eim.deriv1 <- function(mu, jay, w, reverse = FALSE) {
    M <- ncol(mu) - 1
    wz1 <- matrix(0, NROW(mu), M + M-1)  # Tridiagonal
    wz1[, iam(jay, jay, M = M)] <-  1 / (mu[, jay+1])^2 -
                                    1 / (mu[, jay  ])^2
    if (1 <= jay-1)
      wz1[, iam(jay-1, jay-1, M = M)] <-  -1 / (mu[, jay  ])^2
    if (jay+1 <= M)
      wz1[, iam(jay+1, jay+1, M = M)] <-   1 / (mu[, jay+1])^2
    if (1 < M && jay+1 <= M)
      wz1[, iam(jay, jay+1, M = M)] <-    -1 / (mu[, jay+1])^2
    if (1 < M && 1 <= jay-1)
      wz1[, iam(jay-1, jay, M = M)] <-     1 / (mu[, jay  ])^2
    (if (reverse) -c(w) else c(w)) * wz1
  }  # cumulative.eim.deriv1



  cumulative.eim.deriv2 <- function(mu, jay, w) {
    M <- ncol(mu) - 1
    wz2 <- matrix(0, NROW(mu), M + M-1)  # Tridiagonal
    wz2[, iam(jay, jay, M = M)] <-   1 / (mu[, jay+1])^3 +
                                     1 / (mu[, jay  ])^3
    if (1 <= jay-1)
      wz2[, iam(jay-1, jay-1, M = M)] <-  1 / (mu[, jay  ])^3
    if (jay+1 <= M)
      wz2[, iam(jay+1, jay+1, M = M)] <-  1 / (mu[, jay+1])^3
    if (1 < M && jay+1 <= M)
      wz2[, iam(jay, jay+1, M = M)] <-   -1 / (mu[, jay+1])^3
    if (1 < M && 1 <= jay-1)
      wz2[, iam(jay-1, jay, M = M)] <-   -1 / (mu[, jay  ])^3
    2 * c(w) * wz2
  }  # cumulative.eim.deriv2



    probs <- if ( .reverse ) {
        ccump <- cbind(1, eta2theta(eta, .link , earg = .earg ))
        cbind(-tapplymat1(ccump, "diff"), ccump[, ncol(ccump)])
      } else {
        cump <- cbind(eta2theta(eta, .link , earg = .earg ), 1)
        cbind(cump[, 1], tapplymat1(cump, "diff"))
      }
    mu.use <- pmax(probs, .Machine$double.eps * 1.0e-0)

    switch(as.character(deriv),
       "0" =  {
         M <- ncol(eta)
         wz <- c(w) * (1 / mu.use[, 1:M] + 1 / mu.use[, -1])
         if (M > 1)
           wz <- cbind(wz, -c(w) / mu.use[, 2:M])
         wz
       },
       "1" = {
         cumulative.eim.deriv1(mu = mu.use, jay = linpred.index, w = w,
                               reverse = .reverse )
       },
       "2" = {
         cumulative.eim.deriv2(mu = mu.use, jay = linpred.index, w = w)
       },
       stop("argument 'deriv' must be 0 or 1 or 2"))
  }, list( .link = link, .earg = earg,
           .reverse = reverse,
           .multiple.responses = multiple.responses ))),




  validparams = eval(substitute(function(eta, y, extra = NULL) {
    if ( .multiple.responses ) {
      return(TRUE)
    }
    probs <-
      if ( .reverse ) {
        ccump <- cbind(1, eta2theta(eta, .link , earg = .earg ))
        cbind(-tapplymat1(ccump, "diff"), ccump[, ncol(ccump)])
      } else {
        cump <- cbind(eta2theta(eta, .link , earg = .earg ), 1)
        cbind(cump[, 1], tapplymat1(cump, "diff"))
      }
    okay1 <- all(is.finite(probs)) && all(0 < probs & probs < 1)
    if (!okay1)
      warning("It seems that the nonparallelism assumption has ",
              "resulted in intersecting linear/additive predictors. ",
              "Try propodds() or fitting a partial nonproportional ",
              "odds model or choosing ",
              "some other link function, etc.")
    okay1
  }, list( .link = link, .earg = earg,
           .reverse = reverse,
           .multiple.responses = multiple.responses ))),

  deriv = eval(substitute(expression({
    mu.use <- pmax(mu, .Machine$double.eps * 1.0e-0)
    deriv.answer <-
    if ( .multiple.responses ) {
      NOS <- extra$NOS
      Llevels <- extra$Llevels
      dcump.deta <- resmat <- matrix(0, n, NOS * (Llevels-1))
      for (iii in 1:NOS) {
        cindex <- (iii-1)*(Llevels-1) + 1:(Llevels-1)
        aindex <- (iii-1)*(Llevels)   + 1:(Llevels-1)
        cump <- eta2theta(eta[,cindex, drop = FALSE],
                         .link , earg = .earg )
        dcump.deta[,cindex] <- dtheta.deta(cump, .link , earg = .earg )
        resmat[,cindex] <-
          (y[,   aindex, drop = FALSE] /mu.use[,   aindex, drop = FALSE] -
           y[, 1+aindex, drop = FALSE] /mu.use[, 1+aindex, drop = FALSE])
      }
      (if ( .reverse ) -c(w)  else c(w)) * dcump.deta * resmat
    } else {
      cump <- eta2theta(eta, .link , earg = .earg )
      dcump.deta <- dtheta.deta(cump, .link , earg = .earg )
      c(if ( .reverse ) -c(w)  else c(w)) *
      dcump.deta *
      (y[, -(M+1)] / mu.use[, -(M+1)] - y[, -1] / mu.use[, -1])
    }
    deriv.answer
  }), list( .link = link, .earg = earg,
            .reverse = reverse,
            .multiple.responses = multiple.responses ))),
  weight = eval(substitute(expression({
    if ( .multiple.responses ) {
      NOS <- extra$NOS
      Llevels <- extra$Llevels
      wz <- matrix(0, n, NOS*(Llevels-1))  # Diag elts only for a start
      for (iii in 1:NOS) {
        cindex <- (iii-1)*(Llevels-1) + 1:(Llevels-1)
        aindex <- (iii-1)*(Llevels)   + 1:(Llevels-1)
        wz[, cindex] <- c(w) * dcump.deta[, cindex, drop = FALSE]^2 *
                        (1 / mu.use[,   aindex, drop = FALSE] +
                         1 / mu.use[, 1+aindex, drop = FALSE])
      }
      if (Llevels-1 > 1) {
        iii <- 1
        oindex <- (iii-1) * (Llevels-1) + 1:(Llevels-2)
        wz <- cbind(wz, -c(w) *
              dcump.deta[, oindex] * dcump.deta[, 1+oindex])


        if (NOS > 1) {
          cptrwz <- ncol(wz)  # Like a pointer
          wz <- cbind(wz, matrix(0, nrow(wz), (NOS-1) * (Llevels-1)))
          for (iii in 2:NOS) {
            oindex <- (iii-1)*(Llevels-1) + 1:(Llevels-2)
            wz[,cptrwz + 1 + (1:(Llevels-2))] =
                  -c(w) * dcump.deta[,oindex] *
                       dcump.deta[, 1+oindex]
            cptrwz <- cptrwz + Llevels - 1 # Move it along a bit
          }
        }



        }
    } else {
      wz <- c(w) * dcump.deta^2 * (1/mu.use[, 1:M] + 1/mu.use[, -1])
      if (M > 1)
        wz <- cbind(wz,
                    -c(w) * dcump.deta[, -M] *
                    dcump.deta[, 2:M] / mu.use[, 2:M])
    }
    wz
  }), list( .earg = earg, .link = link,
            .multiple.responses = multiple.responses ))))
}  # cumulative()




 propodds <- function(reverse = TRUE, whitespace = FALSE) {
  if (!is.logical(reverse) ||
      length(reverse) != 1)
    stop("argument 'reverse' must be a single logical")

   cumulative(parallel = TRUE, reverse = reverse,
              whitespace = whitespace)
}



 acat <-
  function(link = "loglink", parallel = FALSE,
           reverse = FALSE, zero = NULL, whitespace = FALSE) {


  link <- as.list(substitute(link))
  earg  <- link2list(link)
  link <- attr(earg, "function.name")

  if (!is.logical(reverse) || length(reverse) != 1)
    stop("argument 'reverse' must be a single logical")

  stopifnot(is.logical(whitespace) &&
            length(whitespace) == 1)
  fillerChar <- ifelse(whitespace, " ", "")


  new("vglmff",
  blurb = c("Adjacent-categories model\n\n",
            "Links:    ",
    namesof(if (reverse)
                ifelse(whitespace, "P[Y = j] / P[Y = j + 1]",
                       "P[Y=j]/P[Y=j+1]") else
                ifelse(whitespace, "P[Y = j + 1] / P[Y = j]",
                                   "P[Y=j+1]/P[Y=j]"),
            link, earg = earg), "\n",
            "Variance: ",
            ifelse(whitespace,
            "mu[,j] * (1 - mu[,j]); -mu[,j] * mu[,k]",
            "mu[,j]*(1-mu[,j]); -mu[,j]*mu[,k]")),
  infos = eval(substitute(function(...) {
    list(M1 = NA,  # zz -1?
         Q1 = NA,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = as.character(NA),
         parallel = .parallel ,
         reverse = .reverse ,
         whitespace = .whitespace ,
         zero = .zero ,
         link = .link )
  }, list( .link = link,
           .zero = zero,
           .parallel = parallel,
           .reverse = reverse,
           .whitespace = whitespace ))),

  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints = constraints)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero ,
                     M = M, M1 = M,
                     predictors.names = predictors.names)
  }), list( .parallel = parallel, .zero = zero ))),

  deviance = Deviance.categorical.data.vgam,

  initialize = eval(substitute(expression({

    if (is.factor(y) && !is.ordered(y))
      warning("response should be ordinal---see ordered()")



    delete.zero.colns <- TRUE
    eval(process.categorical.data.VGAM)
    M <- ncol(y) - 1
    mynames <- if ( .reverse )
      paste("P[Y", .fillerChar , "=",
            1:M, "]", .fillerChar , "/", .fillerChar ,
            "P[Y", .fillerChar , "=", .fillerChar , 2:(M+1), "]",
            sep = "") else
      paste("P[Y", .fillerChar , "=", .fillerChar , 2:(M+1), "]",
            .fillerChar , "/", .fillerChar , "P[Y", .fillerChar ,
            "=", .fillerChar , 1:M,     "]", sep = "")

    predictors.names <-
      namesof(mynames, .link , short = TRUE, earg = .earg )
    y.names <- param.names("mu", M+1)

    extra$colnames.y  <- colnames(y)
  }), list( .earg = earg, .link = link, .reverse = reverse,
            .fillerChar = fillerChar,
            .whitespace = whitespace ))),

  linkinv = eval(substitute( function(eta, extra = NULL) {
    if (!is.matrix(eta))
      eta <- as.matrix(eta)
    M <- ncol(eta)
    fv.mat <- if ( .reverse ) {
      zetar <- eta2theta(eta, .link , earg = .earg )
      temp <- tapplymat1(zetar[, M:1], "cumprod")[, M:1, drop = FALSE]
      cbind(temp, 1) / drop(1 + temp %*% rep(1, ncol(temp)))
    } else {
      zeta <- eta2theta(eta, .link , earg = .earg )
      temp <- tapplymat1(zeta, "cumprod")
      cbind(1, temp) / drop(1 + temp %*% rep(1, ncol(temp)))
    }
   label.cols.y(fv.mat, colnames.y = extra$colnames.y, NOS = 1)
  }, list( .earg = earg, .link = link, .reverse = reverse) )),

  last = eval(substitute(expression({
    misc$link <- rep_len( .link , M)
    names(misc$link) <- mynames

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$parameters <- mynames
    misc$reverse <- .reverse
    misc$fillerChar <- .fillerChar
    misc$whitespace <- .whitespace
  }), list( .earg = earg, .link = link, .reverse = reverse,
            .fillerChar = fillerChar,
            .whitespace = whitespace ))),
  linkfun = eval(substitute( function(mu, extra = NULL) {
    M <- ncol(mu) - 1
    theta2eta(if ( .reverse ) mu[, 1:M]  / mu[,  -1] else
                              mu[,  -1]  / mu[, 1:M],
              .link , earg = .earg )
  }, list( .earg = earg, .link = link, .reverse = reverse) )),
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ycounts <- if (is.numeric(extra$orig.w))
                 y * w / extra$orig.w else
                 y * w  # Convert proportions to counts
      nvec <- if (is.numeric(extra$orig.w))
                  round(w / extra$orig.w) else round(w)

      smallno <- 1.0e4 * .Machine$double.eps
      if (max(abs(ycounts - round(ycounts))) > smallno)
        warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dmultinomial(x = ycounts, size = nvec, prob = mu,
                      log = TRUE, dochecking = FALSE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("acat", "VGAMordinal", "VGAMcategorical"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    if (!is.matrix(eta))
      eta <- as.matrix(eta)
    M <- ncol(eta)
    probs <- if ( .reverse ) {
      zetar <- eta2theta(eta, .link , earg = .earg )
      temp <- tapplymat1(zetar[, M:1], "cumprod")[, M:1, drop = FALSE]
      cbind(temp, 1) / drop(1 + temp %*% rep(1, ncol(temp)))
    } else {
      zeta <- eta2theta(eta, .link , earg = .earg )
      temp <- tapplymat1(zeta, "cumprod")
      cbind(1, temp) / drop(1 + temp %*% rep(1, ncol(temp)))
    }
    okay1 <- all(is.finite(probs)) && all(0 < probs & probs < 1)
    okay1
  }, list( .earg = earg, .link = link, .reverse = reverse) )),
  deriv = eval(substitute(expression({
    zeta <- eta2theta(eta, .link , earg = .earg )  # May be zetar

    dzeta.deta <- dtheta.deta(zeta, .link , earg = .earg )

    d1 <- acat.deriv(zeta, M = M, n = n, reverse = .reverse )
    score <- attr(d1, "gradient") / d1


    answer <- if ( .reverse ) {
      cumy <- tapplymat1(y, "cumsum")
      c(w) * dzeta.deta * (cumy[, 1:M] / zeta - score)
    } else {
      ccumy <- tapplymat1(y[, ncol(y):1], "cumsum")[, ncol(y):1]
      c(w) * dzeta.deta * (ccumy[, -1] / zeta - score)
    }


    answer
  }), list( .earg = earg, .link = link, .reverse = reverse) )),
  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, n, dimm(M))

    hess <- attr(d1, "hessian") / d1

    if (M > 1)
      for (jay in 1:(M-1))
        for (kay in (jay+1):M)
          wz[,iam(jay, kay,M)] <-
             (hess[, jay, kay] - score[, jay] * score[, kay]) *
             dzeta.deta[, jay] * dzeta.deta[, kay]
    if ( .reverse ) {
      cump <- tapplymat1(mu, "cumsum")
      wz[, 1:M] <- (cump[, 1:M] / zeta^2 - score^2) * dzeta.deta^2
    } else {
      ccump <- tapplymat1(mu[, ncol(mu):1], "cumsum")[, ncol(mu):1]
      wz[, 1:M] <- (ccump[, -1] / zeta^2 - score^2) * dzeta.deta^2
    }
    c(w) * wz
  }), list( .earg = earg, .link = link, .reverse = reverse ))))
}  # acat()



acat.deriv <- function(zeta, reverse, M, n) {

  alltxt <- NULL
  for (ii in 1:M) {
    index <- if (reverse) ii:M else 1:ii
    vars <- paste("zeta", index, sep = "")
    txt <- paste(vars, collapse = "*")
    alltxt <- c(alltxt, txt)
  }
  alltxt <- paste(alltxt, collapse = " + ")
  alltxt <- paste(" ~ 1 +", alltxt)
  txt <- as.formula(alltxt)

  allvars <- param.names("zeta", M)
  d1 <- deriv3(txt, allvars, hessian = TRUE)

  zeta <- as.matrix(zeta)
  for (ii in 1:M)
    assign(paste("zeta", ii, sep = ""), zeta[, ii])

  ans <- eval(d1)
  ans
}  # acat.deriv()






 brat <- function(refgp = "last",
                  refvalue = 1,
                  ialpha = 1) {
  if (!is.Numeric(ialpha, positive = TRUE))
    stop("'ialpha' must contain positive values only")
  if (!is.Numeric(refvalue, length.arg = 1, positive = TRUE))
    stop("'refvalue' must be a single positive value")

  if (!is.character(refgp) &&
      !is.Numeric(refgp, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE))
    stop("'refgp' must be a single positive integer")


  new("vglmff",
  blurb = c(paste("Bradley-Terry model (without ties)\n\n"),
            "Links:   ",
            namesof("alpha's", "loglink")),
  infos = eval(substitute(function(...) {
    list(M1 = NA,  # zz -1?
         Q1 = NA,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = as.character(NA),
         refvalue = .refvalue ,
         refgp = .refgp ,
         ialpha = .ialpha )
  }, list( .ialpha = ialpha,
           .refgp = refgp,
           .refvalue = refvalue ))),

  initialize = eval(substitute(expression({
    are.ties <- attr(y, "are.ties")  # If Brat() was used
    if (is.logical(are.ties) && are.ties)
        stop("use bratt(), not brat(), when there are ties")

    try.index <- 1:400
    M <- (seq_along(try.index))[(try.index+1)*(try.index) == ncol(y)]
    if (!is.finite(M))
      stop("cannot determine 'M'")
    ialpha <- matrix(rep_len( .ialpha , M), n, M, byrow = TRUE)
    etastart <- matrix(theta2eta(ialpha, "loglink",
                                 earg = list(theta = NULL)),
                       n, M, byrow = TRUE)
    refgp <- .refgp
    if (!intercept.only)
      warning("this function only works with intercept-only models")
    extra$ybrat.indices <- .brat.indices(NCo = M+1, are.ties = FALSE)
    uindex <- if ( .refgp == "last") 1:M else (1:(M+1))[-( .refgp ) ]

    predictors.names <-
      namesof(paste("alpha", uindex, sep = ""), "loglink",
              short = TRUE)

  }), list( .refgp = refgp, .ialpha = ialpha ))),

  linkinv = eval(substitute( function(eta, extra = NULL) {
    probs <- NULL
    eta <- as.matrix(eta)  # in case M = 1; prior to 20171227
    for (ii in 1:nrow(eta)) {
      alpha <- .brat.alpha(eta2theta(eta[ii, ], "loglink",
                                     earg = list(theta = NULL)),
                           .refvalue , .refgp )
      alpha1 <- alpha[extra$ybrat.indices[, "rindex"]]
      alpha2 <- alpha[extra$ybrat.indices[, "cindex"]]
      probs <- rbind(probs, alpha1 / (alpha1 + alpha2))
    }
    if (NROW(probs) == NROW(eta) &&
        NCOL(probs) == NCOL(eta))
      dimnames(probs) <- dimnames(eta)
    probs
  }, list( .refgp = refgp, .refvalue = refvalue) )),

  last = eval(substitute(expression({
    misc$link <- rep_len("loglink", M)
    names(misc$link) <- paste("alpha", uindex, sep = "")

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- list(theta = NULL)

    misc$refgp <- .refgp
    misc$refvalue <- .refvalue
  }), list( .refgp = refgp, .refvalue = refvalue ))),

  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ycounts <- if (is.numeric(extra$orig.w))
                 y * w / extra$orig.w else
                 y * w  # Convert proportions to counts
      nvec <- if (is.numeric(extra$orig.w))
              round(w / extra$orig.w) else round(w)

      smallno <- 1.0e4 * .Machine$double.eps
      if (max(abs(ycounts - round(ycounts))) > smallno)
          warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dmultinomial(x = ycounts, size = nvec, prob = mu,
                      log = TRUE, dochecking = FALSE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("brat", "VGAMcategorical"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    probs <- NULL
    eta <- as.matrix(eta)  # in case M = 1
    for (ii in 1:nrow(eta)) {
      alpha <- .brat.alpha(eta2theta(eta[ii, ], "loglink",
                                     earg = list(theta = NULL)),
                           .refvalue , .refgp )
      alpha1 <- alpha[extra$ybrat.indices[, "rindex"]]
      alpha2 <- alpha[extra$ybrat.indices[, "cindex"]]
      probs <- rbind(probs, alpha1 / (alpha1 + alpha2))
    }


    okay1 <- all(is.finite(probs)) && all(0 < probs & probs < 1)
    okay1
  }, list( .refvalue = refvalue, .refgp = refgp) )),
  deriv = eval(substitute(expression({
    ans <- NULL
    uindex <- if ( .refgp == "last") 1:M else (1:(M+1))[-( .refgp ) ]
    eta <- as.matrix(eta)  # in case M = 1
    for (ii in 1:nrow(eta)) {
      alpha <- .brat.alpha(eta2theta(eta[ii, ], "loglink",
                                     earg = list(theta = NULL)),
                           .refvalue, .refgp )
      ymat <- InverseBrat(y[ii, ], NCo = M+1, diag = 0)
      answer <- rep_len(0, M)
      for (aa in 1:(M+1)) {
        answer <- answer + (1 - (aa == uindex)) *
                  (ymat[uindex, aa] * alpha[aa] - ymat[aa, uindex] *
                  alpha[uindex]) / (alpha[aa] + alpha[uindex])
      }
      ans <- rbind(ans, w[ii] * answer)
    }
    dimnames(ans) <- dimnames(eta)
    ans
  }), list( .refvalue = refvalue, .refgp = refgp) )),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, dimm(M))
    for (ii in 1:nrow(eta)) {
      alpha <- .brat.alpha(eta2theta(eta[ii, ], "loglink",
                                     earg = list(theta = NULL)),
                          .refvalue, .refgp)
      ymat <- InverseBrat(y[ii, ], NCo = M+1, diag = 0)
      for (aa in 1:(M+1)) {
        wz[ii, 1:M] <- wz[ii, 1:M] + (1 - (aa == uindex)) *
                       (ymat[aa, uindex] + ymat[uindex, aa]) *
            alpha[aa] * alpha[uindex] / (
            alpha[aa] + alpha[uindex])^2
      }
      if (M > 1) {
        ind5 <- iam(1, 1, M, both = TRUE, diag = FALSE)
        wz[ii, (M+1):ncol(wz)] <-
          -(ymat[cbind(uindex[ind5$row], uindex[ind5$col])] +
            ymat[cbind(uindex[ind5$col], uindex[ind5$row])]) *
             alpha[uindex[ind5$col]] * alpha[uindex[ind5$row]] / (
             alpha[uindex[ind5$row]] + alpha[uindex[ind5$col]])^2
      }
    }
    wz <- c(w) * wz
    wz
  }), list( .refvalue = refvalue, .refgp = refgp ))))
}  # brat()







 bratt <- function(refgp = "last",
                   refvalue = 1,
                   ialpha = 1,
                   i0 = 0.01) {
  if (!is.Numeric(i0, length.arg = 1, positive = TRUE))
    stop("'i0' must be a single positive value")
  if (!is.Numeric(ialpha, positive = TRUE))
    stop("'ialpha' must contain positive values only")
  if (!is.Numeric(refvalue, length.arg = 1, positive = TRUE))
    stop("'refvalue' must be a single positive value")

  if (!is.character(refgp) &&
     !is.Numeric(refgp, length.arg = 1,
                 integer.valued = TRUE, positive = TRUE))
    stop("'refgp' must be a single positive integer")


  new("vglmff",
  blurb = c(paste("Bradley-Terry model (with ties)\n\n"),
            "Links:   ",
            namesof("alpha's", "loglink"), ", log(alpha0)"),
  infos = eval(substitute(function(...) {
    list(M1 = NA,  # zz -1?
         Q1 = NA,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = as.character(NA),
         refvalue = .refvalue ,
         refgp = .refgp ,
         i0 = .i0 ,
         ialpha = .ialpha )
  }, list( .ialpha = ialpha,
           .i0 = i0,
           .refgp = refgp,
           .refvalue = refvalue ))),

  initialize = eval(substitute(expression({
    try.index <- 1:400
    M <- (seq_along(try.index))[(try.index*(try.index-1)) == ncol(y)]
    if (!is.Numeric(M, length.arg = 1, integer.valued = TRUE))
      stop("cannot determine 'M'")
    NCo <- M  # Number of contestants

    are.ties <- attr(y, "are.ties")  # If Brat() was used
    if (is.logical(are.ties)) {
      if (!are.ties)
        stop("use brat(), not bratt(), when there are no ties")
      ties <- attr(y, "ties")
    } else {
      are.ties <- FALSE
      ties <- 0 * y
    }

    ialpha <- rep_len( .ialpha, NCo-1)
    ialpha0 <- .i0
    etastart <-
      cbind(matrix(theta2eta(ialpha,
                             "loglink",
                             list(theta = NULL)),
                   n, NCo-1, byrow = TRUE),
            theta2eta(rep_len(ialpha0, n), "loglink",
                      list(theta = NULL)))
    refgp <- .refgp
    if (!intercept.only)
      warning("this function only works with intercept-only models")
    extra$ties <- ties # Flat (1-row) matrix
    extra$ybrat.indices <- .brat.indices(NCo = NCo, are.ties = FALSE)
    extra$tbrat.indices <- .brat.indices(NCo = NCo, are.ties = TRUE)
    extra$dnties <- dimnames(ties)
    uindex <- if (refgp == "last") 1:(NCo-1) else (1:(NCo))[-refgp ]

    predictors.names <- c(
      namesof(paste("alpha", uindex, sep = ""), "loglink",
              short = TRUE),
      namesof("alpha0", "loglink", short = TRUE))
  }), list( .refgp = refgp,
           .i0 = i0,
           .ialpha = ialpha ))),

  linkinv = eval(substitute( function(eta, extra = NULL) {
    probs <- qprobs <- NULL
    M <- ncol(eta)
    for (ii in 1:nrow(eta)) {
      alpha <- .brat.alpha(eta2theta(eta[ii, -M],
                                     "loglink"),
                           .refvalue , .refgp )
      alpha0 <- loglink(eta[ii, M], inverse = TRUE)
      alpha1 <- alpha[extra$ybrat.indices[, "rindex"]]
      alpha2 <- alpha[extra$ybrat.indices[, "cindex"]]
       probs <- rbind( probs, alpha1 / (alpha1 + alpha2 + alpha0))  #
      qprobs <- rbind(qprobs, alpha0 / (alpha1 + alpha2 + alpha0))  #
    }
    if (length(extra$dnties))
      dimnames(qprobs) <- extra$dnties
    attr(probs, "probtie") <- qprobs
    probs
  }, list( .refgp = refgp, .refvalue = refvalue) )),
  last = eval(substitute(expression({
    misc$link <- rep_len("loglink", M)
    names(misc$link) <- c(paste("alpha", uindex, sep = ""), "alpha0")


    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- list(theta = NULL)


    misc$refgp <- .refgp
    misc$refvalue <- .refvalue
    misc$alpha  <- alpha
    misc$alpha0 <- alpha0
  }), list( .refgp = refgp, .refvalue = refvalue ))),
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * (y * log(mu) +
                         0.5 * extra$ties * log(attr(mu, "probtie")))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("bratt", "VGAMcategorical"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    probs <- qprobs <- NULL
    M <- ncol(eta)
    for (ii in 1:nrow(eta)) {
      alpha <- .brat.alpha(eta2theta(eta[ii, -M],
                                     "loglink"),
                           .refvalue , .refgp )
      alpha0 <- loglink(eta[ii, M], inverse = TRUE)
      alpha1 <- alpha[extra$ybrat.indices[, "rindex"]]
      alpha2 <- alpha[extra$ybrat.indices[, "cindex"]]
       probs <- rbind( probs, alpha1 / (alpha1 + alpha2 + alpha0))
      qprobs <- rbind(qprobs, alpha0 / (alpha1 + alpha2 + alpha0))
    }
    okay1 <- all(is.finite( probs)) && all(0 <  probs &  probs < 1) &&
             all(is.finite(qprobs)) && all(0 < qprobs & qprobs < 1)
    okay1
  }, list( .refvalue = refvalue, .refgp = refgp) )),
  deriv = eval(substitute(expression({
    ans <- NULL
    ties <- extra$ties
    NCo <- M
    uindex <- if ( .refgp == "last") 1:(M-1) else (1:(M))[-( .refgp )]
    eta <- as.matrix(eta)
    for (ii in 1:nrow(eta)) {
      alpha <- .brat.alpha(eta2theta(eta[ii, -M], "loglink",
                                     earg = list(theta = NULL)),
                           .refvalue, .refgp )
      alpha0 <- loglink(eta[ii, M], inverse = TRUE)
      ymat <- InverseBrat(   y[ii, ], NCo = M, diag = 0)
      tmat <- InverseBrat(ties[ii, ], NCo = M, diag = 0)
      answer <- rep_len(0, NCo-1)  # deriv wrt eta[-M]
      for (aa in 1:NCo) {
        Daj <- alpha[aa] + alpha[uindex] + alpha0
        pja <- alpha[uindex] / Daj
        answer <- answer + alpha[uindex] *
             (-ymat[aa, uindex] + ymat[uindex, aa] * (1 - pja) / pja -
             tmat[uindex, aa]) / Daj
      }
      deriv0 <- 0 # deriv wrt eta[M]
      for (aa in 1:(NCo-1))
        for (bb in (aa+1):NCo) {
          Dab <- alpha[aa] + alpha[bb] + alpha0
          qab <- alpha0 / Dab
          deriv0 <- deriv0 + alpha0 *
                    (-ymat[aa, bb] - ymat[bb,aa] +
                    tmat[aa, bb] * (1 - qab) / qab) / Dab
          }
        ans <- rbind(ans, w[ii] * c(answer, deriv0))
    }
    dimnames(ans) <- dimnames(eta)
    ans
  }), list( .refvalue = refvalue, .refgp = refgp) )),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, dimm(M))   # includes diagonal
    for (ii in 1:nrow(eta)) {
      alpha <- .brat.alpha(eta2theta(eta[ii, -M], "loglink",
                           earg = list(theta = NULL)),
                          .refvalue, .refgp)
      alpha0 <- loglink(eta[ii, M], inverse = TRUE)
      ymat <- InverseBrat(   y[ii, ], NCo = M, diag = 0)
      tmat <- InverseBrat(ties[ii, ], NCo = M, diag = 0)

      for (aa in 1:(NCo)) {
        Daj <- alpha[aa] + alpha[uindex] + alpha0
        pja <- alpha[uindex] / Daj
        nja <- ymat[aa,uindex] + ymat[uindex,aa] + tmat[uindex,aa]
        wz[ii, 1:(NCo-1)] <- wz[ii, 1:(NCo - 1)] +
                             alpha[uindex]^2 * nja *
                             (1 - pja) / (pja * Daj^2)
        if (aa < NCo)
          for (bb in (aa+1):(NCo)) {
            nab <- ymat[aa,bb] + ymat[bb,aa] + tmat[bb,aa]
            Dab <- alpha[aa] + alpha[bb] + alpha0
            qab <- alpha0 / Dab
            wz[ii, NCo] <- wz[ii,NCo] + alpha0^2 * nab *
                     (1-qab) / (qab * Dab^2)
          }
      }

      if (NCo > 2) {
        ind5 <- iam(1, 1, M = NCo, both = TRUE, diag = FALSE)
        alphajunk <- c(alpha, junk = NA)
        mat4 <- cbind(uindex[ind5$row],uindex[ind5$col])
        wz[ii,(M+1):ncol(wz)] <- -(ymat[mat4] + ymat[mat4[, 2:1]] +
           tmat[mat4]) * alphajunk[uindex[ind5$col]] *
           alphajunk[uindex[ind5$row]] / (alpha0 +
           alphajunk[uindex[ind5$row]] +
           alphajunk[uindex[ind5$col]])^2
      }
      for (sss in seq_along(uindex)) {
        jay <- uindex[sss]
        naj <- ymat[, jay] + ymat[jay, ] + tmat[, jay]
        Daj <- alpha[jay] + alpha + alpha0
        wz[ii, iam(sss, NCo, M = NCo, diag = TRUE)] <-
            -alpha[jay] * alpha0 * sum(naj / Daj^2)
      }
    }
    wz <- c(w) * wz
    wz
  }), list( .refvalue = refvalue, .refgp = refgp ))))
}  # bratt()



.brat.alpha <- function(vec, value, posn) {
  if (is.character(posn))
    if (posn != "last")
      stop("can only handle \"last\"") else return(c(vec, value))
  c(if (posn == 1) NULL else vec[1:(posn-1)], value,
    if (posn == length(vec) + 1) NULL else vec[posn:length(vec)])
}


.brat.indices <- function(NCo, are.ties = FALSE) {
  if (!is.Numeric(NCo, length.arg = 1,
                  integer.valued = TRUE) ||
      NCo < 2)
    stop("bad input for 'NCo'")
  m <- diag(NCo)
  if (are.ties) {
    cbind(rindex = row(m)[col(m) <  row(m)],
          cindex = col(m)[col(m) <  row(m)])
  } else
    cbind(rindex = row(m)[col(m) != row(m)],
          cindex = col(m)[col(m) != row(m)])
}



 Brat <- function(mat,
                  ties = 0 * mat,
                  string = c(">", "=="),
                  whitespace = FALSE) {


  stopifnot(is.logical(whitespace) &&
            length(whitespace) == 1)
  fillerChar <- ifelse(whitespace, " ", "")
  string <- paste(fillerChar, string, fillerChar, sep = "")


  allargs <- list(mat)  # ,...
  callit <- if (length(names(allargs))) names(allargs) else
            as.character(seq_along(allargs))
  ans <- ans.ties <- NULL
  for (ii in seq_along(allargs)) {
    m <- allargs[[ii]]
    if (!is.matrix(m) || dim(m)[1] != dim(m)[2])
      stop("m must be a square matrix")

    diag(ties) <- 0
    if (!all(ties == t(ties)))
      stop("ties must be a symmetric matrix")
    are.ties <- any(ties > 0)
    diag(ties) <- NA

    diag(m) <- 0  # Could have been NAs
    if (anyNA(m))
      stop("missing values not allowed (except on the diagonal)")
    diag(m) <- NA

    dm <- as.data.frame.table(m)
    dt <- as.data.frame.table(ties)
    dm <- dm[!is.na(dm$Freq), ]
    dt <- dt[!is.na(dt$Freq), ]
    usethis1 <- paste(dm[, 1], string[1], dm[, 2], sep = "")
    usethis2 <- paste(dm[, 1], string[2], dm[, 2], sep = "")
    ans <- rbind(ans, matrix(dm$Freq, nrow = 1))
    ans.ties <- rbind(ans.ties, matrix(dt$Freq, nrow = 1))
  }
  dimnames(ans) <- list(callit, usethis1)
  dimnames(ans.ties) <- list(callit, usethis2)
  attr(ans, "ties") <- ans.ties
  attr(ans, "are.ties") <- are.ties
  ans
}  # Brat()




InverseBrat <-
  function(yvec,
           NCo = (1:900)[(1:900)*((1:900)-1) == ncol(rbind(yvec))],
           multiplicity = if (is.matrix(yvec)) nrow(yvec) else 1,
           diag = NA,
           string = c(">", "=="),
           whitespace = FALSE) {



  stopifnot(is.logical(whitespace) &&
            length(whitespace) == 1)
  fillerChar <- ifelse(whitespace, " ", "")
  string <- paste(fillerChar, string, fillerChar, sep = "")


  ans <- array(diag, c(NCo, NCo, multiplicity))
  yvec.orig <- yvec
  yvec <- c(yvec)
  ptr <- 1
  for (mul in 1:multiplicity)
    for (i1 in 1:(NCo))
      for (i2 in 1:(NCo))
        if (i1 != i2) {
          ans[i2, i1, mul] <- yvec[ptr]
          ptr <- ptr + 1
        }
  ans <- if (multiplicity > 1) ans else matrix(ans, NCo, NCo)

  if (is.array(yvec.orig) ||
      is.matrix(yvec.orig)) {
    names.yvec <- dimnames(yvec.orig)[[2]]
    ii <- strsplit(names.yvec, string[1])
    cal <- NULL
    for (kk in c(NCo, 1:(NCo-1)))
      cal <- c(cal, (ii[[kk]])[1])
    if (multiplicity>1) {
      dimnames(ans) <- list(cal, cal, dimnames(yvec.orig)[[1]])
    } else {
      dimnames(ans) <- list(cal, cal)
    }
  }
  ans
}  # InverseBrat()






 ordpoisson <- function(cutpoints,
                        countdata = FALSE, NOS = NULL, Levels = NULL,
                        init.mu = NULL, parallel = FALSE,
                        zero = NULL,
                        link = "loglink") {

  link <- as.list(substitute(link))
  earg  <- link2list(link)
  link <- attr(earg, "function.name")




  fcutpoints <- cutpoints[is.finite(cutpoints)]
  if (!is.Numeric(fcutpoints, integer.valued = TRUE) ||
      any(fcutpoints < 0))
    stop("'cutpoints' must have non-negative integer or Inf ",
         "values only")
  if (is.finite(cutpoints[length(cutpoints)]))
    cutpoints <- c(cutpoints, Inf)

  if (!is.logical(countdata) || length(countdata) != 1)
    stop("argument 'countdata' must be a single logical")
  if (countdata) {
    if (!is.Numeric(NOS, integer.valued = TRUE, positive = TRUE))
      stop("'NOS' must have integer values only")
    if (!is.Numeric(Levels, integer.valued = TRUE,
                    positive = TRUE) ||
        any(Levels < 2))
      stop("'Levels' must have integer values (>= 2) only")
    Levels <- rep_len(Levels, NOS)
  }


  new("vglmff",
  blurb = c(paste("Ordinal Poisson model\n\n"),
            "Link:     ", namesof("mu", link, earg = earg)),


  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           apply.int = TRUE,
                           constraints = constraints)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero ,
                     M = M, M1 = 1,
                     predictors.names = predictors.names)
  }), list( .parallel = parallel, .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("mu"),
         lmu = .link ,
         zero = .zero )
  }, list( .zero = zero, .link = link ))),

  initialize = eval(substitute(expression({
    orig.y <- cbind(y)  # Convert y into a matrix if necessary
    if ( .countdata ) {
      extra$NOS <- M <- NOS <- .NOS
      extra$Levels <- Levels <- .Levels
      y.names <- dimnames(y)[[2]]  # Hopefully the user inputted them
    } else {
      if (any(w != 1) || NCOL(w) != 1)
        stop("the 'weights' argument must be a vector of all ones")
      extra$NOS <- M <- NOS <- if (is.Numeric( .NOS )) .NOS else
          ncol(orig.y)
      Levels <- rep_len(if (is.Numeric( .Levels )) .Levels else 0, NOS)
      if (!is.Numeric( .Levels ))
        for (iii in 1:NOS) {
          Levels[iii] <- length(unique(sort(orig.y[,iii])))
        }
      extra$Levels <- Levels
    }


    initmu <- if (is.Numeric( .init.mu ))
      rep_len( .init.mu , NOS) else NULL
    cutpoints <- rep_len( .cutpoints, sum(Levels))
    delete.zero.colns <- FALSE
    use.y <- if ( .countdata ) y else matrix(0, n, sum(Levels))
    use.etastart <- matrix(0, n, M)
    cptr <- 1
    for (iii in 1:NOS) {
        y <- factor(orig.y[,iii], levels=(1:Levels[iii]))
        if ( !( .countdata )) {
            eval(process.categorical.data.VGAM)  # Creates mustart and y
            use.y[,cptr:(cptr+Levels[iii]-1)] <- y
        }
        use.etastart[,iii] <- if (is.Numeric(initmu))
            initmu[iii] else
            median(cutpoints[cptr:(cptr+Levels[iii]-1-1)])
        cptr <- cptr + Levels[iii]
    }
    mustart <- NULL # Overwrite it
    etastart <- theta2eta(use.etastart, .link , earg = .earg )
    y <- use.y  # n x sum(Levels)
    M <- NOS
    for (iii in 1:NOS) {
        mu.names <- paste("mu", iii, ".", sep = "")
    }

    ncoly <- extra$ncoly <- sum(Levels)
    cp.vector <- rep_len( .cutpoints , ncoly)
    extra$countdata <- .countdata
    extra$cutpoints <- cp.vector
    extra$n <- n

    mynames <- param.names("mu", M, skip1 = TRUE)
    predictors.names <-
      namesof(mynames, .link , earg = .earg , tag = FALSE)
  }), list( .link = link, .countdata = countdata, .earg = earg,
            .cutpoints = cutpoints, .NOS = NOS, .Levels = Levels,
            .init.mu = init.mu))),
  linkinv = eval(substitute( function(eta, extra = NULL) {
    mu <- eta2theta(eta, link= .link , earg = .earg )  # Poisson means
    mu <- cbind(mu)
    mu
  }, list( .link = link, .earg = earg, .countdata = countdata ))),
  last = eval(substitute(expression({
    if ( .countdata ) {
      misc$link <- .link
      misc$earg <- list( .earg )
    } else {
      misc$link <- rep_len( .link , M)
      names(misc$link) <- mynames

      misc$earg <- vector("list", M)
      names(misc$earg) <- names(misc$link)
      for (ii in 1:M)
        misc$earg[[ii]] <- .earg
    }
    misc$parameters <- mynames
    misc$countdata <- .countdata
    misc$true.mu = FALSE    # $fitted is not a true mu
  }), list( .link = link, .countdata = countdata, .earg = earg ))),
    loglikelihood =
      function(mu, y, w, residuals = FALSE, eta, extra = NULL,
               summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      probs <- ordpoissonProbs(extra, mu)
      index0 <- y == 0
      probs[index0] <- 1
      pindex0 <- probs == 0
      probs[pindex0] <- 1
      if (summation) {
        sum(pindex0) * (-1.0e+10) + sum(w * y * log(probs))
      } else {
        stop("20140311; 'summation=F' not done yet")
      }
      }
  },
  vfamily = c("ordpoisson", "VGAMcategorical"),
  deriv = eval(substitute(expression({
    probs <- ordpoissonProbs(extra, mu)
    probs.use <- pmax(probs, .Machine$double.eps * 1.0e-0)

    cp.vector <- extra$cutpoints
    NOS <- extra$NOS
    Levels <- extra$Levels
    resmat <- matrix(0, n, M)
    dl.dprob <- y / probs.use
    dmu.deta <- dtheta.deta(mu, .link , earg = .earg )
    dprob.dmu <- ordpoissonProbs(extra, mu, deriv = 1)
    cptr <- 1
    for (iii in 1:NOS) {
      for (kkk in 1:Levels[iii]) {
       resmat[,iii] <- resmat[,iii] +
                      dl.dprob[,cptr] * dprob.dmu[,cptr]
       cptr <- cptr + 1
      }
    }
    resmat <- c(w) * resmat * dmu.deta
    resmat
  }), list( .link = link, .earg = earg, .countdata=countdata ))),
  weight = eval(substitute(expression({
    d2l.dmu2 <- matrix(0, n, M)  # Diagonal matrix
    cptr <- 1
    for (iii in 1:NOS) {
      for (kkk in 1:Levels[iii]) {
        d2l.dmu2[,iii] <- d2l.dmu2[,iii] +
            dprob.dmu[,cptr]^2 / probs.use[,cptr]
        cptr <- cptr + 1
      }
    }
    wz <- c(w) * d2l.dmu2 * dmu.deta^2
    wz
  }), list( .earg = earg, .link = link, .countdata = countdata ))))
}



ordpoissonProbs <- function(extra, mu, deriv = 0) {
  cp.vector <- extra$cutpoints
  NOS <- extra$NOS
  if (deriv == 1) {
    dprob.dmu <- matrix(0, extra$n, extra$ncoly)
  } else {
    probs <- matrix(0, extra$n, extra$ncoly)
  }
  mu <- cbind(mu)
  cptr <- 1
  for (iii in 1:NOS) {
    if (deriv == 1) {
      dprob.dmu[,cptr] <- -dpois(x = cp.vector[cptr], lambda = mu[,iii])
    } else {
      probs[,cptr] <- ppois(q = cp.vector[cptr], lambda = mu[,iii])
    }
    cptr <- cptr + 1
    while (is.finite(cp.vector[cptr])) {
      if (deriv == 1) {
        dprob.dmu[,cptr] <-
                dpois(x = cp.vector[cptr-1], lambda = mu[,iii]) -
                dpois(x = cp.vector[cptr  ], lambda = mu[,iii])
      } else {
        probs[,cptr] <-
                ppois(q = cp.vector[cptr  ], lambda = mu[,iii]) -
                ppois(q = cp.vector[cptr-1], lambda = mu[,iii])
      }
      cptr <- cptr + 1
    }
    if (deriv == 1) {
        dprob.dmu[,cptr] <-
                dpois(x = cp.vector[cptr-1], lambda = mu[,iii]) -
                dpois(x = cp.vector[cptr  ], lambda = mu[,iii])
    } else {
        probs[,cptr] <-
                ppois(q = cp.vector[cptr  ], lambda = mu[,iii]) -
                ppois(q = cp.vector[cptr-1], lambda = mu[,iii])
    }
    cptr <- cptr + 1
  }
  if (deriv == 1) dprob.dmu else probs
}









findFirstMethod <- function(methodsfn, charvec) {
  answer <- NULL
  for (ii in seq_along(charvec)) {
    if (existsMethod(methodsfn, signature(VGAMff = charvec[ii]))) {
      answer <- charvec[ii]
      break
    }
  }
  answer
}



margeff <- function(object, subset = NULL, ...) {


  try.this <- findFirstMethod("margeffS4VGAM", object@family@vfamily)
  if (length(try.this)) {
    margeffS4VGAM(object = object,
            subset = subset,
            VGAMff = new(try.this),
        ...)
  } else {
    stop("Could not find a methods function for 'margeffS4VGAM' ",
         "emanating from '", object@family@vfamily[1], "'")
  }
}





subsetarray3 <- function(array3, subset = NULL) {
  if (is.null(subset)) {
    return(array3)
  } else
    if (is.numeric(subset) && (length(subset) == 1)) {
      return(array3[, , subset])
  } else {
    return(array3[, , subset])
  }
  warning("argument 'subset' unmatched. Doing nothing")
  array3
}





setClass("VGAMcategorical",     contains = "vglmff")

setClass("VGAMordinal",         contains = "VGAMcategorical")
setClass("multinomial",         contains = "VGAMcategorical")

setClass("acat",                contains = "VGAMordinal")
setClass("cumulative",          contains = "VGAMordinal")
setClass("cratio",              contains = "VGAMordinal")
setClass("sratio",              contains = "VGAMordinal")



setMethod("margeffS4VGAM",
          signature(VGAMff = "VGAMcategorical"),
  function(object,
           subset = NULL,
           VGAMff,
           ...) {
  object@post$M <- M   <- object@misc$M
  object@post$n <- nnn <- object@misc$n
  invisible(object)
  })




setMethod("margeffS4VGAM",  signature(VGAMff ="multinomial"),
  function(object,
           subset = NULL,
           VGAMff,
           ...) {

  object <- callNextMethod(VGAMff = VGAMff,
                           object = object,
                           subset = subset,
                           ...)

  M   <- object@misc$M
  nnn <- object@misc$n
    cfit <- coefvlm(object, matrix.out = TRUE)
    rlev <- object@misc$refLevel
    if (!length(rlev))
      relev <- M+1  # Default
    Bmat <- matrix(0, nrow(cfit), 1 + ncol(cfit))
    Bmat[, -rlev] <- cfit
    ppp   <- nrow(Bmat)
    pvec1 <- fitted(object)[1, ]
    rownames(Bmat) <- rownames(cfit)
    colnames(Bmat) <- if (length(names(pvec1))) names(pvec1) else
                      param.names("mu", M+1)


    BB <- array(Bmat, c(ppp, M+1, nnn))
    pvec  <- c(t(fitted(object)))
    pvec  <- rep(pvec, each = ppp)
    temp1 <- array(BB * pvec, c(ppp, M+1, nnn))
    temp2 <- aperm(temp1, c(2, 1, 3))  # (M+1) x ppp x nnn
    temp2 <- colSums(temp2)  # ppp x nnn
    temp2 <- array(rep(temp2, each = M+1), c(M+1, ppp, nnn))
    temp2 <- aperm(temp2, c(2, 1, 3))  # ppp x (M+1) x nnn
    temp3 <- pvec
    ans.mlm <- array((BB - temp2) * temp3, c(ppp, M+1, nnn),
                     dimnames = list(dimnames(Bmat)[[1]],
                                     dimnames(Bmat)[[2]],
                                     dimnames(fitted(object))[[1]]))
    return(subsetarray3(ans.mlm, subset = subset))
  })




setMethod("margeffS4VGAM",  signature(VGAMff = "VGAMordinal"),
  function(object,
           subset = NULL,
           VGAMff,
           ...) {
  M   <- object@misc$M
  nnn <- object@misc$n

  object@post$reverse <- object@misc$reverse
  object@post$linkfunctions <- linkfunctions <- object@misc$link
  object@post$all.eargs <- all.eargs <- object@misc$earg
  object@post$Bmat <- Bmat <- coefvlm(object, matrix.out = TRUE)
  object@post$ppp <- nrow(Bmat)
  etamat <- predict(object)

  hdot <- Thetamat <- etamat
  for (jlocal in 1:M) {
    Thetamat[, jlocal] <- eta2theta(etamat[, jlocal],
                                    link = linkfunctions[jlocal],
                                    earg = all.eargs[[jlocal]])
    hdot[, jlocal] <- dtheta.deta(Thetamat[, jlocal],
                                  link = linkfunctions[jlocal],
                                  earg = all.eargs[[jlocal]])
  }  # jlocal


  object@post$hdot <- hdot
  object@post$Thetamat <- Thetamat
  object
  })





setClass("tobit",          contains = "VGAMcategorical")



setMethod("margeffS4VGAM",  signature(VGAMff = "tobit"),
  function(object,
           subset = NULL,
           VGAMff,
           ...) {


  object <- callNextMethod(VGAMff = VGAMff,
                           object = object,
                           subset = subset,
                           ...)

  etamat <- predict(object)  # nnn x M
  nnn <- nrow(etamat)
  M   <- ncol(etamat)
  Lowmat <- object@misc$Lower
  Uppmat <- object@misc$Upper
  censoL <- object@extra$censoredL
  censoU <- object@extra$censoredU

  mustar <- fitted(object, type.fitted = "uncensored")
  M1 <- npred(object, type = "one.response")
  NOS <- M / M1
  sigmahat <- matrix(NA_real_, nnn, NOS)
  for (jay in seq(NOS))
    sigmahat[, jay] <-
      eta2theta(etamat[, M1 * jay],
                link = object@misc$link[[M1 * jay]],
                earg = object@misc$earg[[M1 * jay]])

  ymat <- depvar(object)
  cnymat <- colnames(ymat)
  rnymat <- rownames(ymat)
  ymat[censoL] <- Lowmat[censoL]
  ymat[censoU] <- Uppmat[censoU]
  zedd <- (ymat - mustar) / sigmahat  # nnn x NOS
  coefmat <- coef(object, matrix = TRUE)[, c(TRUE, FALSE),
                                         drop = FALSE]
  ppp <- nrow(coefmat)

  betamat <- matrix(c(coefmat), nnn, ppp * NOS, byrow = TRUE)
  ans <- betamat  *
         kronecker(dnorm(zedd) / sigmahat, matrix(1, 1, ppp))

  uncens <- !censoL && !censoU
  ans[uncens] <- ans[uncens] * zedd[uncens] / sigmahat[uncens]

  ans[censoL] <- -ans[censoL]


  dim(ans) <- c(nnn, ppp, NOS)
  ans <- aperm(ans, c(2, 3, 1))  # ppp x NOS x nnn
  dimnames(ans) <- list(rownames(coefmat), cnymat, rnymat)

  ans
})











setMethod("margeffS4VGAM",  signature(VGAMff = "cumulative"),
  function(object,
           subset = NULL,
           VGAMff,
           ...) {


  object <- callNextMethod(VGAMff = VGAMff,
                           object = object,
                           subset = subset,
                           ...)
  reverse <- object@post$reverse
  linkfunctions <- object@post$linkfunctions
  all.eargs <- object@post$all.eargs
  Bmat <- cfit <- object@post$Bmat
  ppp <- object@post$ppp
  etamat <- predict(object)  # nnn x M
  fitmat <- fitted(object)   # nnn x (M + 1)
  nnn <- nrow(etamat)
  M   <- ncol(etamat)
  hdot <- object@post$hdot
  Thetamat <- object@post$Thetamat




    hdot.big <- kronecker(hdot, matrix(1, ppp, 1))  # Enlarged
    resmat <- cbind(hdot.big, 1)
    resmat[, 1] <- ifelse(reverse, -1, 1) * hdot.big[, 1] * cfit[, 1]

    if (M > 1) {
      for (jlocal in 2:M) {
        resmat[, jlocal] <- ifelse(reverse, -1, 1) *
          (hdot.big[, jlocal    ] * cfit[, jlocal    ] -
           hdot.big[, jlocal - 1] * cfit[, jlocal - 1])
      }  # jlocal

    }  # if

    resmat[, M+1] <- ifelse(reverse, 1, -1) *
                     hdot.big[, M] * cfit[, M]

    ans.cum <- array(resmat, c(ppp, nnn, M+1),
                     dimnames = list(dimnames(Bmat)[[1]],
                                     dimnames(fitted(object))[[1]],
                                     dimnames(fitted(object))[[2]]))
    ans.cum <- aperm(ans.cum, c(1, 3, 2))  # ppp x (M+1) x nnn

    subsetarray3(ans.cum, subset = subset)
  })







setMethod("margeffS4VGAM",  signature(VGAMff = "acat"),
  function(object,
           subset = NULL,
           VGAMff,
           ...) {


  object <- callNextMethod(VGAMff = VGAMff,
                           object = object,
                           subset = subset,
                           ...)
  reverse <- object@post$reverse
  linkfunctions <- object@post$linkfunctions
  all.eargs <- object@post$all.eargs
  Bmat <- cfit <- object@post$Bmat
  ppp <- object@post$ppp
  etamat <- predict(object)  # nnn x M
  fitmat <- fitted(object)   # nnn x (M + 1)
  nnn <- nrow(etamat)
  M   <- ncol(etamat)
  hdot <- object@post$hdot
  Thetamat <- object@post$Thetamat




    expcs.etamat <- if (reverse)
      exp(tapplymat1(etamat[, M:1, drop = FALSE],
                     "cumsum")[, M:1, drop = FALSE]) else
      exp(tapplymat1(etamat, "cumsum"))
    csexpcs.etavec <- rowSums(expcs.etamat)



    if (!all(object@misc$link == "loglink"))
      stop("currently only the 'loglink' link is supported")


  acat.derivs <- function(jay, tee,
                          M, expcs.etamat, Thetamat,
                          prob1, probMplus1,
                          reverse = FALSE) {

    if (jay > M+1) stop("argument 'jay' out of range")
    if (M   < tee) stop("argument 'tee' out of range")

    if (reverse) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      dpMplus1.detat <- -(probMplus1^2) *
                        rowSums(expcs.etamat[, 1:tee, drop = FALSE])
      if (jay == M+1) {
        return(dpMplus1.detat)
      }
      if (jay <= tee) {
        return((probMplus1 + dpMplus1.detat) * expcs.etamat[, jay])
      }
      if (tee < jay) {
        return(dpMplus1.detat * expcs.etamat[, jay])
      }
    } else {  # reverse = FALSE ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      dp1.detat <-
        -(prob1^2) * rowSums(expcs.etamat[, tee:M, drop = FALSE])
      if (jay == 1) {
        return(dp1.detat)
      }
      if (jay <= tee) {
        return(dp1.detat * expcs.etamat[, jay-1])
      }
      if (tee < jay) {
        return((prob1 + dp1.detat) * expcs.etamat[, jay-1])
      }
    } # reverse ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  }  # acat.derivs



  A        <- array(0, c(i = nnn, vars = ppp, probs = M + 1, etas = M))
  ansarray <- array(0, c(vars = ppp, i = nnn, probs = M + 1))
  if (reverse) {
    probMplus1 <- 1 / (1 +  csexpcs.etavec) # Last level of Y
  } else {
    prob1 <- 1 / (1 +  csexpcs.etavec)  # First level of Y
  }

    for (jlocal in 1:(M+1)) {
      for (tlocal in 1:M) {
        A[, , jlocal, tlocal] <-
          acat.derivs(jay = jlocal, tee = tlocal,
                      M = M, expcs.etamat = expcs.etamat,
                      Thetamat = Thetamat,
                      prob1 = prob1, probMplus1 = probMplus1,
                      reverse = reverse)
      }
    }


    A <- aperm(A, c(2, 1, 3, 4))  # c(ppp, nnn, M+1, M)
    for (jlocal in 1:(M + 1)) {
      for (tlocal in 1:M) {
        ansarray[,, jlocal]  <- ansarray[,, jlocal] +
                                A[,, jlocal, tlocal] * Bmat[, tlocal]
      }
    }
    ans.acat <- aperm(ansarray, c(1, 3, 2))  # c(ppp, M+1, nnn)
    dimnames(ans.acat) <- list(rownames(Bmat),
                               colnames(fitmat),
                               rownames(etamat))
    subsetarray3(ans.acat, subset = subset)
  })





  cratio.derivs <- function(jay, tee,
                            hdot, M, cpThetamat, Thetamat,
                            reverse = FALSE) {

    if (jay >= M+1) stop("argument 'jay' out of range")
    if (M   <  tee) stop("argument 'tee' out of range")

    if (reverse) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      if (jay == 1) {
        return(hdot[, tee] * cpThetamat[, 1] / Thetamat[, tee])
      }

      if (jay-1 == tee) {
        return(-hdot[, jay-1] * cpThetamat[, jay])
      }
      if (jay <= tee) {
        return((1 - Thetamat[, jay-1]) *
                hdot[, tee] * cpThetamat[, jay] / Thetamat[, tee])
      }
      return(rep_len(0, nrow(Thetamat)))  # Since jay-1 > tee
    } else {  # reverse = FALSE ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      if (jay == 1 && tee == 1) {
        return(-hdot[, 1])
      }

      if (jay == tee) {
        return(-hdot[, jay] * cpThetamat[, jay-1])
      }
      if (tee < jay) {
        return((1 - Thetamat[, jay]) *
                hdot[, tee] * cpThetamat[, jay-1] / Thetamat[, tee])
      }
      return(rep_len(0, nrow(Thetamat)))  # Since jay < tee
    } # reverse ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  }  # cratio.derivs






setMethod("margeffS4VGAM",  signature(VGAMff = "cratio"),
  function(object,
           subset = NULL,
           VGAMff,
           ...) {


  object <- callNextMethod(VGAMff = VGAMff,
                           object = object,
                           subset = subset,
                           ...)
  reverse <- object@post$reverse
  linkfunctions <- object@post$linkfunctions
  all.eargs <- object@post$all.eargs
  Bmat <- cfit <- object@post$Bmat
  ppp <- object@post$ppp
  etamat <- predict(object)  # nnn x M
  fitmat <- fitted(object)   # nnn x (M + 1)
  nnn <- nrow(etamat)
  M   <- ncol(etamat)
  hdot <- object@post$hdot
  Thetamat <- object@post$Thetamat







  vfamily <- object@family@vfamily
  c.nots <- any(vfamily == "cratio")

  if (any(vfamily == "cratio")) {
    cpThetamat <- if (reverse)
      tapplymat1(    Thetamat[, M:1, drop = FALSE],
                 "cumprod")[, M:1, drop = FALSE] else
      tapplymat1(    Thetamat, "cumprod")
  }



  A        <- array(0, c(i = nnn, vars = ppp, probs = M + 1, etas = M))
  ansarray <- array(0, c(vars = ppp, i = nnn, probs = M + 1))


  choosemat <- if (c.nots) Thetamat else 1 - Thetamat
  if (min(choosemat) <= 0)
    warning("division by 0 may occur")






    if (reverse) {
      for (tlocal in 1:M) {
        for (jlocal in 1:tlocal) {
          A[, , jlocal, tlocal] <-
            cratio.derivs(jay = jlocal, tee = tlocal,
                          hdot = ifelse(c.nots, 1, -1) * hdot,
                          M = M, cpThetamat = cpThetamat,
                          Thetamat = choosemat,
                          reverse = reverse)
        }
      }
      if (M > 1)
        for (jlocal in 2:M) {
          A[, , jlocal, jlocal-1] <-
            cratio.derivs(jay = jlocal, tee = jlocal-1,
                          hdot = ifelse(c.nots, 1, -1) * hdot,
                          M = M, cpThetamat = cpThetamat,
                          Thetamat = choosemat,
                          reverse = reverse)
        }
    } else {
     for (jlocal in 1:M) {
        for (tlocal in 1:jlocal) {
          A[, , jlocal, tlocal] <-
            cratio.derivs(jay = jlocal, tee = tlocal,
                          hdot = ifelse(c.nots, 1, -1) * hdot,
                          M = M, cpThetamat = cpThetamat,
                          Thetamat = choosemat,
                          reverse = reverse)
        }
      }
    }

    if (reverse) {
      A[, , M+1, M] <- ifelse(c.nots, -1, 1) * hdot[, M]
    } else {
      for (jlocal in 1:M) {
        for (tlocal in 1:jlocal) {
          A[, , M+1, tlocal] <- if (c.nots) {
            A[, , M+1, tlocal] - A[, , jlocal, tlocal]
          } else {
            -hdot[, tlocal] * cpThetamat[, M] / choosemat[, tlocal]
          }
        }
      }
    }

    A <- aperm(A, c(2, 1, 3, 4))  # c(ppp, nnn, M+1, M)
    for (jlocal in 1:(M + 1)) {
      for (tlocal in 1:M) {
        ansarray[,, jlocal]  <- ansarray[,, jlocal] +
                                A[,, jlocal, tlocal] * Bmat[, tlocal]
      }
    }
    ans.csratio <- aperm(ansarray, c(1, 3, 2))  # c(ppp, M+1, nnn)
    dimnames(ans.csratio) <- list(rownames(Bmat),
                                  colnames(fitmat),
                                  rownames(etamat))
    subsetarray3(ans.csratio, subset = subset)  # "cratio" & "sratio"
  })




setMethod("margeffS4VGAM",  signature(VGAMff = "sratio"),
  function(object,
           subset = NULL,
           VGAMff,
           ...) {








  object <- callNextMethod(VGAMff = VGAMff,
                           object = object,
                           subset = subset,
                           ...)
  reverse <- object@post$reverse
  linkfunctions <- object@post$linkfunctions
  all.eargs <- object@post$all.eargs
  Bmat <- cfit <- object@post$Bmat
  ppp <- object@post$ppp
  etamat <- predict(object)  # nnn x M
  fitmat <- fitted(object)   # nnn x (M + 1)
  nnn <- nrow(etamat)
  M   <- ncol(etamat)
  hdot <- object@post$hdot
  Thetamat <- object@post$Thetamat




  vfamily <- object@family@vfamily
  c.nots <- any(vfamily == "cratio")
  if (any(vfamily == "sratio")) {
    cpThetamat <- if (reverse)
      tapplymat1(1 - Thetamat[, M:1, drop = FALSE],
                 "cumprod")[, M:1, drop = FALSE] else
      tapplymat1(1 - Thetamat, "cumprod")
  }



  A        <- array(0, c(i = nnn, vars = ppp, probs = M + 1, etas = M))
  ansarray <- array(0, c(vars = ppp, i = nnn, probs = M + 1))


  choosemat <- if (c.nots) Thetamat else 1 - Thetamat
  if (min(choosemat) <= 0)
    warning("division by 0 may occur")











    if (reverse) {
      for (tlocal in 1:M) {
        for (jlocal in 1:tlocal) {
          A[, , jlocal, tlocal] <-
            cratio.derivs(jay = jlocal, tee = tlocal,
                          hdot = ifelse(c.nots, 1, -1) * hdot,
                          M = M, cpThetamat = cpThetamat,
                          Thetamat = choosemat,
                          reverse = reverse)
        }
      }
      if (M > 1)
        for (jlocal in 2:M) {
          A[, , jlocal, jlocal-1] <-
            cratio.derivs(jay = jlocal, tee = jlocal-1,
                          hdot = ifelse(c.nots, 1, -1) * hdot,
                          M = M, cpThetamat = cpThetamat,
                          Thetamat = choosemat,
                          reverse = reverse)
        }
    } else {
     for (jlocal in 1:M) {
        for (tlocal in 1:jlocal) {
          A[, , jlocal, tlocal] <-
            cratio.derivs(jay = jlocal, tee = tlocal,
                          hdot = ifelse(c.nots, 1, -1) * hdot,
                          M = M, cpThetamat = cpThetamat,
                          Thetamat = choosemat,
                          reverse = reverse)
        }
      }
    }

    if (reverse) {
      A[, , M+1, M] <- ifelse(c.nots, -1, 1) * hdot[, M]
    } else {
      for (jlocal in 1:M) {
        for (tlocal in 1:jlocal) {
          A[, , M+1, tlocal] <- if (c.nots) {
            A[, , M+1, tlocal] - A[, , jlocal, tlocal]
          } else {
            -hdot[, tlocal] * cpThetamat[, M] / choosemat[, tlocal]
          }
        }
      }
    }

    A <- aperm(A, c(2, 1, 3, 4))  # c(ppp, nnn, M+1, M)
    for (jlocal in 1:(M + 1)) {
      for (tlocal in 1:M) {
        ansarray[,, jlocal]  <- ansarray[,, jlocal] +
                                A[,, jlocal, tlocal] * Bmat[, tlocal]
      }
    }
    ans.csratio <- aperm(ansarray, c(1, 3, 2))  # c(ppp, M+1, nnn)
    dimnames(ans.csratio) <- list(rownames(Bmat),
                                  colnames(fitmat),
                                  rownames(etamat))
    subsetarray3(ans.csratio, subset = subset)  # "cratio" & "sratio"
  })








 margefff <- function(object, subset = NULL) {


  ii <- subset
  if (!is(object, "vglm"))
    stop("'object' is not a vglm() object")
  if (!any(temp.logical <-
    is.element(c("multinomial", "cumulative", "acat",
                 "cratio", "sratio"),
               object@family@vfamily)))
    stop("'object' is not a 'multinomial' or 'acat' ",
         "or 'cumulative'  or 'cratio' or 'sratio' VGLM!")
  vfamily <- object@family@vfamily
  if (is(object, "vgam"))
    stop("'object' is a vgam() object")
  if (length(object@control$xij))
    stop("'object' contains 'xij' terms")
  if (length(object@misc$form2))
    stop("'object' contains 'form2' terms")

  oassign <- object@misc$orig.assign
  if (any(unlist(lapply(oassign, length)) > 1))
    warning("some terms in 'object' create more than one column of ",
            "the LM design matrix")

  nnn <- object@misc$n
  M <- object@misc$M  # ncol(B)  # length(pvec) - 1


    if (any(vfamily == "multinomial")) {
    rlev <- object@misc$refLevel
    cfit <- coefvlm(object, matrix.out = TRUE)
    B <- if (!length(rlev)) {
      cbind(cfit, 0)
    } else {
      if (rlev == M+1) {  # Default
        cbind(cfit, 0)
      } else if (rlev == 1) {
        cbind(0, cfit)
      } else {
        cbind(cfit[, 1:(rlev-1)], 0, cfit[, rlev:M])
      }
    }
    ppp   <- nrow(B)
    pvec1 <- fitted(object)[1, ]
    colnames(B) <- if (length(names(pvec1))) names(pvec1) else
                   param.names("mu", M+1)

    if (is.null(ii)) {
      BB <- array(B, c(ppp, M+1, nnn))
      pvec  <- c(t(fitted(object)))
      pvec  <- rep(pvec, each = ppp)
      temp1 <- array(BB * pvec, c(ppp, M+1, nnn))
      temp2 <- aperm(temp1, c(2, 1, 3))  # (M+1) x ppp x nnn
      temp2 <- colSums(temp2)  # ppp x nnn
      temp2 <- array(rep(temp2, each = M+1), c(M+1, ppp, nnn))
      temp2 <- aperm(temp2, c(2, 1, 3))  # ppp x (M+1) x nnn
      temp3 <- pvec
      ans <- array((BB - temp2) * temp3, c(ppp, M+1, nnn),
                   dimnames = list(dimnames(B)[[1]],
                                   dimnames(B)[[2]],
                                   dimnames(fitted(object))[[1]]))
      return(ans)
    } else if (is.numeric(ii) && length(ii) == 1) {
        pvec  <- fitted(object)[ii, ]
        temp1 <- B * matrix(pvec, ppp, M+1, byrow = TRUE)
        temp2 <- matrix(rowSums(temp1), ppp, M+1)
        temp3 <- matrix(pvec, nrow(B), M+1, byrow = TRUE)
        return((B - temp2) * temp3)
    } else {
        if (is.logical(ii))
          ii <- (1:nnn)[ii]

        ans <- array(0, c(ppp, M+1, length(ii)),
                   dimnames = list(dimnames(B)[[1]],
                                   dimnames(B)[[2]],
                                   dimnames(fitted(object)[ii, ])[[1]]))
        for (ilocal in seq_along(ii)) {
          pvec  <- fitted(object)[ii[ilocal], ]
          temp1 <- B * matrix(pvec, ppp, M+1, byrow = TRUE)
          temp2 <- matrix(rowSums(temp1), ppp, M+1)
          temp3 <- matrix(pvec, nrow(B), M+1, byrow = TRUE)
          ans[ , , ilocal] <- (B - temp2) * temp3
        }
        return(ans)
      }
  }  # "multinomial"





  reverse <- object@misc$reverse
  linkfunctions <- object@misc$link
  all.eargs <- object@misc$earg
  B <- cfit <- coefvlm(object, matrix.out = TRUE)
  ppp <- nrow(B)
  etamat <- predict(object)  # nnn x M
  fitmat <- fitted(object)   # nnn x (M + 1)
  nnn <- nrow(etamat)


  hdot <- Thetamat <- etamat
  for (jlocal in 1:M) {
    Thetamat[, jlocal] <- eta2theta(etamat[, jlocal],
                                    link = linkfunctions[jlocal],
                                    earg = all.eargs[[jlocal]])
    hdot[, jlocal] <- dtheta.deta(Thetamat[, jlocal],
                                  link = linkfunctions[jlocal],
                                  earg = all.eargs[[jlocal]])
  }  # jlocal



  if (any(vfamily == "acat")) {
    expcs.etamat <- if (reverse)
      exp(tapplymat1(etamat[, M:1, drop = FALSE],
                     "cumsum")[, M:1, drop = FALSE]) else
      exp(tapplymat1(etamat, "cumsum"))
    csexpcs.etavec <- rowSums(expcs.etamat)
  }
  if (any(vfamily == "cratio")) {
    cpThetamat <- if (reverse)
      tapplymat1(    Thetamat[, M:1, drop = FALSE],
                 "cumprod")[, M:1, drop = FALSE] else
      tapplymat1(    Thetamat, "cumprod")
  }
  if (any(vfamily == "sratio")) {
    cpThetamat <- if (reverse)
      tapplymat1(1 - Thetamat[, M:1, drop = FALSE],
                 "cumprod")[, M:1, drop = FALSE] else
      tapplymat1(1 - Thetamat, "cumprod")
  }





     if (is.logical(is.multivariateY <-
         object@misc$multiple.responses) && is.multivariateY)
    stop("cannot handle cumulative(multiple.responses = TRUE)")




  if (any(vfamily == "cumulative")) {
    hdot.big <- kronecker(hdot, matrix(1, ppp, 1))  # Enlarged
    resmat <- cbind(hdot.big, 1)
    resmat[, 1] <- ifelse(reverse, -1, 1) * hdot.big[, 1] * cfit[, 1]

    if (M > 1) {
      for (jlocal in 2:M)
        resmat[, jlocal] <- ifelse(reverse, -1, 1) *
          (hdot.big[, jlocal    ] * cfit[, jlocal    ] -
           hdot.big[, jlocal - 1] * cfit[, jlocal - 1])

    }  # jlocal

    resmat[, M+1] <- ifelse(reverse, 1, -1) * hdot.big[, M] * cfit[, M]

    temp1 <- array(resmat, c(ppp, nnn, M+1),
                   dimnames = list(dimnames(B)[[1]],
                                   dimnames(fitted(object))[[1]],
                                   dimnames(fitted(object))[[2]]))
    temp1 <- aperm(temp1, c(1, 3, 2))  # ppp x (M+1) x nnn

    if (is.null(ii)) {
      return(temp1)
    } else
    if (is.numeric(ii) && (length(ii) == 1)) {
      return(temp1[, , ii])
    } else {
      return(temp1[, , ii])
    }
  }  # "cumulative"









  if (any(vfamily == "acat")) {
    if (!all(object@misc$link == "loglink"))
      stop("currently only the 'loglink' link is supported")


  acat.derivs <- function(jay, tee,
                          M, expcs.etamat, Thetamat,
                          prob1, probMplus1,
                          reverse = FALSE) {

    if (jay > M+1) stop("argument 'jay' out of range")
    if (M   < tee) stop("argument 'tee' out of range")

    if (reverse) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      dpMplus1.detat <- -(probMplus1^2) *
                        rowSums(expcs.etamat[, 1:tee, drop = FALSE])
      if (jay == M+1) {
        return(dpMplus1.detat)
      }
      if (jay <= tee) {
        return((probMplus1 + dpMplus1.detat) * expcs.etamat[, jay])
      }
      if (tee < jay) {
        return(dpMplus1.detat * expcs.etamat[, jay])
      }
    } else {  # reverse = FALSE ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      dp1.detat <-
        -(prob1^2) * rowSums(expcs.etamat[, tee:M, drop = FALSE])
      if (jay == 1) {
        return(dp1.detat)
      }
      if (jay <= tee) {
        return(dp1.detat * expcs.etamat[, jay-1])
      }
      if (tee < jay) {
        return((prob1 + dp1.detat) * expcs.etamat[, jay-1])
      }
    } # reverse ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  }  # acat.derivs



  A        <- array(0, c(i = nnn, vars = ppp, probs = M + 1, etas = M))
  ansarray <- array(0, c(vars = ppp, i = nnn, probs = M + 1))
  if (reverse) {
    probMplus1 <- 1 / (1 +  csexpcs.etavec) # Last level of Y
  } else {
    prob1 <- 1 / (1 +  csexpcs.etavec)  # First level of Y
  }

    for (jlocal in 1:(M+1)) {
      for (tlocal in 1:M) {
        A[, , jlocal, tlocal] <-
          acat.derivs(jay = jlocal, tee = tlocal,
                      M = M, expcs.etamat = expcs.etamat,
                      Thetamat = Thetamat,
                      prob1 = prob1, probMplus1 = probMplus1,
                      reverse = reverse)
      }
    }


    A <- aperm(A, c(2, 1, 3, 4))  # c(ppp, nnn, M+1, M)
    for (jlocal in 1:(M + 1)) {
      for (tlocal in 1:M) {
        ansarray[,, jlocal]  <- ansarray[,, jlocal] +
                                A[,, jlocal, tlocal] * B[, tlocal]
      }
    }
    ans.acat <- aperm(ansarray, c(1, 3, 2))  # c(ppp, M+1, nnn)
    dimnames(ans.acat) <- list(rownames(B),
                               colnames(fitmat),
                               rownames(etamat))
    return(ans.acat)
  }  # "acat"






  c.nots <- any(vfamily == "cratio")

  cratio.derivs <- function(jay, tee,
                            hdot, M, cpThetamat, Thetamat,
                            reverse = FALSE) {

    if (jay >= M+1) stop("argument 'jay' out of range")
    if (M   <  tee) stop("argument 'tee' out of range")

    if (reverse) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      if (jay == 1) {
        return(hdot[, tee] * cpThetamat[, 1] / Thetamat[, tee])
      }

      if (jay-1 == tee) {
        return(-hdot[, jay-1] * cpThetamat[, jay])
      }
      if (jay <= tee) {
        return((1 - Thetamat[, jay-1]) *
                hdot[, tee] * cpThetamat[, jay] / Thetamat[, tee])
      }
      return(rep_len(0, nrow(Thetamat)))  # Since jay-1 > tee
    } else {  # reverse = FALSE ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      if (jay == 1 && tee == 1) {
        return(-hdot[, 1])
      }

      if (jay == tee) {
        return(-hdot[, jay] * cpThetamat[, jay-1])
      }
      if (tee < jay) {
        return((1 - Thetamat[, jay]) *
                hdot[, tee] * cpThetamat[, jay-1] / Thetamat[, tee])
      }
      return(rep_len(0, nrow(Thetamat)))  # Since jay < tee
    } # reverse ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  }  # cratio.derivs


  A        <- array(0, c(i = nnn, vars = ppp, probs = M + 1, etas = M))
  ansarray <- array(0, c(vars = ppp, i = nnn, probs = M + 1))


  choosemat <- if (c.nots) Thetamat else 1 - Thetamat
  if (min(choosemat) <= 0)
    warning("division by 0 may occur")




  if (any(vfamily == "cratio" | vfamily == "sratio")) {


    if (reverse) {
      for (tlocal in 1:M) {
        for (jlocal in 1:tlocal) {
          A[, , jlocal, tlocal] <-
            cratio.derivs(jay = jlocal, tee = tlocal,
                          hdot = ifelse(c.nots, 1, -1) * hdot,
                          M = M, cpThetamat = cpThetamat,
                          Thetamat = choosemat,
                          reverse = reverse)
        }
      }
      if (M > 1)
        for (jlocal in 2:M) {
          A[, , jlocal, jlocal-1] <-
            cratio.derivs(jay = jlocal, tee = jlocal-1,
                          hdot = ifelse(c.nots, 1, -1) * hdot,
                          M = M, cpThetamat = cpThetamat,
                          Thetamat = choosemat,
                          reverse = reverse)
        }
    } else {
     for (jlocal in 1:M) {
        for (tlocal in 1:jlocal) {
          A[, , jlocal, tlocal] <-
            cratio.derivs(jay = jlocal, tee = tlocal,
                          hdot = ifelse(c.nots, 1, -1) * hdot,
                          M = M, cpThetamat = cpThetamat,
                          Thetamat = choosemat,
                          reverse = reverse)
        }
      }
    }

    if (reverse) {
      A[, , M+1, M] <- ifelse(c.nots, -1, 1) * hdot[, M]
    } else {
      for (jlocal in 1:M) {
        for (tlocal in 1:jlocal) {
          A[, , M+1, tlocal] <- if (c.nots) {
            A[, , M+1, tlocal] - A[, , jlocal, tlocal]
          } else {
            -hdot[, tlocal] * cpThetamat[, M] / choosemat[, tlocal]
          }
        }
      }
    }

    A <- aperm(A, c(2, 1, 3, 4))  # c(ppp, nnn, M+1, M)
    for (jlocal in 1:(M + 1)) {
      for (tlocal in 1:M) {
        ansarray[,, jlocal]  <- ansarray[,, jlocal] +
                                A[,, jlocal, tlocal] * B[, tlocal]
      }
    }
    ans.csratio <- aperm(ansarray, c(1, 3, 2))  # c(ppp, M+1, nnn)
    dimnames(ans.csratio) <- list(rownames(B),
                                  colnames(fitmat),
                                  rownames(etamat))
    return(ans.csratio)
  }  # "cratio" and "sratio"


}  # margefff








prplot <- function(object,
                  control = prplot.control(...), ...) {


  if (!any(slotNames(object) == "family") ||
      !any(object@family@vfamily == "VGAMcategorical"))
    stop("'object' does not seem to be a VGAM categorical ",
         "model object")

  if (!any(object@family@vfamily == "cumulative"))
    stop("'object' does not seem to be a VGAM categorical ",
         "model object")

    control <- prplot.control(...)


  object <- plot.vgam(object, plot.arg = FALSE, raw = FALSE)  # , ...

  if (length(names(object@preplot)) != 1)
      stop("object needs to have only one term")


  MM <- object@misc$M
  use.y <- cbind((object@preplot[[1]])$y)
  Constant <- attr(object@preplot, "Constant")
  if (is.numeric(Constant) && length(Constant) == ncol(use.y))
    use.y <- use.y + matrix(Constant, nrow(use.y), ncol(use.y),
                            byrow = TRUE)
  for (ii in 1:MM) {
    use.y[, ii] <- eta2theta(use.y[, ii],
                             link = object@misc$link[[ii]],
                             earg = object@misc$earg[[ii]])
  }
  if (ncol(use.y) != MM) use.y = use.y[, 1:MM, drop = FALSE]

  use.x <- (object@preplot[[1]])$x
  myxlab <- if (length(control$xlab))
           control$xlab else (object@preplot[[1]])$xlab
  mymain <- if (MM <= 3)
           paste(object@misc$parameters, collapse = ", ") else
           paste(object@misc$parameters[c(1, MM)],
                 collapse = ",...,")
  if (length(control$main)) mymain = control$main
  if (length(control$ylab)) myylab = control$ylab

  matplot(use.x, use.y, type = "l",
          xlab = myxlab, ylab = myylab,
          lty = control$lty, col = control$col, las = control$las,
          xlim = if (is.Numeric(control$xlim))
                 control$xlim else range(use.x),
          ylim = if (is.Numeric(control$ylim))
                 control$ylim else range(use.y),
          main=mymain)
  if (control$rug.arg)
    rug(use.x, col=control$rcol, lwd=control$rlwd)

  invisible(object)
}




 prplot.control <- function(xlab = NULL, ylab = "Probability",
                           main = NULL,
                           xlim = NULL, ylim = NULL,
                           lty = par()$lty,
                           col = par()$col,
                           rcol = par()$col,
                           lwd = par()$lwd,
                           rlwd = par()$lwd,
                           las = par()$las,
                           rug.arg  = FALSE,
                           ...) {

    list(xlab = xlab, ylab = ylab,
         xlim = xlim, ylim = ylim,
         lty = lty, col = col, rcol = rcol,
         lwd = lwd, rlwd = rlwd, rug.arg = rug.arg,
         las = las, main = main)
}










is.parallel.matrix <- function(object, ...)
  is.matrix(object) && all(!is.na(object)) &&
  all(c(object) == 1) && ncol(object) == 1


is.parallel.vglm <- function(object, type = c("term", "lm"), ...) {

  type <- match.arg(type, c("term", "lm"))[1]
  Hlist <- constraints(object, type = type)

  unlist(lapply(Hlist, is.parallel.matrix))
}


if (!isGeneric("is.parallel"))
  setGeneric("is.parallel", function(object, ...)
             standardGeneric("is.parallel"),
             package = "VGAM")


setMethod("is.parallel",  "matrix", function(object, ...)
    is.parallel.matrix(object, ...))


setMethod("is.parallel",  "vglm", function(object, ...)
    is.parallel.vglm(object, ...))






is.zero.matrix <- function(object, ...) {

  rnames <- rownames(object)
  intercept.index <- if (length(rnames)) {
    if (any(rnames == "(Intercept)")) {
      (seq_along(rnames))[rnames == "(Intercept)"]
    } else {
      stop("the matrix does not seem to have an intercept")
      NULL
    }
  } else {
      stop("the matrix does not seem to have an intercept")
      NULL
  }

  if (nrow(object) <= 1)
    stop("the matrix needs to have more than one row, ",
         "i.e., more than ",
         "an intercept on the RHS of the formula")

  cfit <- object[-intercept.index, , drop = FALSE]

  foo <- function(conmat.col)
    all(!is.na(conmat.col)) &&
    all(c(conmat.col) == 0)

  unlist(apply(cfit, 2, foo))
}


is.zero.vglm <- function(object, ...) {
  is.zero.matrix(coef(object, matrix = TRUE))
}


if (!isGeneric("is.zero"))
  setGeneric("is.zero", function(object, ...)
             standardGeneric("is.zero"),
             package = "VGAM")


setMethod("is.zero",  "matrix", function(object, ...)
          is.zero.matrix(object, ...))


setMethod("is.zero",  "vglm", function(object, ...)
          is.zero.vglm(object, ...))






setMethod("showvglmS4VGAM",
          signature(VGAMff = "acat"),
  function(object,
           VGAMff,
           ...) {
  cat("\nThis is an adjacent categories model with",
      1 + object@misc$M, "levels\n")
  invisible(object)
  })


setMethod("showvgamS4VGAM",
          signature(VGAMff = "acat"),
  function(object,
           VGAMff,
           ...) {
  cat("\nThis is an adjacent categories model with",
      1 + object@misc$M, "levels\n")
  invisible(object)
  })



setMethod("showvglmS4VGAM",
          signature(VGAMff = "multinomial"),
  function(object,
           VGAMff,
           ...) {
  cat("\nThis is a multinomial logit model with",
      1 + object@misc$M, "levels\n")
  invisible(object)
  })


setMethod("showvgamS4VGAM",
          signature(VGAMff = "multinomial"),
  function(object,
           VGAMff,
           ...) {
  cat("\nThis is a multinomial logit model with",
      1 + object@misc$M, "levels\n")
  invisible(object)
  })






R2latvar <- function(object) {


  if (!is(object, "vglm"))
      stop("argument 'object' is not a vglm() fit")
 

  vfam <- object@family@vfamily
  fam.permitted <- c("cumulative", "propodds", "binomialff")
  if (!(vfam[1] %in% fam.permitted))
    stop("the family function must be one of ", fam.permitted)


  linkfn <- linkfun(object)[1]
  link.permitted <- c("logitlink", "probitlink",
                      "clogloglink")  # "cauchitlink"
  if (!(linkfn %in% link.permitted))
    stop("allowable link functions supported are ", link.permitted)

  infos <- object@family@infos()
  if (length(infos$parallel) == 1 &&
      is.logical(infos$parallel))
    if (!infos$parallel)
      stop("the linear predictors are not parallel")
  if (!all(unlist(constraints(object)[-1]) == 1))
      stop("the linear predictors are not parallel")



  eta1 <- predict(object)[, 1]
  offset <- switch(linkfn,
                   logitlink   = (pi^2)/3,
                   probitlink  = 1,
                   clogloglink = (pi^2) / 6,
                   stop("link unrecognized"))

  veta1 <- var(eta1)
  veta1 / (veta1 + offset) 
}








ordsup.vglm <-
  function(object, all.vars = FALSE, confint = FALSE, ...) {


  if (!is(object, "vglm"))
    stop("argument 'object' is not a vglm() fit")
 

  vfam <- object@family@vfamily
  fam.permitted <- c("uninormal", "propodds", "cumulative")
  if (!(vfam[1] %in% fam.permitted))
    stop("the family function must be one of ",
         paste(fam.permitted, collapse = ", "))


  linkfns <- linkfun(object)
  link.permitted <- c("identitylink", "logitlink", "probitlink")
  if (!any(linkfns %in% link.permitted))
    stop("allowable link functions supported are ",
         paste("'", link.permitted, "'",
               sep = "", collapse = ", "))


  R <- npred(object) / npred(object, type = "one.response")      
  if (R > 1)
    stop("currently cannot handle multiple responses")


  infos <- object@family@infos()
  cobj <- coef(object)
  cobj.mat <- coef(object, matrix = TRUE)
  x.LM  <- model.matrix(object, type =  "lm")
  is.binary <- apply(if (has.intercept(object))
                     x.LM[, -1, drop = FALSE] else x.LM,
                     2, function(coln) length(unique(coln)) == 2)

  if (all.vars) {
    is.binary[] <- TRUE
  } else {
    if (sum(is.binary) == 0) {
      warning("no binary explanatory variables; returning a NULL")
      return(NULL)
    }
  }


  is.binary <- names(is.binary)[is.binary]


      

  switch(vfam[1],
  uninormal = {
    if (!all(linkfns[c(TRUE, FALSE)] == "identitylink"))
      stop("must use 'identitylink' for the mean parameter")               
    if (!all(cobj.mat[-1, c(FALSE, TRUE)] == 0))
      stop("the sdev or var parameter must be intercept-only")
    sdev <- numeric(ncol(cobj.mat) / 2)
    for (jay in 1:(ncol(cobj.mat) / 2))
      sdev[jay] <- eta2theta(cobj.mat["(Intercept)", 2*jay],
                             link = object@misc$link[2*jay],
                             earg = object@misc$earg[2*jay])
    if (infos$parameters.names[2] == "var")
      sdev <- sqrt(sdev)  # It was actually the variance before
  },
  cumulative = {

  if (length(infos$parallel) == 1 &&
      is.logical(infos$parallel))
    if (!infos$parallel)
      stop("the linear predictors are not parallel")
  if (!all(unlist(constraints(object)[-1]) == 1))
      stop("the linear predictors are not parallel")


  reverse <- infos$reverse
  if (FALSE)
  if (!reverse) {
    stop("the 'reverse' argument must be TRUE")
  }
  },
  {
  stop("family unrecognized")
  })



  gamma.fun <- function(cobj, is.binary, vfam,
                        reverse, linkfns = NULL) {
  switch(vfam[1],
  uninormal = {
    gamma <- pnorm(cobj[ is.binary ] / (sqrt(2) * sdev))
  },
  cumulative = {
    gamma <-
    switch(linkfns[1],
    clogloglink = clogloglink(ifelse(reverse, 1, -1) *  # Not sure
                              cobj[ is.binary ], inverse = TRUE),
    logitlink   =   logitlink(ifelse(reverse, 1, -1) *
                              cobj[ is.binary ] / sqrt(2),
                              inverse = TRUE),
    probitlink  =  probitlink(ifelse(reverse, 1, -1) *
                              cobj[ is.binary ] / sqrt(2),
                              inverse = TRUE))
  },
  zzzz = {
  })
  gamma
  }


      

  gamma <- gamma.fun(cobj, is.binary, vfam, reverse, linkfns)
  Delta <- 2 * gamma - 1

      

  

  if (confint) {
    ans2 <- confint(object, parm = is.binary, ...)
    if (!is.matrix(ans2))
      ans2 <- matrix(ans2, 1, 2,
                     dimnames = list(is.binary, names(ans2)))
    ans2.low <- gamma.fun(ans2[, 1], TRUE, vfam, reverse, linkfns)
    ans2.upp <- gamma.fun(ans2[, 2], TRUE, vfam, reverse, linkfns)
    Delta.low <- 2 * ans2.low - 1
    Delta.upp <- 2 * ans2.upp - 1
    names(ans2.low) <- names(ans2.upp) <-
    names(Delta.low) <- names(Delta.upp) <- is.binary
  }


  c(
  list(gamma = gamma,
       Delta = Delta),
       if (confint) {
         list(lower.gamma = ans2.low,
              upper.gamma = ans2.upp,
              Lower.Delta = Delta.low,
              Upper.Delta = Delta.upp)
       } else NULL)
}  # ordsup.vglm



 if (!isGeneric("ordsup"))
     setGeneric("ordsup",
  function(object, ...) standardGeneric("ordsup"))

setMethod("ordsup",  "vglm",
  function(object, ...) ordsup.vglm(object, ...))















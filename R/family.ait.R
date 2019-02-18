# These functions are
# Copyright (C) 1998-2019 T.W. Yee, University of Auckland.
# All rights reserved.


















dgentpois <-
  function(x, lambda, truncate = 0,  # NULL,
           max.support = Inf,  # 20181018
           log = FALSE) {
  if (!is.numeric(max.support) ||
      length(max.support) != 1 ||
      round(max.support) != max.support ||
      max.support <= max(truncate))
    stop("bad input for argument 'max.support'")

  if (is.null(truncate) && is.infinite(max.support) && 0 < max.support)
    return(dpois(x, lambda, log = log))
    
  if (is.list(truncate))
    truncate <- unlist(truncate)
    
  if (!is.Numeric(truncate, integer.valued = TRUE) ||
      any(truncate < 0))
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")

  LLL <- max(length(x), length(lambda))
  if (length(x)      != LLL) x      <- rep(x,      length = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, length = LLL)
  sump <- rep(0, length = LLL)
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
}  # dgentpois



pgentpois <-
  function(q, lambda, truncate = 0,  # NULL,
           max.support = Inf) {
  if (is.null(truncate))
    return(ppois(q, lambda))
    
  if (is.list(truncate))
    truncate <- unlist(truncate)
    
  if (!is.Numeric(truncate, integer.valued = TRUE) ||
      any(truncate < 0))
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")
  if (!is.numeric(max.support) ||
      length(max.support) != 1 ||
      round(max.support) != max.support ||
      max.support <= max(truncate))
    stop("bad input for argument 'max.support'")

  LLL <- max(length(q), length(lambda))
  if (length(q)      != LLL) q      <- rep(q,      length = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, length = LLL)
  sump <- rep(0, length = LLL)
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
}  # pgentpois




rgentpois <-
  function(n, lambda, truncate = 0,  # NULL
           max.support = Inf,  # 20181124
           maxits = 10000) {
  if (!is.numeric(max.support) ||
      length(max.support) != 1 ||
      round(max.support) != max.support ||
      max.support <= max(truncate))
    stop("bad input for argument 'max.support'")

  use.n <- if ((length.n <- length(n)) > 1) 
    length.n else
      if (!is.Numeric(n, integer.valued = TRUE, 
                      length.arg = 1, positive = TRUE))
        stop("bad input for argument 'n'") else n
    
  if (is.null(truncate) && is.infinite(max.support) && 0 < max.support)
    return(rpois(use.n, lambda))
    
  if (is.list(truncate))
    truncate <- unlist(truncate)
    
  if (!is.Numeric(truncate, integer.valued = TRUE) ||
      any(truncate < 0))
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")

  lambda <- rep(lambda, length = use.n)
  ans <- rpois(use.n, lambda)
  ind.replace <- (ans %in% truncate) | (max.support < ans)

  iter <- 0
  while (any(ind.replace)) {
    ans[ind.replace] <- rpois(sum(ind.replace), lambda[ind.replace])
    ind.replace <- (ans %in% truncate) | (max.support < ans)
    iter <- iter + 1
    if (iter > maxits) {
      warning("reached 'maxits' iterations; breaking and returning NAs")
      ans[ind.replace] <- NA
      break
    }
  }  # while
    
  ans
}  # rgentpois





 gentpoisson <-
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
      max.support <= max(truncate))
    stop("bad input for argument 'max.support'")

  if (length( ilambda) && !is.Numeric(ilambda, positive = TRUE))
    stop("bad input for argument 'ilambda'")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "prob.t"))[1]




  new("vglmff",
  blurb = c("Generally truncated-Poisson distribution\n\n",
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
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "prob.t"))[1]

    lambda <- eta2theta(eta, .link , earg = .earg )
    sump <- matrix(0, NROW(eta), NCOL(eta))
    Sume <- matrix(0, NROW(eta), NCOL(eta))
    for (tval in .truncate ) {
      pmf <- dpois(tval, lambda)
      sump <- sump + pmf
      Sume <- Sume + pmf * tval
    }
    ans <- switch(type.fitted,
       "mean"      = (lambda - Sume) / (1 - sump),
       "lambda"    = lambda,
       "prob.t"    = sump)  # Pr(Y=truncatedvalue) as it were
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
      ll.elts <- c(w) * dgentpois(y, lambda = lambda, log = TRUE,
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
  vfamily = c("gentpoisson"),
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
    rgentpois(nsim * length(lambda), lambda, truncate = .truncate ,
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

    sumderiv1 <- matrix(0, NROW(eta), NCOL(eta))
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
}  # gentpoisson







dgentbinom <-
  function(x, size, prob, truncate = 0,  # NULL,
           log = FALSE) {
  if (is.null(truncate))
    return(dbinom(x, size, prob, log = log))
    
  LLL <- max(length(x), length(size), length(prob))
  if (length(x)      != LLL) x      <- rep(x,      length = LLL)
  if (length(size  ) != LLL) size   <- rep(size,   length = LLL)
  if (length(prob  ) != LLL) prob   <- rep(prob,   length = LLL)

  if (is.list(truncate))
    truncate <- unlist(truncate)
    
  if (!is.Numeric(truncate, integer.valued = TRUE) ||
      any(truncate < 0))  # || max(truncate) > size
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")
    
  sump <- rep(0, length = LLL)
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
}  # dgentbinom



pgentbinom <-
  function(q, size, prob, truncate = 0) {
  if (is.null(truncate))
    return(pbinom(q, size, prob))
    
  if (is.list(truncate))
    truncate <- unlist(truncate)
    
  if (!is.Numeric(truncate, integer.valued = TRUE) ||
      any(truncate < 0) ||
      any(size < truncate))
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")

  LLL <- max(length(q), length(size), length(prob))
  if (length(q)      != LLL) q      <- rep(q,      length = LLL)
  if (length(size  ) != LLL) size   <- rep(size,   length = LLL)
  if (length(prob  ) != LLL) prob   <- rep(prob,   length = LLL)

  sump <- rep(0, length = LLL)
  for (tval in truncate)
    sump <- sump + dbinom(tval, size, prob)
  denom <- 1 - sump

  numer <- pbinom(q, size, prob)
  for (tval in truncate) {
    if (any(vecTF <- tval <= q))
      numer[vecTF] <- numer[vecTF] - dbinom(tval, size[vecTF], prob[vecTF])
  }
  ans <- numer / denom
  ans
}  # pgentbinom




rgentbinom <-
  function(n, size, prob, truncate = 0,  # NULL
           maxits = 10000) {
  use.n <- if ((length.n <- length(n)) > 1) 
    length.n else
      if (!is.Numeric(n, integer.valued = TRUE, 
                      length.arg = 1, positive = TRUE))
        stop("bad input for argument 'n'") else n
    
  if (is.null(truncate))
    return(rbinom(use.n, size, prob))
    
  if (is.list(truncate))
    truncate <- unlist(truncate)
    
  if (!is.Numeric(truncate, integer.valued = TRUE) ||
      any(truncate < 0))  # || any(truncate > size))
    stop("bad input for argument 'truncate'")
  if (!identical(truncate, (unique(truncate))))
    stop("values of argument 'truncate' must be unique")

  size <- rep(size, length = use.n)
  prob <- rep(prob, length = use.n)
  ans <- rbinom(use.n, size, prob)
  ind.replace <- ans %in% truncate

  iter <- 0
  while (any(ind.replace)) {
    ans[ind.replace] <- rbinom(sum(ind.replace), size[ind.replace],
                               prob[ind.replace])
    ind.replace <- ans %in% truncate
    iter <- iter + 1
    if (iter > maxits) {
      warning("reached 'maxits' iterations; breaking and returning NAs")
      ans[ind.replace] <- NA
      break
    }
  }  # while

  ans
}  # rgentbinom






 gentbinomial <-
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



  if (!is.logical(multiple.responses) || length(multiple.responses) != 1)
    stop("bad input for argument 'multiple.responses'")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "prob", "prob.t"))[1]


  new("vglmff",
  blurb = c("Generally truncated binomial distribution\n\n",
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

    mustart.orig <- mustart
    if ( .multiple.responses ) {
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


    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly

      extra$orig.w <- w
      mustart <- matrix(colSums(y) / colSums(w),  # Not colSums(y * w)...
                        n, ncoly, byrow = TRUE)

    } else {
      eval(binomialff(link = .earg ,  # earg = .earg ,
                      earg.link = TRUE)@initialize)
    }
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)


    if ( .multiple.responses ) {

      dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
      dn2 <- if (length(dn2)) {
        paste("E[", dn2, "]", sep = "")
      } else {
        param.names("prob", M)
      }
      predictors.names <-
        namesof(if (M > 1) dn2 else "prob",
                .link , earg = .earg, short = TRUE)

      w <- matrix(w, n, ncoly)
      y <- y / w  # Now sample proportion
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



    nvec <- if (NCOL(y) > 1) {
              NULL
            } else {
              if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)
            }
    extra$tau <- if (length(nvec) && length(unique(nvec) == 1))
                   nvec[1] else NULL
    if (any(round(y * nvec) %in% .truncate ))
      stop("some response values equal values of the 'truncate' ",
           "argument")
  }), list( .link = link,
            .truncate = truncate,
            .type.fitted = type.fitted,
            .earg = earg, .multiple.responses = multiple.responses ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    w <- extra$w
    binprob <- eta2theta(eta, .link , earg = .earg )
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }
    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "prob.t"))[1]
    nvec <- if ( .multiple.responses ) {
             w
           } else {
             if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
               round(w)
           }
    NOS <- NCOL(eta)
    sump <- matrix(0, NROW(eta), NCOL(eta))
    Sume <- matrix(0, NROW(eta), NCOL(eta))
    for (tval in .truncate ) {
      pmf <- dbinom(tval, size = nvec, binprob)
      sump <- sump + pmf
      Sume <- Sume + pmf * (tval / nvec)
    }
    ans <- switch(type.fitted,
       "mean"      = (binprob - Sume) / (1 - sump),
       "prob"      = binprob,
       "prob.t"    = sump)  # Pr(Y=truncatedvalue) as it were
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



    if (length(extra$tau)) {
      R <- tfit$qr$qr[1:ncol.X.vlm, 1:ncol.X.vlm, drop = FALSE]
      R[lower.tri(R)] <- 0
      tmp6 <- N.hat.posbernoulli(eta = eta, link = .link , earg = .earg ,
                             R = R, w = w,
                             X.vlm = X.vlm.save,
                             Hlist = Hlist,  # 20150428; bug fixed here
                             extra = extra, model.type = "0")
      extra$N.hat    <- tmp6$N.hat
      extra$SE.N.hat <- tmp6$SE.N.hat
    }


  }), list( .link = link, .earg = earg,
            .truncate = truncate,
            .multiple.responses = multiple.responses ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

      ycounts <- if ( .multiple.responses ) {
                  round(y * extra$orig.w)
                 } else {
                   if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                   y * w  # Convert proportions to counts
                 }
      nvec <- if ( .multiple.responses ) {
                w
              } else {
                if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                  round(w)
              }
      use.orig.w <- if (is.numeric(extra$orig.w)) extra$orig.w else 1
    binprob <- eta2theta(eta, .link , earg = .earg )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      answer <- c(use.orig.w) * dgentbinom(ycounts, size = nvec,
                                           truncate = .truncate ,
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

  vfamily = c("gentbinomial"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    binprob <- eta2theta(eta, .link , earg = .earg )
    okay1 <- all(is.finite(binprob)) && all(0 < binprob & binprob < 1)
    okay1
  }, list( .link = link,
           .truncate = truncate,
           .earg = earg ))),





  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")

    if ( .multiple.responses )
      stop("cannot run simulate() when 'multiple.responses = TRUE'")

    eta <- predict(object)
    binprob <- eta2theta(eta, .link , earg = .earg )

    extra <- object@extra
    w <- extra$w  # Usual code
    w <- pwts  # 20140101


    nvec <- if ( .multiple.responses ) {
              w
            } else {
              if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                round(w)
            }
    rgentbinom(nsim * length(eta), size = nvec, prob = binprob,
               truncate = .truncate )
  }, list( .link = link, .earg = earg,
           .truncate = truncate,
           .multiple.responses = multiple.responses ))),



  deriv = eval(substitute(expression({
    use.orig.w <- if (is.numeric(extra$orig.w)) extra$orig.w else
                  rep_len(1, n)

    ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
               y * w # Convert proportions to counts
    nvec <- if ( .multiple.responses ) {
              w
            } else {
              if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)
            }
    size <- nvec
    binprob <- eta2theta(eta, .link , earg = .earg )
    dmu.deta <- dtheta.deta(binprob, .link , earg = .earg )

    sump <-
    Sume <- matrix(0, NROW(eta), NCOL(eta))
    for (tval in .truncate ) {
      pmf <- dbinom(tval, size, binprob)
      sump <- sump + pmf
      Sume <- Sume + pmf * (tval / size)
    }

    pmf.deriv1 <- function(ycount, size, prob)
      dbinom(ycount, size, prob) * size *
      (ycount / size - prob) / (prob * (1 - prob))
    pmf.deriv2 <- function(ycount, size, prob)
      dbinom(ycount, size, prob) * size * (
      size * ((ycount/size / prob) -
              (1 - ycount / size) / (1 - prob))^2 -
      ycount / size / prob^2 - (1 - ycount / size) / (1 - prob)^2)

    sumderiv1 <- matrix(0, NROW(eta), NCOL(eta))
    sumderiv2 <- matrix(0, NROW(eta), NCOL(eta))
    for (tval in .truncate ) {
      sumderiv1 <- sumderiv1 + pmf.deriv1(tval, size, binprob)
      sumderiv2 <- sumderiv2 + pmf.deriv2(tval, size, binprob)
    }

    dl.dmu <- size * y / binprob -
              size * (1 - y) / (1 - binprob) +
              sumderiv1 / (1 - sump)  # - (1 - binprob) * temp3 / temp1
    dl.dmu * dmu.deta
  }), list( .link = link, .earg = earg,
            .truncate = truncate,
            .multiple.responses = multiple.responses ))),
  weight = eval(substitute(expression({
    EY <- (binprob - Sume) / (1 - sump) # The expectation of Y
    ned2l.dmu2 <- size * EY / binprob^2 +
                  size * (1 - EY) / (1 - binprob)^2 -
       sumderiv2 / (1 - sump) -
      (sumderiv1 / (1 - sump))^2

    wz <- ned2l.dmu2 * dmu.deta^2  # c(w) already incorporated inside
    wz
  }), list( .link = link, .earg = earg,
            .truncate = truncate,
            .multiple.responses = multiple.responses ))))
}  # gentbinomial






dgenapois <-
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
  if (length(x)      != LLL) x      <- rep(x,      length = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, length = LLL)

  lalter <- length(alter)
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  suma <- rowSums(pobs.a)
  if (any(1 < suma))
    stop("bad input for argument 'pobs.a'")

  if (log.arg) {
    logpmf <- log1p(-suma) +
              dgentpois(x, lambda, truncate = alter, log = TRUE)
  } else {
    pmf <- exp(log1p(-suma)) * dgentpois(x, lambda, truncate = alter)
  }
 for (jay in seq_along(alter)) {
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
}  # dgenapois




pgenapois <-
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
  if (length(q)      != LLL) q      <- rep(q,      length = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, length = LLL)


  lalter <- length(alter)
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  suma <- rowSums(pobs.a)
  if (any(1 < suma))
    stop("bad input for argument 'pobs.a'")
  numer1 <- 1 - suma

  offset <- rep(0, length = LLL)
  for (jay in seq_along(alter)) {
    aval <- alter[jay]
    if (any(vecTF <- aval <= q)) {
      offset[vecTF] <- offset[vecTF] + pobs.a[vecTF, jay]
    }
  }
  ans <- numer1 * pgentpois(q, lambda, truncate = alter) + offset
  ans
}  # pgenapois





rgenapois <-
  function(n, lambda, alter = 0,  # NULL
           pobs.a = 0, byrow.arg = FALSE,
           maxits = 10000) {
  use.n <- if ((length.n <- length(n)) > 1) 
    length.n else
      if (!is.Numeric(n, integer.valued = TRUE, 
                      length.arg = 1, positive = TRUE))
        stop("bad input for argument 'n'") else n
    
  if (is.null(alter))
    return(rpois(use.n, lambda))
    
  if (is.list(alter))
    alter <- unlist(alter)
    
  if (!is.Numeric(alter, integer.valued = TRUE) ||
      any(alter < 0))
    stop("bad input for argument 'alter'")
  if (!identical(alter, (unique(alter))))
    stop("values of argument 'alter' must be unique")
    

  lalter <- length(alter)
  LLL <- use.n
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  suma <- rowSums(pobs.a)
  if (any(1 < suma))
    stop("bad input for argument 'pobs.a'")

  pobs.a <- cbind(pobs.a, 1 - suma)


  lambda <- rep(lambda, length = use.n)


  cpobs.a <- tapplymat1(pobs.a, "cumsum")
  choice <- rep(1, use.n)
  rrr <- runif(use.n)
  if (lalter > 1) {
    for (jay in 2:(lalter + 1)) {
       jayvec <- rep(jay, use.n)
       choice <- ifelse(cpobs.a[, jay-1] < rrr & rrr <= cpobs.a[, jay],
                        jay, choice)
    }
  }

  ans <- alter[choice]
  ind7 <- choice == (lalter + 1)
  if (any(ind7))
    ans[ind7] <- rgentpois(sum(ind7), lambda = lambda[ind7],
                           truncate = alter, maxits = maxits)
  ans
}  # rgenapois









 genapoisson <-
  function(alter = 0,
           zero = NULL,
           llambda = "loglink",
           type.fitted = c("mean", "lambda", "pobs.a", "onempobs.a"),
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
  }
 
  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "pobs.a", "onempobs.a"))[1]
  temp7 <- paste("pobs", alter, sep = "")


  new("vglmff",
  blurb = c("Generally-altered Poisson distribution\n",
            "(multinomial and truncated-Poisson conditional model)\n\n",
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
         link = c(if (length( .alter ))
                  rep("multilogitlink",  length( .alter )) else NULL,
                  .llambda ),
         link1parameter = if (length( .alter ))
                          FALSE else TRUE,  # mulilogitlink is multiparameter
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
           .alter = alter
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
    for (jay in seq_along(alter))
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

    temp7 <- paste("pobs", alter, sep = "")
    mynames1 <- paste("multinomial(", temp7, ")", sep = "")
    mynames2 <- param.names("lambda", ncoly, skip1 = TRUE)
    predictors.names <-
        c(        mynames1,
          namesof(mynames2, .llambda , earg = .elambda , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                             imu = .ilambda ,
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
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "pobs.a", "onempobs.a"))[1]

    alter <- ( .alter )
    M1 <- length(alter) + 1
    NOS <- NCOL(eta) / M1

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
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
      "pobs.a"     =     pobs.a,  # P(Y is altered)
      "onempobs.a" = 1 - pobs.a)  # P(Y is not altered)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .llambda = llambda,
           .elambda = elambda,
           .alter = alter ))),
  last = eval(substitute(expression({
    temp.names <- c(rep_len( "multinomial" , lalter),
                    rep_len( .llambda , NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]
    
    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- NULL  # ( .epobs.0 )
      misc$earg[[M1*ii  ]] <- .elambda
    }
  }), list( .llambda = llambda, .elambda = elambda,
            .alter = alter ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
    lambda <- cbind(eta2theta(eta[,  NCOL(eta), drop = FALSE],
                              .llambda, earg = .elambda ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgenapois(y, lambda = lambda, log = TRUE,
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
  vfamily = c("genapoisson"),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
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
                         inverse = TRUE)  # Am sure this is correct
    lambda <- eta2theta(eta[, ncol(eta)],
                        .llambda , earg = .elambda )
    rgenapois(nsim * length(lambda), lambda = lambda,
              pobs.a = pobs.a, alter = .alter )
  }, list( .llambda = llambda, .elambda = elambda,
           .alter = alter ))),


  deriv = eval(substitute(expression({
    alter <- ( .alter )
    M1 <- length(alter) + 1
    NOS <- ncol(eta) / M1  # extra$NOS
    y0 <- extra$y0
    skip <- extra$skip.these
    is.altered <- rowSums(skip) > 0  # TRUE if (any) y %in% avec

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
    lambda <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                        .llambda , earg = .elambda )


    sump <-
    Sume <- matrix(0, NROW(eta), NOS)
    for (aval in alter) {
      pmf <- dpois(aval, lambda)
      sump <- sump + pmf
      Sume <- Sume + pmf * aval
    }
    pmf.deriv1 <- function(y, lambda)
      dpois(y-1, lambda) - dpois(y, lambda)
    pmf.deriv2 <- function(y, lambda)
      dpois(y-2, lambda) - 2 * dpois(y-1, lambda) + dpois(y, lambda)
    sumderiv1 <- matrix(0, NROW(eta), NOS)
    sumderiv2 <- matrix(0, NROW(eta), NOS)
    for (aval in alter) {
      sumderiv1 <- sumderiv1 + pmf.deriv1(aval, lambda)
      sumderiv2 <- sumderiv2 + pmf.deriv2(aval, lambda)
    }

    lalter <- length(alter)
    pobs.a <- rowSums(phimat[, -NCOL(phimat), drop = FALSE])  # P(altered)
    onempobs.a <- 1 - pobs.a
    temp4 <- 1 / onempobs.a

    dl.deta <- skip - phimat[, -M, drop = FALSE]

    dl.dlambda <- (1 - is.altered) *
      (y/lambda - 1 + sumderiv1 / (1 - sump))
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    ans <- cbind(c(w) * dl.deta,
                 c(w) * dl.dlambda * dlambda.deta)
    ans
  }), list( .llambda = llambda, .elambda = elambda,
            .alter = alter ))),


  weight = eval(substitute(expression({
    MM12 <- M1 * (M1 + 1) / 2
    wz <- matrix(0.0, n, MM12)  # A full matrix

      ned2l.dlambda2 <- onempobs.a *
        ((lambda - Sume) / ((1 - sump) * lambda^2) -
         sumderiv2 / (1 - sump) -
        (sumderiv1 / (1 - sump))^2)
      wz[, iam(M1, M1, M = M1)] <- c(w) * ned2l.dlambda2 * dlambda.deta^2

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
    if (lalter > 0)
      for (jay in seq(lalter))
        for (kay in jay:lalter)
           wz[, iam(jay, kay, M = M1)] <- c(w) *
          wz4[, iam(jay, kay, M = M1-1)]
    wz
  }), list( .alter = alter ))))
}  # End of genapoisson







dgenipois <-
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
  if (length(x)      != LLL) x      <- rep(x,      length = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, length = LLL)

  linflate <- length(inflate)
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)
  sumi <- .rowSums(pstr.i, LLL, linflate)
  if (any(1 < sumi))
    stop("bad input for argument 'pstr.i'")


  numer1 <- 1 - sumi

  pmf <- numer1 * dpois(x, lambda) 
  for (jay in seq_along(inflate)) {
    ival <- inflate[jay]
    if (any(vecTF <- ival == x)) {
      pmf[vecTF] <- pmf[vecTF] + pstr.i[vecTF, jay]
    }
  }

  if (log.arg) log(pmf) else pmf
}  # dgenipois





pgenipois <-
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
  if (length(q)      != LLL) q      <- rep(q,      length = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, length = LLL)

  linflate <- length(inflate)
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)

  sumi <- rowSums(pstr.i)
  if (any(1 < sumi))
    stop("bad input for argument 'pstr.i'")
  numer1 <- 1 - sumi

  offset <- rep(0, length = LLL)
  for (jay in seq_along(inflate)) {
    ival <- inflate[jay]
    if (any(vecTF <- ival <= q)) {
      offset[vecTF] <- offset[vecTF] + pstr.i[vecTF, jay]
    }
  }
  ans <- numer1 * ppois(q, lambda) + offset
  ans
}  # pgenipois










rgenipois <-
  function(n, lambda, inflate = 0,  # NULL
           pstr.i = 0, byrow.arg = FALSE) {
  use.n <- if ((length.n <- length(n)) > 1) 
    length.n else
      if (!is.Numeric(n, integer.valued = TRUE, 
                      length.arg = 1, positive = TRUE))
        stop("bad input for argument 'n'") else n
    
  if (is.null(inflate))
    return(rpois(use.n, lambda))
    
  if (is.list(inflate))
    inflate <- unlist(inflate)
    
  if (!is.Numeric(inflate, integer.valued = TRUE) ||
      any(inflate < 0))
    stop("bad input for argument 'inflate'")
  if (!identical(inflate, (unique(inflate))))
    stop("values of argument 'inflate' must be unique")

  if (any(pstr.i < 0) || any(1 < pstr.i))
    stop("bad input for argument 'pstr.i'")
  
  linflate <- length(inflate)  # >= 1
  LLL <- use.n
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)

  sumi <- rowSums(pstr.i)
  if (any(1 < sumi))
    stop("bad input for argument 'pstr.i'")


  lambda <- rep(lambda, length = use.n)


  cpstr.i <- cbind(tapplymat1(pstr.i, "cumsum"), 1)
  choice <- rep(1, use.n)
  rrr <- runif(use.n)
  for (jay in 2:(linflate + 1)) {
    jayvec <- rep(jay, use.n)
    choice <- ifelse(cpstr.i[, jay-1] < rrr & rrr <= cpstr.i[, jay],
                     jayvec, choice)
  }

  ans <- inflate[choice]
  ind7 <- choice == (linflate + 1)
  if (any(ind7))
    ans[ind7] <- rpois(sum(ind7), lambda = lambda[ind7])
  ans
}  # rgenipois








 genipoisson <-
  function(inflate = 0,
           zero = NULL,
           llambda = "loglink",
           type.fitted = c("mean", "lambda", "pstr.i", "onempstr.i"),
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
  }
 
  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "pstr.i", "onempstr.i"))[1]
  temp7 <- paste("pstr", inflate, sep = "")


  new("vglmff",
  blurb = c("Generally-inflated Poisson distribution\n",
            "(multinomial and Poisson distribution)\n\n",
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
         link = c(if (length( .inflate ))
                  rep("multilogitlink",  length( .inflate )) else NULL,
                  .llambda ),
         link1parameter = if (length( .inflate ))
                          FALSE else TRUE,  # mulilogitlink is multiparameter
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
           .inflate = inflate
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
    for (jay in seq_along(inflate))
      extra$y0[, jay] <- y0[, jay] <- as.numeric(y == inflate[jay])
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)
    if (any((css <- colSums(y0)) == 0))
      stop("some 'inflate' argument values have no response values: ",
           paste(inflate[css == 0], collapse = ", "))          

    
    fillerChar <- " "
    use.refLevel <- M1+1  # Assumes only one response
    allbut.refLevel <- (1:(M+1))[-use.refLevel]
    predictors.names <-
      paste("log(prphi[,", allbut.refLevel,
            "]", fillerChar, "/", fillerChar, "prphi[,",
            use.refLevel, "])", sep = "")

    temp7 <- paste("pstr", inflate, sep = "")
    mynames1 <- paste("multinomial(", temp7, ")", sep = "")
    mynames2 <- param.names("lambda", ncoly, skip1 = TRUE)
    predictors.names <-
        c(        mynames1,
          namesof(mynames2, .llambda , earg = .elambda , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                             imu = .ilambda ,
                             ishrinkage = .ishrinkage ,
                             pos.only = TRUE,
                             probs.y = .probs.y )

      phimat <- colSums(c(w) * y0) / c(colSums(w))  # Weighted by 'w'.
      phimat <- phimat * ( .mux.inflate )
      phimat <- matrix(phimat, n, linflate, byrow = TRUE)
      etastart <-  multilogitlink(cbind(phimat, 1 - rowSums(phimat)))
      etastart <-
        cbind(etastart,
              theta2eta(lambda.init, .llambda , earg = .elambda ))
 #   etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .llambda = llambda,
            .elambda = elambda,
            .ipstr0 = ipstr0, .ilambda = ilambda,
            .mux.inflate = mux.inflate,
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .imethod = imethod,
            .inflate = inflate,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "pstr.i", "onempstr.i"))[1]

    inflate <- ( .inflate )
    M1 <- length(inflate) + 1
    NOS <- NCOL(eta) / M1

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
    lambda <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                        .llambda , earg = .elambda )
    pstr.i <- rowSums(phimat[, -NCOL(eta), drop = FALSE])
    Sume <- matrix(0, NROW(eta), NOS)  # Different defn of Sume
    for (jay in length(inflate)) {
      Sume <- Sume + phimat[, jay] * inflate[jay]
    }

    ans <- switch(type.fitted,
      "mean"       = (1 - pstr.i) * lambda + Sume,
      "lambda"     = lambda,
      "pstr.i"     =     pstr.i,  # P(Y is structurally inflated)
      "onempstr.i" = 1 - pstr.i)  # P(Y is not structurally inflated)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .llambda = llambda,
           .elambda = elambda,
           .inflate = inflate ))),
  last = eval(substitute(expression({
    temp.names <- c(rep_len( "multinomial" , linflate),
                    rep_len( .llambda , NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]
    
    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- NULL  # ( .epstr.0 )
      misc$earg[[M1*ii  ]] <- .elambda
    }
  }), list( .llambda = llambda, .elambda = elambda,
            .inflate = inflate ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pstr.i <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
    lambda <- cbind(eta2theta(eta[, NCOL(eta), drop = FALSE],
                              .llambda, earg = .elambda ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgenipois(y, lambda = lambda, log.arg = TRUE,
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
  vfamily = c("genipoisson"),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pstr.i <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
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
                         inverse = TRUE)  # Am sure this is correct
    lambda <- eta2theta(eta[, ncol(eta)],
                        .llambda , earg = .elambda )
    rgenipois(nsim * length(lambda), lambda = lambda,
              pstr.i = pstr.i, inflate = .inflate )
  }, list( .llambda = llambda, .elambda = elambda,
           .inflate = inflate ))),

  deriv = eval(substitute(expression({
    inflate <- ( .inflate )
    M1 <- length(inflate) + 1
    NOS <- ncol(eta) / M1  # extra$NOS
    y0 <- extra$y0  # Restore the value
    is.inflated <- rowSums(y0) > 0  # TRUE if (any) y %in% ivec

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
    lambda <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                        .llambda , earg = .elambda )
    linflate <- length(inflate)
    pstr.i <- rowSums(phimat[, -NCOL(phimat), drop = FALSE])  #P(inflated)
    onempstr.i <- 1 - pstr.i


    pmf.deriv1 <- function(y, lambda)
      dpois(y-1, lambda) - dpois(y, lambda)
    pmf.deriv2 <- function(y, lambda)
      dpois(y-2, lambda) - 2 * dpois(y-1, lambda) +
      dpois(y  , lambda)

    dl.deta <- -phimat[, -M, drop = FALSE]
    for (jay in seq_along(inflate)) {
      dl.deta[, jay] <- dl.deta[, jay] + y0[, jay] * phimat[, jay] / (
          phimat[, jay] + onempstr.i * dpois(inflate[jay], lambda))
    }  # for jay

    dl.dlambda <- (1 - is.inflated) * (y / lambda - 1)
    for (jay in seq_along(inflate)) {
      dl.dlambda <- dl.dlambda + y0[, jay] *
          onempstr.i * pmf.deriv1(inflate[jay], lambda) / (
          phimat[, jay] + onempstr.i * dpois(inflate[jay], lambda))
    }  # for jay
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    ans <- cbind(c(w) * dl.deta,
                 c(w) * dl.dlambda * dlambda.deta)
    ans
  }), list( .llambda = llambda, .elambda = elambda,
            .inflate = inflate ))),


  weight = eval(substitute(expression({
    MM12 <- M1 * (M1 + 1) / 2
    wz <- matrix(0.0, n, MM12)  # A full matrix

  if (FALSE) { 
      ned2l.dlambda2.old <- onempstr.i / lambda
      for (jay in seq_along(inflate))
        ned2l.dlambda2.old <- ned2l.dlambda2.old -
          (1 - inflate[jay] / lambda)^2 *
          phimat[, jay] * onempstr.i * dpois(inflate[jay], lambda) / (
          phimat[, jay] + onempstr.i * dpois(inflate[jay], lambda))
  }


  if (TRUE) { 
      ned2l.dlambda2 <- 1 / lambda
      for (jay in seq_along(inflate))
        ned2l.dlambda2 <- ned2l.dlambda2 -
          dpois(inflate[jay], lambda) * inflate[jay] / lambda^2
      for (jay in seq_along(inflate))
        ned2l.dlambda2 <- ned2l.dlambda2 -
          pmf.deriv2(inflate[jay], lambda) +
          onempstr.i  * (pmf.deriv1(inflate[jay], lambda)^2) / (
          phimat[, jay] + onempstr.i * dpois(inflate[jay], lambda))
      ned2l.dlambda2 <- ned2l.dlambda2 * onempstr.i
  }

      wz[, iam(M1, M1, M = M1)] <- c(w) * ned2l.dlambda2 * dlambda.deta^2

      for (jay in seq_along(inflate)) {
        wz[, iam(jay, M1, M = M1)] <- c(w) * dlambda.deta *
        phimat[, jay] * onempstr.i * pmf.deriv1(inflate[jay], lambda) / (
        phimat[, jay] + onempstr.i * dpois(inflate[jay], lambda))
      }  # for jay


    use.refLevel <- M
    if (linflate > 0) {
      index <- iam(NA, NA, M = M-1, both = TRUE, diag = TRUE)
      wz4 <- -phimat[, index$row, drop = FALSE] *
              phimat[, index$col, drop = FALSE]
      wz4[, 1:linflate] <- wz4[, 1:linflate] + phimat[, -M]
      for (jay in seq_along(inflate)) {
        wz4[, jay] <- wz4[, jay] -
          phimat[, jay] * onempstr.i * dpois(inflate[jay], lambda) / (
          phimat[, jay] + onempstr.i * dpois(inflate[jay], lambda))
      }  # for jay
    }
    wz4 <- as.matrix(wz4)  #  Needed when linflate == 1.
    if (linflate > 0)
      for (jay in seq(linflate))
        for (kay in jay:linflate)
           wz[, iam(jay, kay, M = M1  )] <- c(w) *
          wz4[, iam(jay, kay, M = M1-1)]
    wz
  }), list( .inflate = inflate ))))
}  # End of genipoisson






dgenibinom <-
  function(x, size, prob, inflate = 0,
           pstr.i = 0, byrow.arg = FALSE,
           log.arg = FALSE) {
  if (is.null(inflate))
    return(dbinom(x, size, prob, log = log.arg))
    
  if (is.list(inflate))
    inflate <- unlist(inflate)
    
  if (!is.Numeric(inflate, integer.valued = TRUE) ||
      any(inflate < 0) ||
      any(size < inflate))
    stop("bad input for argument 'inflate'")
  if (!identical(inflate, (unique(inflate))))
    stop("values of argument 'inflate' must be unique")

  if (any(pstr.i < 0) || any(1 < pstr.i))
    stop("bad input for argument 'pstr.i'")

  LLL <- max(length(x), length(size), length(prob))  # , length(pstr.i)
  if (length(x)      != LLL) x      <- rep(x,      length = LLL)
  if (length(size  ) != LLL) size   <- rep(size,   length = LLL)
  if (length(prob  ) != LLL) prob   <- rep(prob,   length = LLL)

  linflate <- length(inflate)
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)
  sumi <- .rowSums(pstr.i, LLL, linflate)
  if (any(1 < sumi))
    stop("bad input for argument 'pstr.i'")


  numer1 <- 1 - sumi

  pmf <- numer1 * dbinom(x, size, prob) 
  for (jay in seq_along(inflate)) {
    ival <- inflate[jay]
    if (any(vecTF <- ival == x)) {
      pmf[vecTF] <- pmf[vecTF] + pstr.i[vecTF, jay]
    }
  }

  if (log.arg) log(pmf) else pmf
}  # dgenibinom








pgenibinom <-
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
  if (length(q)      != LLL) q      <- rep(q,      length = LLL)
  if (length(size)   != LLL) size   <- rep(size,   length = LLL)
  if (length(prob)   != LLL) prob   <- rep(prob,   length = LLL)

  linflate <- length(inflate)
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)

  sumi <- rowSums(pstr.i)
  if (any(1 < sumi))
    stop("bad input for argument 'pstr.i'")
  numer1 <- 1 - sumi

  offset <- rep(0, length = LLL)
  for (jay in seq_along(inflate)) {
    ival <- inflate[jay]
    if (any(vecTF <- ival <= q)) {
      offset[vecTF] <- offset[vecTF] + pstr.i[vecTF, jay]
    }
  }
  ans <- numer1 * pbinom(q, size, prob) + offset
  ans
}  # pgenibinom







rgenibinom <-
  function(n, size, prob, inflate = 0,  # NULL
           pstr.i = 0, byrow.arg = FALSE) {
  use.n <- if ((length.n <- length(n)) > 1)
    length.n else
      if (!is.Numeric(n, integer.valued = TRUE, 
                      length.arg = 1, positive = TRUE))
        stop("bad input for argument 'n'") else n
    
  if (is.null(inflate))
    return(rbinom(use.n, size, prob))
    
  if (is.list(inflate))
    inflate <- unlist(inflate)
    
  if (!is.Numeric(inflate, integer.valued = TRUE) ||
      any(inflate < 0))
    stop("bad input for argument 'inflate'")
  if (!identical(inflate, (unique(inflate))))
    stop("values of argument 'inflate' must be unique")

  if (any(pstr.i < 0) || any(1 < pstr.i))
    stop("bad input for argument 'pstr.i'")
  
  linflate <- length(inflate)  # >= 1
  LLL <- use.n
  pstr.i <- matrix(pstr.i, LLL, linflate, byrow = byrow.arg)

  sumi <- rowSums(pstr.i)
  if (any(1 < sumi))
    stop("bad input for argument 'pstr.i'")

  size <- rep(size, length = use.n)
  prob <- rep(prob, length = use.n)

  cpstr.i <- cbind(tapplymat1(pstr.i, "cumsum"), 1)
  choice <- rep(1, use.n)
  rrr <- runif(use.n)
  for (jay in 2:(linflate + 1)) {
    jayvec <- rep(jay, use.n)
    choice <- ifelse(cpstr.i[, jay-1] < rrr & rrr <= cpstr.i[, jay],
                     jayvec, choice)
  }

  ans <- inflate[choice]
  ind7 <- choice == (linflate + 1)
  if (any(ind7))
    ans[ind7] <- rbinom(sum(ind7), size[ind7], prob[ind7])
  ans
}  # rgenibinom








dgenabinom <-
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
      any(size < alter))
    stop("bad input for argument 'alter'")
  if (!identical(alter, (unique(alter))))
    stop("values of argument 'alter' must be unique")

  if (any(pobs.a < 0) || any(1 < pobs.a))
    stop("bad input for argument 'pobs.a'")

  LLL <- max(length(x), length(size), length(prob))
  if (length(x)      != LLL) x      <- rep(x,      length = LLL)
  if (length(size  ) != LLL) size   <- rep(size,   length = LLL)
  if (length(prob  ) != LLL) prob   <- rep(prob,   length = LLL)

  lalter <- length(alter)
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  suma <- rowSums(pobs.a)
  if (any(1 < suma))
    stop("bad input for argument 'pobs.a'")

  if (log.arg) {
    logpmf <- log1p(-suma) +
              dgentbinom(x, size, prob, truncate = alter, log = TRUE)
  } else {
    pmf <- exp(log1p(-suma)) * dgentbinom(x, size, prob, truncate = alter)
  }
 for (jay in seq_along(alter)) {
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
}  # dgenabinom







pgenabinom <-
  function(q, size, prob, alter = 0,
           pobs.a = 0, byrow.arg = FALSE) {
  if (is.null(alter))
    return(pbinom(q, size, prob))
    
  if (is.list(alter))
    truncate <- unlist(alter)
    
  if (!is.Numeric(alter, integer.valued = TRUE) ||
      any(alter < 0) ||
      any(size < alter))
    stop("bad input for argument 'alter'")
  if (!identical(alter, (unique(alter))))
    stop("values of argument 'alter' must be unique")

  if (any(pobs.a < 0) || any(1 < pobs.a))
    stop("bad input for argument 'pobs.a'")

  LLL <- max(length(q), length(size), length(prob))
  if (length(q)      != LLL) q      <- rep(q,      length = LLL)
  if (length(size  ) != LLL) size   <- rep(size,   length = LLL)
  if (length(prob  ) != LLL) prob   <- rep(prob,   length = LLL)

  lalter <- length(alter)
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  suma <- rowSums(pobs.a)
  if (any(1 < suma))
    stop("bad input for argument 'pobs.a'")
  numer1 <- 1 - suma

  offset <- rep(0, length = LLL)
  for (jay in seq_along(alter)) {
    aval <- alter[jay]
    if (any(vecTF <- aval <= q)) {
      offset[vecTF] <- offset[vecTF] + pobs.a[vecTF, jay]
    }
  }
  ans <- numer1 * pgentbinom(q, size, prob, truncate = alter) + offset
  ans
}  # pgenabinom





rgenabinom <-
  function(n, size, prob, alter = 0,  # NULL
           pobs.a = 0, byrow.arg = FALSE,
           maxits = 10000) {
  use.n <- if ((length.n <- length(n)) > 1) 
    length.n else
      if (!is.Numeric(n, integer.valued = TRUE, 
                      length.arg = 1, positive = TRUE))
        stop("bad input for argument 'n'") else n
    
  if (is.null(alter))
    return(rbinom(use.n, size, prob))
    
  if (is.list(alter))
    alter <- unlist(alter)
    
  if (!is.Numeric(alter, integer.valued = TRUE) ||
      any(alter < 0) ||
      any(size < alter))
    stop("bad input for argument 'alter'")
  if (!identical(alter, (unique(alter))))
    stop("values of argument 'alter' must be unique")
    

  lalter <- length(alter)
  LLL <- use.n
  pobs.a <- matrix(pobs.a, LLL, lalter, byrow = byrow.arg)

  suma <- rowSums(pobs.a)
  if (any(1 < suma))
    stop("bad input for argument 'pobs.a'")

  pobs.a <- cbind(pobs.a, 1 - suma)


  size <- rep(size, length = use.n)
  prob <- rep(prob, length = use.n)


  cpobs.a <- tapplymat1(pobs.a, "cumsum")
  choice <- rep(1, use.n)
  rrr <- runif(use.n)
  if (lalter > 1) {
    for (jay in 2:(lalter + 1)) {
       jayvec <- rep(jay, use.n)
       choice <- ifelse(cpobs.a[, jay-1] < rrr & rrr <= cpobs.a[, jay],
                        jay, choice)
    }
  }

  ans <- alter[choice]
  ind7 <- choice == (lalter + 1)
  if (any(ind7))
    ans[ind7] <- rgentbinom(sum(ind7), size[ind7], prob[ind7],
                            truncate = alter, maxits = maxits)
  ans
}  # rgenabinom






 genabinomial <-
  function(alter = 0,  # NULL,  # 0  #,
           zero = NULL,  # Was zero = 2 prior to 20130917
           lprob  = "logitlink",
           type.fitted = c("mean", "prob", "pobs.a"),
           imethod = 1,
           iprob = NULL
          ) {

 print("20181203; does not work!")


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
                           c("mean", "prob", "pobs.a"))[1]
  temp7 <- paste("pobs", alter, sep = "")


  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Generally-altered binomial distribution\n",
            "(multinomial and truncated-binomial conditional model)\n\n",
            "Links:    ", if (length(alter))
      paste("multilogitlink(cbind(",
            paste(temp7, collapse = ", "),
            ", ",
            "\n                           ",
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
         Q1 = 1,
         link = c(if (length( .alter ))
                  rep("multilogitlink",  length( .alter )) else NULL,
                  prob = .lprob ),
         link1parameter = if (length( .alter ))
                          FALSE else TRUE,  # mulilogitlink is multiparameter
         alter = .alter ,
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE,
         parameters.names = c(if (length( .alter )) temp7 else NULL,
                              "prob"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted,
           .lprob = lprob, .eprob = eprob,
           .alter = alter
         ))),

  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w



    if (NCOL(y) == 1) {
      if (is.factor(y))
        y <- y != levels(y)[1]
      nn <- rep_len(1, n)
      if (!all(y >= 0 & y <= 1))
        stop("response values must be in [0, 1]")
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + w * y) / (1.0 + w)


      no.successes <- y
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(no.successes - round(no.successes)) > 1.0e-8))
        stop("Number of successes must be integer-valued")

    } else if (NCOL(y) == 2) {
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(y - round(y)) > 1.0e-8))
        stop("Count data must be integer-valued")
      y <- round(y)
      nvec <- y[, 1] + y[, 2]
      y <- ifelse(nvec > 0, y[, 1] / nvec, 0)
      w <- w * nvec
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + nvec * y) / (1 + nvec)
    } else {
      stop("for the binomialff family, response 'y' must be a ",
           "vector of 0 and 1's\n",
           "or a factor ",
           "(first level = fail, other levels = success),\n",
           "or a 2-column matrix where col 1 is the no. of ",
           "successes and col 2 is the no. of failures")
    }
    if (!all(w == 1))
      extra$new.w <- w


    y <- as.matrix(y)




    alter <- ( .alter )
    lalter <- length(alter)
    M1 <- lalter + 1
    NOS <- 1  # Since multiple.responses = FALSE;  NCOL(y)
    M <- NOS * M1
    extra$NOS <- ncoly <- ncol(y)  # Number of species
 print("ncoly")
 print( ncoly )


 print("head(w)")
 print( head(w) )
 print("head(y)")
 print( head(y) )

    extra$y0 <- y0 <- matrix(0, n, lalter)
    for (jay in seq_along(alter))
      extra$y0[, jay] <- y0[, jay] <- as.numeric(y %in% alter[jay])
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, lalter)
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)
    if (any((css <- colSums(skip.these)) == 0))
      stop("some 'alter' argument values have no response values: ",
           paste(alter[css == 0], collapse = ", "))          
 print("head(y0)")
 print( head(y0) )
 print("head(skip.these)")
 print( head(skip.these) )

    
    fillerChar <- " "
    use.refLevel <- M1+1  # Assumes only one response
    allbut.refLevel <- (1:(M+1))[-use.refLevel]
    predictors.names <-
      paste("log(pobs.a[,", allbut.refLevel,
            "]", fillerChar, "/", fillerChar,
            "pobs.a[,",  # "promega[,",
            use.refLevel, "])", sep = "")

    temp7 <- paste("pobs", alter, sep = "")
    mynames1 <- paste("multinomial(", temp7, ")", sep = "")
    mynames2 <- param.names("prob", ncoly, skip1 = TRUE)
    predictors.names <-
        c(        mynames1,
          namesof(mynames2, .lprob , earg = .eprob , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {
      prob.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                             imu = .iprob ,
                             pos.only = TRUE
                           )

      phimat <- colMeans(skip.these)  # yettodo: weight this by 'w'
      phimat <- matrix(phimat, n, lalter, byrow = TRUE)
      etastart <-  multilogitlink(cbind(phimat, 1 - rowSums(phimat)))
      etastart <-
        cbind(etastart,
              theta2eta(prob.init, .lprob , earg = .eprob ))
 #   etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
 print("head(etastart)1")
 print( head(etastart) )
    mustart <- NULL  # 20181203
    }
  }), list( .lprob = lprob,
            .eprob = eprob,
            .iprob = iprob,
            .imethod = imethod,
            .alter = alter,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w
 print("head(Size) in @linkinv")
 print( head(Size) )

   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "pobs.a", "onempobs.a"))[1]

    alter <- ( .alter )
    M1 <- length(alter) + 1
    NOS <- NCOL(eta) / M1

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
    prob <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                        .lprob , earg = .eprob )
    pobs.a <- rowSums(phimat[, -NCOL(eta), drop = FALSE])
    sump <- matrix(0, NROW(eta), NOS)
    Sume <- matrix(0, NROW(eta), NOS)
    for (aval in alter) {
      pmf <- dbinom(aval, size = Size, prob)
      sump <- sump + pmf
      Sume <- Sume + pmf * aval
    }

    ans <- switch(type.fitted,
      "mean"       = (1 - pobs.a) * (prob - Sume) / (1 - sump),
      "prob"       = prob,
      "pobs.a"     =     pobs.a,  # P(Y is altered)
      "onempobs.a" = 1 - pobs.a)  # P(Y is not altered)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lprob = lprob,
           .eprob = eprob,
           .alter = alter ))),
  last = eval(substitute(expression({
    temp.names <- c(rep_len( "multinomial" , lalter),
                    prob = rep_len( .lprob , NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]
    
    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- NULL  # ( .epobs.0 )
      misc$earg[[M1*ii  ]] <- .eprob
    }
  }), list( .lprob = lprob, .eprob = eprob,
            .alter = alter ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w
 print("head(Size) in @loglikelihood")
 print( head(Size) )

    pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
    prob <- cbind(eta2theta(eta[,  NCOL(eta), drop = FALSE],
                              .lprob, earg = .eprob ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgenabinom(round(y * Size),  # Was y, prior to 20181203
                   prob = prob, log = TRUE,
                   alter = .alter ,  # byrow.arg = FALSE,
                   size = Size,
                   pobs.a = pobs.a[, -NCOL(pobs.a), drop = FALSE])

      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .eprob = eprob,
           .alter = alter ))),
  vfamily = c("genabinomial"),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w
    okay2 <- max( .alter ) <= max(Size)
 print("head(Size) in @validparams")
 print( head(Size) )
    pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
    prob <- eta2theta(eta[,  NCOL(eta), drop = FALSE],
                        .lprob , earg = .eprob )
    okay1 <- all(is.finite(prob))   && all(0 < prob   & prob   < 1) &&
             all(is.finite(pobs.a)) && all(0 < pobs.a & pobs.a < 1)
    okay1 && okay2
  }, list( .lprob = lprob, .eprob = eprob,
           .alter = alter ))),


  simslot = eval(substitute(
  function(object, nsim) {
    extra <- object@extra
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs.a <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
    prob <- eta2theta(eta[, ncol(eta)],
                      .lprob , earg = .eprob )
    rgenabinom(nsim * length(prob), prob = prob,
               size = Size,
               pobs.a = pobs.a, alter = .alter )
  }, list( .lprob = lprob, .eprob = eprob,
           .alter = alter ))),


  deriv = eval(substitute(expression({
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w
    
    alter <- ( .alter )
    M1 <- length(alter) + 1
    NOS <- ncol(eta) / M1  # extra$NOS
 print("NOS in @deriv")
 print( NOS )
    y0 <- extra$y0
    skip <- extra$skip.these
    is.altered <- rowSums(skip) > 0  # TRUE if (any) y %in% avec

    phimat <- multilogitlink(eta[, -NCOL(eta), drop = FALSE],
                         refLevel = NCOL(eta),  # Assumes one response
                         inverse = TRUE)  # Am sure this is correct
    prob <- eta2theta(eta[, NCOL(eta), drop = FALSE],
                      .lprob , earg = .eprob )


    sump <-
    Sume <- matrix(0, NROW(eta), NOS)
    for (aval in alter) {
      pmf <- dbinom(aval, size = Size, prob)
      sump <- sump + pmf
      Sume <- Sume + pmf * aval
    }

    pmf.deriv1 <- function(ycount, size, prob)
      dbinom(ycount, size, prob) * size *
      (ycount / size - prob) / (prob * (1 - prob))
    pmf.deriv2 <- function(ycount, size, prob)
      dbinom(ycount, size, prob) * size * (
      size * (((ycount / size) / prob) -
              (1 - ycount / size) / (1 - prob))^2 -
      (ycount / size) / prob^2 -
      (1 - ycount / size) / (1 - prob)^2)
    sumderiv1 <- matrix(0, NROW(eta), NOS)
    sumderiv2 <- matrix(0, NROW(eta), NOS)
    for (aval in alter) {
      sumderiv1 <- sumderiv1 + pmf.deriv1(aval, Size, prob)
      sumderiv2 <- sumderiv2 + pmf.deriv2(aval, Size, prob)
    }

    lalter <- length(alter)
    pobs.a <- rowSums(phimat[, -NCOL(phimat), drop = FALSE])  # P(altered)
    onempobs.a <- 1 - pobs.a
    temp4 <- 1 / onempobs.a

    dl.deta <- skip - phimat[, -M, drop = FALSE]

    dl.dprob <- (1 - is.altered) * Size *
      ((y / Size) / prob -
       (1 - y / Size) / (1 - prob) +
      sumderiv1 / (1 - sump))
    dprob.deta <- dtheta.deta(prob, .lprob , earg = .eprob )
    ans <- cbind(c(w) * dl.deta,
                 c(w) * dl.dprob * dprob.deta)
 print("head(ans) in @deriv")
 print( head(ans) )
    ans
  }), list( .lprob = lprob, .eprob = eprob,
            .alter = alter ))),


  weight = eval(substitute(expression({
    MM12 <- M1 * (M1 + 1) / 2
    wz <- matrix(0.0, n, MM12)  # A full matrix

      ned2l.dprob2 <- onempobs.a *
        ((Size / prob^2) *
         (prob - Sume) / (1 - sump) +
         (1 - (prob - Sume) / (1 - sump)) * Size / (1 - prob)^2 -
         sumderiv2 / (1 - sump) -
        (sumderiv1 / (1 - sump))^2)
      wz[, iam(M1, M1, M = M1)] <- c(w) * ned2l.dprob2 * dprob.deta^2

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
 print("head(wz4)")
 print( head(wz4) )
    if (lalter > 0)
      for (jay in seq(lalter))
        for (kay in jay:lalter)
           wz[, iam(jay, kay, M = M1)] <- c(w) *
          wz4[, iam(jay, kay, M = M1-1)]
 print("head(wz)")
 print( head(wz) )
    wz
  }), list( .alter = alter ))))
}  # genabinomial








if (FALSE)
 zabinomial <-
  function(lpobs0 = "logitlink",
           lprob  = "logitlink",
           type.fitted = c("mean", "prob", "pobs0"),
           ipobs0 = NULL, iprob = NULL,
           imethod = 1,
           zero = NULL  # Was zero = 2 prior to 20130917
          ) {



  lpobs0 <- as.list(substitute(lpobs0))
  epobs0 <- link2list(lpobs0)
  lpobs0 <- attr(epobs0, "function.name")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")


  type.fitted <- match.arg(type.fitted,
                           c("mean", "prob", "pobs0"))[1]

  if (length(ipobs0))
    if (!is.Numeric(ipobs0, positive = TRUE) ||
        ipobs0 >= 1)
      stop("argument 'ipobs0' is out of range")

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")



  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Zero-altered binomial distribution ",
            "(Bernoulli and positive-binomial conditional model)\n\n",
            "P[Y = 0] = pobs0,\n",
            "P[Y = y] = (1 - pobs0) * dposbinom(x = y, size, prob), ",
            "y = 1, 2, ..., size,\n\n",
            "Link:     ",
            namesof("pobs0",   lpobs0, earg = epobs0), ", ",
            namesof("prob" ,   lprob,  earg = eprob),  "\n",
            "Mean:     (1 - pobs0) * prob / (1 - (1 - prob)^size)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = NA,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("pobs0", "prob"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted ))),

  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w



    if (NCOL(y) == 1) {
      if (is.factor(y))
        y <- y != levels(y)[1]
      nn <- rep_len(1, n)
      if (!all(y >= 0 & y <= 1))
        stop("response values must be in [0, 1]")
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + w * y) / (1.0 + w)


      no.successes <- y
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(no.successes - round(no.successes)) > 1.0e-8))
        stop("Number of successes must be integer-valued")

    } else if (NCOL(y) == 2) {
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(y - round(y)) > 1.0e-8))
        stop("Count data must be integer-valued")
      y <- round(y)
      nvec <- y[, 1] + y[, 2]
      y <- ifelse(nvec > 0, y[, 1] / nvec, 0)
      w <- w * nvec
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + nvec * y) / (1 + nvec)
    } else {
      stop("for the binomialff family, response 'y' must be a ",
           "vector of 0 and 1's\n",
           "or a factor ",
           "(first level = fail, other levels = success),\n",
           "or a 2-column matrix where col 1 is the no. of ",
           "successes and col 2 is the no. of failures")
    }
    if (!all(w == 1))
      extra$new.w <- w


    y <- as.matrix(y)
    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)

    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)


    predictors.names <-
        c(namesof("pobs0", .lpobs0 , earg = .epobs0 , tag = FALSE),
          namesof("prob" , .lprob  , earg = .eprob  , tag = FALSE))



    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    phi.init <- if (length( .ipobs0 )) .ipobs0 else {
        prob0.est <- sum(Size[y == 0]) / sum(Size)
        if ( .imethod == 1) {
          (prob0.est - (1 - mustart)^Size) / (1 - (1 - mustart)^Size)
        } else
        if ( .imethod == 2) {
          prob0.est
        } else {
          prob0.est * 0.5
        }
    }

    phi.init[phi.init <= -0.10] <- 0.50  # Lots of sample variation
    phi.init[phi.init <=  0.01] <- 0.05  # Last resort
    phi.init[phi.init >=  0.99] <- 0.95  # Last resort




    if (!length(etastart)) {
      etastart <-
        cbind(theta2eta(phi.init, .lpobs0, earg = .epobs0 ),
              theta2eta( mustart, .lprob,  earg = .eprob  ))


      mustart <- NULL
    }
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0,
            .iprob = iprob, .ipobs0 = ipobs0,
            .imethod = imethod,
            .type.fitted = type.fitted ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- ncol(eta) / c(M1 = 2)
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "prob", "pobs0"))[1]

    phi0  <- eta2theta(eta[, 1], .lpobs0 , earg = .epobs0 )
    prob  <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    ans <- switch(type.fitted,
                  "mean"      = (1 - phi0) * prob / (1 - (1 - prob)^Size),
                  "prob"      = prob,
                  "pobs0"     = phi0)  # P(Y=0)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lprob = lprob, .lpobs0 = lpobs0,
           .eprob = eprob, .epobs0 = epobs0 ))),

  last = eval(substitute(expression({
    misc$link <-    c(prob = .lprob, pobs0 = .lpobs0 )
    misc$earg <- list(prob = .eprob, pobs0 = .epobs0 )

    misc$imethod  <- .imethod
    misc$zero     <- .zero
    misc$expected <- TRUE
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0,
            .zero = zero,
            .imethod = imethod ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w
    pobs0 <- eta2theta(eta[, 1], .lpobs0 , earg = .epobs0 )
    prob  <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        orig.w * dzabinom(x = round(y * Size), size = Size,
                          prob = prob, pobs0 = pobs0,
                          log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .lpobs0 = lpobs0,
           .eprob = eprob, .epobs0 = epobs0 ))),
  vfamily = c("zabinomial"),


  validparams = eval(substitute(function(eta, y, extra = NULL) {
    phi0 <- eta2theta(eta[, 1], .lpobs0 , earg = .epobs0 )
    prob <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )
    okay1 <- all(is.finite(phi0)) && all(0 < phi0 & phi0 < 1) &&
             all(is.finite(prob)) && all(0 < prob & prob < 1)
    okay1
  }, list( .lprob = lprob, .lpobs0 = lpobs0,
           .eprob = eprob, .epobs0 = epobs0 ))),


  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- if (length(extra$NOS)) extra$NOS else 1

    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    phi0 <- eta2theta(eta[, 1], .lpobs0 , earg = .epobs0 )
    prob <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )

    dphi0.deta <- dtheta.deta(phi0, .lpobs0, earg = .epobs0 )
    dprob.deta <- dtheta.deta(prob, .lprob , earg = .eprob  )

    df0.dprob   <- -Size *              (1 -  prob)^(Size - 1)
    df02.dprob2 <-  Size * (Size - 1) * (1 -  prob)^(Size - 2)
    prob0  <- (1 -  prob)^(Size)
    oneminusf0  <- 1 - prob0


    dl.dphi0 <- -1 / (1 - phi0)
    dl.dprob <-  c(w)      * (y / prob - (1 - y) / (1 - prob)) +
                 c(orig.w) * df0.dprob / oneminusf0


    dl.dphi0[y == 0] <- 1 / phi0[y == 0]  # Do it in one line
    skip <- extra$skip.these
    for (spp. in 1:NOS) {
      dl.dprob[skip[, spp.], spp.] <- 0
    }


    ans <- cbind(c(orig.w) * dl.dphi0 * dphi0.deta,
                             dl.dprob * dprob.deta)


    ans
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0 ))),


  weight = eval(substitute(expression({
    wz <- matrix(0.0, n, M1)

    usualmeanY <-  prob
    meanY <- (1 - phi0) * usualmeanY / oneminusf0


    term1 <-  c(Size) * (meanY /      prob^2 -
                         meanY / (1 - prob)^2) +
             c(Size) * (1 - phi0) / (1 - prob)^2

    term2 <-  -(1 - phi0) * df02.dprob2 / oneminusf0
    term3 <-  -(1 - phi0) * (df0.dprob  / oneminusf0)^2
    ned2l.dprob2 <- term1 + term2 + term3
    wz[, iam(2, 2, M)] <- ned2l.dprob2 * dprob.deta^2


    mu.phi0 <- phi0
    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if ( .lpobs0 == "logitlink" && is.empty.list( .epobs0 )) {
      tmp100
    } else {
      (dphi0.deta^2) / tmp100
    }
    wz[, iam(1, 1, M)] <- tmp200


    c(orig.w) * wz
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0 ))))
}  #  zabinomial















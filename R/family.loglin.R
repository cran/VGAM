# These functions are
# Copyright (C) 1998-2025 T.W. Yee, University of Auckland.
# All rights reserved.






 loglinb2 <-
    function(exchangeable = FALSE, zero = "u12") {


  if (!isFALSE(exchangeable) && !isTRUE(exchangeable))
    warning("'exchangeable' is not a single logical")


  new("vglmff",
  blurb = c("Loglinear model for two binary responses\n\n",
            "Links:    ",
            "Identity: u1, u2, u12",
            "\n"),
  constraints = eval(substitute(expression({
    cm.int.default <- diag(3)  # cm.intercept.default

    constraints <-
      cm.VGAM(matrix(c(1,1,0, 0,0,1), 3, 2), x = x,
              constraints = constraints,
              apply.int = TRUE, bool = .exchangeable ,
              cm.default           = cm.int.default,
              cm.intercept.default = cm.int.default)
    constraints <-
    cm.zero.VGAM(constraints, x = x, .zero , M = M,
                 predictors.names = predictors.names,
                 M1 = 3)
  }),
  list( .exchangeable = exchangeable,
       .zero = zero ))),



  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 2,  # ncol(depvar(object))
         expected = TRUE,
         multipleResponses = FALSE,  # TRUE,
         parameters.names = c("u1", "u2", "u12"),
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = expression({
    predictors.names <- c("u1", "u2", "u12")
    Q1 <- 2


    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = Q1,
              out.wy = TRUE,
              colsyperw = Q1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    if (ncol(y) != Q1)
      stop("ncol(y) must be = ", Q1)

    if (length(mustart) + length(etastart) == 0) {
      mustart <- matrix(NA_real_, nrow(y), 4)
      mustart[,1] <- weighted.mean((1-y[,1]) * (1-y[,2]), w)
      mustart[,2] <- weighted.mean((1-y[,1]) *    y[,2] , w)
      mustart[,3] <- weighted.mean(   y[,1]  * (1-y[,2]), w)
      mustart[,4] <- weighted.mean(   y[,1]  *    y[,2] , w)
      if (any(mustart == 0))
        stop("some combinations of the response not realized")
    }
  }),
  linkinv = function(eta, extra = NULL) {
    u1 <-  eta[, 1]
    u2 <-  eta[, 2]
    u12 <- eta[, 3]
    denom <- 1 + exp(u1) + exp(u2) + exp(u1 + u2 + u12)
    cbind("00" = 1/denom,
          "01" = exp(u2) / denom,
          "10" = exp(u1) / denom,
          "11" = exp(u1+u2+u12) / denom)
  },
  last = expression({
    misc$link <- c("u1"  = "identitylink",
                   "u2"  = "identitylink",
                   "u12" = "identitylink")
    misc$earg <- list("u1"  = list(theta = NULL),
                      "u2"  = list(theta = NULL),
                      "u12" = list(theta = NULL))

  }),
  linkfun = function(mu, extra = NULL)  {
    u0 <-  log(mu[, 1])
    u2 <-  log(mu[, 2]) - u0
    u1 <-  log(mu[, 3]) - u0
    u12 <- log(mu[, 4]) - u0 - u1 - u2
    cbind(u1, u2, u12)
  },
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    u1 <-  eta[, 1]
    u2 <-  eta[, 2]
    u12 <- eta[, 3]
    denom <- 1 + exp(u1) + exp(u2) + exp(u1+u2+u12)
    u0 <- -log(denom)
    if (residuals) {
      stop("loglikelihood residuals not implemented")
    } else {
      ll.elts <- c(w) * (u0 + u1*y[, 1] +
                 u2*y[, 2] + u12*y[, 1]*y[, 2])
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("loglinb2"),
  validparams = function(eta, y, extra = NULL) {
    u1 <-  eta[, 1]
    u2 <-  eta[, 2]
    u12 <- eta[, 3]
    okay1 <- all(is.finite(u1 )) &&
             all(is.finite(u2 )) &&
             all(is.finite(u12))
    okay1
  },
  deriv = expression({
    u1 <-  eta[, 1]
    u2 <-  eta[, 2]
    u12 <- eta[, 3]
    denom <- 1 + exp(u1) + exp(u2) + exp(u1 + u2 + u12)
    du0.du1 <- -(exp(u1) + exp(u1 + u2 + u12)) / denom
    du0.du2 <- -(exp(u2) + exp(u1 + u2 + u12)) / denom
    du0.du12 <- -exp(u1 + u2 + u12) / denom
    c(w) * cbind(du0.du1  + y[, 1],
                 du0.du2  + y[, 2],
                 du0.du12 + y[, 1] * y[, 2])
  }),
  weight = expression({
    d2u0.du1.2 <- -(exp(u1) + exp(u1 + u2 + u12)) * (1+exp(u2)) / denom^2
    d2u0.du22 <-  -(exp(u2) + exp(u1 + u2 + u12)) * (1+exp(u1)) / denom^2
    d2u0.du122 <- -exp(u1 + u2 + u12) * (1+exp(u1)+exp(u2)) / denom^2
    d2u0.du1u2 <- -(exp(u1 + u2 + u12) - exp(u1 + u2)) / denom^2
    d2u0.du1u3 <- -(1 + exp(u2)) * exp(u1 + u2 + u12) / denom^2
    d2u0.du2u3 <- -(1 + exp(u1)) * exp(u1 + u2 + u12) / denom^2

    wz <- matrix(NA_real_, n, dimm(M))
    wz[,iam(1, 1, M)] <- -d2u0.du1.2
    wz[,iam(2, 2, M)] <- -d2u0.du22
    wz[,iam(3, 3, M)] <- -d2u0.du122
    wz[,iam(1, 2, M)] <- -d2u0.du1u2
    wz[,iam(1, 3, M)] <- -d2u0.du1u3
    wz[,iam(2, 3, M)] <- -d2u0.du2u3
    c(w) * wz
  }))
}  # loglinb2





 if (FALSE)
 loglinb3.orig <-
    function(exchangeable = FALSE,
             zero = c("u12", "u13", "u23")) {


  if (!isFALSE(exchangeable) && !isTRUE(exchangeable))
    warning("'exchangeable' should be a single logical")


  new("vglmff",
  blurb = c("Loglinear model for three binary responses\n\n",
            "Links:    ",
            "Identity: u1, u2, u3, u12, u13, u23",
            "\n"),
  constraints = eval(substitute(expression({
    cm.intercept.default <- diag(6)

    constraints <-
        cm.VGAM(matrix(c(1,1,1,0,0,0, 0,0,0,1,1,1), 6, 2),
                x = x,
                bool = .exchangeable ,
                constraints = constraints,
                apply.int = TRUE,
                cm.default           = cm.intercept.default,
                cm.intercept.default = cm.intercept.default)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 6)
  }),
  list( .exchangeable = exchangeable,
        .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 6,
         Q1 = 3,  # ncol(depvar(object))
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("u1", "u2", "u3", "u12", "u13", "u23"),
         zero = .zero )
  }, list( .zero = zero
         ))),


  initialize = expression({
    predictors.names <- c("u1", "u2", "u3", "u12", "u13", "u23")
    Q1 <- 3

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = Q1,
              out.wy = TRUE,
              colsyperw = Q1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    if (ncol(y) != Q1)
      stop("ncol(y) must be = ", Q1)


    if (FALSE)
    extra$my.expression <- expression({
      u1  <- eta[, 1]
      u2  <- eta[, 2]
      u3  <- eta[, 3]
      u12 <- eta[, 4]
      u13 <- eta[, 5]
      u23 <- eta[, 6]
      denom <- 1 + exp(u1) + exp(u2) + exp(u3) + exp(u1 + u2 + u12) +
               exp(u1 + u3 + u13) + exp(u2 + u3 + u23) +
               exp(u1 + u2 + u3 + u12 + u13 + u23)
    })



    if (length(mustart) + length(etastart) == 0) {
      mustart <- matrix(NA_real_, nrow(y), 2^3)
      mustart[,1] <- weighted.mean((1-y[,1])*(1-y[,2])*(1-y[,3]), w)
      mustart[,2] <- weighted.mean((1-y[,1])*(1-y[,2])*   y[,3] , w)
      mustart[,3] <- weighted.mean((1-y[,1])*   y[,2] *(1-y[,3]), w)
      mustart[,4] <- weighted.mean((1-y[,1])*   y[,2] *   y[,3] , w)
      mustart[,5] <- weighted.mean(   y[,1] *(1-y[,2])*(1-y[,3]), w)
      mustart[,6] <- weighted.mean(   y[,1] *(1-y[,2])*   y[,3] , w)
      mustart[,7] <- weighted.mean(   y[,1] *   y[,2] *(1-y[,3]), w)
      mustart[,8] <- weighted.mean(   y[,1] *   y[,2] *   y[,3] , w)
      if (any(mustart == 0))
        stop("some combinations of the response not realized")
    }
  }),
  linkinv = function(eta, extra = NULL) {
      u1  <- eta[, 1]
      u2  <- eta[, 2]
      u3  <- eta[, 3]
      u12 <- eta[, 4]
      u13 <- eta[, 5]
      u23 <- eta[, 6]
      denom <- 1 + exp(u1) + exp(u2) + exp(u3) + exp(u1 + u2 + u12) +
               exp(u1 + u3 + u13) + exp(u2 + u3 + u23) +
               exp(u1 + u2 + u3 + u12 + u13 + u23)


    cbind("000" = 1,
          "001" = exp(u3),
          "010" = exp(u2),
          "011" = exp(u2+u3+u23),
          "100" = exp(u1),
          "101" = exp(u1+u3+u13),
          "110" = exp(u1+u2+u12),
          "111" = exp(u1+u2+u3+u12+u13+u23)) / denom
  },
  last = expression({
    misc$link <- rep_len("identitylink", M)
    names(misc$link) <- predictors.names
    misc$earg <- list(u1  = list(theta = NULL),
                      u2  = list(theta = NULL),
                      u3  = list(theta = NULL),
                      u12 = list(theta = NULL),
                      u13 = list(theta = NULL),
                      u23 = list(theta = NULL))


  }),
  linkfun = function(mu, extra = NULL)  {
    u0  <- log(mu[, 1])
    u3  <- log(mu[, 2]) - u0
    u2  <- log(mu[, 3]) - u0
    u23 <- log(mu[, 4]) - u0 - u2 - u3
    u1  <- log(mu[, 5]) - u0
    u13 <- log(mu[, 6]) - u0 - u1 - u3
    u12 <- log(mu[, 7]) - u0 - u1 - u2
    cbind(u1, u2, u3, u12, u13, u23)
  },
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    u1  <- eta[, 1]
    u2  <- eta[, 2]
    u3  <- eta[, 3]
    u12 <- eta[, 4]
    u13 <- eta[, 5]
    u23 <- eta[, 6]
    denom <- 1 + exp(u1) + exp(u2) + exp(u3) + exp(u1 + u2 + u12) +
             exp(u1 + u3 + u13) + exp(u2 + u3 + u23) +
             exp(u1 + u2 + u3 + u12 + u13 + u23)

    u0 <- -log(denom)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * (u0 + u1*y[, 1] + u2*y[, 2] + u3*y[, 3] +
                u12*y[, 1]*y[, 2] +
                u13*y[, 1]*y[, 3] + u23*y[, 2]*y[, 3])
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("loglinb3.orig"),
  validparams = function(eta, y, extra = NULL) {
    okay1 <- all(is.finite(eta))
    okay1
  },
  deriv = expression({
    u1  <- eta[, 1]
    u2  <- eta[, 2]
    u3  <- eta[, 3]
    u12 <- eta[, 4]
    u13 <- eta[, 5]
    u23 <- eta[, 6]
    denom <- 1 + exp(u1) + exp(u2) + exp(u3) + exp(u1 + u2 + u12) +
             exp(u1 + u3 + u13) + exp(u2 + u3 + u23) +
             exp(u1 + u2 + u3 + u12 + u13 + u23)



    allterms <- exp(u1+u2+u3+u12+u13+u23)
    A1 <- exp(u1) + exp(u1 + u2 + u12) + exp(u1 + u3 + u13) +
          allterms
    A2 <- exp(u2) + exp(u1 + u2 + u12) + exp(u2 + u3 + u23) +
          allterms
    A3 <- exp(u3) + exp(u3 + u2 + u23) + exp(u1 + u3 + u13) +
          allterms
    A12 <- exp(u1 + u2 + u12) + allterms
    A13 <- exp(u1 + u3 + u13) + allterms
    A23 <- exp(u2 + u3 + u23) + allterms


    c(w) * cbind(-A1/denom + y[, 1],
                 -A2/denom + y[, 2],
                 -A3/denom + y[, 3],
                 -A12/denom + y[, 1]*y[, 2],
                 -A13/denom + y[, 1]*y[, 3],
                 -A23/denom + y[, 2]*y[, 3])
  }),
  weight = expression({
    u0 <- -log(denom)
    dA2.du1 <- exp(u1 + u2 + u12) + allterms
    dA3.du1 <- exp(u1 + u3 + u13) + allterms
    dA3.du2 <- exp(u2 + u3 + u23) + allterms

    wz <- matrix(NA_real_, n, dimm(6))
    expu0 <- exp(u0)

    wz[,iam(1,1,M)] <- A1 * (1 - expu0 * A1)
    wz[,iam(2,2,M)] <- A2 * (1 - expu0 * A2)
    wz[,iam(3,3,M)] <- A3 * (1 - expu0 * A3)
    wz[,iam(1,2,M)] <- (dA2.du1 - expu0 * A1 * A2)
    wz[,iam(1,3,M)] <- (dA3.du1 - expu0 * A1 * A3)
    wz[,iam(2,3,M)] <- (dA3.du2 - expu0 * A2 * A3)
    wz[,iam(4,4,M)] <- A12 * (1 - expu0 * A12)
    wz[,iam(5,5,M)] <- A13 * (1 - expu0 * A13)
    wz[,iam(6,6,M)] <- A23 * (1 - expu0 * A23)
    wz[,iam(4,6,M)] <- (allterms - expu0 * A12 * A23)
    wz[,iam(5,6,M)] <- (allterms - expu0 * A12 * A23)
    wz[,iam(4,5,M)] <- (allterms - expu0 * A12 * A13)
    wz[,iam(1,4,M)] <- A12 * (1 - expu0 * A1)
    wz[,iam(1,5,M)] <- A13 * (1 - expu0 * A1)
    wz[,iam(1,6,M)] <- (allterms - expu0 * A1 * A23)
    wz[,iam(2,4,M)] <- A12 * (1 - expu0 * A2)
    wz[,iam(2,5,M)] <- (allterms - expu0 * A2 * A13)
    wz[,iam(2,6,M)] <- A23 * (1 - expu0 * A2)
    wz[,iam(3,4,M)] <- (allterms - expu0 * A3 * A12)
    wz[,iam(3,5,M)] <- A13 * (1 - expu0 * A3)
    wz[,iam(3,6,M)] <- A23 * (1 - expu0 * A3)
    wz <- expu0 * wz
    c(w) * wz
  }))
}  # loglinb3.orig





 loglinb3 <-
    function(exchangeable = FALSE,
             zero = c("u12", "u13", "u23", if
                      (u123.arg) "u123" else NULL),
             u123.arg = FALSE) {


  if (!isFALSE(u123.arg) && !isTRUE(u123.arg))
    warning("'u123.arg' not a single logical")
  if (!isFALSE(exchangeable) && !isTRUE(exchangeable))
    warning("'exchangeable' not a single logical")


  new("vglmff",
  blurb = c("Loglinear model: 3 binary responses\n\n",
            "Links:    ",
            "Identity: u1, u2, u3, u12, u13, u23",
            if (u123.arg) ", u123" else NULL,
            "\n"),
  constraints = eval(substitute(expression({
    M <- ifelse( .u123.arg , 7, 6)
    cm.int.default <- diag(M)  # cm.intercept.default 

    use.mat <- matrix(c(1,1,1,0,0,0, 0,0,0,1,1,1),
                      6, 2)
    if (.u123.arg )
        use.mat <- rbind(cbind(use.mat, 0),
                         c(0, 0, 1))

    constraints <-
        cm.VGAM(use.mat, x = x,
                bool = .exchangeable ,
                constraints = constraints,
                apply.int = TRUE,
                cm.default           = cm.int.default,
                cm.intercept.default = cm.int.default)
    constraints <-
    cm.zero.VGAM(constraints, x = x, .zero ,
                 predictors.names = predictors.names,
                 quiet = TRUE, M = M, M1 = M)
  }),
  list( .exchangeable = exchangeable,
        .u123.arg = u123.arg,
        .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 6 + ( .u123.arg ),   # 7 or 6
         Q1 = 3,  # ncol(depvar(object))
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("u1", "u2", "u3",
                              "u12", "u13", "u23",
           if ( .u123.arg ) "u123" else NULL),
         zero = .zero )
  },
  list( .zero = zero,
        .u123.arg = u123.arg ))),


  initialize = eval(substitute(expression({
    predictors.names <- c("u1", "u2", "u3",
               "u12", "u13", "u23",
               if ( .u123.arg ) "u123" else NULL)
    Q1 <- 3

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = Q1,
              out.wy = TRUE,
              colsyperw = Q1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    if (ncol(y) != Q1)
      stop("ncol(y) must be = ", Q1)


    if (length(mustart) + length(etastart) == 0) {
      mustart <- matrix(NA_real_, nrow(y), 2^3)
      mustart[,1] <- weighted.mean((1-y[,1])*(1-y[,2])*(1-y[,3]), w)
      mustart[,2] <- weighted.mean((1-y[,1])*(1-y[,2])*   y[,3] , w)
      mustart[,3] <- weighted.mean((1-y[,1])*   y[,2] *(1-y[,3]), w)
      mustart[,4] <- weighted.mean((1-y[,1])*   y[,2] *   y[,3] , w)
      mustart[,5] <- weighted.mean(   y[,1] *(1-y[,2])*(1-y[,3]), w)
      mustart[,6] <- weighted.mean(   y[,1] *(1-y[,2])*   y[,3] , w)
      mustart[,7] <- weighted.mean(   y[,1] *   y[,2] *(1-y[,3]), w)
      mustart[,8] <- weighted.mean(   y[,1] *   y[,2] *   y[,3] , w)
      if (any(mustart == 0))
        stop("some combinations of the response not realized")
    }
  }), list( .u123.arg = u123.arg )) ),
  
  linkinv = eval(substitute(
      function(eta, extra = NULL) {
      u1   <- eta[, 1]
      u2   <- eta[, 2]
      u3   <- eta[, 3]
      u12  <- eta[, 4]
      u13  <- eta[, 5]
      u23  <- eta[, 6]
      u123 <- if ( .u123.arg ) eta[, 7] else 0
      denom <- 1 + exp(u1) + exp(u2) + exp(u3) +
         exp(u1 + u2 + u12) +
         exp(u1 + u3 + u13) +
         exp(u2 + u3 + u23) +
         exp(u1 + u2 + u3 + u12 + u13 + u23 + u123)

    cbind("000" = 1,
          "001" = exp(u3),
          "010" = exp(u2),
          "011" = exp(u2+u3+u23),
          "100" = exp(u1),
          "101" = exp(u1+u3+u13),
          "110" = exp(u1+u2+u12),
          "111" = exp(u1+u2+u3+u12+u13+u23+
                      u123)) / denom
  }, list( .u123.arg = u123.arg ))),
  last = eval(substitute(expression({
    misc$link <- rep_len("identitylink", M)
    names(misc$link) <- predictors.names
    misc$earg <- list(u1  = list(theta = NULL),
                      u2  = list(theta = NULL),
                      u3  = list(theta = NULL),
                      u12 = list(theta = NULL),
                      u13 = list(theta = NULL),
                      u23 = list(theta = NULL))
    if ( .u123.arg ) 
      misc$earg <- c(misc$earg,
                     list(u123  = list(theta = NULL)))
  }), list( .u123.arg = u123.arg ))),
  linkfun = eval(substitute(
      function(mu, extra = NULL)  {
    u0   <- log(mu[, 1])
    u3   <- log(mu[, 2]) - u0
    u2   <- log(mu[, 3]) - u0
    u23  <- log(mu[, 4]) - u0 - u2 - u3
    u1   <- log(mu[, 5]) - u0
    u13  <- log(mu[, 6]) - u0 - u1 - u3
    u12  <- log(mu[, 7]) - u0 - u1 - u2
    u123 <- if ( .u123.arg ) 
            log(mu[, 8]) - u0 - u1 - u2 - u3 -
            u12 - u13 - u23 else NULL
    cbind(u1, u2, u3, u12, u13, u23, u123)
  }, list( .u123.arg = u123.arg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    u1   <- eta[, 1]
    u2   <- eta[, 2]
    u3   <- eta[, 3]
    u12  <- eta[, 4]
    u13  <- eta[, 5]
    u23  <- eta[, 6]
    u123 <- if ( .u123.arg ) eta[, 7] else 0
    denom <- 1 + exp(u1) + exp(u2) + exp(u3) +
         exp(u1 + u2 + u12) +
         exp(u1 + u3 + u13) +
         exp(u2 + u3 + u23) +
         exp(u1 + u2 + u3 + u12 + u13 + u23 + u123)

    u0 <- -log(denom)
    if (residuals) {
      stop("loglikelihood residuals not implemented")
    } else {
      ll.elts <-
          c(w) * (u0 + u1 * y[, 1] + u2 * y[, 2] +
                       u3 * y[, 3] +
                u12  * y[, 1] * y[, 2] +
                u13  * y[, 1] * y[, 3] +
                u23  * y[, 2] * y[, 3] +
                u123 * y[, 1] * y[, 2] * y[, 3])
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .u123.arg = u123.arg ))),
  vfamily = c("loglinb3"),
 validparams = eval(substitute(
     function(eta, y, extra = NULL) {
    okay1 <- all(is.finite(eta))
    okay1
  }, list( .u123.arg = u123.arg ))),
  deriv = eval(substitute(expression({
  u1   <- eta[, 1]
  u2   <- eta[, 2]
  u3   <- eta[, 3]
  u12  <- eta[, 4]
  u13  <- eta[, 5]
  u23  <- eta[, 6]
  u123 <- if ( .u123.arg ) eta[, 7] else 0
  esigma <- exp(u1 + u2 + u3 + u12 + u13 + u23 + u123)
  denom <- 1 + exp(u1) + exp(u2) + exp(u3) +
    exp(u1 + u2 + u12) +
    exp(u1 + u3 + u13) + exp(u2 + u3 + u23) + esigma
  A1 <- exp(u1) + exp(u1 + u2 + u12) +
      exp(u1 + u3 + u13) + esigma
  A2 <- exp(u2) + exp(u1 + u2 + u12) +
      exp(u2 + u3 + u23) + esigma
  A3 <- exp(u3) + exp(u3 + u2 + u23) +
      exp(u1 + u3 + u13) + esigma
  A12 <- exp(u1 + u2 + u12) + esigma
  A13 <- exp(u1 + u3 + u13) + esigma
  A23 <- exp(u2 + u3 + u23) + esigma
  c(w) *
  cbind(-A1 / denom + y[, 1],
        -A2 / denom + y[, 2],
        -A3 / denom + y[, 3],
        -A12 / denom + y[, 1] * y[, 2],
        -A13 / denom + y[, 1] * y[, 3],
        -A23 / denom + y[, 2] * y[, 3],
        if ( .u123.arg ) -esigma / denom +
         y[, 1] * y[, 2] * y[, 3] else NULL)
   }), list( .u123.arg = u123.arg ))),
  weight = eval(substitute(expression({
  wz <- matrix(NA_real_, n, dimm(M))
  wz[, iam(1, 1, M)] <- A1 * (denom - A1)
  wz[, iam(2, 2, M)] <- A2 * (denom - A2)
  wz[, iam(3, 3, M)] <- A3 * (denom - A3)
  wz[, iam(4, 4, M)] <- A12 * (denom - A12)
  wz[, iam(5, 5, M)] <- A13 * (denom - A13)
  wz[, iam(6, 6, M)] <- A23 * (denom - A23)
  wz[, iam(1, 2, M)] <- denom * A12 - A1 * A2
  wz[, iam(1, 3, M)] <- denom * A13 - A1 * A3
  wz[, iam(2, 3, M)] <- denom * A23 - A2 * A3
  wz[, iam(1, 4, M)] <- A12 * (denom - A1)
  wz[, iam(2, 4, M)] <- A12 * (denom - A2)
  wz[, iam(3, 4, M)] <- denom * esigma - A12 * A3
  wz[, iam(1, 5, M)] <- A13 * (denom - A1)
  wz[, iam(2, 5, M)] <- denom * esigma - A2 * A13
  wz[, iam(3, 5, M)] <- A13 * (denom - A3)
  wz[, iam(4, 5, M)] <- denom * esigma - A12 * A13
  wz[, iam(1, 6, M)] <- denom * esigma - A1 * A23
  wz[, iam(2, 6, M)] <- A23 * (denom - A2)
  wz[, iam(3, 6, M)] <- A23 * (denom - A3)
  wz[, iam(4, 6, M)] <- denom * esigma - A12 * A23
  wz[, iam(5, 6, M)] <- denom * esigma - A13 * A23
  if ( .u123.arg ) {
    wz[, iam(1, 7, M)] <- esigma * (denom - A1)
    wz[, iam(2, 7, M)] <- esigma * (denom - A2)
    wz[, iam(3, 7, M)] <- esigma * (denom - A3)
    wz[, iam(4, 7, M)] <- esigma * (denom - A12)
    wz[, iam(5, 7, M)] <- esigma * (denom - A13)
    wz[, iam(6, 7, M)] <- esigma * (denom - A23)
    wz[, iam(7, 7, M)] <- esigma * (denom - esigma)
  }
  wz <- wz / denom^2
  c(w) * wz
  }), list( .u123.arg = u123.arg ))))
}  # loglinb3






 loglinb4 <-
   function(order4 = 4,  # May be in 2:4, T, F (4-1)
            zero = c("u12", "u13", "u14",
                     "u23", "u24", "u34",
                     if (order4 > 2) c("u123", "u124",
                     "u134", "u234") else NULL,
            if (order4 > 3) "u1234" else NULL),
            exchangeable = FALSE) {

  if (is.numeric(order4))
    stopifnot(length(order4) == 1,
              order4 %in% 2:4) else
  if (!isFALSE(order4) && !isTRUE(order4))
    stop("'order4' not a single logical")
  if (is.logical(order4))
    order4 <- if (order4) 4 else 3

  if (!isFALSE(exchangeable) && !isTRUE(exchangeable))
    stop("'exchangeable' not a single logical")

  par.names <- c("u1", "u2", "u3", "u4",
    "u12", "u13", "u14", "u23", "u24", "u34",
    if (order4 >= 3)
        c("u123", "u124", "u134", "u234") else NULL,
    if (order4 == 4) "u1234" else NULL)
  M1 <- length(par.names)

  new("vglmff",
  blurb = c("Loglinear model: 4 binary responses\n\n",
            "Links:    ",
            "Identity: u1, u2, u3, u4, u12, u13,\n",
            "          u14, u23, u24, u34",
            if (order4 >= 3)
              ", u123, u124, u134, u234" else NULL,
            if (order4 == 4) ", u1234" else NULL,
            "\n"),
  constraints = eval(substitute(expression({
    cm.int.default <- diag(M)  # cm.intercept.default 
    use.mat <-
      matrix(c(rep(1, 4), rep(0, 6),
               rep(0, 4), rep(1, 6)), 10, 2)
    if( .order4 >= 3)
      use.mat <- rbind(cbind(use.mat, 0),
         cbind(rep(0, 4), rep(0, 4), rep(1, 4)))
    if (.order4 == 4)
        use.mat <- rbind(cbind(use.mat, 0),
                         c(0, 0, 0, 1))

    constraints <-
        cm.VGAM(use.mat, x = x,
                bool = .exchangeable ,
                constraints = constraints,
                apply.int = TRUE,
                cm.default           = cm.int.default,
                cm.intercept.default = cm.int.default)
    constraints <-
    cm.zero.VGAM(constraints, x = x, .zero ,
                 predictors.names = predictors.names,
                 quiet = TRUE, M = M, M1 = M)
  }),
  list( .exchangeable = exchangeable,
        .order4 = order4,
        .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = .M1 ,   # 6 + ( .order4 ),   # 7 or 6
         Q1 = 4,  # ncol(depvar(object))
         expected = TRUE,
         multipleResponses = FALSE, M = .M1 ,
         parameters.names = .par.names , 
         zero = .zero )
  },
  list( .zero = zero, .M1 = M1,
        .par.names = par.names,
        .order4 = order4 ))),


  initialize = eval(substitute(expression({
    predictors.names <- ( .par.names )
    Q1 <- 4
    M <- M1 <- ( .M1 )

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = Q1,
              out.wy = TRUE,
              colsyperw = Q1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    if (ncol(y) != Q1)
      stop("ncol(y) must be = ", Q1)

    if (length(mustart) + length(etastart) == 0) {
      mustart <- matrix(NA_real_, nrow(y), 2^4)
      wmf <- weighted.mean
      mustart <- cbind(
      wmf((1-y[,1])*(1-y[,2])*(1-y[,3])*(1-y[,4]), w),
      wmf((1-y[,1])*(1-y[,2])*(1-y[,3])*   y[,4] , w),
      wmf((1-y[,1])*(1-y[,2])*   y[,3] *(1-y[,4]), w),
      wmf((1-y[,1])*(1-y[,2])*   y[,3] *   y[,4] , w),
      wmf((1-y[,1])*   y[,2] *(1-y[,3])*(1-y[,4]), w),
      wmf((1-y[,1])*   y[,2] *(1-y[,3])*   y[,4] , w),
      wmf((1-y[,1])*   y[,2] *   y[,3] *(1-y[,4]), w),
      wmf((1-y[,1])*   y[,2] *   y[,3] *   y[,4] , w),
      wmf(   y[,1] *(1-y[,2])*(1-y[,3])*(1-y[,4]), w),
      wmf(   y[,1] *(1-y[,2])*(1-y[,3])*   y[,4] , w),
      wmf(   y[,1] *(1-y[,2])*   y[,3] *(1-y[,4]), w),
      wmf(   y[,1] *(1-y[,2])*   y[,3] *   y[,4] , w),
      wmf(   y[,1] *   y[,2] *(1-y[,3])*(1-y[,4]), w),
      wmf(   y[,1] *   y[,2] *(1-y[,3])*   y[,4] , w),
      wmf(   y[,1] *   y[,2] *   y[,3] *(1-y[,4]), w),
      wmf(   y[,1] *   y[,2] *   y[,3] *   y[,4] , w))
      if (any(mustart == 0))
          stop("some combinations of the response ",
               "not realized")
      mustart <- matrix(mustart, n, 2^4, byrow = TRUE)
    }
  }),
  list( .order4 = order4, .M1 = M1,
        .par.names = par.names )) ),
  
  linkinv = eval(substitute(
      function(eta, extra = NULL) {
      u1    <- eta[,  1]
      u2    <- eta[,  2]
      u3    <- eta[,  3]
      u4    <- eta[,  4]
      u12   <- eta[,  5]
      u13   <- eta[,  6]
      u14   <- eta[,  7]
      u23   <- eta[,  8]
      u24   <- eta[,  9]
      u34   <- eta[, 10]
      u123  <- if ( .order4 >= 3) eta[, 11] else 0
      u124  <- if ( .order4 >= 3) eta[, 12] else 0
      u134  <- if ( .order4 >= 3) eta[, 13] else 0
      u234  <- if ( .order4 >= 3) eta[, 14] else 0
      u1234 <- if ( .order4 == 4) eta[, 15] else 0

  esigma <- exp(u1 + u2 + u3 + u4 +
                u12 + u13 + u14 + u23 + u24 + u34 +
                u123 + u124 + u134 + u234 + u1234)
  denom <- 1 + exp(u1) + exp(u2) + exp(u3) +
      exp(u4) + 
      exp(u1 + u2 + u12) + exp(u1 + u3 + u13) +
      exp(u1 + u4 + u14) + exp(u2 + u3 + u23) +
      exp(u2 + u4 + u24) + exp(u3 + u4 + u34) + 
      exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
      exp(u1 + u2 + u4 + u12 + u14 + u24 + u124) +
      exp(u1 + u3 + u4 + u13 + u14 + u34 + u134) +
      exp(u2 + u3 + u4 + u23 + u24 + u34 + u234) +
      esigma
      if ( .order4 >= 3) denom <- denom +
         exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
         exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
         exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
         exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
         exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
         exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
         exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
         exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
    if ( .order4 == 4) denom <- denom +
         exp(u1 + u2 + u3 + u4 + u12 + u13 +
             u14 + u23 + u24 + u34 + u1234)

    cbind("0000" = 1,
          "0001" = exp(u4),
          "0010" = exp(u3),
          "0011" = exp(u3+u4+u34),
          "0100" = exp(u2),
          "0101" = exp(u2+u4+u24),
          "0110" = exp(u2+u3+u23),
          "0111" = exp(u2+u3+u4+u23+u24+u34+u234),
          "1000" = exp(u1),
          "1001" = exp(u1+u4+u14),
          "1010" = exp(u1+u3+u13),
          "1011" = exp(u1+u3+u4+u13+u14+u34+u134),
          "1100" = exp(u1+u2+u12),
          "1101" = exp(u1+u2+u4+u12+u14+u24+u124),
          "1110" = exp(u1+u2+u3+u12+u13+u23+u123),
          "1111" = exp(u1+u2+u3+u4+
                       u12+u13+u14+u23+u24+u34+
                       u123+u124+u134+u234+
                       u1234)) / denom
  }, list( .order4 = order4 ))),
  last = eval(substitute(expression({
    misc$link <- rep_len("identitylink", M)
    names(misc$link) <- predictors.names
    misc$earg <- list(u1  = list(theta = NULL),
                      u2  = list(theta = NULL),
                      u3  = list(theta = NULL),
                      u4  = list(theta = NULL),
                      u12 = list(theta = NULL),
                      u13 = list(theta = NULL),
                      u14 = list(theta = NULL),
                      u23 = list(theta = NULL),
                      u24 = list(theta = NULL),
                      u34 = list(theta = NULL))
    if ( .order4 >= 3) 
      misc$earg <- c(misc$earg,
                     list(u123 = list(theta = NULL)),
                     list(u124 = list(theta = NULL)),
                     list(u134 = list(theta = NULL)),
                     list(u234 = list(theta = NULL)))
    if ( .order4 == 4) 
      misc$earg <- c(misc$earg,
                     list(u1234 = list(theta = NULL)))
  }), list( .order4 = order4 ))),
  linkfun = eval(substitute(
      function(mu, extra = NULL)  {
    u0    <- log(mu[,  1])
    u4    <- log(mu[,  2])-u0
    u3    <- log(mu[,  3])-u0
    u34   <- log(mu[,  4])-u0-u3-u4
    u2    <- log(mu[,  5])-u0
    u24   <- log(mu[,  6])-u0-u2-u4
    u23   <- log(mu[,  7])-u0-u2-u3
    u1    <- log(mu[,  9])-u0
    u14   <- log(mu[, 10])-u0-u1-u4
    u13   <- log(mu[, 11])-u0-u1-u3
    u12   <- log(mu[, 13])-u0-u1-u2
    u123  <- u124 <- u134 <- u234 <- NULL
    if ( .order4 >= 3) {
      u234  <- log(mu[,  8])-u0-u2-u3-u4-u23-u24-u34
      u134  <- log(mu[, 12])-u0-u1-u3-u4-u13-u14-u34
      u124  <- log(mu[, 14])-u0-u1-u2-u4-u12-u14-u24
      u123  <- log(mu[, 15])-u0-u1-u2-u3-u12-u13-u23
    }
    u1234 <- if ( .order4 == 4) 
      log(mu[, 16])-u0-u1-u2-u3-u4-
      u12-u13-u14-u23-u24-u34 else NULL
    ans <-
    cbind(u1, u2, u3, u4,
          u12, u13, u14, u23, u24, u34,
          u123, u124, u134, u234, u1234)
    ans
  }, list( .order4 = order4 ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
      u1    <- eta[,  1]
      u2    <- eta[,  2]
      u3    <- eta[,  3]
      u4    <- eta[,  4]
      u12   <- eta[,  5]
      u13   <- eta[,  6]
      u14   <- eta[,  7]
      u23   <- eta[,  8]
      u24   <- eta[,  9]
      u34   <- eta[, 10]
      u123  <- if ( .order4 >= 3) eta[, 11] else 0
      u124  <- if ( .order4 >= 3) eta[, 12] else 0
      u134  <- if ( .order4 >= 3) eta[, 13] else 0
      u234  <- if ( .order4 >= 3) eta[, 14] else 0
      u1234 <- if ( .order4 == 4) eta[, 15] else 0
  esigma <- exp(u1 + u2 + u3 + u4 +
                u12 + u13 + u14 + u23 + u24 + u34 +
                u123 + u124 + u134 + u234 + u1234)
  denom <- 1 + exp(u1) + exp(u2) + exp(u3) +
      exp(u4) + 
      exp(u1 + u2 + u12) + exp(u1 + u3 + u13) +
      exp(u1 + u4 + u14) + exp(u2 + u3 + u23) +
      exp(u2 + u4 + u24) + exp(u3 + u4 + u34) + 
      exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
      exp(u1 + u2 + u4 + u12 + u14 + u24 + u124) +
      exp(u1 + u3 + u4 + u13 + u14 + u34 + u134) +
      exp(u2 + u3 + u4 + u23 + u24 + u34 + u234) +
      esigma
      
    u0 <- -log(denom)
    if (residuals) {
      stop("loglikelihood residuals not implemented")
    } else {
      ll.elts <-
          c(w) * (u0 + u1 * y[, 1] + u2 * y[, 2] +
                       u3 * y[, 3] + u4 * y[, 4] +
                u12  * y[, 1] * y[, 2] +
                u13  * y[, 1] * y[, 3] +
                u14  * y[, 1] * y[, 4] +
                u23  * y[, 2] * y[, 3] +
                u24  * y[, 2] * y[, 4] +
                u34  * y[, 3] * y[, 4] +
                u123 * y[, 1] * y[, 2] * y[, 3] +
                u124 * y[, 1] * y[, 2] * y[, 4] +
                u134 * y[, 1] * y[, 3] * y[, 4] +
                u234 * y[, 2] * y[, 3] * y[, 4] +
               u1234 * y[, 1] * y[, 2] * y[, 3] *
                       y[, 4])
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .order4 = order4 ))),
  vfamily = c("loglinb4"),
 validparams = eval(substitute(
     function(eta, y, extra = NULL) {
    okay1 <- all(is.finite(eta))
    okay1
  }, list( .order4 = order4 ))),
  deriv = eval(substitute(expression({
    u1    <- eta[,  1]
    u2    <- eta[,  2]
    u3    <- eta[,  3]
    u4    <- eta[,  4]
    u12   <- eta[,  5]
    u13   <- eta[,  6]
    u14   <- eta[,  7]
    u23   <- eta[,  8]
    u24   <- eta[,  9]
    u34   <- eta[, 10]
    u123  <- if ( .order4 >= 3) eta[, 11] else 0
    u124  <- if ( .order4 >= 3) eta[, 12] else 0
    u134  <- if ( .order4 >= 3) eta[, 13] else 0
    u234  <- if ( .order4 >= 3) eta[, 14] else 0
    u1234 <- if ( .order4 == 4) eta[, 15] else 0

  esigma <- exp(u1 + u2 + u3 + u4 +
                u12 + u13 + u14 + u23 + u24 + u34 +
                u123 + u124 + u134 + u234 + u1234)
  denom <- 1 + exp(u1) + exp(u2) + exp(u3) +
      exp(u4) + 
      exp(u1 + u2 + u12) + exp(u1 + u3 + u13) +
      exp(u1 + u4 + u14) + exp(u2 + u3 + u23) +
      exp(u2 + u4 + u24) + exp(u3 + u4 + u34) + 
      exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
      exp(u1 + u2 + u4 + u12 + u14 + u24 + u124) +
      exp(u1 + u3 + u4 + u13 + u14 + u34 + u134) +
      exp(u2 + u3 + u4 + u23 + u24 + u34 + u234) +
      esigma
  A1 <- exp(u1) + exp(u1 + u2 + u12) +
      exp(u1 + u3 + u13) + exp(u1 + u4 + u14) +
      exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
      exp(u1 + u2 + u4 + u12 + u14 + u24 + u124) +
      exp(u1 + u3 + u4 + u13 + u14 + u34 + u134) +
      esigma
  A2 <- exp(u2) + exp(u1 + u2 + u12) +
      exp(u2 + u3 + u23) + exp(u2 + u4 + u24) +
      exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
      exp(u1 + u2 + u4 + u12 + u14 + u24 + u124) +
      exp(u2 + u3 + u4 + u23 + u24 + u34 + u234) +
      esigma
  A3 <- exp(u3) + exp(u1 + u3 + u13) +
      exp(u2 + u3 + u23) + exp(u3 + u4 + u34) +
      exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
      exp(u1 + u3 + u4 + u13 + u14 + u34 + u134) +
      exp(u2 + u3 + u4 + u23 + u24 + u34 + u234) +
      esigma
  A4 <- exp(u4) + exp(u1 + u4 + u14) +
      exp(u2 + u4 + u24) + exp(u3 + u4 + u34) +
      exp(u1 + u2 + u4 + u12 + u14 + u24 + u124) +
      exp(u1 + u3 + u4 + u13 + u14 + u34 + u134) +
      exp(u2 + u3 + u4 + u23 + u24 + u34 + u234) +
      esigma
  A12 <- exp(u1 + u2 + u12) + 
      exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
      exp(u1 + u2 + u4 + u12 + u14 + u24 + u124) + 
    esigma
  A13 <- exp(u1 + u3 + u13) + 
      exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
      exp(u1 + u3 + u4 + u13 + u14 + u34 + u134) +
    esigma
  A14 <- exp(u1 + u4 + u14) + 
      exp(u1 + u2 + u4 + u12 + u14 + u24 + u124) +
      exp(u1 + u3 + u4 + u13 + u14 + u34 + u134) +
    esigma
  A23 <- exp(u2 + u3 + u23) + 
      exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
      exp(u2 + u3 + u4 + u23 + u24 + u34 + u234) +
    esigma
  A24 <- exp(u2 + u4 + u24) + 
      exp(u1 + u2 + u4 + u12 + u14 + u24 + u124) +
      exp(u2 + u3 + u4 + u23 + u24 + u34 + u234) +
    esigma
  A34 <- exp(u3 + u4 + u34) + 
      exp(u1 + u3 + u4 + u13 + u14 + u34 + u134) +
      exp(u2 + u3 + u4 + u23 + u24 + u34 + u234) +
      esigma
  A123 <- exp(u1 + u2 + u3 + u12 + u13 + u23 + u123) +
      esigma
  A124 <- exp(u1 + u2 + u4 + u12 + u14 + u24 + u124) +
      esigma
  A134 <- exp(u1 + u3 + u4 + u13 + u14 + u34 + u134) +
      esigma
  A234 <- exp(u2 + u3 + u4 + u23 + u24 + u34 + u234) +
      esigma
  ans <- 
  c(w) * cbind(  # numeric(0) act like NULL
    -A1/denom + y[, 1],
    -A2/denom + y[, 2],
    -A3/denom + y[, 3],
    -A4/denom + y[, 4],
    -A12/denom + y[, 1] * y[, 2],
    -A13/denom + y[, 1] * y[, 3],
    -A14/denom + y[, 1] * y[, 4], 
    -A23/denom + y[, 2] * y[, 3],
    -A24/denom + y[, 2] * y[, 4],
    -A34/denom + y[, 3] * y[, 4])  # Added
  if ( .order4 >= 3)
     ans <- cbind(ans, c(w) * cbind(   
    -A123/denom + y[, 1] * y[, 2] * y[, 3],
    -A124/denom + y[, 1] * y[, 2] * y[, 4],
    -A134/denom + y[, 1] * y[, 3] * y[, 4],
    -A234/denom + y[, 2] * y[, 3] * y[, 4]))
  if ( .order4 >= 4)
     ans <- cbind(ans, c(w) * cbind(-esigma / denom +
                  y[, 1] * y[, 2] * y[, 3] * y[, 4]))
  ans
  }), list( .order4 = order4))),
  weight = eval(substitute(expression({
  wz <- matrix(NA_real_, n, dimm(M))
  wz[, iam( 1,  1, M)] <- A1 * (denom - A1)
  wz[, iam( 2,  2, M)] <- A2 * (denom - A2)
  wz[, iam( 3,  3, M)] <- A3 * (denom - A3)
  wz[, iam( 4,  4, M)] <- A4 * (denom - A4)
  wz[, iam( 5,  5, M)] <- A12 * (denom - A12)
  wz[, iam( 6,  6, M)] <- A13 * (denom - A13)
  wz[, iam( 7,  7, M)] <- A14 * (denom - A14)
  wz[, iam( 8,  8, M)] <- A23 * (denom - A23)
  wz[, iam( 9,  9, M)] <- A24 * (denom - A24)
  wz[, iam(10, 10, M)] <- A34 * (denom - A34)
  if ( .order4 >= 3) {
  wz[, iam(11, 11, M)] <- A123 * (denom - A123)
  wz[, iam(12, 12, M)] <- A124 * (denom - A124)
  wz[, iam(13, 13, M)] <- A134 * (denom - A134)
  wz[, iam(14, 14, M)] <- A234 * (denom - A234)
  if ( .order4 >= 4)
    wz[, iam(15, 15, M)] <- esigma * (denom - esigma)
  }
  
  wz[, iam(1, 2, M)] <- denom * A12 - A1 * A2
  wz[, iam(1, 3, M)] <- denom * A13 - A1 * A3
  wz[, iam(2, 3, M)] <- denom * A23 - A2 * A3
  wz[, iam(1, 4, M)] <- denom * A14 - A1 * A4
  wz[, iam(2, 4, M)] <- denom * A24 - A2 * A4
  wz[, iam(3, 4, M)] <- denom * A34 - A3 * A4
  
  wz[, iam(1, 5, M)] <- A12 * (denom - A1)
  wz[, iam(2, 5, M)] <- A12 * (denom - A2)
  wz[, iam(3, 5, M)] <- denom * A123 - A12 * A3
  wz[, iam(4, 5, M)] <- denom * A124 - A12 * A4
  
  wz[, iam(1, 6, M)] <- A13 * (denom - A1)
  wz[, iam(2, 6, M)] <- denom * A123 - A13 * A2
  wz[, iam(3, 6, M)] <- A13 * (denom - A3)
  wz[, iam(4, 6, M)] <- denom * A134 - A13 * A4
  wz[, iam(5, 6, M)] <- denom * A123 - A13 * A12
  
  wz[, iam(1, 7, M)] <- A14 * (denom - A1)
  wz[, iam(2, 7, M)] <- denom * A124 - A14 * A2
  wz[, iam(3, 7, M)] <- denom * A134 - A14 * A3
  wz[, iam(4, 7, M)] <- A14 * (denom - A4)
  wz[, iam(5, 7, M)] <- denom * A124 - A14 * A12
  wz[, iam(6, 7, M)] <- denom * A134 - A14 * A13
  
  wz[, iam(1, 8, M)] <- denom * A123 - A23 * A1
  wz[, iam(2, 8, M)] <- A23 * (denom - A2)
  wz[, iam(3, 8, M)] <- A23 * (denom - A3)
  wz[, iam(4, 8, M)] <- denom * A234 - A23 * A4
  wz[, iam(5, 8, M)] <- denom * A123 - A23 * A12
  wz[, iam(6, 8, M)] <- denom * A123 - A23 * A13
  wz[, iam(7, 8, M)] <- denom * esigma - A23 * A14
  
  wz[, iam(1, 9, M)] <- denom * A124 - A24 * A1
  wz[, iam(2, 9, M)] <- A24 * (denom - A2)
  wz[, iam(3, 9, M)] <- denom * A234 - A24 * A3
  wz[, iam(4, 9, M)] <- A24 * (denom - A4)
  wz[, iam(5, 9, M)] <- denom * A124 - A24 * A12
  wz[, iam(6, 9, M)] <- denom * esigma - A24 * A13
  wz[, iam(7, 9, M)] <- denom * A124 - A24 * A14
  wz[, iam(8, 9, M)] <- denom * A234 - A24 * A23
  
  wz[, iam(1, 10, M)] <- denom * A134 - A34 * A1
  wz[, iam(2, 10, M)] <- denom * A234 - A34 * A2
  wz[, iam(3, 10, M)] <- A34 * (denom - A3)
  wz[, iam(4, 10, M)] <- A34 * (denom - A4)
  wz[, iam(5, 10, M)] <- denom * esigma - A34 * A12
  wz[, iam(6, 10, M)] <- denom * A134 - A34 * A13
  wz[, iam(7, 10, M)] <- denom * A134 - A34 * A14
  wz[, iam(8, 10, M)] <- denom * A234 - A34 * A23
  wz[, iam(9, 10, M)] <- denom * A234 - A34 * A24

  if ( .order4 >= 3) {
  wz[, iam( 1, 11, M)] <- A123 * (denom - A1)
  wz[, iam( 2, 11, M)] <- A123 * (denom - A2)
  wz[, iam( 3, 11, M)] <- A123 * (denom - A3)
  wz[, iam( 4, 11, M)] <- denom * esigma - A123 * A4
  wz[, iam( 5, 11, M)] <- A123 * (denom - A12)
  wz[, iam( 6, 11, M)] <- A123 * (denom - A13)
  wz[, iam( 7, 11, M)] <- denom * esigma - A123 * A14
  wz[, iam( 8, 11, M)] <- A123 * (denom - A23)
  wz[, iam( 9, 11, M)] <- denom * esigma - A123 * A24
  wz[, iam(10, 11, M)] <- denom * esigma - A123 * A34
  
  wz[, iam( 1, 12, M)] <- A124 * (denom - A1)
  wz[, iam( 2, 12, M)] <- A124 * (denom - A2)
  wz[, iam( 3, 12, M)] <- denom * esigma - A124 * A3
  wz[, iam( 4, 12, M)] <- A124 * (denom - A4)
  wz[, iam( 5, 12, M)] <- A124 * (denom - A12)
  wz[, iam( 6, 12, M)] <- denom * esigma - A124 * A13
  wz[, iam( 7, 12, M)] <- A124 * (denom - A14)
  wz[, iam( 8, 12, M)] <- denom * esigma - A124 * A23
  wz[, iam( 9, 12, M)] <- A124 * (denom - A24)
  wz[, iam(10, 12, M)] <- denom * esigma - A124 * A34
  wz[, iam(11, 12, M)] <- denom * esigma - A124 * A123
  
  wz[, iam( 1, 13, M)] <- A134 * (denom - A1)
  wz[, iam( 2, 13, M)] <- denom * esigma - A134 * A2
  wz[, iam( 3, 13, M)] <- A134 * (denom - A3)
  wz[, iam( 4, 13, M)] <- A134 * (denom - A4)
  wz[, iam( 5, 13, M)] <- denom * esigma - A134 * A12
  wz[, iam( 6, 13, M)] <- A134 * (denom - A13)
  wz[, iam( 7, 13, M)] <- A134 * (denom - A14)
  wz[, iam( 8, 13, M)] <- denom * esigma - A134 * A23
  wz[, iam( 9, 13, M)] <- denom * esigma - A134 * A24
  wz[, iam(10, 13, M)] <- A134 * (denom - A34)
  wz[, iam(11, 13, M)] <- denom * esigma - A134 * A123
  wz[, iam(12, 13, M)] <- denom * esigma - A134 * A124
  
  wz[, iam( 1, 14, M)] <- denom * esigma - A234 * A1
  wz[, iam( 2, 14, M)] <- A234 * (denom - A2)
  wz[, iam( 3, 14, M)] <- A234 * (denom - A3)
  wz[, iam( 4, 14, M)] <- A234 * (denom - A4)
  wz[, iam( 5, 14, M)] <- denom * esigma - A234 * A12
  wz[, iam( 6, 14, M)] <- denom * esigma - A234 * A13
  wz[, iam( 7, 14, M)] <- denom * esigma - A234 * A14
  wz[, iam( 8, 14, M)] <- A234 * (denom - A23)
  wz[, iam( 9, 14, M)] <- A234 * (denom - A24)
  wz[, iam(10, 14, M)] <- A234 * (denom - A34)
  wz[, iam(11, 14, M)] <- denom * esigma - A234*A123
  wz[, iam(12, 14, M)] <- denom * esigma - A234*A124
  wz[, iam(13, 14, M)] <- denom * esigma - A234*A134
  }

  if ( .order4 >= 4) {
  wz[, iam( 1, 15, M)] <- esigma * (denom - A1)
  wz[, iam( 2, 15, M)] <- esigma * (denom - A2)
  wz[, iam( 3, 15, M)] <- esigma * (denom - A3)
  wz[, iam( 4, 15, M)] <- esigma * (denom - A4)
  wz[, iam( 5, 15, M)] <- esigma * (denom - A12)
  wz[, iam( 6, 15, M)] <- esigma * (denom - A13)
  wz[, iam( 7, 15, M)] <- esigma * (denom - A14)
  wz[, iam( 8, 15, M)] <- esigma * (denom - A23)
  wz[, iam( 9, 15, M)] <- esigma * (denom - A24)
  wz[, iam(10, 15, M)] <- esigma * (denom - A34)
  wz[, iam(11, 15, M)] <- esigma * (denom - A123)
  wz[, iam(12, 15, M)] <- esigma * (denom - A124)
  wz[, iam(13, 15, M)] <- esigma * (denom - A134)
  wz[, iam(14, 15, M)] <- esigma * (denom - A234)
  }
 
  wz <- wz / denom^2
  ans <- c(w) * wz
  ans
  }), list( .order4 = order4 ))))
}  # loglinb4























# These functions are
# Copyright (C) 1998-2014 T.W. Yee, University of Auckland.
# All rights reserved.















 abbott <- function(link0 = "logit",
                    link1 = "logit",
                    iprob0 = NULL, iprob1 = NULL,
                    type.fitted = c("observed", "treatment", "control"),
                    mux.offdiagonal = 0.98,
                    zero = 1) {


  type.fitted <- match.arg(type.fitted,
                           c("observed", "treatment", "control"),
                           several.ok = TRUE)


  link0 <- as.list(substitute(link0))
  earg0 <- link2list(link0)
  link0 <- attr(earg0, "function.name")

  link1 <- as.list(substitute(link1))
  earg1 <- link2list(link1)
  link1 <- attr(earg1, "function.name")




  if (!is.Numeric(mux.offdiagonal, length.arg = 1) ||
      mux.offdiagonal >= 1 ||
      mux.offdiagonal <  0)
    stop("argument 'mux.offdiagonal' must be in the interval [0, 1)")


  new("vglmff",
  blurb = c("Abbott's model for binary responses\n",
            "mu = prob0 + (1 - prob0) * prob1\n",
            "where 'prob0' is the 'control' mortality and\n",
            "'prob1' is the 'treatment' mortality and\n",
            "'mu' is the 'observed' mortality\n\n",
            "Links: ",
            namesof("prob0", link0, earg = earg0), ",  ",
            namesof("prob1", link1, earg = earg1)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero # ,
           ))),

  initialize = eval(substitute(expression({
    eval(binomialff(link = .link0 )@initialize)  # w, y, mustart are assigned


    predictors.names <-
      c(namesof("prob0", .link0, earg = .earg0, short =  TRUE),
        namesof("prob1", .link1, earg = .earg1, short =  TRUE))


    if (is.null(etastart)) {
      prob0.init <- if (length( .iprob0 )) {
        rep( .iprob0, length.out = n)
      } else {
        mustart / 2
      }

      prob1.init <- if (length( .iprob1 )) {
        rep( .iprob1, length.out = n)
      } else {
        mustart / 2
      }


      mustart <- NULL


        etastart <-
         cbind(theta2eta(prob0.init, link = .link0 , earg = .earg0 ),
               theta2eta(prob1.init, link = .link1 , earg = .earg1 ))
    }
  }), list( .link0 = link0, .earg0 = earg0,
            .link1 = link1, .earg1 = earg1,
            .iprob0 = iprob0, .iprob1 = iprob1 ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob0 <- eta2theta(eta[, 1], .link0 , earg = .earg0 )
    prob1 <- eta2theta(eta[, 2], .link1 , earg = .earg1 )

    con.fv <- prob0
    trt.fv <- prob1
    obs.fv <- prob0 + (1 - prob0) * prob1



    ans <- cbind("observed"  = obs.fv,
                "treatment" = trt.fv,
                "control"   = con.fv)

                   
    ans[, .type.fitted , drop = FALSE]
  }, list( .link0 = link0, .earg0 = earg0,
           .link1 = link1, .earg1 = earg1,
           .type.fitted = type.fitted ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {


    prob0 <- eta2theta(eta[, 1], .link0 , earg = .earg0 )
    prob1 <- eta2theta(eta[, 2], .link1 , earg = .earg1 )
    mymu <- prob0 + (1 - prob0) * prob1


    if (residuals) {
      w * (y / mymu - (1 - y) / (1 - mymu))
    } else {
      ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                  y * w # Convert proportions to counts
      nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                  round(w)
      smallno <- 1.0e6 * .Machine$double.eps
      smallno <- sqrt(.Machine$double.eps)
      if (max(abs(ycounts - round(ycounts))) > smallno)
        warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
        dbinom(x = ycounts, size = nvec, prob = mymu, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link0 = link0, .earg0 = earg0,
           .link1 = link1, .earg1 = earg1 ))),


  last = eval(substitute(expression({
    misc$link <-    c(prob0 = .link0 , prob1 = .link1 )
    misc$earg <- list(prob0 = .earg0 , prob1 = .earg1 )

    misc$mux.offdiagonal <- .mux.offdiagonal
    misc$type.fitted <- .type.fitted
    misc$true.mu <- ( .type.fitted == "observed")


  }), list( .link0 = link0, .earg0 = earg0,
            .link1 = link1, .earg1 = earg1,
            .mux.offdiagonal = mux.offdiagonal,
            .type.fitted = type.fitted
          ))),
  vfamily = c("abbott", "vquantal"),
  deriv = eval(substitute(expression({
    prob0 <- eta2theta(eta[, 1], .link0, earg = .earg0 )
    prob1 <- eta2theta(eta[, 2], .link1, earg = .earg1 )
    dprob0.deta <- dtheta.deta(prob0, .link0 , earg = .earg0 )
    dprob1.deta <- dtheta.deta(prob1, .link1 , earg = .earg1 )


    mymu <- prob0 + (1 - prob0) * prob1


    dl.dmu <- y / mymu - (1 - y) / (1 - mymu)
    dmu.dprob0 <- 1 - prob1
    dmu.dprob1 <- 1 - prob0
    dl.dprob0 <- dl.dmu * dmu.dprob0 
    dl.dprob1 <- dl.dmu * dmu.dprob1 


    c(w) * cbind(dl.dprob0 * dprob0.deta,
                 dl.dprob1 * dprob1.deta)
  }), list( .link0 = link0, .earg0 = earg0,
            .link1 = link1, .earg1 = earg1 ))),
  weight = eval(substitute(expression({


    ned2l.dmu2 <- 1 / (mymu * (1-mymu))
    ned2l.dprob02     <- ned2l.dmu2 * dmu.dprob0^2
    ned2l.dprob12     <- ned2l.dmu2 * dmu.dprob1^2
    ned2l.dprob1prob2 <-              ( 1)  # seems sort of ok but slow cvgc
    ned2l.dprob1prob2 <-              ( 0)  # kill it
    ned2l.dprob1prob2 <- ned2l.dmu2 * ( 1)  # dont seem to work

    ned2l.dprob1prob2 <- ned2l.dmu2 * dmu.dprob1 * dmu.dprob0 *
                         .mux.offdiagonal

    od2l.dmu2 <- y / mymu^2 + (1 - y) / (1 - mymu)^2
    od2l.dprob02     <- od2l.dmu2 * dmu.dprob0^2
    od2l.dprob12     <- od2l.dmu2 * dmu.dprob1^2
    od2l.dprob1prob2 <- od2l.dmu2 * dmu.dprob1 * dmu.dprob0 + dl.dmu


    wz <- cbind(ned2l.dprob02 * dprob0.deta^2,
                ned2l.dprob12 * dprob1.deta^2,
                ned2l.dprob1prob2 * dprob1.deta * dprob0.deta)




    c(w) * wz
  }), list( .link0 = link0, .earg0 = earg0,
            .link1 = link1, .earg1 = earg1,
             .mux.offdiagonal = mux.offdiagonal ))))
}








if (FALSE)
 Abbott <- function(lprob1 = elogit(min = 0, max = 1),  # For now, that is
                   lprob0 = "logit",
                   iprob0 = NULL, iprob1 = NULL,
                   nointercept = 2,  # NULL,
                   zero = 1) {




 stop("does not work")


  lprob1 <- as.list(substitute(lprob1))
  eprob1 <- link2list(lprob1)
  lprob1 <- attr(eprob1, "function.name")

  lprob0 <- as.list(substitute(lprob0))
  eprob0 <- link2list(lprob0)
  lprob0 <- attr(eprob0, "function.name")



  new("vglmff",
  blurb = c("Abbott's model for binary response\n",
            "mu = prob0 + prob1\n",
            "where 'prob0' is the control mortality and\n",
            "'prob1' is the treatment mortality\n\n",
            "Links: ",
            namesof("prob0", lprob0, earg = eprob0), ",  ",
            namesof("prob1", lprob1, earg = eprob1)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
    constraints <- cm.nointercept.vgam(constraints, x, .nointercept, M)
  }), list( .zero = zero,
            .nointercept = nointercept ))),


  initialize = eval(substitute(expression({
 print("here1")

    eval(binomialff(link = .lprob1)@initialize)  # w, y, mustart are assigned

 print("here2")
 print("summary(mustart)")
 print( summary(mustart) )

    predictors.names <-
      c(namesof("prob0", .lprob0, earg = .eprob0, short =  TRUE),
        namesof("prob1", .lprob1, earg = .eprob1, short =  TRUE))


    if (is.null(etastart)) {
      prob0.init <- if (length( .iprob0 )) {
          rep( .iprob0, len = n)
        } else {
          mustart / 2
        }

      prob1.init <- if (length( .iprob1 )) {
          rep( .iprob1, len = n)
        } else {
          mustart * 1 / 4
        }


      mustart <- NULL


 print("prob0.init ")
 print( sort(prob0.init) )
 print("prob1.init ")
 print( sort(prob1.init) )


        eprob1 <- list(min = prob0.init, max = 1)
        etastart <-
         cbind(theta2eta(prob0.init, link = .lprob0 , earg = .eprob0 ),
               theta2eta(prob1.init, link = .lprob1 , earg =  eprob1 ))
 print("head(etastart)")
 print( head(etastart) )
    }
  }), list( .lprob1 = lprob1, .eprob1 = eprob1,
            .lprob0 = lprob0, .eprob0 = eprob0,
            .iprob0 = iprob0, .iprob1 = iprob1 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob0 <- eta2theta(eta[, 1], .lprob0 , earg = .eprob0)

    eprob1 <- list(min = prob0, max = 1)

    prob1 <- eta2theta(eta[, 2], .lprob1 , earg =  eprob1)
    prob0 + prob1
  }, list( .lprob1 = lprob1, .eprob1 = eprob1,
           .lprob0 = lprob0, .eprob0 = eprob0 ))),
  last = eval(substitute(expression({
    eprob1 <- list(min = prob0, max = 1)
    misc$link <-    c(prob0 = .lprob0, prob1 = .lprob1)

    misc$earg <- list(prob0 = .eprob0, prob1 =  eprob1)

    misc$nointercept = .nointercept
  }), list( .lprob1 = lprob1, .eprob1 = eprob1,
            .lprob0 = lprob0, .eprob0 = eprob0,
            .nointercept = nointercept ))),
  vfamily = c("Abbott", "vquantal"),
  deriv = eval(substitute(expression({
    prob0 <- eta2theta(eta[,1], .lprob0, earg = .eprob0)

    eprob1 <- list(min = prob0, max = 1)
    prob1 <- eta2theta(eta[,2], .lprob1, earg =  eprob1)
    dprob0.deta <- dtheta.deta(prob0, .lprob0 , earg = .eprob0 )
    dprob1.deta <- dtheta.deta(prob1, .lprob1 , earg =  eprob1 )

    dl.dmu <- y / mu - (1 - y) / (1 - mu)
    dmu.dprob0 <- 1 # - prob1
    dmu.dprob1 <- 1 # - prob0
    dl.dprob0 <- dl.dmu * dmu.dprob0 
    dl.dprob1 <- dl.dmu * dmu.dprob1 

    c(w) * cbind(dl.dmu * dmu.dprob0 * dprob0.deta,
                 dl.dmu * dmu.dprob1 * dprob1.deta)
  }), list( .lprob1 = lprob1, .eprob1 = eprob1,
            .lprob0 = lprob0, .eprob0 = eprob0 ))),
    weight = eval(substitute(expression({


    ned2l.dmu2 <- 1 / (mu * (1-mu))
    ned2l.dprob02 <- ned2l.dmu2 * dmu.dprob0^2
    ned2l.dprob12 <- ned2l.dmu2 * dmu.dprob1^2

    wz <- cbind(ned2l.dprob02 * dprob0.deta^2,
                ned2l.dprob12 * dprob1.deta^2)

 print("head(wz)")
 print( head(wz) )
    c(w) * wz
  }), list( .lprob1 = lprob1, .eprob1 = eprob1,
            .lprob0 = lprob0, .eprob0 = eprob0 ))))
}















abbott.EM.control <- function(maxit = 1000, ...) {
  list(maxit = maxit)
}


 abbott.EM <-
  function(link = "probit",
           b1.arg = 0, b2.arg = 0,
           imethod = 1, ilambda = 0.5,
           iprob = NULL) {


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.Numeric(b1.arg,  # length.arg = 1,
                  integer.valued = TRUE) ||
      b1.arg < 0)
    stop("argument 'b1.arg' must be a vector of non-negative integers")


  if (!is.Numeric(b2.arg,  # length.arg = 1,
                  integer.valued = TRUE) ||
      b2.arg < 0)
    stop("argument 'b2.arg' must be a vector of non-negative integers")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  zero <- NULL
  if (length(zero) &&
      !is.Numeric(zero, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'zero'")


  new("vglmff",
  blurb = c("Probit regression with nonzero background (EM algorithm)\n",
            "P[Y=1] = mu = prob0 + (1 - prob0) * linkinv(eta)\n\n",
            "Link:     ",
            namesof("pi", link, earg = earg), "\n",
            "Mean:     mu"),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 1
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         zero = .zero )
  }, list( .zero = zero ))),

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


    if (length(table(y)) != 2 || max(y) > 1)
      stop("response must be a vector of 0s and 1s only")


    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly
    extra$lambda <- matrix( .ilambda , n, M, byrow = TRUE)
    extra$orig.w <- w


    mynames1 <- paste("prob0", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
      namesof(mynames1, .link , earg = .earg , tag = FALSE)


    if (!length(etastart)) {
      prob.init <- if ( .imethod == 2)
                      1 / (1 + y + 1/16) else
                  if ( .imethod == 3)
                      1 / (1 + apply(y, 2, median) + 1/16) else
                      rnorm(n * M, mean = 0.5, sd = 0.01)  # Mean 0.5

      if (!is.matrix(prob.init))
        prob.init <- matrix(prob.init, n, M, byrow = TRUE)


      if (length( .iprob ))
        prob.init <- matrix( .iprob , n, M, byrow = TRUE)  # Mean 0.5


        etastart <- theta2eta(prob.init, .link , earg = .earg )  # Mean 0
    }
  }), list( .link = link, .earg = earg,
            .ilambda = ilambda,
            .imethod = imethod, .iprob = iprob ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob <- eta2theta(eta, .link , earg = .earg )
    mymu <- extra$lambda + (1 - extra$lambda) * prob  # Eqn (3)
    mymu
  }, list( .link = link, .earg = earg ))),

  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <- c(rep( .link , length = ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$M1 <- M1
    misc$multipleResponses <- TRUE
    misc$imethod <- .imethod
    misc$iprob <- .iprob
    misc$b1.arg <- .b1.arg
    misc$b2.arg <- .b2.arg

    extra$lambda <- extra$lambda[1, ]  # Now a vector
  }), list( .link = link, .earg = earg,
            .iprob = iprob,
            .b1.arg = b1.arg, .b2.arg = b2.arg,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    prob <- eta2theta(eta, .link , earg = .earg )
    mymu <- extra$lambda + (1 - extra$lambda) * prob  # Eqn (3)

   if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)
      ll.elts <- c(w) * dbinom(x = y, prob = mymu,
                               size = nvec, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("abbott.EM"),
  deriv = eval(substitute(expression({
    prob <- eta2theta(eta, .link , earg = .earg )

    mymu <- extra$lambda + (1 - extra$lambda) * prob  # Eqn (3)

    wz <- cbind((1 - extra$lambda) * prob / mymu)  # Eqn (4)

    Deriv1 <- ifelse(y == 0, -dnorm(eta) / pnorm(eta, lower.tail = FALSE),
                              dnorm(eta) / pnorm(eta))

    c(w) * wz * Deriv1
  }), list( .link = link, .earg = earg ))),

  weight = eval(substitute(expression({

    extra$lambda <-
      matrix((colSums((1 - wz) * y) + .b1.arg ) / (n + .b1.arg + .b2.arg ),
             n, M, byrow = TRUE)  # Between eqns (6),(7)




    c(w) * wz
  }), list( .link = link, .earg = earg,
            .b1.arg = b1.arg, .b2.arg = b2.arg ))))
}






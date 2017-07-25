# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.

















calibrate.qrrvglm.control <-
  function(object,
           trace = FALSE,  # passed into optim()
           Method.optim = "BFGS",  # passed into optim(method = Method)
           gridSize = ifelse(Rank == 1, 9, 5),
         varI.latvar = FALSE, ...) {

  Rank <- object@control$Rank
  eq.tolerances <- object@control$eq.tolerances
  if (!is.Numeric(gridSize, positive = TRUE,
                  integer.valued = TRUE, length.arg = 1))
    stop("bad input for argument 'gridSize'")
  if (gridSize < 2)
    stop("'gridSize' must be >= 2")

  list(
       trace = as.numeric(trace)[1],
       Method.optim = Method.optim,
       gridSize = gridSize,
       varI.latvar = as.logical(varI.latvar)[1])
}  # calibrate.qrrvglm.control










calibrate.qrrvglm <-
  function(object,
           newdata = NULL,
           type = c("latvar", "predictors", "response", "vcov", "all3or4"),
           initial.vals = NULL, ...) {


  Quadratic <- if (is.logical(object@control$Quadratic))
               object@control$Quadratic else FALSE  # T if CQO, F if CAO


  newdata.orig <- newdata
  if (!length(newdata)) {
    newdata <- data.frame(depvar(object))
  }


  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("latvar", "predictors",
                            "response", "vcov", "all3or4"))[1]

  if (!Quadratic && type == "vcov")
    stop("cannot have 'type=\"vcov\"' when object is ",
         "a \"rrvgam\" object")



  if (!all(weights(object, type = "prior") == 1))
    warning("not all the prior weights of 'object' are 1; assuming ",
            "they all are here")



  if (!is.data.frame(newdata))
    newdata <- data.frame(newdata)
  if (!all(object@misc$ynames %in% colnames(newdata)))
    stop("need the following colnames in 'newdata': ",
         paste(object@misc$ynames, collapse = ", "))
  newdata <- newdata[, object@misc$ynames, drop = FALSE]



  if (!is.matrix(newdata))
    newdata <- as.matrix(newdata)


  nn <- nrow(newdata)




  obfunct <- slot(object@family, object@misc$criterion)
  minimize.obfunct <-
    if (Quadratic) object@control$min.criterion else
    TRUE  # Logical; TRUE for CAO objects because deviance is minimized
  if (!is.logical(minimize.obfunct))
    stop("object@control$min.criterion is not a logical")
  optim.control <- calibrate.qrrvglm.control(object = object, ...)

  use.optim.control <- optim.control
  use.optim.control$Method.optim <-
  use.optim.control$gridSize <-
  use.optim.control$varI.latvar <- NULL


  if ((Rank <- object@control$Rank) > 2)
    stop("currently can only handle Rank = 1 and 2")
  Coefobject <- if (Quadratic) {
    Coef(object, varI.latvar = optim.control$varI.latvar)
  } else {
    Coef(object)
  }



  if (!length(initial.vals)) {
    L <- apply(latvar(object), 2, min)
    U <- apply(latvar(object), 2, max)
    initial.vals <- if (Rank == 1)
        cbind(seq(L, U, length = optim.control$gridSize)) else
        expand.grid(seq(L[1], U[1], length = optim.control$gridSize),
                    seq(L[2], U[2], length = optim.control$gridSize))
  }



  M <- npred(object)
  v.simple <- if (Quadratic) {
    length(object@control$colx1.index) == 1 &&
    names(object@control$colx1.index) == "(Intercept)" &&
    (if (any(names(constraints(object)) == "(Intercept)"))
    trivial.constraints(constraints(object))["(Intercept)"] == 1 else TRUE)
  } else {
    FALSE  # To simplify things for "rrvgam" objects
  }
  B1bix <- if (v.simple) {
    matrix(Coefobject@B1, nn, M, byrow = TRUE)
  } else {
    Xlm <- predict.vlm(as(object, "vglm"),  # object,
                       newdata = newdata.orig,
                       type = "Xlm")
    Xlm[, names(object@control$colx1.index), drop = FALSE] %*%
    (if (Quadratic) Coefobject@B1 else object@coefficients[1:M])
  }  # !v.simple







  objfun1 <- function(lv1val,
                      x = NULL, y, w = 1, extraargs) {
      ans <- sum(c(w) * extraargs$Obfunction(
              bnu = c(lv1val),
              y0 = y,
              extra = extraargs$object.extra,
              objfun = extraargs$obfunct,
              Object = extraargs$Object,  # Needed for "rrvgam" objects
              Coefs = extraargs$Coefobject,
              B1bix = extraargs$B1bix,
              misc.list = extraargs$object.misc,
              everything = FALSE,
              mu.function = extraargs$mu.function))
    ans
  }

  objfun2 <- function(lv1val, lv2val,
                      x = NULL, y, w = 1, extraargs) {
      ans <- sum(c(w) * extraargs$Obfunction(
              bnu = c(lv1val, lv2val),
              y0 = y,
              extra = extraargs$object.extra,
              objfun = extraargs$obfunct,
              Object = extraargs$Object,  # Needed for "rrvgam" objects
              Coefs = extraargs$Coefobject,
              B1bix = extraargs$B1bix,
              misc.list = extraargs$object.misc,
              everything = FALSE,
              mu.function = extraargs$mu.function))
    ans
  }



  mu.function <- slot(object@family, "linkinv")
  wvec <- 1  # zz; Assumed here
  mylist <-
    list(object.extra = object@extra,
         Obfunction   = if (Quadratic) .my.calib.objfunction.qrrvglm else
                                       .my.calib.objfunction.rrvgam,
         Coefobject   = Coefobject,
         B1bix        = NA,  # Will be replaced below
         obfunct      = obfunct,
         object.misc  = object@misc,
         Object       = if (Quadratic) 666 else object,
         mu.function  = mu.function)



  init.vals <- matrix(NA_real_, nn, Rank)
  for (i1 in 1:nn) {
    if (optim.control$trace)
      cat("Grid searching initial values for observation",
          i1, "-----------------\n")


    y0 <- newdata[i1, , drop = FALSE]  # drop added 20150624
    mylist$B1bix <- B1bix[i1, ]
    try.this <- if (Rank == 1)
      grid.search(initial.vals[, 1],
                  objfun = objfun1, y = y0 , w = wvec,
                  ret.objfun = TRUE,
                  extraargs = mylist) else
      grid.search2(initial.vals[, 1], initial.vals[, 2],
                   objfun = objfun2, y = y0, w = wvec,
                   ret.objfun = TRUE,
                   extraargs = mylist)
    lv1.init <- try.this[if (Rank == 1) "Value" else "Value1"]
    lv2.init <- if (Rank >= 2) try.this["Value2"] else NULL

    init.vals[i1, ] <- c(lv1.init, lv2.init)
  }  # for i1





  BestOFpar <- matrix(NA_real_, nn, Rank)
  BestOFvalues <- rep(NA_real_, nn)  # Best OF objective function values

  for (i1 in 1:nn) {
    if (optim.control$trace) {
      cat("\nOptimizing for observation", i1, "-----------------\n")
      flush.console()
    }

    ans <-
      optim(par = init.vals[i1, ],
            fn = if (Quadratic) .my.calib.objfunction.qrrvglm else
                                .my.calib.objfunction.rrvgam,
            method = optim.control$Method.optim,  # "BFGS" or "CG" or...
            control = c(fnscale = ifelse(minimize.obfunct, 1, -1),
                        use.optim.control),  # as.vector() needed
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              extra = object@extra,
              objfun = obfunct,
              Object = if (Quadratic) 666 else object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              everything = FALSE,  # Differs from below
              mu.function = mu.function)


        if (optim.control$trace) {
          if (ans$convergence == 0)
            cat("Successful convergence\n") else
            cat("Unsuccessful convergence\n")
          flush.console()
        }
        if (ans$convergence == 0) {
          BestOFpar[i1, ] <- ans$par
          BestOFvalues[i1] <- ans$value
        }  # else do nothing since NA_real_ signifies convergence failure
  }  # for i1





  pretty <- function(BestOFpar, newdata, Rank) {
    if (Rank == 1) {
      if (!is.null(dimnames(newdata)[[1]])) {
        BestOFpar <- c(BestOFpar)
        names(BestOFpar) <- dimnames(newdata)[[1]]
      }
    } else {
      dimnames(BestOFpar) <-
        list(dimnames(newdata)[[1]],
             if (Rank == 1) "latvar" else
                            paste("latvar", 1:Rank, sep = ""))
    }
    BestOFpar
  }  # pretty

  BestOFpar <- pretty(BestOFpar, newdata, Rank)
  attr(BestOFpar,"objectiveFunction") <-
    pretty(BestOFvalues, newdata, Rank = 1)
  if (type == "latvar")
    return(BestOFpar)






    choose.fun <- if (Quadratic) .my.calib.objfunction.qrrvglm else
                                 .my.calib.objfunction.rrvgam



  etaValues <- matrix(NA_real_, nn, M)
  muValues <- matrix(NA_real_, nn, ncol(fitted(object)))
  vcValues <- if (Quadratic) array(0, c(Rank, Rank, nn)) else NULL



  if (optim.control$trace)
    cat("\n")

  for (i1 in 1:nn) {
    if (optim.control$trace) {
      cat("Evaluating quantities for observation", i1,
          "-----------------\n")
      flush.console()
    }
    ans <- choose.fun(
              bnu = if (Rank == 1) BestOFpar[i1] else BestOFpar[i1, ],
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              extra = object@extra,
              objfun = obfunct,
              Object = if (Quadratic) 666 else object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              everything = TRUE,
              mu.function = mu.function)

    muValues[i1, ] <- ans$mu
    etaValues[i1, ] <- ans$eta
    if (Quadratic)
      vcValues[, , i1] <- ans$vcmat  # Might be NULL, e.g., "rvgam"
  }  # for i1



  dimnames(muValues) <- dimnames(newdata)
  dimnames(etaValues) <- list(dimnames(newdata)[[1]],
                              dimnames(object@predictors)[[2]])
  if (Quadratic)
    dimnames(vcValues) <- list(as.character(1:Rank),
                               as.character(1:Rank),
                               dimnames(newdata)[[1]])

  switch(type,
         response   = muValues,
         predictors = etaValues,
         vcov       = vcValues,
         all3or4 = list(latvar     = BestOFpar,
                        predictors = etaValues,
                        response   = muValues,
                        vcov       = if (Quadratic) vcValues else NULL))
}  # calibrate.qrrvglm







.my.calib.objfunction.qrrvglm <-
  function(bnu, y0, extra = NULL,
           objfun, Coefs, Object,
           B1bix,
           misc.list,
           everything = TRUE,
           mu.function) {

  bnumat <- cbind(bnu)
  Rank <- length(bnu)
  eta <- cbind(c(B1bix)) + Coefs@A %*% bnumat
  M <- misc.list$M
  check.eta <- matrix(0, M, 1)
  for (ss in 1:M) {
    temp <- Coefs@D[, , ss, drop = FALSE]
    dim(temp) <- dim(temp)[1:2]  # c(Rank, Rank)
    eta[ss, 1] <- eta[ss, 1] + t(bnumat) %*% temp %*% bnumat

    if (FALSE) {
    warning("this line is wrong:")
    alf <- loge(Coefs@Maximum[ss])  # zz get the link function
    tolmat <- Coefs@Tolerance[, , ss, drop = FALSE]
    check.eta[ss, 1] <- alf - 0.5 * t(bnumat) %*%
                        solve(tolmat) %*% bnumat  
    }  # FALSE
  }  # for ss
  eta <- matrix(eta, 1, M, byrow = TRUE)
  mu <- rbind(mu.function(eta, extra))  # Make sure it has one row
  value <- objfun(mu = mu, y = y0,
                  w = 1,  # ignore prior.weights on the object
                  residuals = FALSE, eta = eta, extra = extra)
  if (everything) {
    vcmat <- matrix(0, Rank, Rank)
    for (ss in 1:M) {
      vec1 <- cbind(Coefs@A[ss, ]) +
                    2 * matrix(Coefs@D[, , ss], Rank, Rank) %*% bnumat
      vcmat <- vcmat + mu[1, ss] * vec1 %*% t(vec1)
    }
    vcmat <- solve(vcmat)
  } else {
    vcmat <- NULL
  }
  if (everything)
    list(eta = eta,
         mu = mu,
         value = value,
         vcmat = vcmat) else
    value
}  # .my.calib.objfunction.qrrvglm






.my.calib.objfunction.rrvgam <-
  function(bnu, y0, extra = NULL,
           objfun,
           Object,  # Needed for "rrvgam" objects
           Coefs,
           B1bix,  # Actually not needed here
           misc.list,
           everything = TRUE,
           mu.function) {
    Rank <- length(bnu)
    NOS <- Coefs@NOS
    eta <- matrix(NA_real_, 1, NOS)
    for (jlocal in 1:NOS) {
      eta[1, jlocal] <- predictrrvgam(Object, grid = bnu, sppno = jlocal,
                                      Rank = Rank, deriv = 0)$yvals
    }
    mu <- rbind(mu.function(eta, extra))  # Make sure it has one row
    value <- objfun(mu = mu, y = y0,
                   w = 1,  # ignore prior.weights on the object
                   residuals = FALSE, eta = eta, extra = extra)
    vcmat <- NULL  # No theory as of yet to compute the vcmat
  if (everything)
    list(eta = eta,
         mu = mu,
         value = value,
         vcmat = vcmat) else
    value
}  # .my.calib.objfunction.rrvgam

















calibrate.rrvglm <-
  function(object,
           newdata = NULL,
           type = c("latvar", "predictors", "response", "vcov", "all3or4"),
           initial.vals = NULL,  # For one observation only
           ...) {


  Quadratic <- FALSE  # Because this function was adapted from CQO code.
  newdata.orig <- newdata
  if (!length(newdata)) {
    newdata <- data.frame(depvar(object))
  }

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("latvar", "predictors",
                            "response", "vcov", "all3or4"))[1]

  if (!all(weights(object, type = "prior") == 1))
    warning("not all the prior weights of 'object' are 1; assuming ",
            "they all are here")


  if (!is.data.frame(newdata))
    newdata <- data.frame(newdata)
  if (!all(object@misc$ynames %in% colnames(newdata)))
    stop("need the following colnames in 'newdata': ",
         paste(object@misc$ynames, collapse = ", "))
  newdata <- newdata[, object@misc$ynames, drop = FALSE]


  if (!is.matrix(newdata))
    newdata <- as.matrix(newdata)


  nn <- nrow(newdata)


  obfunct <- slot(object@family, object@misc$criterion)
  minimize.obfunct <- object@control$min.criterion  # zz
  if (!is.logical(minimize.obfunct))
    stop("object@control$min.criterion is not a logical")
  minimize.obfunct <- as.vector(minimize.obfunct)
  optim.control <- calibrate.rrvglm.control(object = object, ...)

  use.optim.control <- optim.control
  use.optim.control$Method.optim <-
  use.optim.control$gridSize <-
  use.optim.control$varI.latvar <- NULL


  if ((Rank <- object@control$Rank) > 3)
    stop("currently can only handle Rank = 1, 2 and 3")
  Coefobject <- if (Quadratic) {
    Coef(object, varI.latvar = optim.control$varI.latvar)
  } else {
    Coef(object)
  }



 
  if (!length(initial.vals)) {
    L <- apply(latvar(object), 2, min)
    U <- apply(latvar(object), 2, max)
    initial.vals <- if (Rank == 1)
      cbind(seq(L, U, length = optim.control$gridSize)) else
      if (Rank == 2)
      expand.grid(seq(L[1], U[1], length = optim.control$gridSize),
                  seq(L[2], U[2], length = optim.control$gridSize)) else
      expand.grid(seq(L[1], U[1], length = optim.control$gridSize),
                  seq(L[2], U[2], length = optim.control$gridSize),
                  seq(L[3], U[3], length = optim.control$gridSize))
  }  # !length(initial.vals)


  v.simple <- length(object@control$colx1.index) == 1 &&
              names(object@control$colx1.index) == "(Intercept)" &&
              (if (any(names(constraints(object)) == "(Intercept)"))
    trivial.constraints(constraints(object))["(Intercept)"] == 1 else TRUE)
  B1bix <- if (v.simple) {
    matrix(Coefobject@B1, nn, M, byrow = TRUE)
  } else {
    Xlm <- predict.vlm(as(object, "vglm"),  # object,
                       newdata = newdata.orig,
                       type = "Xlm")
    Xlm[, names(object@control$colx1.index), drop = FALSE] %*%
    Coefobject@B1
  }  # !v.simple




  objfun1 <- function(lv1val,
                      x = NULL, y, w = 1, extraargs) {
      ans <- sum(c(w) * extraargs$Obfunction(
              bnu = c(lv1val),
              y0 = y,
              extra = extraargs$object.extra,
              objfun = extraargs$obfunct,
              Object = extraargs$Object,
              Coefs = extraargs$Coefobject,
              B1bix = extraargs$B1bix,
              misc.list = extraargs$object.misc,
              everything = FALSE,
              mu.function = extraargs$mu.function))
    ans
  }

  objfun2 <- function(lv1val, lv2val,
                      x = NULL, y, w = 1, extraargs) {
      ans <- sum(c(w) * extraargs$Obfunction(
              bnu = c(lv1val, lv2val),
              y0 = y,
              extra = extraargs$object.extra,
              objfun = extraargs$obfunct,
              Object = extraargs$Object,
              Coefs = extraargs$Coefobject,
              B1bix = extraargs$B1bix,
              misc.list = extraargs$object.misc,
              everything = FALSE,
              mu.function = extraargs$mu.function))
    ans
  }

  objfun3 <- function(lv1val, lv2val, lv3val,
                      x = NULL, y, w = 1, extraargs) {
      ans <- sum(c(w) * extraargs$Obfunction(
              bnu = c(lv1val, lv2val, lv3val),
              y0 = y,
              extra = extraargs$object.extra,
              objfun = extraargs$obfunct,
              Object = extraargs$Object,
              Coefs = extraargs$Coefobject,
              B1bix = extraargs$B1bix,
              misc.list = extraargs$object.misc,
              everything = FALSE,
              mu.function = extraargs$mu.function))
    ans
  }


  mu.function <- slot(object@family, "linkinv")
  wvec <- 1  # zz; Assumed here
  mylist <-
    list(object.extra = object@extra,
         Obfunction  = .my.calib.objfunction.rrvglm,
         Coefobject   = Coefobject,
         B1bix        = NA,  # Will be replaced below
         obfunct      = obfunct,
         object.misc  = object@misc,
         Object       = 666,  # object,
         mu.function  = mu.function)


  M <- npred(object)
  init.vals <- matrix(NA_real_, nn, Rank)
  for (i1 in 1:nn) {
    if (optim.control$trace)
      cat("Grid searching initial values for observation",
          i1, "-----------------\n")


    y0 <- newdata[i1, , drop = FALSE]  # drop added 20150624
    mylist$B1bix <- B1bix[i1, ]
    try.this <- if (Rank == 1)
      grid.search(initial.vals[, 1],
                  objfun = objfun1, y = y0 , w = wvec,
                  ret.objfun = TRUE,
                  extraargs = mylist) else if (Rank == 2)
      grid.search2(initial.vals[, 1], initial.vals[, 2],
                   objfun = objfun2, y = y0, w = wvec,
                   ret.objfun = TRUE,
                   extraargs = mylist) else
      grid.search3(initial.vals[, 1], initial.vals[, 2], initial.vals[, 3], 
                   objfun = objfun3, y = y0, w = wvec,
                   ret.objfun = TRUE,
                   extraargs = mylist)
    lv1.init <- try.this[if (Rank == 1) "Value" else "Value1"]
    lv2.init <- if (Rank >= 2) try.this["Value2"] else NULL
    lv3.init <- if (Rank >= 3) try.this["Value3"] else NULL

    init.vals[i1, ] <- c(lv1.init, lv2.init, lv3.init)
  }  # for i1





  BestOFpar <- matrix(NA_real_, nn, Rank)
  BestOFvalues <- rep(NA_real_, nn)  # Best OF objective function values

  for (i1 in 1:nn) {
    if (optim.control$trace) {
      cat("\nOptimizing for observation", i1, "-----------------\n")
      flush.console()
    }
    ans <-
      optim(par = init.vals[i1, ],
            fn = .my.calib.objfunction.rrvglm,
            method = optim.control$Method.optim,  # "BFGS" or "CG" or...
            control = c(fnscale = ifelse(minimize.obfunct, 1, -1),
                        use.optim.control),  # as.vector() needed
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              extra = object@extra,
              objfun = obfunct,
              Object = 666,  # object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              everything = FALSE,  # Differs from below
              mu.function = mu.function)

        if (optim.control$trace) {
          if (ans$convergence == 0)
            cat("Successful convergence\n") else
            cat("Unsuccessful convergence\n")
          flush.console()
        }
        if (ans$convergence == 0) {
          BestOFpar[i1, ] <- ans$par
          BestOFvalues[i1] <- ans$value
        }  # else do nothing since NA_real_ signifies convergence failure
  }  # for i1





  pretty <- function(BestOFpar, newdata, Rank) {
    if (Rank == 1) {
      BestOFpar <- c(BestOFpar)
      names(BestOFpar) <- dimnames(newdata)[[1]]
    } else
      dimnames(BestOFpar) <-
        list(dimnames(newdata)[[1]],
             if (Rank == 1) "latvar" else
                            paste("latvar", 1:Rank, sep = ""))
    BestOFpar
  }  # pretty


  BestOFpar <- pretty(BestOFpar, newdata, Rank)
  attr(BestOFpar,"objectiveFunction") <-
    pretty(BestOFvalues, newdata, Rank = 1)
  if (type == "latvar")
    return(BestOFpar)





  etaValues <- matrix(NA_real_, nn, M)
  muValues <- matrix(NA_real_, nn, ncol(fitted(object)))
  vcValues <- if (Quadratic) array(0, c(Rank, Rank, nn)) else NULL

  if (optim.control$trace)
    cat("\n")

  for (i1 in 1:nn) {
    if (optim.control$trace) {
      cat("Evaluating quantities for observation", i1,
          "-----------------\n")
      flush.console()
    }
    ans <- .my.calib.objfunction.rrvglm(
              bnu = if (Rank == 1) BestOFpar[i1] else BestOFpar[i1, ],
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              extra = object@extra,
              objfun = obfunct,
              Object = 666,  # object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              everything = TRUE,
              mu.function = mu.function)

    muValues[i1, ] <- ans$mu
    etaValues[i1, ] <- ans$eta
    if (Quadratic)
      vcValues[, , i1] <- ans$vcmat  # Might be NULL, e.g., "rvgam"
  }  # for i1


  dimnames(muValues) <- dimnames(newdata)
  dimnames(etaValues) <- list(dimnames(newdata)[[1]],
                              dimnames(object@predictors)[[2]])
  if (Quadratic)
    dimnames(vcValues) <- list(as.character(1:Rank),
                               as.character(1:Rank),
                               dimnames(newdata)[[1]])

  switch(type,
         response   = muValues,
         predictors = etaValues,
         vcov       = vcValues,
         all3or4 = list(latvar     = BestOFpar,
                        predictors = etaValues,
                        response   = muValues,
                        vcov       = if (Quadratic) vcValues else NULL))
}  # calibrate.rrvglm





.my.calib.objfunction.rrvglm <-
  function(bnu,
           y0, extra = NULL,
           objfun, Coefs,
           Object,  # Not needed
           B1bix,
           misc.list,
           everything = TRUE,
           mu.function) {



  bnumat <- cbind(bnu)
  Rank <- length(bnu)
  eta <- cbind(c(B1bix)) + Coefs@A %*% bnumat
  M <- misc.list$M
  eta <- matrix(eta, 1, M, byrow = TRUE)
  mu <- rbind(mu.function(eta, extra = extra))  # Make sure it has 1 row
  value <- objfun(mu = mu, y = y0,
                  w = 1,  # ignore prior.weights on the object zz
                  residuals = FALSE, eta = eta, extra = extra)
  if (everything) {
    vcmat <- matrix(0, Rank, Rank)
    for (ss in 1:M) {
      vec1 <- cbind(Coefs@A[ss, ])
      vcmat <- vcmat + mu[1, ss] * vec1 %*% t(vec1)
    }
    vcmat <- solve(vcmat)
  } else {
    vcmat <- NULL
  }
  if (everything)
    list(eta   = eta,
         mu    = mu,
         value = value,
         vcmat = vcmat) else
    value
}  # .my.calib.objfunction.rrvglm




calibrate.rrvglm.control <-
  function(object,
           trace = FALSE,  # passed into optim()
           Method.optim = "BFGS",  # passed into optim(method = Method)
           gridSize = if (Rank == 1) 7 else 5,
           ...) {


  Rank <- object@control$Rank
  if (!is.Numeric(gridSize, positive = TRUE,
                  integer.valued = TRUE, length.arg = 1))
    stop("bad input for argument 'gridSize'")
  if (gridSize < 2)
    stop("argument 'gridSize' must be >= 2")

  list(
       trace = as.numeric(trace)[1],
       Method.optim = Method.optim,
       gridSize = gridSize
      )
}  # calibrate.rrvglm.control








# These functions are
# Copyright (C) 1998-2018 T.W. Yee, University of Auckland.
# All rights reserved.


















fitCQOagain <-
  function(Object,
           which.spp = 1,
           error.check = TRUE) {

  lv1 <- latvar(Object)

  X.lm1 <- model.matrix(Object, type = "vlm")  # Has latvar vars too.
  etol <- Object@control$eq.tolerances
  Itol <- Object@control$I.tolerances
  colx1.index <- Object@control$colx1.index
  colx2.index <- Object@control$colx2.index
  nice31 <- !etol || Itol
  NOS <- npred(Object) / npred(Object, type = "one.response")


  if (length(colx1.index) != 1 ||
      names(colx1.index) != "(Intercept)")
    stop("calibrate() can only operate when B1 is intercept-only")

  if (!etol) {


 print("could fit all species separately ,,,,,,,,,,,,,,,,,,,,,")




  if (FALSE) {
    X.lm <-  # Same order as in the C code:
    cbind(lv1,
          lv1^2,  # Always here even when I.tol is TRUE
          X.lm1[, colx1.index])

    dev.total <- 0
    index.spp <- c(3, 1, 2)  # Order becomes: Intercept, \nu, and \nu^2.
    Y.all <- depvar(Object)
    loop.spp. <- if (error.check) 1:NOS else which.spp
    for (spp. in loop.spp.) {
      y.use <- Y.all[, spp.]  # One species response
      fit1.tmp <- vglm(y.use ~ X.lm - 1,
                       family = Object@family, trace = FALSE)
      dev.total <- dev.total + deviance(fit1.tmp)
      if (spp. == which.spp) {
        vcov.spp <- vcov(fit1.tmp)[index.spp, index.spp]
        cofs.spp <- coef(fit1.tmp)[index.spp]
      }
    }
  }




 print("attempt2")
    index.spp <- which.spp + (0:2) * NOS
    lv12 <- lv1^2
    clist12 <- list("(Intercept)" = diag(NOS),
                    "lv1"         = diag(NOS),
                    "lv12"        = diag(NOS))
    fit1 <-
      vglm(depvar(Object) ~ lv1 + lv12,
           family = Object@family,  # Handles link = "cloglog", etc
           etastart = predict(Object),  # Gives immediate convergence
           constraints = clist12,
           trace = TRUE)

    dev.total <- deviance(fit1)
    vcov.spp <- vcov(fit1)[index.spp, index.spp]
    cofs.spp <- coef(fit1)[index.spp]





  } else {
 print("fitting all species simultaneously ,,,,,,,,,,,,,,,,,,,,,")


 print("attempt2")
    index.spp <- c(which.spp, NOS + which.spp, 2 * NOS + 1)
    lv12 <- lv1^2
    clist12 <- list("(Intercept)" = diag(NOS),
                    "lv1"         = diag(NOS),
                    "lv12"        = matrix(1, NOS, 1))
    fit1 <-
      vglm(depvar(Object) ~ lv1 + lv12,
           family = Object@family,  # Handles link = "cloglog", etc
           etastart = predict(Object),  # Gives immediate convergence
           constraints = clist12,
           trace = TRUE)
 print("coef(fit1)")
 print( coef(fit1) )

    dev.total <- deviance(fit1)
    vcov.spp <- vcov(fit1)[index.spp, index.spp]
    cofs.spp <- coef(fit1)[index.spp]
  }


  if (error.check) {
    check.diff <- deviance(Object) - dev.total  # Should be 0
    if (abs(check.diff) / (1 + deviance(Object)) > 0.01) {
      warning("the two models do not appear to be the same")
    } else {
      cat("Both models appear to be the same (a good thing).\n")
    }
  }


  if (cofs.spp[3] >= 0)
    warning("the 3rd coefficient should be negative!!")


  list(Coefficients = cofs.spp,
       VarCov       = vcov.spp)
}  # fitCQOagain







fitCLOagain <-
  function(Object,
           which.spp = 1,
           error.check = TRUE) {

 print("in fitCLOagain()")
  lvmat <- latvar(Object)

  X.lm1 <- model.matrix(Object, type = "lm")  # Has latvar vars too.
  colx1.index <- Object@control$colx1.index
  colx2.index <- Object@control$colx2.index
  NOS <- npred(Object) / npred(Object, type = "one.response")


  if (length(colx1.index) != 1 ||
      names(colx1.index) != "(Intercept)")
    warning("calibrate() can only operate when B1 is intercept-only")

  X.lm.use <- X.lm1
 print("head(X.lm.use)")
 print( head(X.lm.use) )

 print("attempt2")
  index.spp <- TRUE  # which.spp + (0:2) * NOS
  clist12 <- constraints(Object, type = "term")


  data.frame.use <- if(exists(Object@misc$dataname))
      get(Object@misc$dataname) else list()

  fit1 <-
    vglm(
         formula = formula(Object@call),
         data = data.frame.use,
         family = Object@family,  # Handles link = "cloglog", etc
         etastart = predict(Object),  # Gives immediate convergence
         constraints = clist12,
         trace = TRUE)

  dev.total <- deviance(fit1)
  vcov.spp <- vcov(fit1)[index.spp, index.spp]
  cofs.spp <- coef(fit1)[index.spp]



  if (error.check) {
    check.diff <- deviance(Object) - dev.total  # Should be 0
    if (abs(check.diff) / (1 + deviance(Object)) > 0.01) {
      warning("the two models do not appear to be the same")
    } else {
      cat("Both models appear to be the same (a good thing).\n")
    }
  }



  list(Coefficients = cofs.spp,
       VarCov       = vcov.spp)
}  # fitCLOagain





getMaxOptTol <-
  function(coeff3, vc, varcov.arg = FALSE) {

  if (length(coeff3) != 3 || coeff3[3] > 0)
    stop("bad input for argument 'coeff3'")
  beta0 <- coeff3[1]
  beta1 <- coeff3[2]
  beta2 <- coeff3[3]
  uj <- -0.5 * beta1 / beta2
  tolj <- 1 / sqrt(-2 * beta2)  # Like a sdev, it's sqrt(Tol) really.
  alphaj <- beta0 - 0.25 * (beta1^2) / beta2
  dimn <- c("Maximum", "Optimum", "Tolerance")
  ans <- c(alphaj, uj, tolj)
  names(ans) <- dimn
  if (varcov.arg) {
    dmax.dbeta <- c(1, uj, uj^2)
    dopt.dbeta <- c(0, tolj^2, 2 * uj * tolj)
    dtol.dbeta <- c(0, 0, tolj^3)
    dTheta.dbeta <- cbind(dmax.dbeta, dopt.dbeta, dtol.dbeta)
    varcov <- t(dTheta.dbeta) %*% vc %*% dTheta.dbeta
    dimnames(varcov) <- list(dimn, dimn)
    list(coeffs = ans, varcov = varcov)
  } else {   
    ans
  }
}  # getMaxOptTol











getG.calib <-
  function(nu0, coeff3 = NULL, 
           uj = NULL, tolj = NULL, alphaj = NULL,
           numerator = c("xi", "eta"), ...) {

  if (mode(numerator) != "character" && mode(numerator) != "name")
    numerator <- as.character(substitute(numerator))
  numerator <- match.arg(numerator, c("xi", "eta"))[1]


  if (!is.null(coeff3)) {
    if (length(coeff3) != 3 || coeff3[3] > 0)
      stop("bad input for argument 'coeff3'")
    beta1 <- coeff3[2]  # beta0 <- coeff3[1] is unneeded
    beta2 <- coeff3[3]
  }
  if (is.null(uj))     uj <- -0.5 * beta1 / beta2
  if (is.null(tolj))   tolj <- 1 / sqrt(-2 * beta2)
  switch(numerator,
         "xi"  = c(0, 1, 2 * nu0),
         "eta" = c(1, nu0,
                   nu0^2 - 2 * uj * (nu0 - uj) * (1 - 1 / tolj)))
}  # getG.calib







dzwald.qrrvglm <-
  function(bnu0, y0, extra = NULL,
           objfun, Coefs, Object,
           B1bix,
           misc.list,
           everything = TRUE,
           mu.function,
           link.function, earg = list(),  # These go together
           which.spp = 1) {
 print("y0")
 print( y0)

  
  M <- misc.list$M
  if ((Rank <- length(bnu0)) != 1)
    stop("can only handle rank-1 objects")

  linkfun <- Object@family@linkfun
  if (!is.function(linkfun))
    stop("could not obtain linkfun")



  

  aaa <- fitCQOagain(Object, which.spp = which.spp)
  cofs.spp <- aaa$Coefficients  # 3-vector

  bbb <- getMaxOptTol(cofs.spp,
                      vc = aaa$VarCov, varcov.arg = TRUE)
  Omegamat <- bbb$varcov





  NOS <- M  # TRUE for binomial and Poisson
  Dmat <- Gmat <- 0
  for (spp. in 1:NOS) {
    uj.jay <- Coefs@Optimum[, spp.]
    tolj.jay <- sqrt(Coefs@Tolerance[Rank, Rank, spp.])  # care!!
    alphaj.jay <- linkfun(Coefs@Maximum,
                          extra = Object@extra)[spp.]

    eta0.jay <- alphaj.jay - 0.5 * ((bnu0 - uj.jay) / tolj.jay)^2
    eta0.jay <- rbind(eta0.jay)



    fv0.jay <- mu.function(eta0.jay, extra = Object@extra)
    fv0.jay <- c(fv0.jay)  # Remove array attributes
    dTheta.deta0j <- dtheta.deta(fv0.jay, link.function, earg = earg)
    dTheta.deta0j <- c(dTheta.deta0j)  # Remove array attributes


    xi0.jay <- -(bnu0 - uj.jay)  / tolj.jay^2
    Dmat <- Dmat + dTheta.deta0j * xi0.jay^2  # More general
    dxi0.dtheta <-
      getG.calib(nu0 = bnu0, coeff3 = cofs.spp,  # NULL, 
                 numerator = "xi",
                 which.spp = spp.)
    deta0.dtheta <-
      getG.calib(nu0 = bnu0, coeff3 = cofs.spp,  # NULL, 
                 numerator = "eta",
                 which.spp = spp.)
    Gmat <- Gmat + (y0[1, spp.] - fv0.jay) * dxi0.dtheta -
            xi0.jay * dTheta.deta0j * deta0.dtheta
  }

  1 / Dmat + (rbind(Gmat) %*% Omegamat %*% cbind(Gmat)) / Dmat^2
}  # dzwald.qrrvglm












calibrate.qrrvglm.control <-
  function(object,
           trace = FALSE,  # passed into optim()
           Method.optim = "BFGS",  # passed into optim(method = Method)
           gridSize = ifelse(Rank == 1, 21, 9),
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
           type = c("latvar", "predictors", "response",
                    "vcov", "all3or4"),
           se.type = c("asbefore", "wald"),
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

  get.SEs <- Quadratic && type %in% c("vcov", "all3or4")
  if (mode(se.type) != "character" && mode(se.type) != "name")
    se.type <- as.character(substitute(se.type))
  se.type <- match.arg(se.type, c("asbefore", "wald"))[1]

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


  nn <- nrow(newdata)  # Number of sites to calibrate




  obfunct <- slot(object@family, object@misc$criterion)  # deviance
 print("object@misc$criterion")
 print( object@misc$criterion )
 print("obfunct7")
 print( obfunct )
  minimize.obfunct <-
    if (Quadratic) object@control$min.criterion else
    TRUE  # Logical; TRUE for CAO objects because deviance is minimized
  if (!is.logical(minimize.obfunct))
    stop("object@control$min.criterion is not a logical")
 print("minimize.obfunct")
 print( minimize.obfunct )
  optim.control <- calibrate.qrrvglm.control(object = object, ...)

  use.optim.control <- optim.control
  use.optim.control$Method.optim <-
  use.optim.control$gridSize     <-
  use.optim.control$varI.latvar  <- NULL


  if ((Rank <- object@control$Rank) > 2)
    stop("currently can only handle Rank = 1 and 2")
  Coefobject <- if (Quadratic) {
    Coef(object, varI.latvar = optim.control$varI.latvar)
  } else {
    Coef(object)
  }



  if (!length(initial.vals)) {
    Lvec <- apply(latvar(object), 2, min)
    Uvec <- apply(latvar(object), 2, max)
    initial.vals <- if (Rank == 1)
      cbind(seq(Lvec, Uvec, length = optim.control$gridSize)) else
      expand.grid(seq(Lvec[1], Uvec[1], length = optim.control$gridSize),
                  seq(Lvec[2], Uvec[2], length = optim.control$gridSize))
  }

 print("head(initial.vals)")
 print( head(initial.vals) )
 print("tail(initial.vals)")
 print( tail(initial.vals) )





  M <- npred(object)
  v.simple <- if (Quadratic) {
    length(object@control$colx1.index) == 1 &&
    names(object@control$colx1.index) == "(Intercept)" &&
    (if (any(names(constraints(object)) == "(Intercept)"))
    trivial.constraints(constraints(object))["(Intercept)"] == 1 else
    TRUE)
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
 print("head(B1bix)")
 print( head(B1bix) )
 print("dim(B1bix)")
 print( dim(B1bix) )







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


  choose.fun <- if (Quadratic) .my.calib.objfunction.qrrvglm else
                               .my.calib.objfunction.rrvgam
  mu.function <- slot(object@family, "linkinv")
  wvec <- 1  # zz; Assumed here
  mylist <-
    list(object.extra = object@extra,
         Obfunction   = choose.fun,  # e.g. .my.calib.objfunction.qrrvglm
         Coefobject   = Coefobject,
         B1bix        = NA,  # Will be replaced below
         obfunct      = obfunct,  # deviance
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
                  maximize = !minimize.obfunct,  # Most general.
                  extraargs = mylist) else
      grid.search2(initial.vals[, 1], initial.vals[, 2],
                   objfun = objfun2, y = y0, w = wvec,
                   ret.objfun = TRUE,
                   maximize = !minimize.obfunct,  # Most general.
                   extraargs = mylist)
    lv1.init <- try.this[if (Rank == 1) "Value" else "Value1"]
    lv2.init <- if (Rank >= 2) try.this["Value2"] else NULL

    init.vals[i1, ] <- c(lv1.init, lv2.init)
  }  # for i1

 print("hi34")
 print("head(init.vals)")
 print( head(init.vals) )




  BestOFpar <- matrix(NA_real_, nn, Rank)
  BestOFvalues <- rep(NA_real_, nn)  # Best OF objective function values

  for (i1 in 1:nn) {
    if (optim.control$trace) {
      cat("\nOptimizing for observation", i1, "-----------------\n")
      flush.console()
    }

    ans <-
      optim(par = init.vals[i1, ],
            fn  = choose.fun,
            method = optim.control$Method.optim,  # "BFGS" or "CG" or...
            control = c(fnscale = ifelse(minimize.obfunct, 1, -1),
                        use.optim.control),  # as.vector() needed
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              extra = object@extra,
              objfun = obfunct,  # deviance
              Object = if (Quadratic) 666 else object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              everything = FALSE,  # Differs from Step 4. below
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

 print("hi43")
 print("BestOFvalues")
 print( BestOFvalues )
 print("BestOFpar")
 print( BestOFpar )
 print("dimnames(newdata)")
 print( dimnames(newdata) )




  pretty <- function(BestOFpar, newdata, Rank) {
    if (Rank == 1) {
      if (!is.null(dimnames(newdata)[[1]])) {
        BestOFpar <- c(BestOFpar)
        names(BestOFpar) <- dimnames(newdata)[[1]]
      }
    } else {
      dimnames(BestOFpar) <- list(dimnames(newdata)[[1]],
                                  param.names("latvar", Rank, skip1 = TRUE))
    }
    BestOFpar
  }  # pretty

  BestOFpar <- pretty(BestOFpar, newdata, Rank)
  attr(BestOFpar,"objectiveFunction") <-
    pretty(BestOFvalues, newdata, Rank = 1)
  if (type == "latvar")
    return(BestOFpar)






  etaValues <- matrix(NA_real_, nn, M)
  muValues <- matrix(NA_real_, nn, ncol(fitted(object)))
  if (get.SEs)
    vcValues <- array(0, c(Rank, Rank, nn))


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
              objfun = obfunct,  # deviance
              Object = if (Quadratic) 666 else object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              everything = get.SEs,  # Differs from Step 3.
              mu.function = mu.function)

    muValues[i1, ] <- ans$mu
    etaValues[i1, ] <- ans$eta


    if (get.SEs) {
      vcValues[, , i1] <-
        if (se.type == "wald")
        dzwald.qrrvglm(
              bnu0 = if (Rank == 1) BestOFpar[i1] else BestOFpar[i1, ],
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              extra = object@extra,
              objfun = obfunct,  # deviance
              Object = object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              everything = TRUE,
              mu.function = mu.function,
              link.function = linkfun(object)[i1],
              earg          = object@misc$earg[[i1]], 
              which.spp = i1) else
        ans$vcmat  # Might be NULL, e.g., "rvgam"
    }  # if (get.SEs)
  }  # for i1



  dimnames(muValues) <- dimnames(newdata)
  dimnames(etaValues) <- list(dimnames(newdata)[[1]],
                              dimnames(object@predictors)[[2]])
  if (get.SEs)
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
                        vcov       = if (get.SEs) vcValues else NULL))
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
  mu <- rbind(mu.function(eta, extra))  # Make sure it has 1 row
  obvalue <- objfun(mu = mu, y = y0,
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
         obvalue = obvalue,
         vcmat = vcmat) else
    obvalue
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
      eta[1, jlocal] <- predictrrvgam(Object, grid = bnu,
                          sppno = jlocal, Rank = Rank, deriv = 0)$yvals
    }
    mu <- rbind(mu.function(eta, extra))  # Make sure it has one row
    obvalue <- objfun(mu = mu, y = y0,
                      w = 1,  # ignore prior.weights on the object
                      residuals = FALSE, eta = eta, extra = extra)
    vcmat <- NULL  # No theory as of yet to compute the vcmat
  if (everything)
    list(eta = eta,
         mu = mu,
         obvalue = obvalue,
         vcmat = vcmat) else
    obvalue
}  # .my.calib.objfunction.rrvgam















dzwald.rrvglm <-
  function(bnu0, y0, extra = NULL,
           objfun, Coefs, Object,
           B1bix,
           misc.list,
           everything = TRUE,
           mu.function,
           link.function, earg = list(),  # These go together
           which.spp = 1) {
 print("in dzwald.rrvglm()")
 warning("dzwald.rrvglm() is not working yet")
 print("y0")
 print( y0)

  
  M <- misc.list$M
  if ((Rank <- length(bnu0)) != 1)
    stop("can only handle rank-1 objects")

  linkfun <- Object@family@linkfun
  if (!is.function(linkfun))
    stop("could not obtain linkfun")



  

  aaa <- fitCQOagain(Object, which.spp = which.spp)
  cofs.spp <- aaa$Coefficients  # 3-vector

  bbb <- getMaxOptTol(cofs.spp,
                      vc = aaa$VarCov, varcov.arg = TRUE)
  Omegamat <- bbb$varcov





  NOS <- M  # TRUE for binomial and Poisson
  Dmat <- Gmat <- 0
  for (spp. in 1:NOS) {
    uj.jay <- Coefs@Optimum[, spp.]
    tolj.jay <- sqrt(Coefs@Tolerance[Rank, Rank, spp.])  # care!!
    alphaj.jay <- linkfun(Coefs@Maximum,
                          extra = Object@extra)[spp.]

    eta0.jay <- alphaj.jay - 0.5 * ((bnu0 - uj.jay) / tolj.jay)^2
    eta0.jay <- rbind(eta0.jay)



    fv0.jay <- mu.function(eta0.jay, extra = Object@extra)
    fv0.jay <- c(fv0.jay)  # Remove array attributes
    dTheta.deta0j <- dtheta.deta(fv0.jay, link.function, earg = earg)
    dTheta.deta0j <- c(dTheta.deta0j)  # Remove array attributes


    xi0.jay <- -(bnu0 - uj.jay)  / tolj.jay^2
    Dmat <- Dmat + dTheta.deta0j * xi0.jay^2  # More general
    dxi0.dtheta <-
      getG.calib(nu0 = bnu0, coeff3 = cofs.spp,  # NULL, 
                 numerator = "xi",
                 which.spp = spp.)
    deta0.dtheta <-
      getG.calib(nu0 = bnu0, coeff3 = cofs.spp,  # NULL, 
                 numerator = "eta",
                 which.spp = spp.)
    Gmat <- Gmat + (y0[1, spp.] - fv0.jay) * dxi0.dtheta -
            xi0.jay * dTheta.deta0j * deta0.dtheta
  }

  1 / Dmat + (rbind(Gmat) %*% Omegamat %*% cbind(Gmat)) / Dmat^2
}  # dzwald.rrvglm






calibrate.rrvglm <-
  function(object,
           newdata = NULL,
           type = c("latvar", "predictors", "response", "vcov",
                    "all3or4"),
           se.type = c("asbefore", "wald"),
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
  get.SEs <- type %in% c("vcov", "all3or4")
  if (mode(se.type) != "character" && mode(se.type) != "name")
    se.type <- as.character(substitute(se.type))
  se.type <- match.arg(se.type, c("asbefore", "wald"))[1]

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


  nn <- nrow(newdata)  # Number of sites to calibrate


  obfunct <- slot(object@family, object@misc$criterion)
 print("object@misc$criterion")
 print( object@misc$criterion )
 print("obfunct7")
 print( obfunct )
  minimize.obfunct <- object@control$min.criterion  # deviance
  if (!is.logical(minimize.obfunct))
    stop("object@control$min.criterion is not a logical")
  minimize.obfunct <- as.vector(minimize.obfunct)
 print("minimize.obfunct")
 print( minimize.obfunct )
  optim.control <- calibrate.rrvglm.control(object = object, ...)
 print("optim.control")
 print( optim.control )

  use.optim.control <- optim.control
  use.optim.control$Method.optim <-
  use.optim.control$gridSize     <-
  use.optim.control$varI.latvar  <- NULL


  if ((Rank <- object@control$Rank) > 3)
    stop("currently can only handle Rank = 1, 2 and 3")
  Coefobject <- if (Quadratic) {
    Coef(object, varI.latvar = optim.control$varI.latvar)
  } else {
    Coef(object)
  }



 
  if (!length(initial.vals)) {
    Lvec <- apply(latvar(object), 2, min)
    Uvec <- apply(latvar(object), 2, max)
    initial.vals <- if (Rank == 1)
      cbind(seq(Lvec, Uvec, length = optim.control$gridSize)) else
      if (Rank == 2)
  expand.grid(seq(Lvec[1], Uvec[1], length = optim.control$gridSize),
              seq(Lvec[2], Uvec[2], length = optim.control$gridSize)) else
  expand.grid(seq(Lvec[1], Uvec[1], length = optim.control$gridSize),
              seq(Lvec[2], Uvec[2], length = optim.control$gridSize),
              seq(Lvec[3], Uvec[3], length = optim.control$gridSize))
 print("head(initial.vals)")
 print( head(initial.vals) )
 print("tail(initial.vals)")
 print( tail(initial.vals) )
  }  # !length(initial.vals)




  M <- npred(object)
  v.simple <- length(object@control$colx1.index) == 1 &&
              names(object@control$colx1.index) == "(Intercept)" &&
              (if (any(names(constraints(object)) == "(Intercept)"))
      trivial.constraints(constraints(object))["(Intercept)"] == 1 else
      TRUE)
  B1bix <- if (v.simple) {
    matrix(Coefobject@B1, nn, M, byrow = TRUE)
  } else {
 print("not v.simple")
    Xlm <- predict.vlm(as(object, "vglm"),  # object,
                       newdata = newdata.orig,
                       type = "Xlm")
 print("head(Xlm)0")
 print( head(Xlm) )
    Xlm[, names(object@control$colx1.index), drop = FALSE] %*%
    Coefobject@B1
  }  # !v.simple
 print("head(B1bix)")
 print( head(B1bix) )
 print("dim(B1bix)")
 print( dim(B1bix) )




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
                  maximize = !minimize.obfunct,  # Most general.
                  extraargs = mylist) else if (Rank == 2)
      grid.search2(initial.vals[, 1], initial.vals[, 2],
                   objfun = objfun2, y = y0, w = wvec,
                   ret.objfun = TRUE,
                   maximize = !minimize.obfunct,  # Most general.
                   extraargs = mylist) else
      grid.search3(initial.vals[, 1], initial.vals[, 2],
                   initial.vals[, 3], 
                   objfun = objfun3, y = y0, w = wvec,
                   ret.objfun = TRUE,
                   maximize = !minimize.obfunct,  # Most general.
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
            fn  = .my.calib.objfunction.rrvglm,
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


 print("hi43")
 print("BestOFvalues")
 print( BestOFvalues )
 print("BestOFpar")
 print( BestOFpar )
 print("dimnames(newdata)")
 print( dimnames(newdata) )





  pretty <- function(BestOFpar, newdata, Rank) {
    if (Rank == 1) {
      BestOFpar <- c(BestOFpar)
      names(BestOFpar) <- dimnames(newdata)[[1]]
    } else
      dimnames(BestOFpar) <- list(dimnames(newdata)[[1]],
                                  param.names("latvar", Rank, skip1 = TRUE))
    BestOFpar
  }  # pretty


  BestOFpar <- pretty(BestOFpar, newdata, Rank)
  attr(BestOFpar,"objectiveFunction") <-
    pretty(BestOFvalues, newdata, Rank = 1)
  if (type == "latvar")
    return(BestOFpar)





  etaValues <- matrix(NA_real_, nn, M)
  muValues <- matrix(NA_real_, nn, ncol(fitted(object)))
  if (get.SEs)
    vcValues <- array(0, c(Rank, Rank, nn))


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
              objfun = obfunct,  # deviance
              Object = 666,  # object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              everything = get.SEs,  # Differs from Step 3.
              mu.function = mu.function)

    muValues[i1, ] <- ans$mu
    etaValues[i1, ] <- ans$eta


    if (get.SEs)
      vcValues[, , i1] <- if (se.type == "wald")
        dzwald.rrvglm(
              bnu0 = if (Rank == 1) BestOFpar[i1] else BestOFpar[i1, ],
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              extra = object@extra,
              objfun = obfunct,  # deviance
              Object = object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              everything = TRUE,
              mu.function = mu.function,
              link.function = linkfun(object)[i1],
              earg          = object@misc$earg[[i1]], 
              which.spp = i1) else
         ans$vcmat


  }  # for i1


  dimnames(muValues) <- dimnames(newdata)
  dimnames(etaValues) <- list(dimnames(newdata)[[1]],
                              dimnames(object@predictors)[[2]])
  if (get.SEs)
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
                        vcov       = if (get.SEs) vcValues else NULL))
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
  obvalue <- objfun(mu = mu, y = y0,
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
         obvalue = obvalue,
         vcmat = vcmat) else
    obvalue
}  # .my.calib.objfunction.rrvglm




calibrate.rrvglm.control <-
  function(object,
           trace = FALSE,  # passed into optim()
           Method.optim = "BFGS",  # passed into optim(method = Method)
           gridSize = ifelse(Rank == 1, 17, 9),
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








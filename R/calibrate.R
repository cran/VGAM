# These functions are
# Copyright (C) 1998-2023 T.W. Yee, University of Auckland.
# All rights reserved.


























phifun.integrand <- function(charfun, muvec, parlink, earg,
                             cc, xi, ee, ss, sigma2.fuzz = 1e-8) {
  cnew <- cc / ss

  eta.use <- theta2eta(muvec, link = parlink, earg = earg)  # Dangerous
  chfw <- charfun(x = xi * cnew, eta = eta.use, extra = list())
  chfw <- prod(chfw)
  phi <- chfw * exp(-1i * ee * cnew) * exp(-sigma2.fuzz * (cnew^2) / 2)
  phi
}  # phifun.integrand 



cdffun.integrand <-
  function(Object, A.qnorm = 4, wtil, lam, xi, ee, ss,
           sigma2.fuzz = 1e-8,
           nodes = sqrt(2) * ghn100,
           qwts  = sqrt(2) * ghw100 / (2 * pi)) {
  if (!(is(Object, "qrrvglm") || is(Object,  "rrvglm")))
    stop("argument 'Object' must be a 'qrrvglm' or 'rrvglm' Object")
  famfun <- Object@family
  infos <- famfun@infos()
  charfun <- if (is.logical(infos$charfun) && infos$charfun)
    famfun@charfun else stop("the family function has no charfun slot")
  parlink <- linkfun(Object)[1]  # All the same, choose the first
  Earg <- Object@misc$earg[[1]]  # Dangerous

  N <- length(nodes)
  intgrnd <- numeric(N)
  for (k in 1:N) {
    curnode <- nodes[k]
    exptrm <- (exp( 1i * curnode * A.qnorm) -
               exp(-1i * curnode * wtil)) / (1i * curnode) 
    intgrnd[k] <- phifun.integrand(charfun = charfun,
        muvec = lam, parlink = parlink, earg = Earg,  # Dangerous
        cc = curnode, xi = xi,
        ee = ee, ss = ss, sigma2.fuzz = sigma2.fuzz) *
                  exptrm * exp(0.5 * curnode^2)
  }
  sum(qwts * intgrnd)
}  # cdffun.integrand 




charfun.cqo.cdf <-
  function(
                bnu,
                y0,
                extra ,
                objfun,
                Object,
                Coefs,
                B1bix,
                misc.list = list(),
                Everything,
                mu.function,


           A.qnorm = 4,
           sigma2.fuzz = 1e-8, lower.latvar = NULL, upper.latvar = NULL,
           nodes = sqrt(2) * ghn100,
           qwts  = sqrt(2) * ghw100 / (2 * pi)) {
  wt.mean <- function(x, y, alpha = 0.5) (1 - alpha) * x + alpha * y
  site.scores <- latvar(Object)
  range.ss <- range(site.scores)
  if (is.null(lower.latvar))
    lower.latvar <- wt.mean(range.ss[1], range.ss[2],  0.0)
  if (is.null(upper.latvar))
    upper.latvar <- wt.mean(range.ss[1], range.ss[2],  1.0)

  famfun <- Object@family
  infos <- famfun@infos()
  charfun <- if (is.logical(infos$charfun) && infos$charfun)
    famfun@charfun else stop("the family function has no charfun slot")
  Earg <- Object@misc$earg[[1]]  # Dangerous
  
  uopt <- Opt(Object)
  tol <- Tol(Object)[1, 1, ]  # Interpreted as a variance, not SDs
  parlink <- linkfun(Object)[1]  # All the same, choose the first
  alf <- theta2eta(Max(Object), parlink, earg = list())




    nu0 <- bnu
    qtrm <- (nu0 - uopt) / sqrt(tol)  # sqrt(tol), not tol
    eta <- alf - 0.5 * qtrm^2

    muvec <- eta2theta(eta, parlink, earg = Earg)  # Dangerous
    lam <- muvec

    xi <- -qtrm / sqrt(tol)  # sqrt(tol), not tol

    nxi <- sqrt(sum(xi^2))
    xi <- xi / nxi

    ee <- sum(lam * xi)
    varfun <- charfun(eta = eta,  # log(muvec),
                      extra = list(),  # Object@extra,  # For now
                      varfun = TRUE)
    ss <- sqrt(sum(varfun * xi^2))  # For both poissonff and binomialff
    w.obs <- sum(y0 * xi)
    wtil <- (w.obs - ee) / ss
    if (is.na(wtil)) {
      prb <- 1  # zz
      prb <- 0  # zz
    } else {
      nrm.prb <- pnorm(wtil)
      prb <- if (wtil < -A.qnorm) 0 else
             if ( A.qnorm < wtil) 1 else {
        prbc <- cdffun.integrand(Object, A.qnorm = A.qnorm,
                                 wtil = wtil, lam = lam, xi = xi,
                                 ee = ee, ss = ss,
                                 sigma2.fuzz = sigma2.fuzz)
        Re(prbc)
      }
    }
    prb
}  # charfun.cqo.cdf








charfun.clo.cdf <-
  function(
           bnu,
           y0,
           extra ,
           objfun,
           Object,
           Coefs,
           B1bix,
           misc.list = list(),
           Everything,
           mu.function,


           A.qnorm = 4,
           sigma2.fuzz = 1e-8, lower.latvar = NULL, upper.latvar = NULL,
           nodes = sqrt(2) * ghn100,
           qwts  = sqrt(2) * ghw100 / (2 * pi)) {
  if (length(bnu) > 1)
    stop("can only handle rank-1 objects")
  vfam <- intersect(Object@family@vfamily, c("binomialff", "poissonff"))
  if (!length(vfam))
    stop("only 'poissonff' and 'binomialff' families allowed")
  all.links <- linkfun(Object)
  parlink <- all.links[1]  # All the same, choose the first
  canon.link <- switch(vfam,
                       binomialff = all(all.links == "logitlink"),
                       poissonff  = all(all.links == "loglink"))
  if (!canon.link) stop("model does not use the canonical link")  # else
  A.mat <- Coefs@A
  B1.mat <- Coefs@B1
  Index.corner <- Object@control$Index.corner  # Corner constraints


  wt.mean <- function(x, y, alpha = 0.5) (1 - alpha) * x + alpha * y
  site.scores <- latvar(Object)
  range.ss <- range(site.scores)
  if (is.null(lower.latvar))
    lower.latvar <- wt.mean(range.ss[1], range.ss[2],  0.0)
  if (is.null(upper.latvar))
    upper.latvar <- wt.mean(range.ss[1], range.ss[2],  1.0)

  famfun <- Object@family
  infos <- famfun@infos()
  charfun <- if (is.logical(infos$charfun) && infos$charfun)
    famfun@charfun else stop("the family function has no charfun slot")
  Earg <- Object@misc$earg[[1]]  # Dangerous


  
if (FALSE) {
  uopt <- Opt(Object)
  tol <- Tol(Object)[1, 1, ]  # Interpreted as a variance, not SDs
  parlink <- linkfun(Object)[1]  # All the same, choose the first
  alf <- theta2eta(Max(Object), parlink, earg = list())
}


  







    nu0 <- bnu


    eta <- Coefs@B1["(Intercept)", ] + A.mat %*% nu0
    eta <- rbind(c(eta))  # Make sure it is 1 x M


    muvec <- eta2theta(eta, parlink, earg = Earg)  # Dangerous
    lam <- muvec

    xi <- A.mat[, 1]



    nxi <- sqrt(sum(xi^2))
    xi <- xi / nxi



    ee <- sum(lam * xi)
    varfun <- charfun(eta = eta,  # log(muvec),
                      extra = list(),  # Object@extra,  # For now
                      varfun = TRUE)
    ss <- sqrt(sum(varfun * xi^2))  # For both poissonff and binomialff
    w.obs <- sum(y0 * xi)
    wtil <- (w.obs - ee) / ss
    if (is.na(wtil)) {
      prb <- 1  # zz
      prb <- 0  # zz
    } else {
      nrm.prb <- pnorm(wtil)
      prb <- if (wtil < -A.qnorm) 0 else
             if ( A.qnorm < wtil) 1 else {
        prbc <- cdffun.integrand(Object, A.qnorm = A.qnorm,
                                 wtil = wtil, lam = lam, xi = xi,
                                 ee = ee, ss = ss,
                                 sigma2.fuzz = sigma2.fuzz)
        Re(prbc)
      }
    }
    prb
}  # charfun.clo.cdf
















fnumat2R <- 
  function(object,
           refit.model = FALSE) {



  numat <- latvar(object)  # After scaling
  Rank <- Rank(object)
  control <- object@control

  M <- npred(object)
  M1 <- npred(object, type = "one.response")
  if (M1 > 1)
    stop("this function only works with M1==1 models")
  nsmall <- nobs(object)

  X.lm <- model.matrix(object, type = "lm")  # Has latvar vars too.
  etol <- control$eq.tolerances
  Itol <- control$I.tolerances
  colx1.index <- control$colx1.index
  colx2.index <- control$colx2.index
  NOS <- M / M1

  if (!etol) {
    index.rc <- iam(NA, NA, M = Rank, both = TRUE)
    numat2 <- numat[, index.rc$row.index, drop = FALSE] *
              numat[, index.rc$col.index, drop = FALSE]
    if (Rank > 1)
      numat2[, -(1:Rank)] <- numat2[, -(1:Rank)] * 2
  }  # !etol
  if (etol && !Itol) {
    numat2 <- numat[, 1:Rank, drop = FALSE] *
              numat[, 1:Rank, drop = FALSE]  # Same as Itol
  }
  if (Itol) {
    numat2 <- numat[, 1:Rank, drop = FALSE] *
              numat[, 1:Rank, drop = FALSE]
  }  # Itol



  if (Rank >= 2 && (!etol || (etol && !Itol)))
    stop("cannot currently handle the given scalings")




  ansA <- kronecker(numat, diag(M))
  colnames(ansA) <- param.names("A", NCOL(ansA), skip1 = FALSE)


  ansD <- kronecker(numat2, if (Itol || etol) matrix(1, M, 1) else
                            diag(M)) * (-0.5)
  colnames(ansD) <- param.names("D", NCOL(ansD), skip1 = FALSE)
  ansx1 <- if (length(colx1.index))
           kronecker(X.lm[, colx1.index, drop = FALSE], diag(M)) else
           NULL  # Trivial constraints are assumed.
  if (length(ansx1))
    colnames(ansx1) <- param.names("x1.", NCOL(ansx1), skip1 = FALSE)
  X.vlm <- cbind(ansA, ansD, ansx1)
  if (!refit.model)
    return(X.vlm)


  mf <- model.frame(object)
  mt <- attr(mf, "terms")
  clist <- vector("list", ncol(X.vlm))
  names(clist) <- colnames(X.vlm)
  for (ii in seq(length(clist)))
    clist[[ii]] <- diag(M)  # Wrong but doesnt matter; trivial constraints

  somejunk <- clist
  inci <- 0
  for (ii in seq(length(somejunk))) {
    inci <- inci + 1
    somejunk[[ii]] <- inci
  }
  attr(numat, "assign") <- somejunk
  attr(X.vlm, "assign") <- somejunk
  control$Check.cm.rank <- FALSE
  control$stepsize <- 1
  control$half.stepsizing <- FALSE
  control$Check.rank <- FALSE

  eta.mat.orig <- predict(object)
  OOO.orig <- object@offset
  if (!length(OOO.orig) || all(OOO.orig == 0))
    OOO.orig <- matrix(0, nsmall, M)

  fit1 <-
    vglm.fit(x = numat, y = depvar(object),
             w = c(weights(object, type = "prior")),
             X.vlm.arg = X.vlm,
             Xm2 = NULL, Terms = mt,
             constraints = clist, extra = object@extra,
             etastart = eta.mat.orig,
             offset = OOO.orig, family = object@family,
             control = control)
  if (fit1$iter > 3)
    warning("refitted model took an unusually large number of ",
            "iterations to converge")  # Usually 1 iteration is needed
  return(fit1)
}  # fnumat2R






forG.calib.qrrvglm <-
  function(bnu0,  # coeff3 = NULL, 
           numerator = c("xi", "eta"),
           lenb1bix = 1) {



  numerator <- match.arg(numerator, c("xi", "eta"))[1]

  if (FALSE) {
    if (!is.null(coeff3)) {
      if (length(coeff3) != 3 || coeff3[3] > 0)
        stop("bad input for argument 'coeff3'")
      beta1 <- coeff3[2]  # beta0 <- coeff3[1] is unneeded
      beta2 <- coeff3[3]
    }
    if (is.null(uj))     uj <- -0.5 * beta1 / beta2
    if (is.null(tolj))   tolj <- 1 / sqrt(-2 * beta2)
  }  # FALSE





  Rank <- length(bnu0)  # Usually 1, sometimes 2.
  switch(numerator,
    "xi"  = cbind(rep_len(1, Rank), 2 * bnu0, matrix(0, Rank, lenb1bix)),
    "eta" = cbind(bnu0,               bnu0^2, matrix(1, Rank, lenb1bix)))
}  # forG.calib.qrrvglm








dzwald.qrrvglm <-
  function(bnu0, y0, Object, CoefsObject,
           B1bix, mu.function) {




  M <- npred(Object)
  if ((M1 <- npred(Object, type = "one.response")) > 1)
    stop("this function only works with M1==1 models")
  nsmall <- nobs(Object)
  NOS <- M / M1
  etol <- Object@control$eq.tolerances
  Itol <- Object@control$I.tolerances
  Earg <- Object@misc$earg
  all.links <- linkfun(Object)

  if ((Rank <- Rank(Object)) >= 2 && !Itol)
    warning("can only handle rank-1 (or rank-2 Itol) objects")

  linkfun <- Object@family@linkfun
  if (!is.function(linkfun))
    stop("could not obtain @linkfun")


  vfam <- intersect(Object@family@vfamily, c("binomialff", "poissonff"))
  if (!length(vfam))
    stop("only 'poissonff' and 'binomialff' families allowed")

  canon.link <- switch(vfam,
                       binomialff = all(all.links == "logitlink"),
                       poissonff  = all(all.links == "loglink"))
  if (!canon.link) stop("model does not use the canonical link")  # else
 


  Omegamat <- vcov(Object)
  dimn <- colnames(Omegamat)


  DD <- 0  # Later becomes a matrix.
  Gmat <- matrix(0, ncol(Omegamat), Rank)
  for (spp. in 1:NOS) {
    index.A.spp. <- spp. + (seq(Rank) - 1) * M  # Rank-vector
    index.D.spp. <- if (Itol || etol) (Rank * M) + seq(Rank) else
    if (!etol) (Rank * M) +
      spp. + (seq(Rank * (Rank+1) / 2) - 1) * M  # Rank-vector
    dd <- min(which(substr(dimn, 1, 3) == "x1."))
    index.B1.spp. <- dd - 1 + spp. + (seq(1000) - 1) * M
    index.B1.spp. <- index.B1.spp.[index.B1.spp. <= ncol(Omegamat)]
    all.index <- c(index.A.spp., index.D.spp., index.B1.spp.)






  if (!all(should.be.TRUE <- substr(dimn[index.D.spp.], 1, 1) == "D"))
    stop("encountered a bookkeeping error; failed to pick up the ",
         "D array coefficients")





    uj.jay <- CoefsObject@Optimum[, spp.]  # Rank-vector
    tolj.jay <- if (Rank == 1) sqrt(CoefsObject@Tolerance[, , spp.]) else
                diag(sqrt(CoefsObject@Tolerance[, , spp.]))  # Ditto
    alphaj.jay <- linkfun(CoefsObject@Maximum,
                          extra = Object@extra)[spp.]  # Scalar


    eta0.jay <- alphaj.jay - 0.5 * sum(((bnu0 - uj.jay) / tolj.jay)^2)
    eta0.jay <- rbind(eta0.jay)
    fv0.jay <- mu.function(eta0.jay, extra = Object@extra)
    fv0.jay <- c(fv0.jay)  # Remove array attributes
    dTheta.deta0j <- dtheta.deta(fv0.jay,
                                 all.links[spp.],
                                 earg = Earg[[spp.]])  # Scalar
    dTheta.deta0j <- c(dTheta.deta0j)  # Remove array attributes


    xi0.jay <- -(bnu0 - uj.jay) / tolj.jay^2  # Rank-vector
    DD <- DD + dTheta.deta0j *
        (cbind(xi0.jay) %*% rbind(xi0.jay))  # More general
    dxi0.dtheta <- forG.calib.qrrvglm(bnu0, numerator = "xi")  # R x dim11
    deta0.dtheta <- forG.calib.qrrvglm(bnu0, numerator = "eta")  # Rxdim11
    Index.All <- cbind(index.A.spp., index.D.spp.,
                       rep_len(index.B1.spp., Rank))  # Just to make sure
    for (rlocal in seq(Rank)) {
      Gmat[Index.All[rlocal, ], rlocal] <-
      Gmat[Index.All[rlocal, ], rlocal] +
      (y0[1, spp.] - fv0.jay) * dxi0.dtheta[rlocal, ] -
      xi0.jay[rlocal] * dTheta.deta0j * deta0.dtheta[rlocal, ]
    }  # rlocal

  }  #  for spp.




  DDinv <- solve(DD)
  muxf <- diag(Rank) + t(Gmat) %*% Omegamat %*% Gmat %*% DDinv
  Vmat <- DDinv %*% muxf
  Vmat
}  # dzwald.qrrvglm












calibrate.qrrvglm.control <-
  function(object,
           trace = FALSE,  # passed into optim()
           method.optim = "BFGS",  # passed into optim(method = Method)
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
       method.optim = method.optim,
       gridSize = gridSize,
       varI.latvar = as.logical(varI.latvar)[1])
}  # calibrate.qrrvglm.control





 if (!isGeneric("calibrate"))
     setGeneric("calibrate",
                function(object, ...) standardGeneric("calibrate"))





calibrate.qrrvglm <-
  function(object,
           newdata = NULL,
           type = c("latvar", "predictors", "response", "vcov",
                    "everything"),
           lr.confint = FALSE,  # 20180430
           cf.confint = FALSE,  # 20180602
           level = 0.95,        # 20180430
           initial.vals = NULL,
           ...) {




      

  se.type <- c("dzwald", "asbefore")  # Effectively only the 1st one used
  Quadratic <- if (is.logical(object@control$Quadratic))
               object@control$Quadratic else FALSE  # T if CQO, F if CAO


  newdata.orig <- newdata
  if (!length(newdata)) {
    newdata <- data.frame(depvar(object))
  }


  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("latvar", "predictors",
                            "response", "vcov", "everything"))[1]




  get.SEs <- Quadratic && type %in% c("vcov", "everything")



  if (mode(se.type) != "character" && mode(se.type) != "name")
    se.type <- as.character(substitute(se.type))
  se.type <- match.arg(se.type, c("dzwald", "asbefore"))[1]
  if (se.type == "asbefore") warning("'asbefore' is buggy")

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
  minimize.obfunct <-
    if (Quadratic) object@control$min.criterion else
    TRUE  # Logical; TRUE for CAO objects because deviance is minimized
  if (!is.logical(minimize.obfunct))
    stop("object@control$min.criterion is not a logical")
  optim.control <- calibrate.qrrvglm.control(object = object, ...)

  use.optim.control <- optim.control
  use.optim.control$method.optim <-
  use.optim.control$gridSize     <-
  use.optim.control$varI.latvar  <- NULL


  if ((Rank <- object@control$Rank) > 2)
    stop("currently can only handle Rank = 1 and Rank = 2")
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
    if (se.type == "dzwald" && (type == "everything" || type == "vcov"))
      stop("only noRRR = ~ 1 models are handled for ",
           "type = 'everything' or type = 'vcov'")

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
              Everything = FALSE,
              mu.function = extraargs$mu.function))
    ans
  }  # objfun1

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
              Everything = FALSE,
              mu.function = extraargs$mu.function))
    ans
  }  # objfun2


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
         Object       = if (Quadratic) 777 else object,
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
            method = optim.control$method.optim,  # "BFGS" or "CG" or...
            control = c(fnscale = ifelse(minimize.obfunct, 1, -1),
                        use.optim.control),  # as.vector() needed
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              extra = object@extra,
              objfun = obfunct,  # deviance
              Object = if (Quadratic) 777 else object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              Everything = FALSE,  # Differs from Step 6. below
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







  prettyCQO <- function(BestOFpar, newdata, Rank) {
    if (Rank == 1) {
      BestOFpar <- c(BestOFpar)  # Drop the dimension
      if (!is.null(dimnames(newdata)[[1]])) {
        names(BestOFpar) <- dimnames(newdata)[[1]]
      }
    } else {
      dimnames(BestOFpar) <- list(dimnames(newdata)[[1]],
                             param.names("latvar", Rank, skip1 = TRUE))
    }
    BestOFpar
  }  # prettyCQO

  BestOFpar <- prettyCQO(BestOFpar, newdata, Rank)
  attr(BestOFpar,"objectiveFunction") <-
    prettyCQO(BestOFvalues, newdata, Rank = 1)
  if (type == "latvar" && (!cf.confint && !lr.confint)) {
    return(BestOFpar)
  }



  
   


  if (lr.confint && Rank > 1) {
    warning("argument 'lr.confint' should only be TRUE if Rank == 1. ",
            "Setting 'lr.confint = FALSE'.")
    lr.confint <- FALSE
  }


  if (lr.confint && !(type %in% c("latvar", "everything"))) {
    warning("argument 'lr.confint' should only be TRUE if ",
            "'type = \"latvar\"' or 'type = \"everything\"'. ",
            "Setting 'lr.confint = FALSE'.")
    lr.confint <- FALSE
  }
  if (lr.confint && Rank == 1) {

    format.perc <- function(probs, digits)
      paste(format(100 * probs, trim = TRUE, scientific = FALSE,
            digits = digits), "%")
    aa <- (1 - level) / 2
    aa <- c(aa, 1 - aa)

    pct <- format.perc(aa, 3)
    cimat <- array(NA, dim = c(nn, 2L),
                   dimnames = list(dimnames(newdata)[[1]], pct))

    for (i1 in 1:nn) {
      if (optim.control$trace) {
        cat("\nSolving for the roots for obsn", i1, "---------------\n")
        flush.console()
      }


      foo1.lhs.rhs <-
      function(bnu,
               y0, extra = NULL,
               objfun, Coefs,
               Object,  # Not needed
               B1bix,
               misc.list,
               Everything = FALSE,
               mu.function,
               BestOFvalues = NA,
               level = 0.95,
               criterion.arg = "loglikelihood") {
        if (!(criterion.arg %in% c("deviance", "loglikelihood")))
          stop("'criterion.arg' must be 'deviance' or 'loglikelihood'")
        ifelse(criterion.arg == "deviance", 1, -2) *
        (-BestOFvalues +
             .my.calib.objfunction.qrrvglm(bnu = bnu,
                y0 = y0,
                extra = extra,
                objfun = objfun,
                Object = Object,
                Coefs = Coefs,
                B1bix = B1bix,
                misc.list = misc.list,
                Everything = Everything,
                mu.function = mu.function)) -
        qchisq(level, df = 1)
      }


      for (Side in 1:2) {
        ans.lhs.rhs <-
        uniroot(f = foo1.lhs.rhs,
                interval = if (Side == 1) c(Lvec, BestOFpar[i1]) else
                                          c(BestOFpar[i1], Uvec),
                extendInt = ifelse(Side == 1, "downX", "upX"),
                y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
                extra = object@extra,
                objfun = obfunct,
                Coefs = Coefobject,
                Object = 777,  # object,
                B1bix = B1bix[i1, , drop = FALSE],
                misc.list = object@misc,
                Everything = FALSE,
                mu.function = mu.function,
                BestOFvalues = BestOFvalues[i1], level = level,
                criterion.arg = object@misc$criterion)
        cimat[i1, Side] <- ans.lhs.rhs$root
      }  # Side
    }  # for i1

    if (type == "latvar")
      return(cbind(estimate = BestOFpar,
                   objfun   = BestOFvalues,
                   cimat))
  }  # if (lr.confint && Rank == 1)




  
   



  if (cf.confint && Rank > 1) {
    warning("argument 'cf.confint' should only be TRUE if Rank == 1. ",
            "Setting 'cf.confint = FALSE'.")
    cf.confint <- FALSE
  }


  if (cf.confint && !(type %in% c("latvar", "everything"))) {
    warning("argument 'cf.confint' should only be TRUE if ",
            "'type = \"latvar\"' or 'type = \"everything\"'. ",
            "Setting 'cf.confint = FALSE'.")
    cf.confint <- FALSE
  }
  if (cf.confint && Rank == 1) {

    format.perc <- function(probs, digits)
      paste(format(100 * probs, trim = TRUE, scientific = FALSE,
            digits = digits), "%")
    aa <- (1 - level) / 2
    aa <- c(aa, 1 - aa)
    pct <- format.perc(aa, 3)
    cimat2 <- array(NA, dim = c(nn, 2L),
                    dimnames = list(dimnames(newdata)[[1]], pct))



    for (i1 in 1:nn) {
      if (optim.control$trace) {
        cat("\nSolving for the roots for obsn", i1, "---------------\n")
        flush.console()
      }


      foo2.lhs.rhs <-
      function(bnu,
               y0, extra = NULL,
               objfun, Coefs,
               Object,  # Not needed
               B1bix,
               misc.list,
               Everything = FALSE,
               mu.function,
               BestOFvalues = NA,
               pr.level = c(0.05, 0.95)[1]
               ) {
    charfun.cqo.cdf(bnu = bnu,
                y0 = y0,
                extra = extra,
                objfun = objfun,
                Coefs = Coefs,
                Object = Object,
                B1bix = B1bix,
                misc.list = misc.list,
                Everything = Everything,
                mu.function = mu.function
                ) -
        pr.level
      }





      

      for (Side in 1:2) {
        ans.lhs.rhs <-
        uniroot(f = foo2.lhs.rhs,
                interval = if (Side == 1) c(Lvec, BestOFpar[i1]) else
                                          c(BestOFpar[i1], Uvec),
                extendInt = "yes",  # Might be worse than above.
                y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
                extra = object@extra,
                objfun = obfunct,
                Coefs = Coefobject,
                Object = object,
                B1bix = B1bix[i1, , drop = FALSE],
                misc.list = object@misc,
                Everything = FALSE,
                mu.function = mu.function,
                BestOFvalues = BestOFvalues[i1],
                pr.level = ifelse(Side == 1, aa[1], aa[2])
               )
        cimat2[i1, Side] <- ans.lhs.rhs$root
      }  # Side
    }  # for i1

    vecTF <- cimat2[, 2] < cimat2[, 1]
    if (any(vecTF)) {
      temp <- cimat2[vecTF, 1]
      cimat2[vecTF, 1] <- cimat2[vecTF, 2]
      cimat2[vecTF, 2] <- temp
    }
    if (type == "latvar")
      return(cbind(estimate = BestOFpar,
                   objfun   = BestOFvalues,
                   cimat2))

  }  # if (cf.confint && Rank == 1)




      



  etaValues <- matrix(NA_real_, nn, M)
  muValues <- matrix(NA_real_, nn, ncol(fitted(object)))
  vcValues <- if (get.SEs) array(0, c(Rank, Rank, nn)) else NULL


  if (optim.control$trace)
    cat("\n")

  for (i1 in 1:nn) {
    if (optim.control$trace) {
      cat("Evaluating quantities for observation", i1,
          "-----------------\n")
      flush.console()
    }
    ans5 <- choose.fun(
              bnu = if (Rank == 1) BestOFpar[i1] else BestOFpar[i1, ],
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              extra = object@extra,
              objfun = obfunct,  # deviance
              Object = if (Quadratic) 777 else object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              Everything = TRUE,  # Differs from Step 3.
              mu.function = mu.function)

     muValues[i1, ] <- ans5$mu
    etaValues[i1, ] <- ans5$eta


    if (get.SEs) {
      vcValues[, , i1] <- if (se.type == "dzwald")
        dzwald.qrrvglm(
              bnu0 = if (Rank == 1) BestOFpar[i1] else BestOFpar[i1, ],
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              Object = object,
              CoefsObject = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              mu.function = mu.function) else
        ans5$vcmat  # Might be NULL, e.g., "rrvgam"
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
         latvar     = BestOFpar,  # Done already, so not really needed
         predictors = etaValues,
         response   = muValues,
         vcov       = vcValues,
         everything = list(
             latvar     = BestOFpar,
             predictors = etaValues,
             response   = muValues,
             vcov       = vcValues,
             lr.confint = if (lr.confint)
                            cbind(estimate = BestOFpar,
                                  objfun   = BestOFvalues,
                                  cimat) else NULL,
             cf.confint = if (cf.confint)
                            cbind(estimate = BestOFpar,
                                  objfun   = BestOFvalues,
                                  cimat2) else NULL)
        )
}  # calibrate.qrrvglm











.my.calib.objfunction.qrrvglm <-
  function(bnu, y0, extra = NULL,
           objfun, Coefs, Object,
           B1bix,
           misc.list,
           Everything = TRUE,
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
      alf <- loglink(Coefs@Maximum[ss])  # zz get the link function
      tolmat <- Coefs@Tolerance[, , ss, drop = FALSE]
      check.eta[ss, 1] <- alf - 0.5 * t(bnumat) %*%
                          solve(tolmat) %*% bnumat  
    }  # FALSE
    }  # for ss
  eta <- matrix(eta, 1, M, byrow = TRUE)
  mu <- matrix(mu.function(eta, extra = extra), nrow = 1)
  obvalue <- objfun(mu = mu, y = y0,
                    w = 1,  # ignore prior.weights on the object
                    residuals = FALSE, eta = eta, extra = extra)
  if (Everything) {
    vcmat <- diag(Rank)
    if (FALSE && M == NCOL(mu)) {
      for (ss in 1:M) {
        vec1 <- cbind(Coefs@A[ss, ]) +
                      2 * matrix(Coefs@D[, , ss], Rank, Rank) %*% bnumat
        vcmat <- vcmat + mu[1, ss] * vec1 %*% t(vec1)
      }  # ss
    }  # if (M == NCOL(mu))
    vcmat <- solve(vcmat)
  } else {
    vcmat <- NULL
  }
  if (Everything)
    list(eta     = eta,
         mu      = mu,
         obvalue = obvalue,
         vcmat   = vcmat) else
    obvalue
}  # .my.calib.objfunction.qrrvglm






.my.calib.objfunction.rrvgam <-
  function(bnu, y0, extra = NULL,
           objfun,
           Object,  # Needed for "rrvgam" objects
           Coefs,
           B1bix,  # Actually not needed here
           misc.list,
           Everything = TRUE,
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
  if (Everything)
    list(eta = eta,
         mu = mu,
         obvalue = obvalue,
         vcmat = vcmat) else
    obvalue
}  # .my.calib.objfunction.rrvgam













forG.calib.rrvglm <-
  function(bnu0,
           numerator = c("xi", "eta"),
           lenb1bix = 1) {

  numerator <- match.arg(numerator, c("xi", "eta"))[1]

  Rank <- length(bnu0)  # Usually 1, sometimes 2.
  switch(numerator,
    "xi"  = cbind(rep_len(1, Rank), matrix(0, Rank, lenb1bix)),
    "eta" = cbind(bnu0,             matrix(1, Rank, lenb1bix)))
}  # forG.calib.rrvglm






dzwald.rrvglm <-
  function(bnu0, y0,  # extra = NULL, objfun,
           Object,  CoefsObject,
           B1bix,
           mu.function
           ) {

  
  M <- npred(Object)
  if ((M1 <- npred(Object, type = "one.response")) > 1)
    stop("this function only works with M1==1 models")
  NOS <- M / M1
  Earg <- Object@misc$earg
  all.links <- linkfun(Object)




  Rank <- Rank(Object)

  linkfun <- Object@family@linkfun
  if (!is.function(linkfun))
    stop("could not obtain @linkfun")


  vfam <- intersect(Object@family@vfamily, c("binomialff", "poissonff"))
  if (!length(vfam))
    stop("only 'poissonff' and 'binomialff' families allowed")

  canon.link <- switch(vfam,
                       binomialff = all(all.links == "logitlink"),
                       poissonff  = all(all.links == "loglink"))
  if (!canon.link) stop("model does not use the canonical link")  # else
 


  Omegamat <- vcov(Object)  # Numerical problems might occur to get this.
  dimn <- colnames(Omegamat)
  A.mat <- CoefsObject@A
  B1.mat <- CoefsObject@B1
  Index.corner <- Object@control$Index.corner  # Corner constraints

  DD <- 0  # Later becomes a matrix.
  Gmat <- matrix(0, ncol(Omegamat), Rank)
  icounter <- 0  # Number of rows of \bigtilde{\bA}.
  for (spp. in 1:NOS) {
    index.A.spp. <- if (any(spp. == Index.corner)) NULL else {
      icounter <- icounter + 1
      icounter + (M - Rank) * (seq(Rank) - 1)
    }
    index.D.spp. <- NULL
    dd <- max(which(substr(dimn, 1, 8) == "I(latvar"))
    index.B1.spp. <- dd + spp. + (seq(nrow(B1.mat)) - 1) * M


    all.index <- c(index.A.spp., index.D.spp., index.B1.spp.)







    alphaj.jay <- CoefsObject@B1["(Intercept)", spp.]  # Scalar


    eta0.jay <- alphaj.jay + A.mat[spp., , drop = FALSE] %*% bnu0
    eta0.jay <- rbind(eta0.jay)
    fv0.jay <- mu.function(eta0.jay, extra = Object@extra)
    fv0.jay <- c(fv0.jay)  # Remove array attributes
    dTheta.deta0j <- dtheta.deta(fv0.jay,
                                 all.links[spp.],
                                 earg = Earg[[spp.]])  # Scalar
    dTheta.deta0j <- c(dTheta.deta0j)  # Remove array attributes


    xi0.jay <- A.mat[spp., ]  # Rank-vector
    DD <- DD + dTheta.deta0j *
        (cbind(xi0.jay) %*% rbind(xi0.jay))  # More general
    dxi0.dtheta <- forG.calib.rrvglm(bnu0, numerator = "xi")  # R x dim11
    deta0.dtheta <- forG.calib.rrvglm(bnu0, numerator = "eta")  # Rxdim11
    Index.All <- cbind(index.A.spp., index.D.spp.,
                       rep_len(index.B1.spp., Rank))  # Just to make sure

    if (!is.null(index.A.spp.))
    for (rlocal in seq(Rank)) {
      Gmat[Index.All[rlocal, ], rlocal] <-
      Gmat[Index.All[rlocal, ], rlocal] +
      (y0[1, spp.] - fv0.jay) * dxi0.dtheta[rlocal, ] -
      xi0.jay[rlocal] * dTheta.deta0j * deta0.dtheta[rlocal, ]
    }  # rlocal

  }  #  for spp.




  DDinv <- solve(DD)
  muxf <- diag(Rank) + t(Gmat) %*% Omegamat %*% Gmat %*% DDinv
  Vmat <- DDinv %*% muxf
  Vmat
}  # dzwald.rrvglm









calibrate.rrvglm <-
  function(object,
           newdata = NULL,
           type = c("latvar", "predictors", "response", "vcov",
                    "everything"),
           lr.confint = FALSE,  # 20180427
           cf.confint = FALSE,  # 20180604
           level = 0.95,        # 20180428
           initial.vals = NULL,  # For one observation only
           ...) {
  se.type <- c("dzwald", "asbefore")  # Effectively only the 1st one used


  Quadratic <- FALSE  # Because this function was adapted from CQO code.
  newdata.orig <- newdata
  if (!length(newdata)) {
    newdata <- data.frame(depvar(object))
  }

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("latvar", "predictors",
                            "response", "vcov", "everything"))[1]
  get.SEs <- type %in% c("vcov", "everything")
  if (mode(se.type) != "character" && mode(se.type) != "name")
    se.type <- as.character(substitute(se.type))
  se.type <- match.arg(se.type, c("dzwald", "asbefore"))[1]

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
  minimize.obfunct <- object@control$min.criterion  # deviance
  if (!is.logical(minimize.obfunct))
    stop("object@control$min.criterion is not a logical")
  minimize.obfunct <- as.vector(minimize.obfunct)
  optim.control <- calibrate.rrvglm.control(object = object, ...)

  use.optim.control <- optim.control
  use.optim.control$method.optim <-
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
    Xlm <- predict.vlm(as(object, "vglm"),  # object,
                       newdata = newdata.orig,
                       type = "Xlm")
    if (NROW(Xlm) != nn)
      warning("NROW(Xlm) and ", nn, " are unequal")

    if (se.type == "dzwald" && (type == "everything" || type == "vcov"))
      stop("only noRRR = ~ 1 models are handled for ",
           "type = 'everything' or type = 'vcov'")

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
              Everything = FALSE,
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
              Everything = FALSE,
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
              Everything = FALSE,
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
         Object       = 777,  # object,
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
            method = optim.control$method.optim,  # "BFGS" or "CG" or...
            control = c(fnscale = ifelse(minimize.obfunct, 1, -1),
                        use.optim.control),  # as.vector() needed
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              extra = object@extra,
              objfun = obfunct,
              Object = 777,  # object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              Everything = FALSE,  # Differs from Step 5 below
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







  prettyCLO <- function(BestOFpar, newdata, Rank) {
    if (Rank == 1) {
      BestOFpar <- c(BestOFpar)  # Drop the dimension
      if (!is.null(dimnames(newdata)[[1]])) {
        names(BestOFpar) <- dimnames(newdata)[[1]]
      }
    } else
      dimnames(BestOFpar) <- list(dimnames(newdata)[[1]],
                             param.names("latvar", Rank, skip1 = TRUE))
    BestOFpar
  }  # prettyCLO

  BestOFpar <- prettyCLO(BestOFpar, newdata, Rank)  # Dimension may drop.
  attr(BestOFpar,"objectiveFunction") <-
    prettyCLO(BestOFvalues, newdata, Rank = 1)
  if (type == "latvar" && (!cf.confint && !lr.confint)) {
    return(BestOFpar)
  }






  if (lr.confint && Rank > 1) {
    warning("argument 'lr.confint' should only be TRUE if Rank == 1. ",
            "Setting 'lr.confint = FALSE'.")
    lr.confint <- FALSE
  }


  if (lr.confint && !(type %in% c("latvar", "everything"))) {
    warning("argument 'lr.confint' should only be TRUE if ",
            "'type = \"latvar\"' or 'type = \"everything\"'. ",
            "Setting 'lr.confint = FALSE'.")
    lr.confint <- FALSE
  }
  if (lr.confint && Rank == 1) {

    format.perc <- function(probs, digits)
      paste(format(100 * probs, trim = TRUE, scientific = FALSE,
            digits = digits), "%")
    aa <- (1 - level) / 2
    aa <- c(aa, 1 - aa)

    pct <- format.perc(aa, 3)
    cimat <- array(NA, dim = c(nn, 2L),
                   dimnames = list(dimnames(newdata)[[1]], pct))

    for (i1 in 1:nn) {
      if (optim.control$trace) {
        cat("Solving for the roots for obsn", i1, "---------------\n")
        flush.console()
      }


      foo3.lhs.rhs <-
      function(bnu,
               y0, extra = NULL,
               objfun, Coefs,
               Object,  # Not needed
               B1bix,
               misc.list,
               Everything = FALSE,
               mu.function,
               BestOFvalues = NA,
               level = 0.95,
               criterion.arg = "loglikelihood") {
        if (!(criterion.arg %in% c("deviance", "loglikelihood")))
          stop("'criterion.arg' must be 'deviance' or 'loglikelihood'")
        ifelse(criterion.arg == "deviance", 1, -2) *
        (-BestOFvalues +
             .my.calib.objfunction.rrvglm(bnu = bnu,
                y0 = y0,
                extra = extra,
                objfun = objfun,
                Object = Object,
                Coefs = Coefs,
                B1bix = B1bix,
                misc.list = misc.list,
                Everything = Everything,
                mu.function = mu.function)) -
        qchisq(level, df = 1)
      }


      for (Side in 1:2) {
        ans.lhs.rhs <-
        uniroot(f = foo3.lhs.rhs,
                interval = if (Side == 1) c(Lvec, BestOFpar[i1]) else
                                          c(BestOFpar[i1], Uvec),
                extendInt = ifelse(Side == 1, "downX", "upX"),
                y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
                extra = object@extra,
                objfun = obfunct,
                Object = 777,  # object,
                Coefs = Coefobject,
                B1bix = B1bix[i1, , drop = FALSE],
                misc.list = object@misc,
                Everything = FALSE,
                mu.function = mu.function,
                BestOFvalues = BestOFvalues[i1], level = level,
                criterion.arg = object@misc$criterion)
        cimat[i1, Side] <- ans.lhs.rhs$root
      }  # Side
    }  # for i1

    if (type == "latvar")
      return(cbind(estimate = BestOFpar,
                   objfun   = BestOFvalues,
                   cimat))
  }  # if (lr.confint && Rank == 1)


  
   



  if (cf.confint && Rank > 1) {
    warning("argument 'cf.confint' should only be TRUE if Rank == 1. ",
            "Setting 'cf.confint = FALSE'.")
    cf.confint <- FALSE
  }


  if (cf.confint && !(type %in% c("latvar", "everything"))) {
    warning("argument 'cf.confint' should only be TRUE if ",
            "'type = \"latvar\"' or 'type = \"everything\"'. ",
            "Setting 'cf.confint = FALSE'.")
    cf.confint <- FALSE
  }
  if (cf.confint && Rank == 1) {

    format.perc <- function(probs, digits)
      paste(format(100 * probs, trim = TRUE, scientific = FALSE,
            digits = digits), "%")
    aa <- (1 - level) / 2
    aa <- c(aa, 1 - aa)
    pct <- format.perc(aa, 3)
    cimat2 <- array(NA, dim = c(nn, 2L),
                    dimnames = list(dimnames(newdata)[[1]], pct))



    for (i1 in 1:nn) {
      if (optim.control$trace) {
        cat("\nSolving for the roots for obsn", i1, "---------------\n")
        flush.console()
      }


      foo4.lhs.rhs <-
      function(bnu,
               y0, extra = NULL,
               objfun, Coefs,
               Object,  # Not needed
               B1bix,
               misc.list,
               Everything = FALSE,
               mu.function,
               BestOFvalues = NA,
               pr.level = c(0.05, 0.95)[1]
               ) {
    charfun.clo.cdf(bnu = bnu,
                y0 = y0,
                extra = extra,
                objfun = objfun,
                Object = Object,
                Coefs = Coefs,
                B1bix = B1bix,
                misc.list = misc.list,
                Everything = Everything,
                mu.function = mu.function
                ) -
        pr.level
      }





      

      for (Side in 1:2) {
        ans.lhs.rhs <-
        uniroot(f = foo4.lhs.rhs,
                interval = if (Side == 1) c(Lvec, BestOFpar[i1]) else
                                          c(BestOFpar[i1], Uvec),
                extendInt = "yes",  # Might be worse than above.
                y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
                extra = object@extra,
                objfun = obfunct,
                Object = object,
                Coefs = Coefobject,
                B1bix = B1bix[i1, , drop = FALSE],
                misc.list = object@misc,
                Everything = FALSE,
                mu.function = mu.function,
                BestOFvalues = BestOFvalues[i1],
                pr.level = ifelse(Side == 1, aa[1], aa[2])
               )
        cimat2[i1, Side] <- ans.lhs.rhs$root
      }  # Side
    }  # for i1

    vecTF <- cimat2[, 2] < cimat2[, 1]
    if (any(vecTF)) {
      temp <- cimat2[vecTF, 1]
      cimat2[vecTF, 1] <- cimat2[vecTF, 2]
      cimat2[vecTF, 2] <- temp
    }
    if (type == "latvar")
      return(cbind(estimate = BestOFpar,
                   objfun   = BestOFvalues,
                   cimat2))

  }  # if (cf.confint && Rank == 1)








  etaValues <- matrix(NA_real_, nn, M)
  muValues  <- matrix(NA_real_, nn, ncol(fitted(object)))
  vcValues <- if (get.SEs) array(0, c(Rank, Rank, nn)) else NULL


  if (optim.control$trace)
    cat("\n")

  for (i1 in 1:nn) {
    if (optim.control$trace) {
      cat("Evaluating quantities for observation", i1,
          "-----------------\n")
      flush.console()
    }
    ans5 <- .my.calib.objfunction.rrvglm(
              bnu = if (Rank == 1) BestOFpar[i1] else BestOFpar[i1, ],
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              extra = object@extra,
              objfun = obfunct,  # deviance
              Object = 777,  # object,
              Coefs = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              misc.list = object@misc,
              Everything = TRUE,  # Differs from Step 3.
              mu.function = mu.function)


    muValues[i1, ] <- ans5$mu
    etaValues[i1, ] <- ans5$eta


    if (get.SEs)
      vcValues[, , i1] <- if (se.type == "dzwald")
        dzwald.rrvglm(
              bnu0 = if (Rank == 1) BestOFpar[i1] else BestOFpar[i1, ],
              y0 = newdata[i1, , drop = FALSE],  # drop added 20150624
              Object = object,
              CoefsObject = Coefobject,
              B1bix = B1bix[i1, , drop = FALSE],
              mu.function = mu.function
            ) else
        ans5$vcmat



  }  # for i1







  dimnames(muValues) <- dimnames(newdata)
  dimnames(etaValues) <- list(dimnames(newdata)[[1]],
                              dimnames(object@predictors)[[2]])
  if (get.SEs)
    dimnames(vcValues) <- list(as.character(1:Rank),
                               as.character(1:Rank),
                               dimnames(newdata)[[1]])

  switch(type,
         latvar     = BestOFpar,  # Done already, so not really needed
         predictors = etaValues,
         response   = muValues,
         vcov       = vcValues,
         everything = list(
             latvar     = BestOFpar,
             predictors = etaValues,
             response   = muValues,
             vcov       = vcValues,
             lr.confint = if (lr.confint)
                            cbind(estimate = BestOFpar,
                                  objfun   = BestOFvalues,
                                  cimat) else NULL,
             cf.confint = if (cf.confint)
                            cbind(estimate = BestOFpar,
                                  objfun   = BestOFvalues,
                                  cimat2) else NULL)
        )
}  # calibrate.rrvglm






.my.calib.objfunction.rrvglm <-
  function(bnu,
           y0, extra = NULL,
           objfun, Coefs,
           Object,  # Not needed
           B1bix,
           misc.list,
           Everything = TRUE,
           mu.function) {



  bnumat <- cbind(bnu)
  Rank <- length(bnu)
  eta <- cbind(c(B1bix)) + Coefs@A %*% bnumat
  M <- misc.list$M
  eta <- matrix(eta, 1, M, byrow = TRUE)

  mu <- matrix(mu.function(eta, extra = extra), nrow = 1)
  obvalue <- objfun(mu = mu, y = y0,
                    w = 1,  # ignore prior.weights on the object zz
                    residuals = FALSE, eta = eta, extra = extra)
  if (Everything) {
    vcmat <- diag(Rank)
    vcmat <- solve(vcmat)
  } else {
    vcmat <- NULL
  }
  if (Everything)
    list(eta     = eta,
         mu      = mu,
         obvalue = obvalue,
         vcmat   = vcmat) else
    obvalue
}  # .my.calib.objfunction.rrvglm






calibrate.rrvglm.control <-
  function(object,
           trace = FALSE,  # passed into optim()
           method.optim = "BFGS",  # passed into optim(method = Method)
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
       method.optim = method.optim,
       gridSize = gridSize
      )
}  # calibrate.rrvglm.control





setMethod("calibrate", "rrvglm", function(object, ...)
          calibrate.rrvglm(object, ...))



  setMethod("calibrate", "qrrvglm", function(object, ...)
           calibrate.qrrvglm(object, ...))




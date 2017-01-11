# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.












cao.fit <-
  function(x, y, w = rep_len(1, length(x[, 1])),
           etastart = NULL, mustart = NULL, coefstart = NULL,
           offset = 0, family,
           control = cao.control(...), criterion = "coefficients",
           qr.arg = FALSE, constraints = NULL, extra = NULL,
           Terms = Terms, function.name = "cao", ...) {


  maxitl <- NULL
  fv <- NULL


  eff.n <- nrow(x)  # + sum(abs(w[1:nrow(x)]))

  specialCM <- NULL
  post <- list()
  check.rank <- TRUE
  nonparametric <- TRUE
  optim.maxit <- control$optim.maxit
  save.weights <- control$save.weights
  trace <- control$trace
  minimize.criterion <- control$min.criterion

  n <- dim(x)[1]


  copy.X.vlm <- FALSE  # May be overwritten in @initialize

  X.vlm.save <- NULL

  intercept.only <- ncol(x) == 1 && dimnames(x)[[2]] == "(Intercept)"
  y.names <- predictors.names <- NULL # May be overwritten in @initialize


  n.save <- n


  Rank <- control$Rank
  rrcontrol <- control  #

  if (length(family@initialize))
    eval(family@initialize)   # Initialize mu and M (and optionally w)
  n <- n.save

  modelno <- switch(family@vfamily[1], "poissonff" = 2,
                    "binomialff" = 1, "quasipoissonff" = 0,
                    "quasibinomialff" = 0, "negbinomial" = 3,
                    "gamma2" = 5, "gaussianff" = 8,
                    0)  # stop("cannot fit this model using fast algorithm")
  if (!modelno)
    stop("the family function does not work with cao()")
  if (modelno == 1)
    modelno <- get("modelno", envir = VGAMenv)

  eval(rrr.init.expression)

  if (length(etastart)) {
    eta <- etastart
    mu <- if (length(mustart)) mustart else
          if (length(body(slot(family, "linkinv"))))
            slot(family, "linkinv")(eta, extra) else
            warning("argument 'etastart' assigned a value ",
                    "but there is no 'linkinv' slot to use it")
  }

  if (length(mustart)) {
    mu <- mustart
    if (length(body(slot(family, "linkfun")))) {
      eta <- slot(family, "linkfun")(mu, extra)
    } else {
      warning("argument 'mustart' assigned a value ",
              "but there is no 'link' slot to use it")
    }
  }


  M <- if (is.matrix(eta)) ncol(eta) else 1



  if (length(family@constraints))
    eval(family@constraints)


  special.matrix <- matrix(-34956.125, M, M)  # An unlikely used matrix
  just.testing <- cm.VGAM(special.matrix, x, rrcontrol$noRRR, constraints)
  findex <- trivial.constraints(just.testing, special.matrix)
  tc1 <- trivial.constraints(constraints)


  if (all(findex == 1))
    stop("No covariates to form latent variables from.")
  colx1.index <- names.colx1.index <- NULL
  dx2 <- dimnames(x)[[2]]
  if (sum(findex)) {
    asx <- attr(x, "assign")
    for (ii in names(findex))
      if (findex[ii]) {
        names.colx1.index <- c(names.colx1.index, dx2[asx[[ii]]])
        colx1.index <- c(colx1.index, asx[[ii]])
      }
    names(colx1.index) <- names.colx1.index
  }
  rrcontrol$colx1.index <- control$colx1.index <- colx1.index
  colx2.index <- 1:ncol(x)
  names(colx2.index) <- dx2
  colx2.index <- colx2.index[-colx1.index]
  p1 <- length(colx1.index)
  p2 <- length(colx2.index)
  rrcontrol$colx2.index <- control$colx2.index <- colx2.index



  Cmat <- if (length(rrcontrol$Cinit)) {
            matrix(rrcontrol$Cinit, p2, Rank)
          } else {
            if (!rrcontrol$Use.Init.Poisson.QO) {
              matrix(rnorm(p2 * Rank, sd = rrcontrol$sd.Cinit), p2, Rank)
            } else {
                .Init.Poisson.QO(ymat = as.matrix(y),
                          X1 = x[, colx1.index, drop = FALSE],
                          X2 = x[, colx2.index, drop = FALSE],
                          Rank = rrcontrol$Rank, trace = rrcontrol$trace,
                          max.ncol.etamat = rrcontrol$Etamat.colmax,
                          Crow1positive = rrcontrol$Crow1positive,
                          constwt = any(family@vfamily[1] ==
                          c("negbinomial", "gamma2", "gaussianff")),
                  takelog = any(family@vfamily[1] != c("gaussianff")))
            }
          }


  rrcontrol$Cinit <- control$Cinit <- Cmat  # Good for valt()

  Hlist <- process.constraints(constraints, x, M, specialCM = specialCM)

  nice31 <- checkCMCO(Hlist, control = control, modelno = modelno)
  if (nice31 != 1)
    stop("not nice")

  ncolHlist <- unlist(lapply(Hlist, ncol))
  latvar.mat <- x[, colx2.index, drop = FALSE] %*% Cmat


  rmfromVGAMenv(c("etamat", "beta"), prefix = ".VGAM.CAO.")

  Nice21 <- length(names.colx1.index) == 1 &&
            names.colx1.index == "(Intercept)"
  if (!Nice21)
    stop("'noRRR = ~ 1' is supported only, without constraints")
  NOS <- ifelse(modelno %in% c(3, 5), M/2, M)
  p1star. <- if (Nice21) ifelse(modelno %in% c(3, 5), 2, 1) else M
  p2star. <- if (Nice21) Rank else stop("not Nice21")
  pstar. <- p1star. + p2star.
  nstar <- if (Nice21) ifelse(modelno %in% c(3, 5), n * 2, n) else n * M
  lenbeta <- pstar. * ifelse(Nice21, NOS, 1)

  othint <-
        c(Rank, control$eq.tol, pstar. ,
                 dim2wz = 1, inited = 0,  # w(, dimw) cols
          modelno, maxitl = control$maxitl,
          actnits = 0, twice = 0, p1star. ,
          p2star. , Nice21, lenbeta,
          controlI.tolerances = 0, control$trace,
          p1, p2 = p2, imethod = control$imethod, bchat = 0)
  othdbl <- c(small = control$SmallNo, fseps = control$epsilon,
              .Machine$double.eps,
              iKvector = rep_len(control$iKvector, NOS),
              iShape   = rep_len(control$iShape,   NOS),
              resss = 0, bfeps = control$bf.epsilon, hstep = 0.1)

  for (iter in 1:optim.maxit) {
    if (control$trace) {
      cat("\nIteration", iter, "\n")
      flush.console()
    }

      conjgrad <- optim(par = c(Cmat), fn = callcaoc,
                   gr = if (control$GradientFunction) calldcaoc else NULL,
                   method = "BFGS",
                   control = list(fnscale = 1,
                                  trace = as.integer(control$trace),
                                  maxit = control$Maxit.optim,
                                  REPORT = 10),
                   etamat = eta, xmat = x, ymat = y,  # as.matrix(y),
                   wvec = w, modelno = modelno,
                   Control = control,
                   Nice21 = Nice21,
                   p1star. = p1star. , p2star. = p2star. ,
                   n = n, M = M,
                   othint = othint, othdbl = othdbl,
                   alldump = FALSE)


      Cmat <- matrix(conjgrad$par, p2, Rank)  # old becoz of scale(cmatrix)


    if (converged <- (conjgrad$convergence == 0))
      break
  }

  if (!converged) {


    if (control$maxitl > 1) {
      warning("convergence not obtained in ", control$maxitl,
              " iterations.")
    } else {
      warning("convergence not obtained")
    }
  } else {
  }
  Cmat <- crow1C(Cmat, control$Crow1positive)  # Make sure signs are right

  flush.console()
  temp9 <- callcaoc(cmatrix = Cmat,
                    etamat = eta, xmat = x, ymat = y,
                    wvec = w, modelno = modelno,
                    Control = control,
                    Nice21 = Nice21,
                    p1star. = p1star. , p2star. = p2star. ,
                    n = n, M = M,
                    othint = othint, othdbl = othdbl,
                    alldump = TRUE)
  if (!is.list(extra))
    extra <- list()
  extra$Cmat <- temp9$Cmat

  ynames <- dimnames(y)[[2]]
  extra$df1.nl <- temp9$df1.nl
  extra$lambda1 <- temp9$lambda1
  extra$spar1 <- temp9$spar1
  names(extra$df1.nl) <-
  names(extra$lambda1) <-
  names(extra$spar1) <- ynames
  if (Rank == 2) {
    extra$spar2 <- temp9$spar2
    extra$lambda2 <- temp9$lambda2
    extra$df2.nl <- temp9$df2.nl
    names(extra$df2.nl) <-
    names(extra$lambda2) <-
    names(extra$spar2) <- ynames
  }

  extra$alldeviance <- temp9$alldeviance
  names(extra$alldeviance) <- ynames

  mu <- matrix(temp9$fitted, n, NOS, byrow = TRUE)











  dn <- labels(x)
  yn <- dn[[1]]


  if (is.matrix(mu)) {
    if (length(dimnames(y)[[2]])) {
      y.names <- dimnames(y)[[2]]
    }
    if (length(dimnames(mu)[[2]])) {
      y.names <- dimnames(mu)[[2]]
    }
    dimnames(mu) <- list(yn, y.names)
  } else {
    names(mu) <- names(fv)
  }


  fit <- list(
              fitted.values = mu,
              Cmatrix = Cmat,
              terms = Terms)  # terms: This used to be done in vglm()




  misc <- list(
      criterion = criterion,
      predictors.names = predictors.names,
      M = M,
      n = n,
      nonparametric = nonparametric,
      p = ncol(x),
      ynames = ynames)

  crit.list <- list()
  crit.list$deviance <- temp9$deviance





  if (w[1] != 1 || any(w != w[1]))
    fit$prior.weights <- w

  if (length(family@last))
    eval(family@last)

  structure(c(fit,
      temp9,
      list(
      contrasts = attr(x, "contrasts"),
      control = control,
      crit.list = crit.list,
      extra = extra,
      family = family,
      iter = iter,
      misc = misc,
      post = post,
      x = x,
      y = y)),
    vclass = family@vfamily)
}





cao.control <- function(Rank = 1,
          all.knots = FALSE,
          criterion = "deviance",
          Cinit = NULL,
          Crow1positive = TRUE,
          epsilon = 1.0e-05,
          Etamat.colmax = 10,
          GradientFunction = FALSE,  # For now 20041224
          iKvector = 0.1,
          iShape = 0.1,
          noRRR = ~ 1,
          Norrr = NA,
          SmallNo = 5.0e-13,
          Use.Init.Poisson.QO = TRUE,

          Bestof = if (length(Cinit)) 1 else 10,
          maxitl = 10,  # was 40 prior to 20100420
          imethod = 1,
          bf.epsilon = 1.0e-7,
          bf.maxit = 10,  # was 40 prior to 20100420
          Maxit.optim = 250,
          optim.maxit = 20,
          sd.sitescores = 1.0,
          sd.Cinit = 0.02,
          suppress.warnings = TRUE,
          trace = TRUE,
          df1.nl = 2.5,  # About 1.5--2.5 gives the flexibility of a quadratic
          df2.nl = 2.5,  # About 1.5--2.5 gives the flexibility of a quadratic
          spar1 = 0,  # 0 means df1.nl is used
          spar2 = 0,  # 0 means df2.nl is used
          ...) {


  if (length(Norrr) != 1 || !is.na(Norrr)) {
    warning("argument 'Norrr' has been replaced by 'noRRR'. ",
            "Assigning the latter but using 'Norrr' will become an error in ",
            "the next VGAM version soon.")
    noRRR <- Norrr
  }



  if (!is.Numeric(iShape, positive = TRUE))
    stop("bad input for argument 'iShape'")
  if (!is.Numeric(iKvector, positive = TRUE))
    stop("bad input for argument 'iKvector'")
  if (!is.Numeric(imethod, positive = TRUE, length.arg = 1,
                  integer.valued = TRUE))
    stop("bad input for argument 'imethod'")

  if (criterion != "deviance")
    stop("'criterion' must be 'deviance'")
  if (GradientFunction)
    stop("20050114; GradientFunction = TRUE not working yet")

  se.fit <- as.logical(FALSE)
  if (se.fit)
    stop("se.fit = FALSE handled only")

  if (length(Cinit) && !is.Numeric(Cinit))
    stop("Bad input for argument 'Cinit'")
  if (!is.Numeric(Bestof, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE))
    stop("Bad input for argument 'Bestof'")
  if (!is.Numeric(maxitl, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE))
    stop("Bad input for argument 'maxitl'")
  if (!is.Numeric(bf.epsilon, length.arg = 1,
                  positive = TRUE))
    stop("Bad input for argument 'bf.epsilon'")
  if (!is.Numeric(bf.maxit, integer.valued = TRUE,
                  positive = TRUE, length.arg = 1))
    stop("Bad input for argument 'bf.maxit'")

  if (!is.Numeric(Etamat.colmax, positive = TRUE,
                  length.arg = 1) ||
      Etamat.colmax < Rank)
    stop("bad input for argument 'Etamat.colmax'")

  if (!is.Numeric(Maxit.optim, integer.valued = TRUE,
                  positive = TRUE, length.arg = 1))
    stop("Bad input for argument 'Maxit.optim'")
  if (!is.Numeric(optim.maxit, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE))
    stop("Bad input for argument 'optim.maxit'")
  if (!is.Numeric(sd.sitescores, length.arg = 1,
                  positive = TRUE))
    stop("Bad input for argument 'sd.sitescores'")
  if (!is.Numeric(sd.Cinit, length.arg = 1,
                  positive = TRUE))
    stop("Bad input for argument 'sd.Cinit'")
  if (!is.Numeric(df1.nl) || any(df1.nl < 0))
    stop("Bad input for argument 'df1.nl'")
  if (any(df1.nl >= 0 & df1.nl < 0.05)) {
    warning("'df1.nl' values between 0 and 0.05 converted to 0.05")
    df1.nl[df1.nl < 0.05] <- 0.05
  }
  if (any(df1.nl > 3.5)) {
    warning("'df1.nl' values > 3.5 are excessive")
  }
  if (!is.Numeric(df2.nl) || any(df2.nl < 0))
    stop("Bad input for argument 'df2.nl'")
  if (any(df2.nl >= 0 & df2.nl < 0.05)) {
    warning("'df2.nl' values between 0 and 0.05 converted to 0.05")
    df2.nl[df2.nl < 0.05] <- 0.05
  }
  if (!is.Numeric(spar1) || any(spar1 < 0))
    stop("Bad input for argument 'spar1'")
  if (!is.Numeric(spar2) || any(spar2 < 0))
    stop("Bad input for argument 'spar2'")
  if (!is.Numeric(epsilon, positive = TRUE, length.arg = 1))
    stop("Bad input for argument 'epsilon'")

  if (!is.Numeric(SmallNo, positive = TRUE, length.arg = 1))
    stop("Bad input for argument 'SmallNo'")
  if ((SmallNo < .Machine$double.eps) ||
      (SmallNo > .0001))
    stop("'SmallNo' is out of range")

    ans <- list(
     Corner = FALSE,  # A constant, not a control parameter; unneeded?
     eq.tolerances = FALSE,  # A constant, not a control parameter; needed
     I.tolerances = FALSE,  # A constant, not a control parameter; unneeded?
     Quadratic = FALSE,  # A constant, not a control parameter; unneeded?
        all.knots = as.logical(all.knots)[1],
        Bestof = Bestof,
        Cinit = Cinit,
        ConstrainedO = TRUE,  # A constant, not a control parameter
        criterion = criterion,
        Crow1positive = as.logical(rep_len(Crow1positive, Rank)),
        epsilon = epsilon,
        Etamat.colmax = Etamat.colmax,
        FastAlgorithm = TRUE,  # A constant, not a control parameter
        GradientFunction = as.logical(GradientFunction),
        maxitl = maxitl,
        bf.epsilon = bf.epsilon,
        bf.maxit = bf.maxit,
        imethod = imethod,
        Maxit.optim = Maxit.optim,
        optim.maxit = optim.maxit,
        noRRR = noRRR,
        Rank = Rank,
        sd.sitescores = sd.sitescores,
        sd.Cinit = sd.Cinit,
        se.fit = se.fit,  # If TRUE, then would need storage for S QR fits
        SmallNo = SmallNo,
        suppress.warnings = as.logical(suppress.warnings),
        trace = as.integer(trace),
        Use.Init.Poisson.QO = Use.Init.Poisson.QO,
        iKvector = as.numeric(iKvector),
        iShape = as.numeric(iShape),
        DF1 = 2.5,    # Used as Default value if df1.nl has no default
        DF2 = 2.5,    # Used as Default value if df2.nl has no default
        SPAR1 = 0,    # Used as Default value if spar1 has no default
        SPAR2 = 0,    # Used as Default value if spar2 has no default
        df1.nl = df1.nl,
        df2.nl = df2.nl,
        spar1 = spar1,
        spar2 = spar2)
    ans
}


create.cms <- function(Rank = 1, M, MSratio = 1, which, p1 = 1) {
  if (!is.Numeric(p1, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'p1'")
  Hlist. <- vector("list", p1 + Rank)
  for (rr in 1:(p1+Rank))
    Hlist.[[rr]] <- diag(M)
  names(Hlist.) <- if (p1 == 1) c("(Intercept)", names(which)) else stop()
  if (MSratio == 2) {
    for (r in 1:Rank)
      Hlist.[[p1+r]] <- eijfun(1, M)
  }
  Hlist.
}




callcaoc <- function(cmatrix,
                    etamat, xmat, ymat, wvec, modelno,
                    Control, Nice21 = TRUE,
                    p1star. = if (modelno %in% c(3, 5)) 2 else 1,
                    p2star. = Rank,
                    n, M,
                    othint, othdbl,
                    alldump = FALSE) {
  flush.console()

  control <- Control
  Rank <- control$Rank
  p1 <- length(control$colx1.index)
  p2 <- length(control$colx2.index)
  yn <- dimnames(ymat)[[2]]
  if (length(yn) != ncol(ymat))
    stop("the column names of argument 'ymat' must be given")

  queue <- qbig <- Rank  # 20051019; number of smooths per species

  NOS <- if (modelno %in% c(3, 5)) M/2 else M
  df1.nl <- procVec(control$df1.nl, yn = yn , Default = control$DF1)
  spar1  <- procVec(control$spar1,  yn = yn , Default = control$SPAR1)
  df2.nl <- procVec(control$df2.nl, yn = yn , Default = control$DF2)
  spar2  <- procVec(control$spar2,  yn = yn , Default = control$SPAR2)
  if (any(c(length(spar1), length(spar2), length(df1.nl),
            length(df2.nl)) != NOS))
    stop("wrong length in at least one of arguments ",
         "'df1.nl', 'df2.nl', 'spar1', 'spar2'")

  cmatrix <- matrix(cmatrix, p2, Rank)  # crow1C() needs a matrix as input
  cmatrix <- crow1C(cmatrix, crow1positive = control$Crow1positive)
  numat <- xmat[, control$colx2.index, drop = FALSE] %*% cmatrix
  evnu <- eigen(var(numat), symmetric = TRUE)
  temp7 <- if (Rank > 1) evnu$vector %*% diag(evnu$value^(-0.5)) else
           evnu$vector %*% evnu$value^(-0.5)
  cmatrix <- cmatrix %*% temp7
  cmatrix <- crow1C(cmatrix, crow1positive = control$Crow1positive)
  numat <- xmat[, control$colx2.index, drop = FALSE] %*% cmatrix


  dim(numat) <- c(n, Rank)
  mynames5 <- if (Rank == 1) "latvar" else paste("latvar", 1:Rank, sep = "")
  nu1mat <- cbind("(Intercept)" = 1, latvar = numat)
  dimnames(nu1mat) <- list(dimnames(xmat)[[1]],
                           c("(Intercept)", mynames5))

  temp.smooth.frame <- vector("list", p1+Rank)  # Temporary makeshift frame
  names(temp.smooth.frame) <- c(names(control$colx1.index), mynames5)
  for (uu in 1:(p1+Rank)) {
    temp.smooth.frame[[uu]] <- nu1mat[, uu]
  }
  temp.smooth.frame <- data.frame(temp.smooth.frame)
  for (uu in 1:Rank) {
    attr(temp.smooth.frame[,uu+p1], "spar") <- 0  # this value unused
    attr(temp.smooth.frame[,uu+p1], "df") <- 4    # this value unused
  }

  pstar. <- p1star. + p2star.  # = Mdot + Rank
  nstar <- if (Nice21) ifelse(modelno %in% c(3, 5), n * 2, n) else n * M
  lenbeta <- pstar. * ifelse(Nice21, NOS, 1)  # Holds the linear coeffs

  inited <- if (exists(".VGAM.CAO.etamat", envir = VGAMenv)) 1 else 0
  usethiseta <-
      if (inited == 1)
        getfromVGAMenv("etamat", prefix = ".VGAM.CAO.") else
        t(etamat)

  if (anyNA(usethiseta)) {
    usethiseta <- t(etamat)  # So that dim(usethiseta) == c(M,n)
    rmfromVGAMenv("etamat", prefix = ".VGAM.CAO.")
  }

  usethisbeta <- if (inited == 2)
                 getfromVGAMenv("beta", prefix = ".VGAM.CAO.") else
                 double(lenbeta)
  othint[5] <- inited   # Refine initialization within C
  pstar <- NOS * pstar.
  bnumat <- if (Nice21) matrix(0, nstar, pstar.) else
                        stop("code not written here")

  M. <- MSratio <- M / NOS  # 1 or 2 usually
  which <- p1 + (1:Rank)  # These columns are smoothed
  nwhich <- names(which) <- mynames5

  origHlist <-
  Hlist. <- create.cms(Rank = Rank, M = M., MSratio = MSratio,
                       which = which, p1 = p1)  # For 1 species only
  ncolHlist. <- unlist(lapply(Hlist. , ncol))
  smooth.frame <- s.vam(x = nu1mat, zedd = NULL,
                        wz = NULL, smomat = NULL,
                        which = which,
                        smooth.frame = temp.smooth.frame,
                        bf.maxit = control$bf.maxit,
                        bf.epsilon = control$bf.epsilon,
                        trace = FALSE, se.fit = control$se.fit,
                        X.vlm.save = bnumat, Hlist = Hlist. ,
                        ncolHlist = ncolHlist. ,
                        M =  M. ,
                        qbig = NULL, Umat = NULL,  # NULL ==> unneeded
                        all.knots = control$all.knots, nk = NULL,
                        sf.only = TRUE)

  ldk <- 3 * max(ncolHlist.[nwhich]) + 1   # 20020711

  dimw. <- M.   # Smoothing one spp. at a time
  dim1U. <- M.
  wz. <- matrix(0, n, dimw. )
  if (names(Hlist.)[1] != "(Intercept)")
    stop("something wrong here")
  Hlist.[[1]] <- NULL

  trivc <- rep_len(2 - M. , queue)
  ncbvec <- ncolHlist.[nwhich]
  ncolb <- max(ncbvec)

  qbig. <- NOS * qbig    # == NOS * Rank; holds all the smooths
  if (!all(as.vector(ncbvec) == rep_len(1, queue)))
    stop("'ncbvec' not right---should be a queue-vector of ones")
  pbig <- pstar. #



  contr.sp <- list(low = -1.5,  ## low = 0.      was default till R 1.3.x
                   high = 1.5,
                   tol = 1e-4,  ## tol = 0.001   was default till R 1.3.x
                   eps = 2e-8,  ## eps = 0.00244 was default till R 1.3.x
                   maxit = 500 )


  npetc <-
    c(n = nrow(nu1mat), p. = ncol(nu1mat), q = length(which),
        se.fit = control$se.fit, 0,
      control$bf.maxit, qrank = 0, M = M. , nbig = nstar, pbig = pbig,
      qbig = qbig, dim2wz = dimw. , dim1U = dim1U. ,
        ierror = 0, ldk = ldk,
      contr.sp$maxit, iinfo = 0)




  if (Rank == 2) {
    smopar <- (c(spar1, spar2))[interleave.VGAM(4 * NOS, M1 = 2)]
    dofvec <- (1.0 + c(df1.nl, df2.nl))[interleave.VGAM(4 * NOS, M1 = 2)]
    lamvec <- 0 * dofvec
    stop("20100414; havent got Rank = 2 going yet")
  } else {
    smopar <- c(spar1, spar2)
    dofvec <- c(df1.nl, df2.nl) + 1.0
    lamvec <- 0 * dofvec
  }

  ans1 <- .C("vcao6",
  numat = as.double(numat), ymat = as.double(ymat), wvec = as.double(wvec),
  etamat = as.double(usethiseta), fv = double(NOS*n), zedd = double(n*M),
  wz = double(n*M), U = double(M*n),  # bnumat = as.double(bnumat),
  qr = double(nstar*pstar.), qraux = double(pstar.),
    qpivot = integer(pstar.),
  n = as.integer(n), M = as.integer(M), NOS = as.integer(NOS),
      nstar = as.integer(nstar), dim1U = as.integer( M ),  # for U, not U.
  errcode = integer(1), othint = as.integer(othint),
  deviance = double(1 + NOS),  # NOS more elts added 20100413
  beta = as.double(usethisbeta),
  othdbl = as.double(othdbl),
      npetc = as.integer(npetc), M. = as.integer( M. ),
  dofvec = as.double(dofvec),
  lamvec = as.double(lamvec),
  smopar = as.double(smopar),
      match = as.integer(smooth.frame$matcho),
      as.integer(smooth.frame$nef),
      which = as.integer(which),
      smomat = as.double(matrix(0, n, qbig. )),
      nu1mat = as.double(nu1mat),
  Hlist = as.double(unlist( Hlist. )),
  as.integer(ncbvec),
      smap = as.integer(1:(Rank+1)),  #
      trivc = as.integer(trivc),




  levmat = double(NOS * sum(smooth.frame$neffec * ncbvec)),



      bcoefficients = double(NOS * sum(smooth.frame$nknots*ncbvec)),
      xknots = as.double(unlist(smooth.frame$knots)),
  bindex = as.integer(smooth.frame$bindex),
  lindex = as.integer(smooth.frame$lindex),
      nknots = as.integer(smooth.frame$nknots),
      kindex = as.integer(smooth.frame$kindex))
flush.console()


    if (ans1$errcode == 0) {
      assign2VGAMenv(c("etamat", "beta"), ans1, prefix = ".VGAM.CAO.")
      assign(".VGAM.CAO.cmatrix", matrix(cmatrix, p2, Rank),
             envir = VGAMenv)
    } else {
      if (!control$suppress.warnings) {
        cat("warning in callcaoc: error code  = ", ans1$errcode, "\n")
        cat("warning in callcaoc: npetc[14]   = ", ans1$npetc[14], "\n")
        flush.console()
      }
      rmfromVGAMenv(c("etamat", "beta"), prefix = ".VGAM.CAO.")
    }

  returnans <- if (alldump) {
      bindex <- ans1$bindex
      ncolHlist <- ncbvec
      Bspline2 <- vector("list", NOS)
      names(Bspline2) <- dimnames(ymat)[[2]]
      Bspline <- vector("list", length(nwhich))
      names(Bspline) <- nwhich
      ind9 <- 0   # moving index
      for (sppno in 1:NOS) {
        for (ii in seq_along(nwhich)) {
          ind7 <- (smooth.frame$bindex[ii]):(smooth.frame$bindex[ii+1]-1)
          ans <- ans1$bcoeff[ind9+ind7]
          ans <- matrix(ans, ncol = ncolHlist[nwhich[ii]])
          Bspline[[ii]] <-
            new(Class = "vsmooth.spline.fit",
                "Bcoefficients" = ans,
                "xmax"          = smooth.frame$xmax[ii],
                "xmin"          = smooth.frame$xmin[ii],
                "knots"         = as.vector(smooth.frame$knots[[ii]]))
        }
        ind9 <- ind9 + smooth.frame$bindex[length(nwhich)+1] - 1
        Bspline2[[sppno]] <- Bspline
      }

      qrank <- npetc[7]  # Assume all species have the same qrank value
      dim(ans1$etamat) <- c(M, n)  # was c(n, M) prior to 20060822



      df1.nl  <- ans1$dofvec[1:NOS] - 1.0
      lambda1 <- ans1$lamvec[1:NOS]
      spar1   <- ans1$smopar[1:NOS]
      if (Rank == 2) {
 stop("20100414; this is not working yet")
        df2.nl  <- ans1$dofvec[NOS + (1:NOS)] - 1.0
        lambda2 <- ans1$lamvec[NOS + (1:NOS)]
        spar2   <- ans1$smopar[NOS + (1:NOS)]
      }

      list(deviance = ans1$deviance[1],
           alldeviance = ans1$deviance[-1],
           bcoefficients = ans1$bcoefficients,
           bindex = ans1$bindex,
           Bspline = Bspline2,
           Cmat = matrix(cmatrix, p2, Rank, dimnames = list(
                         names(control$colx2.index), mynames5)),
           coefficients = ans1$beta,
           df1.nl = df1.nl,
           df2.nl = if (Rank == 2) df2.nl else NULL,
           df.residual = n*M - qrank - sum(ans1$df - 1),
           fitted = ans1$fv,  # NOS x n
           kindex = ans1$kindex,
           lambda1 = lambda1,
           lambda2 = if (Rank == 2) lambda2 else NULL,
           predictors = matrix(ans1$etamat, n, M, byrow = TRUE),
           wresiduals = ans1$zedd - t(ans1$etamat),  # n x M
           spar1 = spar1,
           spar2 = if (Rank == 2) spar2 else NULL)
    } else {
      ans1$deviance[1]
    }
  flush.console()
  returnans
}



calldcaoc <- function(cmatrix,
                     etamat, xmat, ymat, wvec, modelno,
                     Control, Nice21 = TRUE,
                     p1star. = if (modelno %in% c(3, 5)) 2 else 1,
                     p2star. = Rank,
                     n, M,
                     othint, othdbl,
                     alldump = FALSE) {


  if (alldump)
    stop("really used?")
  flush.console()



  U <- NULL





  if (!Nice21)
    stop("'Nice21' must be TRUE")
  control <- Control
  Rank <- control$Rank
  p2 <- length(control$colx2.index)
  yn <- dimnames(ymat)[[2]]
  if (!length( yn ))
    yn <- paste("Y", 1:ncol(ymat), sep = "")


  cmatrix <- scale(cmatrix)

  xmat2 <- xmat[, control$colx2.index, drop = FALSE]   #ccc
  numat <- xmat2 %*% matrix(cmatrix, p2, Rank)
  dim(numat) <- c(nrow(xmat), Rank)
  temp.smooth.frame <- vector("list", 1+Rank)  # Temporary makeshift frame
  mynames5 <- if (Rank == 1) "latvar" else paste("latvar", 1:Rank, sep = "")
  names(temp.smooth.frame) <- c("(Intercept)", mynames5)
  temp.smooth.frame[[1]] <- rep_len(1, n)
  for (uu in 1:Rank) {
    temp.smooth.frame[[uu+1]] <- numat[, uu]
  }
  temp.smooth.frame <- data.frame(temp.smooth.frame)
  for (uu in 1:Rank) {
    attr(temp.smooth.frame[,uu+1], "spar") <- 0  # any old value
    attr(temp.smooth.frame[,uu+1], "df") <- 4    # any old value
  }
  pstar.  <- p1star.  + p2star.
  nstar <- if (Nice21) ifelse(modelno %in% c(3, 5), n * 2, n) else n * M
  NOS <- ifelse(modelno %in% c(3, 5), M / 2, M)
  lenbeta <- pstar. * ifelse(Nice21, NOS, 1)

  if (TRUE) {
    inited <- if (exists(".VGAM.CAO.etamat", envir = VGAMenv))
              1 else 0
    usethiseta <- if (inited == 1)
                  get(".VGAM.CAO.etamat", envir = VGAMenv) else
                  t(etamat)
  }
  usethisbeta <- if (inited == 2)
                 get(".VGAM.CAO.beta", envir = VGAMenv) else
                 double(lenbeta)





  pstar <- NOS * pstar.
  bnumat <- if (Nice21)
            matrix(0, nstar, pstar) else stop("need 'Nice21'")

  M. <- MSratio <- M / NOS  # 1 or 2 usually


  p1 <- 1

  which <- p1 + (1:Rank)  # The first 1 is the intercept term
  nwhich <- names(which) <- mynames5

  origHlist <- Hlist. <-
    create.cms(Rank = Rank, M = M., MSratio = MSratio,
               which = which, p1 = p1)  # For 1 species
  ncolHlist. <- unlist(lapply(Hlist. , ncol))
    nu1mat <- cbind("(Intercept)" = 1, latvar = numat)
    dimnames(nu1mat) <- list(dimnames(xmat)[[1]],
                             c("(Intercept)", "latvar"))

    smooth.frame <- s.vam(x = nu1mat, zedd = NULL, wz = NULL,
                          smomat = NULL, which = which,
                          smooth.frame = temp.smooth.frame,
                          bf.maxit = control$bf.maxit,
                          bf.epsilon = control$bf.epsilon,
                          trace = FALSE, se.fit = control$se.fit,
                          X.vlm.save = bnumat, Hlist = Hlist.,
                          ncolHlist = ncolHlist. ,
                          M = M. , qbig = NULL,

                          Umat = U,  # NULL value ==> not needed
                          all.knots = control$all.knots, nk = NULL,
                          sf.only = TRUE)

    ldk <- 4 * max(ncolHlist.[nwhich])   # was M;     # Prior to 20020711
    ldk <- 3 * max(ncolHlist.[nwhich]) + 1   # 20020711



    wz. <- matrix(0, n, M. )  # not sure
    dimw. <- if (is.matrix( wz. )) ncol( wz. ) else 1


    dim1U. <- M.  # 20100410




    queue <- qbig <- Rank  # 20051019; number of smooths per species



    Hlist.[[1]] <- NULL
    trivc <- rep_len(2 - M. , queue)
    ncbvec <- ncolHlist.[nwhich]
    ncolb <- max(ncbvec)


    qbig. <- NOS * qbig    # == NOS * Rank
    pbig <- pstar. # Not sure
    if (FALSE) {
      df1.nl <- rep_len(control$df1.nl, NOS)  # This is used
      df2.nl <- rep_len(control$df2.nl, NOS)  # This is used
      spar1  <- rep_len(control$spar1,  NOS)  # This is used
      spar2  <- rep_len(control$spar2,  NOS)  # This is used
    } else {
      df1.nl <- procVec(control$df1.nl, yn = yn , Default = control$DF1)
      df2.nl <- df1.nl  # 20100417; stopgap
      spar1  <- procVec(control$spar1,  yn = yn , Default = control$SPAR1)
      spar2  <- spar1  # 20100417; stopgap
      dofvec <- c(df1.nl, df2.nl)
      lamvec <- 0 * dofvec
      smopar <- c(spar1, spar2)
    }





    contr.sp <- list(low = -1.5,  ## low = 0.      was default till R 1.3.x
                     high = 1.5,
                     tol = 1e-4,  ## tol = 0.001   was default till R 1.3.x
                     eps = 2e-8,  ## eps = 0.00244 was default till R 1.3.x
                     maxit = 500 )



warning("20100405; this is old:")
    npetc <-
      c(n = n, p = 1+Rank, length(which), se.fit = control$se.fit, 0,
        maxitl = control$maxitl, qrank = 0, M =  M. , n.M = n* M. ,
          pbig = sum( ncolHlist.),
        qbig = qbig, dimw =  dimw. , dim1U =  dim1U. ,
          ierror = 0, ldk = ldk)

warning("20100405; this is new:")
    npetc <- c(n = nrow(nu1mat), p.  = ncol(nu1mat),
               q = length(which),
               se.fit = control$se.fit, 0,
    control$bf.maxit, qrank = 0, M =  M. , nbig = nstar, pbig = pbig,
    qbig = qbig, dim2wz =  dimw. , dim1U =  dim1U. , ierror = 0, ldk = ldk,
    contr.sp$maxit, iinfo = 0)

    flush.console()

    if (!Nice21)
      stop("need 'Nice21'")

    ans1 <- .C("vdcao6",
    numat = as.double(numat), as.double(ymat), as.double(wvec),
    etamat = as.double(usethiseta), fv = double(NOS*n),
      zedd = double(n*M),
    wz = double(n*M), U = double(M*n),  # bnumat = as.double(bnumat),
    qr = double(nstar*pstar.), qraux = double(pstar.),
      qpivot = integer(pstar.),
    as.integer(n), as.integer(M), NOS = as.integer(NOS),
        as.integer(nstar), dim1U = as.integer(M),
    errcode = integer(1), othint = as.integer(othint),
    deviance  =  double(1 + NOS), beta = as.double(usethisbeta),
    othdbl = as.double(othdbl),
    as.double(xmat2),
    cmat = as.double(cmatrix),
    p2 = as.integer(p2), deriv = double(p2 * Rank),
    betasave = double(lenbeta),
    npetc = as.integer(npetc), M. = as.integer( M. ),
    dofvec = as.double(dofvec + 1.0),
    lamvec = as.double(0 * dofvec),
    smopar = as.double(smopar),
    match = as.integer(smooth.frame$matcho),
    as.integer(smooth.frame$nef),
    as.integer(which),
    smomat = as.double(matrix(0, n, qbig. )),
        nu1mat = as.double(nu1mat),
    as.double(unlist( Hlist. )),
    as.integer(ncbvec), smap = as.integer(1:(Rank+1)),
    trivc = as.integer(trivc),




  levmat = double(NOS * sum(smooth.frame$neffec * ncbvec)),



    bcoefficients = double(NOS * sum(smooth.frame$nknots * ncbvec)),
    xknots = as.double(unlist(smooth.frame$knots)),
    bindex = as.integer(smooth.frame$bindex),
    lindex = as.integer(smooth.frame$lindex),
    nknots = as.integer(smooth.frame$nknots),
    kindex = as.integer(smooth.frame$kindex))
        flush.console()

         assign(".VGAM.CAO.etamat", ans1$etamat, envir = VGAMenv)
         assign(".VGAM.CAO.z", ans1$zedd, envir = VGAMenv)
         assign(".VGAM.CAO.U", ans1$U, envir = VGAMenv)  # U
       if (ans1$errcode == 0) {
       } else {
         cat("warning in calldcaoc: error code  = ", ans1$errcode, "\n")
         flush.console()
       }

  returnans <- if (alldump) {
    bindex <- ans1$bindex
    ncolHlist <- ncbvec
    Bspline2 <- vector("list", NOS)
    names(Bspline2) <- dimnames(ymat)[[2]]
    Bspline <- vector("list", length(nwhich))
    names(Bspline) <- nwhich
    ind9 <- 0   # moving index
    for (jay in 1:NOS) {
      for (ii in seq_along(nwhich)) {
        ind9 <- ind9[length(ind9)] + (bindex[ii]):(bindex[ii+1]-1)
        ans <- ans1$bcoeff[ind9]
        ans <- matrix(ans, ncol = ncolHlist[nwhich[ii]])
        Bspline[[ii]] <-
          new(Class = "vsmooth.spline.fit",
              "Bcoefficients" = ans,
              "xmax"          = smooth.frame$xmax[ii],
              "xmin"          = smooth.frame$xmin[ii],
              "knots"         = as.vector(smooth.frame$knots[[ii]]))
      }
      Bspline2[[jay]] <- Bspline
    }

    qrank <- npetc[7]  # Assume all species have the same qrank value
    dim(ans1$etamat) <- c(M,n)   # bug: was c(n,M) prior to 20060822
    list(deviance    = ans1$deviance[1],
         alldeviance = ans1$deviance[-1],
         bcoefficients = ans1$bcoefficients,
         bindex = ans1$bindex,
         Bspline = Bspline2,
         Cmat = matrix(cmatrix, p2, Rank, dimnames = list(
                     names(control$colx2.index), mynames5)),
         coefficients = ans1$beta,
         df1.nl = ans1$dofvec[1:NOS] - 1,
         df2.nl = if (Rank == 2) ans1$dofvec[2 * (1:NOS) - 1] - 1 else NULL,
         lambda1 = ans1$lambda[1:NOS],
         lambda2 = if (Rank == 2) ans1$lambda[2 * (1:NOS) - 1] else NULL,
         df.residual = n * M - qrank - sum(ans1$df - 1),
         fitted = ans1$fv,
         kindex = ans1$kindex,
         predictors=matrix(ans1$etamat, n, M, byrow = TRUE),
         wresiduals = ans1$zedd - t(ans1$etamat),  # n x M
         spar1 = ans1$smopar[1:NOS],
         spar2 = if (Rank == 2) ans1$smopar[2 * (1:NOS) - 1] else NULL)
  } else {
    ans1$deriv
  }
  flush.console()
  returnans
}






setClass(Class = "Coef.rrvgam", representation(
      "Bspline"      = "list",
      "C"            = "matrix",
      "Constrained"  = "logical",
      "df1.nl"       = "numeric",
      "df2.nl"       = "numeric",
      "dispersion"   = "numeric",
      "eta2"         = "matrix",
      "latvar"       = "matrix",
      "latvar.order" = "matrix",
      "M"            = "numeric",
      "Maximum"      = "numeric",
      "NOS"          = "numeric",
      "Optimum"      = "matrix",
      "Optimum.order"= "matrix",
      "Rank"         = "numeric",
      "spar1"        = "numeric",
      "spar2"        = "numeric"))







Coef.rrvgam <- function(object,
    epsOptimum = 0.00001,  # Determines how accurately Optimum is estimated
    gridlen = 40,      # Number of points on the grid (one level at a time)
    maxgriditer = 10,  # Maximum number of iters allowed for grid search
    smallno = 0.05, ...) {

  if (!is.Numeric(epsOptimum, positive = TRUE, length.arg = 1))
    stop("bad input for argument 'epsOptimum'")
  if (!is.Numeric(gridlen, positive = TRUE, integer.valued = TRUE) ||
      gridlen < 5)
    stop("bad input for argument 'gridlen'")
  if (!is.Numeric(maxgriditer, positive = TRUE,
                  length.arg = 1, integer.valued = TRUE) ||
      maxgriditer < 3)
    stop("bad input for argument 'maxgriditer'")
  if (!is.logical(ConstrainedO <- object@control$ConstrainedO))
    stop("cannot determine whether the model is constrained or not")
  if (!is.Numeric(smallno, positive = TRUE, length.arg = 1) ||
     smallno > 0.5 || smallno < 0.0001)
    stop("bad input for argument 'smallno'")


  ocontrol <- object@control
  if ((Rank <- ocontrol$Rank) > 2) stop("'Rank' must be 1 or 2")
  gridlen <- rep_len(gridlen, Rank)
  M <- if (any(slotNames(object) == "predictors") &&
           is.matrix(object@predictors))
       ncol(object@predictors) else
       object@misc$M
  NOS <- if (length(object@y)) ncol(object@y) else M
    MSratio <- M / NOS  # 1 or 2; First value is g(mean)=quadratic form in latvar
    nice21 <- (length(ocontrol$colx1.index) == 1) &&
              (names(ocontrol$colx1.index) == "(Intercept)")
    if (!nice21)
      stop("Can only handle 'noRRR = ~ 1'")

    p1 <- length(ocontrol$colx1.index)
    p2 <- length(ocontrol$colx2.index)
    modelno <- object@control$modelno  # 1,2,3,... or 0
    ynames <- object@misc$ynames
    if (!length(ynames)) ynames <- object@misc$predictors.names
    if (!length(ynames)) ynames <- object@misc$ynames
    if (!length(ynames)) ynames <- paste("Y", 1:NOS, sep = "")
    lp.names <- object@misc$predictors.names
    if (!length(lp.names)) lp.names <- NULL

    latvar.names <-
      if (Rank == 1) "latvar" else paste("latvar", 1:Rank, sep = "")
    Cmat <- object@extra$Cmat  # p2 x Rank (provided maxitl > 1)
    if (ConstrainedO)
      dimnames(Cmat) <- list(names(ocontrol$colx2.index), latvar.names)
    latvar.mat <- if (ConstrainedO) {
      object@x[, ocontrol$colx2.index, drop = FALSE] %*% Cmat
    } else {
      object@latvar
    }

    optimum <- matrix(NA_real_, Rank, NOS,
                      dimnames = list(latvar.names, ynames))
    extents <- apply(latvar.mat, 2, range)  # 2 by R

    maximum <- rep_len(NA_real_, NOS)

    which.species <- 1:NOS  # Do it for all species
    if (Rank == 1) {
      gridd <- cbind(seq(extents[1, 1], extents[2, 1], len = gridlen))
      eta2matrix <- matrix(0, NOS, 1)  # Added 20160716
    } else {
      gridd <-
        expand.grid(seq(extents[1, 1], extents[2, 1], len = gridlen[1]),
                    seq(extents[1, 2], extents[2, 2], len = gridlen[2]))
      eta2matrix <- matrix(0, NOS, 1)
    }
    gridd.orig <- gridd
    for (sppno in seq_along(which.species)) {
      gridd <- gridd.orig
      gridres1 <- gridd[2, 1] - gridd[1, 1]
      gridres2 <- if (Rank == 2) gridd[2, 2] - gridd[1, 2] else 0
      griditer <- 1

      thisSpecies <- which.species[sppno]
      indexSpecies <- if (is.character(which.species))
          match(which.species[sppno], ynames) else which.species[sppno]

      if (is.na(indexSpecies))
        stop("mismatch found in 'which.species'")

      while (griditer == 1 ||
             ((griditer <= maxgriditer) &&
             ((gridres1 > epsOptimum) ||
              (gridres2 > epsOptimum)))) {
        temp <- predictrrvgam(object, grid = gridd, sppno = thisSpecies,
                           Rank = Rank, deriv = 0, MSratio = MSratio)
        yvals <- temp$yvals  # gridlen-vector
        xvals <- temp$xvals  # gridlen x Rank; gridd
        if (length(temp$eta2))
          eta2matrix[sppno, 1] <- temp$eta2

        nnn <- length(yvals)
        index <- (1:nnn)[yvals == max(yvals)]
        if (length(index) != 1)
          warning("could not find a single maximum")
        if (Rank == 2) {
          initvalue <- rep_len(xvals[index,], Rank)  # for optim()
          if (abs(initvalue[1] - extents[1, 1]) < smallno)
            initvalue[1] <- extents[1, 1] + smallno
          if (abs(initvalue[1] - extents[2, 1]) < smallno)
            initvalue[1] <- extents[2, 1] - smallno
          if (abs(initvalue[2] - extents[1, 2]) < smallno)
            initvalue[2] <- extents[1, 2] + smallno
          if (abs(initvalue[2] - extents[2, 2]) < smallno)
            initvalue[2] <- extents[2, 2] - smallno
          break
        }
        if (index == 1 || index == nnn) {
          maximum[sppno] <- optimum[1, sppno] <- NA
          gridres1 <- epsOptimum + 1  # equivalent to a break
          break  # just in case
        } else {
          maximum[sppno] <- yvals[index]  # On the eta scale
          optimum[1, sppno] <- xvals[index, 1]
          gridd[, 1] <- seq(
                  max(extents[1, 1], optimum[1, sppno] - gridres1),
                  min(extents[2, 1], optimum[1, sppno] + gridres1),
                  len = gridlen)
          gridres1 <- gridd[2, 1] - gridd[1, 1]
          griditer <- griditer + 1
        }
      }  # of while

      if (Rank == 2) {
        myfun <- function(x, object, sppno, Rank = 1,
                          deriv = 0, MSratio = 1) {
          x <- matrix(x, 1, length(x))
          temp <- predictrrvgam(object, grid = x, sppno = sppno,
                             Rank = Rank, deriv = deriv, MSratio = MSratio)
          temp$yval
        }
        answer <- optim(initvalue, myfun, gr = NULL, method = "L-BFGS-B",
                        lower = extents[1, ], upper = extents[2, ],
                        control = list(fnscale = -1),  # maximize!
                        object = object, sppno = sppno, Rank = Rank,
                        deriv = 0, MSratio = MSratio)
        for (rindex in 1:Rank)
          if (abs(answer$par[rindex] - extents[1, rindex]) > smallno &&
              abs(answer$par[rindex] - extents[2, rindex]) > smallno) {
            optimum[rindex,sppno] <- answer$par[rindex]
             maximum[sppno] <- answer$value
          }
        }  # end of Rank = 2
    }  # end of sppno
    myetamat <- rbind(maximum)
    if (MSratio == 2)
      myetamat <- kronecker(myetamat, matrix(1:0, 1, 2))
    maximum <- object@family@linkinv(eta = myetamat, extra = object@extra)
    maximum <- c(maximum)  # Convert from matrix to vector
    names(maximum) <- ynames

    ans <- new(Class = "Coef.rrvgam",
               Bspline = object@Bspline,
               Constrained = ConstrainedO,
               df1.nl = object@extra$df1.nl,
               latvar = latvar.mat,
               latvar.order = latvar.mat,
               Maximum = maximum,
               M = M,
               NOS = NOS,
               Optimum = optimum,
               Optimum.order = optimum,
               Rank = Rank,
               spar1 = object@extra$spar1)
    if (ConstrainedO) {
      ans@C <- Cmat
    } else {
      Cmat <- NULL
    }
    if (Rank == 2) {
      dimnames(eta2matrix) <-
        list(object@misc$predictors.names[c(FALSE, TRUE)], " ")
      ans@eta2 <- eta2matrix
      ans@df2.nl <- object@extra$df2.nl
      ans@spar2  <- object@extra$spar2
    }

    for (rindex in 1:Rank) {
      ans@Optimum.order[rindex, ] <- order(ans@Optimum[rindex, ])
      ans@latvar.order[, rindex]  <- order(ans@latvar[, rindex])
    }

  if (length(object@misc$estimated.dispersion) &&
      object@misc$estimated.dispersion) {
    p <- length(object@coefficients)
    n <- object@misc$n
    M <- object@misc$M
    NOS <- if (length(object@y)) ncol(object@y) else M
    pstar <- p + length(Cmat)  # Adjustment
    adjusted.dispersion <- object@misc$dispersion *
                           (n * M - p) / (n * M - pstar)
    ans@dispersion <- adjusted.dispersion
  }
  if (MSratio == 2) {
    lcoef <- object@coefficients
    temp <- lcoef[((1:NOS)-1) * (2+Rank)+2]
    names(temp) <- object@misc$predictors.names[2 * (1:NOS)]
    ans@dispersion <- temp
  }
  dimnames(ans@Optimum) <- list(latvar.names, ynames)
  ans
}



show.Coef.rrvgam <- function(object,
                          digits = max(2, options()$digits-2), ...) {
  Rank <- object@Rank
  NOS <- object@NOS
  M <- object@M

  Maximum <- if (length(object@Maximum))
             cbind(Maximum = object@Maximum) else NULL
  optmat <- cbind(t(object@Optimum))
  dimnames(optmat) <- list(dimnames(optmat)[[1]],
                           if (Rank > 1)
                           paste("Optimum",
                                 dimnames(optmat)[[2]], sep = ".") else
                                 "Optimum")

  if ( object@Constrained ) {
    cat("\nC matrix (constrained/canonical coefficients)\n")
    print(object@C, digits = digits, ...)
  }
  cat("\nOptimums and maximums\n")
  print(cbind(Optimum = optmat,
              Maximum), digits = max(1, digits-1))
  cat("\nNonlinear degrees of freedom\n")
  if (Rank == 1) {
    print(cbind(df1.nl = object@df1.nl), digits = max(2, digits-1), ...)
  } else {
    print(cbind(df1.nl = object@df1.nl,
                df2.nl = object@df2.nl), digits = max(2, digits-1), ...)
  }
  invisible(object)
}





setMethod("show", "Coef.rrvgam", function(object)
  show.Coef.rrvgam(object))





setMethod("coef", "rrvgam", function(object, ...) Coef.rrvgam(object, ...))
setMethod("coefficients", "rrvgam", function(object, ...)
    Coef.rrvgam(object, ...))
setMethod("Coef", "rrvgam", function(object, ...) Coef.rrvgam(object, ...))




lvplot.rrvgam <- function(object,
          add = FALSE, show.plot = TRUE, rugplot = TRUE, y = FALSE,
          type = c("fitted.values", "predictors"),
          xlab = paste("Latent Variable",
                       if (Rank == 1) "" else " 1", sep = ""),
          ylab = if (Rank == 1) switch(type, predictors = "Predictors",
              fitted.values = "Fitted values") else "Latent Variable 2",
          pcex = par()$cex, pcol = par()$col, pch = par()$pch,
          llty = par()$lty, lcol = par()$col, llwd = par()$lwd,
          label.arg= FALSE, adj.arg=-0.5,
          sites= FALSE, spch = NULL, scol = par()$col, scex = par()$cex,
          sfont = par()$font,
          which.species = NULL,
          check.ok = TRUE, ...) {
    type <- match.arg(type, c("fitted.values", "predictors"))[1]

    if ((Rank <- object@control$Rank) > 2)
      stop("can only handle 'Rank' = 1 or 2 models")
    M <- if (any(slotNames(object) == "predictors") &&
             is.matrix(object@predictors))
         ncol(object@predictors) else
         object@misc$M
    NOS <- ncol(object@y)
    MSratio <- M / NOS  # First value is g(mean) = quadratic form in latvar
    n <- object@misc$n
    colx2.index <- object@control$colx2.index
    cx1i <- object@control$colx1.index
    if (!length(which.species))
      which.species <- 1:NOS
    if (check.ok)
      if (!(length(cx1i) == 1 && names(cx1i) == "(Intercept)"))
          stop("latent variable plots allowable only ",
               "for 'noRRR = ~ 1' models")

    Coeflist <- Coef(object)
    Cmat <- Coeflist@C
    latvarmat <- Coeflist@latvar  # n x Rank

    if (!show.plot)
      return(latvarmat)

    r.curves <- slot(object, type)

    if (MSratio != 1 && type == "predictors")
      stop("can only plot the predictors if M == S")
    MorS <- ncol(r.curves)  # Actually, here, the value is S always.
    if (!add) {
      if (Rank == 1) {
        matplot(latvarmat,
                if ( y && type == "fitted.values")
                    object@y[, which.species, drop = FALSE] else
                    r.curves[, which.species, drop = FALSE],
                type = "n", xlab = xlab, ylab = ylab, ...)
      } else { # Rank == 2
        matplot(c(Coeflist@Optimum[1, which.species], latvarmat[, 1]),
                c(Coeflist@Optimum[2, which.species], latvarmat[, 2]),
                type = "n", xlab = xlab, ylab = ylab, ...)
      }
    }


    pch     <- rep_len(pch,     length(which.species))
    pcol    <- rep_len(pcol,    length(which.species))
    pcex    <- rep_len(pcex,    length(which.species))
    llty    <- rep_len(llty,    length(which.species))
    lcol    <- rep_len(lcol,    length(which.species))
    llwd    <- rep_len(llwd,    length(which.species))
    adj.arg <- rep_len(adj.arg, length(which.species))

    sppnames <- if (type == "predictors") dimnames(r.curves)[[2]] else
                                          dimnames(object@y)[[2]]
    if (Rank == 1) {
      for (sppno in seq_along(which.species)) {
        thisSpecies <- which.species[sppno]
        indexSpecies <- if (is.character(which.species))
           match(which.species[sppno], sppnames) else which.species[sppno]
        if (is.na(indexSpecies))
          stop("mismatch found in 'which.species'")
        xx <- latvarmat
        yy <- r.curves[, indexSpecies]
        ooo <- sort.list(xx)
        xx <- xx[ooo]
        yy <- yy[ooo]
        lines(xx, yy, col = lcol[sppno],
              lwd = llwd[sppno], lty = llty[sppno])
        if (y && type == "fitted.values") {
          ypts <- object@y
          if (ncol(as.matrix(ypts)) == ncol(r.curves))
            points(xx, ypts[ooo, sppno], col = pcol[sppno],
                   cex = pcex[sppno], pch = pch[sppno])
        }
      }
      if (rugplot) rug(xx)
    } else {
      if (sites) {
        text(latvarmat[,1], latvarmat[,2], adj = 0.5,
             labels = if (is.null(spch)) dimnames(latvarmat)[[1]] else
             rep_len(spch, nrow(latvarmat)),
             col = scol, cex = scex, font=sfont)
      }
      for (sppno in seq_along(which.species)) {
          thisSpecies <- which.species[sppno]
          indexSpecies <- if (is.character(which.species))
               match(which.species[sppno], sppnames) else
               which.species[sppno]
          if (is.na(indexSpecies))
            stop("mismatch found in 'which.species'")
          points(Coeflist@Optimum[1, indexSpecies],
                 Coeflist@Optimum[2, indexSpecies],
                 col = pcol[sppno], cex = pcex[sppno], pch = pch[sppno])
      }
      if (label.arg) {
        for (sppno in seq_along(which.species)) {
          thisSpecies <- which.species[sppno]
          indexSpecies <- if (is.character(which.species))
             match(which.species[sppno], sppnames) else
                   which.species[sppno]
          text(Coeflist@Optimum[1, indexSpecies],
               Coeflist@Optimum[2, indexSpecies],
               labels = (dimnames(Coeflist@Optimum)[[2]])[indexSpecies],
               adj = adj.arg[sppno], col = pcol[sppno],
               cex = pcex[sppno])
        }
      }
    }
    invisible(latvarmat)
}


setMethod("lvplot", "rrvgam",
           function(object, ...) {
           invisible(lvplot.rrvgam(object, ...))})



predict.rrvgam <- function (object, newdata = NULL,
                         type = c("link", "response", "terms"),
                         deriv = 0, ...) {
  type <- match.arg(type, c("link", "response", "terms"))[1]
  if (type != "link" && deriv != 0)
    stop("Setting deriv = <positive integer> requires type='link'")
  na.act <- object@na.action
  object@na.action <- list()
  ocontrol <- object@control
  nice21 <- (length(ocontrol$colx1.index) == 1) &&
            (names(ocontrol$colx1.index) == "(Intercept)")
  if (!nice21)
    stop("Can only handle 'noRRR = ~ 1'")

  if (!length(newdata) && type == "response" &&
       length(object@fitted.values)) {
    if (length(na.act)) {
      return(napredict(na.act[[1]], object@fitted.values))
    } else {
      return(object@fitted.values)
    }
  }

  if (!length(newdata)) {
    X <- model.matrixvlm(object, type = "lm", ...)
    offset <- object@offset
    tt <- terms(object)
    if (!length(object@x))
      attr(X, "assign") <- attrassignlm(X, tt)
  } else {
    if (is.smart(object) && length(object@smart.prediction)) {
      setup.smart("read", smart.prediction = object@smart.prediction)
    }

    tt <- terms(object)  # 20030811; object@terms$terms
    X <- model.matrix(delete.response(tt), newdata,
                      contrasts = if (length(object@contrasts))
                                  object@contrasts else NULL,
                      xlev = object@xlevels)

    if (nice21 && nrow(X) != nrow(newdata)) {
      as.save <- attr(X, "assign")
      X <- X[rep_len(1, nrow(newdata)),, drop = FALSE]
      dimnames(X) <- list(dimnames(newdata)[[1]], "(Intercept)")
      attr(X, "assign") <- as.save  # Restored
    }

    offset <- if (!is.null(off.num <- attr(tt, "offset"))) {
                eval(attr(tt, "variables")[[off.num+1]], newdata)
              } else if (!is.null(object@offset))
                eval(object@call$offset, newdata)

    if (is.smart(object) && length(object@smart.prediction)) {
      wrapup.smart()
    }

    attr(X, "assign") <- attrassigndefault(X, tt)
  }

    cancoefs <- concoef(object)

    latvarmat <- X[, ocontrol$colx2.index, drop = FALSE] %*% cancoefs

    Rank <- ocontrol$Rank
    NOS <- ncol(object@y)
    sppnames <- dimnames(object@y)[[2]]
    modelno <- ocontrol$modelno  # 1,2,3,5 or 0
    M <- if (any(slotNames(object) == "predictors") &&
             is.matrix(object@predictors))
         ncol(object@predictors) else
         object@misc$M
    MSratio <- M / NOS  # First value is g(mean) = quadratic form in latvar
    if (type == "terms") {
      terms.mat <- matrix(0, nrow(X), Rank*NOS)  # 1st R cols for spp.1, etc.
      interceptvector <- rep_len(0, NOS)
    } else {
      etamat <- matrix(0, nrow(X), M)  # Could contain derivatives
    }
    ind8 <- 1:Rank
    which.species <- 1:NOS  # Do it all for all species
    for (sppno in seq_along(which.species)) {
      thisSpecies <- which.species[sppno]
      indexSpecies <- if (is.character(which.species))
        match(which.species[sppno], sppnames) else which.species[sppno]
      if (is.na(indexSpecies))
        stop("mismatch found in 'which.species'")

     temp345 <-
       predictrrvgam(object, grid = latvarmat, sppno = thisSpecies,
                  Rank = Rank, deriv = deriv, MSratio = MSratio,
                  type = ifelse(type == "response", "link", type))
     if (MSratio == 2) {
       if (any(type == c("link", "response"))) {
         etamat[, 2*sppno-1] <- temp345$yvals
         etamat[, 2*sppno  ] <- temp345$eta2
       } else {
         terms.mat[, ind8] <- temp345
         interceptvector[sppno] <- attr(temp345, "constant")
       }
     } else {
       if (any(type == c("link", "response"))) {
         etamat[, sppno] <- temp345$yvals
       } else {
         terms.mat[, ind8] <- temp345
         interceptvector[sppno] <- attr(temp345, "constant")
       }
     }
     ind8 <- ind8 + Rank
    }

  if (length(offset) && any(offset != 0))
    etamat <- etamat + offset

  if (type == "link") {
    dimnames(etamat) <-
        list(dimnames(X)[[1]],
             if (deriv == 0)
               object@misc$predictors.names else NULL)
    return(etamat)
  } else if (type == "response") {
    fv <- object@family@linkinv(etamat, extra = object@extra)
    dimnames(fv) <- list(dimnames(fv)[[1]],
                         dimnames(object@fitted.values)[[2]])
    return(fv)
  } else {
    attr(terms.mat, "constant") <- interceptvector
    terms.mat
  }
}



setMethod("predict", "rrvgam", function(object, ...)
           predict.rrvgam(object, ...))




predictrrvgam <- function(object, grid, sppno, Rank = 1,
                          deriv = 0, MSratio = 1, type = "link") {
  if (type != "link" && type != "terms")
    stop("'link' must be \"link\" or \"terms\"")
  if (ncol(grid <- as.matrix(grid)) != Rank)
    stop("'grid' must have ", Rank, " columns")
  if (!is.Numeric(1 + deriv, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE))
    stop("'deriv' must be a non-negative integer")
  if (type == "terms" && deriv != 0)
    stop("'deriv' must be 0 when type=\"terms\"")

  temp.b <- object@Bspline[[sppno]]
  if (type == "terms") {
    meanlatvar <- colMeans(grid)
    answer <- matrix(0, nrow(grid), Rank)
  } else {
    nlfunvalues <- 0
  }

  for (rindex in 1:Rank) {
    temp <- temp.b[[rindex]]  # temp is of class "vsmooth.spline.fit"
    nlpart <- predict(temp, grid[, rindex], deriv = deriv)
    yvals <- nlpart$y
    if (type == "terms") {
      answer[, rindex] <- yvals
    } else {
      nlfunvalues <- nlfunvalues + yvals
    }
  }

  lcoef <- object@coefficients  # linear coefs; dont use coef() (== Coef)
  llcoef <- lcoef[(1+(sppno-1)*(MSratio+Rank)):(sppno*(MSratio+Rank))]
  if (type == "terms") {
    interceptvector <- llcoef[1]
    for (rindex in 1:Rank) {
      answer[, rindex] <- answer[, rindex] + (grid[, rindex] -
                          meanlatvar[rindex]) * llcoef[MSratio+rindex]
      interceptvector <- interceptvector +
          meanlatvar[rindex] * llcoef[MSratio+rindex]
    }
  } else {
    linpar <- if (deriv == 0) {
                llcoef[1] + grid %*% llcoef[-(1:MSratio)]
              } else {
                if (deriv == 1) llcoef[MSratio + rindex] else 0
              }
            nlfunvalues <- nlfunvalues + linpar  # Now complete
  }
  if (type == "terms") {
    attr(answer, "constant") <- interceptvector
    answer
  } else {
    list(xvals = grid,
         yvals = c(nlfunvalues),
         eta2  = if (MSratio == 2) llcoef[MSratio] else NULL)
    }
}






plot.rrvgam <- function(x,
                     xlab = if (Rank == 1) "Latent Variable" else
                            paste("Latent Variable", 1:Rank),
                     ylab = NULL, residuals.arg = FALSE,
                     pcol = par()$col, pcex = par()$cex, pch = par()$pch,
                     lcol = par()$col, lwd = par()$lwd, lty = par()$lty,
                     add = FALSE,
                     main = NULL,
                     center.cf = Rank > 1,
                     WhichRank = 1:Rank,
                     which.species = NULL,  # a numeric or character vector
                     rugplot = TRUE, se.arg = FALSE, deriv = 0,
                     scale = 0, ylim = NULL,
                     overlay = FALSE, ...) {
  Rank <- x@control$Rank
  if (!is.logical(center.cf) || length(center.cf) != 1)
    stop("bad input for argument 'center.cf'")
  if (Rank > 1 &&  !center.cf)
    stop("center.cf = TRUE is needed for models with Rank > 1")
  NOS <- ncol(x@y)
  sppnames <- dimnames(x@y)[[2]]
  modelno <- x@control$modelno  # 1,2,3, or 0
  M <- if (any(slotNames(x) == "predictors") &&
         is.matrix(x@predictors)) ncol(x@predictors) else x@misc$M
  if (all((MSratio <- M / NOS) != c(1,2)))
    stop("bad value for 'MSratio'")

  pcol <- rep_len(pcol, Rank*NOS)
  pcex <- rep_len(pcex, Rank*NOS)
  pch  <- rep_len(pch,  Rank*NOS)
  lcol <- rep_len(lcol, Rank*NOS)
  lwd  <- rep_len(lwd,  Rank*NOS)
  lty  <- rep_len(lty,  Rank*NOS)
  xlab <- rep_len(xlab, Rank)

  if (!length(which.species)) which.species <- 1:NOS
  if (length(ylab))
    ylab <- rep_len(ylab, length(which.species))  # Too long if overlay
  if (length(main))
    main <- rep_len(main, length(which.species))  # Too long if overlay
  latvarmat <- latvar(x)
  nice21 <- length(x@control$colx1.index) == 1 &&
                   names(x@control$colx1.index) == "(Intercept)"
  if (!nice21)
    stop("can only handle intercept-only models")

  counter <- 0
  for (sppno in seq_along(which.species)) {
    thisSpecies <- which.species[sppno]
    indexSpecies <- if (is.character(which.species))
        match(which.species[sppno], sppnames) else which.species[sppno]
    if (is.na(indexSpecies))
      stop("mismatch found in 'which.species'")
    terms.mat <- predictrrvgam(object = x, grid = latvarmat,
                               type = "terms",
                               sppno = indexSpecies, Rank = Rank,
                               deriv = deriv, MSratio = MSratio)
    for (rindex in WhichRank) {
      xvals <- latvarmat[, rindex]
      yvals <- terms.mat[, rindex]
      ooo <- sort.list(xvals)
      xvals <- xvals[ooo]
      yvals <- yvals[ooo]
      if (!center.cf)
        yvals <- yvals + attr(terms.mat, "constant")
      if (!add)
      if (sppno == 1 || !overlay) {
        ylim.use <- if (length(ylim)) ylim else
                    ylim.scale(range(yvals), scale)
        matplot(xvals, yvals, type = "n",
                xlab = xlab[rindex],
                ylab = if (length(ylab)) ylab[sppno] else
                       ifelse(overlay, "Fitted functions",
                                       "Fitted function"),
                main = if (length(main)) main[sppno] else
                       ifelse(overlay, "", sppnames[thisSpecies]),
                ylim = ylim.use,
                ...)
      }
      if (residuals.arg) {
        stop("cannot handle residuals = TRUE yet")
      }
      counter <- counter + 1
      lines(xvals, yvals,
            col = lcol[counter], lwd = lwd[counter], lty = lty[counter])
      if (rugplot) rug(xvals)
    }
  }
  invisible(x)
}




setMethod("plot", "rrvgam",
           function(x, y, ...) {
           if (!missing(y)) stop("cannot process the 'y' argument")
           invisible(plot.rrvgam(x, ...))})



persp.rrvgam <-
  function(x,
           show.plot = TRUE,
           xlim = NULL, ylim = NULL, zlim = NULL,  # zlim ignored if Rank == 1
           gridlength = if (Rank == 1) 301 else c(51, 51),
           which.species = NULL,
           xlab = if (Rank == 1) "Latent Variable" else "Latent Variable 1",
           ylab = if (Rank == 1) "Expected Value"  else "Latent Variable 2",
           zlab = "Expected value",
           labelSpecies = FALSE,   # For Rank == 1 only
           stretch = 1.05,  # quick and dirty, Rank == 1 only
           main = "",
           ticktype = "detailed",
           col = if (Rank == 1) par()$col else "white",
           lty = par()$lty,
           lwd = par()$lwd,
           rugplot = FALSE,
           ...) {
  object <- x  # don't like x as the primary argument
  coefobj <- Coef(object)
  if ((Rank <- coefobj@Rank) > 2)
    stop("object must be a rank-1 or rank-2 model")
  fvmat <- fitted(object)
  NOS <- ncol(fvmat)    # Number of species
  M <- if (any(slotNames(object) == "predictors") &&
         is.matrix(object@predictors)) ncol(object@predictors) else
         object@misc$M
  MSratio <- M / NOS  # First value is g(mean) = quadratic form in latvar

  xlim <- if (length(xlim))
            xlim else
            range(coefobj@latvar[, 1])
  if (!length(ylim.orig <- ylim)) {
    ylim <- if (Rank == 1)
              c(0, max(fvmat)*stretch) else
              range(coefobj@latvar[,2])
  }
  xlim <- rep_len(xlim, 2)
  ylim <- rep_len(ylim, 2)
  gridlength <- rep_len(gridlength, Rank)
  latvar1 <- seq(xlim[1], xlim[2], length = gridlength[1])
  latvar2 <- if (Rank == 2)
               seq(ylim[1], ylim[2], len = gridlength[2]) else
               NULL
    latvarmat <- if (Rank == 2)
                   expand.grid(latvar1, latvar2) else
                   cbind(latvar1)

  sppNames <- dimnames(object@y)[[2]]
  if (!length(which.species)) {
    which.species <- sppNames[1:NOS]
    which.species.numer <- 1:NOS
  } else
  if (is.numeric(which.species)) {
    which.species.numer <- which.species
    which.species <- sppNames[which.species.numer]  # Convert to character
  } else {
    which.species.numer <- match(which.species, sppNames)
  }

  LP <- matrix(NA_real_, nrow(latvarmat), NOS)
  for (sppno in 1:NOS) {
    temp <- predictrrvgam(object = object, grid = latvarmat, sppno = sppno,
                          Rank = Rank, deriv = 0, MSratio = MSratio)
    LP[, sppno] <- temp$yval
  }
  if (MSratio == 2) {
    LP <- kronecker(LP, matrix(1:0, 1, 2))  # n x M
  }
  fitvals <- object@family@linkinv(LP, extra = object@extra)  # n by NOS
  dimnames(fitvals) <- list(NULL, dimnames(fvmat)[[2]])

  if (Rank == 1) {
    if (show.plot) {
      if (!length(ylim.orig))
        ylim <- c(0, max(fitvals[,which.species.numer]) * stretch)  # A revision
      col <- rep_len(col, length(which.species.numer))
      lty <- rep_len(lty, length(which.species.numer))
      lwd <- rep_len(lwd, length(which.species.numer))
      matplot(latvar1, fitvals, xlab = xlab, ylab = ylab,
              type = "n", main = main, xlim = xlim, ylim = ylim, ...)
      if (rugplot) rug(latvar(object))
      for (sppno in seq_along(which.species.numer)) {
        ptr2 <- which.species.numer[sppno]  # points to species column
        lines(latvar1, fitvals[,ptr2], col = col[sppno],
              lty = lty[sppno], lwd = lwd [sppno], ...)
        if (labelSpecies) {
          ptr1 <- (1:nrow(fitvals))[max(fitvals[, ptr2]) ==
                                        fitvals[, ptr2]]
          ptr1 <- ptr1[1]
          text(latvar1[ptr1], fitvals[ptr1, ptr2] + (stretch-1) *
               diff(range(ylim)), label = sppNames[sppno],
               col = col[sppno], ...)
        }
      }
    }
  } else {
    max.fitted <- matrix(fitvals[,which.species[1]],
                         length(latvar1), length(latvar2))
    if (length(which.species) > 1)
      for (sppno in which.species[-1]) {
        max.fitted <- pmax(max.fitted,
                           matrix(fitvals[, sppno],
                                  length(latvar1), length(latvar2)))
    }
    if (!length(zlim))
      zlim <- range(max.fitted, na.rm = TRUE)


    perspdefault <- getS3method("persp", "default")
    if (show.plot)
      perspdefault(latvar1, latvar2, max.fitted,
                   zlim = zlim,
                   xlab = xlab, ylab = ylab, zlab = zlab,
                   ticktype = ticktype, col = col, main = main, ...)


  }

  invisible(list(fitted      = fitvals,
                 latvar1grid = latvar1,
                 latvar2grid = if (Rank == 2) latvar2 else NULL,
                 max.fitted  = if (Rank == 2) max.fitted else NULL))
}


if (!isGeneric("persp"))
  setGeneric("persp", function(x, ...) standardGeneric("persp"))
setMethod("persp", "rrvgam", function(x, ...) persp.rrvgam(x = x, ...))



latvar.rrvgam <- function(object, ...) {
  Coef(object, ...)@latvar
}




if (!isGeneric("lv"))
  setGeneric("lv",
             function(object, ...) {
    .Deprecated("latvar")

               standardGeneric("lv")
             },
             package = "VGAM")

 setMethod("lv", "rrvgam",
           function(object, ...) latvar.rrvgam(object, ...))



 if (!isGeneric("latvar"))
    setGeneric("latvar",
  function(object, ...) standardGeneric("latvar"))

setMethod("latvar", "rrvgam",
  function(object, ...) latvar.rrvgam(object, ...))












setClass(Class = "summary.rrvgam",
         representation("misc" = "list",
                        "call" = "call"),
         contains = "Coef.rrvgam")





summary.rrvgam <- function(object, ...) {
  answer <- Coef(object, ...)


  answer <- as(answer, "summary.rrvgam")


  answer@misc <- object@misc
  answer@call <- object@call
  answer
}


setMethod("summary", "rrvgam", function(object, ...)
  summary.rrvgam(object, ...))




show.summary.rrvgam <- function(x, ...) {
  cat("\nCall:\n")
  dput(x@call)

  show.Coef.rrvgam(x, ...)

  cat("\nNumber of species: ", x@NOS, "\n")

  if (length(x@misc$dispersion) == 1) {
    cat("\nDispersion parameter(s): ", x@misc$dispersion, "\n")
  } else if (is.Numeric(x@dispersion)) {
    cat("\nDispersion parameter(s)\n")
    print( x@dispersion, ... )
  }
  invisible(x)
}





setMethod("show", "summary.rrvgam",
          function(object)
          show.summary.rrvgam(object))




concoef.rrvgam <- function(object, ...) {
  Coef(object, ...)@C
}


concoef.Coef.rrvgam <- function(object, ...) {
  if (length(list(...)))
    warning("Too late! Ignoring the extra arguments")
  object@C
}


if (FALSE) {
 if (!isGeneric("ccoef"))
     setGeneric("ccoef", function(object, ...) {
    .Deprecated("concoef")

    standardGeneric("ccoef")
    })


setMethod("ccoef", "rrvgam", function(object, ...)
    concoef.rrvgam(object, ...))
setMethod("ccoef", "Coef.rrvgam", function(object, ...)
    concoef.Coef.rrvgam(object, ...))
}



setMethod("concoef", "rrvgam", function(object, ...)
    concoef.rrvgam(object, ...))
setMethod("concoef", "Coef.rrvgam", function(object, ...)
    concoef.Coef.rrvgam(object, ...))










if (!isGeneric("calibrate"))
  setGeneric("calibrate", function(object, ...)
  standardGeneric("calibrate"))


setMethod("calibrate", "rrvgam", function(object, ...)
          calibrate.qrrvglm(object, ...))


setMethod("calibrate", "qrrvglm", function(object, ...)
          calibrate.qrrvglm(object, ...))


Tol.rrvgam <- function(object, ...) {
  stop("The tolerance for a 'rrvgam' object is undefined")
}

if (!isGeneric("Tol"))
  setGeneric("Tol", function(object, ...) standardGeneric("Tol"))
setMethod("Tol", "rrvgam", function(object, ...)
          Tol.rrvgam(object, ...))








setMethod("show",  "rrvgam", function(object) show.vgam(object))












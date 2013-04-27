# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.





callcqoc <- function(cmatrix, etamat, xmat, ymat, wvec,
                    X_vlm_1save, modelno, Control,
                    n, M, p1star, p2star, nice31, allofit = FALSE) {
    ocmatrix <- cmatrix
    control <- Control
    Rank <- control$Rank
    p1 <- length(control$colx1.index); p2 <- length(control$colx2.index)
    dim(cmatrix) <- c(p2, Rank)  # for crow1C
    pstar <- p1star + p2star
    maxMr <- max(M, Rank)
    nstar <- if (nice31) ifelse(modelno == 3 || modelno == 5,n*2,n) else n*M
    NOS <- ifelse(modelno == 3 || modelno==5, M/2, M)
    lenbeta <- pstar * ifelse(nice31, NOS, 1)

    if (itol <- control$ITolerances) {
        if (Rank > 1) {
            numat <- xmat[,control$colx2.index,drop = FALSE] %*% cmatrix
            evnu <- eigen(var(numat))
            cmatrix <- cmatrix %*% evnu$vector
        }

        cmatrix <- crow1C(cmatrix, control$Crow1positive)
        numat <- xmat[, control$colx2.index, drop = FALSE] %*% cmatrix
        sdnumat <- apply(numat, 2, sd)
        for(lookat in 1:Rank)
            if (sdnumat[lookat] >
                control$MUXfactor[lookat] * control$isdlv[lookat]) {
                muxer <- control$isdlv[lookat] *
                        control$MUXfactor[lookat] / sdnumat[lookat]
                numat[,lookat] <- numat[,lookat] * muxer
                cmatrix[,lookat] <- cmatrix[,lookat] * muxer
                if (control$trace) {
                    cat(paste("Taking evasive action for latent variable ",
                              lookat, ".\n", sep = ""))
                    flush.console()
                }
                rmfromVGAMenv(c("etamat", "z", "U", "beta", "deviance",
                                "cmatrix", "ocmatrix"), prefix = ".VGAM.CQO.")
            }
    } else {
        numat <- xmat[, control$colx2.index, drop = FALSE] %*% cmatrix
        evnu <- eigen(var(numat))
        temp7 <- if (Rank > 1) evnu$vector %*% diag(evnu$value^(-0.5)) else
                evnu$vector %*% evnu$value^(-0.5)
        cmatrix <- cmatrix %*% temp7
        cmatrix <- crow1C(cmatrix, control$Crow1positive)
        numat <- xmat[, control$colx2.index, drop = FALSE] %*% cmatrix
    }

    inited <- if (is.R()) {
        if (exists(".VGAM.CQO.etamat", envir = VGAM:::VGAMenv)) 1 else 0
    } else 0


    usethiseta <- if (inited == 1) 
        getfromVGAMenv("etamat", prefix = ".VGAM.CQO.") else t(etamat)
    usethisbeta <- if (inited == 2) 
        getfromVGAMenv("beta", prefix = ".VGAM.CQO.") else double(lenbeta)

    othint <- c(Rank = Rank, control$EqualTol, pstar = pstar,
               dimw = 1, inited = inited, modelno = modelno,
               maxitl = control$maxitl, actnits = 0, twice = 0,
               p1star = p1star, p2star = p2star, nice31 = nice31,
               lenbeta = lenbeta, itol = itol, control$trace,
               p1 = p1, p2 = p2, control$imethod)
    bnumat <- if (nice31) matrix(0,nstar,pstar) else
             cbind(matrix(0,nstar,p2star), X_vlm_1save)

 

    ans1 <- 
    dotC(name = if (nice31) "cqo_1" else "cqo_2",
       numat = as.double(numat), as.double(ymat), 
       as.double(if (p1) xmat[,control$colx1.index] else 999),
       as.double(wvec), etamat = as.double(usethiseta),
       moff = double(if (itol) n else 1),
       fv = double(NOS*n), z = double(n*M), wz = double(n*M),
       U = double(M*n), bnumat = as.double(bnumat),
       qr = double(nstar*pstar), qraux = double(pstar),
           qpivot = integer(pstar),
       as.integer(n), as.integer(M), NOS = as.integer(NOS),
           as.integer(nstar), dim1U = as.integer(M),
           errcode = integer(1 + NOS), othint = as.integer(othint),
       deviance = double(1+NOS), beta = as.double(usethisbeta),
           othdbl = as.double(c(small = control$SmallNo,
                  epsilon = control$epsilon, .Machine$double.eps,
                  iKvector = rep(control$iKvector, len = NOS),
                  iShape = rep(control$iShape, len = NOS))))


    if (ans1$errcode[1] == 0) {
        assign2VGAMenv(c("etamat", "z", "U", "beta", "deviance"),
                            ans1, prefix = ".VGAM.CQO.")
        if (is.R()) {
            assign(".VGAM.CQO.cmatrix",   cmatrix, envir = VGAM:::VGAMenv)
            assign(".VGAM.CQO.ocmatrix", ocmatrix, envir = VGAM:::VGAMenv)
        } else {
        }
    } else {
 print("hi 88 20100402; all the species did not converge in callcqo")
        warning("error code in callcqoc = ", ans1$errcode[1])
    if (nice31) {
 print("ans1$errcode[-1]") # Only if (nice31)
 print( ans1$errcode[-1] )
    }
        rmfromVGAMenv(c("etamat", "z", "U", "beta", "deviance",
                        "cmatrix", "ocmatrix"), prefix = ".VGAM.CQO.")
    }
    if (control$trace)
        flush.console()
    if (allofit) list(deviance     = ans1$deviance[1],
                      alldeviance  = ans1$deviance[-1],
                      coefficients = ans1$beta) else ans1$deviance[1]
}



calldcqo <- function(cmatrix, etamat, xmat, ymat, wvec,
                     X_vlm_1save, modelno, Control,
                     n, M, p1star, p2star, nice31, allofit = FALSE) {
    control <- Control
    Rank <- control$Rank
    p1 <- length(control$colx1.index); p2 <- length(control$colx2.index)
    dim(cmatrix) <- c(p2, Rank)  # for crow1C

    xmat2 <- xmat[, control$colx2.index, drop = FALSE]   #ccc
    numat <- double(n*Rank)  #ccc
    pstar <- p1star + p2star
    maxMr <- max(M, Rank)
    nstar <- if (nice31)
             ifelse(modelno == 3 || modelno == 5,n*2,n) else n*M
    NOS <- ifelse(modelno == 3 || modelno == 5, M/2, M)
    lenbeta <- pstar * ifelse(nice31, NOS, 1)

    if (itol <- control$ITolerances) {
        if (Rank > 1) {
            numat <- xmat[, control$colx2.index, drop=FALSE] %*% cmatrix
            evnu <- eigen(var(numat))
            cmatrix <- cmatrix %*% evnu$vector
        }

        cmatrix <- crow1C(cmatrix, control$Crow1positive)
        numat <- xmat[,control$colx2.index,drop=FALSE] %*% cmatrix
        sdnumat <- apply(numat, 2, sd)
        for(lookat in 1:Rank)
          if (sdnumat[lookat] > control$MUXfactor[lookat] *
                                control$isdlv[lookat]) {
                muxer <- control$isdlv[lookat] *
                        control$MUXfactor[lookat] / sdnumat[lookat]
                cmatrix[,lookat] <- cmatrix[,lookat] * muxer
                if (control$trace) {
                    cat(paste("Taking evasive action for latent variable ",
                              lookat, ".\n", sep=""))
                    flush.console()
                }
                rmfromVGAMenv(c("etamat", "z", "U", "beta", "deviance",
                                "cmatrix", "ocmatrix"), prefix = ".VGAM.CQO.")
          }
    } else {
        numat <- xmat[,control$colx2.index,drop=FALSE] %*% cmatrix
        evnu <- eigen(var(numat))
        temp7 <- if (Rank > 1) evnu$vector %*% diag(evnu$value^(-0.5)) else
                              evnu$vector %*% evnu$value^(-0.5)
        cmatrix <- cmatrix %*% temp7
        cmatrix <- crow1C(cmatrix, control$Crow1positive)
        numat <- xmat[, control$colx2.index, drop = FALSE] %*% cmatrix
    }

    inited <- if (is.R()) {
        if (exists(".VGAM.CQO.etamat", envir = VGAM:::VGAMenv)) 1 else 0
    } else 0


    usethiseta <- if (inited == 1) 
        getfromVGAMenv("etamat", prefix = ".VGAM.CQO.") else t(etamat)
    usethisbeta <- if (inited == 2) 
        getfromVGAMenv("beta", prefix = ".VGAM.CQO.") else double(lenbeta)

    othint <- c(Rank, control$EqualTol, pstar, dimw = 1, inited = inited,
               modelno, maxitl = control$maxitl, actnits = 0, twice = 0, 
               p1star = p1star, p2star = p2star, nice31 = nice31, lenbeta,
               itol = itol, control$trace,
               p1, p2, control$imethod) # other ints
    bnumat <- if (nice31) matrix(0,nstar,pstar) else
             cbind(matrix(0,nstar,p2star), X_vlm_1save)
    flush.console()

    ans1 <- 
    dotC(name = "dcqo1", numat = as.double(numat), as.double(ymat), 
       as.double(if (p1) xmat[,control$colx1.index] else 999),
       as.double(wvec), etamat = as.double(usethiseta),
           moff = double(if (itol) n else 1),
           fv = double(NOS*n), z = double(n*M), wz = double(n*M),
           U = double(M*n), bnumat = as.double(bnumat),
       qr = double(nstar * pstar), qraux = double(pstar), qpivot = integer(pstar),
       as.integer(n), as.integer(M), NOS = as.integer(NOS),
       as.integer(nstar), dim1U = as.integer(M),
           errcode = integer(1 + NOS), othint = as.integer(othint),
       deviance = double(1 + NOS), beta = as.double(usethisbeta),
       othdbl = as.double(c(small = control$SmallNo,
                epsilon = control$epsilon, .Machine$double.eps,
                iKvector = rep(control$iKvector, len = NOS),
                iShape = rep(control$iShape, len = NOS))),
       xmat2 = as.double(xmat2),
           cmat = as.double(cmatrix),
       p2 = as.integer(p2), deriv = double(p2*Rank),
           hstep = as.double(control$Hstep))

    if (ans1$errcode[1] != 0) {
        warning("error code in calldcqo = ", ans1$errcode[1])
 print("hi 88 20100402; all the species did not converge in calldcqo")
 print("ans1$errcode[]")
 print( ans1$errcode[] )
    }

    flush.console()
    ans1$deriv
}


checkCMCO <- function(Blist, control, modelno) {

  p1 <- length(colx1.index <- control$colx1.index)
  p2 <- length(colx2.index <- control$colx2.index)
  if (p1 + p2 != length(Blist))
    stop("'Blist' is the wrong length")
  if (p1 == 0 || p2 == 0)
    stop("Some variables are needed in noRRR and non-noRRR arguments")
  if (all(names(colx1.index) != "(Intercept)"))
    stop("an intercept term must be in the argument 'noRRR' formula")
  Blist1 <- vector("list", p1) 
  Blist2 <- vector("list", p2)
  for(kk in 1:p1)
    Blist1[[kk]] <- Blist[[(colx1.index[kk])]]
  for(kk in 1:p2)
    Blist2[[kk]] <- Blist[[(colx2.index[kk])]]

  if (modelno == 3 || modelno == 5) {
    if (p1 > 1)
      for(kk in 2:p1)
        Blist1[[kk]] <- (Blist1[[kk]])[c(TRUE,FALSE),,drop = FALSE]
    for(kk in 1:p2)
      Blist2[[kk]] <- (Blist2[[kk]])[c(TRUE,FALSE),,drop = FALSE]
  }

  if (!all(trivial.constraints(Blist2) == 1))
      stop("the constraint matrices for the non-noRRR terms ",
           "are not trivial")
    if (!trivial.constraints(Blist1[[1]]))
        stop("the constraint matrices for intercept term is ",
             "not trivial")
    if (p1 > 1)
        for(kk in 2:p1)
            if (!trivial.constraints(list(Blist1[[kk]])))
                stop("the constraint matrices for some 'noRRR' ",
                     "terms is not trivial")
            
  nice31 <- if (control$Quadratic)
              (!control$EqualTol || control$ITolerances) else TRUE
  as.numeric(nice31)
}



cqo.fit <- function(x, y, w = rep(1, length(x[, 1])),
    etastart = NULL, mustart = NULL, coefstart = NULL,
    offset = 0, family,
    control = qrrvglm.control(...),
    constraints = NULL,
    extra = NULL,
    Terms = Terms, function.name = "cqo", ...) {




  modelno <- quasi.newton <- NOS <- z <- fv <- NULL





  if (!all(offset == 0))
    stop("cqo.fit() cannot handle offsets")
  specialCM <- NULL
  post <- list()
  nonparametric <- FALSE
  epsilon <- control$epsilon
  maxitl <- control$maxitl
  save.weight <- control$save.weight
  trace <- control$trace
  orig.stepsize <- control$stepsize


  n <- dim(x)[1]


  intercept.only <- ncol(x) == 1 && dimnames(x)[[2]] == "(Intercept)"
  y.names <- predictors.names <- NULL  # May be overwritten in @initialize

 
    n.save <- n 



    Rank <- control$Rank
    rrcontrol <- control  #

    if (length(family@initialize))
        eval(family@initialize)     # Initialize mu and M (and optionally w)
    n <- n.save 

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

    if (is.character(rrcontrol$Dzero)) {
        index <- match(rrcontrol$Dzero, dimnames(as.matrix(y))[[2]]) 
        if (any(is.na(index)))
            stop("Dzero argument didn't fully match y-names")
        if (length(index) == M)
            stop("all linear predictors are linear in the",
                 " latent variable(s); so set 'Quadratic=FALSE'")
        rrcontrol$Dzero <- control$Dzero <- index
    }

    if (length(family@constraints))
        eval(family@constraints)


    special.matrix <- matrix(-34956.125, M, M)    # An unlikely used matrix
    just.testing <- cm.vgam(special.matrix, x, rrcontrol$noRRR, constraints)
    findex <- trivial.constraints(just.testing, special.matrix)
    tc1 <- trivial.constraints(constraints)

    if (!control$Quadratic && sum(!tc1)) {
        for(ii in names(tc1))
            if (!tc1[ii] && !any(ii == names(findex)[findex == 1]))
              warning("'", ii, "' is a non-trivial constraint that will ",
                      "be overwritten by reduced-rank regression")
    }

    if (all(findex == 1))
        stop("use vglm(), not rrvglm()!")
    colx1.index <- names.colx1.index <- NULL
    dx2 <- dimnames(x)[[2]]
    if (sum(findex)) {
        asx <- attr(x, "assign")
        for(ii in names(findex))
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
    p1 <- length(colx1.index); p2 <- length(colx2.index)
    rrcontrol$colx2.index <- control$colx2.index <- colx2.index




    Amat <- if (length(rrcontrol$Ainit)) rrcontrol$Ainit else
            matrix(rnorm(M * Rank, sd = rrcontrol$SD.Cinit), M, Rank)

    Cmat <- if (length(rrcontrol$Cinit)) {
               matrix(rrcontrol$Cinit, p2, Rank)
           } else {
                if (!rrcontrol$Use.Init.Poisson.QO) {
                  matrix(rnorm(p2 * Rank, sd = rrcontrol$SD.Cinit), p2, Rank)
                } else
                  .Init.Poisson.QO(ymat = as.matrix(y), 
                      X1 = x[, colx1.index, drop = FALSE],
                      X2 = x[, colx2.index, drop = FALSE],
                      Rank = rrcontrol$Rank, trace = rrcontrol$trace,
                      max.ncol.etamat = rrcontrol$Etamat.colmax,
                      Crow1positive = rrcontrol$Crow1positive,
                      isdlv = rrcontrol$isdlv,
                      constwt = any(family@vfamily[1] ==
                      c("negbinomial","gamma2","gaussianff")),
                      takelog = any(family@vfamily[1] != c("gaussianff")))
            }

    if (rrcontrol$ITolerances) {
        lvmat <- x[, rrcontrol$colx2.index, drop = FALSE] %*% Cmat
        lvmatmeans <- t(lvmat) %*% matrix(1/n, n, 1)
        if (!all(abs(lvmatmeans) < 4))
            warning("ITolerances=TRUE but the variables making up the ",
                    "latent variable(s) do not appear to be centered.")
    }
    if (modelno == 3 || modelno == 5) 
        Amat[c(FALSE,TRUE),] <- 0  # Intercept only for log(k)

    if (length(control$szero))
        Amat[control$szero,] <- 0

    rrcontrol$Ainit <- control$Ainit <- Amat   # Good for valt()
    rrcontrol$Cinit <- control$Cinit <- Cmat   # Good for valt()

    Blist <- process.constraints(constraints, x, M, specialCM = specialCM)
    nice31 <- checkCMCO(Blist, control = control, modelno = modelno)
    ncolBlist <- unlist(lapply(Blist, ncol))
    dimB <- sum(ncolBlist)

    X_vlm_save <- if (nice31) {
        NULL 
    } else {
        tmp500 <- lm2qrrvlm.model.matrix(x = x, Blist = Blist,
                                        C = Cmat, control = control)
        xsmall.qrr <- tmp500$new.lv.model.matrix 
        B.list <- tmp500$constraints
        lv.mat <- tmp500$lv.mat
        if (length(tmp500$offset)) {
            offset <- tmp500$offset 
        }
        lm2vlm.model.matrix(xsmall.qrr, B.list, xij = control$xij)
    }

    if (length(coefstart) && length(X_vlm_save)) {
        eta <- if (ncol(X_vlm_save) > 1) X_vlm_save %*% coefstart +
                   offset else X_vlm_save * coefstart + offset
        eta <- if (M > 1) matrix(eta, ncol = M, byrow = TRUE) else c(eta) 
        mu <- family@linkinv(eta, extra)
    }

    rmfromVGAMenv(c("etamat", "z", "U", "beta", "deviance",
                    "cmatrix", "ocmatrix"), prefix = ".VGAM.CQO.")

    eval(cqo.init.derivative.expression)
    for(iter in 1:control$optim.maxit) {
        eval(cqo.derivative.expression)
        if (!quasi.newton$convergence) break
    }
    if (maxitl > 1 && iter >= maxitl && quasi.newton$convergence)
        warning("convergence not obtained in", maxitl, "iterations.")

    if (length(family@fini))
        eval(family@fini)

    asgn <- attr(x, "assign")
    coefs <- getfromVGAMenv("beta", prefix = ".VGAM.CQO.")
    if (control$ITolerances) {
        if (NOS == M) {
            coefs <- c(t(matrix(coefs, ncol = M))) # Get into right order
        } else {
            coefs <- coefs
        }
    }

    dn <- labels(x)
    yn <- dn[[1]]
    xn <- dn[[2]]
    residuals <- z - fv
    if (M == 1) {
        residuals <- as.vector(residuals)
        names(residuals) <- yn
    } else {
        dimnames(residuals) <- list(yn, predictors.names)
    }

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
        y.names <- NULL
    }

    df.residual <- 55 - 8 - Rank*p2
    fit <- list(assign = asgn,
                coefficients = coefs,
                constraints = Blist,
                df.residual = df.residual,
                df.total = n*M,
                fitted.values = mu,
                offset = offset, 
                residuals = residuals,
                terms = Terms) # terms: This used to be done in vglm() 

    if (M == 1) {
        wz <- as.vector(wz)  # Convert wz into a vector
    }
    fit$weights <- if (save.weight) wz else NULL

    misc <- list(
        colnames.x = xn,
        criterion = "deviance",
        function.name = function.name, 
        intercept.only=intercept.only,
        predictors.names = predictors.names,
        M = M,
        n = n,
        nonparametric = nonparametric,
        orig.assign = attr(x, "assign"),
        p = ncol(x),
        ynames = dimnames(y)[[2]])

    if (w[1] != 1 || any(w != w[1]))
        fit$prior.weights <- w

    if (length(family@last))
        eval(family@last)

    edeviance <- getfromVGAMenv("deviance", prefix = ".VGAM.CQO.")
    crit.list <- list(deviance = edeviance[1], alldeviance = edeviance[-1])
    if (is.character(y.names) &&
        length(y.names) == length(crit.list$alldeviance))
            names(crit.list$alldeviance) = y.names
    structure(c(fit, list(predictors = matrix(eta, n, M),
        contrasts = attr(x, "contrasts"),
        control = control,
        crit.list = crit.list,
        extra = extra,
        family = family,
        iter = iter,
        misc = misc,
        post = post,
        rss = 000,
        x = x,
        y = y)),
        vclass = family@vfamily)
}



.Init.Poisson.QO <-
  function(ymat, X1, X2, Rank = 1, epsilon = 1/32,
           max.ncol.etamat = 10,
           trace = FALSE, Crow1positive = rep(TRUE, len = Rank),
           isdlv = rep(1, lengt = Rank),
           constwt = FALSE, takelog = TRUE) {

    print.CQO.expression = expression({
        if (trace && length(X2)) {
            cat("\nUsing initial values\n")
            dimnames(ans) <- list(dimnames(X2)[[2]],
                            if (Rank == 1) "lv" else
                            paste("lv", 1:Rank, sep = ""))
            if (p2 > 5) print(ans, dig = 3) else  print(t(ans), dig = 3)
        }
        flush.console()
    })

    sd.scale.X2.expression <- expression({
        if (length(isdlv)) {
            actualSD <- c( sqrt(diag(var(X2 %*% ans))) )
            for(ii in 1:Rank)
                ans[,ii] <- ans[,ii] * isdlv[ii] / actualSD[ii]
        }
    })

    Crow1positive <- if (length(Crow1positive))
        rep(Crow1positive, len = Rank) else
        rep(TRUE, len = Rank)
    if (epsilon <= 0) 
        stop("epsilon > 0 is required")
    ymat <- cbind(ymat) + epsilon  # ymat == 0 cause problems
    NOS <- ncol(ymat)
    p2 <- ncol(X2)
    if (NOS < 2*Rank) {
      ans <- crow1C(matrix(rnorm(p2 * Rank, sd = 0.02), p2, Rank),
                    Crow1positive)
      eval(sd.scale.X2.expression)
      if (NOS == 1) {
        eval(print.CQO.expression) 
        return(ans)
      } else {
        ans.save <- ans;   # ans.save contains scaled guesses
      }
    }

    calS <- 1:NOS  # Set of all species available for the approximation
    effrank <- min(Rank, floor(NOS/2))  # effective rank
    ncol.etamat <- min(if (length(X2)) floor(NOS/2) else effrank,
                      max.ncol.etamat)
    etamat <-
    wts <- matrix(0, nrow = nrow(ymat), ncol = ncol.etamat) # has >=1 coln
    rr <- 1
    for(ii in 1:floor(NOS/2)) {
        if (length(calS) < 2) break
        index <- sample(calS, size = 2)   # Randomness here
        etamat[,rr] <- etamat[,rr] + (if (takelog)
                      log(ymat[,index[1]] / ymat[,index[2]]) else
                          ymat[,index[1]] - ymat[,index[2]])
        wts[,rr] <- wts[,rr] +
                   (if (constwt) 1 else ymat[,index[1]] + ymat[,index[2]])
        calS <- setdiff(calS, index)
        rr <- (rr %% ncol.etamat) + 1
    }
    if (trace)
        cat("\nObtaining initial values\n")

    if (length(X2)) {
        alt <- valt(x = cbind(X1, X2), z = etamat,
                    U = sqrt(t(wts)), Rank = effrank,
                    Blist = NULL, Cinit = NULL, trace = FALSE,
                    colx1.index = 1:ncol(X1), Criterion = "rss")
        temp.control <- list(Rank = effrank, colx1.index = 1:ncol(X1),
                             Alpha = 0.5,
                           colx2.index = (ncol(X1)+1):(ncol(X1) + ncol(X2)),
                             Corner = FALSE, Svd.arg = TRUE,
                             Uncorrelated.lv = TRUE, Quadratic = FALSE)
        
        ans2 <- if (Rank > 1)
               rrr.normalize(rrcontrol = temp.control, A = alt$A, 
                             C = alt$C, x = cbind(X1, X2)) else alt
        ans <- crow1C(ans2$C, rep(Crow1positive, length.out = effrank))

        Rank.save <- Rank
        Rank <- effrank
        eval(sd.scale.X2.expression)
        Rank <- Rank.save 

        if (effrank < Rank) {
            ans <- cbind(ans, ans.save[,-(1:effrank)]) # ans is better
        }
        eval(print.CQO.expression)
    } else {
        xij <- NULL # temporary measure
        U <- t(sqrt(wts))
        tmp <- vlm.wfit(xmat = X1, zmat = etamat, Blist = NULL, U = U,
                       matrix.out = TRUE,
                       is.vlmX = FALSE, rss = TRUE, qr = FALSE, xij = xij)
        ans <- crow1C(as.matrix(tmp$resid),
                     rep(Crow1positive, length.out = effrank))
        if (effrank < Rank) {
            ans <- cbind(ans, ans.save[,-(1:effrank)]) # ans is better
        }

        if (Rank > 1) {
            evnu <- eigen(var(ans))
            ans <- ans %*% evnu$vector
        }

        if (length(isdlv)) {
            actualSD <- apply(cbind(ans), 2, sd)
            for(ii in 1:Rank)
                ans[,ii] <- ans[,ii] * isdlv[ii] / actualSD[ii]
        }
        ans <- crow1C(ans, rep(Crow1positive, length.out = Rank))
        dimnames(ans) <- list(dimnames(X1)[[1]],
                       if (Rank == 1) "lv" else paste("lv", 1:Rank, sep = ""))
        if (trace) {
          if (nrow(ans) > 10) print(t(ans), dig = 3) else
                              print(ans, dig = 3)
        }
    }
    ans
}



cqo.init.derivative.expression <- expression({
    which.optimizer <- if (is.R()) {
        if (control$Quadratic && control$FastAlgorithm) {
          "BFGS" 
        } else {
          if (iter <= rrcontrol$Switch.optimizer) "Nelder-Mead" else "BFGS"
        }
    } else "Quasi-Newton" 
    if (trace && control$OptimizeWrtC) {
        cat("\nUsing", which.optimizer, "algorithm\n")
        flush.console()
    } 


 if (FALSE) {
    constraints <- replace.constraints(constraints, diag(M),
                                      rrcontrol$colx2.index)

    nice31 <- (!control$EqualTol || control$ITolerances) &&
             all(trivial.constraints(constraints) == 1)
}

    NOS <- ifelse(modelno == 3 || modelno == 5, M/2, M)
    canfitok <- if (is.R()) 
        (exists("CQO.FastAlgorithm", envir=VGAM:::VGAMenv) &&
        get("CQO.FastAlgorithm", envir = VGAM:::VGAMenv)) else
    (exists("CQO.FastAlgorithm",inherits=TRUE) && CQO.FastAlgorithm)
    if (!canfitok)
        stop("cannot fit this model using fast algorithm")

    p2star <- if (nice31) 
      ifelse(control$IToleran, Rank, Rank+0.5*Rank*(Rank+1)) else
      (NOS*Rank + Rank*(Rank+1)/2 * ifelse(control$EqualTol,1,NOS))
    p1star <- if (nice31) ifelse(modelno == 3 || modelno == 5,1+p1,p1) else
             (ncol(X_vlm_save)-p2star)
    X_vlm_1save <- if (p1star > 0) X_vlm_save[,-(1:p2star)] else NULL
})
    



cqo.derivative.expression <- expression({


    if (iter == 1 || quasi.newton$convergence) {
        quasi.newton <- optim(par = Cmat, fn = callcqoc,
                gr = if (control$GradientFunction) calldcqo else NULL,
                method = which.optimizer,
                control = list(fnscale = 1,trace = as.integer(control$trace),
                    parscale = rep(control$Parscale, len = length(Cmat)),
                    maxit = control$Maxit.optim),
                etamat = eta, xmat = x, ymat = y, wvec = w,
                X_vlm_1save  =  X_vlm_1save,
                modelno = modelno, Control = control,
                n = n, M = M, p1star = p1star,
                p2star = p2star, nice31 = nice31)

        z <- matrix(getfromVGAMenv("z", prefix = ".VGAM.CQO."), n, M)
        U <- matrix(getfromVGAMenv("U", prefix = ".VGAM.CQO."), M, n)
    }


    ocmatrix <- getfromVGAMenv("ocmatrix", prefix = ".VGAM.CQO.")
    maxdiff <- max(abs(c(ocmatrix) - c(quasi.newton$par)) / (1 +
              abs(c(ocmatrix))))
    if (maxdiff < 1.0e-4) {
        Cmat <- getfromVGAMenv("cmatrix", prefix = ".VGAM.CQO.")
    } else {
        warning("solution does not correspond to .VGAM.CQO.cmatrix")
    }

    alt <- valt.1iter(x = x, z = z, U = U, Blist = Blist,
                      C = Cmat, nice31 = nice31,
                     control = rrcontrol, lp.names = predictors.names,
                     MSratio = M/NOS)

    if (length(alt$offset))
        offset <- alt$offset

    B1.save <- alt$B1 # Put later into extra  
    tmp.fitted <- alt$fitted  # contains \bI_{Rank} \bnu if Corner

    if (trace && control$OptimizeWrtC) {
       cat("\n")
       cat(which.optimizer, "using", if (is.R()) "optim():" else
           "nlminb():", "\n")
       cat("Objective =", quasi.newton$value, "\n")
       cat("Parameters (= c(C)) = ", if (length(quasi.newton$par) < 5) ""
           else "\n")
       cat( if (is.R()) alt$Cmat else format(alt$Cmat), fill=TRUE)
       cat("\n")
       if (!is.R()) {
           cat("Gradient norm =", format(quasi.newton$grad.norm), "\n")
           cat("Number of gradient evaluations =", quasi.newton$g.evals,
               "\n")
       }
       cat("Number of function evaluations =", if (is.R()) 
           quasi.newton$count[1] else quasi.newton$f.evals, "\n")
       if (length(quasi.newton$message))
           cat("Message =", quasi.newton$message, "\n")
       cat("\n")
       flush.console()
    }

    Amat <- alt$Amat  # 
    Cmat <- alt$Cmat  # 
    Dmat <- alt$Dmat  # 

    eval(cqo.end.expression) #
})



cqo.end.expression <- expression({

    rmfromVGAMenv(c("etamat"), prefix = ".VGAM.CQO.")


    if (control$Quadratic) {
        if (!length(extra)) extra =list()
        extra$Amat <- Amat      # Not the latest iteration ??
        extra$Cmat <- Cmat      # Saves the latest iteration 
        extra$Dmat <- Dmat      # Not the latest iteration
        extra$B1   <- B1.save   # Not the latest iteration (not good)
    } else {
        Blist <- replace.constraints(Blist.save, Amat, colx2.index)
    }


    fv <- tmp.fitted            # Contains \bI \bnu
    eta <- fv + offset
    mu <- family@linkinv(eta, extra)

    if (any(is.na(mu)))
        warning("there are NAs in mu") 

    deriv.mu <- eval(family@deriv)
    wz <- eval(family@weight)
    if (control$checkwz)
        wz <- checkwz(wz, M = M, trace = trace, wzeps = control$wzepsilon)
    U <- vchol(wz, M = M, n = n, silent = !trace)
    tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
    z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset # Contains \bI \bnu




})

crow1C <- function(cmat,
                  crow1positive = rep(TRUE, length.out = ncol(cmat)),
                  amat = NULL) {
  if (!is.logical(crow1positive) || length(crow1positive) != ncol(cmat))
    stop("bad input in crow1C")

  for(LV in 1:ncol(cmat))
    if (( crow1positive[LV] && cmat[1,LV] < 0) ||
       (!crow1positive[LV] && cmat[1,LV] > 0)) {
          cmat[,LV] <- -cmat[,LV]
          if (length(amat)) amat[,LV] <- -amat[,LV]
    }
  if (length(amat)) list(cmat = cmat, amat = amat) else cmat
}





printqrrvglm <- function(x, ...)
{
    if (!is.null(cl <- x@call)) {
            cat("Call:\n")
            dput(cl)
    }

    if (FALSE) {
    }

    if (FALSE) {
        nobs <- if (length(x@df.total)) x@df.total else length(x@residuals)
        rdf <- x@df.residual
        if (!length(rdf))
            rdf <- nobs - Rank
    }
    cat("\n")

    if (length(deviance(x)))
        cat("Residual deviance:", format(deviance(x)), "\n")

    if (FALSE && length(x@criterion)) {
      ncrit <- names(x@criterion)
      for(ii in ncrit)
        if (ii != "loglikelihood" && ii != "deviance")
          cat(paste(ii, ":", sep=""), format(x@criterion[[ii]]), "\n")
    }

    invisible(x)
}


setMethod("Coef", "qrrvglm", function(object, ...)
          Coef.qrrvglm(object, ...))


setMethod("coef",         "qrrvglm", function(object, ...)
          Coef.qrrvglm(object, ...))
setMethod("coefficients", "qrrvglm", function(object, ...)
          Coef.qrrvglm(object, ...))


if (!isGeneric("deviance"))
    setGeneric("deviance", function(object, ...)
    standardGeneric("deviance"))
setMethod("deviance", "qrrvglm", function(object,...)
          object@criterion$deviance)


setMethod("fitted",        "qrrvglm", function(object, ...)
          fittedvlm(object))
setMethod("fitted.values", "qrrvglm", function(object, ...)
          fittedvlm(object))









setMethod("show",  "qrrvglm", function(object) printqrrvglm(object))










# These functions are
# Copyright (C) 1998-2011 T.W. Yee, University of Auckland.
# All rights reserved.






uqo.control = function(Rank=1,
          Bestof = if (length(lvstart) && !jitter.sitescores) 1 else 10,
          CA1 = FALSE,
          Crow1positive = TRUE,
          epsilon = 1.0e-07,
          EqualTolerances = ITolerances,
          Etamat.colmax = 10,
          GradientFunction=TRUE,
          Hstep = 0.001,
          isdlv = rep(c(2, 1, rep(0.5, len=Rank)), len=Rank),
          ITolerances = FALSE,
          lvstart = NULL,
          jitter.sitescores = FALSE,
          maxitl = 40,
          Maxit.optim = 250,
          MUXfactor = rep(3, length=Rank),
          optim.maxit = 20,
          nRmax = 250,
          SD.sitescores = 1.0,
          SmallNo = 5.0e-13,
          trace = TRUE,
          Use.Init.Poisson.QO=TRUE,
          ...)
{

    Kinit = 0.001
    if (!is.Numeric(MUXfactor, posit=TRUE))
        stop("bad input for \"MUXfactor\"")
    if (any(MUXfactor < 1 | MUXfactor > 10))
        stop("MUXfactor values must lie between 1 and 10")
    if (!is.Numeric(isdlv, posit=TRUE)) stop("bad input for \"isdlv\"")
    if (any(isdlv < 0.2 | isdlv > 10))
        stop("isdlv values must lie between 0.2 and 10")
    if (length(isdlv) > 1 && any(diff(isdlv) > 0))
        stop("successive isdlv values must not increase")
    if (!is.Numeric(Rank, allow=1, integ=TRUE, posit=TRUE))
        stop("Bad input for \"Rank\"")
    if (!is.Numeric(Bestof, allow=1, integ=TRUE, posit=TRUE))
        stop("Bad input for \"Bestof\"")
    if (!is.Numeric(Etamat.colmax, posit=TRUE, allow=1) || Etamat.colmax < Rank)
        stop("bad input for \"Etamat.colmax\"")
    if (!is.Numeric(maxitl, allow=1, integ=TRUE, posit=TRUE))
        stop("Bad input for \"maxitl\"")
    if (!is.Numeric(Maxit.optim, integ=TRUE, posit=TRUE, allow=1))
        stop("Bad input for \"Maxit.optim\"")
    if (!is.Numeric(optim.maxit, allow=1, integ=TRUE, posit=TRUE))
        stop("Bad input for \"optim.maxit\"")
    if (!is.Numeric(nRmax, allow=1, integ=TRUE, posit=TRUE))
        stop("Bad input for \"nRmax\"")
    if (!is.Numeric(Hstep, allow=1, posit=TRUE))
        stop("Bad input for \"Hstep\"")
    if (!is.Numeric(epsilon, allow=1, posit=TRUE))
        stop("Bad input for \"epsilon\"")
    if (!is.Numeric(SmallNo, allow=1, posit=TRUE))
        stop("Bad input for \"SmallNo\"")

    if ((SmallNo < .Machine$double.eps) || (SmallNo > .0001))
        stop("SmallNo is out of range") 

    if (Use.Init.Poisson.QO && CA1)
        stop("cannot have both Use.Init.Poisson.QO=TRUE and CA1=TRUE")

    ans = list(
           Bestof = Bestof,
           CA1 = CA1,
           ConstrainedQO = FALSE, # A constant, not a control parameter
           Corner = FALSE, # Needed for valt.1iter()
           Crow1positive=as.logical(rep(Crow1positive, len=Rank)),
           epsilon = epsilon,
           EqualTolerances = as.logical(EqualTolerances)[1],
           Etamat.colmax = Etamat.colmax,
           FastAlgorithm = TRUE, # A constant, not a control parameter
           GradientFunction = GradientFunction,
           Hstep = Hstep,
           isdlv = rep(isdlv, len=Rank),
           ITolerances = as.logical(ITolerances)[1],
           lvstart = lvstart,
           jitter.sitescores = as.logical(jitter.sitescores),
           Kinit = Kinit,
           maxitl= maxitl,
           Maxit.optim = Maxit.optim,
           MUXfactor = rep(MUXfactor, length=Rank),
           nRmax = nRmax,
           optim.maxit = optim.maxit,
           OptimizeWrtC = FALSE,
           Quadratic = TRUE,
           Rank = Rank,
           SD.sitescores = SD.sitescores,
           SmallNo = SmallNo,
           trace = as.logical(trace),
           Use.Init.Poisson.QO=as.logical(Use.Init.Poisson.QO)[1])
    ans
}




uqo  <- function(formula,
                 family, data=list(), 
                 weights=NULL, subset=NULL, na.action=na.fail,
                 etastart=NULL, mustart=NULL, coefstart=NULL,
                 control=uqo.control(...), 
                 offset=NULL, 
                 method="uqo.fit",
                 model=FALSE, x.arg=TRUE, y.arg=TRUE,
                 contrasts=NULL, 
                 constraints=NULL,
                 extra=NULL, 
                 qr.arg=FALSE, ...)
{
    dataname <- as.character(substitute(data))  # "list" if no data=
    function.name <- "uqo"

    ocall <- match.call()

    mt <- terms(formula, data = data)
    if (missing(data)) 
        data <- environment(formula)

    mf <- match.call(expand=FALSE)
    mf$family <- mf$method <- mf$model <- mf$x.arg <- mf$y.arg <- mf$control <-
        mf$contrasts <- mf$constraints <- mf$extra <- mf$qr.arg <- NULL
    mf$coefstart <- mf$etastart <- mf$... <- NULL
    mf$drop.unused.levels <- TRUE 
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame()) 
    if (method == "model.frame")
        return(mf)
    na.act <- attr(mf, "na.action")

    xvars <- as.character(attr(mt, "variables"))[-1]
    if ((yvar <- attr(mt, "response")) > 0)
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(mf[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }

    y <- model.response(mf, "numeric") # model.extract(mf, "response")
    x <- model.matrix(mt, mf, contrasts)
    attr(x, "assign") = attrassigndefault(x, mt)
    offset <- model.offset(mf)
    if (is.null(offset)) 
        offset <- 0 # yyy ???
    w <- model.weights(mf)
    if (!length(w))
        w <- rep(1, nrow(mf))
    else if (ncol(as.matrix(w))==1 && any(w < 0))
        stop("negative weights not allowed")

    if (is.character(family))
        family <- get(family)
    if (is.function(family))
        family <- family()
    if (!inherits(family, "vglmff")) {
        stop("'family=", family, "' is not a VGAM family function")
    }

    if (!is.null(family@first))
        eval(family@first)

    uqo.fitter <- get(method)
    if (ncol(x) != 1 || dimnames(x)[[2]] != "(Intercept)")
        stop("uqo()'s formula must have ~ 1 on the RHS") 

    if (control$FastAlgorithm &&
       length(as.list(family@deviance)) <= 1)
        stop("The fast algorithm requires the family ",
             "function to have a deviance slot")
    deviance.Bestof = rep(as.numeric(NA), len=control$Bestof)
    for(tries in 1:control$Bestof) {
         if (control$trace && (control$Bestof>1))
         cat(paste("\n========================= Fitting model", tries,
                     "=========================\n"))
         it <- uqo.fitter(x=x, y=y, w=w, offset=offset,
                   etastart=etastart, mustart=mustart, coefstart=coefstart,
                   family=family, control=control,
                   constraints=constraints, extra=extra,
                   qr.arg = qr.arg, Terms=mt, function.name=function.name,
                   ca1 = control$CA1 && tries==1, ...)
        deviance.Bestof[tries] = it$crit.list$deviance
        if (tries==1||min(deviance.Bestof[1:(tries-1)]) > deviance.Bestof[tries])
            fit = it
    }
    fit$misc$deviance.Bestof = deviance.Bestof
    fit$misc$criterion = "deviance"  # Needed for calibrate; 21/1/05

    fit$misc$dataname <- dataname

    answer <-
    new("uqo",
      "call"         = ocall,
      "coefficients" = fit$coefficients,
      "constraints"  = fit$constraints,
      "control"      = fit$control,
      "criterion"    = fit$crit.list,
      "lv"           = fit$sitescores,
      "family"       = fit$family,
      "fitted.values"= as.matrix(fit$fitted.values),
      "iter"         = fit$iter,
      "misc"         = fit$misc,
      "model"        = if (model) mf else data.frame(),
      "na.action"    = if (length(na.act)) list(na.act) else list(),
      "predictors"  = as.matrix(fit$predictors))

    answer@control$min.criterion = TRUE # Needed for calibrate; 21/1/05

    if (length(fit$weights))
        slot(answer, "weights") = as.matrix(fit$weights)
    if (x.arg)
        slot(answer, "x") = x
    if (y.arg)
        slot(answer, "y") = as.matrix(fit$y)
    slot(answer, "extra") = if (length(fit$extra)) {
        if (is.list(fit$extra)) fit$extra else {
            warning("'extra' is not a list, therefore ",
                    "placing 'extra' into a list")
            list(fit$extra)
        }
    } else list() # R-1.5.0
    if (length(fit$prior.weights))
        slot(answer, "prior.weights") = as.matrix(fit$prior.weights)

    answer
}


calluqof = function(sitescores, etamat, ymat, wvec, modelno, nice31, xmat,
                    Control,
                    n, M, maxMr5, othint, othdbl, bnumat, Hstep=NA, alldump) {
    control = Control
    Rank = control$Rank
    itol = othint[14]
    inited = if (is.R()) {
        as.numeric(existsinVGAMenv("etamat", prefix=".VGAM.UQO."))
    } else 0
    othint[5] = inited  # Replacement
    usethiseta = if (inited==1)
        getfromVGAMenv("etamat", prefix = ".VGAM.UQO.") else t(etamat)
    usethisbeta = double(othint[13])
    pstar = othint[3]
    nstar = if (nice31) ifelse(modelno==3 || modelno==5,n*2,n) else n*M
    NOS = ifelse(modelno == 3 || modelno==5, M/2, M)

    sitescores = matrix(sitescores, ncol=Rank)
    sitescores = scale(sitescores, center = TRUE, scale = FALSE)
    if (itol) {
        numat = matrix(sitescores, ncol=Rank)
        if (Rank > 1) {
            evnu = eigen(var(numat))
            numat = numat %*% evnu$vector
        }



        sdnumat = apply(numat, 2, sd)
        for(lookat in 1:Rank)
            if (sdnumat[lookat]>control$MUXfactor[lookat]*control$isdlv[lookat]){
                muxer = control$isdlv[lookat] * control$MUXfactor[lookat] / 
                        sdnumat[lookat]
                numat[,lookat] = numat[,lookat] * muxer
                if (control$trace) {
                }
            }
    } else {
        numat = matrix(sitescores, ncol=Rank)
        evnu = eigen(var(numat))
        temp7 = if (Rank > 1) evnu$vector %*% diag(evnu$value^(-0.5)) else
                evnu$vector %*% evnu$value^(-0.5)
        numat = numat %*% temp7
    }

    ans1 <- 
    dotFortran(name = if (nice31) "cqo1f" else "cqo2f",
       numat=as.double(numat), as.double(ymat), 
             as.double(xmat),
       as.double(wvec), etamat=as.double(usethiseta),
           moff=double(if(itol) n else 1),
           fv=double(NOS*n), z=double(n*M), wz=double(n*M),
           U=double(M*n), bnumat=as.double(bnumat),
       qr=double(nstar*pstar), qraux=double(pstar), qpivot=integer(pstar),
       as.integer(n), as.integer(M), NOS=as.integer(NOS),
       as.integer(nstar), dimu=as.integer(M),
           errcode=integer(1), othint=as.integer(othint),
           rowind=integer(maxMr5), colind=integer(maxMr5),
       deviance=double(1), beta=as.double(usethisbeta),
       twk=double(if(nice31) nstar*3 else M*n*2), wkmm=double(M*(M+pstar)),
           othdbl=as.double(othdbl))

       if (ans1$errcode == 0) {
            assign2VGAMenv(c("etamat","numat"), ans1, prefix=".VGAM.UQO.")
            if (alldump) {
                ans1$fv = matrix(ans1$fv,n,M,byrow=TRUE,dimnames=dimnames(ymat))
                assign2VGAMenv(c("beta","fv"), ans1, prefix=".VGAM.UQO.")
                assign2VGAMenv(c("z","U"), ans1, prefix=".VGAM.UQO.")
            }
       } else {
           cat("warning in calluqof: error code =", ans1$errcode, "\n")
           rmfromVGAMenv(c("etamat"), prefix=".VGAM.UQO.")
       }
    ans1$deviance
}

callduqof = function(sitescores, etamat, ymat, wvec, modelno, nice31, xmat,
                     Control,
                     n, M, maxMr5, othint, othdbl, bnumat, Hstep, alldump) {
    control = Control
    itol = othint[14]
    inited = if (is.R()) {
        if (exists(".VGAM.UQO.etamat", envir = VGAM:::VGAMenv)) 1 else 0
    } else 0 # 0 means fortran initializes the etamat
    othint[5] = inited  # Replacement
    usethiseta = if (inited==1)
        getfromVGAMenv("etamat", prefix = ".VGAM.UQO.") else t(etamat)
    usethisbeta = double(othint[13])
    pstar = othint[3]
    nstar = if (nice31) ifelse(modelno==3 || modelno==5,n*2,n) else n*M
    NOS = ifelse(modelno == 3 || modelno==5, M/2, M)
    Rank = othint[1]

    sitescores = matrix(sitescores, ncol=Rank)
    sitescores = scale(sitescores, center = TRUE, scale = FALSE)
    if (itol) {
        numat = matrix(sitescores, ncol=Rank)
        if (Rank > 1) {
            evnu = eigen(var(numat))
            numat = numat %*% evnu$vector
        }

        sdnumat = apply(numat, 2, sd)
        for(lookat in 1:Rank)
            if (sdnumat[lookat]>control$MUXfactor[lookat]*control$isdlv[lookat]){
                muxer = control$isdlv[lookat] * control$MUXfactor[lookat] / 
                        sdnumat[lookat]
                numat[,lookat] = numat[,lookat] * muxer
                if (control$trace) {
                }
            }
    } else {
        numat = matrix(sitescores, ncol=Rank)
        evnu = eigen(var(numat))
        temp7 = if (Rank > 1) evnu$vector %*% diag(evnu$value^(-0.5)) else
                evnu$vector %*% evnu$value^(-0.5)
        numat = numat %*% temp7
    }



    ans1 <- 
    dotFortran(name = "duqof", numat=as.double(numat), as.double(ymat),
           as.double(xmat),
           as.double(wvec), etamat=as.double(usethiseta),
           moff=double(if(itol) n else 1),
           fv=double(NOS*n), z=double(n*M), wz=double(n*M),
           U=double(M*n), bnumat=as.double(bnumat),
       qr=double(nstar*pstar), qraux=double(pstar), qpivot=integer(pstar),
       as.integer(n), as.integer(M), NOS=as.integer(M), 
       as.integer(nstar), dimu=as.integer(M),
           errcode=integer(1), othint=as.integer(othint),
           rowind=integer(maxMr5), colind=integer(maxMr5),
       deviance=double(1), beta=as.double(usethisbeta),
       twk=double(if(nice31) nstar*3 else M*n*2), wkmm=double(M*(M+pstar)),
           othdbl=as.double(othdbl),
       onumat=as.double(numat),
       deriv=double(n*Rank), hstep=as.double(Hstep),
       betasave=usethisbeta)

       if (ans1$errcode == 0) {
           assign2VGAMenv(c("etamat"), ans1, prefix=".VGAM.UQO.")
       } else {
           cat("warning in callduqof: error code =", ans1$errcode, "\n")
           rmfromVGAMenv(c("etamat"), prefix=".VGAM.UQO.")
       }
    ans1$deriv
}




uqo.fit <- function(x, y, w=rep(1, len=nrow(x)),
    etastart=NULL, mustart=NULL, coefstart=NULL,
    offset=0, family, control=uqo.control(...),
    qr.arg=FALSE, constraints=NULL, extra=NULL,
    Terms=Terms, function.name="uqo", ca1=TRUE, ...)
{
    if (!all(offset == 0)) stop("cqo.fit() cannot handle offsets")
    nonparametric <- FALSE
    epsilon <- control$epsilon
    optim.maxit <- control$optim.maxit
    save.weight <- control$save.weight
    trace <- control$trace
    orig.stepsize <- control$stepsize


    n <- dim(x)[1]


    copy_X_vlm <- FALSE    # May be overwritten in @initialize
    stepsize <- orig.stepsize
    old.coeffs <- coefstart

    intercept.only <- ncol(x) == 1 && dimnames(x)[[2]] == "(Intercept)"
    y.names <- predictors.names <- NULL    # May be overwritten in @initialize

    n.save <- n 


    Rank <- control$Rank
    rrcontrol <- control  #

    if (length(family@initialize))
        eval(family@initialize)       # Initialize mu and M (and optionally w)
    n <- n.save 


    eval(rrr.init.expression)

    if (length(etastart)) {
        eta <- etastart
        mu <- if (length(mustart)) mustart else family@linkinv(eta, extra)
    } else {
        if (length(mustart))
            mu <- mustart
        eta <- family@link(mu, extra)
    }

    M <- if (is.matrix(eta)) ncol(eta) else 1

    if (is.character(rrcontrol$Dzero)) {
        index = match(rrcontrol$Dzero, dimnames(as.matrix(y))[[2]]) 
        if (any(is.na(index)))
            stop("Dzero argument didn't fully match y-names")
        if (length(index) == M)
            stop("all linear predictors are linear in the latent variable(s)")
        rrcontrol$Dzero = control$Dzero = index
    }



    if (length(family@constraints))
        eval(family@constraints)


    colx1.index = 1:ncol(x)
    names(colx1.index) = names.colx1.index = dimnames(x)[[2]]

    rrcontrol$colx1.index=control$colx1.index=colx1.index #Save it on the object
    colx2.index = NULL
    p1 = length(colx1.index); p2 = 0
    rrcontrol$colx2.index=control$colx2.index=colx2.index #Save it on the object
    rrcontrol$Quadratic = control$Quadratic = TRUE



    sitescores <- if (length(rrcontrol$lvstart)) {
        matrix(rrcontrol$lvstart, n, Rank)
    } else {
        if (rrcontrol$Use.Init.Poisson) {
            .Init.Poisson.QO(ymat=as.matrix(y),
                             X1=x, X2=NULL,
                             Rank=rrcontrol$Rank, trace=rrcontrol$trace,
                             max.ncol.etamat = rrcontrol$Etamat.colmax,
                             Crow1positive=rrcontrol$Crow1positive,
                             isdlv=rrcontrol$isdlv,
               constwt= any(family@vfamily[1] ==
                            c("negbinomial","gamma2","gaussianff")),
               takelog= any(family@vfamily[1] != c("gaussianff")))
        } else if (ca1) {
            if (Rank == 1) .VGAM.UQO.CA(y)[,1:Rank]  else {
                temp = .VGAM.UQO.CA(y)[,1:Rank]
                temp %*% solve(chol(var(temp)))
            }
        } else {
            matrix((runif(n*Rank)-0.5)*rrcontrol$SD.sitescores,n,Rank)
        }
    }
    if (rrcontrol$jitter.sitescores)
        sitescores <- jitteruqo(sitescores)


    Blist <- process.constraints(constraints, x, M)
    ncolBlist <- unlist(lapply(Blist, ncol))
    dimB <- sum(ncolBlist)


    modelno = switch(family@vfamily[1], "poissonff"=2,
              "binomialff"=1, "quasipoissonff"=0, "quasibinomialff"=0,
              "negbinomial"=0,
              "gamma2"=5,
              0)  # stop("can't fit this model using fast algorithm")
    if (!modelno) stop("the family function does not work with uqo()")
    if (modelno == 1) modelno = get("modelno", envir = VGAM:::VGAMenv)
    rmfromVGAMenv(c("etamat", "beta"), prefix=".VGAM.UQO.")

    cqofastok = if (is.R()) (exists("CQO.FastAlgorithm", envir = VGAM:::VGAMenv) &&
                  get("CQO.FastAlgorithm", envir = VGAM:::VGAMenv)) else
              (exists("CQO.FastAlgorithm", inherits=TRUE) && CQO.FastAlgorithm)
    if (!cqofastok)
        stop("can't fit this model using fast algorithm")

    nice31 = (!control$EqualTol || control$ITolerances) && control$Quadratic &&
             all(trivial.constraints(Blist))


    X_vlm_1save <- if (nice31) {
        NULL
    } else {
        lm2vlm.model.matrix(x, Blist, xij=control$xij)
    }

    NOS = ifelse(modelno==3 || modelno==5, M/2, M)
    p1star = if (nice31) p1*ifelse(modelno==3 || modelno==5,2,1) else ncol(X_vlm_1save)
    p2star = if (nice31)
      ifelse(control$IToleran, Rank, Rank+0.5*Rank*(Rank+1)) else
      (NOS*Rank + Rank*(Rank+1)/2 * ifelse(control$EqualTol,1,NOS))

    pstar = p1star + p2star
    nstar = if (nice31) ifelse(modelno==3 || modelno==5,n*2,n) else n*M
    maxMr = max(M, Rank)
    maxMr5 = maxMr*(maxMr+1)/2
    lenbeta = pstar * ifelse(nice31, NOS, 1)

    othint = c(Rank, control$EqualTol, pstar, dimw=1, inited=290, # other ints
               modelno, maxitl=control$maxitl, actnits=0, twice=0, p1star,
               p2star, nice31, lenbeta, control$ITolerances, control$trace,
               p1, p2, control$imethod)
    othdbl = c(small=control$SmallNo, fseps=control$epsilon,
               .Machine$double.eps,
               kinit=rep(control$Kinit, len=NOS),
               shapeinit=rep(control$shapeinit, len=NOS))
    bnumat = if (nice31) matrix(0,nstar,pstar) else
             cbind(matrix(0,nstar,p2star), X_vlm_1save)

    rmfromVGAMenv(c("etamat", "z", "U", "beta", "deviance", "fv",
                         "cmatrix", "ocmatrix"), prefix=".VGAM.UQO.")


    for(iter in 1:optim.maxit) {
        if (control$trace)
            cat("\nIteration", iter, "\n")
        conjgrad <- optim(par=sitescores, fn=calluqof, 
                     gr = if (control$GradientFunction) callduqof else NULL,
                     method = if (n*Rank>control$nRmax) "CG" else "BFGS",
                     control=list(fnscale=1, trace=as.integer(control$trace),
                                  maxit=control$Maxit.optim),
                     etamat=eta, ymat=y, wvec=w, modelno=modelno,
                     Control=rrcontrol,
                     nice31=nice31, xmat = x,
                     n=n, M=M, maxMr5=maxMr5, othint=othint, othdbl=othdbl,
                     bnumat=bnumat, Hstep=control$Hstep, alldump=FALSE)

        sitescores = getfromVGAMenv("numat", prefix = ".VGAM.UQO.")
        dim(sitescores) = c(n, Rank)
        sitescores = scale(sitescores, center = TRUE, scale = FALSE)
        sitescores = crow1C(sitescores, rrcontrol$Crow1positive)
        dimnames(sitescores) = list(dimnames(y)[[1]], if (Rank==1) "lv" else
                                    paste("lv", 1:Rank, sep=""))

        if (converged <- (conjgrad$convergence == 0)) break
    }

    if (!converged && optim.maxit>1)
        warning("convergence not obtained")


    temp9 = 
    calluqof(sitescores, etamat=eta, ymat=y, wvec=w, modelno=modelno,
             nice31=nice31, xmat = x,
             Control=rrcontrol,
             n=n, M=M, maxMr5=maxMr5, othint=othint, othdbl=othdbl,
             bnumat=bnumat, Hstep=NA, alldump=TRUE)

    coefs = getfromVGAMenv("beta", prefix = ".VGAM.UQO.")
    VGAM.fv = getfromVGAMenv("fv", prefix = ".VGAM.UQO.")
    etamat = getfromVGAMenv("etamat", prefix = ".VGAM.UQO.")
    dim(etamat) = c(M,n)
    etamat = t(etamat)
    wresids = getfromVGAMenv("z", prefix = ".VGAM.UQO.") - etamat
    dim(wresids) = c(n,M)


    if (!intercept.only)
        stop("can only handle intercept.only==TRUE currently")
    if (nice31) {
        coefs = c(t(matrix(coefs, ncol=M))) # Get into right order
        coefs = matrix(coefs, nrow=M)
        Amat = coefs[,1:Rank,drop=FALSE]
        if (rrcontrol$IToleran) {
            B1 = coefs[,-(1:Rank),drop=FALSE]
            Dmat = matrix(0, M, Rank*(Rank+1)/2)
            Dmat[,1:Rank] = -0.5
        } else {
            Dmat = coefs[,(Rank+1):(Rank + Rank*(Rank+1)/2),drop=FALSE]
            B1 = coefs[,(1+(Rank + Rank*(Rank+1)/2)):ncol(coefs),drop=FALSE]
        }
    } else {
        Amat = t(matrix(coefs[1:(Rank*M)], Rank, M))
        cptr1 = (Rank*M)
        Dmat = coefs[(cptr1+1):(cptr1+Rank*(Rank+1)/2)]
        Dmat = matrix(Dmat, M, Rank*(Rank+1)/2, byrow=TRUE)
        cptr1 = (Rank*M) + Rank*(Rank+1)/2
        B1 = coefs[(cptr1+1):length(coefs)]
    }

    lv.names = if (Rank==1) "lv" else paste("lv", 1:Rank, sep="") 
    lp.names = predictors.names
    if (!length(lp.names)) lp.names = NULL
    extra$Amat = matrix(Amat, M, Rank, dimnames = list(lp.names, lv.names))
    extra$B1   = matrix(B1, ncol=M, dimnames = 
                        list(names.colx1.index, predictors.names))
    extra$Dmat = matrix(Dmat, M, Rank*(Rank+1)/2)
    extra$Cmat = NULL  # This is UQO!!

    VGAM.etamat = getfromVGAMenv("etamat", prefix = ".VGAM.UQO.") 
    VGAM.etamat = matrix(VGAM.etamat, n, M, byrow=TRUE,
                         dimnames = list(dimnames(y)[[1]], predictors.names))

    coefficients = c(coefs) # Make a vector because of class "numeric"

    rmfromVGAMenv(c("etamat", "beta", "fv"), prefix=".VGAM.UQO.")

    if (length(family@fini))
        eval(family@fini)

    misc <- list(function.name = function.name, 
        intercept.only=intercept.only,
        predictors.names = predictors.names,
        modelno = modelno,
        M = M, n = n,
        nstar = nstar, nice31 = nice31,
        p = ncol(x),
        pstar = pstar, p1star = p1star, p2star = p2star,
        ynames = dimnames(y)[[2]])

    crit.list <- list(deviance=conjgrad$value)


    structure(c(list(
        coefficients = coefficients,
        constraints = Blist,
        sitescores = sitescores,
        crit.list = crit.list,
        control=control,
        extra=extra,
        family=family,
        fitted.values=VGAM.fv,
        iter=iter,
        misc=misc,
        predictors=VGAM.etamat,
        prior.weights = w,
        x=x,
        y=y)),
        vclass=family@vfamily)
}




printuqo <- function(x, ...)
{
    if (!is.null(cl <- x@call)) {
        cat("Call:\n")
        dput(cl)
    }

    cat("\n")
    cat(x@misc$n, "sites and", x@misc$M, "responses/species\n")
    cat("Rank", x@control$Rank)
    cat(",", ifelse(x@control$EqualToler, "equal-tolerances", 
        "unequal-tolerances"), "\n")

    if (length(deviance(x)))
        cat("\nResidual Deviance:", format(deviance(x)), "\n")

    invisible(x)
}


setMethod("print",  "uqo", function(x, ...)  printuqo(x, ...))

    setMethod("show",  "uqo", function(object)  printuqo(object))



deviance.uqo <- function(object, ...)
    object@criterion$deviance

setMethod("deviance", "uqo", function(object, ...)
           deviance.uqo(object, ...))


setMethod("coefficients", "uqo", function(object, ...)
          Coef.qrrvglm(object, ...))
setMethod("coef", "uqo", function(object, ...)
          Coef.qrrvglm(object, ...))
setMethod("Coef", "uqo", function(object, ...)
          Coef.qrrvglm(object, ...))




setMethod("show", "Coef.uqo", function(object)
    printCoef.qrrvglm(object, C = FALSE))
setMethod("print", "Coef.uqo", function(x, ...)
    printCoef.qrrvglm(x, ...))




residualsuqo  <- function(object,
              type = c("deviance", "pearson", "working", "response"),
              matrix.arg= TRUE) {

    if (mode(type) != "character" && mode(type) != "name")
        type = as.character(substitute(type))
    type = match.arg(type, c("deviance", "pearson", "working", "response"))[1]

    switch(type,
        response = object@y - fitted(object),
        stop("this type of residual hasn't been implemented yet")
    )
}

setMethod("resid", "uqo", function(object, ...) 
          residualsuqo(object, ...))
setMethod("residuals", "uqo", function(object, ...)
          residualsuqo(object, ...))

fitted.values.uqo  <- function(object, ...)
    object@fitted.values

setMethod("fitted", "uqo", function(object, ...) 
          fitted.values.uqo(object, ...))
setMethod("fitted.values", "uqo", function(object, ...)
          fitted.values.uqo(object, ...))



predict.uqo  <- function(object, newdata = NULL, ...) {
    if (length(newdata) > 0)
        stop("can't handle newdata argument yet")
    object@predictors
}

setMethod("predict", "uqo", function(object, ...) 
          predict.uqo(object, ...))


setMethod("persp", "uqo", function(x, ...) 
          perspqrrvglm(x, ...))

setMethod("trplot", "uqo", function(object, ...) 
          trplot.qrrvglm(object, check.ok=FALSE, ...))


setMethod("plot", "uqo", function(x, y, ...) 
         invisible(plotqrrvglm(object=x, ...)))


setMethod("lvplot", "uqo", function(object, ...) 
         invisible(lvplot.qrrvglm(object, C=FALSE, check.ok=FALSE, ...)))



.VGAM.UQO.CA = function(Y) {
    Y = as.matrix(Y) / sum(Y)
    rowsum = c(Y %*% rep(1, len=ncol(Y)))
    colsum = c(t(Y) %*% rep(1, len=nrow(Y)))
    rc = outer(rowsum, colsum)
    Ybar = (Y - rc) / sqrt(rc)
    Q = qr(Ybar)
    if (Q$rank > 0) {
        temp = svd(Ybar)
        colnames(temp$u) = paste("CA", 1:length(temp$d), sep = "")
        rownames(temp$u) = dimnames(Y)[[1]]
        sweep(as.matrix(temp$u[,1:Q$rank, drop=FALSE]),
            1, 1/sqrt(rowsum), "*")
    } else stop("Null rank")
}



if (FALSE) {
    scores.uqo <- function (x, type = c("sites", "species"), ...) {
        if (mode(type) != "character" && mode(type) != "name")
            type = as.character(substitute(type))
        type = match.arg(type, c("sites", "species"))[1]
    
        switch(type,
            sites = if (any(slotNames(x)=="lv")) x@lv else Coef(x)@lv,
            species = if (any(slotNames(x)=="Optimum")) x@Optimum else Coef(x)@Optimum
        )
    }

    setMethod("scores", "uqo", function(x, ...) scores.uqo(x, ...))
}

jitteruqo = function(mat) {
    mat * ifelse(runif(length(mat)) < 0.5, -1, 1)
}

setMethod("Opt", "uqo", function(object, ...) Opt.qrrvglm(object, ...))
setMethod("Max", "uqo", function(object, ...) Max.qrrvglm(object, ...))
setMethod("lv",  "uqo", function(object, ...)  lv.qrrvglm(object, ...))



if (!isGeneric("calibrate"))
    setGeneric("calibrate", function(object, ...) standardGeneric("calibrate"))
setMethod("calibrate", "uqo", function(object, ...)
          calibrate.qrrvglm(object, ...))



summary.uqo = function(object, ...) {
    answer = Coef(object, ...)
    class(answer) = "summary.uqo"
    answer@call = object@call
    answer@misc = object@misc
    answer
}

printsummary.uqo = function(x, ...) {

    cat("\nCall:\n")
    dput(x@call)

    printCoef.qrrvglm(x, ...)

    cat("\nNumber of responses/species: ", x@NOS, "\n")

    if (length(x@misc$dispersion) == 1) 
    cat("\nDispersion parameter(s): ", x@misc$dispersion, "\n")
    invisible(x)
}

setClass("summary.uqo", representation("Coef.uqo",
         "misc" = "list",
         "call" = "call"))

setMethod("summary", "uqo", function(object, ...)
    summary.uqo(object, ...))

setMethod("print", "summary.uqo",
          function(x, ...)
          invisible(printsummary.uqo(x, ...)))

setMethod("show", "summary.uqo",
          function(object)
          invisible(printsummary.uqo(object)))



Tol.uqo = function(object, ...) {
    Coef(object, ...)@Tolerance
}

Tol.Coef.uqo = function(object, ...) {
    if (length(list(...))) warning("Too late! Ignoring the extra arguments")
    Coef(object, ...)@Tolerance
}

if (!isGeneric("Tol"))
    setGeneric("Tol", function(object, ...) standardGeneric("Tol"))
setMethod("Tol", "uqo", function(object, ...) Tol.uqo(object, ...))
setMethod("Tol", "Coef.uqo", function(object, ...) Tol.Coef.uqo(object, ...))






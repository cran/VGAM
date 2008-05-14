# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.







cao.fit <- function(x, y, w=rep(1, length(x[, 1])),
    etastart=NULL, mustart=NULL, coefstart=NULL,
    offset=0, family,
    control=cao.control(...), criterion="coefficients",
    qr.arg=FALSE, constraints=NULL, extra=NULL,
    Terms=Terms, function.name="cao", ...)
{
    specialCM = NULL
    post = list()
    check.rank = TRUE # 
    nonparametric <- TRUE
    optim.maxit <- control$optim.maxit
    backchat <- FALSE # control$backchat && !control$Quadratic  # rrr;
    save.weight <- control$save.weight
    trace <- control$trace
    minimize.criterion <- control$min.criterion

    n <- dim(x)[1]


    copyxbig <- FALSE    # May be overwritten in @initialize

    xbig.save <- NULL

    intercept.only <- ncol(x) == 1 && dimnames(x)[[2]] == "(Intercept)"
    y.names <- predictors.names <- NULL    # May be overwritten in @initialize

 
    n.save <- n 


    Rank <- control$Rank
    rrcontrol <- control  #

    if(length(family@initialize))
        eval(family@initialize)       # Initialize mu and M (and optionally w)
    n <- n.save 

    modelno = switch(family@vfamily[1], "poissonff"=2,
              "binomialff"=1, "quasipoissonff"=0, "quasibinomialff"=0,
              "negbinomial"=3,
              "gamma2"=5, "gaussianff"=8,
              0)  # stop("can't fit this model using fast algorithm")
    if(!modelno) stop("the family function does not work with cao()")
    if(modelno == 1) modelno = get("modelno", envir = VGAMenv)

    eval(rrr.init.expression)

    if(length(etastart)) {
        eta <- etastart
        mu <- if(length(mustart)) mustart else family@inverse(eta, extra)
    } else {
        if(length(mustart))
            mu <- mustart
        eta <- family@link(mu, extra)
    }

    M <- if(is.matrix(eta)) ncol(eta) else 1



    if(length(family@constraints))
        eval(family@constraints)


    special.matrix = matrix(-34956.125, M, M)    # An unlikely used matrix 
    just.testing <- cm.vgam(special.matrix, x, rrcontrol$Norrr, constraints)
    findex = trivial.constraints(just.testing, special.matrix)
    tc1 = trivial.constraints(constraints)


    if(all(findex == 1))
        stop("No covariates to form latent variables from.")
    colx1.index = names.colx1.index = NULL
    dx2 = dimnames(x)[[2]]
    if(sum(findex)) {
        asx = attr(x, "assign")
        for(ii in names(findex))
            if(findex[ii]) {
                names.colx1.index = c(names.colx1.index, dx2[asx[[ii]]])
                colx1.index = c(colx1.index, asx[[ii]])
        }
        names(colx1.index) = names.colx1.index
    }
    rrcontrol$colx1.index=control$colx1.index=colx1.index #Save it on the object
    colx2.index = 1:ncol(x)
    names(colx2.index) = dx2
    colx2.index = colx2.index[-colx1.index]
    p1 = length(colx1.index); p2 = length(colx2.index)
    rrcontrol$colx2.index=control$colx2.index=colx2.index #Save it on the object



    Cmat = if(length(rrcontrol$Cinit)) matrix(rrcontrol$Cinit,p2,Rank) else {
                if(!rrcontrol$Use.Init.Poisson.QO) {
                    matrix(rnorm(p2 * Rank, sd=rrcontrol$SD.Cinit), p2, Rank)
                } else {
                    .Init.Poisson.QO(ymat=as.matrix(y),
                                     X1=x[,colx1.index,drop=FALSE],
                                     X2=x[,colx2.index,drop=FALSE],
                                     Rank=rrcontrol$Rank, trace=rrcontrol$trace,
                                     max.ncol.etamat = rrcontrol$Etamat.colmax,
                                     Crow1positive=rrcontrol$Crow1positive,
                      constwt= any(family@vfamily[1] ==
                               c("negbinomial","gamma2","gaussianff")),
                      takelog= any(family@vfamily[1] != c("gaussianff")))
                }
            }


    rrcontrol$Cinit = control$Cinit = Cmat   # Good for valt()

    Blist <- process.constraints(constraints, x, M, specialCM=specialCM)

    nice31 = checkCMCO(Blist, control=control, modelno=modelno)
    if(nice31 != 1) stop("not nice")

    ncolBlist <- unlist(lapply(Blist, ncol))
    lv.mat = x[,colx2.index,drop=FALSE] %*% Cmat 


    rmfromVGAMenv(c("etamat", "beta"), prefix=".VGAM.CAO.")

    Nice21 = length(names.colx1.index)==1 && names.colx1.index == "(Intercept)"
    if(!Nice21) stop("Norrr = ~ 1 is supported only, without constraints")
    NOS = ifelse(modelno==3 || modelno==5, M/2, M)
    p1star. = if(Nice21) ifelse(modelno==3 || modelno==5,2,1) else M
    p2star. = if(Nice21) Rank else stop("not Nice21")
    pstar. = p1star. + p2star. 
    nstar = if(Nice21) ifelse(modelno==3 || modelno==5,n*2,n) else n*M
    lenbeta = pstar. * ifelse(Nice21, NOS, 1)

    othint = c(Rank,control$EqualTol, pstar. , dimw=1, inited=0, # w(,dimw) cols
               modelno, maxitl=control$maxitl, actnits=0, twice=0, p1star. ,
               p2star. , Nice21, lenbeta, controlITolerances=0, control$trace,
               p1, p2=p2, imethod=control$method.init, bchat=0)
    othdbl = c(small=control$SmallNo, fseps=control$epsilon,
               .Machine$double.eps,
               iKvector=rep(control$iKvector, len=NOS),
               iShape=rep(control$iShape, len=NOS),
               resss=0, bfeps=control$bf.epsilon, hstep=0.1)

    for(iter in 1:optim.maxit) {
        if(control$trace) {
            cat("\nIteration", iter, "\n")
            if(exists("flush.console"))
                flush.console()
        }
flush.console()

        conjgrad = optim(par=c(Cmat), fn=callcaof, 
                     gr=if(control$GradientFunction) calldcaof else NULL,
                     method="BFGS",
                     control=list(fnscale=1, trace=as.integer(control$trace),
                                  maxit=control$Maxit.optim, REPORT=10),
                     etamat=eta, xmat=x, ymat=y, # as.matrix(y), 
                     wvec=w, modelno=modelno,
                     Control=control,
                     Nice21=Nice21,
                     p1star. = p1star. , p2star. = p2star. ,
                     n=n, M=M, 
                     othint=othint, othdbl=othdbl,
                     alldump=FALSE)


        Cmat = matrix(conjgrad$par, p2, Rank) # old because of scale(cmatrix)

   #    Cmat <- Cmat %*% Ut  # Normalized

        if(converged <- (conjgrad$convergence == 0)) break
    }

    if(!converged) {
        if(maxitl > 1) {
            warning(paste("convergence not obtained in", maxitl, "iterations."))
        } else {
            warning(paste("convergence not obtained"))
        }
    } else {
    }
    Cmat = crow1C(Cmat, control$Crow1positive)   # Make sure the signs are right

flush.console()
    temp9 = 
    callcaof(cmatrix=Cmat,
             etamat=eta, xmat=x, ymat=y, wvec=w, modelno=modelno,
             Control=control,
             Nice21=Nice21,
             p1star. = p1star. , p2star. = p2star. ,
             n=n, M=M, 
             othint=othint, othdbl=othdbl,
             alldump=TRUE)
    if(!is.list(extra))
        extra = list()
    extra$Cmat = temp9$Cmat
    ynames = dimnames(y)[[2]]
    extra$df1.nl = temp9$df1.nl
    extra$spar1 = temp9$spar1
    names(extra$spar1) = ynames
    names(extra$df1.nl) = ynames
    if(Rank == 2) {
        extra$spar2 = temp9$spar2
        extra$df2.nl = temp9$df2.nl
        names(extra$spar2) = ynames
        names(extra$df2.nl) = ynames
    }

    mu = matrix(temp9$fitted, n, NOS, byrow=TRUE)











    dn <- labels(x)
    yn <- dn[[1]]


    if(is.matrix(mu)) {
          if(length(dimnames(y)[[2]])) {
              y.names <- dimnames(y)[[2]]
          }
          if(length(dimnames(mu)[[2]])) {
              y.names <- dimnames(mu)[[2]]
          }
          dimnames(mu) <- list(yn, y.names)
    } else {
        names(mu) <- names(fv)
    }


    fit <- list(
                fitted.values=mu,
                Cmatrix = Cmat,
                terms=Terms) # terms: This used to be done in vglm() 




    misc <- list(
        criterion = criterion,
        predictors.names = predictors.names,
        M = M,
        n = n,
        nonparametric = nonparametric,
        p = ncol(x),
        ynames = ynames)

    crit.list <- list()
    crit.list$deviance = temp9$deviance


                                    


    if(w[1] != 1 || any(w != w[1]))
        fit$prior.weights <- w

    if(length(family@last))
        eval(family@last)

    structure(c(fit, 
        temp9,
        list(
        contrasts=attr(x, "contrasts"),
        control=control,
        crit.list=crit.list,
        extra=extra,
        family=family,
        iter=iter,
        misc=misc,
        post = post,
        x=x,
        y=y)),
        vclass=family@vfamily)
}





cao.control = function(Rank=1,
          all.knots = FALSE,
          criterion="deviance",
          Cinit=NULL,
          Crow1positive=TRUE,
          epsilon = 1.0e-05,
          Etamat.colmax = 10,
          GradientFunction=FALSE,  # For now 24/12/04
          iKvector = 0.1,
          iShape = 0.1,
          Norrr = ~ 1,
          SmallNo = 5.0e-13,
          Use.Init.Poisson.QO=TRUE,

          Bestof = if(length(Cinit)) 1 else 10,
          maxitl = 40,
          method.init = 1,
          bf.epsilon = 1.0e-7,
          bf.maxit = 40,
          Maxit.optim = 250,
          optim.maxit = 20,
          SD.sitescores = 1.0,
          SD.Cinit = 0.02,
          trace = TRUE,
          df1.nl = 2.5, # About 1.5-2.5 gives the flexibility of a quadratic
          df2.nl = 2.5, # About 1.5-2.5 gives the flexibility of a quadratic
          spar1 = 0,    # 0 means df1.nl is used
          spar2 = 0,    # 0 means df2.nl is used
          ...)
{
    if(!is.Numeric(iShape, posit=TRUE)) stop("bad input for \"iShape\"")
    if(!is.Numeric(iKvector, posit=TRUE)) stop("bad input for \"iKvector\"")
    if(!is.Numeric(method.init, posit=TRUE, allow=1, integer=TRUE))
        stop("bad input for \"method.init\"")
    if(criterion != "deviance") stop("\"criterion\" must be \"deviance\"")
    if(GradientFunction) stop("14/1/05; GradientFunction=TRUE not working yet")
    se.fit = as.logical(FALSE)
    if(se.fit) stop("se.fit = FALSE handled only")

    if(length(Cinit) && !is.Numeric(Cinit))
        stop("Bad input for \"Cinit\"")
    if(!is.Numeric(Bestof, allow=1, integ=TRUE, posit=TRUE))
        stop("Bad input for \"Bestof\"")
    if(!is.Numeric(maxitl, allow=1, integ=TRUE, posit=TRUE))
        stop("Bad input for \"maxitl\"")
    if(!is.Numeric(bf.epsilon, allow=1, posit=TRUE))
        stop("Bad input for \"bf.epsilon\"")
    if(!is.Numeric(bf.maxit, integ=TRUE, posit=TRUE, allow=1))
        stop("Bad input for \"bf.maxit\"")
    if(!is.Numeric(Etamat.colmax, posit=TRUE, allow=1) || Etamat.colmax < Rank)
        stop("bad input for \"Etamat.colmax\"")
    if(!is.Numeric(Maxit.optim, integ=TRUE, posit=TRUE, allow=1))
        stop("Bad input for \"Maxit.optim\"")
    if(!is.Numeric(optim.maxit, allow=1, integ=TRUE, posit=TRUE))
        stop("Bad input for \"optim.maxit\"")
    if(!is.Numeric(SD.sitescores, allow=1, posit=TRUE))
        stop("Bad input for \"SD.sitescores\"")
    if(!is.Numeric(SD.Cinit, allow=1, posit=TRUE))
        stop("Bad input for \"SD.Cinit\"")
    if(!is.Numeric(df1.nl) || any(df1.nl < 0))
        stop("Bad input for \"df1.nl\"")
    if(any(df1.nl >= 0 & df1.nl < 0.05)) {
        warning("df1.nl values between 0 and 0.05 converted to 0.05")
        df1.nl[df1.nl < 0.05] = 0.05
    }
    if(!is.Numeric(df2.nl) || any(df2.nl < 0))
        stop("Bad input for \"df2.nl\"")
    if(any(df2.nl >= 0 & df2.nl < 0.05)) {
        warning("df2.nl values between 0 and 0.05 converted to 0.05")
        df2.nl[df2.nl < 0.05] = 0.05
    }
    if(!is.Numeric(spar1) || any(spar1 < 0))
        stop("Bad input for \"spar1\"")
    if(!is.Numeric(spar2) || any(spar2 < 0))
        stop("Bad input for \"spar2\"")
    if(!is.Numeric(epsilon, posit=TRUE, allow=1))
        stop("Bad input for \"epsilon\"")

    if(!is.Numeric(SmallNo, posit=TRUE, allow=1))
        stop("Bad input for \"SmallNo\"")
    if((SmallNo < .Machine$double.eps) ||
       (SmallNo > .0001)) stop("SmallNo is out of range") 

    ans = list(
        Corner=FALSE, # A constant, not a control parameter; unneeded?
        EqualTolerances=FALSE, # A constant, not a control parameter; needed
        ITolerances=FALSE, # A constant, not a control parameter; unneeded?
        Quadratic=FALSE, # A constant, not a control parameter; unneeded?
           all.knots = as.logical(all.knots)[1],
           Bestof = Bestof,
           Cinit=Cinit,
           ConstrainedO = TRUE, # A constant, not a control parameter
           criterion=criterion,
           Crow1positive=as.logical(rep(Crow1positive, len=Rank)),
           epsilon = epsilon,
           Etamat.colmax = Etamat.colmax,
           FastAlgorithm = TRUE, # A constant, not a control parameter
           GradientFunction = as.logical(GradientFunction),
           maxitl = maxitl,
           bf.epsilon = bf.epsilon,
           bf.maxit = bf.maxit,
           method.init = method.init,
           Maxit.optim = Maxit.optim,
           optim.maxit = optim.maxit,
           Norrr=Norrr,
           Rank = Rank,
           SD.sitescores = SD.sitescores,
           SD.Cinit = SD.Cinit,
           se.fit = se.fit, # If TRUE, then would need storage for S QR fits
           SmallNo = SmallNo,
           trace = as.integer(trace),
           Use.Init.Poisson.QO=Use.Init.Poisson.QO,
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


create.cms <- function(Rank=1, M, MSratio=1, which, p1=1) {
    if(!is.Numeric(p1, allow=1, integ=TRUE, pos=TRUE)) stop("bad input for p1")
    Blist. = vector("list", p1+Rank)
    for(r in 1:(p1+Rank))
        Blist.[[r]] = diag( M )
    names(Blist.) = if(p1 == 1) c("(Intercept)", names(which)) else stop()
    if(MSratio == 2) {
        for(r in 1:Rank) 
            Blist.[[p1+r]] = eij(1, M)
    }
    Blist.
}


callcaof = function(cmatrix,
                    etamat, xmat, ymat, wvec, modelno, 
                    Control, Nice21=TRUE,
                    p1star. = if(any(modelno==c(3,5))) 2 else 1, p2star. =Rank,
                    n, M, 
                    othint, othdbl,
                    alldump=FALSE) {
if(exists("flush.console")) flush.console()

    control = Control
    Rank = control$Rank
    p1 = length(control$colx1.index)
    p2 = length(control$colx2.index)
    yn = dimnames(ymat)[[2]]
    if(length(yn) != ncol(ymat)) stop("the column names of ymat must be given")
    queue = qbig = Rank # 19/10/05; number of smooths per species
    NOS = if(any(modelno==c(3,5))) M/2 else M
    df1.nl = procVec(control$df1.nl, yn= yn , Def=control$DF1)
    spar1 = procVec(control$spar1, yn= yn , Def= control$SPAR1)
    df2.nl = procVec(control$df2.nl, yn= yn , Def=control$DF2)
    spar2 = procVec(control$spar2, yn= yn , Def= control$SPAR2)
    if(any(c(length(spar1),length(spar2),length(df1.nl),length(df2.nl)) != NOS))
        stop("wrong length in at least one of df1.nl, df2.nl, spar1, spar2")

    cmatrix = matrix(cmatrix, p2, Rank)  # crow1C() needs a matrix as input
        cmatrix = crow1C(cmatrix, crow=control$Crow1positive)
        numat = xmat[,control$colx2.index,drop=FALSE] %*% cmatrix
        evnu = eigen(var(numat))
        temp7 = if(Rank > 1) evnu$vector %*% diag(evnu$value^(-0.5)) else
                evnu$vector %*% evnu$value^(-0.5)
        cmatrix = cmatrix %*% temp7
        cmatrix = crow1C(cmatrix, crow=control$Crow1positive)
        numat = xmat[,control$colx2.index,drop=FALSE] %*% cmatrix


    dim(numat) = c(n, Rank)
    mynames5 = if(Rank==1) "lv" else paste("lv", 1:Rank, sep="")
    nu1mat = cbind("(Intercept)"=1, lv=numat)
    dimnames(nu1mat) = list(dimnames(xmat)[[1]], c("(Intercept)", mynames5))

    temp.smooth.frame = vector("list", p1+Rank)  # A temporary makeshift frame
    names(temp.smooth.frame) = c(names(control$colx1.index), mynames5)
    for(uu in 1:(p1+Rank)) {
        temp.smooth.frame[[uu]] = nu1mat[,uu]
    }
    temp.smooth.frame = data.frame(temp.smooth.frame)
    for(uu in 1:Rank) {
        attr(temp.smooth.frame[,uu+p1], "spar") = 0  # this value unused
        attr(temp.smooth.frame[,uu+p1], "df") = 4    # this value unused
    }

    pstar.  = p1star.  + p2star.   # = Mdot + Rank
    nstar = if(Nice21) ifelse(modelno==3 || modelno==5,n*2,n) else n*M
    lenbeta = pstar. * ifelse(Nice21, NOS, 1) # Holds the linear coeffs

    inited = if(is.R()) {
        if(exists(".VGAM.CAO.etamat", envir=VGAMenv)) 1 else 0
    } else 0
    usethiseta = if(inited==1)
        getfromVGAMenv("etamat", prefix = ".VGAM.CAO.") else t(etamat)

    if(any(is.na(usethiseta))) {
        usethiseta = t(etamat)  # So that dim(usethiseta)==c(M,n)
        rmfromVGAMenv("etamat", prefix=".VGAM.CAO.")
    }

    usethisbeta = if(inited==2)
        getfromVGAMenv("beta", prefix = ".VGAM.CAO.") else double(lenbeta)
    othint[5] = inited   # Refine initialization within FORTRAN
    pstar = NOS * pstar. 
    bnumat=if(Nice21) matrix(0,nstar,pstar.) else stop("code not written here")

    M. = MSratio = M / NOS     # 1 or 2 usually
    which = p1 + (1:Rank) # These columns are smoothed
    nwhich = names(which) = mynames5

    origBlist = Blist. = create.cms(Rank=Rank, M=M., MSratio=MSratio, 
                                    which=which, p1=p1) # For 1 species only
    ncolBlist. <- unlist(lapply(Blist. , ncol))
    smooth.frame = s.vam(x=nu1mat, z=NULL, wz=NULL, s=NULL,
                         which=which,
                         smooth.frame=temp.smooth.frame,
                         bf.maxit=control$bf.maxit,
                         bf.epsilon=control$bf.epsilon,
                         trace=FALSE, se.fit=control$se.fit,
                         xbig.save=bnumat, Blist=Blist. ,
                         ncolBlist=ncolBlist. ,
                         M= M. , qbig=NULL, U=NULL, # NULL implies not needed
                         backchat=FALSE,
                         all.knots=control$all.knots, nk=NULL,
                         sf.only=TRUE)

    ldk <- 3 * max(ncolBlist.[nwhich]) + 1   # 11/7/02

    dimw. = M.   # Smoothing one spp. at a time
    dimu. = M.
    wz. = matrix(0, n, dimw. )
    U. = matrix(0, dimu. , n)
    if(names(Blist.)[1] != "(Intercept)") stop("something wrong here")
    Blist.[[1]] <- NULL

    trivc = rep(2 - M. , len=queue)   # All of queue smooths are basic smooths
    ncbvec <- ncolBlist.[nwhich]
    ncolb <- max(ncbvec)
    pmax.mwk <- pmax(ncbvec*(ncbvec+1)/2, dimw. )
    size.twk <- max((4+4*smooth.frame$nef)*ncbvec + dimu. * smooth.frame$nef)
    size.twk <- max(size.twk, M*smooth.frame$n)

    qbig. = NOS * qbig    # == NOS * Rank; holds all the smooths
    if(!all.equal(as.vector(ncbvec), rep(1, len=queue)))
        stop("ncbvec not right---should be a queue-vector of ones")
    pbig = pstar. #
    backchat = FALSE


    npetc = c(n=nrow(nu1mat), p. =ncol(nu1mat), q=Rank, # q=length(which),
                  se.fit=control$se.fit, 0,
        control$bf.maxit, qrank=0, M= M. , nbig=nstar, pbig=pbig,
        qbig=qbig, dimw= dimw. , dimu= dimu. , ier=0, ldk=ldk)





    if(Rank == 2) {
        spardf = (c(spar1,1+df1.nl,spar2,1+df2.nl))[interleave.VGAM(4*NOS,M=2)]
    } else {
        spardf = c(spar1, 1.0+df1.nl)
    }

    ans1 <- dotFortran(name = "vcao6f",
       numat=as.double(numat),
           ymat=as.double(ymat), wvec=as.double(wvec),
       etamat=as.double(usethiseta),
           fv=double(NOS*n), z=double(n*M), wz=double(n*M),
           U=double(M*n), # bnumat=as.double(bnumat),
       qr=double(nstar*pstar.), qraux=double(pstar.), qpivot=integer(pstar.),
       n=as.integer(n), M=as.integer(M), NOS=as.integer(NOS),
       nstar=as.integer(nstar), dimu=as.integer( M ), # for U, not U. 
           errcode=integer(1), othint=as.integer(othint),
       deviance=double(1), beta=as.double(usethisbeta),
       twk=double(if(Nice21) nstar*3 else M*n*2), wkmm=double(M*(M+pstar)),
           othdbl=as.double(othdbl),
            npetc = as.integer(npetc), M. = as.integer( M. ),
            spardf = as.double(spardf),
        match=as.integer(smooth.frame$o), as.integer(smooth.frame$nef), 
            which=as.integer(which),
            etal = double( M. * n ),
            smomat = as.double(matrix(0, n, qbig. )),
            s0 = double((2* M. )*(2* M. )*2),
            U. = as.double( U. ), etapmat = as.double( U. ),
                nu1mat=as.double(nu1mat),
            blist=as.double(unlist( Blist. )), as.integer(ncbvec), 
            smap=as.integer(1:(Rank+1)), # 
            rcind = integer( M. *( M. +1)), trivc = as.integer(trivc),
        work1 = double(3*qbig + (9+2*4+max(smooth.frame$nknots))*
                     max(smooth.frame$nknots)),
            wk2 = double(n* M. *3),
           wwkmm = double( M. * M. *16 + M. * pbig),
            work3 = double(max(max(2 * smooth.frame$nef * ncbvec^2),
                           max(smooth.frame$nknots * ncbvec * (4*ncbvec+1)))),
        sgdub = double(max(smooth.frame$nknots) * max(4,ncolb)),
            bmb = double( M. * M. ),
            lev = double(NOS * max(smooth.frame$nef * ncbvec)),
        mwk = double(max(smooth.frame$nef * (1 + 2* M. + pmax.mwk)) ),
           ttwk = double(size.twk),
        bcoefficients = double(NOS * sum(smooth.frame$nknots*ncbvec)),
            knots = as.double(unlist(smooth.frame$knots)),
        bindex = as.integer(smooth.frame$bindex),
            nknots = as.integer(smooth.frame$nknots),
            itwk = integer(2 * M. ),
            kindex = as.integer(smooth.frame$kindex))
if(exists("flush.console")) flush.console()


    if(ans1$errcode == 0) {
        assign2VGAMenv(c("etamat", "beta"), ans1, prefix=".VGAM.CAO.")
        if(is.R()) {
            assign(".VGAM.CAO.cmatrix", matrix(cmatrix,p2,Rank), envir=VGAMenv)
        } else {
            .VGAM.CAO.cmatrix <<- matrix(cmatrix,p2,Rank) # matrix reqd for R=2
        }

    } else {
        cat("warning in callcaof: error code =", ans1$errcode, "\n")
        cat("warning in callcaof: npetc[14] =", ans1$npetc[14], "\n")
        if(exists("flush.console")) flush.console()
        rmfromVGAMenv(c("etamat", "beta"), prefix=".VGAM.CAO.")
    }

    returnans = if(alldump) {
        bindex = ans1$bindex
        ncolBlist = ncbvec
        Bspline2 <- vector("list", NOS)
        names(Bspline2) <- dimnames(ymat)[[2]]
        Bspline <- vector("list", length(nwhich))
        names(Bspline) <- nwhich
        ind9 = 0   # moving index
        for(sppno in 1:NOS) {
            for(ii in 1:length(nwhich)) {
                ind7 = (smooth.frame$bindex[ii]):(smooth.frame$bindex[ii+1]-1)
                ans = ans1$bcoeff[ind9+ind7]
                ans = matrix(ans, ncol=ncolBlist[nwhich[ii]])
                Bspline[[ii]] = new(Class="vsmooth.spline.fit",
                    "Bcoefficients" = ans,
                    "xmax"          = smooth.frame$xmax[ii],
                    "xmin"          = smooth.frame$xmin[ii],
                    "knots"         = as.vector(smooth.frame$knots[[ii]]))
            }
            ind9 = ind9 + smooth.frame$bindex[length(nwhich)+1]-1
            Bspline2[[sppno]] = Bspline
        }

        qrank = npetc[7]  # Assume all species have the same qrank value
        dim(ans1$etamat) = c(M,n)    # was c(n,M) prior to 22/8/06
        if(Rank == 2) {
             spardf = array(ans1$spardf, c(Rank,NOS,2))
             df1.nl = spardf[1,,2] - 1
             df2.nl = spardf[2,,2] - 1
             spar1 = spardf[1,,1]
             spar2 = spardf[2,,1]
        } else {
             df1.nl = ans1$spardf[  NOS+(1:NOS)] - 1
             spar1 = ans1$spardf[1:NOS]
        }
        list(deviance=ans1$deviance, 
             bcoefficients = ans1$bcoefficients,
             bindex = ans1$bindex,
             Bspline = Bspline2,
             Cmat = matrix(cmatrix, p2, Rank, dimnames=list(
                           names(control$colx2.index), mynames5)),
             coefficients = ans1$beta,
             df1.nl = df1.nl,
             df2.nl = if(Rank == 2) df2.nl else NULL,
             df.residual = n*M - qrank - sum(ans1$df - 1),
             fitted = ans1$fv,  # NOS x n
             kindex = ans1$kindex,
             predictors = matrix(ans1$etamat, n, M, byrow=TRUE),
             wresiduals = ans1$z - t(ans1$etamat),   # n x M
             spar1=spar1,
             spar2=if(Rank == 2) spar2 else NULL)
    } else
        ans1$deviance
    if(exists("flush.console")) flush.console()
    returnans
}



calldcaof = function(cmatrix,
                     etamat, xmat, ymat, wvec, modelno, 
                     Control, Nice21=TRUE,
                     p1star. = if(any(modelno==c(3,5))) 2 else 1, p2star. =Rank,
                     n, M, 
                     othint, othdbl,
                     alldump=FALSE) {


    if(alldump) stop("really used?")
if(exists("flush.console")) flush.console()

    if(!Nice21) stop("Nice21 must be TRUE")
    control = Control
    Rank = control$Rank
    p2 = length(control$colx2.index)
    yn = dimnames(ymat)[[2]]
    if(!length( yn )) yn = paste("Y", 1:ncol(ymat), sep="")


    cmatrix = scale(cmatrix)

    xmat2 <- xmat[,control$colx2.index,drop=FALSE]   #ccc
    numat <- xmat2 %*% matrix(cmatrix, p2, Rank)
    dim(numat) <- c(nrow(xmat), Rank)
    temp.smooth.frame = vector("list", 1+Rank)  # A temporary makeshift frame
    mynames5 = if(Rank==1) "lv" else paste("lv",1:Rank,sep="")
    names(temp.smooth.frame) = c("(Intercept)", mynames5)
    temp.smooth.frame[[1]] = rep(1, len=n)
    for(uu in 1:Rank) {
        temp.smooth.frame[[uu+1]] = numat[,uu]
    }
    temp.smooth.frame = data.frame(temp.smooth.frame)
    for(uu in 1:Rank) {
        attr(temp.smooth.frame[,uu+1], "spar") = 0 # any old value
        attr(temp.smooth.frame[,uu+1], "df") = 4 # any old value
    }
    pstar.  = p1star.  + p2star. 
    nstar = if(Nice21) ifelse(modelno==3 || modelno==5,n*2,n) else n*M
    NOS = ifelse(modelno == 3 || modelno==5, M/2, M)
    lenbeta = pstar. * ifelse(Nice21, NOS, 1)

    if(TRUE) {
        inited = if(is.R()) {
            if(exists(".VGAM.CAO.etamat", envir = VGAMenv)) 1 else 0
        } else 0
        usethiseta = if(inited==1) {if(is.R()) get(".VGAM.CAO.etamat",
            envir = VGAMenv) else .VGAM.CAO.etamat} else t(etamat)
    }
    usethisbeta = if(inited==2) {if(is.R()) get(".VGAM.CAO.beta",
        envir = VGAMenv) else .VGAM.CAO.beta} else double(lenbeta)





 pstar = NOS * pstar. 
    bnumat = if(Nice21) matrix(0,nstar,pstar) else stop("need Nice21")

    M. = MSratio = M / NOS     # 1 or 2 usually
    which = p1 + (1:Rank)   # The first 1 is the intercept term
    nwhich = names(which) = mynames5

    origBlist = Blist. = create.cms(Rank=Rank, M=M., MSratio=MSratio,
                                    which=which,p1=p1) # For 1 species
    ncolBlist. <- unlist(lapply(Blist. , ncol))
    nu1mat = cbind("(Intercept)"=1, lv=numat)
    dimnames(nu1mat) = list(dimnames(xmat)[[1]], c("(Intercept)","lv"))

    smooth.frame = s.vam(x=nu1mat, z=NULL, wz=NULL, s=NULL,
                         which=which,
                         smooth.frame=temp.smooth.frame,
                         bf.maxit=control$bf.maxit,
                         bf.epsilon=control$bf.epsilon,
                         trace=FALSE, se.fit=control$se.fit,
                         xbig.save=bnumat, Blist=Blist.,
                         ncolBlist=ncolBlist. ,
                         M= M. , qbig=NULL, U=U, # NULL value implies not needed
                         backchat=FALSE,
                         all.knots=control$all.knots, nk=NULL,
                         sf.only=TRUE)

    ldk <- 4 * max(ncolBlist.[nwhich])   # was M;     # Prior to 11/7/02
    ldk <- 3 * max(ncolBlist.[nwhich]) + 1   # 11/7/02



    wz. = matrix(0, n, M. )  # not sure
    U. = matrix(0, M. , n)
    dimw. = if(is.matrix( wz. )) ncol( wz. ) else 1
    dimu. <- if(is.matrix( U. )) nrow( U. ) else 1
    Blist.[[1]] <- NULL
    trivc = rep(2 - M. , len=queue)   # All of queue smooths are basic smooths
    ncbvec <- ncolBlist.[nwhich]
    ncolb <- max(ncbvec)
    pmax.mwk <- rep( dimw. , length(trivc))
    pmax.mwk <- pmax(ncbvec*(ncbvec+1)/2, dimw. )
    size.twk <- max((4+4*smooth.frame$nef)*ncbvec + dimu. *smooth.frame$nef)
    size.twk <- max(size.twk, M*smooth.frame$n)


    qbig. = NOS * qbig    # == NOS * Rank
    pbig = pstar. # Not sure
    if(FALSE) {
        df1.nl = rep(control$df1.nl, len=NOS)  # This is used
        spar1 = rep(control$spar1, len=NOS)   # This is used
    } else {
        # This is used
        df1.nl = procVec(control$df1.nl, yn= yn , Def=control$DF1)
        spar1 = procVec(control$spar1, yn= yn , Def= control$SPAR1)
    }
    backchat = FALSE


    npetc = c(n=n, p=1+Rank, length(which), se.fit=control$se.fit, 0,
        maxitl=control$maxitl, qrank=0, M= M. , n.M = n* M. ,
            pbig=sum( ncolBlist.),
        qbig=qbig, dimw= dimw. , dimu= dimu. , ier=0, ldk=ldk)

if(exists("flush.console")) flush.console()

    ans1 <- 
  dotFortran(name = if(Nice21) "vdcaof" else stop("need Nice21"),
       numat=as.double(numat),
           as.double(ymat), as.double(wvec),
       etamat=as.double(usethiseta),
           fv=double(NOS*n), z=double(n*M), wz=double(n*M),
           U=double(M*n), # bnumat=as.double(bnumat),
       qr=double(nstar*pstar.), qraux=double(pstar.), qpivot=integer(pstar.),
       as.integer(n), as.integer(M), NOS=as.integer(NOS),
       as.integer(nstar), dimu=as.integer(M),
           errcode=integer(1), othint=as.integer(othint),
       deviance=double(1), beta=as.double(usethisbeta),
       twk=double(if(Nice21) nstar*3 else M*n*2), wkmm=double(M*(M+pstar)),
           othdbl=as.double(othdbl),
       xmat2=as.double(xmat2), onumat=as.double(numat), cmat=as.double(cmatrix),
       p2=as.integer(p2), deriv=double(p2*Rank),
       betasave=double(lenbeta), 
            npetc = as.integer(npetc), M. = as.integer( M. ),
            spardf = as.double(c(spar1, 1.0+df1.nl, spar2, 1.0+df2.nl)),
        match=as.integer(smooth.frame$o), as.integer(smooth.frame$nef), 
            as.integer(which),
        etal = double( M. * n ),
            smomat = as.double(matrix(0, n, qbig. )),
            s0 = double((2* M. )*(2* M. )*2),
            U. = as.double( U. ), etapmat = as.double( U. ),
                nu1mat=as.double(nu1mat),
            as.double(unlist( Blist. )),
        as.integer(ncbvec), smap=as.integer(1:(Rank+1)),
            rcind = integer( M. *( M. +1)), trivc = as.integer(trivc),
        work1 = double(3*qbig + (9+2*4+max(smooth.frame$nknots))*
                     max(smooth.frame$nknots)),
            wk2 = double(n* M. *3),
           wwkmm = double( M. * M. *16 + M. *pbig),
            work3 = double(max(max(2 * smooth.frame$nef * ncbvec^2),
                           max(smooth.frame$nknots * ncbvec * (4*ncbvec+1)))),
        sgdub = double(max(smooth.frame$nknots) * max(4,ncolb)),
            bmb = double( M. * M. ),
            lev = double(NOS * max(smooth.frame$nef * ncbvec)),
        mwk = double(max(smooth.frame$nef * (1 + 2* M. + pmax.mwk)) ),
           ttwk = double(size.twk),
        bcoefficients = double(NOS * sum(smooth.frame$nknots*ncbvec)),
            knots = as.double(unlist(smooth.frame$knots)),
        bindex = as.integer(smooth.frame$bindex),
            nknots = as.integer(smooth.frame$nknots),
            itwk = integer(2* M. ),
            kindex = as.integer(smooth.frame$kindex))
if(exists("flush.console")) flush.console()

           if(is.R()) {
               assign(".VGAM.CAO.etamat", ans1$etamat, envir = VGAMenv)
           assign(".VGAM.CAO.z",ans1$z,envir=VGAMenv)# z; minus any offset
               assign(".VGAM.CAO.U", ans1$U, envir=VGAMenv)  # U
           } else {
               .VGAM.CAO.etamat <<- ans1$etamat
               .VGAM.CAO.z <<- ans1$z
               .VGAM.CAO.U <<- ans1$U
           }
       if(ans1$errcode == 0) {
       } else {
           cat("warning in calldcaof: error code =", ans1$errcode, "\n")
            if(exists("flush.console"))
                flush.console()
       }

    returnans = if(alldump) {
        bindex = ans1$bindex
        ncolBlist = ncbvec
        Bspline2 <- vector("list", NOS)
        names(Bspline2) <- dimnames(ymat)[[2]]
        Bspline <- vector("list", length(nwhich))
        names(Bspline) <- nwhich
        ind9 = 0   # moving index
        for(j in 1:NOS) {
            for(i in 1:length(nwhich)) {
                ind9 = ind9[length(ind9)] + (bindex[i]):(bindex[i+1]-1)
                ans = ans1$bcoeff[ind9]
                ans = matrix(ans, ncol=ncolBlist[nwhich[i]])
                Bspline[[i]] = new(Class="vsmooth.spline.fit",
                    "Bcoefficients" = ans,
                    "xmax"          = smooth.frame$xmax[i],
                    "xmin"          = smooth.frame$xmin[i],
                    "knots"         = as.vector(smooth.frame$knots[[i]]))
            }
            Bspline2[[j]] = Bspline
        }

        qrank = npetc[7]  # Assume all species have the same qrank value
        dim(ans1$etamat) = c(M,n)   # bug: was c(n,M) prior to 22/8/06
        list(deviance=ans1$deviance, 
             bcoefficients = ans1$bcoefficients,
             bindex = ans1$bindex,
             Bspline = Bspline2,
             Cmat=matrix(cmatrix, p2, Rank, dimnames=list(
                         names(control$colx2.index), mynames5)),
             coefficients=ans1$beta,
             df1.nl=ans1$spardf[  NOS+(1:NOS)] - 1,
             df2.nl=ans1$spardf[3*NOS+(1:NOS)] - 1,
             df.residual = n*M - qrank - sum(ans1$df - 1),
             fitted=ans1$fv,
             kindex = ans1$kindex,
             predictors=matrix(ans1$etamat, n, M, byrow=TRUE),
             wresiduals = ans1$z - t(ans1$etamat),   # n x M
             spar1=ans1$spardf[1:NOS],
             spar2=ans1$spardf[2*NOS+(1:NOS)])
    } else
        ans1$deriv
    if(exists("flush.console")) flush.console()
    returnans 
}






setClass(Class="Coef.cao", representation(
      "Bspline"      = "list",
      "C"            = "matrix",
      "Constrained"  = "logical",
      "df1.nl"       = "numeric",
      "df2.nl"       = "numeric",
      "dispersion"   = "numeric",
      "eta2"         = "matrix",
      "lv"           = "matrix",
      "lvOrder"      = "matrix",
      "M"            = "numeric",
      "Maximum"      = "numeric",
      "NOS"          = "numeric",
      "Optimum"      = "matrix",
      "OptimumOrder" = "matrix",
      "Rank"         = "numeric",
      "spar1"        = "numeric",
      "spar2"        = "numeric"))


Coef.cao = function(object,
    epsOptimum = 0.00001, # determines how accurately Optimum is estimated
    gridlen = 40,  # Number of points on the grid (one level at a time)
    maxgriditer = 10,    # Maximum number of iterations allowed for grid search
    smallno = 0.05,
    ...) {

    if(!is.Numeric(epsOptimum, posit=TRUE, allow=1))
        stop("bad input for argument 'epsOptimum'")
    if(!is.Numeric(gridlen, posit=TRUE, integer=TRUE) || gridlen < 5)
        stop("bad input for argument 'gridlen'")
    if(!is.Numeric(maxgriditer, posit=TRUE, allow=1, int=TRUE) || maxgriditer<3)
        stop("bad input for argument 'maxgriditer'")
    if(!is.logical(ConstrainedO <- object@control$ConstrainedO))
        stop("can't determine whether the model is constrained or not")
    if(!is.Numeric(smallno, posit=TRUE, allow=1) ||
       smallno > 0.5 || smallno < 0.0001)
        stop("bad input for argument 'smallno'")
    ocontrol = object@control
    if((Rank <- ocontrol$Rank) > 2) stop("Rank must be 1 or 2") 
    gridlen = rep(gridlen, length=Rank)
    M = if(any(slotNames(object) == "predictors") &&
           is.matrix(object@predictors)) ncol(object@predictors) else
           object@misc$M
    NOS = if(length(object@y)) ncol(object@y) else M
    MSratio = M / NOS # 1 or 2; First value is g(mean) = quadratic form in lv
    nice21 = (length(ocontrol$colx1.index) == 1) &&
             (names(ocontrol$colx1.index) == "(Intercept)")
    if(!nice21) stop("Can only handle Norrr = ~ 1")

    p1 = length(ocontrol$colx1.index)
    p2 = length(ocontrol$colx2.index)
    modelno = object@control$modelno  # 1,2,3,... or 0
    ynames = object@misc$ynames
    if(!length(ynames)) ynames = object@misc$predictors.names
    if(!length(ynames)) ynames = object@misc$ynames
    if(!length(ynames)) ynames = paste("Y", 1:NOS, sep="")
    lp.names = object@misc$predictors.names
    if(!length(lp.names)) lp.names = NULL 

    lv.names = if(Rank==1) "lv" else paste("lv", 1:Rank, sep="")
    Cmat = object@extra$Cmat   # p2 x Rank (provided maxitl > 1)
    if(ConstrainedO)
        dimnames(Cmat) = list(names(ocontrol$colx2.index), lv.names)
    lv.mat = if(ConstrainedO) {
        object@x[,ocontrol$colx2.index,drop=FALSE] %*% Cmat 
    } else {
        object@lv
    }

    optimum = matrix(as.numeric(NA), Rank, NOS, dimnames=list(lv.names, ynames))
    extents = apply(lv.mat, 2, range)  # 2 by R

    maximum = rep(as.numeric(NA), len=NOS)

    whichSpecies = 1:NOS  # Do it for all species
    if(Rank == 1) {
        gridd = cbind(seq(extents[1,1], extents[2,1], len=gridlen))
    } else {
        gridd = expand.grid(seq(extents[1,1], extents[2,1], len=gridlen[1]),
                            seq(extents[1,2], extents[2,2], len=gridlen[2]))
        eta2matrix = matrix(0, NOS, 1)
    }
    gridd.orig = gridd
    # if(Rank == 2) then this is for initial values
    for(sppno in 1:length(whichSpecies)) {
        gridd = gridd.orig 
        gridres1 = gridd[2,1] - gridd[1,1]
        gridres2 = if(Rank==2) gridd[2,2] - gridd[1,2] else 0
        griditer = 1

        thisSpecies = whichSpecies[sppno]
        indexSpecies = if(is.character(whichSpecies))
            match(whichSpecies[sppno], ynames) else whichSpecies[sppno]

        if(is.na(indexSpecies))
            stop("mismatch found in \"whichSpecies\"")

        while(griditer == 1 ||
              ((griditer <= maxgriditer) &&
              ((gridres1 > epsOptimum) || (gridres2 > epsOptimum)))) {
            temp = predictcao(object, grid=gridd, sppno=thisSpecies,
                              Rank=Rank, deriv=0, MSratio=MSratio)
            yvals = temp$yvals  # gridlen-vector
            xvals = temp$xvals  # gridlen x Rank; gridd
            if(length(temp$eta2)) eta2matrix[sppno,1] = temp$eta2

            nnn = length(yvals)
            index = (1:nnn)[yvals==max(yvals)]
            if(length(index)!=1) warning("could not find a single maximum")
            if(Rank == 2) {
                initvalue = rep(xvals[index,], length=Rank) # for optim()
                # Make sure initvalue is in the interior
                if(abs(initvalue[1] - extents[1,1]) < smallno)
                    initvalue[1] = extents[1,1] + smallno
                if(abs(initvalue[1] - extents[2,1]) < smallno)
                    initvalue[1] = extents[2,1] - smallno
                if(abs(initvalue[2] - extents[1,2]) < smallno)
                    initvalue[2] = extents[1,2] + smallno
                if(abs(initvalue[2] - extents[2,2]) < smallno)
                    initvalue[2] = extents[2,2] - smallno
                break
            }
            if(index == 1 || index == nnn) {
                maximum[sppno] = optimum[1,sppno] = NA
                gridres1 = epsOptimum + 1 # equivalent to a break
                break          # just in case
            } else {
                maximum[sppno] = yvals[index] # on the eta scale
                optimum[1,sppno] = xvals[index,1]
                gridd[,1] = seq(
                    max(extents[1,1], optimum[1,sppno]-gridres1),
                    min(extents[2,1], optimum[1,sppno]+gridres1),
                    len=gridlen)
                gridres1 = gridd[2,1] - gridd[1,1]
                griditer = griditer + 1
            }
        } # of while 

        if(Rank == 2) {
            # Rank = 2, so use optim(). The above was to get initial values.
            myfun = function(x, object, sppno, Rank=1, deriv=0, MSratio=1) {
                # x is a 2-vector
                x = matrix(x, 1, length(x))
                temp = predictcao(object, grid=x, sppno=sppno,
                                  Rank=Rank, deriv=deriv, MSratio=MSratio)
                temp$yval
            }
            answer = optim(initvalue, myfun, gr=NULL, method="L-BFGS-B",
                           lower=extents[1,], upper=extents[2,],
                           control=list(fnscale = -1),  # maximize!
                           object=object, sppno=sppno, Rank=Rank,
                           deriv=0, MSratio=MSratio)
            # Check to see if the soln is at the boundary. If not, assign it.
            for(rindex in 1:Rank)
                if(abs(answer$par[rindex] - extents[1,rindex]) > smallno &&
                   abs(answer$par[rindex] - extents[2,rindex]) > smallno) {
                    optimum[rindex,sppno] = answer$par[rindex]
                    maximum[sppno] = answer$value
                }
        } # end of Rank=2
    } # end of sppno 
    myetamat = rbind(maximum)
    if(MSratio == 2) myetamat = kronecker(myetamat, matrix(1:0, 1, 2))
    maximum = object@family@inverse(eta=myetamat, extra=object@extra)
    maximum = c(maximum)  # Convert from matrix to vector 
    names(maximum) = ynames

    ans = new(Class="Coef.cao",
              Bspline = object@Bspline,
              Constrained=ConstrainedO,
              df1.nl = object@extra$df1.nl,
              lv = lv.mat,
              lvOrder = lv.mat,
              Maximum = maximum,
              M = M,
              NOS = NOS, 
              Optimum=optimum, 
              OptimumOrder=optimum, 
              Rank = Rank,
              spar1 = object@extra$spar1)
    if(ConstrainedO) {ans@C = Cmat} else {Cmat = NULL}
    if(Rank == 2) {
        dimnames(eta2matrix) = list(
            object@misc$predictors.names[c(FALSE,TRUE)], " ")
        ans@eta2 = eta2matrix
        ans@df2.nl = object@extra$df2.nl 
        ans@spar2  = object@extra$spar2
    }

    for(rindex in 1:Rank) {
        ans@OptimumOrder[rindex,] = order(ans@Optimum[rindex,])
        ans@lvOrder[,rindex] = order(ans@lv[,rindex])
    }

    if(length(object@misc$estimated.dispersion) &&
       object@misc$estimated.dispersion) {
        p = length(object@coefficients)
        n = object@misc$n
        M = object@misc$M
        NOS = if(length(object@y)) ncol(object@y) else M
        pstar = p + length(Cmat) # Adjustment 
        adjusted.dispersion = object@misc$dispersion * (n*M - p) /
                (n*M - pstar)
        ans@dispersion = adjusted.dispersion 
    }
    if(MSratio == 2) {
        lcoef = object@coefficients
        temp = lcoef[((1:NOS)-1)*(2+Rank)+2]
        names(temp) = object@misc$predictors.names[2*(1:NOS)]
        ans@dispersion = temp
    }
    dimnames(ans@Optimum) = list(lv.names, ynames)
    ans 
}


printCoef.cao = function(object, digits = max(2, options()$digits-2), ...) {
    Rank = object@Rank
    NOS = object@NOS
    M = object@M

    Maximum = if(length(object@Maximum)) cbind(Maximum=object@Maximum) else NULL
    optmat = cbind(t(object@Optimum))
    dimnames(optmat) = list(dimnames(optmat)[[1]],
        if(Rank > 1) paste("Optimum", dimnames(optmat)[[2]], sep=".")
        else "Optimum")

    if( object@Constrained ) {
        cat("\nC matrix (constrained/canonical coefficients)\n")
        print(object@C, digits=digits, ...)
    }
    cat("\nOptima and maxima\n")
    print(cbind(Optimum=optmat,
                Maximum), digits = max(1, digits-1))
    cat("\nNonlinear degrees of freedom\n")
    if(Rank == 1) {
        print(cbind(df1.nl = object@df1.nl), digits=max(2, digits-1), ...)
    } else {
        print(cbind(df1.nl = object@df1.nl,
                    df2.nl = object@df2.nl), digits=max(2, digits-1), ...)
    }
    invisible(object)
}





    setMethod("show", "Coef.cao", function(object)
        printCoef.cao(object))
    setMethod("print", "Coef.cao", function(x, ...)
        printCoef.cao(object=x, ...))


setMethod("coef", "cao", function(object, ...) Coef.cao(object, ...))
setMethod("coefficients", "cao", function(object, ...) Coef.cao(object, ...))
setMethod("Coef", "cao", function(object, ...) Coef.cao(object, ...))




lvplot.cao = function(object,
          add= FALSE, plot.it= TRUE, rugplot = TRUE, y = FALSE, 
          type=c("fitted.values", "predictors"),
          xlab=paste("Latent Variable", if(Rank==1) "" else " 1", sep=""),
          ylab=if(Rank==1) switch(type, predictors="Predictors", 
              fitted.values="Fitted values") else "Latent Variable 2",
          pcex=par()$cex, pcol=par()$col, pch=par()$pch, 
          llty=par()$lty, lcol=par()$col, llwd=par()$lwd,
          label.arg= FALSE, adj.arg=-0.5, 
          sites= FALSE, spch=NULL, scol=par()$col, scex=par()$cex,
          sfont=par()$font,
          whichSpecies = NULL,
          check.ok = TRUE, ...)
{
    type <- match.arg(type, c("fitted.values", "predictors"))[1]

    if((Rank <- object@control$Rank) > 2)
        stop("can only handle rank 1 or 2 models")
    M = if(any(slotNames(object) == "predictors") &&
           is.matrix(object@predictors)) ncol(object@predictors) else
           object@misc$M
    NOS = ncol(object@y)
    MSratio = M / NOS  # First value is g(mean) = quadratic form in lv
    n = object@misc$n
    colx2.index = object@control$colx2.index
    cx1i = object@control$colx1.index
    if(!length(whichSpecies)) whichSpecies = 1:NOS
    if(check.ok)
    if(!(length(cx1i)==1 && names(cx1i)=="(Intercept)"))
        stop("latent variable plots allowable only for Norrr = ~ 1 models")

    Coeflist = Coef(object)
    Cmat = Coeflist@C
    lvmat = Coeflist@lv # n x Rank 

    if(!plot.it) return(lvmat)

    r.curves = slot(object, type) # n times (M or S) (\boldeta or \boldmu) 
    if(MSratio != 1 && type == "predictors")
        stop("can only plot the predictors if M = S")
    MorS = ncol(r.curves) # Actually, here, the value is S always.
    if(!add) {
        if(Rank==1) {
            matplot(lvmat,
                    if( y && type=="fitted.values")
                        object@y[,whichSpecies,drop=FALSE] else
                        r.curves[,whichSpecies,drop=FALSE],
                    type="n", xlab=xlab, ylab=ylab, ...)
        } else { # Rank==2
            matplot(c(Coeflist@Optimum[1,whichSpecies], lvmat[,1]),
                    c(Coeflist@Optimum[2,whichSpecies], lvmat[,2]),
                    type="n", xlab=xlab, ylab=ylab, ...)
        }
    }


    pch  <- rep(pch,  leng=length(whichSpecies))
    pcol <- rep(pcol, leng=length(whichSpecies))
    pcex <- rep(pcex, leng=length(whichSpecies))
    llty <- rep(llty, leng=length(whichSpecies))
    lcol <- rep(lcol, leng=length(whichSpecies))
    llwd <- rep(llwd, leng=length(whichSpecies))
    adj.arg <- rep(adj.arg, leng=length(whichSpecies))

    sppnames = if(type=="predictors") dimnames(r.curves)[[2]] else
        dimnames(object@y)[[2]]
    if(Rank==1) {
        for(sppno in 1:length(whichSpecies)) {
            thisSpecies = whichSpecies[sppno]
            indexSpecies = if(is.character(whichSpecies))
                 match(whichSpecies[sppno], sppnames) else whichSpecies[sppno]
            if(is.na(indexSpecies))
                stop("mismatch found in \"whichSpecies\"")
            xx = lvmat 
            yy = r.curves[,indexSpecies]
            o = sort.list(xx)
            xx = xx[ o ]
            yy = yy[ o ]
            lines(xx, yy, col=lcol[sppno], lwd=llwd[sppno], lty=llty[sppno])
            if( y && type=="fitted.values") {
                ypts = object@y
                if(ncol(as.matrix(ypts)) == ncol(r.curves))
                    points(xx, ypts[o,sppno], col=pcol[sppno],
                           cex=pcex[sppno], pch=pch[sppno])
            } 
        } 
        if(rugplot) rug(xx) 
    } else {
        if(sites) {
            text(lvmat[,1], lvmat[,2], adj=0.5,
                 labels=if(is.null(spch)) dimnames(lvmat)[[1]] else 
                 rep(spch, length=nrow(lvmat)), col=scol, cex=scex, font=sfont)
        }
        for(sppno in 1:length(whichSpecies)) {
            thisSpecies = whichSpecies[sppno]
            indexSpecies = if(is.character(whichSpecies))
                 match(whichSpecies[sppno], sppnames) else whichSpecies[sppno]
            if(is.na(indexSpecies))
                stop("mismatch found in \"whichSpecies\"")
            points(Coeflist@Optimum[1,indexSpecies],
                   Coeflist@Optimum[2,indexSpecies],
                   col=pcol[sppno], cex=pcex[sppno], pch=pch[sppno])
        }
        if(label.arg) {
            for(sppno in 1:length(whichSpecies)) {
                thisSpecies = whichSpecies[sppno]
                indexSpecies = if(is.character(whichSpecies))
                   match(whichSpecies[sppno], sppnames) else whichSpecies[sppno]
                text(Coeflist@Optimum[1,indexSpecies],
                     Coeflist@Optimum[2,indexSpecies],
                     labels=(dimnames(Coeflist@Optimum)[[2]])[indexSpecies], 
                     adj=adj.arg[sppno], col=pcol[sppno], cex=pcex[sppno])
            }
        }
    }
    invisible(lvmat)
}


setMethod("lvplot", "cao",
           function(object, ...) {
           invisible(lvplot.cao(object, ...))})



predict.cao <- function (object, newdata=NULL,
                         type = c("link", "response", "terms"), 
                         deriv = 0, ...) {
    type <- match.arg(type, c("link", "response", "terms"))[1]
    if(type != "link" && deriv != 0)
        stop("Setting deriv=<positive integer> requires type=\"link\"")
    na.act = object@na.action
    object@na.action = list()
    ocontrol = object@control
    nice21 = (length(ocontrol$colx1.index) == 1) &&
             (names(ocontrol$colx1.index) == "(Intercept)")
    if(!nice21) stop("Can only handle Norrr = ~ 1")

    if(!length(newdata) && type=="response" && length(object@fitted.values)) {
        if(length(na.act)) {
            return(napredict(na.act[[1]], object@fitted.values))
        } else {
            return(object@fitted.values)
        }
    }

    attrassignlm <- function(object, ...) 
        attrassigndefault(model.matrix(object), object@terms)

    attrassigndefault <- function(mmat, tt) {
      if (!inherits(tt, "terms"))
        stop("need terms object")
      aa <- attr(mmat, "assign")
      if (is.null(aa))
        stop("argument is not really a model matrix")
      ll <- attr(tt, "term.labels")
      if (attr(tt, "intercept") > 0)
        ll <- c("(Intercept)", ll)
      aaa <- factor(aa, labels = ll)
      split(order(aa), aaa)
    }

    if(!length(newdata)) {
        X <- model.matrixvlm(object, type="lm", ...)
        offset <- object@offset
        tt <- terms(object)
        if(is.R() && !length(object@x))
            attr(X, "assign") <- attrassignlm(X, tt)
    } else {
        if(is.smart(object) && length(object@smart.prediction)) {
            setup.smart("read", smart.prediction=object@smart.prediction)
        }

        tt <- terms(object)  # 11/8/03; object@terms$terms 
        X <- model.matrix(delete.response(tt), newdata, contrasts = 
                      if(length(object@contrasts)) object@contrasts else NULL,
                      xlev = object@xlevels)

        if(is.R() && nice21 && nrow(X)!=nrow(newdata)) {
            as.save = attr(X, "assign")
            X = X[rep(1, nrow(newdata)),,drop=FALSE]
            dimnames(X) = list(dimnames(newdata)[[1]], "(Intercept)")
            attr(X, "assign") = as.save  # Restored 
        }

        offset <- if (!is.null(off.num<-attr(tt,"offset"))) {
            eval(attr(tt,"variables")[[off.num+1]], newdata)
        } else if (!is.null(object@offset))
            eval(object@call$offset, newdata)

        if(is.smart(object) && length(object@smart.prediction)) {
            wrapup.smart() 
        }

        if(is.R())
            attr(X, "assign") <- attrassigndefault(X, tt)
    }

    cancoefs = ccoef(object)

    lvmat = X[,ocontrol$colx2.index,drop=FALSE] %*% cancoefs   # n x Rank

    Rank = ocontrol$Rank
    NOS = ncol(object@y)
    sppnames = dimnames(object@y)[[2]]
    modelno = ocontrol$modelno  # 1,2,3,5 or 0
    M = if(any(slotNames(object) == "predictors") &&
           is.matrix(object@predictors)) ncol(object@predictors) else
           object@misc$M
    MSratio = M / NOS  # First value is g(mean) = quadratic form in lv
    if(type == "terms") {
        terms.mat = matrix(0, nrow(X), Rank*NOS) # 1st R colns for spp.1, etc.
        interceptvector = rep(0, len=NOS)
    } else {
        etamat = matrix(0, nrow(X), M)  # Could contain derivatives
    }
    ind8 = 1:Rank
    whichSpecies = 1:NOS  # Do it all for all species
    for(sppno in 1:length(whichSpecies)) {
        thisSpecies = whichSpecies[sppno]
        indexSpecies = if(is.character(whichSpecies))
            match(whichSpecies[sppno], sppnames) else whichSpecies[sppno]
        if(is.na(indexSpecies))
            stop("mismatch found in \"whichSpecies\"")

        temp345 = predictcao(object, grid=lvmat, sppno=thisSpecies,
                             Rank=Rank, deriv=deriv, MSratio=MSratio,
                             type=ifelse(type=="response", "link", type))
        if(MSratio == 2) {
            if(any(type == c("link", "response"))) {
                etamat[,2*sppno-1] = temp345$yvals 
                etamat[,2*sppno  ] = temp345$eta2 
            } else {
                terms.mat[,ind8] = temp345
                interceptvector[sppno] = attr(temp345, "constant")
            }
        } else {
            if(any(type == c("link", "response"))) {
                etamat[,sppno] = temp345$yvals 
            } else {
                terms.mat[,ind8] = temp345
                interceptvector[sppno] = attr(temp345, "constant")
            }
        }
        ind8 = ind8 + Rank
    }

    if(length(offset) && any(offset != 0))
        etamat <- etamat + offset

    if(type == "link") {
        dimnames(etamat) = list(dimnames(X)[[1]], if(deriv == 0) 
                                object@misc$predictors.names else NULL)
        return(etamat)
    } else if(type == "response") {
        fv <- object@family@inverse(etamat, extra=object@extra)
        dimnames(fv) = list(dimnames(fv)[[1]],
                            dimnames(object@fitted.values)[[2]])
        return(fv)
    } else {
        attr(terms.mat, "constant") = interceptvector
        terms.mat
    }
}



setMethod("predict", "cao", function(object, ...)
           predict.cao(object, ...))


predictcao <- function(object, grid, sppno, Rank=1, deriv=0, MSratio=1,
                       type="link") {
    if(type != "link" && type != "terms")
        stop("'link' must be \"link\" or \"terms\"")
    if(ncol(grid <- as.matrix(grid)) != Rank)
        stop(paste("'grid' must have", Rank, "columns"))
    if(!is.Numeric(1+deriv, allow=1, positive=TRUE, integ=TRUE))
        stop("'deriv' must be a non-negative integer")
    if(type == "terms" && deriv != 0)
        stop("'deriv' must be 0 when type=\"terms\"")

    temp.b = object@Bspline[[sppno]]
    if(type == "terms") {
        meanlv = apply(grid, 2, mean)
        answer = matrix(0, nrow(grid), Rank)
    } else {
        nlfunvalues = 0
    }
    for(rindex in 1:Rank) {
        temp = temp.b[[rindex]]  # temp is of class "vsmooth.spline.fit"
        nlpart = predict(temp, grid[,rindex], deriv=deriv)
        yvals = nlpart$y
        if(type == "terms") {
            answer[,rindex] = yvals
        } else {
            nlfunvalues = nlfunvalues + yvals
        }
    }

    # Get the linear part of the additive predictor (intercept and slopes)
        lcoef = object@coefficients # linear coeffs; dont use coef() (==Coef)
        llcoef = lcoef[(1+(sppno-1)*(MSratio+Rank)):(sppno*(MSratio+Rank))]
        if(type == "terms") {
            interceptvector = llcoef[1]
            for(rindex in 1:Rank) {
                answer[,rindex] = answer[,rindex] +
                    (grid[,rindex] - meanlv[rindex]) * llcoef[MSratio+rindex]
                interceptvector = interceptvector +
                    meanlv[rindex] * llcoef[MSratio+rindex]
            }
        } else {
            linpar = if(deriv==0) {llcoef[1]+grid %*% llcoef[-(1:MSratio)]} else
                {if(deriv==1) llcoef[MSratio+rindex] else 0}
            nlfunvalues = nlfunvalues + linpar # Now complete
        }
    if(type == "terms") {
        attr(answer, "constant") = interceptvector
        answer
    } else {
        list(xvals = grid,
             yvals = c(nlfunvalues),
             eta2 = if(MSratio == 2) llcoef[MSratio] else NULL)
    }
}




plot.cao = function(x,
                    xlab=if(Rank==1) "Latent Variable" else 
                         paste("Latent Variable", 1:Rank),
                    ylab=NULL, residuals.arg=FALSE,
                    pcol=par()$col, pcex=par()$cex, pch=par()$pch,
                    lcol=par()$col, lwd=par()$lwd, lty=par()$lty, 
                    add=FALSE, 
                    main=NULL,
                    center.cf = Rank > 1,
                    WhichRank = 1:Rank, 
                    whichSpecies = NULL, # a numeric or character vector
                    rugplot=TRUE, se.arg=FALSE, deriv=0,
                    scale=0, ylim=NULL,
                    overlay = FALSE, ...)
{
    Rank = x@control$Rank
    if(!is.logical(center.cf) || length(center.cf) != 1)
        stop("bad input for argument 'center.cf'")
    if(Rank > 1 &&  !center.cf)
        stop("center.cf=TRUE is needed for models with Rank>1")
    NOS = ncol(x@y)
    sppnames = dimnames(x@y)[[2]]
    modelno = x@control$modelno  # 1,2,3, or 0
    M = if(any(slotNames(x) == "predictors") &&
           is.matrix(x@predictors)) ncol(x@predictors) else x@misc$M
    if(all((MSratio <- M / NOS) != c(1,2))) stop("bad value for 'MSratio'")
    pcol = rep(pcol, length=Rank*NOS)
    pcex = rep(pcex, length=Rank*NOS)
    pch  = rep(pch,  length=Rank*NOS)
    lcol = rep(lcol, length=Rank*NOS)
    lwd  = rep(lwd,  length=Rank*NOS)
    lty  = rep(lty,  length=Rank*NOS)
    xlab = rep(xlab, length=Rank)
    if(!length(whichSpecies)) whichSpecies = 1:NOS
    if(length(ylab)) 
        ylab = rep(ylab, len=length(whichSpecies)) # Too long if overlay
    if(length(main))
         main = rep(main, len=length(whichSpecies)) # Too long if overlay
    lvmat = lv(x)
    nice21 = length(x@control$colx1.index) == 1 &&
                    names(x@control$colx1.index) == "(Intercept)"
    if(!nice21)
        stop("can only handle intercept-only models")
    counter = 0
    for(sppno in 1:length(whichSpecies)) {
        thisSpecies = whichSpecies[sppno]
        indexSpecies = if(is.character(whichSpecies))
            match(whichSpecies[sppno], sppnames) else whichSpecies[sppno]
        if(is.na(indexSpecies))
            stop("mismatch found in \"whichSpecies\"")
        terms.mat = predictcao(object=x, grid=lvmat, type="terms",
                               sppno=indexSpecies, Rank=Rank,
                               deriv=deriv, MSratio=MSratio)
        for(rindex in WhichRank) {
            xvals = lvmat[,rindex]
            yvals = terms.mat[,rindex]
            o = sort.list(xvals)
            xvals = xvals[ o ]
            yvals = yvals[ o ]
            if(!center.cf) yvals = yvals + attr(terms.mat, "constant")
            if(!add)
            if(sppno==1 || !overlay) {
                ylim.use = if(length(ylim)) ylim else
                    ylim.scale(range(yvals), scale)
                matplot(xvals, yvals, type="n", 
                        xlab=xlab[rindex], 
                        ylab=if(length(ylab)) ylab[sppno] else 
                        ifelse(overlay, "Fitted functions", "Fitted function"),
                        main=if(length(main)) main[sppno] else 
                             ifelse(overlay, "", sppnames[thisSpecies]),
                        ylim=ylim.use,
                        ...)
            }
            if(residuals.arg) {
                stop("can't handle residuals=TRUE yet")
            } 
            counter = counter + 1
            lines(xvals, yvals,
                  col=lcol[counter], lwd=lwd[counter], lty=lty[counter])
            if(rugplot) rug(xvals)
        }
    }
    invisible(x)
}




setMethod("plot", "cao",
           function(x, y, ...) {
           if(!missing(y)) stop("can't process the \"y\" argument")
           invisible(plot.cao(x, ...))})



persp.cao = function(x,
                     plot.it=TRUE,
                     xlim=NULL, ylim=NULL, zlim=NULL, # zlim ignored if Rank==1
                     gridlength=if(Rank==1) 301 else c(51,51),
                     whichSpecies = NULL,
                    xlab=if(Rank==1) "Latent Variable" else "Latent Variable 1",
                    ylab=if(Rank==1) "Expected Value" else "Latent Variable 2",
                     zlab="Expected value",
                     labelSpecies = FALSE,   # For Rank==1 only
                     stretch = 1.05,  # quick and dirty, Rank==1 only
                     main="",
                     ticktype = "detailed",
                     col = if(Rank==1) par()$col else "white",
                     lty=par()$lty,
                     lwd=par()$lwd,
                     rugplot=FALSE,
                     ...) {
    object = x  # don't like x as the primary argument 
    coefobj = Coef(object) 
    if((Rank <- coefobj@Rank) > 2)
        stop("object must be a rank-1 or rank-2 model")
    fvmat = fitted(object)
    NOS = ncol(fvmat)    # Number of species
    M = if(any(slotNames(object) == "predictors") &&
           is.matrix(object@predictors)) ncol(object@predictors) else
           object@misc$M
    MSratio = M / NOS  # First value is g(mean) = quadratic form in lv

    xlim = if(length(xlim)) xlim else range(coefobj@lv[,1])
    if(!length(ylim)) {
        ylim = if(Rank==1) c(0, max(fvmat)*stretch) else range(coefobj@lv[,2])
    }
    xlim = rep(xlim, length=2)
    ylim = rep(ylim, length=2)
    gridlength = rep(gridlength, length=Rank)
    lv1 = seq(xlim[1], xlim[2], length=gridlength[1])
    lv2 = if(Rank == 2) seq(ylim[1], ylim[2], len=gridlength[2]) else NULL
    lvmat = if(Rank == 2) expand.grid(lv1, lv2) else cbind(lv1)

    sppNames = dimnames(object@y)[[2]]
    if(!length(whichSpecies)) {
        whichSpecies = sppNames[1:NOS]
        whichSpecies.numer = 1:NOS
    } else
    if(is.numeric(whichSpecies)) {
        whichSpecies.numer = whichSpecies
        whichSpecies = sppNames[whichSpecies.numer]  # Convert to character
    } else
        whichSpecies.numer = match(whichSpecies, sppNames)

    LP = matrix(as.numeric(NA), nrow(lvmat), NOS) # For first eta for each spp.
    for(sppno in 1:NOS) {
        temp = predictcao(object=object, grid=lvmat, sppno=sppno, 
                          Rank=Rank, deriv=0, MSratio=MSratio)
        LP[,sppno] = temp$yval
    }
    if(MSratio == 2) {
        LP = kronecker(LP, matrix(1:0, 1, 2))  # n x M
    }
    fitvals = object@family@inverse(LP, extra=object@extra)   # n by NOS
    dimnames(fitvals) = list(NULL, dimnames(fvmat)[[2]])

    if(Rank==1) {
        if(plot.it) {
            ylim = c(0, max(fitvals[,whichSpecies.numer])*stretch) # A revision
            col = rep(col, len=length(whichSpecies.numer))
            lty = rep(lty, len=length(whichSpecies.numer))
            lwd = rep(lwd, len=length(whichSpecies.numer))
            matplot(lv1, fitvals, xlab=xlab, ylab=ylab,
                    type="n", main=main, xlim=xlim, ylim=ylim, ...)
            if(rugplot) rug(lv(object)) 
            for(sppno in 1:length(whichSpecies.numer)) {
                ptr2 = whichSpecies.numer[sppno]  # points to species column
                lines(lv1, fitvals[,ptr2], col=col[sppno], 
                      lty=lty[sppno], lwd=lwd [sppno], ...)
                if(labelSpecies) {
                    ptr1=(1:nrow(fitvals))[max(fitvals[,ptr2])==fitvals[,ptr2]]
                    ptr1 = ptr1[1]
                    text(lv1[ptr1], fitvals[ptr1,ptr2]+(stretch-1) *
                         diff(range(ylim)), label=sppNames[sppno],
                         col=col[sppno], ...)
                }
            }
        }
    } else {
        maxfitted = matrix(fitvals[,whichSpecies[1]], length(lv1), length(lv2))
        if(length(whichSpecies) > 1)
        for(sppno in whichSpecies[-1]) {
            maxfitted = pmax(maxfitted, matrix(fitvals[,sppno], 
                                               length(lv1), length(lv2)))
        }
        if(!length(zlim))
            zlim = range(maxfitted, na.rm = TRUE)
        if(plot.it)
            graphics:::persp.default(lv1, lv2, maxfitted,
                  zlim=zlim,
                  xlab=xlab, ylab=ylab, zlab=zlab,
                  ticktype = ticktype, col = col, main=main, ...) 
    }

    invisible(list(fitted=fitvals,
                   lv1grid=lv1,
                   lv2grid=if(Rank==2) lv2 else NULL,
                   maxfitted=if(Rank==2) maxfitted else NULL))
}


if(!isGeneric("persp"))
    setGeneric("persp", function(x, ...) standardGeneric("persp"))
setMethod("persp", "cao", function(x, ...) persp.cao(x=x, ...))



lv.cao = function(object, ...) {
    Coef(object, ...)@lv
}
lv.Coef.cao = function(object, ...) {
    if(length(list(...))) warning("Too late! Ignoring the extra arguments")
    object@lv
}


if(!isGeneric("lv"))
    setGeneric("lv", function(object, ...) standardGeneric("lv"))
setMethod("lv",  "cao", function(object, ...) lv.cao(object, ...))
setMethod("lv", "Coef.cao", function(object, ...) lv.Coef.cao(object, ...))





setClass(Class="summary.cao", representation("Coef.cao",
         "misc" = "list",
         "call" = "call"))

summary.cao = function(object, ...) {
    answer = Coef(object, ...)
    class(answer) = "summary.cao"
    answer@misc = object@misc
    answer@call = object@call
    answer
}

printsummary.cao = function(x, ...) {
    cat("\nCall:\n")
    dput(x@call)

    printCoef.cao(x, ...)

    cat("\nNumber of species: ", x@NOS, "\n")

    if(length(x@misc$dispersion) == 1) {
        cat("\nDispersion parameter(s): ", x@misc$dispersion, "\n")
    } else if(is.Numeric(x@dispersion)) {
        cat("\nDispersion parameter(s)\n")
        print( x@dispersion, ... )
    }
    invisible(x)
}

setMethod("summary", "cao", function(object, ...)
    summary.cao(object, ...))

setMethod("print", "summary.cao",
          function(x, ...)
          invisible(printsummary.cao(x, ...)))

setMethod("show", "summary.cao",
          function(object)
          invisible(printsummary.cao(object)))




ccoef.cao = function(object, ...) {
    Coef(object, ...)@C
}

ccoef.Coef.cao = function(object, ...) {
    if(length(list(...))) warning("Too late! Ignoring the extra arguments")
    object@C
}


if(!isGeneric("ccoef"))
    setGeneric("ccoef", function(object, ...) standardGeneric("ccoef"))
setMethod("ccoef", "cao", function(object, ...) ccoef.cao(object, ...))
setMethod("ccoef", "Coef.cao", function(object, ...) ccoef.Coef.cao(object, ...))


if(!isGeneric("calibrate"))
    setGeneric("calibrate", function(object, ...) standardGeneric("calibrate"))
setMethod("calibrate", "cao", function(object, ...)
          calibrate.qrrvglm(object, ...))

    
setMethod("calibrate", "qrrvglm", function(object, ...)
          calibrate.qrrvglm(object, ...))


Tol.cao = function(object, ...) {
    stop("The tolerance for a \"cao\" object is undefined")
}

if(!isGeneric("Tol"))
    setGeneric("Tol", function(object, ...) standardGeneric("Tol"))
setMethod("Tol", "cao", function(object, ...)
          Tol.cao(object, ...))






setMethod("show",  "cao", function(object) print.vgam(object))
setMethod("print", "cao", function(x, ...) print.vgam(x, ...))




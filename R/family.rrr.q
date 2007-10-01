# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.





replace.constraints <- function(Blist, cm, index)
{

    for(i in index)
        Blist[[i]] = cm
    Blist
}


valt.control <- function(
                 Alphavec=c(2, 4, 6, 9, 12, 16, 20, 25, 30, 40, 50,
                            60, 80, 100, 125, 2^(8:12)),
                 Criterion = c("rss", "coefficients"),
                 Linesearch= FALSE, Maxit=7,
                 Suppress.warning=TRUE,
                 Tolerance=1e-7, ...)
{

    if(mode(Criterion) != "character" && mode(Criterion) != "name")
        Criterion <- as.character(substitute(Criterion))
    Criterion <- match.arg(Criterion, c("rss", "coefficients"))[1]

    list(Alphavec=Alphavec,
         Criterion = Criterion, 
         Linesearch=Linesearch,
         Maxit=Maxit,
         Suppress.warning=Suppress.warning,
         Tolerance=Tolerance)

}


qrrvglm.xprod = function(numat, Aoffset, Quadratic, ITolerances) {
    Rank = ncol(numat)
    moff = NULL
    ans = if(Quadratic) {
            index = iam(NA, NA, M=Rank, diagonal=TRUE, both=TRUE) 
            temp1 = cbind(numat[,index$row] * numat[,index$col])
            if(ITolerances) {
                moff = 0
                for(ii in 1:Rank)
                    moff = moff - 0.5 * temp1[,ii]
            }
            cbind(numat, if(ITolerances) NULL else temp1)
    } else 
        as.matrix(numat)
    list(matrix=if(Aoffset>0) ans else ans[,-(1:Rank),drop=FALSE],
         offset = moff)
}


valt <- function(x, z, U, Rank=1,
                 Blist=NULL, 
                 Cinit=NULL,
                 Alphavec=c(2, 4, 6, 9, 12, 16, 20, 25, 30, 40, 50,
                           60, 80, 100, 125, 2^(8:12)),
                 Criterion=c("rss", "coefficients"),
                 Crow1positive = rep(TRUE, len=Rank),
                 colx1.index,
                 Linesearch=FALSE,
                 Maxit=20, 
                 Structural.zero=NULL,
                 SD.Cinit=0.02,
                 Suppress.warning= FALSE,
                 Tolerance=1e-6, 
                 trace= FALSE,
                 xij=NULL)
{





    if(mode(Criterion) != "character" && mode(Criterion) != "name")
        Criterion <- as.character(substitute(Criterion))
    Criterion <- match.arg(Criterion, c("rss", "coefficients"))[1]

    if(any(diff(Alphavec)) <= 0)
        stop("Alphavec must be an increasing sequence") 

    if(!is.matrix(z))
        z <- as.matrix(z)
    n <- nrow(z)
    M <- ncol(z)
    if(!is.matrix(x))
        x <- as.matrix(x)

    colx2.index = (1:ncol(x))[-colx1.index]
    p1 = length(colx1.index)
    p2 = length(colx2.index)
    p  = p1 + p2
    if(!p2) stop("p2, the dimension of vars for reduced-rank regn, must be > 0")

    if(!length(Blist)) {
        Blist = replace.constraints(vector("list", p), diag(M), 1:p)
    }

    dU <- dim(U)
    if(dU[2] != n)
        stop("input unconformable")

    cmat2 = replace.constraints(vector("list", Rank+p1),
                 if(length(Structural.zero))
                 diag(M)[,-Structural.zero,drop=FALSE] else diag(M), 1:Rank)
    if(p1)
        for(k in 1:p1)
            cmat2[[Rank+k]] <- Blist[[colx1.index[k]]]

    if(is.null(Cinit))
        Cinit <- matrix(rnorm(p2*Rank, sd=SD.Cinit), p2, Rank)

    fit <- list(rss=0)  # Only for initial old.crit below

    C <- Cinit # This is input for the main iter loop
    old.crit <- switch(Criterion, coefficients=C, rss=fit$rss)

    recover = 0  # Allow a few iterations between different line searches 
    for(iter in 1:Maxit) {
        iter.save <- iter

        lv.mat <- x[,colx2.index,drop=FALSE] %*% C
        new.lv.model.matrix = cbind(lv.mat, if(p1) x[,colx1.index] else NULL)
        fit = vlm.wfit(new.lv.model.matrix, z, Blist=cmat2, U=U, 
              matrix.out=TRUE, XBIG=FALSE, rss=FALSE, qr=FALSE, xij=xij)
        A <- t(fit$mat.coef[1:Rank,,drop=FALSE])

        cmat1 = replace.constraints(Blist, A, colx2.index)
        fit = vlm.wfit(x, z, Blist=cmat1, U=U, 
              matrix.out=TRUE, XBIG= FALSE, rss=TRUE, qr= FALSE, xij=xij)
        C = fit$mat.coef[colx2.index,,drop=FALSE] %*% A %*% solve(t(A) %*% A)

        numat = x[,colx2.index,drop=FALSE] %*% C
        evnu = eigen(var(numat))
        temp7 = if(Rank > 1) evnu$vector %*% diag(evnu$value^(-0.5)) else
                evnu$vector %*% evnu$value^(-0.5)
        C = C %*% temp7
        A = A %*% t(solve(temp7))
        temp8 = crow1C(cmat=C, Crow1positive, amat=A)
        C = temp8$cmat
        A = temp8$amat


        ratio=switch(Criterion,
                     coefficients=max(abs(C - old.crit) / (Tolerance+abs(C))),
                     rss=max(abs(fit$rss - old.crit) / (Tolerance+fit$rss)))

        if(trace) {
            cat("    Alternating iteration", iter,
                "    ratio =", format(ratio), "\n")
            if(!is.null(fit$rss))
                cat("    rss =", fit$rss, "\n")
            if(exists("flush.console"))
                flush.console()
        }

        if(ratio < Tolerance) {
            if(!Linesearch || (Linesearch && iter >= 3)) break
        } else if(iter == Maxit && !Suppress.warning) {
            warning("did not converge")
        }

        fini.linesearch = FALSE
        if(Linesearch && iter - recover >= 2) {
            xnew <- C

            direction1 <- (xnew-xold) # / sqrt(1 + sum((xnew-xold)^2))
            ftemp <- fit$rss  # Most recent objective function 
            use.alpha <- 0   # The current step relative to (xold, yold)
            for(itter in 1:length(Alphavec)) {
                CC <- xold + Alphavec[itter] * direction1

                try.lv.mat <- x[,colx2.index,drop=FALSE] %*% CC
                try.new.lv.model.matrix = cbind(try.lv.mat,
                                   if(p1) x[,colx1.index] else NULL)

                try = vlm.wfit(try.new.lv.model.matrix, z, Blist=cmat2, U=U, 
                               matrix.out=TRUE, XBIG= FALSE, rss=TRUE, qr= FALSE,
                               xij=xij)
                if(try$rss < ftemp) {
                    use.alpha <- Alphavec[itter]
                    fit <- try 
                    ftemp <- try$rss
                    C <- CC 
                    A = t(fit$mat.coef[1:Rank,,drop=FALSE])
                    lv.mat <- x[,colx2.index,drop=FALSE] %*% C
                    recover = iter # Give it some alt'g iterations to recover
                } else {
                    if(trace && use.alpha>0) {
                        cat("    Finished line search using Alpha =", 
                            use.alpha, "\n")
                        if(exists("flush.console"))
                            flush.console()
                    }
                    fini.linesearch = TRUE
                }
                if(fini.linesearch) break 
            } # End of itter loop 
        }

        xold <- C # Don't take care of drift
        old.crit <- switch(Criterion, coefficients=C, rss=fit$rss)
    } # End of iter loop

    list(A=A, C=C, fitted=fit$fitted, new.coeffs = fit$coef, rss=fit$rss)
}



lm2qrrvlm.model.matrix = function(x, Blist, C, control, assign=TRUE,
                                  no.thrills=FALSE)
{

    Rank = control$Rank
    colx1.index = control$colx1.index
    Quadratic = control$Quadratic
    Dzero = control$Dzero
    Corner = control$Corner
    ITolerances = control$ITolerances

    M = nrow(Blist[[1]])
    p1 = length(colx1.index)
    combine2 = c(control$Structural.zero,
                 if(Corner) control$Index.corner else NULL)

    Qoffset = if(Quadratic) ifelse(ITolerances, 0, sum(1:Rank)) else 0
    NoA = length(combine2) == M    # No unknown parameters in A
    cmat2 = if(NoA) {
        Aoffset = 0
        vector("list", Aoffset+Qoffset+p1)
    } else {
        Aoffset = Rank
        replace.constraints(vector("list", Aoffset+Qoffset+p1),
           if(length(combine2)) diag(M)[,-combine2,drop=FALSE] else diag(M),
           1:Rank) # If Corner then doesn't contain \bI_{Rank}
    }
    if(Quadratic && !ITolerances)
        cmat2 = replace.constraints(cmat2,
            if(control$EqualTolerances)
                matrix(1, M, 1) - eij(Dzero, M) else {
            if(length(Dzero)) diag(M)[,-Dzero,drop=FALSE] else diag(M)},
            Aoffset + (1:Qoffset))
    if(p1)
        for(k in 1:p1)
            cmat2[[Aoffset+Qoffset+k]] <- Blist[[colx1.index[k]]]
    if(!no.thrills) {
        i63 = iam(NA, NA, M=Rank, both=TRUE)
        names(cmat2) = c(
               if(NoA) NULL else paste("(lv", 1:Rank, ")", sep=""), 
               if(Quadratic && Rank==1 && !ITolerances)
                   "(lv^2)" else 
               if(Quadratic && Rank>1 && !ITolerances)
                   paste("(lv", i63$row, ifelse(i63$row==i63$col, "^2",
                   paste("*lv", i63$col, sep="")), ")", sep="") else NULL,
               if(p1) names(colx1.index) else NULL)
    }

    lv.mat = x[,control$colx2.index,drop=FALSE] %*% C


    tmp900 = qrrvglm.xprod(lv.mat, Aoffset, Quadratic, ITolerances)
    new.lv.model.matrix = cbind(tmp900$matrix,
                                if(p1) x[,colx1.index] else NULL)
    if(!no.thrills)
        dimnames(new.lv.model.matrix) = list(dimnames(x)[[1]], names(cmat2))

    if(assign) {
        asx = attr(x, "assign")
        asx = vector("list", ncol(new.lv.model.matrix))
        names(asx) = names(cmat2)   # wrong zz 
        for(i in 1:length(names(asx))) {
            asx[[i]] = i
        }
        attr(new.lv.model.matrix, "assign") = asx
    }

    if(no.thrills)
        list(new.lv.model.matrix = new.lv.model.matrix, constraints = cmat2,
             offset = tmp900$offset) else
        list(new.lv.model.matrix = new.lv.model.matrix, constraints = cmat2,
             NoA = NoA, Aoffset = Aoffset, lv.mat = lv.mat,
             offset = tmp900$offset)
}

valt.2iter <- function(x, z, U, Blist, A, control) {


    cmat1 = replace.constraints(Blist, A, control$colx2.index)
    fit <- vlm.wfit(x, z, Blist=cmat1, U=U, matrix.out=TRUE, 
                    XBIG= FALSE, rss=TRUE, qr= FALSE, xij=control$xij)
    C = fit$mat.coef[control$colx2.index,,drop=FALSE] %*% A %*% solve(t(A) %*% A)

    list(A=A, C=C, fitted=fit$fitted, new.coeffs = fit$coef,
         Blist=cmat1, rss=fit$rss)
}



valt.1iter = function(x, z, U, Blist, C, control, lp.names=NULL, nice31=FALSE,
                      MSratio = 1) {

    Rank = control$Rank
    Quadratic = control$Quadratic
    Index.corner = control$Index.corner
    p1 = length(control$colx1.index)
    M = ncol(zedd <- as.matrix(z))
    NOS = M / MSratio
    Corner = control$Corner
    ITolerances = control$ITolerances

    Qoffset = if(Quadratic) ifelse(ITolerances, 0, sum(1:Rank)) else 0
    tmp833 = lm2qrrvlm.model.matrix(x=x, Blist=Blist, C=C, control=control)
    new.lv.model.matrix = tmp833$new.lv.model.matrix 
    cmat2.save = cmat2 = tmp833$constraints     # Doesn't contain \bI_{Rank}
    lv.mat = tmp833$lv.mat
    if(Corner)
        zedd[,Index.corner] = zedd[,Index.corner] - lv.mat

    if(nice31 && MSratio == 1) {
        fit = list(mat.coef = NULL, fitted.values = NULL, rss = 0)

        cmat2 = NULL # for vlm.wfit

        i5 = rep(0, len=MSratio)
        for(ii in 1:NOS) {
            i5 = i5 + 1:MSratio

            tmp100 = vlm.wfit(new.lv.model.matrix, zedd[,i5,drop=FALSE], 
                              Blist=cmat2, U=U[i5,,drop=FALSE],
                              matrix.out=TRUE, XBIG=FALSE, rss=TRUE, qr=FALSE,
                              Eta.range = control$Eta.range,
                              xij=control$xij, lp.names=lp.names[i5])
            fit$rss = fit$rss + tmp100$rss
            fit$mat.coef = cbind(fit$mat.coef, tmp100$mat.coef)
            fit$fitted.values = cbind(fit$fitted.values, tmp100$fitted.values)
        }
    } else {
        fit = vlm.wfit(new.lv.model.matrix, zedd, Blist=cmat2, U=U, 
                       matrix.out=TRUE, XBIG= FALSE, rss=TRUE, qr= FALSE,
                       Eta.range = control$Eta.range,
                       xij=control$xij, lp.names=lp.names)
    }
    A = if(tmp833$NoA) matrix(0, M, Rank) else
        t(fit$mat.coef[1:Rank,,drop=FALSE])
    if(Corner)
        A[Index.corner,] = diag(Rank)     

    B1=if(p1) fit$mat.coef[-(1:(tmp833$Aoffset+Qoffset)),,drop=FALSE] else NULL
    fv = as.matrix(fit$fitted.values)
    if(Corner)
        fv[,Index.corner] = fv[,Index.corner] + lv.mat
    Dmat = if(Quadratic) {
            if(ITolerances) {
                tmp800 = matrix(0, M, Rank*(Rank+1)/2)
                tmp800[if(MSratio==2) c(TRUE,FALSE) else TRUE,1:Rank] = -0.5
                tmp800
            } else 
                t(fit$mat.coef[(tmp833$Aoffset+1):
                  (tmp833$Aoffset+Qoffset),,drop=FALSE])
    } else
        NULL

    list(Amat=A, B1=B1, Cmat=C, Dmat=Dmat, fitted=if(M==1) c(fv) else fv,
         new.coeffs = fit$coef, constraints=cmat2, rss=fit$rss,
         offset = if(length(tmp833$offset)) tmp833$offset else NULL)
}





rrr.init.expression <- expression({
    if(backchat || control$Quadratic) 
        copyxbig <- TRUE

    modelno = switch(family@vfamily[1], "poissonff"=2,
              "quasipoissonff"=2, "quasipoisson"=2,
              "binomialff"=1, "quasibinomialff"=1,
              "quasibinomial"=1, "negbinomial"=3,
              "gamma2"=5, "gaussianff"=8,
              0)  # stop("can't fit this model using fast algorithm")
    if(modelno == 1) modelno = get("modelno", envir = VGAMenv)
    rrcontrol$modelno = control$modelno = modelno
    if(modelno==3 || modelno==5) {


        M = 2 * ifelse(is.matrix(y), ncol(y), 1)
        control$Structural.zero =
        rrcontrol$Structural.zero = seq(from=2, to=M, by=2)  # Handles A
        control$Dzero =
        rrcontrol$Dzero = seq(from=2, to=M, by=2)  # Handles D

    }


})



rrr.alternating.expression <- expression({

    alt <- valt(x, z, U, Rank=Rank,
                Blist=Blist,
                Cinit=rrcontrol$Cinit,
                Criterion=rrcontrol$Criterion,
                colx1.index=rrcontrol$colx1.index,
                Linesearch=rrcontrol$Linesearch,
                Maxit=rrcontrol$Maxit,
                Structural.zero=rrcontrol$Structural.zero,
                SD.Cinit=rrcontrol$SD.Cinit,
                Suppress.warning=rrcontrol$Suppress.warning,
                Tolerance=rrcontrol$Tolerance,
                trace=trace,
                xij=control$xij) # This is subject to drift in A and C

    ans2 = rrr.normalize(rrcontrol=rrcontrol, A=alt$A, C=alt$C, x=x)

    Amat = ans2$A              # Fed into Blist below (in rrr.end.expression)
    tmp.fitted = alt$fitted    # Also fed; was alt2$fitted 

    rrcontrol$Cinit <- ans2$C   # For next valt() call

    eval(rrr.end.expression)    # Put Amat into Blist, and create new z
})

    adjust.Dmat.expression = expression({
    if(length(Dmat)) {
        ind0 = iam(NA, NA, both= TRUE, M=Rank)
        for(kay in 1:M) {
            elts = Dmat[kay,,drop=FALSE] # Manual recycling
            if(length(elts) < Rank)
                elts = matrix(elts, 1, Rank)
            Dk = m2adefault(elts, M=Rank)[,,1]
            Dk = matrix(Dk, Rank, Rank)
            Dk = t(Mmat) %*% Dk  %*% Mmat # 22/8/03; Not diagonal in general
            Dmat[kay,] = Dk[cbind(ind0$row.index[1:ncol(Dmat)],
                                  ind0$col.index[1:ncol(Dmat)])] 
        }
    }})

rrr.normalize = function(rrcontrol, A, C, x, Dmat=NULL) {



    colx2.index = rrcontrol$colx2.index
    Rank = rrcontrol$Rank
    Index.corner = rrcontrol$Index.corner
    M = nrow(A)
    C.old = C

    if(rrcontrol$Corner) {
        tmp87 = A[Index.corner,,drop=FALSE]
        Mmat <- solve(tmp87) # The normalizing matrix
        C <- C %*% t(tmp87)
        A <- A %*% Mmat
        A[Index.corner,] <- diag(Rank)  # Make sure 
        eval(adjust.Dmat.expression)
    }

    if(rrcontrol$Svd.arg) {
        temp = svd(C %*% t(A))
        if(!is.matrix(temp$v))
            temp$v = as.matrix(temp$v) 
        C = temp$u[,1:Rank,drop=FALSE] %*%
            diag(temp$d[1:Rank]^(1-rrcontrol$Alpha), nrow=Rank)
        A = diag(temp$d[1:Rank]^(rrcontrol$Alpha), nrow=Rank) %*%
            t(temp$v[,1:Rank,drop=FALSE])
        A = t(A)
        Mmat = t(C.old)  %*% C.old %*% solve(t(C) %*% C.old)
        eval(adjust.Dmat.expression)
    }

    if(rrcontrol$Uncor) {
        lv.mat <- x[,colx2.index,drop=FALSE] %*% C
        var.lv.mat <- var(lv.mat)
        UU = chol(var.lv.mat)
        Ut <- solve(UU)
        Mmat <- t(UU)
        C <- C %*% Ut
        A <- A %*% t(UU)
        eval(adjust.Dmat.expression)
    }


    if(rrcontrol$Quadratic) {
        Mmat = diag(Rank)
        for(LV in 1:Rank)
            if(( rrcontrol$Crow1positive[LV] && C[1,LV] < 0) ||
               (!rrcontrol$Crow1positive[LV] && C[1,LV] > 0)) {
                C[,LV] = -C[,LV]
                A[,LV] = -A[,LV]
                Mmat[LV,LV] = -1
            }
        eval(adjust.Dmat.expression) # Using Mmat above 
    }


    list(Amat=A, Cmat=C, Dmat=Dmat)
}


rrr.end.expression = expression({

    if(is.R()) {
        if(exists(".VGAM.etamat", envir = VGAMenv))
            rm(".VGAM.etamat", envir = VGAMenv)
    } else {
        while(exists(".VGAM.etamat", inherits=TRUE))
            rm(".VGAM.etamat", inherits=TRUE)
    }


    if(control$Quadratic) {
        if(!length(extra)) extra=list()
        extra$Cmat = Cmat      # Saves the latest iteration 
        extra$Dmat = Dmat      # Not the latest iteration
        extra$B1   = B1.save   # Not the latest iteration (not good)
    } else {
        Blist = replace.constraints(Blist.save, Amat, colx2.index)
    }

    xbig.save = if(control$Quadratic) {
        tmp300 = lm2qrrvlm.model.matrix(x=x, Blist=Blist.save,
                                        C=Cmat, control=control)
        lv.mat = tmp300$lv.mat  # Needed at the top of new.s.call

        lm2vlm.model.matrix(tmp300$new.lv.model.matrix,B.list,xij=control$xij)
    } else {
        lm2vlm.model.matrix(x, Blist, xij=control$xij)
    }


    fv <- tmp.fitted            # Contains \bI \bnu
    eta <- fv + offset
    if(FALSE && control$Rank == 1) {
        ooo = order(lv.mat[,1])
    }
    mu <- family@inverse(eta, extra)

    if(any(is.na(mu)))
        warning("there are NAs in mu") 

    deriv.mu <- eval(family@deriv)
    wz <- eval(family@weight)
    if(control$checkwz)
        wz = checkwz(wz, M=M, trace=trace, wzeps=control$wzepsilon)
    U <- vchol(wz, M=M, n=n, silent=!trace)
    tvfor <- vforsub(U, as.matrix(deriv.mu), M=M, n=n)
    z <- eta + vbacksub(U, tvfor, M=M, n=n) - offset # Contains \bI \bnu



})



rrr.derivative.expression <- expression({






    which.optimizer = if(is.R()) {
        if(control$Quadratic && control$FastAlgorithm) {
            "BFGS" 
        } else {
            if(iter <= rrcontrol$Switch.optimizer) "Nelder-Mead" else "BFGS"
        }
    } else "Quasi-Newton" 
    if(trace && control$OptimizeWrtC) {
        cat("\n\n")
        cat("Using", which.optimizer, "\n")
        if(exists("flush.console"))
            flush.console()
    } 

    constraints=replace.constraints(constraints,diag(M),rrcontrol$colx2.index)
    nice31 = (!control$EqualTol || control$ITolerances) &&
             all(trivial.constraints(constraints))

    theta0 <- c(Cmat) # zz; Possibly bad because of normalization?
    if(is.R()) assign(".VGAM.dot.counter", 0, envir = VGAMenv) else
        .VGAM.dot.counter <<- 0 
if(control$OptimizeWrtC) {
    if(is.R()) {
        if(control$Quadratic && control$FastAlgorithm) {
            if(iter == 2) {
               if(is.R()) {
                   if(exists(".VGAM.etamat", envir = VGAMenv))
                       rm(".VGAM.etamat", envir = VGAMenv)
               } else {
                   if(exists(".VGAM.etamat", inherits=TRUE))
                       rm(".VGAM.etamat", inherits=TRUE)
               }
            }
            if(iter > 2 && !quasi.newton$convergence) {
                if(is.R()) {
                    if(zthere <- exists(".VGAM.z", envir = VGAMenv)) {
                        ..VGAM.z = get(".VGAM.z", envir = VGAMenv)
                        ..VGAM.U = get(".VGAM.U", envir = VGAMenv)
                        ..VGAM.beta = get(".VGAM.beta", envir = VGAMenv)
                    }
                } else {
                    if(zthere <- exists(".VGAM.z")) {
                        ..VGAM.z = .VGAM.z
                        ..VGAM.U = .VGAM.U
                        ..VGAM.beta = .VGAM.beta
                    }
                }
                if(zthere) {
                    z = matrix(..VGAM.z, n, M)  # minus any offset
                    U = matrix(..VGAM.U, M, n)
                }

            }
    
            if(iter == 2 || quasi.newton$convergence) {
                NOS = ifelse(modelno==3 || modelno==5, M/2, M)

                canfitok = if(is.R()) 
                    (exists("CQO.FastAlgorithm", envir=VGAMenv) &&
                    get("CQO.FastAlgorithm", envir = VGAMenv)) else
                (exists("CQO.FastAlgorithm",inherits=TRUE) && CQO.FastAlgorithm)
                if(!canfitok)
                    stop("can't fit this model using fast algorithm")
                p2star = if(nice31) 
                  ifelse(control$IToleran, Rank, Rank+0.5*Rank*(Rank+1)) else
                  (NOS*Rank + Rank*(Rank+1)/2 * ifelse(control$EqualTol,1,NOS))
                p1star = if(nice31) p1 * ifelse(modelno==3 || modelno==5,2,1) else
                         (ncol(xbig.save) - p2star)
                xbig.save1 = if(p1star > 0) xbig.save[,-(1:p2star)] else NULL
                quasi.newton = optim(par=Cmat, fn=callcqof, 
                        gr=if(control$GradientFunction) calldcqof else NULL,
                        method=which.optimizer,
                        control=list(fnscale=1,trace=as.integer(control$trace),
                            parscale=rep(control$Parscale, len=length(Cmat)),
                            maxit=250),
                        etamat=eta, xmat=x, ymat=y, wvec=w,
                        xbig.save1 = if(nice31) NULL else xbig.save1,
                        modelno=modelno, Control=control,
                        n=n, M=M, p1star=p1star, p2star=p2star, nice31=nice31)


                if(is.R()) {
                    if(zthere <- exists(".VGAM.z", envir = VGAMenv)) {
                        ..VGAM.z = get(".VGAM.z", envir = VGAMenv)
                        ..VGAM.U = get(".VGAM.U", envir = VGAMenv)
                        ..VGAM.beta = get(".VGAM.beta", envir = VGAMenv)
                    }
                } else {
                    if(zthere <- exists(".VGAM.z")) {
                        ..VGAM.z = .VGAM.z
                        ..VGAM.U = .VGAM.U
                        ..VGAM.beta = .VGAM.beta
                    }
                }
                if(zthere) {
                    z = matrix(..VGAM.z, n, M)  # minus any offset
                    U = matrix(..VGAM.U, M, n)
                }
            } else {
                if(is.R()) {
                    if(exists(".VGAM.offset", envir = VGAMenv))
                        rm(".VGAM.offset", envir = VGAMenv)
                } else {
                    while(exists(".VGAM.offset", inherits=TRUE))
                        rm(".VGAM.offset", inherits=TRUE)
                }
            }
        } else {
            use.reltol = if(length(rrcontrol$Reltol) >= iter) 
                rrcontrol$Reltol[iter] else rev(rrcontrol$Reltol)[1]
            quasi.newton <-
            optim(par=theta0,
                  fn=rrr.derivC.rss, 
                  method=which.optimizer,
                  control=list(fnscale=rrcontrol$Fnscale, 
                               maxit=rrcontrol$Maxit,
                               abstol=rrcontrol$Abstol,
                               reltol=use.reltol),
                  U=U, z= if(control$ITolerances) z+offset else z,
                  M=M, xmat=x,    # varbix2=varbix2,
                  Blist=Blist, rrcontrol=rrcontrol)
        }
    } else {
        quasi.newton <-
        nlminb(start=theta0,
               objective=rrr.derivC.rss, 
               control = nlminb.control(x.tol = rrcontrol$X.tol,
                                  eval.max = rrcontrol$Eval.max,
                                  iter.max = rrcontrol$Iter.max,
                                  abs.tol = rrcontrol$Abs.tol,
                                  rel.tol = rrcontrol$Rel.tol,
                                  step.min = rrcontrol$Step.min,
                                  rel.err = rrcontrol$Rel.err),
               U=U, z= if(control$ITolerances) z+offset else z,
               M=M, xmat=x, # varbix2=varbix2,
               Blist=Blist, rrcontrol=rrcontrol)
    }



    Cmat = matrix(quasi.newton$par, p2, Rank, byrow=FALSE)

    if(Rank > 1 && rrcontrol$ITolerances) {
            numat = x[,rrcontrol$colx2.index,drop=FALSE] %*% Cmat
            evnu = eigen(var(numat))
            Cmat = Cmat %*% evnu$vector
            numat = x[,rrcontrol$colx2.index,drop=FALSE] %*% Cmat
            offset = if(Rank > 1) -0.5*apply(numat^2, 1, sum) else -0.5*numat^2
    }
}


    alt = valt.1iter(x=x, z=z, U=U, Blist=Blist, C=Cmat, nice31=nice31,
                     control=rrcontrol, lp.names=predictors.names)


    if(length(alt$offset))
        offset = alt$offset

    B1.save = alt$B1 # Put later into extra  
    tmp.fitted = alt$fitted  # contains \bI_{Rank} \bnu if Corner

    if(modelno!=33 && control$OptimizeWrtC)
        alt = rrr.normalize(rrc=rrcontrol, A=alt$Amat, C=alt$Cmat, 
                            x=x, Dmat=alt$Dmat)

    if(trace && control$OptimizeWrtC) {
        cat("\n")
        cat(which.optimizer, "using",
            if(is.R()) "optim():" else "nlminb():", "\n")
        cat("Objective =", if(is.R()) 
             quasi.newton$value else format(quasi.newton$objective), "\n")
        cat("Parameters (= c(C)) = ", if(length(quasi.newton$par) < 5)
            "" else "\n")
        cat(if(is.R()) alt$Cmat else format(alt$Cmat), fill=TRUE)
        cat("\n")
        if(!is.R())
            cat("Gradient norm =", format(quasi.newton$grad.norm), "\n")
        cat("Number of function evaluations =", if(is.R()) 
             quasi.newton$count[1] else quasi.newton$f.evals, "\n")
        if(!is.R())
            cat("Number of gradient evaluations =", quasi.newton$g.evals, "\n")
        if(length(quasi.newton$message))
            cat("Message =", quasi.newton$message, "\n")
        cat("\n")
        if(exists("flush.console"))
            flush.console()
    }



    Amat = alt$Amat  # Needed in rrr.end.expression 
    Cmat = alt$Cmat  # Needed in rrr.end.expression if Quadratic 
    Dmat = alt$Dmat  # Put later into extra  

    eval(rrr.end.expression)    # Put Amat into Blist, and create new z
})


rrr.derivC.rss = function(theta, U, z, M, xmat, Blist, rrcontrol,
                          omit.these=NULL) {

    if(rrcontrol$trace) {
        cat(".")
        if(exists("flush.console"))
            flush.console()
    }
    alreadyThere = if(is.R())
        exists(".VGAM.dot.counter", envir = VGAMenv) else
        exists(".VGAM.dot.counter")
    if(alreadyThere) {
        if(is.R()) {
            VGAM.dot.counter = get(".VGAM.dot.counter", envir = VGAMenv)
            VGAM.dot.counter = VGAM.dot.counter + 1 
            assign(".VGAM.dot.counter", VGAM.dot.counter, envir = VGAMenv)
        } else {
            .VGAM.dot.counter <<- .VGAM.dot.counter + 1
        }
        if(VGAM.dot.counter > max(50, options()$width - 5)) {
            if(rrcontrol$trace) {
                cat("\n")
                if(exists("flush.console"))
                    flush.console()
            }
            if(is.R()) assign(".VGAM.dot.counter", 0, envir = VGAMenv) else
                .VGAM.dot.counter <<- 0
        }
    }

    Cmat = matrix(theta, length(rrcontrol$colx2.index), rrcontrol$Rank)


    tmp700 = lm2qrrvlm.model.matrix(x=xmat, Blist=Blist,
                   no.thrills = !rrcontrol$Corner,
                   C=Cmat, control=rrcontrol, assign= FALSE)
    Blist = tmp700$constraints # Doesn't contain \bI_{Rank} \bnu

    if(rrcontrol$Corner) {
        z = as.matrix(z) # should actually call this zedd
        z[,rrcontrol$Index.corner] = z[,rrcontrol$Index.corner] - tmp700$lv.mat
    }

    if(length(tmp700$offset)) z = z - tmp700$offset


    vlm.wfit(x=tmp700$new.lv.model.matrix, z=z,
             Blist=Blist, ncolx=ncol(xmat), U=U, only.rss=TRUE,
             matrix.out= FALSE, XBIG= FALSE, rss= TRUE, qr= FALSE,
             Eta.range = rrcontrol$Eta.range,
             xij=rrcontrol$xij)$rss
}



rrvglm.optim.control = function(Fnscale=1,
                                Maxit=100, 
                                Switch.optimizer=3,
                                Abstol= -Inf, 
                                Reltol=sqrt(.Machine$double.eps),
                                ...)
{




    list(Fnscale=Fnscale, 
         Maxit=Maxit,
         Switch.optimizer=Switch.optimizer,
         Abstol=Abstol,
         Reltol=Reltol)
}



if(is.R())
nlminbcontrol = function(Abs.tol = 10^(-6),
                       Eval.max=91,
                       Iter.max=91,
                       Rel.err = 10^(-6),
                       Rel.tol = 10^(-6),
                       Step.min = 10^(-6),
                       X.tol = 10^(-6),
                       ...)
{


    list(Abs.tol = Abs.tol,
         Eval.max=Eval.max,
         Iter.max = Iter.max,
         Rel.err=Rel.err,
         Rel.tol=Rel.tol,
         Step.min=Step.min,
         X.tol=X.tol)
}




Coef.qrrvglm <- function(object, varlvI = FALSE, reference = NULL, ...) {


    if(length(varlvI) != 1 || !is.logical(varlvI)) 
        stop("\"varlvI\" must be TRUE or FALSE")
    if(length(reference) > 1) stop("\"reference\" must be of length 0 or 1")
    if(length(reference) && is.Numeric(reference))
        if(!is.Numeric(reference, allow=1, integ=TRUE))
            stop("bad input for argument \"reference\"")
    if(!is.logical(ConstrainedQO <- object@control$ConstrainedQO))
        stop("can't determine whether the model is constrained or not")
    ocontrol = object@control
    coef.object = object@coefficients 
    Rank = ocontrol$Rank 
    M = object@misc$M
    NOS = if(length(object@y)) ncol(object@y) else M
    MSratio = M / NOS  # First value is g(mean) = quadratic form in lv
    Quadratic = if(ConstrainedQO) ocontrol$Quadratic else TRUE
    if(!Quadratic) stop("object is not a quadratic ordination object")
    p1 = length(ocontrol$colx1.index)
    p2 = length(ocontrol$colx2.index)
    Index.corner = ocontrol$Index.corner
    Structural.zero = ocontrol$Structural.zero
    EqualTolerances = ocontrol$EqualTolerances
    Dzero = ocontrol$Dzero
    Corner = if(ConstrainedQO) ocontrol$Corner else FALSE
    estITol = if(ConstrainedQO) object@control$ITolerances else FALSE
    modelno = object@control$modelno  # 1,2,3,4,5,6,7 or 0
    combine2 = c(Structural.zero, if(Corner) Index.corner else NULL)
    NoA = length(combine2) == M # A is fully known # doesn't handle !Corner yet

    Qoffset = if(Quadratic) ifelse(estITol, 0, sum(1:Rank)) else 0

    ynames = object@misc$ynames
    if(!length(ynames)) ynames = object@misc$predictors.names
    if(!length(ynames)) ynames = object@misc$ynames
    if(!length(ynames)) ynames = paste("Y", 1:NOS, sep="")
    lp.names = object@misc$predictors.names
    if(!length(lp.names)) lp.names = NULL 

    dzero.vector = rep(FALSE, length=M)
    if(length(Dzero))
        dzero.vector[Dzero] = TRUE
    names(dzero.vector) = ynames 
    lv.names = if(Rank==1) "lv" else paste("lv", 1:Rank, sep="")

    td.expression = expression({
        Tolerance = Darray = m2adefault(Dmat, M=Rank)
        for(i in 1:M)
            if(length(Dzero) && any(Dzero == i)) {
                Tolerance[,,i] = NA   # Darray[,,i] == O 
                bellshaped[i] = FALSE 
            } else {
                Tolerance[,,i] = -0.5 * solve(Darray[,,i])
                bellshaped[i] = all(eigen(Tolerance[,,i])$values > 0)
            }
        optimum = matrix(as.numeric(NA),Rank,M) # dimnames=list(lv.names,ynames)
        for(i in 1:M)
            if(bellshaped[i])
                optimum[,i] = Tolerance[,,i] %*% cbind(Amat[i,])
    })
    Amat = object@extra$Amat   # M  x Rank
    Cmat = object@extra$Cmat   # p2 x Rank
    Dmat = object@extra$Dmat   #
    B1   = object@extra$B1     #
    bellshaped = rep(FALSE, length=M)

    if(is.character(reference)) {
        reference = (1:NOS)[reference == ynames]
        if(length(reference) != 1)
           stop("could not match argument \"reference\" with any response")
    }
    ptr1 = 1
    candidates = if(length(reference)) reference else {
        if(length(ocontrol$Dzero)) (1:M)[-ocontrol$Dzero] else (1:M)}
    repeat {
        if(ptr1 > 0) {
            this.spp = candidates[ptr1]
        }
        elts = Dmat[this.spp,,drop=FALSE]
        if(length(elts) < Rank)
            elts = matrix(elts, 1, Rank)
        Dk = m2adefault(elts, M=Rank)[,,1]    # Hopefully negative-def 
        temp400 = eigen(Dk)
        ptr1 = ptr1 + 1 
        if(all(temp400$value < 0)) break
        if(ptr1 > length(candidates)) break
    }
    if(all(temp400$value < 0)) {
        temp1tol = -0.5 * solve(Dk)
        dim(temp1tol) = c(Rank,Rank)
        Mmat = t(chol(temp1tol))
        if(ConstrainedQO) {
            temp900 = solve(t(Mmat))
            Cmat = Cmat %*% temp900
            Amat = Amat %*% Mmat
        }
        if(length(Cmat)) {
            temp800 = crow1C(Cmat, ocontrol$Crow1positive, amat=Amat)
            Cmat = temp800$cmat
            Amat = temp800$amat
        }
        eval(adjust.Dmat.expression)
        eval(td.expression)
    } else {
        if(length(reference) == 1) 
            stop(paste("tolerance matrix specified by \"reference\"",
                       "is not positive-definite")) else
            warning(paste("could not find any positive-definite",
                          "tolerance matrix"))
    }


    if(ConstrainedQO)
    if(Rank > 1) {
        if(!length(xmat <- object@x)) stop("cannot obtain the model matrix")
        numat = xmat[,ocontrol$colx2.index,drop=FALSE] %*% Cmat
        evnu = eigen(var(numat))
        Mmat = solve(t(evnu$vector))
        Cmat = Cmat %*% evnu$vector  # == Cmat %*% solve(t(Mmat))
        Amat = Amat %*% Mmat
        temp800 = crow1C(Cmat, ocontrol$Crow1positive, amat=Amat)
        Cmat = temp800$cmat
        Amat = temp800$amat
        eval(adjust.Dmat.expression)
        eval(td.expression)
    }


    if(ConstrainedQO)
    if(varlvI) {
        if(!length(xmat <- object@x)) stop("cannot obtain the model matrix")
        numat = xmat[,ocontrol$colx2.index,drop=FALSE] %*% Cmat
        sdnumat = sd(numat)
        Mmat = if(Rank > 1) diag(sdnumat) else matrix(sdnumat, 1, 1)
        Cmat = Cmat %*% solve(t(Mmat))
        Amat = Amat %*% Mmat
        temp800 = crow1C(Cmat, ocontrol$Crow1positive, amat=Amat)
        Cmat = temp800$cmat
        Amat = temp800$amat
        eval(adjust.Dmat.expression)
        eval(td.expression)
    }


    cx1i = ocontrol$colx1.index
    maximum = if(length(cx1i)==1 && names(cx1i)=="(Intercept)") {
        eta.temp = B1
        for(i in 1:M)
            eta.temp[i] = eta.temp[i] + 
                Amat[i,,drop=FALSE] %*% optimum[,i,drop=FALSE] +
                t(optimum[,i,drop=FALSE]) %*%
                Darray[,,i,drop= TRUE] %*% optimum[,i,drop=FALSE]
        mymax = object@family@inverse(rbind(eta.temp), extra=object@extra)  
        c(mymax)  # Convert from matrix to vector 
    } else {
        5 * rep(as.numeric(NA), len=M)  # Make "numeric"
    }
    names(maximum) = ynames
    
    lv.mat = if(ConstrainedQO) {
        object@x[,ocontrol$colx2.index,drop=FALSE] %*% Cmat 
    } else {
        object@lv
    }

    dimnames(Amat) = list(lp.names, lv.names)
    if(ConstrainedQO)
        dimnames(Cmat) = list(names(ocontrol$colx2.index), lv.names)
    if(!length(xmat <- object@x)) stop("cannot obtain the model matrix")
    dimnames(lv.mat) = list(dimnames(xmat)[[1]], lv.names)

    ans = 
    new(Class=if(ConstrainedQO) "Coef.qrrvglm" else "Coef.uqo",
         A=Amat, B1=B1, Constrained=ConstrainedQO, D=Darray,
         NOS = NOS, Rank = Rank,
         lv = lv.mat,
         lvOrder = lv.mat,
         Optimum=optimum, 
         OptimumOrder=optimum, 
         bellshaped=bellshaped,
         Dzero=dzero.vector,
         Maximum = maximum,
         Tolerance=Tolerance)
    if(ConstrainedQO) {ans@C = Cmat} else {Cmat = NULL}

    for(r in 1:Rank)
        ans@OptimumOrder[r,] = order(ans@Optimum[r,])
    for(r in 1:Rank)
        ans@lvOrder[,r] = order(ans@lv[,r])

    if(length(object@misc$estimated.dispersion) &&
       object@misc$estimated.dispersion) {
        p = length(object@coefficients)
        n = object@misc$n
        M = object@misc$M
        NOS = if(length(object@y)) ncol(object@y) else M
        pstar = if(ConstrainedQO) (p + length(Cmat)) else
                p + n*Rank # Adjustment; not sure about UQO 
        adjusted.dispersion = object@misc$dispersion * (n*M - p) /
                (n*M - pstar)
        ans@dispersion = adjusted.dispersion 
    }

    if(MSratio > 1) {
        keepIndex = seq(from=1, to=M, by=MSratio)
        ans@Dzero = ans@Dzero[keepIndex]
        ans@Optimum = ans@Optimum[,keepIndex,drop=FALSE]
        ans@Tolerance = ans@Tolerance[,,keepIndex,drop=FALSE]
        ans@bellshaped = ans@bellshaped[keepIndex]
        names(ans@Dzero) = ynames
    } else {
        dimnames(ans@D) = list(lv.names, lv.names, ynames)
    }
    names(ans@bellshaped) = ynames 
    dimnames(ans@Optimum) = list(lv.names, ynames)
    dimnames(ans@Tolerance) = list(lv.names, lv.names, ynames)
    ans 
}


setClass(Class="Coef.rrvglm", representation(
      "A"            = "matrix",
      "B1"           = "matrix",
      "C"            = "matrix",
      "Rank"         = "numeric",
      "colx1.index"  = "numeric",
      "colx2.index"  = "numeric",
      "Atilde"       = "matrix"))

setClass(Class="Coef.uqo", representation(
      "A"            = "matrix",
      "B1"           = "matrix",
      "Constrained"  = "logical",
      "D"            = "array",
      "NOS"          = "numeric",
      "Rank"         = "numeric",
      "lv"           = "matrix",
      "lvOrder"      = "matrix",
      "Maximum"      = "numeric",
      "Optimum"      = "matrix",
      "OptimumOrder" = "matrix",
      "bellshaped"   = "logical",
      "dispersion"   = "numeric",
      "Dzero"        = "logical",
      "Tolerance"    = "array"))

setClass(Class="Coef.qrrvglm", representation("Coef.uqo",
      "C"            = "matrix"))

printCoef.qrrvglm = function(x, ...) {

    object = x 
    Rank = object@Rank
    M = nrow(object@A)
    NOS = object@NOS
    iii = matrix(as.numeric(NA), NOS, Rank)
    if(Rank == 1) {  # || object@Diagonal
        for(i in 1:NOS) {
            fred = if(Rank>1) diag(object@Tolerance[,,i,drop=F]) else
                   object@Tolerance[,,i]
            if(all(fred > 0))
                iii[i,] = sqrt(fred)
        }
        dimnames(iii) = list(dimnames(object@Tolerance)[[3]],
                             if(Rank==1) "lv" else 
                             paste("Tolerance", dimnames(iii)[[2]], sep=""))
    } else {
        for(i in 1:NOS) {
            fred = eigen(object@Tolerance[,,i])
            if(all(fred$value > 0))
                iii[i,] = sqrt(fred$value)
        }
        dimnames(iii) = list(dimnames(object@Tolerance)[[3]],
                             paste("tol", 1:Rank, sep=""))
    }

    dimnames(object@A) = list(dimnames(object@A)[[1]],
        if(Rank > 1) paste("A", dimnames(object@A)[[2]], sep=".") else "A")

    Maximum = if(length(object@Maximum)) cbind(Maximum=object@Maximum) else NULL
    if(length(Maximum) && length(iii) && Rank==1)
        Maximum[is.na(iii),] = NA

    optmat = cbind(t(object@Optimum))
    dimnames(optmat) = list(dimnames(optmat)[[1]],
        if(Rank > 1) paste("Optimum", dimnames(optmat)[[2]], sep=".")
        else "Optimum")
    if(length(optmat) && length(iii) && Rank==1)
        optmat[is.na(iii),] = NA

    if( object@Constrained ) {
        cat("\nC matrix (constrained/canonical coefficients)\n")
        print(object@C, ...)
    }
    cat("\nB1 and A matrices\n")
    print(cbind(t(object@B1),
                A=object@A), ...)
    cat("\nOptima and maxima\n")
    print(cbind(Optimum=optmat,
                Maximum), ...)
    if(Rank > 1) { # !object@Diagonal && Rank > 1
        cat("\nTolerances\n") } else
        cat("\nTolerance\n")
    print(iii, ...)

    cat("\nStandard deviation of the latent variables (site scores)\n")
    print(sd(object@lv))
    invisible(object)
}


    setMethod("show", "Coef.qrrvglm", function(object)
        printCoef.qrrvglm(object))
    setMethod("print", "Coef.qrrvglm", function(x, ...)
        printCoef.qrrvglm(x, ...))
    setMethod("summary", "qrrvglm", function(object, ...)
        summary.qrrvglm(object, ...))

predict.qrrvglm <- function(object,
                         newdata=NULL,
                         type=c("link", "response", "lv", "terms"),
                         se.fit=FALSE,
                         deriv=0,
                         dispersion=NULL,
                         extra=object@extra, 
                         varlvI = FALSE, reference = NULL, ...)
{
    if(se.fit)
        stop("can't handle se.fit==TRUE yet")
    if(deriv != 0)
        stop("derivative is not equal to 0")

    if(mode(type) != "character" && mode(type) != "name")
        type <- as.character(substitute(type))
    type <- match.arg(type, c("link", "response", "lv", "terms"))[1]
    if(type=="lv")
        stop("can't handle type='lv' yet")
    if(type=="terms")
        stop("can't handle type='terms' yet")

    na.act = object@na.action
    object@na.action = list()
    M = object@misc$M
    Rank  = object@control$Rank

    if(length(newdata)) {
        p1 = length(object@control$colx1.index)
        Coefs = Coef(object, varlvI = varlvI, reference = reference)
        temptype = ifelse(type=="link", "response", type[1]) 

        conmat =  replace.constraints(vector("list", object@misc$p),
                      diag(M), object@control$colx1.index)
        conmat =  replace.constraints(conmat, Coefs@A,
                                      object@control$colx2.index)
        names(conmat) = object@misc$colnames.x
        object@constraints = conmat # To get RR-VGLM type eta

        newcoefs = lm2vlm.model.matrix(object@x, conmat, xij=object@control$xij)
        newcoefs = dimnames(newcoefs)[[2]]
        newcoefs = coef(object)[newcoefs] # Fixes up order of coefficients
        c.index = is.na(newcoefs) # Indices corresponding to C matrix 
        newcoefs[c.index] = t(Coefs@C) # Fill in the rest (C) of the outerprod

        object@coefficients = newcoefs 

        pvlm =  predict.vlm(object, newdata=newdata, type=temptype,
                            se.fit=se.fit, deriv=deriv,
                            dispersion=dispersion, extra=extra, ...)

        newcoefs[c.index] = 0  # Trick; set C==0 here to give B1 terms only
        object@coefficients = newcoefs
        pvlm1 =  if(!length(object@control$colx1.index)) 0 else 
                 predict.vlm(object, newdata=newdata, type=temptype,
                             se.fit=se.fit, deriv=deriv,
                             dispersion=dispersion, extra=extra, ...)

        lvmat = if(object@control$Corner)
            (pvlm - pvlm1)[,object@control$Index.corner,drop=FALSE] else
            stop("corner constraints needed") 

        for(j in 1:M)
            pvlm[,j] = pvlm[,j] + (if(Rank==1) (lvmat^2 * Coefs@D[,,j]) else
                (((lvmat %*% Coefs@D[,,j]) * lvmat) %*% rep(1, Rank)))
    } else {
        pvlm =  predict.vglm(object, type=type, se.fit=se.fit,
                             deriv=deriv, dispersion=dispersion,
                             extra=extra, ...)
    }

    pred = switch(type,
    response={ fv = if(length(newdata)) object@family@inverse(pvlm, extra) else
                    pvlm
        if(M > 1 && is.matrix(fv)) {
            dimnames(fv) <- list(dimnames(fv)[[1]],
                                 dimnames(object@fitted.values)[[2]])
        }
        fv
    },
    link = pvlm,
    terms=stop("failure here"))

    if(!length(newdata) && length(na.act)) {
        if(se.fit) {
            pred$fitted.values = napredict(na.act[[1]], pred$fitted.values)
            pred$se.fit = napredict(na.act[[1]], pred$se.fit)
        } else {
            pred = napredict(na.act[[1]], pred)
        }
    }
    pred
}

setMethod("predict", "qrrvglm", function(object, ...)
    predict.qrrvglm(object, ...))

coefqrrvglm = function(object, matrix.out = FALSE,
                        label = TRUE, compress = TRUE) {
    if(matrix.out)
        stop("currently can't handle matrix.out=TRUE")
    coefvlm(object, matrix.out = matrix.out, label = label, compress = compress)
}



residualsqrrvglm  <- function(object,
              type = c("deviance", "pearson", "working", "response", "ldot"),
              matrix.arg= TRUE) {
    stop("this function hasn't been written yet")

}

setMethod("residuals",  "qrrvglm", function(object, ...)
          residualsqrrvglm(object, ...))




printrrvglm <- function(x, ...)
{
    if(!is.null(cl <- x@call)) {
            cat("Call:\n")
            dput(cl)
    }
    coef <- x@coefficients
    if(any(nas <- is.na(coef))) {
            if(is.null(names(coef)))
                    names(coef) <- paste("b", 1:length(coef), sep = "")
            cat("\nCoefficients: (", sum(nas),
                    " not defined because of singularities)\n", sep = "")
    } else 
        cat("\nCoefficients:\n")

    if(FALSE) {
    Rank <- x@Rank
    if(!length(Rank))
        Rank <- sum(!nas)
    }

    if(FALSE) {
        nobs <- if(length(x@df.total)) x@df.total else length(x@residuals)
        rdf <- x@df.residual
        if(!length(rdf))
            rdf <- nobs - Rank
    }
    cat("\n")

    if(length(deviance(x)))
        cat("Residual Deviance:", format(deviance(x)), "\n")
    if(length(logLik(x)))
        cat("Log-likelihood:", format(logLik(x)), "\n")

    if(length(x@criterion)) {
        ncrit <- names(x@criterion)
        for(i in ncrit)
            if(i!="loglikelihood" && i!="deviance")
                cat(paste(i, ":", sep=""), format(x@criterion[[i]]), "\n")
    }

    invisible(x)
}




setMethod("print", "rrvglm", function(x, ...) printrrvglm(x, ...))

    setMethod("show", "rrvglm", function(object) printrrvglm(object))




rrvglm.control.Gaussian <- function(backchat= FALSE, half.stepsizing= FALSE,
                                    save.weight= TRUE, ...)
{

    list(backchat= FALSE, half.stepsizing= FALSE,
         save.weight=as.logical(save.weight)[1])
}



summary.rrvglm <- function(object, correlation= FALSE,
                           dispersion=NULL, digits=NULL, 
                           numerical= TRUE,
                           h.step = 0.0001, 
                           kill.all= FALSE, omit13= FALSE, fixA= FALSE, ...)
{





    if(!is.Numeric(h.step, allow=1) || abs(h.step)>1)
        stop("bad input for \"h.step\"")

    if(!object@control$Corner)
        stop("this function works with corner constraints only")

    if(is.null(dispersion))
        dispersion <- object@misc$dispersion

    newobject <- object
    class(newobject) <- "vglm"   # 6/2/02; For Splus6
    stuff <- summaryvglm(newobject, correlation=correlation,
                          dispersion=dispersion)

    answer <-
    new(Class="summary.rrvglm",
        object,
        call=stuff@call,
        coef3=stuff@coef3,
        cov.unscaled=stuff@cov.unscaled,
        correlation=stuff@correlation,
        df=stuff@df,
        pearson.resid=stuff@pearson.resid,
        sigma=stuff@sigma)


    if(is.numeric(stuff@dispersion))
        slot(answer, "dispersion") = stuff@dispersion



    tmp5 <- get.rrvglm.se1(object, omit13=omit13,
                           numerical=numerical, h.step=h.step,
                           kill.all=kill.all, fixA=fixA, ...) 
    if(any(diag(tmp5$cov.unscaled) <= 0) ||
       any(eigen(tmp5$cov.unscaled)$value <= 0)) {
        warning("cov.unscaled is not positive definite") 
    }

    answer@cov.unscaled <- tmp5$cov.unscaled 

    od <- if(is.numeric(object@misc$disper)) object@misc$disper else
        object@misc$default.disper
    if(is.numeric(dispersion)) {
        if(is.numeric(od) && dispersion!=od)
            warning("dispersion != object@misc$dispersion; using the former")
    } else {
        dispersion <- if(is.numeric(od)) od else 1
    }

    tmp8 = object@misc$M - object@control$Rank - 
           length(object@control$Structural.zero)
    answer@df[1] <- answer@df[1] + tmp8 * object@control$Rank
    answer@df[2] <- answer@df[2] - tmp8 * object@control$Rank
    if(dispersion==0) {
        dispersion <- tmp5$rss / answer@df[2]  # Estimate 
    }

    answer@coef3 <- get.rrvglm.se2(answer@cov.unscaled, dispersion=dispersion,
                                   coef=tmp5$coefficients)

    answer@dispersion <- dispersion
    answer@sigma <- dispersion^0.5


    answer
}





printsummary.rrvglm <- function(x, digits=NULL, quote= TRUE, prefix="")
{


    printsummary.vglm(x, digits = NULL, quote = TRUE, prefix = "")


    invisible(x)
}



get.rrvglm.se1 <- function(fit, omit13= FALSE, kill.all= FALSE,
                           numerical= TRUE,
                           fixA= FALSE, h.step=0.0001,
                           trace.arg= FALSE, ...) {




    if(length(fit@control$Nested) && fit@control$Nested)
        stop("sorry, can't handle nested models yet")

    Structural.zero = fit@control$Structural.zero


    if(!length(fit@x))
        stop("fix@x is empty. Run rrvglm(... , x= TRUE)")

    colx1.index = fit@control$colx1.index 
    colx2.index = fit@control$colx2.index 
    Blist <- fit@constraints
    ncolBlist <- unlist(lapply(Blist, ncol))

    p1 = length(colx1.index)
    p2 = length(colx2.index)

    Rank <- fit@control$Rank  # fit@misc$Nested.Rank   

    Amat <- fit@constraints[[colx2.index[1]]]
    Bmat <- if(p1) coef(fit, mat= TRUE)[colx1.index,,drop=FALSE] else NULL
    C.try <- coef(fit, mat= TRUE)[colx2.index,,drop=FALSE]
    Cmat <- C.try %*% Amat %*% solve(t(Amat) %*% Amat)

    x1mat <- if(p1) fit@x[,colx1.index,drop=FALSE] else NULL
    x2mat <- fit@x[,colx2.index,drop=FALSE]
 
    wz <- weights(fit, type="w")  # old: wweights(fit)  #fit@weights
    if(!length(wz))
        stop("can't get fit@weights")

    M <- fit@misc$M
    n <- fit@misc$n
    Index.corner <- fit@control$Index.corner   # used to be (1:Rank);
    zmat <- fit@predictors + fit@residuals
    theta <- c(Amat[-c(Index.corner,Structural.zero),])
    if(fit@control$checkwz)
        wz = checkwz(wz, M=M, trace=trace, wzeps=fit@control$wzepsilon)
    U <- vchol(wz, M=M, n=n, silent= TRUE)

    if(numerical) {
        delct.da <- num.deriv.rrr(fit, M=M, r=Rank,
                                  x1mat=x1mat, x2mat=x2mat, p2=p2, 
                                  Index.corner, Aimat=Amat, Bmat=Bmat, Cimat=Cmat,
                                  h.step=h.step, colx2.index=colx2.index,
                                  xij=fit@control$xij,
                                  Structural.zero=Structural.zero)
    } else {
        delct.da <- dctda.fast.only(theta=theta, wz=wz, U=U, zmat, M=M, r=Rank,
                                    x1mat=x1mat, x2mat=x2mat,
                                    p2=p2, Index.corner, Aimat=Amat,
                                    Bmat=Bmat, Cimat=Cmat,
                                    xij=fit@control$xij,
                                    Structural.zero=Structural.zero)
    }


    newobject <- fit
    class(newobject) <- "vglm"   # 6/2/02; For Splus6
    sfit2233 <- summaryvglm(newobject) 
    d8 <-  dimnames(sfit2233@cov.unscaled)[[1]]
    cov2233 <- solve(sfit2233@cov.unscaled) # Includes any intercepts
    dimnames(cov2233) = list(d8, d8)

    log.vec33 = NULL 
    nassign = names(fit@constraints) 
    choose.from =  varassign(fit@constraints, nassign)
    for(ii in nassign)
        if(any(ii== names(colx2.index))) {
            log.vec33 = c(log.vec33, choose.from[[ii]])
        }
    cov33 = cov2233[ log.vec33, log.vec33, drop=FALSE]   # r*p2 by r*p2
    cov23 = cov2233[-log.vec33, log.vec33, drop=FALSE]
    cov22 = cov2233[-log.vec33,-log.vec33, drop=FALSE]


    lv.mat <- x2mat %*% Cmat
    offs = matrix(0, n, M)     # The "0" handles Structural.zero's 
    offs[,Index.corner] = lv.mat
    if(M == (Rank+length(Structural.zero)))
        stop("can't handle full-rank models yet")
    cm = matrix(0, M, M-Rank-length(Structural.zero))
    cm[-c(Index.corner,Structural.zero),] = diag(M-Rank-length(Structural.zero))

    Blist = vector("list", length(colx1.index)+1) 
    names(Blist) = c(names(colx1.index), "I(lv.mat)")
    for(ii in names(colx1.index))
        Blist[[ii]]  = fit@constraints[[ii]]
    Blist[["I(lv.mat)"]] = cm


    if(p1) {
        ooo = fit@assign
        bb = NULL 
        for(ii in 1:length(ooo)) {
            if(any(ooo[[ii]][1] == colx1.index))
                bb = c(bb, names(ooo)[ii])
        }

        has.intercept = any(bb=="(Intercept)")
        bb[bb=="(Intercept)"] = "1"
        if(p1>1)
            bb = paste(bb, collapse="+")
        if(has.intercept) {
            bb = paste("zmat - offs ~ ", bb, " + I(lv.mat)", collapse=" ")
        } else {
            bb = paste("zmat - offs ~ -1 + ", bb, " + I(lv.mat)", collapse=" ")
        }
        bb = as.formula(bb)
    } else {
        bb = as.formula("zmat - offs ~ -1 + I(lv.mat)")
    }


    if(( is.R() && fit@misc$dataname == "list") ||
       (!is.R() && fit@misc$dataname == "sys.parent")) {
        dspec = FALSE
    } else {
        if(is.R()) {
            mytext1 = "exists(x=fit@misc$dataname, envir = VGAMenv)"
            myexp1 = parse(text=mytext1)
            is.there = eval(myexp1)
            bbdata= if(is.there) get(fit@misc$dataname, envir=VGAMenv) else
                    get(fit@misc$dataname)
        } else {
            bbdata = get(fit@misc$dataname)
        }
        dspec = TRUE
    }

    if(!is.R()) {

        stop("26-9-2007: uncomment out the following lines to run it in Splus")
    }

    fit1122 <- if(dspec) vlm(bb,
                  constraint=Blist, crit="d", weight=wz, data=bbdata, 
                  save.weight= TRUE, smart= FALSE, trace=trace.arg, x= TRUE) else 
              vlm(bb,
                  constraint=Blist, crit="d", weight=wz,
                  save.weight= TRUE, smart= FALSE, trace=trace.arg, x= TRUE)



    sfit1122 <- summaryvlm(fit1122)
    d8 <-  dimnames(sfit1122@cov.unscaled)[[1]]
    cov1122 <- solve(sfit1122@cov.unscaled)
    dimnames(cov1122) = list(d8, d8)

    lcs = length(coef(sfit1122))
    log.vec11 = (lcs-(M-Rank-length(Structural.zero))*Rank+1):lcs
    cov11 = cov1122[log.vec11,  log.vec11, drop=FALSE]
    cov12 = cov1122[ log.vec11, -log.vec11, drop=FALSE]
    cov22 = cov1122[-log.vec11, -log.vec11, drop=FALSE]
    cov13 = delct.da %*% cov33    # zz; this always seems to be negative 


    if(omit13) 
        cov13 = cov13 * 0   # zero it

    if(kill.all) {
        cov13 = cov13 * 0   # zero it
        if(fixA) {
            cov12 = cov12 * 0   # zero it
        } else {
            cov23 = cov23 * 0   # zero it
        }
    }

 cov13 = -cov13   # Richards (1961)

    if(fixA) {
        cov.unscaled <- rbind(cbind(cov1122, rbind(cov13, cov23)),
                              cbind(t(cov13), t(cov23), cov33))
    } else {
        cov.unscaled <- rbind(cbind(cov11, cov12, cov13),
                              cbind(rbind(t(cov12), t(cov13)), cov2233))
    }

    ans <- solve(cov.unscaled)

    # Get all the coefficients 
    acoefs <- c(fit1122@coefficients[log.vec11], fit@coefficients)
    dimnames(ans) = list(names(acoefs), names(acoefs))
    list(cov.unscaled=ans, coefficients=acoefs, rss=sfit1122@rss)
}



get.rrvglm.se2 <- function(cov.unscaled, dispersion=1, coefficients) {

    d8 <-  dimnames(cov.unscaled)[[1]]
    ans <- matrix(coefficients, length(coefficients), 3) 
    ans[,2] <- sqrt(dispersion) * sqrt(diag(cov.unscaled))
    ans[,3] <- ans[,1] / ans[,2]
    dimnames(ans) <- list(d8, c("Value", "Std. Error", "t value"))
    ans
}



num.deriv.rrr <- function(fit, M, r, x1mat, x2mat,
                          p2, Index.corner, Aimat, Bmat, Cimat, 
                          h.step=0.0001, colx2.index,
                          xij=NULL, Structural.zero=NULL)
{

    nn <- nrow(x2mat)
    if(nrow(Cimat)!=p2 || ncol(Cimat)!=r)
        stop("Cimat wrong shape")

    dct.da <- matrix(as.numeric(NA), (M-r-length(Structural.zero))*r, r*p2)

    if((length(Index.corner) + length(Structural.zero)) == M)
        stop("can't handle full rank models yet")
    cbindex = (1:M)[-c(Index.corner, Structural.zero)]

    ptr = 1
    for(s in 1:r)
        for(tt in cbindex) {
            small.Blist = vector("list", p2)
            pAmat = Aimat
            pAmat[tt,s] = pAmat[tt,s] + h.step   # Perturb it
            for(ii in 1:p2)
                small.Blist[[ii]] = pAmat

            offset = if(length(fit@offset)) fit@offset else 0
            if(all(offset==0)) offset = 0
            neweta = x1mat %*% Bmat + x2mat %*% Cimat %*% t(pAmat)
            fit@predictors = neweta


            newmu <- fit@family@inverse(neweta, fit@extra) 
            fit@fitted.values = newmu

            fred = weights(fit, type="w", deriv= TRUE, ignore.slot= TRUE)
            if(!length(fred))
                stop("can't get @weights and $deriv from object")
            wz = fred$weights
            deriv.mu <- fred$deriv

            U <- vchol(wz, M=M, n=nn, silent= TRUE)
            tvfor <- vforsub(U, as.matrix(deriv.mu), M=M, n=nn)
            newzmat <- neweta + vbacksub(U, tvfor, M=M, n=nn) - offset

            newfit = vlm.wfit(x=x2mat, z=newzmat - x1mat %*% Bmat,
                              Blist=small.Blist, U = U, 
                              matrix.out = FALSE, XBIG = FALSE,
                              rss = TRUE, qr = FALSE, x.ret = FALSE, offset = NULL,
                              xij=xij)
            dct.da[ptr,] <- (newfit$coef - t(Cimat)) / h.step
            ptr = ptr + 1
        }

    dct.da
}




dctda.fast.only <- function(theta, wz, U, zmat, M, r, x1mat, x2mat,
                            p2, Index.corner, Aimat, Bmat, Cimat,
                            xij=NULL,
                            Structural.zero=NULL)
{


    if(length(Structural.zero))
        stop("can't handle Structural.zero in dctda.fast.only()")

    nn <- nrow(x2mat)
    if(nrow(Cimat)!=p2 || ncol(Cimat)!=r)
        stop("Cimat wrong shape")

    fred <- kronecker(matrix(1,1,r), x2mat)
    fred <- kronecker(fred, matrix(1,M,1))
    barney <- kronecker(Aimat, matrix(1,1,p2))
    barney <- kronecker(matrix(1,nn,1), barney)

    temp <- array(t(barney*fred), c(p2*r, M, nn))
    temp <- aperm(temp, c(2,1,3))     # M by p2*r by nn
    temp <- mux5(wz, temp, M=M, matrix.arg= TRUE)
    temp <- m2adefault(temp, M=p2*r)         # Note M != M here!
    G <- solve(apply(temp,1:2,sum))   # p2*r by p2*r 

    dc.da <- array(NA, c(p2, r, M, r))  # different from other functions
    if(length(Index.corner) == M)
        stop("can't handle full rank models yet")
    cbindex <- (1:M)[-Index.corner]    # complement of Index.corner 
    resid2 <- if(length(x1mat))
        mux22(t(wz), zmat - x1mat %*% Bmat, M=M, upper= FALSE, as.mat= TRUE) else 
        mux22(t(wz), zmat                 , M=M, upper= FALSE, as.mat= TRUE)

    for(s in 1:r)
        for(tt in cbindex) {
            fred <- t(x2mat) * matrix(resid2[,tt], p2, nn, byrow= TRUE)  # p2 * nn
            temp2 <- kronecker(ei(s,r), apply(fred,1,sum))
            for(k in 1:r) {
                Wiak <- mux22(t(wz), matrix(Aimat[,k], nn, M, byrow= TRUE), 
                              M=M, upper= FALSE, as.mat= TRUE)  # nn * M
                wxx <- Wiak[,tt] * x2mat
                blocki <- t(x2mat) %*% wxx 
                temp4a <- blocki %*% Cimat[,k]
                if(k==1) {
                    temp4b <- blocki %*% Cimat[,s]
                }
                temp2 = temp2 - kronecker(ei(s,r), temp4a) -
                                kronecker(ei(k,r), temp4b)
            }
            dc.da[,,tt,s] <- G %*% temp2 
        }
    ans1 <- dc.da[,,cbindex,,drop=FALSE]  # p2 x r x (M-r) x r 
    ans1 <- aperm(ans1, c(2,1,3,4))   # r x p2 x (M-r) x r 

    ans1 <- matrix(c(ans1), r*p2, (M-r)*r)
    ans1 <- t(ans1)
    ans1
}



dcda.fast <- function(theta, wz, U, z, M, r, xmat, pp, Index.corner,
                      intercept= TRUE, xij=NULL)
{



    nn <- nrow(xmat)

    Aimat <- matrix(as.numeric(NA), M, r)
    Aimat[Index.corner,] <- diag(r)
    Aimat[-Index.corner,] <- theta    # [-(1:M)]

    if(intercept) {
        Blist <- vector("list", pp+1)
        Blist[[1]] <- diag(M)
        for(i in 2:(pp+1))
            Blist[[i]] <- Aimat
    } else {
        Blist <- vector("list", pp)
        for(i in 1:(pp))
            Blist[[i]] <- Aimat
    }

    coeffs <- vlm.wfit(xmat, z, Blist, U=U, matrix.out= TRUE,
                        xij=xij)$mat.coef
    c3 <- coeffs <- t(coeffs)  # transpose to make M x (pp+1)


    int.vec <- if(intercept) c3[,1] else 0  # \boldeta_0
    Cimat <- if(intercept) t(c3[Index.corner,-1,drop=FALSE]) else 
             t(c3[Index.corner,,drop=FALSE])
    if(nrow(Cimat)!=pp || ncol(Cimat)!=r)
        stop("Cimat wrong shape")

    fred <- kronecker(matrix(1,1,r), if(intercept) xmat[,-1,drop=FALSE] else xmat)
    fred <- kronecker(fred, matrix(1,M,1))
    barney <- kronecker(Aimat, matrix(1,1,pp))
    barney <- kronecker(matrix(1,nn,1), barney)

    temp <- array(t(barney*fred), c(r*pp,M,nn))
    temp <- aperm(temp, c(2,1,3))
    temp <- mux5(wz, temp, M=M, matrix.arg= TRUE)
    temp <- m2adefault(temp, M=r*pp)     # Note M != M here!
    G <- solve(apply(temp,1:2,sum))

    dc.da <- array(NA, c(pp,r,M,r))  # different from other functions
    cbindex <- (1:M)[-Index.corner]
    resid2 <- mux22(t(wz), z - matrix(int.vec, nn, M, byrow= TRUE), M=M,
                    upper= FALSE, as.mat= TRUE)  # mat= TRUE,

    for(s in 1:r)
        for(tt in cbindex) {
            fred <- (if(intercept) t(xmat[,-1,drop=FALSE]) else
                     t(xmat)) * matrix(resid2[,tt],pp,nn,byrow= TRUE) 
            temp2 <- kronecker(ei(s,r), apply(fred,1,sum))

            temp4 <- rep(0,pp)
            for(k in 1:r) {
                Wiak <- mux22(t(wz), matrix(Aimat[,k],nn,M,byrow= TRUE), 
                              M=M, upper= FALSE, as.mat= TRUE)  # mat= TRUE, 
                wxx <- Wiak[,tt] * (if(intercept) xmat[,-1,drop=FALSE] else xmat)
                blocki <- (if(intercept) t(xmat[,-1,drop=FALSE]) else t(xmat)) %*% wxx 
                temp4 <- temp4 + blocki %*% Cimat[,k]
            }
            dc.da[,,tt,s] <- G %*% (temp2 - 2 * kronecker(ei(s,r),temp4))
        }
    ans1 <- dc.da[,,cbindex,,drop=FALSE]  # pp x r x (M-r) x r 
    ans1 <- aperm(ans1, c(2,1,3,4))   # r x pp x (M-r) x r 

    ans1 <- matrix(c(ans1), (M-r)*r, r*pp, byrow= TRUE)


    detastar.da <- array(0,c(M,r,r,nn))
    for(s in 1:r)
        for(j in 1:r) {
            t1 <- t(dc.da[,j,,s])
            t1 <- matrix(t1, M, pp)
            detastar.da[,j,s,] <- t1 %*% (if(intercept)
                                  t(xmat[,-1,drop=FALSE]) else t(xmat))
        }

    etastar <- (if(intercept) xmat[,-1,drop=FALSE] else xmat) %*% Cimat
    eta <- matrix(int.vec, nn, M, byrow= TRUE) + etastar %*% t(Aimat)

    sumWinv <- solve((m2adefault(t(apply(wz, 2, sum)), M=M))[,,1])

    deta0.da <- array(0,c(M,M,r))
    AtWi <- kronecker(matrix(1,nn,1), Aimat)
    AtWi <- mux111(t(wz), AtWi, M=M, upper= FALSE)  # matrix.arg= TRUE, 
    AtWi <- array(t(AtWi), c(r,M,nn))
    for(ss in 1:r) {
        temp90 <- (m2adefault(t(apply(etastar[,ss]*wz,2,sum)), M=M))[,,1] #MxM
        temp92 <- array(detastar.da[,,ss,], c(M,r,nn))
        temp93 <- mux7(temp92, AtWi)
        temp91 <- apply(temp93, 1:2, sum)   # M x M
        deta0.da[,,ss] <- -(temp90 + temp91) %*% sumWinv
    }
    ans2 <- deta0.da[-(1:r),,,drop=FALSE]   # (M-r) x M x r
    ans2 <- aperm(ans2, c(1,3,2))       # (M-r) x r x M
    ans2 <- matrix(c(ans2), (M-r)*r, M) 

    list(dc.da=ans1, dint.da=ans2)
}



rrr.deriv.rss <- function(theta, wz, U, z, M, r, xmat,
                          pp, Index.corner, intercept= TRUE,
                          xij=NULL)
{

    Amat <- matrix(as.numeric(NA), M, r)
    Amat[Index.corner,] <- diag(r)
    Amat[-Index.corner,] <- theta    # [-(1:M)]

    if(intercept) {
        Blist <- vector("list", pp+1)
        Blist[[1]] <- diag(M)
        for(i in 2:(pp+1))
            Blist[[i]] <- Amat
    } else {
        Blist <- vector("list", pp)
        for(i in 1:(pp))
            Blist[[i]] <- Amat
    }

    vlm.wfit(xmat, z, Blist, U=U, matrix.out= FALSE, rss= TRUE, xij=xij)$rss
}




rrr.deriv.gradient.fast <- function(theta, wz, U, z, M, r, xmat,
                                    pp, Index.corner, intercept= TRUE)
{




    nn <- nrow(xmat)

    Aimat <- matrix(as.numeric(NA), M, r)
    Aimat[Index.corner,] <- diag(r)
    Aimat[-Index.corner,] <- theta    # [-(1:M)]

    if(intercept) {
        Blist <- vector("list", pp+1)
        Blist[[1]] <- diag(M)
        for(i in 2:(pp+1))
            Blist[[i]] <- Aimat
    } else {
        Blist <- vector("list", pp)
        for(i in 1:(pp))
            Blist[[i]] <- Aimat
    }

    coeffs <- vlm.wfit(xmat, z, Blist, U=U, matrix.out= TRUE,
                       xij=NULL)$mat.coef
    c3 <- coeffs <- t(coeffs)  # transpose to make M x (pp+1)


    int.vec <- if(intercept) c3[,1] else 0  # \boldeta_0
    Cimat <- if(intercept) t(c3[Index.corner,-1,drop=FALSE]) else
             t(c3[Index.corner,,drop=FALSE])
    if(nrow(Cimat)!=pp || ncol(Cimat)!=r)
        stop("Cimat wrong shape")

    fred = kronecker(matrix(1,1,r), if(intercept) xmat[,-1,drop=FALSE] else xmat)
    fred <- kronecker(fred, matrix(1,M,1))
    barney <- kronecker(Aimat, matrix(1,1,pp))
    barney <- kronecker(matrix(1,nn,1), barney)

    temp <- array(t(barney*fred), c(r*pp,M,nn))
    temp <- aperm(temp, c(2,1,3))
    temp <- mux5(wz, temp, M=M, matrix.arg= TRUE)
    temp <- m2adefault(temp, M=r*pp)     # Note M != M here!
    G <- solve(apply(temp,1:2,sum))

    dc.da <- array(NA,c(pp,r,r,M))
    cbindex <- (1:M)[-Index.corner]
    resid2 <- mux22(t(wz), z - matrix(int.vec,nn,M,byrow= TRUE), M=M,
                    upper= FALSE, as.mat= TRUE)  # mat= TRUE,

    for(s in 1:r)
        for(tt in cbindex) {
            fred <- (if(intercept) t(xmat[,-1,drop=FALSE]) else
                     t(xmat)) * matrix(resid2[,tt],pp,nn,byrow= TRUE) 
            temp2 <- kronecker(ei(s,r), apply(fred,1,sum))

            temp4 <- rep(0,pp)
            for(k in 1:r) {
                Wiak <- mux22(t(wz), matrix(Aimat[,k],nn,M,byrow= TRUE), 
                              M=M, upper= FALSE, as.mat= TRUE)  # mat= TRUE, 
                wxx <- Wiak[,tt] * (if(intercept) xmat[,-1,drop=FALSE] else xmat)
                blocki <- (if(intercept) t(xmat[,-1,drop=FALSE]) else t(xmat)) %*% wxx 
                temp4 <- temp4 + blocki %*% Cimat[,k]
            }
            dc.da[,,s,tt] <- G %*% (temp2 - 2 * kronecker(ei(s,r),temp4))
        }

    detastar.da <- array(0,c(M,r,r,nn))
    for(s in 1:r)
        for(j in 1:r) {
            t1 <- t(dc.da[,j,s,])
            t1 <- matrix(t1, M, pp)
            detastar.da[,j,s,] <- t1 %*% (if(intercept)
                                  t(xmat[,-1,drop=FALSE]) else t(xmat))
        }

    etastar <- (if(intercept) xmat[,-1,drop=FALSE] else xmat) %*% Cimat
    eta <- matrix(int.vec, nn, M, byrow= TRUE) + etastar %*% t(Aimat)

    sumWinv <- solve((m2adefault(t(apply(wz, 2, sum)), M=M))[,,1])

    deta0.da <- array(0,c(M,M,r))

    AtWi <- kronecker(matrix(1,nn,1), Aimat)
    AtWi <- mux111(t(wz), AtWi, M=M, upper= FALSE)  # matrix.arg= TRUE, 
    AtWi <- array(t(AtWi), c(r,M,nn))

    for(ss in 1:r) {
        temp90 <- (m2adefault(t(apply(etastar[,ss]*wz,2,sum)), M=M))[,,1]   # M x M
        temp92 <- array(detastar.da[,,ss,],c(M,r,nn))
        temp93 <- mux7(temp92,AtWi)
        temp91 <- apply(temp93,1:2,sum)   # M x M
        deta0.da[,,ss] <- -(temp90 + temp91) %*% sumWinv
    }

    ans <- matrix(0,M,r)
    fred <- mux22(t(wz), z-eta, M=M, upper= FALSE, as.mat= TRUE) # mat= TRUE, 
    fred.array <- array(t(fred %*% Aimat),c(r,1,nn))
    for(s in 1:r) {
        a1 <- apply(fred %*% t(deta0.da[,,s]),2,sum)
        a2 <- apply(fred * etastar[,s],2,sum)
        temp92 <- array(detastar.da[,,s,],c(M,r,nn))
        temp93 <- mux7(temp92, fred.array)
        a3 <- apply(temp93,1:2,sum)
        ans[,s] <- a1 + a2 + a3
    }

    ans <- -2 * c(ans[cbindex,])

    ans
}




vellipse = function(R, ratio=1, orientation=0, center=c(0,0), N=300) {
    if(length(center) != 2) stop("center must be of length 2")
    theta =       2*pi*(0:N)/N
    x1 =       R*cos(theta)
    y1 = ratio*R*sin(theta)
    x = center[1] + cos(orientation)*x1 - sin(orientation)*y1
    y = center[2] + sin(orientation)*x1 + cos(orientation)*y1
    cbind(x, y)
}


biplot.qrrvglm = function(x, ...) {
    stop("biplot.qrrvglm has been replaced by the function lvplot.qrrvglm")
}


lvplot.qrrvglm = function(object, varlvI = FALSE, reference = NULL,
          add= FALSE, plot.it= TRUE, rug= TRUE, y = FALSE, 
          type=c("fitted.values", "predictors"),
          xlab=paste("Latent Variable", if(Rank==1) "" else " 1", sep=""),
          ylab=if(Rank==1) switch(type, predictors="Predictors", 
              fitted.values="Fitted values") else "Latent Variable 2",
          pcex=par()$cex, pcol=par()$col, pch=par()$pch, 
          llty=par()$lty, lcol=par()$col, llwd=par()$lwd,
          label.arg= FALSE, adj.arg=-0.1, 
          ellipse=0.95, Absolute= FALSE, 
              elty=par()$lty, ecol=par()$col, elwd=par()$lwd, egrid=200,
          chull.arg= FALSE, clty=2, ccol=par()$col, clwd=par()$lwd,
              cpch = "   ",
          C = FALSE,
              OriginC = c("origin","mean"),
              Clty=par()$lty, Ccol=par()$col, Clwd=par()$lwd,
              Ccex=par()$cex, Cadj.arg=-0.1, stretchC=1, 
          sites= FALSE, spch=NULL, scol=par()$col, scex=par()$cex,
          sfont=par()$font,
          check.ok = TRUE, ...)
{
    if(mode(type) != "character" && mode(type) != "name")
        type <- as.character(substitute(type))
    type <- match.arg(type, c("fitted.values", "predictors"))[1]

    if(is.numeric(OriginC)) OriginC = rep(OriginC, len=2) else {
        if(mode(OriginC) != "character" && mode(OriginC) != "name")
            OriginC <- as.character(substitute(OriginC))
        OriginC <- match.arg(OriginC, c("origin","mean"))[1]
    }

    if(length(ellipse) > 1) stop("ellipse must be of length 1 or 0")
    if(is.logical(ellipse)) {ellipse = if(ellipse) 0.95 else NULL}

    Rank <- object@control$Rank
    if(Rank > 2)
        stop("can only handle rank 1 or 2 models")
    M = object@misc$M
    NOS = ncol(object@y)
    MSratio = M / NOS  # First value is g(mean) = quadratic form in lv
    n = object@misc$n
    colx2.index = object@control$colx2.index
    cx1i = object@control$colx1.index
    if(check.ok)
        if(!(length(cx1i)==1 && names(cx1i)=="(Intercept)"))
            stop("latent variable plots allowable only for Norrr = ~ 1 models")

    Coef.list = Coef(object, varlvI = varlvI, reference = reference)
    if( C) Cmat = Coef.list@C
    nustar = Coef.list@lv # n x Rank 

    if(!plot.it) return(nustar)

    r.curves = slot(object, type)   # n times M (\boldeta or \boldmu) 
    if(!add) {
        if(Rank==1) {
            matplot(nustar,
                    if( y && type=="fitted.values") object@y else r.curves,
                    type="n", xlab=xlab, ylab=ylab, ...)
        } else { # Rank==2
            matplot(c(Coef.list@Optimum[1,], nustar[,1]),
                    c(Coef.list@Optimum[2,], nustar[,2]),
                    type="n", xlab=xlab, ylab=ylab, ...)
        }
    }

    if((length(pch)  != 1 && length(pch)  != ncol(r.curves)) ||
       (length(pcol) != 1 && length(pcol) != ncol(r.curves)) ||
       (length(pcex) != 1 && length(pcex) != ncol(r.curves)))
        stop("pch, pcol and pcex must be of length 1 or ncol(r.curves)")

    pch  <- rep(pch,  leng=ncol(r.curves))
    pcol <- rep(pcol, leng=ncol(r.curves))
    pcex <- rep(pcex, leng=ncol(r.curves))
    llty <- rep(llty, leng=ncol(r.curves))
    lcol <- rep(lcol, leng=ncol(r.curves))
    llwd <- rep(llwd, leng=ncol(r.curves))
    elty <- rep(elty, leng=ncol(r.curves))
    ecol <- rep(ecol, leng=ncol(r.curves))
    elwd <- rep(elwd, leng=ncol(r.curves))
    adj.arg <- rep(adj.arg, leng=ncol(r.curves))
    if( C ) {
        Clwd <- rep(Clwd, leng=nrow(Cmat))
        Clty <- rep(Clty, leng=nrow(Cmat))
        Ccol <- rep(Ccol, leng=nrow(Cmat))
        Cadj.arg <- rep(Cadj.arg, leng=nrow(Cmat))
        Ccex <- rep(Ccex, leng=nrow(Cmat))
    }

    if(Rank==1) {
        for(i in 1:ncol(r.curves)) {
            xx = nustar 
            yy = r.curves[,i]
            o = sort.list(xx)
            xx = xx[o]
            yy = yy[o]
            lines(xx, yy, col=lcol[i], lwd=llwd[i], lty=llty[i])
            if( y && type=="fitted.values") {
                ypts = object@y
                if(ncol(as.matrix(ypts)) == ncol(r.curves))
                    points(xx, ypts[o,i], col=pcol[i], cex=pcex[i], pch=pch[i])
            } 
        } 
        if(rug) rug(xx) 
    } else {
        for(i in 1:ncol(r.curves))
            points(Coef.list@Optimum[1,i], Coef.list@Optimum[2,i],
                   col=pcol[i], cex=pcex[i], pch=pch[i])
        if(label.arg) {
            for(i in 1:ncol(r.curves))
                text(Coef.list@Optimum[1,i], Coef.list@Optimum[2,i],
                     labels=(dimnames(Coef.list@Optimum)[[2]])[i], 
                     adj=adj.arg[i], col=pcol[i], cex=pcex[i])
        }
        if(chull.arg) {
            hull = chull(nustar[,1], nustar[,2])
            hull = c(hull, hull[1])
            lines(nustar[hull,1], nustar[hull,2], type="b", pch=cpch,
                  lty=clty, col=ccol, lwd=clwd)
        }
        if(length(ellipse)) {
            ellipse.temp = if(ellipse > 0) ellipse else 0.95
            if(ellipse < 0 && (!object@control$EqualTolerances || varlvI))
                stop(paste("an equal-tolerances assumption and varlvI=FALSE",
                     "is needed for \"ellipse\" < 0"))
            if( check.ok ) {
                colx1.index = object@control$colx1.index
                if(!(length(colx1.index)==1 &&
                     names(colx1.index)=="(Intercept)"))
                     stop("can only plot ellipses for intercept models only")
            }
            for(i in 1:ncol(r.curves)) {
                cutpoint = object@family@link(if(Absolute) ellipse.temp
                                else Coef.list@Maximum[i] * ellipse.temp,
                                extra=object@extra)
                if(MSratio > 1) 
                    cutpoint = cutpoint[1,1]

                cutpoint = object@family@link(Coef.list@Maximum[i],
                               extra=object@extra) - cutpoint
                if(is.finite(cutpoint) && cutpoint > 0) {
                    Mmat = diag(rep(ifelse(object@control$Crow1positive, 1, -1),
                                    len=Rank))
                    etoli = eigen(t(Mmat) %*% Coef.list@Tolerance[,,i] %*% Mmat)
                    A=ifelse(etoli$val[1]>0,sqrt(2*cutpoint*etoli$val[1]),Inf)
                    B=ifelse(etoli$val[2]>0,sqrt(2*cutpoint*etoli$val[2]),Inf)
                    if(ellipse < 0) A = B = -ellipse / 2

                    theta.angle = asin(etoli$vector[2,1]) *
                        ifelse(object@control$Crow1positive[2], 1, -1)
                    if(object@control$Crow1positive[1])
                        theta.angle = pi - theta.angle
                    if(all(is.finite(c(A,B))))
                        lines(vellipse(R=2*A, ratio=B/A, orient=theta.angle,
                                       center=Coef.list@Optimum[,i], N=egrid),
                              lwd=elwd[i], col=ecol[i], lty=elty[i])
                }
            }
        }

        if( C ) {
            if(is.character(OriginC) && OriginC=="mean")
                OriginC = c(mean(nustar[,1]), mean(nustar[,2]))
            if(is.character(OriginC) && OriginC=="origin")
                OriginC = c(0,0)
            for(i in 1:nrow(Cmat))
                arrows(x0=OriginC[1], y0=OriginC[2],
                       x1=OriginC[1] + stretchC*Cmat[i,1],
                       y1=OriginC[2] + stretchC*Cmat[i,2],
                       lty=Clty[i], col=Ccol[i], lwd=Clwd[i])
            if(label.arg) {
                temp200 = dimnames(Cmat)[[1]]
                for(i in 1:nrow(Cmat))
                    text(OriginC[1] + stretchC*Cmat[i,1],
                         OriginC[2] + stretchC*Cmat[i,2], col=Ccol[i],
                         labels=temp200[i], adj=Cadj.arg[i], cex=Ccex[i])
            }
        }
        if(sites) {
            text(nustar[,1], nustar[,2], adj=0.5,
                 labels=if(is.null(spch)) dimnames(nustar)[[1]] else 
                 rep(spch, length=nrow(nustar)), col=scol, cex=scex, font=sfont)
        }
    }
    invisible(nustar)
}



lvplot.rrvglm = function(object,
                         A=TRUE,
                         C=TRUE,
                         scores=FALSE, plot.it= TRUE,
                         groups=rep(1,n),
                         gapC=sqrt(sum(par()$cxy^2)), scaleA=1,
                         xlab="Latent Variable 1",
                         ylab="Latent Variable 2",
                         Alabels= if(length(object@misc$predictors.names))
                   object@misc$predictors.names else paste("LP", 1:M, sep=""),
                         Aadj=par()$adj,
                         Acex=par()$cex,
                         Acol=par()$col,
                         Apch=NULL,
                         Clabels=dimnames(Cmat)[[1]],
                         Cadj=par()$adj,
                         Ccex=par()$cex,
                         Ccol=par()$col, 
                         Clty=par()$lty, 
                         Clwd=par()$lwd, 
                         chull.arg=FALSE,
                         ccex=par()$cex,
                         ccol=par()$col,
                         clty=par()$lty,
                         clwd=par()$lwd,
                         spch=NULL,
                         scex=par()$cex,
                         scol=par()$col,
                         slabels=dimnames(x2mat)[[1]],
                         ...)
{


    if(object@control$Rank != 2 && plot.it)
        stop("can only handle rank-2 models")
    M = object@misc$M
    n = object@misc$n
    colx1.index = object@control$colx1.index
    colx2.index = object@control$colx2.index
    p1 = length(colx1.index)
    Coef.list = Coef(object)
    Amat = Coef.list@A
    Cmat = Coef.list@C

    Amat = Amat * scaleA
    dimnames(Amat) = list(object@misc$predictors.names, NULL) 
    Cmat = Cmat / scaleA

    if(!length(object@x)) {
        object@x = model.matrixvlm(object, type="lm")
    }
    x2mat = object@x[,colx2.index,drop=FALSE]
    nuhat = x2mat %*% Cmat
    if(!plot.it) return(as.matrix(nuhat))

    index.nosz = 1:M   # index of no structural zeros; zz  
    allmat = rbind(if(A) Amat else NULL, 
                   if(C) Cmat else NULL, 
                   if(scores) nuhat else NULL)

    plot(allmat[,1], allmat[,2], type="n",
         xlab=xlab, ylab=ylab, ...) # xlim etc. supplied through ...

    if(A) {
        Aadj = rep(Aadj, len=length(index.nosz))
        Acex = rep(Acex, len=length(index.nosz))
        Acol = rep(Acol, len=length(index.nosz))
        if(length(Alabels) != M) stop(paste("Alabels must be of length", M))
        if(length(Apch)) {
            Apch = rep(Apch, len=length(index.nosz))
            for(i in index.nosz)
                points(Amat[i,1],Amat[i,2],pch=Apch[i],cex=Acex[i],col=Acol[i])
        } else {
            for(i in index.nosz)
                text(Amat[i,1], Amat[i,2], Alabels[i], cex=Acex[i],
                     col=Acol[i], adj=Aadj[i])
        }
    }

    if(C) {
        p2 = nrow(Cmat)
        gapC = rep(gapC, len=p2)
        Cadj = rep(Cadj, len=p2)
        Ccex = rep(Ccex, len=p2)
        Ccol = rep(Ccol, len=p2)
        Clwd = rep(Clwd, len=p2)
        Clty = rep(Clty, len=p2)
        if(length(Clabels) != p2)
            stop(paste("length(Clabels) must be equal to", p2))
        for(i in 1:p2) {
            if(is.R()) arrows(0, 0, Cmat[i,1], Cmat[i,2],
                              lwd=Clwd[i], lty=Clty[i], col=Ccol[i]) else
                       arrows(0,0,Cmat[i,1],Cmat[i,2],open=TRUE,
                              lwd=Clwd[i], lty=Clty[i], col=Ccol[i])
            const = 1 + gapC[i] / sqrt(Cmat[i,1]^2 + Cmat[i,2]^2)
            text(const*Cmat[i,1], const*Cmat[i,2], Clabels[i], cex=Ccex[i],
                 adj=Cadj[i], col=Ccol[i])
        }
    }

    if(scores) {
        ugrp = unique(groups)
        nlev = length(ugrp)  # number of groups
        clty = rep(clty, len=nlev)
        clwd = rep(clwd, len=nlev)
        ccol = rep(ccol, len=nlev)
        if(length(spch))
            spch = rep(spch, len=n)
        scol = rep(scol, len=n)
        scex = rep(scex, len=n)
        for(i in ugrp) {
            gp = groups==i
            if(nlev > 1 && (length(unique(spch[gp])) != 1 ||
               length(unique(scol[gp])) != 1 ||
               length(unique(scex[gp])) != 1))
               warning(paste("spch/scol/scex is different for individuals",
                             "from the same group"))

            temp = nuhat[gp,,drop=FALSE]
            if(length(spch)) {
                points(temp[,1], temp[,2], cex=scex[gp], pch=spch[gp],
                       col=scol[gp])
            } else {
                text(temp[,1], temp[,2], label=slabels, cex=scex[gp],
                     col=scol[gp])
            }
            if(chull.arg) {
                hull = chull(temp[,1],temp[,2])
                hull = c(hull, hull[1])
                lines(temp[hull,1], temp[hull,2], type="b", lty=clty[i],
                      col=ccol[i], lwd=clwd[i], pch="  ")
            }
        }
    }

    invisible(nuhat)
}






Coef.rrvglm <- function(object, ...) {
    M <- object@misc$M
    n <- object@misc$n
    colx1.index = object@control$colx1.index
    colx2.index = object@control$colx2.index
    p1 = length(colx1.index)
    Amat <- object@constraints[[colx2.index[1]]]
    B1mat <- if(p1) coef(object, mat= TRUE)[colx1.index,,drop=FALSE] else NULL
    C.try <- coef(object, mat= TRUE)[colx2.index,,drop=FALSE]
    Cmat <- C.try %*% Amat %*% solve(t(Amat) %*% Amat)


    Rank = object@control$Rank
    lv.names = if(Rank>1) paste("lv", 1:Rank, sep="") else "lv"
    dimnames(Amat) = list(object@misc$predictors.names, lv.names)
    dimnames(Cmat) = list(dimnames(Cmat)[[1]], lv.names)

    ans = new(Class="Coef.rrvglm",
      A            = Amat,
      B1           = B1mat,
      C            = Cmat,
      Rank         = Rank,
      colx1.index  = colx1.index,
      colx2.index  = colx2.index)
    if(object@control$Corner)
        ans@Atilde = Amat[-c(object@control$Index.corner,
                         object@control$Structural.zero),,drop=FALSE]
    ans
}

setMethod("Coef", "rrvglm", function(object, ...) Coef.rrvglm(object, ...))

printCoef.rrvglm = function(x, ...) {

    object = x

    cat("\nA matrix:\n")
    print(object@A, ...)
    cat("\n")

    cat("\nC matrix:\n")
    print(object@C, ...)
    cat("\n")

    cat("\nB1 matrix:\n")
    print(object@B1, ...)
    cat("\n")

    invisible(object)
} 


if(is.R()) {
    if(!isGeneric("biplot"))
    setGeneric("biplot", function(x, ...) standardGeneric("biplot")) 
}


setMethod("Coef", "qrrvglm", function(object, ...) Coef.qrrvglm(object, ...))



setMethod("biplot", "qrrvglm",
           function(x, ...) {
           biplot.qrrvglm(x, ...)})

setMethod("lvplot", "qrrvglm",
           function(object, ...) {
           invisible(lvplot.qrrvglm(object, ...))})

setMethod("lvplot", "rrvglm",
           function(object, ...) {
           invisible(lvplot.rrvglm(object, ...))})


biplot.rrvglm = function(x, ...)
    lvplot(object=x, ...)

setMethod("biplot",  "rrvglm", function(x, ...)
           invisible(biplot.rrvglm(x, ...)))




summary.qrrvglm = function(object,
                           varlvI = FALSE, reference = NULL, ...) {
    answer = object
    answer@post$Coef = Coef(object, varlvI = varlvI, reference = reference, 
                            ...) # Store it here; non-elegant

    if(length((answer@post$Coef)@dispersion) &&
       length(object@misc$estimated.dispersion) &&
       object@misc$estimated.dispersion)
        answer@dispersion = 
        answer@misc$dispersion = (answer@post$Coef)@dispersion

    class(answer) = "summary.qrrvglm"
    answer
}

printsummary.qrrvglm = function(x, ...) {



    cat("\nCall:\n")
    dput(x@call)

    print(x@post$Coef, ...) # non-elegant programming

    if(length(x@dispersion) > 1) {
        cat("\nDispersion parameters:\n")
        if(length(x@misc$ynames)) {
            names(x@dispersion) = x@misc$ynames 
            print(x@dispersion, ...)
        } else
            cat(x@dispersion, fill=TRUE)
        cat("\n")
    } else if(length(x@dispersion) == 1) {
        cat("\nDispersion parameter:  ", x@dispersion, "\n")
    }

}

setClass(Class="summary.qrrvglm", representation("qrrvglm"))

setMethod("summary", "qrrvglm",
          function(object, ...)
          summary.qrrvglm(object, ...))

setMethod("print", "summary.qrrvglm",
          function(x, ...)
          invisible(printsummary.qrrvglm(x, ...)))

setMethod("show", "summary.qrrvglm",
          function(object)
          invisible(printsummary.qrrvglm(object)))

setMethod("print", "Coef.rrvglm", function(x, ...)
          invisible(printCoef.rrvglm(x, ...)))

setMethod("show", "Coef.rrvglm", function(object)
          invisible(printCoef.rrvglm(object)))





grc = function(y, Rank=1, Index.corner=2:(1+Rank), Structural.zero=1,
               summary.arg= FALSE, h.step=0.0001, ...) {
                           


    myrrcontrol = rrvglm.control(Rank=Rank, Index.corner=Index.corner,
                                 Structural.zero = Structural.zero, ...)
    object.save = y
    if(is(y, "rrvglm")) {
        y = object.save@y
    } else {
        y = as.matrix(y)
        class(y) = "matrix"   # Needed in R 
    }
    if(length(dim(y)) != 2 || nrow(y) < 3 || ncol(y) < 3)
     stop("y must be a matrix with >= 3 rows & columns, or a rrvglm() object")

    ei = function(i, n) diag(n)[,i,drop=FALSE]
    .grc.df = data.frame(Row2 = ei(2, nrow(y)))

    yn1 = if(length(dimnames(y)[[1]])) dimnames(y)[[1]] else
              paste("x2", 1:nrow(y), sep="")
    warn.save = options()$warn
    options(warn=-3)    # Suppress the warnings (hopefully, temporarily)
    if(any(!is.na(as.numeric(substring(yn1, 1, 1)))))
        yn1 = paste("x2", 1:nrow(y), sep="")
    options(warn=warn.save)

    Row = factor(1:nrow(y))
    modmat.row = if(is.R()) model.matrix(~ Row) else {
        tmp3 = contrasts(Row) 
        dimnames(tmp3) = list(dimnames(tmp3)[[1]], paste("Row", 2:nrow(y), sep=""))
        cbind("(Intercept)"=1, tmp3) 
    }
    Col = factor(1:ncol(y))
    modmat.col = if(is.R()) model.matrix(~ Col) else {
        tmp3 = contrasts(Col)
        dimnames(tmp3) = list(dimnames(tmp3)[[1]], paste("Col", 2:ncol(y), sep=""))
        cbind("(Intercept)"=1, tmp3) 
    }

    cms = list("(Intercept)" = matrix(1, ncol(y), 1))
    for(i in 2:nrow(y)) {
        cms[[paste("Row", i, sep="")]] = matrix(1, ncol(y), 1)
        .grc.df[[paste("Row", i, sep="")]] = modmat.row[,i]
    }
    for(i in 2:ncol(y)) {
        cms[[paste("Col", i, sep="")]] = modmat.col[,i,drop=FALSE]
        .grc.df[[paste("Col", i, sep="")]] = rep(1, nrow(y))
    }
    for(i in 2:nrow(y)) {
        cms[[yn1[i]]] = diag(ncol(y))
        .grc.df[[yn1[i]]] = ei(i, nrow(y))
    }

    dimnames(.grc.df) = list(if(length(dimnames(y)[[1]])) dimnames(y)[[1]] else 
                             as.character(1:nrow(y)),
                             dimnames(.grc.df)[[2]])

    str1 = "~ Row2"
    if(nrow(y)>2) 
    for(i in 3:nrow(y))
        str1 = paste(str1, paste("Row", i, sep=""), sep=" + ")
    for(i in 2:ncol(y))
        str1 = paste(str1, paste("Col", i, sep=""), sep=" + ")
    str2 = paste("y ", str1)
    for(i in 2:nrow(y))
        str2 = paste(str2, yn1[i], sep=" + ")
    myrrcontrol$Norrr = as.formula(str1)  # Overwrite this

    if(is.R()) assign(".grc.df", .grc.df, envir = VGAMenv) else
        .grc.df <<- .grc.df

    warn.save = options()$warn
    options(warn=-3)    # Suppress the warnings (hopefully, temporarily)
    answer = if(is(object.save, "rrvglm")) object.save else 
             rrvglm(as.formula(str2), fam=poissonff,
                    constraints=cms, control=myrrcontrol, data=.grc.df)
    options(warn=warn.save)

    if(summary.arg) {
        class(answer) = "rrvglm"
        answer = summary.rrvglm(answer, h.step=h.step)
    } else { 
        class(answer) = "grc"
    }

    if(is.R()) {
        if(exists(".grc.df", envir = VGAMenv))
            rm(".grc.df", envir = VGAMenv)
    } else {
        remove(".grc.df")
    }

    answer
}

summary.grc = function(object, ...) {
    grc(object, summary.arg= TRUE, ...)
}





trplot.qrrvglm = function(object,
                          whichSpecies=NULL,
                          add=FALSE, plot.it=TRUE,
                          label.sites=FALSE, 
                          sitenames = dimnames(object@y)[[1]],
                          axes.equal = TRUE,
                          cex=par()$cex,
                          col=1:(nos*(nos-1)/2),
                          log="", 
                          lty = rep(par()$lty, len=nos*(nos-1)/2),
                          lwd = rep(par()$lwd, len=nos*(nos-1)/2),
                          tcol= rep(par()$col, len=nos*(nos-1)/2),
                          xlab = NULL, ylab = NULL, 
                          main="",   # "Trajectory plot",
                          type="b", check.ok=TRUE, ...) {
    coef.obj = Coef(object)  # use defaults for those two arguments
    if(coef.obj@Rank != 1) stop("object must be a rank-1 model")
    fv = fitted(object)
    modelno = object@control$modelno  # 1,2,3, or 0
    NOS = ncol(fv)   # Number of species
    M = object@misc$M # 
    nn = nrow(fv)  # Number of sites 
    if(length(sitenames))
        sitenames = rep(sitenames, len=nn)
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
        nos = length(whichSpecies) # nos = number of species to be plotted

    if(length(whichSpecies.numer) <= 1)
        stop("must have at least 2 species to be plotted")
    cx1i = object@control$colx1.index
    if(check.ok)
    if(!(length(cx1i)==1 && names(cx1i)=="(Intercept)"))
        stop("trajectory plots allowable only for Norrr = ~ 1 models")

    first.spp  = iam(1,1,M=M,both=TRUE,diag=FALSE)$row.index
    second.spp = iam(1,1,M=M,both=TRUE,diag=FALSE)$col.index
    myxlab = if(length(whichSpecies.numer)==2) {
                paste("Fitted value for",
                if(is.character(whichSpecies.numer)) whichSpecies.numer[1] else
                    sppNames[whichSpecies.numer[1]])
                 } else "Fitted value for 'first' species"
    myxlab = if(length(xlab)) xlab else myxlab
    myylab = if(length(whichSpecies.numer)==2) {
                paste("Fitted value for",
                if(is.character(whichSpecies.numer)) whichSpecies.numer[2] else
                    sppNames[whichSpecies.numer[2]])
                 } else "Fitted value for 'second' species"
    myylab = if(length(ylab)) ylab else myylab
    if(!add) {
        xxx = if(axes.equal) fv[,whichSpecies.numer] else
              fv[,whichSpecies.numer[first.spp]]
        yyy = if(axes.equal) fv[,whichSpecies.numer] else
              fv[,whichSpecies.numer[second.spp]]
        matplot(xxx, yyy, type="n", log=log, xlab=myxlab,
                ylab=myylab, main=main, ...)
    }

    lwd = rep(lwd, len=nos*(nos-1)/2)
    col = rep(col, len=nos*(nos-1)/2)
    lty = rep(lty, len=nos*(nos-1)/2)
    tcol = rep(tcol, len=nos*(nos-1)/2)

    oo = order(coef.obj@lv)   # Sort by the latent variable
    ii = 0
    col = rep(col, length=nos*(nos-1)/2)
    species.names = NULL
    if(plot.it)
    for(i1 in seq(whichSpecies.numer)) {
        for(i2 in seq(whichSpecies.numer))
            if(i1 < i2) {
                ii = ii + 1
                species.names = rbind(species.names,
                                      cbind(sppNames[i1], sppNames[i2]))
                matplot(fv[oo,whichSpecies.numer[i1]],
                        fv[oo,whichSpecies.numer[i2]],
                        type=type, add=TRUE,
                        lty=lty[ii], lwd=lwd[ii], col=col[ii],
                        pch = if(label.sites) "   " else "*" )
                if(label.sites && length(sitenames))
                    text(fv[oo,whichSpecies.numer[i1]],
                         fv[oo,whichSpecies.numer[i2]],
                         labels=sitenames[oo], cex=cex, col=tcol[ii])
            }
    }
    invisible(list(species.names=species.names, 
                   sitenames=sitenames[oo]))
}

if(!isGeneric("trplot"))
    setGeneric("trplot", function(object, ...) standardGeneric("trplot")) 
setMethod("trplot", "qrrvglm", function(object, ...) trplot.qrrvglm(object, ...))




vcovrrvglm = function(object, ...) {
    summary.rrvglm(object, ...)@cov.unscaled 
}



vcovqrrvglm = function(object,
                       ITolerances = object@control$EqualTolerances,
                       MaxScale = c("predictors", "response"),
           dispersion = rep(if(length(sobj@dispersion)) sobj@dispersion else 1,
                            len=M), ...) {
    stop("this function is not yet completed")

    if(mode(MaxScale) != "character" && mode(MaxScale) != "name")
        MaxScale <- as.character(substitute(MaxScale))
    MaxScale <- match.arg(MaxScale, c("predictors", "response"))[1]
    if(MaxScale != "predictors")
        stop("can currently only handle MaxScale=\"predictors\"")

    sobj = summary(object)
    cobj = Coef(object, ITolerances = ITolerances, ...)
    M = nrow(cobj@A)
    dispersion = rep(dispersion, len=M)
    if(cobj@Rank != 1)
        stop("object must be a rank 1 model")

    dvecMax = cbind(1, -0.5 * cobj@A / c(cobj@D), (cobj@A / c(2*cobj@D))^2)
    dvecTol = cbind(0, 0, 1 / c(-2 * cobj@D)^1.5)
    dvecOpt = cbind(0, -0.5 / c(cobj@D), 0.5 * cobj@A / c(cobj@D^2))

    if((length(object@control$colx1.index) != 1) ||
       (names(object@control$colx1.index) != "(Intercept)"))
        stop("Can only handle Norrr=~1 models")
    okvals=c(3*M,2*M+1) # Tries to correspond to EqualTol==c(FALSE,TRUE) resp.
    if(all(length(coef(object)) != okvals))
        stop("Can only handle intercepts-only model with EqualTolerances=FALSE")

    answer = NULL
    Cov.unscaled = array(NA, c(3,3,M), dimnames=list(
        c("(Intercept)", "lv", "lv^2"),
        c("(Intercept)", "lv", "lv^2"), dimnames(cobj@D)[[3]]))
    for(spp in 1:M) {
        index = c(M+ifelse(object@control$EqualTolerances, 1, M) + spp,
                  spp,
                  M+ifelse(object@control$EqualTolerances, 1, spp))
        vcov = Cov.unscaled[,,spp] =
            sobj@cov.unscaled[index,index] # Order is A, D, B1
        se2Max = dvecMax[spp,,drop=FALSE] %*% vcov %*% cbind(dvecMax[spp,])
        se2Tol = dvecTol[spp,,drop=FALSE] %*% vcov %*% cbind(dvecTol[spp,])
        se2Opt = dvecOpt[spp,,drop=FALSE] %*% vcov %*% cbind(dvecOpt[spp,])
        answer = rbind(answer, dispersion[spp]^0.5 *
                       c(se2Opt=se2Opt, se2Tol=se2Tol, se2Max=se2Max))
    }

    link.function = if(MaxScale=="predictors")
        remove.arg(object@misc$predictors.names[1]) else ""
    dimnames(answer) = list(dimnames(cobj@D)[[3]], c("Optimum", "Tolerance",
        if(nchar(link.function)) paste(link.function,"(Maximum)",sep="") else
        "Maximum"))
    NAthere = is.na(answer %*% rep(1, len=3))
    answer[NAthere,] = NA # NA in tolerance means NA everywhere else
    new(Class="vcov.qrrvglm",
        Cov.unscaled=Cov.unscaled,
        dispersion=dispersion,
        se=sqrt(answer))
}


setMethod("vcov", "rrvglm", function(object, ...)
    vcovrrvglm(object, ...))

setMethod("vcov", "qrrvglm", function(object, ...)
    vcovqrrvglm(object, ...))

setClass(Class="vcov.qrrvglm", representation(
         Cov.unscaled="array",  # permuted cov.unscaled
         dispersion="numeric",
         se="matrix"))



model.matrix.qrrvglm <- function(object, type=c("lv", "vlm"), ...) {

    if(mode(type) != "character" && mode(type) != "name")
    type = as.character(substitute(type))
    type = match.arg(type, c("lv","vlm"))[1]

    switch(type, lv=Coef(object, ...)@lv, vlm=object@x) 
}

setMethod("model.matrix",  "qrrvglm", function(object, ...)
           model.matrix.qrrvglm(object, ...))







persp.qrrvglm = function(x, varlvI = FALSE, reference = NULL,
                  plot.it=TRUE,
                  xlim=NULL, ylim=NULL, zlim=NULL, # zlim ignored if Rank==1
                  gridlength=if(Rank==1) 301 else c(51,51),
                  whichSpecies = NULL,
                  xlab = if(Rank==1) "Latent Variable" else "Latent Variable 1",
                  ylab = if(Rank==1) "Expected Value" else "Latent Variable 2",
                  zlab="Expected value",
                  labelSpecies = FALSE,   # For Rank==1 only
                  stretch = 1.05,  # quick and dirty, Rank==1 only
                  main="",
                  ticktype = "detailed", 
                  col = if(Rank==1) par()$col else "white",
                  add1 = FALSE,
                  ...) {
    oylim = ylim
    object = x  # don't like x as the primary argument 
    coef.obj = Coef(object, varlvI = varlvI, reference = reference)
    if((Rank <- coef.obj@Rank) > 2)
        stop("object must be a rank-1 or rank-2 model")
    fv = fitted(object)
    NOS = ncol(fv)    # Number of species
    M = object@misc$M # 

    xlim = rep(if(length(xlim)) xlim else range(coef.obj@lv[,1]), length=2)
    if(!length(oylim)) {
        ylim = if(Rank==1) c(0, max(fv)*stretch) else
            rep(range(coef.obj@lv[,2]), length=2)
    }
    gridlength = rep(gridlength, length=Rank)
    lv1 = seq(xlim[1], xlim[2], length=gridlength[1])
    if(Rank==1) {
        m = cbind(lv1)
    } else {
        lv2 = seq(ylim[1], ylim[2], length=gridlength[2])
        m = expand.grid(lv1,lv2)
    }

    if(dim(coef.obj@B1)[1] != 1 || dimnames(coef.obj@B1)[[1]] != "(Intercept)")
        stop("Norrr = ~ 1 is needed")
    LP = coef.obj@A %*% t(cbind(m))   # M by n
    LP = LP + c(coef.obj@B1) # Assumes \bix_1 = 1 (intercept only)

    mm = as.matrix(m)
    N = ncol(LP)
    for(j in 1:M) {
        for(i in 1:N) {
            LP[j,i] = LP[j,i] + mm[i,,drop=FALSE] %*% coef.obj@D[,,j] %*%
                                t(mm[i,,drop=FALSE])
        }
    }
    LP = t(LP)   # n by M


    fitvals = object@family@inverse(LP)   # n by NOS
    dimnames(fitvals) = list(NULL, dimnames(fv)[[2]])
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
    if(Rank==1) {
        if(plot.it) {
            if(!length(oylim))
            ylim = c(0, max(fitvals[,whichSpecies.numer])*stretch) # A revision
            col = rep(col, len=length(whichSpecies.numer))
            if(!add1)
            matplot(lv1, fitvals, xlab=xlab, ylab=ylab, type="n", 
                    main=main, xlim=xlim, ylim=ylim, ...) 
            for(j in 1:length(whichSpecies.numer)) {
                ptr2 = whichSpecies.numer[j]  # points to species column
                lines(lv1, fitvals[,ptr2], col=col[j], ...)
                if(labelSpecies) {
                    ptr1=(1:nrow(fitvals))[max(fitvals[,ptr2])==fitvals[,ptr2]]
                    ptr1 = ptr1[1]
                    text(lv1[ptr1], fitvals[ptr1,ptr2]+
                         (stretch-1)*diff(range(ylim)),
                         label=sppNames[j], col=col[j], ...)
                }
            }
        }
    } else {
        maxfitted = matrix(fitvals[,whichSpecies[1]], length(lv1), length(lv2))
        if(length(whichSpecies) > 1)
        for(j in whichSpecies[-1]) {
            maxfitted = pmax(maxfitted, matrix(fitvals[,j], 
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
setMethod("persp", "qrrvglm", function(x, ...) persp.qrrvglm(x=x, ...))




ccoef.qrrvglm = function(object, varlvI = FALSE, reference = NULL, ...) {
    Coef(object, varlvI = varlvI, reference = reference, ...)@C
}

ccoef.Coef.qrrvglm = function(object, ...) {
    if(length(list(...))) warning("Too late! Ignoring the extra arguments")
    object@C
}

lv.qrrvglm <- function(object, varlvI = FALSE, reference = NULL, ...) {
    Coef(object, varlvI = varlvI, reference = reference, ...)@lv
}

lv.rrvglm = function(object, ...) {
    ans = lvplot(object, plot.it=FALSE)
    if(ncol(ans) == 1) dimnames(ans) = list(dimnames(ans)[[1]], "lv")
    ans
}

lv.Coef.qrrvglm = function(object, ...) {
    if(length(list(...))) warning("Too late! Ignoring the extra arguments")
    object@lv
}

Max.qrrvglm = function(object, varlvI = FALSE, reference = NULL, ...) {
    Coef(object, varlvI = varlvI, reference = reference, ...)@Maximum
}

Max.Coef.qrrvglm = function(object, ...) {
    if(length(list(...))) warning("Too late! Ignoring the extra arguments")
    if(any(slotNames(object) == "Maximum")) object@Maximum else
    Max(object, ...)
}

Opt.qrrvglm = function(object, varlvI = FALSE, reference = NULL, ...) {
    Coef(object, varlvI = varlvI, reference = reference, ...)@Optimum
}

Opt.Coef.qrrvglm = function(object, ...) {
    if(length(list(...))) warning("Too late! Ignoring the extra arguments")
    Coef(object, ...)@Optimum
}

Tol.qrrvglm = function(object, varlvI = FALSE, reference = NULL, ...) {
    Coef(object, varlvI = varlvI, reference = reference, ...)@Tolerance
}

Tol.Coef.qrrvglm = function(object, ...) {
    if(length(list(...))) warning("Too late! Ignoring the extra arguments")
    if(any(slotNames(object) == "Tolerance")) object@Tolerance else
    Tol(object, ...)
}


if(!isGeneric("ccoef"))
    setGeneric("ccoef", function(object, ...) standardGeneric("ccoef")) 
setMethod("ccoef",  "rrvglm", function(object, ...) ccoef.qrrvglm(object, ...))
setMethod("ccoef", "qrrvglm", function(object, ...) ccoef.qrrvglm(object, ...))
setMethod("ccoef",  "Coef.rrvglm", function(object, ...) ccoef.Coef.qrrvglm(object, ...))
setMethod("ccoef", "Coef.qrrvglm", function(object, ...) ccoef.Coef.qrrvglm(object, ...))

setMethod("coef", "qrrvglm", function(object, ...) Coef.qrrvglm(object, ...))
setMethod("coefficients", "qrrvglm", function(object, ...) Coef.qrrvglm(object, ...))

if(!isGeneric("lv"))
    setGeneric("lv", function(object, ...) standardGeneric("lv")) 
setMethod("lv",  "rrvglm", function(object, ...) lv.rrvglm(object, ...))
setMethod("lv", "qrrvglm", function(object, ...) lv.qrrvglm(object, ...))
setMethod("lv",  "Coef.rrvglm", function(object, ...) lv.Coef.qrrvglm(object, ...))
setMethod("lv", "Coef.qrrvglm", function(object, ...) lv.Coef.qrrvglm(object, ...))

if(!isGeneric("Max"))
    setGeneric("Max", function(object, ...) standardGeneric("Max")) 
setMethod("Max", "qrrvglm", function(object, ...) Max.qrrvglm(object, ...))
setMethod("Max", "Coef.qrrvglm", function(object, ...) Max.Coef.qrrvglm(object, ...))

if(!isGeneric("Opt"))
    setGeneric("Opt", function(object, ...) standardGeneric("Opt"))
setMethod("Opt", "qrrvglm", function(object, ...) Opt.qrrvglm(object, ...))
setMethod("Opt", "Coef.qrrvglm", function(object, ...) Opt.Coef.qrrvglm(object, ...))

if(!isGeneric("Tol"))
    setGeneric("Tol", function(object, ...) standardGeneric("Tol")) 
setMethod("Tol", "qrrvglm", function(object, ...) Tol.qrrvglm(object, ...))
setMethod("Tol", "Coef.qrrvglm", function(object, ...) Tol.Coef.qrrvglm(object, ...))



cgo <- function(...) {
    stop("The function \"cgo\" has been renamed \"cqo\". Ouch! Sorry!")
}

clo <- function(...) {
    stop("Constrained linear ordination is fitted with the function \"rrvglm\"")
}




is.bell.vlm <-
is.bell.rrvglm <- function(object, ...) {
    M = object@misc$M
    ynames = object@misc$ynames
    ans = rep(FALSE, len=M)
    if(length(ynames)) names(ans) = ynames
    ans
}

is.bell.uqo <-
is.bell.qrrvglm <- function(object, ...) {
    is.finite(Max(object, ...))
}

is.bell.cao <- function(object, ...) {
    NA * Max(object, ...)
}

if(!isGeneric("is.bell"))
    setGeneric("is.bell", function(object, ...) standardGeneric("is.bell"))
setMethod("is.bell","uqo", function(object, ...) is.bell.uqo(object, ...))
setMethod("is.bell","qrrvglm", function(object,...) is.bell.qrrvglm(object,...))
setMethod("is.bell","rrvglm", function(object, ...) is.bell.rrvglm(object, ...))
setMethod("is.bell","vlm", function(object, ...) is.bell.vlm(object, ...))
setMethod("is.bell","cao", function(object, ...) is.bell.cao(object, ...))
setMethod("is.bell","Coef.qrrvglm", function(object,...) is.bell.qrrvglm(object,...))






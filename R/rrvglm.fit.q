# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.







rrvglm.fit <- function(x, y, w=rep(1, length(x[, 1])),
    etastart=NULL, mustart=NULL, coefstart=NULL,
    offset=0, family,
    control=rrvglm.control(...),
    criterion="coefficients",
    qr.arg=FALSE,
    constraints=NULL,
    extra=NULL,
    Terms=Terms, function.name="rrvglm", ...)
{
    post = list()
    check.rank = TRUE # !control$Quadratic
    nonparametric <- FALSE
    epsilon <- control$epsilon
    maxit <- control$maxit
    backchat <- control$backchat && !control$Quadratic  # rrr;
    save.weight <- control$save.weight
    trace <- control$trace
    orig.stepsize <- control$stepsize
    minimize.criterion <- control$min.criterion


    n <- dim(x)[1]

    new.s.call <- expression({
        if(c.list$one.more)
        {
            fv <- if(backchat) {
                      if(M>1) matrix(c.list$fit,n,M,byrow=TRUE) else c.list$fit
                  } else c.list$fit
            new.coeffs <- c.list$coeff

            if(length(family@middle))
                eval(family@middle)

            eta <- fv + offset

            mu <- family@inverse(eta, extra)

            if(length(family@middle2))
                eval(family@middle2)

            old.crit <- new.crit
            new.crit <- 
                switch(criterion,
                    coefficients=new.coeffs,
                    tfun(mu=mu, y=y, w=w, res=FALSE, eta=eta, extra))



            if(trace && orig.stepsize==1) {
                cat(if(control$Quadratic) "QRR-VGLM" else "RR-VGLM",
                    "   linear loop ", iter, ": ", criterion, "= ")
                uuuu =
                    switch(criterion,
                    coefficients= if(is.R())
                        format(new.crit, dig=round(2-log10(epsilon))) else
                        format(round(new.crit, round(2-log10(epsilon)))),
                    format(round(new.crit, 4)))
                switch(criterion,
                    coefficients={if(length(new.crit) > 2) cat("\n");
                       cat(uuuu, fill=TRUE, sep=", ")},
                    cat(uuuu, fill=TRUE, sep=", "))
           }

            {
                take.half.step <- (control$half.stepsizing && length(old.coeffs)) && 
                             !control$Quadratic &&
                             ((orig.stepsize!=1) ||
                              (criterion!="coefficients" &&
                             (if(minimize.criterion) new.crit > old.crit else
                             new.crit < old.crit)))
                if(take.half.step) {
                    stepsize <- 2 * min(orig.stepsize, 2*stepsize)
                    new.coeffs.save <- new.coeffs
                    if(trace) 
                        cat("Taking a modified step")
                    repeat {
                        if(trace) {
                            cat(".")
                            if(exists("flush.console"))
                                flush.console()
                        }
                        stepsize <- stepsize / 2
                        if(too.small <- stepsize < 0.001)
                            break
                        new.coeffs <- (1-stepsize)*old.coeffs +
                                       stepsize*new.coeffs.save

                        if(length(family@middle))
                            eval(family@middle)

                        fv <- xbig.save %*% new.coeffs
                        if(M > 1)
                            fv <- matrix(fv, n, M, byrow=TRUE)

                        eta <- fv + offset

                        mu <- family@inverse(eta, extra)

                        if(length(family@middle2))
                            eval(family@middle2)


                        new.crit <- 
                            switch(criterion,
                                coefficients=new.coeffs,
                                tfun(mu=mu,y=y,w=w,res=FALSE,eta=eta,extra))

                        if((criterion=="coefficients") || 
                           ( minimize.criterion && new.crit < old.crit) ||
                           (!minimize.criterion && new.crit > old.crit))
                            break
                    }

                    if(trace) 
                        cat("\n")
                    if(too.small) {
                        warning(paste("iterations terminated because",
                              "half-step sizes are very small"))
                        one.more <- FALSE
                    } else {
                        if(trace) {
                       cat(if(control$Quadratic) "QRR-VGLM" else "RR-VGLM",
                    "   linear loop ", iter, ": ", criterion, "= ")
                            uuuu =
                            switch(criterion,
                            coefficients=if(is.R())
                            format(new.crit, dig=round(2-log10(epsilon))) else
                            format(round(new.crit, round(2-log10(epsilon)))),
                            format(round(new.crit, 4)))

                            switch(criterion,
                            coefficients={if(length(new.crit) > 2) cat("\n");
                               cat(uuuu, fill=TRUE, sep=", ")},
                            cat(uuuu, fill=TRUE, sep=", "))
                        }

                        one.more <- eval(control$convergence)
                    }
                } else {
                    one.more <- eval(control$convergence)
                }
            }
            if(exists("flush.console"))
                flush.console()

            if(one.more) {
                iter <- iter + 1
                deriv.mu <- eval(family@deriv)
                wz <- eval(family@weight)
                if(control$checkwz)
                    wz = checkwz(wz, M=M, trace=trace, wzeps=control$wzepsilon)


                wz = matrix(wz, nrow=n)   # zz 3/10/05
                U <- vchol(wz, M=M, n=n, silent=!trace)
                tvfor <- vforsub(U, as.matrix(deriv.mu), M=M, n=n)
                z = eta + vbacksub(U, tvfor, M, n) - offset # Contains \bI \bnu

                rrr.expression = paste("rrr", control$Algorithm,
                                       "expression", sep=".")
                rrr.expression = get(rrr.expression)
                eval(rrr.expression)

                c.list$z <- z  # contains \bI_{Rank} \bnu
                c.list$U <- U
                if(copyxbig) c.list$xbig <- xbig.save
            }

            c.list$one.more <- one.more
            c.list$coeff <- runif(length(new.coeffs)) # 12/3/03; twist needed!
            old.coeffs <- new.coeffs
        }
        c.list
    }) # end of new.s.call





    copyxbig <- FALSE    # May be overwritten in @initialize
    stepsize <- orig.stepsize
    old.coeffs <- coefstart

    intercept.only <- ncol(x) == 1 && dimnames(x)[[2]] == "(Intercept)"
    y.names <- predictors.names <- NULL    # May be overwritten in @initialize

 
    n.save <- n 



    Rank <- control$Rank
    rrcontrol <- control  #

    if(length(slot(family, "initialize")))
        eval(slot(family, "initialize")) # Initialize mu & M (and optionally w)


    eval(rrr.init.expression)  

    if(length(etastart)) {
        eta <- etastart
        mu <- if(length(mustart)) mustart else
              slot(family, "inverse")(eta, extra)
    } else {
        if(length(mustart))
            mu <- mustart
        eta <- slot(family, "link")(mu, extra)
    }

    M <- if(is.matrix(eta)) ncol(eta) else 1

    if(is.character(rrcontrol$Dzero)) {
        index = match(rrcontrol$Dzero, dimnames(as.matrix(y))[[2]]) 
        if(any(is.na(index)))
            stop("Dzero argument didn't fully match y-names")
        if(length(index) == M)
            stop(paste("all linear predictors are linear in the",
                       "latent variable(s); so set Quadratic=FALSE"))
        rrcontrol$Dzero = control$Dzero = index
    }





    if(length(family@constraints))
        eval(family@constraints)


    special.matrix = matrix(-34956.125, M, M)    # An unlikely used matrix 
    just.testing <- cm.vgam(special.matrix, x, rrcontrol$Norrr, constraints)
    findex = trivial.constraints(just.testing, special.matrix)
    tc1 = trivial.constraints(constraints)

    if(!control$Quadratic && sum(!tc1)) {
        for(ii in names(tc1))
            if(!tc1[ii] && !any(ii == names(findex)[findex==1]))
                warning(paste("\"", ii, "\"", " is a non-trivial constraint",
           " that will be overwritten by reduced-rank regression", sep=""))
    }

    if(all(findex))
        stop("use vglm(), not rrvglm()!")
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
    Index.corner = control$Index.corner




    Amat <- if(length(rrcontrol$Ainit)) rrcontrol$Ainit else
            matrix(rnorm(M * Rank, sd=rrcontrol$SD.Cinit), M, Rank)
    Cmat <- if(length(rrcontrol$Cinit)) rrcontrol$Cinit else {
                if(!rrcontrol$Use.Init.Poisson.QO) {
                    matrix(rnorm(p2 * Rank, sd=rrcontrol$SD.Cinit), p2, Rank)
                } else
                .Init.Poisson.QO(ymat=as.matrix(y), 
                                 X1=x[,colx1.index,drop=FALSE],
                                 X2=x[,colx2.index,drop=FALSE],
                                 Rank=rrcontrol$Rank, trace=rrcontrol$trace,
                                 max.ncol.etamat = rrcontrol$Etamat.colmax,
                                 Crow1positive=rrcontrol$Crow1positive,
                                 isdlv=rrcontrol$isdlv)
            }

    if(modelno==3) 
        Amat[c(FALSE,TRUE),] <- 0  # Intercept only for log(k)


    if(control$Corner)
        Amat[control$Index.corner,] = diag(Rank)
    if(length(control$Structural.zero))
        Amat[control$Structural.zero,] = 0

    rrcontrol$Ainit = control$Ainit = Amat   # Good for valt()
    rrcontrol$Cinit = control$Cinit = Cmat   # Good for valt()

    Blist <- process.constraints(constraints, x, M)

    nice31 = control$Quadratic && (!control$EqualTol || control$ITolerances) &&
             all(trivial.constraints(Blist))

    Blist = Blist.save = replace.constraints(Blist, Amat, colx2.index)


    ncolBlist <- unlist(lapply(Blist, ncol))
    dimB <- sum(ncolBlist)


    xbig.save <- if(control$Quadratic) {
        tmp500 = lm2qrrvlm.model.matrix(x=x, Blist=Blist,
                       C=Cmat, control=control)
        xsmall.qrr = tmp500$new.lv.model.matrix 
        B.list = tmp500$constraints # Doesn't change or contain \bI_{Rank} \bnu
        if(modelno==3 && FALSE) {
            B.list[[1]] = (B.list[[1]])[,c(TRUE,FALSE),drop=FALSE] # Amat
            B.list[[2]] = (B.list[[2]])[,c(TRUE,FALSE),drop=FALSE] # D
        }

        lv.mat = tmp500$lv.mat
        if(length(tmp500$offset)) {
            offset = tmp500$offset 
        }
        lm2vlm.model.matrix(xsmall.qrr, B.list, xij=control$xij)
    } else {
        lv.mat = x[,colx2.index,drop=FALSE] %*% Cmat 
        lm2vlm.model.matrix(x, Blist, xij=control$xij)
    }




    if(length(coefstart)) {
        eta <- if(ncol(xbig.save)>1) xbig.save %*% coefstart +
                   offset else xbig.save * coefstart + offset
        eta <- if(M > 1) matrix(eta, ncol=M, byrow=TRUE) else c(eta) 


        mu <- family@inverse(eta, extra)
    }

    if(criterion != "coefficients") {
        tfun <- slot(family, criterion)   # family[[criterion]]
    }

    iter <- 1
    new.crit <- switch(criterion,
                      coefficients=1,
                      tfun(mu=mu, y=y, w=w, res=FALSE, eta=eta, extra))
    old.crit <- if(minimize.criterion) 10*new.crit+10 else -10*new.crit-10


    deriv.mu <- eval(family@deriv)

    wz <- eval(family@weight)
    if(control$checkwz)
        wz = checkwz(wz, M=M, trace=trace, wzeps=control$wzepsilon)

    U <- vchol(wz, M=M, n=n, silent=!trace)
    tvfor <- vforsub(U, as.matrix(deriv.mu), M=M, n=n)
    z <- eta + vbacksub(U, tvfor, M=M, n=n) - offset

    c.list <- list(z=as.double(z), fit=as.double(t(eta)), one.more=TRUE,
                   coeff=as.double(rep(1,ncol(xbig.save))), U=as.double(U),
                   copyxbig=copyxbig,
                   xbig=if(copyxbig) as.double(xbig.save) else double(3))



    dxbig <- as.integer(dim(xbig.save))
    n.big <- dxbig[[1]]
    p.big <- dxbig[[2]]

    if(n.big < p.big)
        stop(paste(p.big, "parameters but only", n.big, "observations"))

    if(backchat) {
        nframe <- sys.nframe()
        dotC(name="v_init_call", as.integer(nframe), new.s.call)
    }


    if(backchat) {
        if(control$Quadratic)
            stop("29/11/02; can't handle zedd") 
        index.vglm <- iam(NA, NA, M, both=TRUE)
        tfit <- dotFortran(name="vglmf", as.double(xbig.save), n.big, p.big,
                    backchat=as.integer(TRUE),
                    as.integer(n), as.double(z),   # not zedd 
                    coefficients=double(p.big),
                    predictors=double(n.big), effects=double(n.big),
                    qr=as.double(xbig.save), qraux=double(p.big),
                    rank=as.integer(0), pivot=as.integer(seq(p.big)),
                    work=double(max(n.big, 2 * p.big)),
                    wkmm=double(M*M*5 + M*p.big),
                    as.double(U), as.integer(M),
                    dimu=as.integer(if(is.matrix(U)) nrow(U) else 1),
                    dimm=as.integer(if(is.matrix(wz)) ncol(wz) else 1),
                    as.integer(index.vglm$row), as.integer(index.vglm$col),
                    copyxbig=as.integer(copyxbig),
                    rss=double(1))
    } else {
        bf.call = expression(vlm.wfit(xbig.save, zedd, 
            Blist=if(control$Quadratic) B.list else Blist,
            ncolx=ncol(x), U=U,
            Eta.range = control$Eta.range,
            matrix.out=if(control$Quadratic) FALSE else TRUE, XBIG=TRUE,
            qr=qr.arg, xij=control$xij))

        while(c.list$one.more) {
            if(control$Quadratic) {
                zedd = as.matrix(z)
                if(control$Corner)
                    zedd[,Index.corner] = zedd[,Index.corner] - lv.mat 
            } else {
                zedd = z 
            }

            if(!nice31)
                tfit <- eval(bf.call)   # tfit$fitted.values is n x M

            if(!control$Quadratic) {
                Cmat = tfit$mat.coef[colx2.index,,drop=FALSE] %*%
                       Amat %*% solve(t(Amat) %*% Amat)
                rrcontrol$Ainit = control$Ainit = Amat  # Good for valt()
                rrcontrol$Cinit = control$Cinit = Cmat  # Good for valt()
            }
    
            if(!nice31) c.list$coeff <- tfit$coefficients 
    
            if(control$Quadratic) {
                if(control$Corner)
                    tfit$fitted.values[,Index.corner] =
                        tfit$fitted.values[,Index.corner] + lv.mat 
            }

            if(!nice31)
                tfit$predictors = tfit$fitted.values # Doesn't contain the offset
            if(!nice31) c.list$fit = tfit$fitted.values
            c.list <- eval(new.s.call)
            NULL
        }
    }

    if(maxit>1 && iter>=maxit)
        warning(paste("convergence not obtained in", maxit, "iterations."))


    dn.big <- labels(xbig.save)
    xn.big <- dn.big[[2]]
    yn.big <- dn.big[[1]]

    if(length(family@fini))
        eval(family@fini)

    if(M>1 && !nice31)
        tfit$predictors <- matrix(tfit$predictors, n, M, byrow=backchat)

    asgn <- attr(xbig.save, "assign")
    if(nice31) {
        coefs <- rep(0, len=length(xn.big)) # zz
        rank <- p.big  # zz 3/10/05
    } else {
        coefs <- tfit$coefficients
        names(coefs) <- xn.big
        rank <- tfit$rank
    }

    cnames <- xn.big

    if(check.rank && rank < p.big)
        stop("rrvglm only handles full-rank models (currently)")

    if(nice31) {
        R <- matrix(as.numeric(NA), 5, 5) # zz 3/10/05 
    } else {
        R <- if(is.R()) tfit$qr$qr[1:p.big, 1:p.big, drop=FALSE] else {
                 if(backchat) tfit$qr[1:p.big, 1:p.big, drop=FALSE] else 
                              tfit$qr$qr[1:p.big, 1:p.big, drop=FALSE]
             }
        R[lower.tri(R)] <- 0
        attributes(R) <- if(is.R()) list(dim=c(p.big, p.big),
                         dimnames=list(cnames, cnames), rank=rank) else 
                      list(dim=c(p.big, p.big),
                         dimnames=list(cnames, cnames), rank=rank, class="upper")
    }

    if(nice31) {
        effects <- rep(0, len=77) # zz 3/10/05 
    } else {
        effects <- tfit$effects
        neff <- rep("", n.big)
        neff[seq(p.big)] <- cnames
        names(effects) <- neff

        dim(tfit$predictors) <- c(n, M)
    }
    dn <- labels(x)
    yn <- dn[[1]]
    xn <- dn[[2]]

    if(nice31) {
        residuals <- z - fv  # zz - offset ??   # not sure 3/10/05
        if(M==1) {
            residuals <- as.vector(residuals)
            names(residuals) <- yn
        } else {
            dimnames(residuals) <- list(yn, predictors.names)
        }
    } else {
        residuals <- z - tfit$predictors   # zz - offset ??
        if(M==1) {
            tfit$predictors <- as.vector(tfit$predictors)
            residuals <- as.vector(residuals)
            names(residuals) <- names(tfit$predictors) <- yn
        } else {
            dimnames(residuals) <- dimnames(tfit$predictors) <- list(yn, predictors.names)
        }
    }

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



    df.residual <- n.big - rank - (if(control$Quadratic) Rank*p2 else 0)
    fit <- list(assign=asgn,
                coefficients=coefs,
                constraints=if(control$Quadratic) B.list else Blist,
                df.residual=df.residual,
                df.total=n*M,
                effects=effects, 
                fitted.values=mu,
                offset=offset, 
                rank=rank,
                residuals=residuals,
                R=R,
                terms=Terms) # terms: This used to be done in vglm() 

    if(qr.arg && !nice31) {
        fit$qr <- if(is.R())
                      fit$qr <- tfit$qr  else {
                      if(backchat) tfit[c("qr", "rank", "pivot", "qraux")] else
                                   tfit$qr
                  }
        dimnames(fit$qr$qr) <- dn.big
    }

    if(M==1) {
        wz <- as.vector(wz)  # Convert wz into a vector
    } # else
    fit$weights <- if(save.weight) wz else NULL


    misc <- list(
        colnames.x = xn,
        colnames.xbig = xn.big,
        criterion = criterion,
        function.name = function.name, 
        intercept.only=intercept.only,
        predictors.names = predictors.names,
        M = M,
        n = n,
        nonparametric = nonparametric,
        n.big = n.big,
        orig.assign = attr(x, "assign"),
        p = ncol(x),
        p.big = p.big,
        ynames = dimnames(y)[[2]])

    if(one.more)
        misc$rrr.expression = rrr.expression # 


    crit.list <- list()
    if(criterion != "coefficients")
        crit.list[[criterion]] <- fit[[criterion]] <- new.crit

    for(i in names(.min.criterion.VGAM)) {
        if(i != criterion &&
            any(slotNames(family) == i) &&
            (( is.R() && length(body(slot(family, i)))) ||
            ((!is.R() && length(slot(family, i)) > 1)))) {
                fit[[i]] <- crit.list[[i]] <-
                (slot(family, i))(mu=mu, y=y, w=w, res=FALSE, eta=eta, extra)
        }
    }

                                    



    if(w[1] != 1 || any(w != w[1]))
        fit$prior.weights <- w

    if(length(family@last))
        eval(family@last)


    structure(c(fit, list(predictors=if(nice31) matrix(eta,n,M) else tfit$predictors,
        contrasts=attr(x, "contrasts"),
        control=control,
        crit.list=crit.list,
        extra=extra,
        family=family,
        iter=iter,
        misc=misc,
        post = post,
        rss=if(nice31) 000 else tfit$rss,
        x=x,
        y=y)),
        vclass=family@vfamily)
}



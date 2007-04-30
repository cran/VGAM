# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.




vgam.fit <- function(x, y, w, mf,
        etastart, mustart, coefstart,
        offset, family, control, criterion="coefficients",
        constraints=NULL, extra, qr.arg,
        Terms,
        nonparametric, smooth.labels,
        function.name="vgam", ...)
{
    post = list()
    check.Rank <- TRUE # Set this to false for family functions vppr() etc.
    epsilon <- control$epsilon
    maxit <- control$maxit
    save.weight <- control$save.weight
    trace <- control$trace
 
    bf.maxit <- control$bf.maxit
    bf.epsilon <- control$bf.epsilon
    se.fit <- control$se.fit
    minimize.criterion <- control$min.criterion

    n <- dim(x)[1]




    # --------------------------------------------------------------
    new.s.call <- expression({
        if(c.list$one.more)
        {
            fv <- if(backchat && M>1 && nonparametric)
                matrix(c.list$fit,n,M,byrow=TRUE) else c.list$fit
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
            if(trace) {
                cat("VGAM ", bf, " loop ", iter, ": ", criterion, "= ")

                uuuu = switch(criterion,
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

            if(exists("flush.console"))
                flush.console()

            if(!is.finite(one.more) || !is.logical(one.more)) one.more = FALSE
            if(one.more)
            {
                iter <- iter + 1
                deriv.mu <- eval(family@deriv)
                wz <- eval(family@weight)
                if(control$checkwz)
                    wz = checkwz(wz, M=M, trace=trace, wzeps=control$wzepsilon)

                U <- vchol(wz, M=M, n=n, silent=!trace)
                tvfor <- vforsub(U, as.matrix(deriv.mu), M=M, n=n)
                z <- eta + vbacksub(U, tvfor, M=M, n=n) - offset

                c.list$z <- z
                c.list$wz <- wz
                c.list$U <- U
            }

            c.list$one.more <- one.more
            c.list$coeff <- runif(length(new.coeffs)) # 12/3/03; twist needed!
            old.coeffs <- new.coeffs
        }
        c.list
    })





    backchat <-  control$backchat     # if(is.R()) FALSE else TRUE 
    old.coeffs <- coefstart

    intercept.only <- ncol(x) == 1 && dimnames(x)[[2]] == "(Intercept)"
    y.names <- predictors.names <- NULL    # May be overwritten in @initialize

    n.save <- n
    if(length(slot(family, "initialize")))
        eval(slot(family, "initialize")) # Initialize mu & M (and optionally w)

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


    if(length(family@constraints))
        eval(family@constraints)
    Blist <- process.constraints(constraints, x, M)

    ncolBlist <- unlist(lapply(Blist, ncol))
    dimB <- sum(ncolBlist)


    if(nonparametric) {



        smooth.frame <- mf   # zz
        assignx <- attr(x, "assign")
        which <- assignx[smooth.labels]

        bf <- "s.vam"
        bf.call <- parse(text=paste(
                "s.vam(x, z, wz, tfit$smooth, which, tfit$smooth.frame,",
                "bf.maxit, bf.epsilon, trace, se=se.fit, xbig.save, ",
                "Blist, ncolBlist, M, qbig, U, backchat, ",
                "all.knots=control$all.knots, nk=control$nk)",
                sep=""))[[1]]

        qbig <- sum(ncolBlist[smooth.labels])  # Number of component funs
        s <- matrix(0, n, qbig)
        dy <- if(is.matrix(y)) dimnames(y)[[1]] else names(y)
        d2 <- if(is.null(predictors.names))
            paste("(Additive predictor ",1:M,")", sep="") else
            predictors.names
        dimnames(s) <- list(dy, vlabel(smooth.labels,
              ncolBlist[smooth.labels], M))

        tfit <- list(smooth=s, smooth.frame=smooth.frame)
    } else {
        bf.call <- if(is.R()) expression(vlm.wfit(xbig.save, z, Blist=NULL,
            U=U, matrix.out=FALSE, XBIG=TRUE, qr=qr.arg, xij=NULL)) else 
                              expression(vlm.wfit(xbig.save, z, Blist=NULL,
            U=U, matrix.out=FALSE, XBIG=TRUE, singular.ok=TRUE, qr=qr.arg,
            xij=NULL))
        bf <- "vlm.wfit"
    }

    xbig.save <- lm2vlm.model.matrix(x, Blist, xij=control$xij)


    if(length(coefstart)) {
        eta <- if(ncol(xbig.save)>1) xbig.save %*% coefstart +
                   offset else xbig.save * coefstart + offset
        eta <- if(M > 1) matrix(eta, ncol=M, byrow=TRUE) else c(eta)
        mu <- family@inverse(eta, extra)
    }


    if(criterion != "coefficients") {
        tfun <- slot(family, criterion) # Needed for R, so have to follow suit
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

    c.list <- list(wz=as.double(wz), z=as.double(z),
                   fit=as.double(t(eta)),
                   one.more=TRUE, U=as.double(U),
                   coeff=as.double(rep(1,ncol(xbig.save))))


    dxbig <- as.integer(dim(xbig.save))
    n.big <- dxbig[[1]]
    p.big <- dxbig[[2]]
    if(n.big < p.big)
        stop(paste(p.big, "parameters but only", n.big,
             "observations"))

    if(backchat) {
        nframe <- sys.nframe()
        dotC(name="v_init_call", as.integer(nframe), new.s.call)
    }

    while(c.list$one.more) {
        tfit <- eval(bf.call)   # fit$smooth.frame is new

            c.list$coeff <- tfit$coefficients

        tfit$predictors <- tfit$fitted.values + offset

        c.list$fit <- tfit$fitted.values
        c.list <- eval(new.s.call)
        NULL
    }

    if(maxit>1 && iter>=maxit)
        warning(paste("convergence not obtained in ", maxit, " iterations."))


    dn.big <- labels(xbig.save)
    xn.big <- dn.big[[2]]
    yn.big <- dn.big[[1]]

    if(length(family@fini))
        eval(family@fini)

    coefs <- tfit$coefficients
    asgn <- attr(xbig.save, "assign")    # 29/11/01 was x 

    names(coefs) <- xn.big
    cnames <- xn.big

    if(!is.null(tfit$rank)) {
        rank <- tfit$rank
        if(rank < ncol(x)) 
            stop("rank < ncol(x) is bad")
    } else rank <- ncol(x)    # zz 8/12/01 I think rank is all wrong

    R <- if(is.R()) tfit$qr$qr[1:p.big, 1:p.big, drop=FALSE] else {
             if(backchat) tfit$qr[1:p.big, 1:p.big, drop=FALSE] else
                          tfit$qr$qr[1:p.big, 1:p.big, drop=FALSE]
         }
    R[lower.tri(R)] <- 0
    attributes(R) <- if(is.R()) list(dim=c(p.big, p.big),
                     dimnames=list(cnames, cnames), rank=rank) else
                  list(dim=c(p.big, p.big),
                     dimnames=list(cnames, cnames), rank=rank, class="upper")



    dn <- labels(x)
    yn <- dn[[1]]
    xn <- dn[[2]]



    if(is.matrix(mu)) {
          if(length(dimnames(mu)[[2]])) {
              y.names <- dimnames(mu)[[2]]
          } else
          if(length(dimnames(y)[[2]])) {
              y.names <- dimnames(y)[[2]]
          }
          dimnames(mu) <- list(yn, y.names)
    } else {
        names(mu) <- names(fv)
    }

    tfit$fitted.values <- NULL      # Have to kill it off  3/12/01
    fit <- structure(c(tfit, list(
                assign=asgn,
                constraints=Blist,
                control=control,
                fitted.values=mu,
                formula=as.vector(attr(Terms, "formula")),
                iter=iter,
                offset=offset,
                rank=rank,
                R=R,
                terms=Terms)))

    df.residual <- n.big - rank 
    if(backchat) {
        fit$coefficients <- coefs
        fit$df.residual <- df.residual
    }

    if(!se.fit) {
        fit$var <- NULL
    }

    if(M==1) {
        wz <- as.vector(wz)  # Convert wz into a vector
    } # else
    fit$weights <- if(save.weight) wz else NULL



    if(M==1) {
        fit$predictors <- as.vector(fit$predictors)
        fit$residuals <- as.vector(fit$residuals)
        names(fit$residuals) <- names(fit$predictors) <- yn
    } else
        dimnames(fit$residuals) <- dimnames(fit$predictors) <-
            list(yn, predictors.names)

    NewBlist <- process.constraints(constraints, x, M, by.col=FALSE)

    misc <- list(
        colnames.x = xn,
        colnames.xbig = xn.big,
        criterion = criterion,
        function.name = function.name,
        intercept.only=intercept.only,
        predictors.names = predictors.names,
        M = M,
        n = n,
        new.assign = new.assign(x, NewBlist),
        nonparametric = nonparametric,
        n.big = n.big,
        orig.assign = attr(x, "assign"),
        p = ncol(x),
        p.big = p.big,
        ynames = dimnames(y)[[2]])


    if(criterion != "coefficients")
        fit[[criterion]] <- new.crit



    if(se.fit && length(fit$s.xargument)) {
        misc$varassign <- 
            varassign(Blist, names(fit$s.xargument)) # zz or constraints?
    }



    if(nonparametric) {
        misc$smooth.labels <- smooth.labels
    }


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





    if(M==1) {
        fit$predictors <- as.vector(fit$predictors)
        fit$residuals <- as.vector(fit$residuals)
        names(fit$residuals) <- names(fit$predictors) <- yn
    } else
        dimnames(fit$residuals) <- dimnames(fit$predictors) <-
            list(yn, predictors.names)



 
    if(w[1] != 1 || any(w != w[1]))
        fit$prior.weights <- w

    if(length(family@last))
        eval(family@last)


    if(!is.null(fit$smooth)) {
        fit$nl.chisq <- vgam.nlchisq(fit$qr, fit$resid, wz=wz,
                                     s=fit$smooth, deriv=deriv.mu, U=U,
                                     smooth.labels, attr(x, "assign"),
                                     M=M, n=n, constraints=Blist)
    }


    if(!qr.arg) { 
        fit$qr <- NULL
    }




    fit$misc = NULL # 8/6/02; It's necessary to kill it as it exists in vgam
    structure(c(fit, list(
        contrasts=attr(x, "contrasts"),
        control=control,
        crit.list=crit.list,
        extra=extra,
        family=family,
        iter=iter,
        misc=misc,
        post=post,
        x=x,
        y=y)),
        vclass=family@vfamily)
}





new.assign <- function(X, Blist)
{

    M <- nrow(Blist[[1]])
    dn <- labels(X)
    xn <- dn[[2]]

    asgn <- attr(X, "assign")
    nasgn <- names(asgn)
    lasgn <- unlist(lapply(asgn, length))

    ncolBlist <- unlist(lapply(Blist, ncol))
    names(ncolBlist) <- NULL    # This is necessary for below to work 

    temp2 <- vlabel(nasgn, ncolBlist, M)
    L <- length(temp2)
    newasgn <- vector("list", L)

    k <- 0
    low <- 1
    for(i in 1:length(asgn)) {
        len <- low:(low+ncolBlist[i]*lasgn[i]-1)
        temp <- matrix(len, ncolBlist[i], lasgn[i])
        for(m in 1:ncolBlist[i])
            newasgn[[k+m]] <- temp[m,]
        low <- low + ncolBlist[i]*lasgn[i]
        k <- k + ncolBlist[i]
    }

    names(newasgn) <- temp2
    newasgn
}

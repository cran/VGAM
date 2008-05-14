# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.




vglm.fit <- function(x, y, w=rep(1, length(x[, 1])),
    etastart=NULL, mustart=NULL, coefstart=NULL,
    offset=0, family,
    control=vglm.control(),
    criterion="coefficients",
    qr.arg=FALSE,
    constraints=NULL,
    extra=NULL,
    Terms=Terms, function.name="vglm", ...)
{
    specialCM = NULL
    post = list()
    check.rank <- TRUE # Set this to false for family functions vppr() etc.
    nonparametric <- FALSE
    epsilon <- control$epsilon
    maxit <- control$maxit
    backchat <- control$backchat # FALSE 
    save.weight <- control$save.weight
    trace <- control$trace
    orig.stepsize <- control$stepsize
    minimize.criterion <- control$min.criterion



    n <- dim(x)[1]

    new.s.call <- expression({
        if(c.list$one.more) {
            fv <- if(backchat) {
                      if(M>1) matrix(c.list$fit,n,M,byrow=TRUE) else c.list$fit
                  } else c.list$fit
            new.coeffs <- c.list$coeff

            if(length(slot(family, "middle")))
                eval(slot(family, "middle"))

            eta <- fv + offset
            mu <- slot(family, "inverse")(eta, extra)

            if(length(slot(family, "middle2")))
                eval(slot(family, "middle2"))

            old.crit <- new.crit
            new.crit <- 
                switch(criterion,
                    coefficients=new.coeffs,
                    tfun(mu=mu, y=y, w=w, res=FALSE, eta=eta, extra))


            if(trace && orig.stepsize==1) {
                cat("VGLM    linear loop ", iter, ": ", criterion, "= ")
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


            {
                take.half.step=(control$half.stepsizing && length(old.coeffs))&&
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

                        if(length(slot(family, "middle")))
                            eval(slot(family, "middle"))

                        fv <- xbig.save %*% new.coeffs
                        if(M > 1)
                            fv <- matrix(fv, n, M, byrow=TRUE)

                        eta <- fv + offset
                        mu <- slot(family, "inverse")(eta, extra)

                        if(length(slot(family, "middle2")))
                            eval(slot(family, "middle2"))


                        new.crit <- 
                            switch(criterion,
                                coefficients=new.coeffs,
                                tfun(mu=mu,y=y,w=w,res=FALSE,eta=eta,extra))

                        if((criterion=="coefficients") || 
                           ( minimize.criterion && new.crit < old.crit) ||
                           (!minimize.criterion && new.crit > old.crit))
                            break
                    } # of repeat

                    if(trace) 
                        cat("\n")
                    if(too.small) {
                        warning(paste("iterations terminated because",
                              "half-step sizes are very small"))
                        one.more <- FALSE
                    } else {
                        if(trace) {
                            cat("VGLM    linear loop ",
                                iter, ": ", criterion, "= ")

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
                    }
                } else {
                    one.more <- eval(control$convergence)
                }
            }
            if(exists("flush.console"))
                flush.console()

            if(!is.logical(one.more)) one.more = FALSE
            if(one.more) {
                iter <- iter + 1
                deriv.mu <- eval(slot(family, "deriv"))
                wz <- eval(slot(family, "weight"))
                if(control$checkwz)
                    wz = checkwz(wz, M=M, trace=trace, wzeps=control$wzepsilon)

                U <- vchol(wz, M=M, n=n, silent=!trace)
                tvfor <- vforsub(U, as.matrix(deriv.mu), M=M, n=n)
                z <- eta + vbacksub(U, tvfor, M=M, n=n) - offset

                c.list$z <- z
                c.list$U <- U
                if(copyxbig) c.list$xbig <- xbig.save
            }

            c.list$one.more <- one.more
            c.list$coeff <- runif(length(new.coeffs)) # 12/3/03; twist needed!
            old.coeffs <- new.coeffs
        }
        c.list
    })





    copyxbig <- FALSE    # May be overwritten in @initialize
    stepsize <- orig.stepsize
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



    if(length(slot(family, "constraints")))
        eval(slot(family, "constraints"))


    Blist <- process.constraints(constraints, x, M, specialCM=specialCM)


    ncolBlist <- unlist(lapply(Blist, ncol))
    dimB <- sum(ncolBlist)


    xbig.save <- lm2vlm.model.matrix(x, Blist, xij=control$xij)


    if(length(coefstart)) {
        eta <- if(ncol(xbig.save)>1) xbig.save %*% coefstart +
                   offset else xbig.save * coefstart + offset
        eta <- if(M > 1) matrix(eta, ncol=M, byrow=TRUE) else c(eta) 
        mu <- slot(family, "inverse")(eta, extra)
    }


    if(criterion != "coefficients") {
        tfun <- slot(family, criterion)   # family[[criterion]]
    }

    iter <- 1
    new.crit <- switch(criterion,
                      coefficients=1,
                      tfun(mu=mu, y=y, w=w, res=FALSE, eta=eta, extra))
    old.crit <- if(minimize.criterion) 10*new.crit+10 else -10*new.crit-10

    deriv.mu <- eval(slot(family, "deriv"))
    wz <- eval(slot(family, "weight"))
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
        index.vglm <- iam(NA, NA, M, both=TRUE)
        tfit <- dotFortran(name="vglmf", as.double(xbig.save), n.big, p.big,
                    backchat=as.integer(T),
                    as.integer(n), as.double(z), 
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

        bf.call <- if(is.R()) expression(vlm.wfit(xbig.save, z, Blist=NULL,
            U=U, matrix.out=FALSE, XBIG=TRUE, qr=qr.arg, xij=NULL)) else
                              expression(vlm.wfit(xbig.save, z, Blist=NULL,
            U=U, matrix.out=FALSE, XBIG=TRUE, qr=qr.arg, xij=NULL))

        while(c.list$one.more) {
            tfit <- eval(bf.call)   # fit$smooth.frame is new
    
                c.list$coeff <- tfit$coefficients 
    
            tfit$predictors <- tfit$fitted.values
    
            c.list$fit <- tfit$fitted.values
            c.list <- eval(new.s.call)
            NULL
        }
    }

    if(maxit>1 && iter>=maxit)
        warning(paste("convergence not obtained in", maxit, "iterations."))



    dn.big <- labels(xbig.save)
    xn.big <- dn.big[[2]]
    yn.big <- dn.big[[1]]

    if(length(slot(family, "fini")))
        eval(slot(family, "fini"))

    if(M>1) 
        tfit$predictors <- matrix(tfit$predictors, n, M,
                                         byrow=backchat)

    coefs <- tfit$coefficients
    asgn <- attr(xbig.save, "assign")

    names(coefs) <- xn.big

    rank <- tfit$rank
    cnames <- xn.big

    if(check.rank && rank < p.big)
        stop("vglm only handles full-rank models (currently)")

    R <- if(is.R()) tfit$qr$qr[1:p.big, 1:p.big, drop=FALSE] else {
             if(backchat) tfit$qr[1:p.big, 1:p.big, drop=FALSE] else 
                          tfit$qr$qr[1:p.big, 1:p.big, drop=FALSE]
         }
    R[lower.tri(R)] <- 0
    attributes(R) <- if(is.R()) list(dim=c(p.big, p.big),
                     dimnames=list(cnames, cnames), rank=rank) else 
                  list(dim=c(p.big, p.big),
                     dimnames=list(cnames, cnames), rank=rank, class="upper")

    effects <- tfit$effects
    neff <- rep("", n.big)
    neff[seq(p.big)] <- cnames
    names(effects) <- neff

    dim(tfit$predictors) <- c(n, M)
    dn <- labels(x)
    yn <- dn[[1]]
    xn <- dn[[2]]


    residuals <- z - tfit$predictors
    if(M==1) {
        tfit$predictors <- as.vector(tfit$predictors)
        residuals <- as.vector(residuals)
        names(residuals) <- names(tfit$predictors) <- yn
    } else {
        dimnames(residuals) <- dimnames(tfit$predictors) <-
                               list(yn, predictors.names)
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


    df.residual <- n.big - rank
    fit <- list(assign=asgn,
                coefficients=coefs,
                constraints=Blist, 
                df.residual=df.residual,
                df.total=n*M,
                effects=effects, 
                fitted.values=mu,
                offset=offset, 
                rank=rank,
                residuals=residuals,
                R=R,
                terms=Terms) # terms: This used to be done in vglm() 

    if(qr.arg) {
        fit$qr <- if(is.R()) {
            fit$qr <- tfit$qr
        } else {
            if(backchat) tfit[c("qr", "rank", "pivot", "qraux")] else tfit$qr
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

    if(length(slot(family, "last")))
        eval(slot(family, "last"))

    structure(c(fit, list(predictors=tfit$predictors,
        contrasts=attr(x, "contrasts"),
        control=control,
        crit.list=crit.list,
        extra=extra,
        family=family,
        iter=iter,
        misc=misc,
        post=post,
        rss=tfit$rss,
        x=x,
        y=y)),
        vclass=slot(family, "vfamily"))
}


# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.




process.categorical.data.vgam = expression({

 
    if(!is.matrix(y)) {
        yf = as.factor(y)
        lev = levels(yf)
        llev = length(lev)
        nn = length(yf)
        y = matrix(0, nn, llev)
        y[cbind(1:nn,as.vector(unclass(yf)))] = 1
        dimnames(y) = list(names(yf), lev)

        if(llev <= 1)
            stop("the response matrix does not have 2 or more columns")
    } else 
        nn = nrow(y)

    nvec = drop(y %*% rep(1,ncol(y)))

    if(!exists("delete.zero.colns") ||
      (exists("delete.zero.colns") && delete.zero.colns)) {
        sumy2 = as.vector(t(y) %*% rep(1,nrow(y)))
        if(any(index <- sumy2==0))
        {
            y = y[,!index,drop = FALSE]
            sumy2 = sumy2[!index]
            if(all(index) || ncol(y)<=1) stop("y matrix has 0 or 1 columns")
            warning(paste("Deleted", sum(!index),
                          "columns of the response matrix (zero counts)"))
        }
    }

    if(any(miss <- nvec==0))
    {
        warning(paste(sm <- sum(miss),
                      "observation(s) deleted (zero counts)"))
        x = x[!miss,,drop = FALSE]
        y = y[!miss,,drop = FALSE]
        nvec = nvec[!miss]
        nn = nn - sm
    }

    w = w * nvec
    nvec[nvec == 0] = 1
    y = y / nvec             # Convert to proportions

    if(!length(mustart)) {
        mustart = (nvec * y + 1/ncol(y)) / (nvec+1)   # This may be better 
        mustart = y + (1/ncol(y) - y)/nvec   # This may be wrong 
    }
})




Deviance.categorical.data.vgam <-
    function(mu, y, w, residuals = FALSE, eta, extra=NULL)
{



    if(ncol(y) == 1 || ncol(mu) == 1)
        stop("y and mu must have at least 2 columns")

    double.eps = .Machine$double.eps  # ^0.75
    devy = y
    nz = y != 0
    devy[nz] = y[nz] * log(y[nz])
    devmu = y * log(mu)
    if(any(small <- mu * (1 - mu) < double.eps)) {
        warning("fitted values close to 0 or 1")
        smu = mu[small]
        sy = y[small]
        smu = ifelse(smu < double.eps, double.eps, smu)
        devmu[small] = sy * log(smu)
    }
    devi = 2 * (devy - devmu)

    if(residuals) {
        M = if(is.matrix(eta)) ncol(eta) else 1
        if(M > 1)
            return(NULL)
        devi = devi %*% rep(1, ncol(devi))   # deviance = \sum_i devi[i]
        return(c(sign(y[,1] - mu[,1]) * sqrt(abs(devi) * w)))
    } else
        sum(w * devi)
}








sratio = function(link="logit", earg=list(),
                  parallel=FALSE, reverse=FALSE, zero=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Stopping Ratio model\n\n", 
           "Links:    ",
           namesof(if(reverse) "P[Y=j+1|Y<=j+1]" else "P[Y=j|Y>=j]", 
                   link, earg=earg),
           "\n",
           "Variance: mu[,j]*(1-mu[,j]); -mu[,j]*mu[,k]"),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero=zero ))),
    deviance=Deviance.categorical.data.vgam,
    initialize=eval(substitute(expression({
        delete.zero.colns = TRUE 
        eval(process.categorical.data.vgam)
        M = ncol(y) - 1 
        mynames = if( .reverse)
                 paste("P[Y=",2:(M+1),"|Y<=",2:(M+1),"]", sep="") else
                 paste("P[Y=",1:M,"|Y>=",1:M,"]", sep="")
        predictors.names = namesof(mynames, .link, short=TRUE, earg= .earg)
        y.names = paste("mu", 1:(M+1), sep="")
        extra$mymat = if( .reverse ) tapplymat1(y, "cumsum") else
                      tapplymat1(y[,ncol(y):1], "cumsum")[,ncol(y):1]
        if(length(dimnames(y)))
            extra$dimnamesy2 = dimnames(y)[[2]]
    }), list( .earg=earg, .link=link, .reverse=reverse ))),
    inverse=eval(substitute( function(eta, extra=NULL) {
        if(!is.matrix(eta))
            eta = as.matrix(eta)
        fv.matrix =
        if( .reverse ) {
            M = ncol(eta)
            djr = eta2theta(eta, .link, earg= .earg )
            temp = tapplymat1(1-djr[,M:1], "cumprod")[,M:1]
            cbind(1,djr) * cbind(temp,1)
        } else {
            dj = eta2theta(eta, .link, earg= .earg )
            temp = tapplymat1(1-dj, "cumprod")
            cbind(dj,1) * cbind(1, temp)
        }
        if(length(extra$dimnamesy2))
            dimnames(fv.matrix) = list(dimnames(eta)[[1]], extra$dimnamesy2)
        fv.matrix
    }, list( .earg=earg, .link=link, .reverse=reverse) )),
    last=eval(substitute(expression({
        misc$link = rep( .link, length=M)
        names(misc$link) = mynames

        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg

        misc$parameters = mynames
        misc$reverse = .reverse
        extra = list()   # kill what was used 
    }), list( .earg=earg, .link=link, .reverse=reverse ))),
    link=eval(substitute( function(mu, extra=NULL) {
        cump = tapplymat1(mu, "cumsum")
        if( .reverse ) {
            djr = mu[,-1] / cump[,-1]
            theta2eta(djr, .link, earg= .earg )
        } else {
            M = ncol(mu) - 1
            dj = if(M==1) mu[,1] else mu[,1:M]/(1-cbind(0,cump[,1:(M-1)]))
            theta2eta(dj, .link, earg= .earg )
        }
    }, list( .earg=earg, .link=link, .reverse=reverse) )),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * y * log(mu)),
    vfamily=c("sratio", "vcategorical"),
    deriv=eval(substitute(expression({
        if(!length(extra$mymat)) {
            extra$mymat = if( .reverse ) tapplymat1(y, "cumsum") else
                          tapplymat1(y[,ncol(y):1], "cumsum")[,ncol(y):1]
        }
        if( .reverse ) {
            djr = eta2theta(eta, .link, earg= .earg )
            Mp1 = ncol(extra$mymat)
            w * (y[,-1]/djr - extra$mymat[,-Mp1]/(1-djr)) *
              dtheta.deta(djr, .link, earg= .earg )
        } else {
            dj = eta2theta(eta, .link, earg= .earg )
            w * (y[,-ncol(y)]/dj - extra$mymat[,-1]/(1-dj)) *
              dtheta.deta(dj, .link, earg= .earg )
        }
    }), list( .earg=earg, .link=link, .reverse=reverse) )),
    weight= eval(substitute(expression({
        if( .reverse ) {
            cump = tapplymat1(mu, "cumsum")
            ddjr.deta = dtheta.deta(djr, .link, earg= .earg )
            wz = w * ddjr.deta^2 * (mu[,-1]/djr^2 + cump[,1:M]/(1-djr)^2)
        } else {
            ccump = tapplymat1(mu[,ncol(mu):1], "cumsum")[,ncol(mu):1]
            ddj.deta = dtheta.deta(dj, .link, earg= .earg )
            wz = w * ddj.deta^2 * (mu[,1:M]/dj^2 + ccump[,-1]/(1-dj)^2)
        }

        wz
    }), list( .earg=earg, .link=link, .reverse=reverse ))))
}




cratio = function(link="logit", earg=list(),
                  parallel=FALSE, reverse=FALSE, zero=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Continuation Ratio model\n\n", 
           "Links:    ",
           namesof(if(reverse) "P[Y<j+1|Y<=j+1]" else "P[Y>j|Y>=j]", 
                   link, earg=earg),
           "\n",
           "Variance: mu[,j]*(1-mu[,j]); -mu[,j]*mu[,k]"),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero=zero ))),
    deviance=Deviance.categorical.data.vgam,
    initialize=eval(substitute(expression({
        delete.zero.colns = TRUE 
        eval(process.categorical.data.vgam)
        M = ncol(y) - 1 
        mynames = if( .reverse )
            paste("P[Y<",2:(M+1),"|Y<=",2:(M+1),"]", sep="") else
            paste("P[Y>",1:M,"|Y>=",1:M,"]", sep="")
        predictors.names = namesof(mynames, .link, short=TRUE, earg= .earg)
        y.names = paste("mu", 1:(M+1), sep="")
        extra$mymat = if( .reverse ) tapplymat1(y, "cumsum") else
                      tapplymat1(y[,ncol(y):1], "cumsum")[,ncol(y):1]
        if(length(dimnames(y)))
            extra$dimnamesy2 = dimnames(y)[[2]]
    }), list( .earg=earg, .link=link, .reverse=reverse ))),
    inverse=eval(substitute( function(eta, extra=NULL) {
        if(!is.matrix(eta))
            eta = as.matrix(eta)
        fv.matrix =
        if( .reverse ) {
            M = ncol(eta)
            djrs = eta2theta(eta, .link, earg= .earg )
            temp = tapplymat1(djrs[,M:1], "cumprod")[,M:1]
            cbind(1,1-djrs) * cbind(temp,1)
        } else {
            djs = eta2theta(eta, .link, earg= .earg )
            temp = tapplymat1(djs, "cumprod")
            cbind(1-djs,1) * cbind(1, temp)
        }
        if(length(extra$dimnamesy2))
            dimnames(fv.matrix) = list(dimnames(eta)[[1]], extra$dimnamesy2)
        fv.matrix
    }, list( .earg=earg, .link=link, .reverse=reverse) )),
    last=eval(substitute(expression({
        misc$link = rep( .link, length=M)
        names(misc$link) = mynames
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg
        misc$parameters = mynames
        misc$reverse = .reverse
        extra = list()   # kill what was used 
    }), list( .earg=earg, .link=link, .reverse=reverse ))),
    link=eval(substitute( function(mu, extra=NULL) {
        cump = tapplymat1(mu, "cumsum")
        if( .reverse ) {
            djrs = 1 - mu[,-1] / cump[,-1]
            theta2eta(djrs, .link, earg= .earg )
        } else {
            M = ncol(mu) - 1
            djs = if(M==1) 1-mu[,1] else 1-mu[,1:M]/(1-cbind(0,cump[,1:(M-1)]))
            theta2eta(djs, .link, earg= .earg )
        }
    }, list( .earg=earg, .link=link, .reverse=reverse) )),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * y * log(mu)), 
    vfamily=c("cratio", "vcategorical"),
    deriv=eval(substitute(expression({
        if(!length(extra$mymat)) {
            extra$mymat = if( .reverse ) tapplymat1(y, "cumsum") else
                          tapplymat1(y[,ncol(y):1], "cumsum")[,ncol(y):1]
        }
        if( .reverse ) {
            djrs = eta2theta(eta, .link, earg= .earg )
            Mp1 = ncol(extra$mymat)
            -w * (y[,-1]/(1-djrs) - extra$mymat[,-Mp1]/djrs) *
              dtheta.deta(djrs, .link, earg= .earg )
        } else {
            djs = eta2theta(eta, .link, earg= .earg )
            -w * (y[,-ncol(y)]/(1-djs) - extra$mymat[,-1]/djs) *
              dtheta.deta(djs, .link, earg= .earg )
        }
    }), list( .earg=earg, .link=link, .reverse=reverse) )),
    weight= eval(substitute(expression({
        if( .reverse ) {
            cump = tapplymat1(mu, "cumsum")
            ddjrs.deta = dtheta.deta(djrs, .link, earg= .earg )
            wz = w * ddjrs.deta^2 * (mu[,-1]/(1-djrs)^2 + cump[,1:M]/djrs^2)
        } else {
            ccump = tapplymat1(mu[,ncol(mu):1], "cumsum")[,ncol(mu):1]
            ddjs.deta = dtheta.deta(djs, .link, earg= .earg )
            wz = w * ddjs.deta^2 * (mu[,1:M]/(1-djs)^2 + ccump[,-1]/djs^2)
        }

        wz
    }), list( .earg=earg, .link=link, .reverse=reverse ))))
}




vglm.multinomial.deviance.control = function(maxit=21, panic=FALSE, ...)
{
    if(maxit < 1) {
        warning("bad value of maxit; using 21 instead")
        maxit = 21
    }
    list(maxit=maxit, panic=as.logical(panic)[1])
}

vglm.multinomial.control = function(maxit=21, panic=FALSE, 
      criterion=c("aic1", "aic2", names( .min.criterion.VGAM )), ...)
{

    if(mode(criterion) != "character" && mode(criterion) != "name")
        criterion = as.character(substitute(criterion))
    criterion = match.arg(criterion,
        c("aic1", "aic2", names( .min.criterion.VGAM )))[1]

    if(maxit < 1) {
        warning("bad value of maxit; using 21 instead")
        maxit = 21
    }
    list(maxit=maxit, panic=as.logical(panic)[1],
         criterion=criterion,
         min.criterion=c("aic1"=FALSE, "aic2"=TRUE, .min.criterion.VGAM))
}


vglm.vcategorical.control = function(maxit=30, trace=FALSE, panic=TRUE, ...)
{
    if(maxit < 1) {
        warning("bad value of maxit; using 200 instead")
        maxit = 200
    }
    list(maxit=maxit, trace=as.logical(trace)[1], panic=as.logical(panic)[1])
}



multinomial = function(zero=NULL, parallel=FALSE, nointercept=NULL)
{

    new("vglmff",
    blurb=c("Multinomial logit model\n\n", 
           "Links:    log(mu[,j]/mu[,M+1])\n",
           "Variance: mu[,j]*(1-mu[,j]); -mu[,j]*mu[,k]"),
    constraints=eval(substitute(expression({





        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints,
                               intercept.apply=FALSE)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
        constraints = cm.nointercept.vgam(constraints, x, .nointercept, M)
    }), list( .parallel=parallel, .zero=zero, .nointercept=nointercept ))),
    deviance=Deviance.categorical.data.vgam,
    initialize=expression({
        delete.zero.colns = TRUE 
        eval(process.categorical.data.vgam)
        M = ncol(y)-1
        predictors.names = paste("log(mu[,",1:M,"]/mu[,",M+1,"])",sep="")
        y.names = paste("mu", 1:(M+1), sep="")
    }),
    inverse=function(eta, extra=NULL) {
        if(any(is.na(eta)))
            warning("there are NAs in eta in slot inverse")
        phat = cbind(exp(eta), 1) 
        ans = phat / as.vector(phat %*% rep(1, ncol(phat)))
        if(any(is.na(ans)))
            warning("there are NAs here in slot inverse")
        ans
    }, 
    last=expression({
        misc$link = "mlogit"
        misc$earg = list(mlogit = list()) # vector("list", M)

        dy = dimnames(y)
        if(!is.null(dy[[2]]))
            dimnames(fit$fitted.values) = dy
    }),
    link=function(mu, extra=NULL)
        log(mu[,-ncol(mu)]/mu[,ncol(mu)]),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * y * log(mu)), 
    vfamily=c("multinomial", "vcategorical"),
    deriv=expression({
        w * (y[,-ncol(y)] - mu[,-ncol(y)])
    }),
    weight= expression({
        tiny = (mu < .Machine$double.eps^.5) | 
                (mu > 1 - .Machine$double.eps^.5)

        if(M==1) wz = mu[,1] * (1-mu[,1]) else {
            index = iam(NA, NA, M, both=TRUE, diag=TRUE)
            wz = -mu[,index$row] * mu[,index$col]
            wz[,1:M] = wz[,1:M] + mu[,1:M]
        }

        atiny = (tiny %*% rep(1, ncol(mu))) > 0 # apply(tiny, 1, any) 
        if(any(atiny)) {
            if(M==1) wz[atiny] = wz[atiny] * (1 + .Machine$double.eps^0.5) +
                           .Machine$double.eps else 
            wz[atiny,1:M] = wz[atiny,1:M] * (1 + .Machine$double.eps^0.5) +
                             .Machine$double.eps
        }
        w * wz
    }))
}



cumulative = function(link="logit", earg = list(),
                      parallel=FALSE, reverse=FALSE, 
                      mv=FALSE,
                      intercept.apply = FALSE)
{
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.logical(mv) || length(mv)!=1) stop("\"mv\" must be a single logical")
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=if( mv ) c(paste("Multivariate cumulative", link, "model\n\n"), 
           "Links:   ",
           namesof(if(reverse) "P[Y1>=j+1]" else "P[Y1<=j]", link, earg=earg),
           ", ...") else
           c(paste("Cumulative", link, "model\n\n"),
           "Links:   ",
           namesof(if(reverse) "P[Y>=j+1]" else "P[Y<=j]", link, earg=earg)),
    constraints=eval(substitute(expression({
        if( .mv ) {
            if( !length(constraints) ) {
                Llevels = extra$Llevels
                NOS = extra$NOS
                Hk.matrix = kronecker(diag(NOS), matrix(1,Llevels-1,1))
                constraints = cm.vgam(Hk.matrix, x, .parallel, constraints,
                                      intercept.apply = .intercept.apply)
            }
        } else {
            constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints,
                                  intercept.apply = .intercept.apply)
        }
    }), list( .parallel=parallel, .mv=mv, .intercept.apply=intercept.apply ))),
    deviance=eval(substitute(
        function(mu, y, w, residuals=FALSE, eta, extra=NULL) {
        answer =
        if( .mv ) {
            totdev = 0
            NOS = extra$NOS
            Llevels = extra$Llevels
            for(iii in 1:NOS) {
                cindex = (iii-1)*(Llevels-1) + 1:(Llevels-1)
                aindex = (iii-1)*(Llevels) + 1:(Llevels)
                totdev = totdev + Deviance.categorical.data.vgam(
                         mu=mu[,aindex,drop=FALSE],
                         y=y[,aindex,drop=FALSE], w=w, residuals=residuals,
                         eta=eta[,cindex,drop=FALSE], extra=extra)
            }
            totdev
        } else {
            Deviance.categorical.data.vgam(mu=mu, y=y, w=w, residuals=residuals,
                                           eta=eta, extra=extra)
        }
        answer
    }, list( .earg=earg, .link=link, .mv = mv ) )),
    initialize=eval(substitute(expression({
        extra$mv = .mv
        if( .mv ) {
            checkCut(y)  # Check the input; stops if there is an error.
            if(any(w != 1) || ncol(cbind(w)) != 1)
                stop("the 'weights' argument must be a vector of all ones")
            Llevels = max(y)
            delete.zero.colns = FALSE 
            orig.y = cbind(y) # Convert y into a matrix if necessary
            NOS = ncol(cbind(orig.y))
            use.y = use.mustart = NULL
            for(iii in 1:NOS) {
                y = as.factor(orig.y[,iii])
                eval(process.categorical.data.vgam)
                use.y = cbind(use.y, y)
                use.mustart = cbind(use.mustart, mustart)
            }
            mustart = use.mustart
            y = use.y  # n x (Llevels*NOS)
            M = NOS * (Llevels-1)
            mynames = y.names = NULL
            for(iii in 1:NOS) {
                Y.names = paste("Y", iii, sep="")
                mu.names = paste("mu", iii, ".", sep="")
                mynames = c(mynames, if( .reverse )
                    paste("P[",Y.names,">=",2:Llevels,"]", sep="") else
                    paste("P[",Y.names,"<=",1:(Llevels-1),"]", sep=""))
                y.names = c(y.names, paste(mu.names, 1:Llevels, sep=""))
            }
            predictors.names = namesof(mynames, .link, short=TRUE, earg= .earg)
            extra$NOS = NOS
            extra$Llevels = Llevels
        } else {
            delete.zero.colns=TRUE # Cannot have FALSE since then prob(Y=jay)=0
            eval(process.categorical.data.vgam)
            M = ncol(y)-1
            mynames = if( .reverse ) paste("P[Y>=",2:(1+M),"]", sep="") else
                paste("P[Y<=",1:M,"]", sep="")
            predictors.names = namesof(mynames, .link, short=TRUE, earg= .earg)
            y.names = paste("mu", 1:(M+1), sep="")

            if(length(dimnames(y)))
                extra$dimnamesy2 = dimnames(y)[[2]]
        }
    }), list( .link=link, .reverse=reverse, .mv = mv, .earg = earg ))),
    inverse=eval(substitute( function(eta, extra=NULL) {
        answer =
        if( .mv ) {
            NOS = extra$NOS
            Llevels = extra$Llevels
            fv.matrix = matrix(0, nrow(eta), NOS*Llevels)
            for(iii in 1:NOS) {
                cindex = (iii-1)*(Llevels-1) + 1:(Llevels-1)
                aindex = (iii-1)*(Llevels) + 1:(Llevels)
                if( .reverse ) {
                    ccump = cbind(1,eta2theta(eta[,cindex,drop=FALSE], .link,
                                              earg= .earg))
                    fv.matrix[,aindex] =
                        cbind(-tapplymat1(ccump, "diff"), ccump[,ncol(ccump)])
                } else {
                    cump = cbind(eta2theta(eta[,cindex,drop=FALSE], .link,
                                           earg= .earg), 1)
                    fv.matrix[,aindex] =
                        cbind(cump[,1], tapplymat1(cump, "diff"))
                }
            }
            fv.matrix
        } else {
            fv.matrix =
            if( .reverse ) {
                ccump = cbind(1, eta2theta(eta, .link, earg= .earg))
                cbind(-tapplymat1(ccump, "diff"), ccump[,ncol(ccump)])
            } else {
                cump = cbind(eta2theta(eta, .link, earg= .earg), 1)
                cbind(cump[,1], tapplymat1(cump, "diff"))
            }
            if(length(extra$dimnamesy2))
                dimnames(fv.matrix) = list(dimnames(eta)[[1]], extra$dimnamesy2)
            fv.matrix
        }
        answer
    }, list( .link=link, .reverse=reverse, .earg= earg, .mv = mv ))),
    last=eval(substitute(expression({
        if( .mv ) {
            misc$link = .link
            misc$earg = list( .earg )
        } else {
            misc$link = rep( .link, length=M)
            names(misc$link) = mynames
            misc$earg = vector("list", M)
            names(misc$earg) = names(misc$link)
            for(ii in 1:M) misc$earg[[ii]] = .earg
        }

        misc$parameters = mynames
        misc$reverse = .reverse
        misc$parallel = .parallel
        misc$mv = .mv
    }), list( .link=link, .reverse=reverse, .parallel=parallel,
              .mv = mv, .earg= earg ))),
    link=eval(substitute( function(mu, extra=NULL) {
        answer = 
        if( .mv ) {
            NOS = extra$NOS
            Llevels = extra$Llevels
            eta.matrix = matrix(0, nrow(mu), NOS*(Llevels-1))
            for(iii in 1:NOS) {
                cindex = (iii-1)*(Llevels-1) + 1:(Llevels-1)
                aindex = (iii-1)*(Llevels) + 1:(Llevels)
                cump = tapplymat1(as.matrix(mu[,aindex]), "cumsum")
                eta.matrix[,cindex] =
                    theta2eta(if( .reverse) 1-cump[,1:(Llevels-1)] else
                          cump[,1:(Llevels-1)], .link, earg= .earg)
            }
            eta.matrix
        } else {
            cump = tapplymat1(as.matrix(mu), "cumsum")
            M = ncol(as.matrix(mu)) - 1
            theta2eta(if( .reverse ) 1-cump[,1:M] else cump[,1:M], .link,
                      earg= .earg)
        }
        answer
    }, list( .link=link, .reverse=reverse, .earg=earg, .mv=mv ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
       sum(w * y * log(mu)), 
    vfamily=c("cumulative", "vcategorical"),
    deriv=eval(substitute(expression({
        mu.use = pmax(mu, .Machine$double.eps * 1.0e-0)
        deriv.answer = 
        if( .mv ) {
            NOS = extra$NOS
            Llevels = extra$Llevels
            dcump.deta = answer.matrix = matrix(0, n, NOS * (Llevels-1))
            for(iii in 1:NOS) {
                cindex = (iii-1)*(Llevels-1) + 1:(Llevels-1)
                aindex = (iii-1)*(Llevels)   + 1:(Llevels-1)
                cump = eta2theta(eta[,cindex,drop=FALSE], .link, earg= .earg)
                dcump.deta[,cindex] = dtheta.deta(cump, .link, earg= .earg)
                answer.matrix[,cindex] =
                    (y[,aindex,drop=FALSE]/mu.use[,aindex,drop=FALSE] -
                   y[,1+aindex,drop=FALSE]/mu.use[,1+aindex,drop=FALSE])
            }
            (if( .reverse) -w  else w) * dcump.deta * answer.matrix 
        } else {
            cump = eta2theta(eta, .link, earg= .earg)
            dcump.deta = dtheta.deta(cump, .link, earg= .earg)
            (if( .reverse) -w  else w) * dcump.deta *
                (y[,-(M+1)]/mu.use[,-(M+1)] - y[,-1]/mu.use[,-1])
        }
        deriv.answer
    }), list( .link=link, .reverse=reverse, .earg= earg, .mv=mv ))),
    weight= eval(substitute(expression({
        if( .mv ) {
            NOS = extra$NOS
            Llevels = extra$Llevels
            wz = matrix(0, n, NOS*(Llevels-1)) # Diagonal elts only for a start
            for(iii in 1:NOS) {
                cindex = (iii-1)*(Llevels-1) + 1:(Llevels-1)
                aindex = (iii-1)*(Llevels)   + 1:(Llevels-1)
                wz[,cindex] = w * dcump.deta[,cindex,drop=FALSE]^2 *
                (1/mu.use[,aindex,drop=FALSE] + 1/mu.use[,1+aindex,drop=FALSE])
            }
            if(Llevels-1 > 1) {
                iii = 1
                oindex = (iii-1)*(Llevels-1) + 1:(Llevels-2)
                wz = cbind(wz, -w * dcump.deta[,oindex] * dcump.deta[,1+oindex])


                if(NOS > 1) {
                    cptrwz = ncol(wz)  # Like a pointer
                    wz = cbind(wz, matrix(0, nrow(wz), (NOS-1) * (Llevels-1)))
                    for(iii in 2:NOS) {
                        oindex = (iii-1)*(Llevels-1) + 1:(Llevels-2)
                        wz[,cptrwz + 1 + (1:(Llevels-2))] =
                              -w * dcump.deta[,oindex] *
                                   dcump.deta[,1+oindex]
                        cptrwz = cptrwz + Llevels - 1 # Move it along a bit
                    }
                }



            }
        } else {
            wz = w * dcump.deta^2 * (1/mu.use[,1:M] + 1/mu.use[,-1])
            if(M > 1)
                wz = cbind(wz, -w * dcump.deta[,-M] *
                            dcump.deta[,2:M] / mu.use[,2:M])
        }
        wz
    }), list( .earg=earg, .link=link, .mv=mv ))))
}



acat = function(link="loge", earg = list(),
                parallel=FALSE, reverse=FALSE, zero=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Adjacent-categories model\n\n",
           "Links:    ",
         namesof(if(reverse) "P[Y=j]/P[Y=j+1]" else "P[Y=j+1]/P[Y=j]",
                    link, earg=earg),
           "\n",
           "Variance: mu[,j]*(1-mu[,j]); -mu[,j]*mu[,k]"),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero=zero ))),

    deviance=Deviance.categorical.data.vgam,
    initialize=eval(substitute(expression({
        delete.zero.colns = TRUE 
        eval(process.categorical.data.vgam)
        M = ncol(y) - 1
        mynames = if( .reverse )
            paste("P[Y=",1:M,"]/P[Y=",2:(M+1),"]", sep="") else
            paste("P[Y=",2:(M+1),"]/P[Y=",1:M,"]", sep="")
        predictors.names = namesof(mynames, .link, short=TRUE, earg= .earg)
        y.names = paste("mu", 1:(M+1), sep="")
        if(length(dimnames(y)))
            extra$dimnamesy2 = dimnames(y)[[2]]
    }), list( .earg=earg, .link=link, .reverse=reverse ))),
    inverse=eval(substitute( function(eta, extra=NULL) {
        if(!is.matrix(eta))
            eta = as.matrix(eta)
        M = ncol(eta)
        fv.matrix =
        if( .reverse ) {
            zetar = eta2theta(eta, .link, earg= .earg )
            temp = tapplymat1(zetar[,M:1], "cumprod")[,M:1,drop=FALSE]
            cbind(temp,1) / drop(1 + temp %*% rep(1,ncol(temp)))
        } else {
            zeta = eta2theta(eta, .link, earg= .earg )
            temp = tapplymat1(zeta, "cumprod")
            cbind(1,temp) / drop(1 + temp %*% rep(1,ncol(temp)))
        }
        if(length(extra$dimnamesy2))
            dimnames(fv.matrix) = list(dimnames(eta)[[1]], extra$dimnamesy2)
        fv.matrix
    }, list( .earg=earg, .link=link, .reverse=reverse) )),
    last=eval(substitute(expression({
        misc$link = rep( .link, length=M)
        names(misc$link) = mynames

        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg

        misc$parameters = mynames
        misc$reverse = .reverse
    }), list( .earg=earg, .link=link, .reverse=reverse ))),
    link=eval(substitute( function(mu, extra=NULL) {
        M = ncol(mu) - 1
        theta2eta(if( .reverse ) mu[,1:M]/mu[,-1] else mu[,-1]/mu[,1:M],
                  .link, earg= .earg )
    }, list( .earg=earg, .link=link, .reverse=reverse) )),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * y * log(mu)), 
    vfamily=c("acat", "vcategorical"),
    deriv=eval(substitute(expression({
        zeta = eta2theta(eta, .link, earg= .earg )    # May be zetar
        d1 = acat.deriv(zeta, M=M, n=n, reverse=.reverse)
        score = attr(d1, "gradient") / d1
        dzeta.deta = dtheta.deta(zeta, .link, earg= .earg )
        if( .reverse ) {
            cumy = tapplymat1(y, "cumsum")
            w * dzeta.deta * (cumy[,1:M]/zeta - score)
        } else {
            ccumy = tapplymat1(y[,ncol(y):1], "cumsum")[,ncol(y):1]
            w * dzeta.deta * (ccumy[,-1]/zeta - score)
        }
    }), list( .earg=earg, .link=link, .reverse=reverse) )),
    weight= eval(substitute(expression({
        wz = matrix(as.numeric(NA), n, dimm(M)) 

        hess = attr(d1, "hessian") / d1

        if(M>1)
        for(j in 1:(M-1)) 
            for(k in (j+1):M) 
                wz[,iam(j,k,M)] = (hess[,j,k] - score[,j] * score[,k]) * 
                    dzeta.deta[,j] * dzeta.deta[,k]
        if( .reverse ) {
            cump = tapplymat1(mu, "cumsum")
            wz[,1:M] = (cump[,1:M]/zeta^2 - score^2) * dzeta.deta^2 
        } else {
            ccump = tapplymat1(mu[,ncol(mu):1], "cumsum")[,ncol(mu):1]
            wz[,1:M] = (ccump[,-1]/zeta^2 - score^2) * dzeta.deta^2 
        }
        w * wz
    }), list( .earg=earg, .link=link, .reverse=reverse ))))
}

acat.deriv = function(zeta, reverse, M, n)
{

    alltxt = NULL
    for(i in 1:M) {
        index = if(reverse) i:M else 1:i
        vars = paste("zeta", index, sep="")
        txt = paste(vars, collapse="*")
        alltxt = c(alltxt, txt) 
    }
    alltxt = paste(alltxt, collapse=" + ")
    alltxt = paste(" ~ 1 +", alltxt)
    txt = as.formula(alltxt) 

    allvars = paste("zeta", 1:M, sep="")
    d1 = deriv3(txt, allvars, hessian=TRUE)  # deriv3() computes the Hessian

    zeta = as.matrix(zeta)
    for(i in 1:M)
        assign(paste("zeta", i, sep=""), zeta[,i])

    ans = eval(d1)
    ans
}




brat = function(refgp="last",
                 refvalue = 1,
                 init.alpha = 1)
{
    if(!is.Numeric(init.alpha, posit=TRUE))
        stop("\"init.alpha\" must contain positive values only")
    if(!is.Numeric(refvalue, allow=1, posit=TRUE))
        stop("\"refvalue\" must be a single positive value")
    if(!is.character(refgp) &&
       !is.Numeric(refgp, allow=1, integer=TRUE, posit=TRUE))
        stop("\"refgp\" must be a single positive integer")
    new("vglmff",
    blurb=c(paste("Bradley-Terry model (without ties)\n\n"), 
           "Links:   ",
           namesof("alpha's", "loge")),
    initialize=eval(substitute(expression({
        are.ties = attr(y, "are.ties")  # If Brat() was used
        if(is.logical(are.ties) && are.ties)
            stop("use bratt(), not brat(), when there are ties")

        try.index = 1:400
        M = (1:length(try.index))[(try.index+1)*(try.index) == ncol(y)]
        if(!is.finite(M)) stop("can't determine M")
        init.alpha = matrix( rep( .init.alpha, len=M), n, M, byrow=TRUE)
        etastart = matrix(theta2eta(init.alpha, "loge", earg=list()), n, M, byrow=TRUE)
        refgp = .refgp
        if(!intercept.only)
            warning("this function only works with intercept-only models")
        extra$ybrat.indices = .brat.indices(NCo=M+1, are.ties=FALSE)
        uindex = if( .refgp =="last") 1:M else (1:(M+1))[-( .refgp ) ]

        predictors.names=namesof(paste("alpha",uindex,sep=""),"loge",short=TRUE)
    }), list( .refgp=refgp, .init.alpha=init.alpha ))),
    inverse=eval(substitute( function(eta, extra=NULL) {
        probs = NULL
        eta = as.matrix(eta)   # in case M=1
        for(ii in 1:nrow(eta)) {
            alpha = .brat.alpha(eta2theta(eta[ii,], "loge", earg=list()),
                                .refvalue, .refgp)
            alpha1 = alpha[extra$ybrat.indices[,"rindex"]]
            alpha2 = alpha[extra$ybrat.indices[,"cindex"]]
            probs = rbind(probs, alpha1/(alpha1+alpha2))
        }
        dimnames(probs) = dimnames(eta)
        probs
    }, list( .refgp=refgp, .refvalue=refvalue) )),
    last=eval(substitute(expression({
        misc$link = rep( "loge", length=M)
        names(misc$link) = paste("alpha",uindex,sep="")
        misc$refgp = .refgp
        misc$refvalue = .refvalue
    }), list( .refgp=refgp, .refvalue=refvalue ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * y * log(mu)), 
    vfamily=c("brat"),
    deriv=eval(substitute(expression({
        ans = NULL
        uindex = if( .refgp =="last") 1:M else (1:(M+1))[-( .refgp ) ]
        eta = as.matrix(eta)   # in case M=1
        for(ii in 1:nrow(eta)) {
            alpha = .brat.alpha(eta2theta(eta[ii,], "loge", earg=list()),
                                .refvalue, .refgp)
            ymat = InverseBrat(y[ii,], NCo=M+1, diag=0)
            answer = rep(0, len=M)
            for(aa in 1:(M+1)) {
                answer = answer + (1-(aa==uindex)) *
                (ymat[uindex,aa] * alpha[aa] - ymat[aa,uindex] *
                alpha[uindex]) / (alpha[aa] + alpha[uindex])
            }
            ans = rbind(ans, w[ii] * answer)
        }
        dimnames(ans) = dimnames(eta)
        ans
    }), list( .refvalue=refvalue, .refgp=refgp) )),
    weight= eval(substitute(expression({
        wz = matrix(0, n, dimm(M))
        for(ii in 1:nrow(eta)) {
            alpha = .brat.alpha(eta2theta(eta[ii,], "loge", earg=list()),
                                .refvalue, .refgp)
            ymat = InverseBrat(y[ii,], NCo=M+1, diag=0)
            for(aa in 1:(M+1)) {
                wz[ii,1:M] = wz[ii,1:M] + (1-(aa==uindex)) *
                (ymat[aa,uindex] + ymat[uindex,aa]) * alpha[aa] *
                alpha[uindex] / (alpha[aa] + alpha[uindex])^2
            }
            if(M > 1) {
                ind5 = iam(1,1,M, both=TRUE, diag=FALSE)
                wz[ii,(M+1):ncol(wz)] =
                  -(ymat[cbind(uindex[ind5$row],uindex[ind5$col])] +
                    ymat[cbind(uindex[ind5$col],uindex[ind5$row])]) *
                    alpha[uindex[ind5$col]] * alpha[uindex[ind5$row]] /
                    (alpha[uindex[ind5$row]] + alpha[uindex[ind5$col]])^2
            }
        }
        wz = wz * w
        wz
    }), list( .refvalue=refvalue, .refgp=refgp ))))
}




bratt = function(refgp="last",
                  refvalue = 1,
                  init.alpha = 1,
                  i0 = 0.01)
{
    if(!is.Numeric(i0, allow=1, positi=TRUE))
        stop("\"i0\" must be a single positive value")
    if(!is.Numeric(init.alpha, positi=TRUE))
        stop("\"init.alpha\" must contain positive values only")
    if(!is.Numeric(refvalue, allow=1, positi=TRUE))
        stop("\"refvalue\" must be a single positive value")
    if(!is.character(refgp) && 
       !is.Numeric(refgp, allow=1, integer=TRUE, positi=TRUE))
        stop("\"refgp\" must be a single positive integer")
    new("vglmff",
    blurb=c(paste("Bradley-Terry model (with ties)\n\n"), 
           "Links:   ",
           namesof("alpha's", "loge"), ", log(alpha0)"),
    initialize=eval(substitute(expression({
        try.index = 1:400
        M = (1:length(try.index))[(try.index*(try.index-1)) == ncol(y)]
        if(!is.Numeric(M, allow=1, integ=TRUE)) stop("can't determine M")
        NCo = M  # number of contestants

        are.ties = attr(y, "are.ties")  # If Brat() was used
        if(is.logical(are.ties)) {
            if(!are.ties)
                stop("use brat(), not bratt(), when there are no ties")
            ties = attr(y, "ties")
        } else {
            are.ties = FALSE
            ties = 0 * y
        }

        init.alpha = rep( .init.alpha, len=NCo-1)
        ialpha0 = .i0
        etastart = cbind(matrix(theta2eta(init.alpha, "loge"),
                                n, NCo-1, byrow=TRUE),
                         theta2eta( rep(ialpha0, len=n), "loge"))
        refgp = .refgp
        if(!intercept.only)
            warning("this function only works with intercept-only models")
        extra$ties = ties  # Flat (1-row) matrix
        extra$ybrat.indices = .brat.indices(NCo=NCo, are.ties=FALSE)
        extra$tbrat.indices = .brat.indices(NCo=NCo, are.ties=TRUE) # unused
        extra$dnties = dimnames(ties)
        uindex = if(refgp=="last") 1:(NCo-1) else (1:(NCo))[-refgp ]

        predictors.names=c(
            namesof(paste("alpha",uindex,sep=""),"loge",short=TRUE),
            namesof("alpha0", "loge", short=TRUE))
    }), list( .refgp=refgp,
             .i0 = i0,
             .init.alpha=init.alpha ))),
    inverse=eval(substitute( function(eta, extra=NULL) {
        probs = qprobs = NULL
        M = ncol(eta)
        for(ii in 1:nrow(eta)) {
            alpha = .brat.alpha(eta2theta(eta[ii,-M], "loge"), .refvalue, .refgp)
            alpha0 = eta2theta(eta[ii,M], "loge")
            alpha1 = alpha[extra$ybrat.indices[,"rindex"]]
            alpha2 = alpha[extra$ybrat.indices[,"cindex"]]
            probs = rbind(probs, alpha1/(alpha1+alpha2+alpha0)) #
            qprobs = rbind(qprobs, alpha0/(alpha1+alpha2+alpha0)) #
        }
        if(length(extra$dnties))
            dimnames(qprobs) = extra$dnties
        attr(probs, "probtie") = qprobs
        probs
    }, list( .refgp=refgp, .refvalue=refvalue) )),
    last=eval(substitute(expression({
        misc$link = rep( "loge", length=M)
        names(misc$link) = c(paste("alpha",uindex,sep=""), "alpha0")
        misc$refgp = .refgp
        misc$refvalue = .refvalue
        misc$alpha  = alpha   # zz; last one saved (ok since n=1)
        misc$alpha0 = alpha0  # zz; last one saved (ok since n=1)
    }), list( .refgp=refgp, .refvalue=refvalue ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (y * log(mu) + 0.5 * extra$ties * log(attr(mu, "probtie")))),
    vfamily=c("bratt"),
    deriv=eval(substitute(expression({
        ans = NULL
        ties = extra$ties
        NCo = M
        uindex = if( .refgp =="last") 1:(M-1) else (1:(M))[-( .refgp )]
        eta = as.matrix(eta)
        for(ii in 1:nrow(eta)) {
            alpha = .brat.alpha(eta2theta(eta[ii,-M], "loge"), .refvalue, .refgp)
            alpha0 = eta2theta(eta[ii,M], "loge") # M == ncol(eta)
            ymat = InverseBrat(y[ii,], NCo=M, diag=0)
            tmat = InverseBrat(ties[ii,], NCo=M, diag=0)
            answer = rep(0, len=NCo-1) # deriv wrt eta[-M]
            for(aa in 1:NCo) {
                Daj = alpha[aa] + alpha[uindex] + alpha0
                pja = alpha[uindex] / Daj
                answer = answer + alpha[uindex] *
                         (-ymat[aa,uindex] + ymat[uindex,aa]*(1-pja)/pja -
                         tmat[uindex,aa]) / Daj
            }
            deriv0 = 0 # deriv wrt eta[M]
            for(aa in 1:(NCo-1)) 
                for(bb in (aa+1):NCo) {
                        Dab = alpha[aa] + alpha[bb] + alpha0
                        qab = alpha0 / Dab
                        deriv0 = deriv0 + alpha0 *
                                 (-ymat[aa,bb] - ymat[bb,aa] +
                                 tmat[aa,bb]*(1-qab)/qab) / Dab
                }
            ans = rbind(ans, w[ii] * c(answer, deriv0))
        }
        dimnames(ans) = dimnames(eta)
        ans
    }), list( .refvalue=refvalue, .refgp=refgp) )),
    weight= eval(substitute(expression({
        wz = matrix(0, n, dimm(M))   # includes diagonal
        for(ii in 1:nrow(eta)) {
            alpha = .brat.alpha(eta2theta(eta[ii,-M], "loge"), .refvalue, .refgp)
            alpha0 = eta2theta(eta[ii,M], "loge") # M == ncol(eta)
            ymat = InverseBrat(y[ii,], NCo=M, diag=0)
            tmat = InverseBrat(ties[ii,], NCo=M, diag=0)
            for(aa in 1:(NCo)) {
                Daj = alpha[aa] + alpha[uindex] + alpha0
                pja = alpha[uindex] / Daj
                nja = ymat[aa,uindex] + ymat[uindex,aa] + tmat[uindex,aa]
                wz[ii,1:(NCo-1)] = wz[ii,1:(NCo-1)] +
                    alpha[uindex]^2 * nja * (1-pja)/(pja * Daj^2)
                if(aa < NCo)
                    for(bb in (aa+1):(NCo)) {
                        nab = ymat[aa,bb] + ymat[bb,aa] + tmat[bb,aa]
                        Dab = alpha[aa] + alpha[bb] + alpha0
                        qab = alpha0 / Dab
                        wz[ii,NCo] = wz[ii,NCo] + alpha0^2 * nab *
                                 (1-qab) / (qab * Dab^2)
                    }
            }
            if(NCo > 2) {
                ind5 = iam(1,1, M=NCo, both=TRUE, diag=FALSE)
                alphajunk = c(alpha, junk=NA)
                mat4 = cbind(uindex[ind5$row],uindex[ind5$col])
                wz[ii,(M+1):ncol(wz)] =
                  -(ymat[mat4] + ymat[mat4[,2:1]] + tmat[mat4]) *
                    alphajunk[uindex[ind5$col]] * alphajunk[uindex[ind5$row]] /
                    (alpha0 + alphajunk[uindex[ind5$row]] +
                     alphajunk[uindex[ind5$col]])^2
            }
            for(j in 1:length(uindex)) {
                jay = uindex[j]
                naj = ymat[,jay] + ymat[jay,] + tmat[,jay]
                Daj = alpha[jay] + alpha + alpha0
                wz[ii,iam(j, NCo, M=NCo, diag=TRUE)] = 
                    -alpha[jay] * alpha0 * sum(naj / Daj^2)
            }
        }
        wz = wz * w
        wz
    }), list( .refvalue=refvalue, .refgp=refgp ))))
}


.brat.alpha = function(vec, value, posn) {
    if(is.character(posn))
        if(posn!="last")
            stop("can only handle \"last\"") else return(c(vec, value))
    c(if(posn==1) NULL else vec[1:(posn-1)], value,
      if(posn==length(vec)+1) NULL else vec[posn:length(vec)])
}

.brat.indices = function(NCo, are.ties=FALSE) {
    if(!is.Numeric(NCo, allow=1, integ=TRUE) || NCo < 2)
        stop("bad input for NCo")
    m = diag(NCo)
    if(are.ties) {
        cbind(rindex=row(m)[col(m) < row(m)], cindex=col(m)[col(m) < row(m)])
    } else
        cbind(rindex=row(m)[col(m) != row(m)], cindex=col(m)[col(m) != row(m)])
}


Brat = function(mat, ties=0*mat, string=c(" > "," == ")) {
    allargs = list(mat)  # ,...
    callit = if(length(names(allargs))) names(allargs) else
        as.character(1:length(allargs))
    ans = ans.ties = NULL
    for(i in 1:length(allargs)) {
        m = allargs[[i]]
        if(!is.matrix(m) || dim(m)[1] != dim(m)[2]) 
            stop("m must be a square matrix")

        diag(ties) = 0
        if(!all(ties == t(ties)))
            stop("ties must be a symmetric matrix")
        are.ties = any(ties > 0)
        diag(ties) = NA

        diag(m) = 0   # Could have been NAs
        if(any(is.na(m)))
            stop("missing values not allowed (except on the diagonal)")
        diag(m) = NA

        dm = as.data.frame.table(m)
        dt = as.data.frame.table(ties)
        dm = dm[!is.na(dm$Freq),]
        dt = dt[!is.na(dt$Freq),]
        usethis1 = paste(dm[,1], string[1], dm[,2], sep="")
        usethis2 = paste(dm[,1], string[2], dm[,2], sep="")
        ans = rbind(ans, matrix(dm$Freq, nrow=1))
        ans.ties = rbind(ans.ties, matrix(dt$Freq, nrow=1))
    }
    dimnames(ans) = list(callit, usethis1)
    dimnames(ans.ties) = list(callit, usethis2)
    attr(ans, "ties") = ans.ties 
    attr(ans, "are.ties") = are.ties 
    ans
}


InverseBrat = function(yvec, NCo=
                      (1:900)[(1:900)*((1:900)-1)==ncol(rbind(yvec))],
                       multiplicity=if(is.matrix(yvec)) nrow(yvec) else 1,
                       diag=NA, string=c(" > "," == ")) {
    ans = array(diag, c(NCo, NCo, multiplicity))
    yvec.orig = yvec
    yvec = c(yvec)
    ptr = 1
    for(mul in 1:multiplicity)
    for(i1 in 1:(NCo))
        for(i2 in 1:(NCo))
            if(i1 != i2) {
                ans[i2,i1,mul] = yvec[ptr]
                ptr = ptr + 1
            }
    ans = if(multiplicity>1) ans else matrix(ans, NCo, NCo)

    if(is.array(yvec.orig) || is.matrix(yvec.orig)) {
        names.yvec = dimnames(yvec.orig)[[2]]
        i = strsplit(names.yvec, string[1])
        cal = NULL
        for(k in c(NCo, 1:(NCo-1))) cal = c(cal, (i[[k]])[1])
        if(multiplicity>1) {
            dimnames(ans) = list(cal, cal, dimnames(yvec.orig)[[1]])
        } else 
            dimnames(ans) = list(cal, cal)
    } 
    ans
}




tapplymat1 <- function(mat, function.arg=c("cumsum", "diff", "cumprod"))
{


    if(!missing(function.arg))
        function.arg <- as.character(substitute(function.arg))
    function.arg <- match.arg(function.arg, c("cumsum", "diff", "cumprod"))[1]

    type <- switch(function.arg,
        cumsum=1,
        diff=2,
        cumprod=3,
        stop("function.arg not matched"))

    if(!is.matrix(mat))
        mat <- as.matrix(mat)
    nr <- nrow(mat)
    nc <- ncol(mat)
    fred <- dotC(name="tapplymat1", mat=as.double(mat),
        as.integer(nr), as.integer(nc), as.integer(type))

    dim(fred$mat) <- c(nr, nc)
    dimnames(fred$mat) <- dimnames(mat)
    switch(function.arg,
        cumsum=fred$mat,
        diff=fred$mat[,-1,drop=FALSE],
        cumprod=fred$mat)
}




ordpoisson = function(cutpoints,
                      countdata=FALSE, NOS=NULL, Levels=NULL,
                      init.mu = NULL, parallel=FALSE, zero=NULL,
                      link="loge", earg = list()) {
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()
    fcutpoints = cutpoints[is.finite(cutpoints)]
    if(!is.Numeric(fcutpoints, integ=TRUE) || any(fcutpoints < 0))
        stop("\"cutpoints\" must have non-negative integer or Inf values only")
    if(is.finite(cutpoints[length(cutpoints)]))
        cutpoints = c(cutpoints, Inf)

    if(!is.logical(countdata) || length(countdata)!=1)
        stop("\"countdata\" must be a single logical")
    if(countdata) {
        if(!is.Numeric(NOS, integ=TRUE, posit=TRUE))
            stop("\"NOS\" must have integer values only")
        if(!is.Numeric(Levels, integ=TRUE, posit=TRUE)  || any(Levels < 2))
            stop("\"Levels\" must have integer values (>= 2) only")
        Levels = rep(Levels, length=NOS)
    }

    new("vglmff",
    blurb=c(paste("Ordinal Poisson model\n\n"), 
           "Link:     ", namesof("mu", link, earg= earg)),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints,
                              intercept.apply=TRUE)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero=zero ))),
    initialize=eval(substitute(expression({
        orig.y = cbind(y) # Convert y into a matrix if necessary
        if( .countdata ) {
            extra$NOS = M = NOS = .NOS
            extra$Levels = Levels = .Levels
            y.names = dimnames(y)[[2]]  # Hopefully the user inputted them
        } else {
            if(any(w != 1) || ncol(cbind(w)) != 1)
                stop("the 'weights' argument must be a vector of all ones")
            extra$NOS = M = NOS = if(is.Numeric( .NOS )) .NOS else
                ncol(orig.y)
            Levels = rep( if(is.Numeric( .Levels )) .Levels else 0, len=NOS)
            if(!is.Numeric( .Levels ))
                for(iii in 1:NOS) {
                    Levels[iii] = length(unique(sort(orig.y[,iii])))
                }
            extra$Levels = Levels
        }


        initmu = if(is.Numeric( .init.mu )) rep( .init.mu, len=NOS) else NULL
        cutpoints = rep( .cutpoints, len=sum(Levels))
        delete.zero.colns = FALSE 
        use.y = if( .countdata ) y else matrix(0, n, sum(Levels))
        use.etastart = matrix(0, n, M)
        cptr = 1
        for(iii in 1:NOS) {
            y = factor(orig.y[,iii], levels=(1:Levels[iii]))
            if( !( .countdata )) {
                eval(process.categorical.data.vgam)  # Creates mustart and y
                use.y[,cptr:(cptr+Levels[iii]-1)] = y
            }
            use.etastart[,iii] = if(is.Numeric(initmu)) initmu[iii] else
                median(cutpoints[cptr:(cptr+Levels[iii]-1-1)])
            cptr = cptr + Levels[iii]
        }
        mustart = NULL  # Overwrite it
        etastart = theta2eta(use.etastart, .link, earg= .earg)
        y = use.y  # n x sum(Levels)
        M = NOS
        for(iii in 1:NOS) {
            mu.names = paste("mu", iii, ".", sep="")
        }

        ncoly = extra$ncoly = sum(Levels)
        cp.vector = rep( .cutpoints, length=ncoly)
        extra$countdata = .countdata
        extra$cutpoints = cp.vector
        extra$n = n
        mynames = if(M > 1) paste("mu",1:M,sep="") else "mu"
        predictors.names = namesof(mynames, .link, short=TRUE, earg= .earg)
    }), list( .link=link, .countdata = countdata, .earg = earg,
              .cutpoints=cutpoints, .NOS=NOS, .Levels=Levels,
              .init.mu = init.mu
            ))),
    inverse=eval(substitute( function(eta, extra=NULL) {
        mu = eta2theta(eta, link= .link, earg= .earg) # Poisson means
        mu = cbind(mu)
        mu
    }, list( .link=link, .earg= earg, .countdata = countdata ))),
    last=eval(substitute(expression({
        if( .countdata ) {
            misc$link = .link
            misc$earg = list( .earg )
        } else {
            misc$link = rep( .link, length=M)
            names(misc$link) = mynames
            misc$earg = vector("list", M)
            names(misc$earg) = names(misc$link)
            for(ii in 1:M) misc$earg[[ii]] = .earg
        }
        misc$parameters = mynames
        misc$countdata = .countdata
        misc$true.mu <- FALSE    # $fitted is not a true mu
    }), list( .link=link, .countdata = countdata, .earg= earg ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
        probs = ordpoissonProbs(extra, mu)
        index0 <- y == 0
        probs[index0] = 1
        pindex0 <- probs == 0
        probs[pindex0] = 1
        sum(pindex0) * (-1.0e+10) + sum(w * y * log(probs))}},
    vfamily=c("ordpoisson", "vcategorical"),
    deriv=eval(substitute(expression({
        probs = ordpoissonProbs(extra, mu)
        probs.use = pmax(probs, .Machine$double.eps * 1.0e-0)

        cp.vector = extra$cutpoints
        NOS = extra$NOS
        Levels = extra$Levels
        answer.matrix = matrix(0, n, M)
        dl.dprob = y / probs.use
        dmu.deta = dtheta.deta(mu, .link, earg=.earg)
        dprob.dmu = ordpoissonProbs(extra, mu, deriv=1)
        cptr = 1
        for(iii in 1:NOS) {
            for(kkk in 1:Levels[iii]) {
                answer.matrix[,iii] = answer.matrix[,iii] + 
                    dl.dprob[,cptr] * dprob.dmu[,cptr]
                cptr = cptr + 1
            }
        }
        answer.matrix = w * answer.matrix * dmu.deta
        answer.matrix
    }), list( .link=link, .earg= earg, .countdata=countdata ))),
    weight= eval(substitute(expression({
        d2l.dmu2 = matrix(0, n, M)  # Diagonal matrix
        cptr = 1
        for(iii in 1:NOS) {
            for(kkk in 1:Levels[iii]) {
                d2l.dmu2[,iii] = d2l.dmu2[,iii] + 
                    dprob.dmu[,cptr]^2 / probs.use[,cptr]
                cptr = cptr + 1
            }
        }
        wz = w * d2l.dmu2 * dmu.deta^2
        wz
    }), list( .earg=earg, .link=link, .countdata=countdata ))))
}

ordpoissonProbs = function(extra, mu, deriv=0) {
    cp.vector = extra$cutpoints
    NOS = extra$NOS
    if(deriv == 1) {
        dprob.dmu = matrix(0, extra$n, extra$ncoly)
    } else {
        probs = matrix(0, extra$n, extra$ncoly)
    }
    mu = cbind(mu)
    cptr = 1
    for(iii in 1:NOS) {
        if(deriv == 1) {
            dprob.dmu[,cptr] = -dpois(x=cp.vector[cptr], lamb=mu[,iii])
        } else {
            probs[,cptr] = ppois(q=cp.vector[cptr], lambda=mu[,iii])
        }
        cptr = cptr + 1
        while(is.finite(cp.vector[cptr])) {
            if(deriv == 1) {
                dprob.dmu[,cptr] = dpois(x=cp.vector[cptr-1], lamb=mu[,iii]) -
                        dpois(x=cp.vector[cptr], lambda=mu[,iii])
            } else {
                probs[,cptr] = ppois(q=cp.vector[cptr], lambda=mu[,iii]) -
                        ppois(q=cp.vector[cptr-1], lambda=mu[,iii])
            }
            cptr = cptr + 1
        }
        if(deriv == 1) {
            dprob.dmu[,cptr] = dpois(x=cp.vector[cptr-1], lamb=mu[,iii]) -
                    dpois(x=cp.vector[cptr], lambda=mu[,iii])
        } else {
            probs[,cptr] = ppois(q=cp.vector[cptr], lamb=mu[,iii]) -
                    ppois(q=cp.vector[cptr-1], lambda=mu[,iii])
        }
        cptr = cptr + 1
    }
    if(deriv == 1) dprob.dmu else probs
}






scumulative = function(link="logit", earg = list(),
                       lscale="loge", escale = list(),
                       parallel=FALSE, sparallel=TRUE, reverse=FALSE,
                       iscale = 1)
{
    stop("sorry, not working yet")
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(!is.list(escale)) escale = list()
    if(!is.Numeric(iscale, posit=TRUE))
        stop("bad input for argument \"iscale\"")


    new("vglmff",
    blurb=c(paste("Scaled cumulative", link, "model\n\n"),
           "Links:   ",
           namesof(if(reverse) "P[Y>=j+1]" else "P[Y<=j]", link, earg=earg),
           ", ",
           namesof("scale_j", lscale, escale)),
    constraints=eval(substitute(expression({
        J = M / 2
        constraints = cm.vgam(matrix(1,J,1), x, .parallel, constraints,
                              intercept.apply = FALSE)
        constraints[["(Intercept)"]] = rbind(constraints[["(Intercept)"]],
            matrix(0, J, ncol(constraints[["(Intercept)"]])))

        cm2 = cm.vgam(matrix(1,J,1), x, .sparallel, constraints=NULL,
                      intercept.apply = FALSE)

        for(ii in 2:length(constraints))
            constraints[[ii]] =
                cbind(rbind(constraints[[ii]],
                            matrix(0, J, ncol(constraints[[ii]]))),
                      rbind(matrix(0, J, ncol(cm2[[ii]])), cm2[[ii]]))

        for(ii in 1:length(constraints))
            constraints[[ii]] =
                (constraints[[ii]])[interleave.VGAM(M, M=2),,drop=FALSE]
    }), list( .parallel=parallel, .sparallel=sparallel ))),
    deviance=eval(substitute(
        function(mu, y, w, residuals=FALSE, eta, extra=NULL) {
        answer =
            Deviance.categorical.data.vgam(mu=mu, y=y, w=w, residuals=residuals,
                                           eta=eta, extra=extra)
        answer
    }, list( .earg=earg, .link=link ) )),
    initialize=eval(substitute(expression({
        if(intercept.only)
            stop("use cumulative() for intercept-only models")
        delete.zero.colns = TRUE # Cannot have FALSE since then prob(Y=jay)=0
        eval(process.categorical.data.vgam)
        M = 2*(ncol(y)-1)
        J = M / 2
        extra$J = J
        mynames = if( .reverse ) paste("P[Y>=",2:(1+J),"]", sep="") else
            paste("P[Y<=",1:J,"]", sep="")
        predictors.names = c(
            namesof(mynames, .link, short=TRUE, earg= .earg),
            namesof(paste("scale_", 1:J, sep=""),
                    .lscale, short=TRUE, earg= .escale))
        y.names = paste("mu", 1:(J+1), sep="")

        if(length(dimnames(y)))
            extra$dimnamesy2 = dimnames(y)[[2]]

        predictors.names = predictors.names[interleave.VGAM(M, M=2)]

    }), list( .link=link, .lscale=lscale, .reverse=reverse,
              .earg= earg, .escale=escale ))),
    inverse=eval(substitute( function(eta, extra=NULL) {
        J = extra$J
        M = 2*J
        etamat1 = eta[,2*(1:J)-1,drop=FALSE]
        etamat2 = eta[,2*(1:J),  drop=FALSE]
        scalemat = eta2theta(etamat2, .lscale, earg= .escale)
        fv.matrix =
        if( .reverse ) {
            ccump = cbind(1, eta2theta(etamat1/scalemat, .link, earg=.earg))
            cbind(-tapplymat1(ccump, "diff"), ccump[,ncol(ccump)])
        } else {
            cump = cbind(eta2theta(etamat1/scalemat, .link, earg= .earg), 1)
            cbind(cump[,1], tapplymat1(cump, "diff"))
        }
        if(length(extra$dimnamesy2))
            dimnames(fv.matrix) = list(dimnames(eta)[[1]], extra$dimnamesy2)
        fv.matrix
    }, list( .link=link, .lscale=lscale, .reverse=reverse,
             .earg= earg, .escale=escale ))),
    last=eval(substitute(expression({
        J = extra$J
        misc$link = c(rep( .link, length=J),
                      rep( .lscale, length=J))[interleave.VGAM(M, M=2)]
        names(misc$link) = predictors.names  # zz mynames
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:J) misc$earg[[2*ii-1]] = .earg
        for(ii in 1:J) misc$earg[[2*ii  ]] = .escale
        misc$parameters = mynames
        misc$reverse = .reverse
        misc$parallel = .parallel
        misc$sparallel = .sparallel
    }), list( .link=link, .lscale=lscale,
              .reverse=reverse, .parallel=parallel, .sparallel=sparallel,
              .earg=earg, .escale=escale ))),
    link=eval(substitute( function(mu, extra=NULL) {
        cump = tapplymat1(as.matrix(mu), "cumsum")
        J = ncol(as.matrix(mu)) - 1
        M = 2 * J
        answer =  cbind(
            theta2eta(if( .reverse ) 1-cump[,1:J] else cump[,1:J], .link,
                      earg= .earg),
            matrix(theta2eta( .iscale, .lscale, earg = .escale),
                   nrow(as.matrix(mu)), J, byrow=TRUE))
        answer = answer[,interleave.VGAM(M, M=2)]
        answer
    }, list( .link=link, .lscale=lscale, .reverse=reverse,
             .iscale=iscale, .earg=earg, .escale=escale ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
       sum(w * y * log(mu)), 
    vfamily=c("scumulative", "vcategorical"),
    deriv=eval(substitute(expression({
        ooz = iter %% 2
 print("ooz")
 print( ooz )

        J = extra$J
        mu.use = pmax(mu, .Machine$double.eps * 1.0e-0)

        etamat1 = eta[,2*(1:J)-1,drop=FALSE]
        etamat2 = eta[,2*(1:J),  drop=FALSE]
        scalemat = eta2theta(etamat2, .lscale, earg= .escale)

        cump = eta2theta(etamat1 / scalemat, .link, earg= .earg)
        dcump.deta = dtheta.deta(cump, .link, earg= .earg)
        dscale.deta = dtheta.deta(scalemat, .lscale, earg= .escale)
        dl.dcump = (if( .reverse) -w  else w) * 
                (y[,1:J]/mu.use[,1:J] - y[,-1]/mu.use[,-1])
        dcump.dscale = -dcump.deta * etamat1 / scalemat^2
        ans = cbind(dl.dcump * dcump.deta / scalemat,
                    dl.dcump * dcump.dscale * dscale.deta)
        ans = ans[,interleave.VGAM(M, M=2)]
        if(ooz) ans[,c(TRUE,FALSE)] = 0 else ans[,c(FALSE,TRUE)] = 0
        ans
    }), list( .link=link, .lscale=lscale, .reverse=reverse,
              .earg= earg, .escale=escale ))),
    weight= eval(substitute(expression({

        wz = matrix(0, n, 2*(2*M-3))

        wz[,2*(1:J)-1] = if(ooz) w * (dcump.deta / scalemat)^2 *
                         (1/mu.use[,1:J] + 1/mu.use[,-1]) else 1
        wz[,2*(1:J)] = if(ooz) 1 else w * (dcump.dscale * dscale.deta)^2 *
                       (1/mu.use[,1:J] + 1/mu.use[,-1])
        wz0 = w * (dcump.deta / scalemat) * 
                  (dcump.dscale * dscale.deta) *
                  (1/mu.use[,1:J] + 1/mu.use[,-1])
        wz0 = as.matrix(wz0)
        for(ii in 1:J)
            wz[,iam(2*ii-1,2*ii,M=M)] = if(ooz) wz0[,ii] else 0

        if(J > 1) {
            wz0 = -w * (dcump.deta[,-J] / scalemat[,-J]) *
                       (dcump.deta[,-1]  / scalemat[,-1]) / mu.use[,2:J]
            wz0 = as.matrix(wz0) # Just in case J=2
            for(ii in 1:(J-1))
                wz[,iam(2*ii-1,2*ii+1,M=M)] = if(ooz) wz0[,ii] else 0
            wz0 = -w * (dcump.dscale[,-1] * dscale.deta[,-1]) *
                       (dcump.dscale[,-J] * dscale.deta[,-J]) / mu.use[,2:J]
            wz0 = as.matrix(wz0)
            for(ii in 1:(J-1))
                wz[,iam(2*ii,2*ii+2,M=M)] = if(ooz) wz0[,ii] else 0



            wz0 = -w * (dcump.deta[,-J] / scalemat[,-J]) *
                       (dcump.dscale[,-1] * dscale.deta[,-1]) / mu.use[,2:J]
            wz0 = as.matrix(wz0)
            for(ii in 1:(J-1))
                wz[,iam(2*ii-1,2*ii+2,M=M)] = if(ooz) wz0[,ii] else 0
            wz0 = -w * (dcump.deta[,-1] / scalemat[,-1]) *
                       (dcump.dscale[,-J] * dscale.deta[,-J]) / mu.use[,2:J]
            wz0 = as.matrix(wz0)
            for(ii in 1:(J-1))
                wz[,iam(2*ii,2*ii+1,M=M)] = if(ooz) wz0[,ii] else 0
        }
        wz
    }), list( .link=link, .lscale=lscale, .earg=earg, .escale=escale ))))
}






# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.






rcqo <- function(n, p, S,
                 Rank = 1,
                 family = c("poisson", "negbinomial", "binomial-poisson", 
                            "Binomial-negbinomial", "ordinal-poisson",
                            "Ordinal-negbinomial","gamma2"),
                 EqualMaxima = FALSE,
                 EqualTolerances = TRUE,
                 ESOptima = FALSE,
                 loabundance = if(EqualMaxima) hiabundance else 10,
                 hiabundance = 100,
                 sdlv = c(1.5/2^(0:3))[1:Rank],
                 sdOptima = ifelse(ESOptima, 1.5/Rank, 1) *
                            ifelse(scalelv, sdlv, 1),
                 sdTolerances = 0.25,
                 Kvector = 1,
                 Shape = 1,
                 sqrt = FALSE,
                 Log = FALSE,
                 rhox = 0.5,
                 breaks = 4, # ignored unless family="ordinal"
                 seed = NULL,
                 Crow1positive=TRUE,
                 xmat = NULL, # Can be input
                 scalelv = TRUE 
                 ) {
    family = match.arg(family, c("poisson","negbinomial", "binomial-poisson",
         "Binomial-negbinomial", "ordinal-poisson",
         "Ordinal-negbinomial","gamma2"))[1]
    if(!is.Numeric(n, integer=TRUE, posit=TRUE, allow=1))
        stop("bad input for argument 'n'")
    if(!is.Numeric(p, integer=TRUE, posit=TRUE, allow=1) || p < 1 + Rank)
        stop("bad input for argument 'p'")
    if(!is.Numeric(S, integer=TRUE, posit=TRUE, allow=1))
        stop("bad input for argument 'S'")
    if(!is.Numeric(Rank, integer=TRUE, posit=TRUE, allow=1) || Rank > 4)
        stop("bad input for argument 'Rank'")
    if(!is.Numeric(Kvector, posit=TRUE))
        stop("bad input for argument 'Kvector'")
    if(!is.Numeric(rhox) || abs(rhox) >= 1)
        stop("bad input for argument 'rhox'")
    if(length(seed) && !is.Numeric(seed, integer=TRUE, posit=TRUE))
        stop("bad input for argument 'seed'")
    if(!is.logical(EqualTolerances) || length(EqualTolerances)>1)
        stop("bad input for argument 'EqualTolerances)'")
    if(!is.logical(sqrt) || length(sqrt)>1)
        stop("bad input for argument 'sqrt)'")
    if(family != "negbinomial" && sqrt)
        warning("argument 'sqrt' is used only with family='negbinomial'")
    if(!EqualTolerances && !is.Numeric(sdTolerances, posit=TRUE))
        stop("bad input for argument 'sdTolerances'")
    if(!is.Numeric(loabundance, posit=TRUE))
        stop("bad input for argument 'loabundance'")
    if(!is.Numeric(sdlv, posit=TRUE))
        stop("bad input for argument 'sdlv'")
    if(!is.Numeric(sdOptima, posit=TRUE))
        stop("bad input for argument 'sdOptima'")
    if(EqualMaxima && loabundance != hiabundance)
        stop(paste("arguments 'loabundance)' and 'hiabundance)' must",
                   "be equal when EqualTolerances=TRUE"))
    if(any(loabundance > hiabundance))
        stop("loabundance > hiabundance is not allowed")
    if(!is.logical(Crow1positive)) {
        stop("bad input for argument 'Crow1positive)'")
    } else {
        Crow1positive = rep(Crow1positive, len=Rank)
    }
    Shape = rep(Shape, len=S)
    sdlv = rep(sdlv, len=Rank)
    sdOptima = rep(sdOptima, len=Rank)
    sdTolerances = rep(sdTolerances, len=Rank)
    AA = sdOptima / 3^0.5
    if(Rank > 1 && any(diff(sdlv) > 0))
     stop("argument 'sdlv)' must be a vector with decreasing values")

    if(FALSE)
    change.seed.expression = expression({
        if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            runif(1)                       # initialize the RNG if necessary
        }
        if(is.null(seed)) {
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        } else {
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }
    })
    change.seed.expression = expression({
        if(length(seed)) set.seed(seed)
    })
    eval(change.seed.expression)

    V = matrix(rhox, p-1, p-1)
    diag(V) = 1
    L = chol(V)
    if(length(xmat)) {
        xnames = colnames(xmat)
    } else {
        eval(change.seed.expression)
        xmat = matrix(rnorm(n*(p-1)), n, p-1) %*% L
        xmat = scale(xmat, center=TRUE)
        xnames = paste("x", 2:p, sep="")
        dimnames(xmat) = list(as.character(1:n), xnames)
    }
    eval(change.seed.expression)
    ccoefs = matrix(rnorm((p-1)*Rank), p-1, Rank)
    lvmat = cbind(xmat %*% ccoefs)
    if(Rank > 1) {
        Rmat = chol(var(lvmat))
        iRmat = solve(Rmat)
        lvmat = lvmat %*% iRmat  # var(lvmat) == diag(Rank)
        ccoefs = ccoefs %*% iRmat
    }
    for(r in 1:Rank)
        if(( Crow1positive[r] && ccoefs[1,r] < 0) ||
           (!Crow1positive[r] && ccoefs[1,r] > 0)) {
                ccoefs[,r] = -ccoefs[,r]
                lvmat[,r] = -lvmat[,r]
        }

    if(scalelv) {
        for(r in 1:Rank) {
            sdlvr = sd(lvmat[,r])
            lvmat[,r] = lvmat[,r] * sdlv[r] / sdlvr
            ccoefs[,r]  = ccoefs[,r] * sdlv[r] / sdlvr
        }
    } else {
        sdlvr = NULL
        for(r in 1:Rank) {
            sdlvr = c(sdlvr, sd(lvmat[,r]))
        }
    }
    if(ESOptima) {
        if(!is.Numeric(S^(1/Rank), integ=TRUE) || S^(1/Rank) < 2)
            stop("S^(1/Rank) must be an integer greater or equal to 2")
        if(Rank == 1) {
            optima = matrix(as.numeric(NA), S, Rank)
            for(r in 1:Rank) {
                optima[,r] = seq(-AA, AA, len=S^(1/Rank))
            }
        } else if(Rank == 2) {
            optima = expand.grid(lv1=seq(-AA[1], AA[1], len=S^(1/Rank)),
                                 lv2=seq(-AA[2], AA[2], len=S^(1/Rank)))
        } else if(Rank == 3) {
            optima = expand.grid(lv1=seq(-AA[1], AA[1], len=S^(1/Rank)),
                                 lv2=seq(-AA[2], AA[2], len=S^(1/Rank)),
                                 lv3=seq(-AA[3], AA[3], len=S^(1/Rank)))
        } else {
            optima = expand.grid(lv1=seq(-AA[1], AA[1], len=S^(1/Rank)),
                                 lv2=seq(-AA[2], AA[2], len=S^(1/Rank)),
                                 lv3=seq(-AA[3], AA[3], len=S^(1/Rank)),
                                 lv4=seq(-AA[4], AA[4], len=S^(1/Rank)))
        }
        if(Rank > 1)
            optima = matrix(unlist(optima), S, Rank)  # Make sure it is a matrix
    } else {
        optima = matrix(1, S, Rank)
        eval(change.seed.expression)
        for(r in 1:Rank) {
            optima[,r] = rnorm(n=S, sd=sdOptima[r])
        }
    }
    for(r in 1:Rank)
        optima[,r] = optima[,r] * sdOptima[r] / sd(optima[,r])

    ynames = paste("y", 1:S, sep="")
    Kvector = rep(Kvector, len=S)
    names(Kvector) = ynames
    lvnames = if(Rank==1) "lv" else paste("lv", 1:Rank, sep="")
    Tols = if(EqualTolerances) matrix(1, S, Rank) else {
               eval(change.seed.expression)
               temp = matrix(1, S, Rank)
               if(S > 1)
               for(r in 1:Rank) {
                   temp[-1,r] = rnorm(S-1, mean=1, sd=sdTolerances[r])
                   if(any(temp[,r] <= 0)) stop("negative tolerances!")
                   temp[,r] = temp[,r]^2 # Tolerance matrix  = var-cov matrix)
               }
               temp
           }

    dimnames(Tols) = list(ynames, lvnames)
    dimnames(ccoefs) = list(xnames, lvnames)
    dimnames(optima) = list(ynames, lvnames)
    loeta = log(loabundance)  # May be a vector
    hieta = log(hiabundance)
    eval(change.seed.expression)
    logmaxima = runif(S, min=loeta, max=hieta)  # loeta and hieta may be vector
    names(logmaxima) = ynames
    etamat = matrix(logmaxima,n,S,byrow=TRUE) # eta=log(mu) only; intercept term
    for(jay in 1:S) {
        optmat = matrix(optima[jay,], nrow=n, ncol=Rank, byrow=TRUE)
        tolmat = matrix(Tols[jay,], nrow=n, ncol=Rank, byrow=TRUE)
        temp = cbind((lvmat - optmat) / tolmat)
        for(r in 1:Rank)
            etamat[,jay]=etamat[,jay]-0.5*(lvmat[,r] - optmat[jay,r])*temp[,r]
    }

    rootdist = switch(family,
        "poisson"=1, "binomial-poisson"=1, "ordinal-poisson"=1,
        "negbinomial"=2, "Binomial-negbinomial"=2, "Ordinal-negbinomial"=2,
        "gamma2"=3)
    eval(change.seed.expression)
    if(rootdist == 1) {
        ymat = matrix(rpois(n*S, lam=exp(etamat)), n, S)
    } else if(rootdist == 2) {
        mKvector = matrix(Kvector, n, S, byrow=TRUE)
        ymat = matrix(rnbinom(n=n*S, mu=exp(etamat), size=mKvector),n,S)
        if(sqrt) ymat = ymat^0.5
    } else if(rootdist == 3) {
        Shape = matrix(Shape, n, S, byrow=TRUE)
        ymat = matrix(rgamma(n*S, shape=Shape, scale=exp(etamat)/Shape),n,S)
        if(Log) ymat = log(ymat)
    } else stop("argument 'rootdist' unmatched")

    tmp1 = NULL
    if(any(family == c("ordinal-poisson","Ordinal-negbinomial"))) {
        tmp1 = cut(c(ymat), breaks=breaks, labels=NULL) #To get attributes(tmp1)
        ymat = cut(c(ymat), breaks=breaks, labels=FALSE)
        dim(ymat) = c(n,S)
    }
    if(any(family == c("binomial-poisson","Binomial-negbinomial")))
        ymat = 0 + (ymat > 0)

    myform = as.formula(paste(paste("cbind(",
             paste(paste("y",1:S,sep=""), collapse=","),
             ") ~ ", sep=""),
             paste(paste("x",2:p,sep=""), collapse="+"), sep=""))

    dimnames(ymat) = list(as.character(1:n), ynames)
    ans = data.frame(xmat, ymat)
    attr(ans, "ccoefficients") = ccoefs
    attr(ans, "Crow1positive") = Crow1positive
    attr(ans, "family") = family
    attr(ans, "formula") = myform # Useful for running cqo() on the data
    attr(ans, "Rank") = Rank
    attr(ans, "family") = family
    attr(ans, "Kvector") = Kvector
    attr(ans, "logmaxima") = logmaxima
    attr(ans, "loabundance") = loabundance
    attr(ans, "hiabundance") = hiabundance
    attr(ans, "optima") = optima
    attr(ans, "Log") = Log
    attr(ans, "lv") = lvmat
    attr(ans, "eta") = etamat
    attr(ans, "EqualTolerances") = EqualTolerances
    attr(ans, "EqualMaxima") = EqualMaxima || all(loabundance == hiabundance)
    attr(ans, "ESOptima") = ESOptima
    attr(ans, "seed") = seed # RNGstate
    attr(ans, "sdTolerances") = sdTolerances
    attr(ans, "sdlv") =  if(scalelv) sdlv else sdlvr
    attr(ans, "sdOptima") = sdOptima
    attr(ans, "Shape") = Shape
    attr(ans, "sqrt") = sqrt
    attr(ans, "tolerances") = Tols^0.5  # Like a standard deviation
        attr(ans, "breaks") = if(length(tmp1)) attributes(tmp1) else breaks
    ans
}




if(FALSE)
dcqo <- function(x, p, S,
                 family = c("poisson", "binomial", "negbinomial", "ordinal"),
                 Rank = 1,
                 EqualTolerances = TRUE,
                 EqualMaxima = FALSE,
                 EquallySpacedOptima = FALSE,
                 loabundance = if(EqualMaxima) 100 else 10,
                 hiabundance = 100,
                 sdTolerances = 1,
                 sdOptima = 1,
                 nlevels = 4, # ignored unless family="ordinal"
                 seed = NULL
                 ) {
 warning("12/6/06; needs a lot of work based on rcqo()")


    if(mode(family) != "character" && mode(family) != "name")
        family = as.character(substitute(family))
    family = match.arg(family, c("poisson", "binomial",
                                 "negbinomial", "ordinal"))[1]
    if(!is.Numeric(p, integer=TRUE, posit=TRUE, allow=1) || p < 2)
        stop("bad input for argument 'p'")
    if(!is.Numeric(S, integer=TRUE, posit=TRUE, allow=1))
        stop("bad input for argument 'S'")
    if(!is.Numeric(Rank, integer=TRUE, posit=TRUE, allow=1))
        stop("bad input for argument 'Rank'")
    if(length(seed) && !is.Numeric(seed, integer=TRUE, posit=TRUE))
        stop("bad input for argument 'seed'")
    if(!is.logical(EqualTolerances) || length(EqualTolerances)>1)
        stop("bad input for argument 'EqualTolerances)'")
    if(EqualMaxima && loabundance != hiabundance)
        stop(paste("'loabundance)' and 'hiabundance)' must",
                   "be equal when EqualTolerances=TRUE"))
    if(length(seed)) set.seed(seed)

    xmat = matrix(rnorm(n*(p-1)), n, p-1, dimnames=list(as.character(1:n),
                  paste("x", 2:p, sep="")))
    ccoefs = matrix(rnorm((p-1)*Rank), p-1, Rank)
    lvmat = xmat %*% ccoefs
    optima = matrix(rnorm(Rank*S, sd=sdOptima), S, Rank)
    Tols = if(EqualTolerances) matrix(1, S, Rank) else
           matrix(rnorm(Rank*S, mean=1, sd=1), S, Rank)
    loeta = log(loabundance)
    hieta = log(hiabundance)
    logmaxima = runif(S, min=loeta, max=hieta)

    etamat = matrix(logmaxima,n,S,byrow=TRUE) # eta=log(mu) only; intercept term
    for(jay in 1:S) {
        optmat = matrix(optima[jay,], n, Rank, byrow=TRUE)
        tolmat = matrix(Tols[jay,], n, Rank, byrow=TRUE)
        temp = cbind((lvmat - optmat) * tolmat)
        for(r in 1:Rank)
            etamat[,jay] = etamat[,jay] - 0.5 * temp[,r] *
                           (lvmat[,r] - optmat[jay,r])
    }

    ymat = if(family == "negbinomial") {



    } else {
           matrix(rpois(n*S, lam=exp(etamat)), n, S)
    }
    if(family == "binomial")
        ymat = 0 + (ymat > 0)

    dimnames(ymat) = list(as.character(1:n), paste("y", 1:S, sep=""))
    ans = data.frame(xmat, ymat)
    attr(ans, "ccoefficients") = ccoefs
    attr(ans, "family") = family
    ans
}





getInitVals = function(gvals, llfun, ...) {
    LLFUN = match.fun(llfun)
    ff = function(myx, ...) LLFUN(myx, ...)
    objFun = gvals
    for(ii in 1:length(gvals))
        objFun[ii] = ff(myx=gvals[ii], ...) 
    try.this = gvals[objFun == max(objFun)]  # Usually scalar, maybe vector
    try.this
}









campp = function(q, size, prob, mu) {
    if (!missing(mu)) {
        if(!missing(prob))
            stop("'prob' and 'mu' both specified")
        prob <- size/(size + mu)
    }
    K = (1/3) * ((9*q+8)/(q+1) - ((9*size-1)/size) * (mu/(q+1))^(1/3)) /
        sqrt( (1/size) * (mu/(q+1))^(2/3) + 1 / (q+1)) # Note the +, not -
    pnorm(K)
}











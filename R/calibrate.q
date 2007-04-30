# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.








calibrate.qrrvglm.control = function(object,
        trace=FALSE,  # passed into optim()
        Method.optim="BFGS",   # passed into optim(method=Method)
        gridSize = if(Rank==1) 9 else 5,
        varlvI = FALSE, ...) {
    Rank = object@control$Rank
    EqualTolerances = object@control$EqualTolerances
    if(!is.Numeric(gridSize, positive=TRUE, integer=TRUE, allow=1))
        stop("bad input for \"gridSize\"")
    if(gridSize < 2)
        stop("gridSize must be >= 2")
    list(# maxit=Maxit.optim,   # Note the name change
         trace=as.numeric(trace)[1],
         Method.optim=Method.optim,
         gridSize=gridSize,
         varlvI = as.logical(varlvI)[1])
} 

if(!isGeneric("calibrate"))
    setGeneric("calibrate", function(object, ...) standardGeneric("calibrate"))


calibrate.qrrvglm = function(object, 
                             newdata=NULL,
                        type=c("lv","predictors","response","vcov","all3or4"),
                             initial.vals=NULL, ...) {

    Quadratic = if(is.logical(object@control$Quadratic))
                object@control$Quadratic else FALSE  # T if CQO, F if CAO

    if(!length(newdata)) {
        if(!length(object@y)) stop("no newdata") else
        newdata = data.frame(object@y)
    }

    if(mode(type) != "character" && mode(type) != "name")
        type <- as.character(substitute(type))
    type <- match.arg(type, c("lv","predictors","response","vcov","all3or4"))[1]

    if(!Quadratic && type=="vcov")
        stop("cannot have type=\"vcov\" when object is a \"cao\" object")

    if(is.vector(newdata))
        newdata = rbind(newdata)
    if(!is.matrix(newdata))
        newdata = as.matrix(newdata)
    newdata = newdata[,object@misc$ynames,drop=FALSE]

    obfunct = slot(object@family, object@misc$criterion) # Objective function
    minimize.obfunct = if(Quadratic) object@control$min.criterion else
        TRUE  # Logical; TRUE for CAO objects because deviance is minimized
    if(!is.logical(minimize.obfunct)) 
        stop("object@control$min.criterion is not a logical")
    optim.control = calibrate.qrrvglm.control(object=object, ...) # For cao too

    if((Rank <- object@control$Rank) > 2)
        stop("currently can only handle Rank=1 and 2")
    Coefobject = if(Quadratic) {
        Coef(object, varlvI=optim.control$varlvI)
    } else {
        Coef(object)
    }
    if(!length(initial.vals)) {
        L = apply(Coefobject@lv, 2, min)
        U = apply(Coefobject@lv, 2, max)
        initial.vals = if(Rank==1)
            cbind(seq(L, U, length=optim.control$gridSize)) else
            expand.grid(seq(L[1], U[1], length=optim.control$gridSize),
                        seq(L[2], U[2], length=optim.control$gridSize))
    }
    ok = length(object@control$colx1.index)==1 &&
         names(object@control$colx1.index) == "(Intercept)"
    if(!ok) stop("The x1 vector must be an intercept only")

    nn = nrow(newdata)
    BestOFpar = NULL   # It may be more efficient not to append 
    BestOFvalues = NULL   # Best OF objective function values
    for(i1 in 1:nn) {
        if(optim.control$trace)
            cat("\nOptimizing for observation", i1, "-----------------\n")
        OFvalues = OFpar = NULL   # OF means objective function
        for(ii in 1:nrow(initial.vals)) {
            if(optim.control$trace) {
                cat("Starting from grid-point", ii, ":")
                if(exists("flush.console"))
                    flush.console()
            }
            ans = if(is.R()) {
                if(Quadratic)
                optim(par=initial.vals[ii,],
                      fn=.my.calib.objfunction.qrrvglm,
                      method=optim.control$Method.optim,  # "BFGS", or "CG" or ...
                      control=c(fnscale=ifelse(minimize.obfunct,1,-1),
                                optim.control),
                      y=newdata[i1,],
                      extra=object@extra,
                      objfun=obfunct,
                      Coefs=Coefobject,
                      misc.list = object@misc,
                      everything = FALSE,
                      mu.function = slot(object@family, "inverse")) else
                optim(par=initial.vals[ii,],
                      fn=.my.calib.objfunction.cao,
                      method=optim.control$Method.optim,  # "BFGS", or "CG" or ...
                      control=c(fnscale=ifelse(minimize.obfunct,1,-1),
                                optim.control),
                      y=newdata[i1,],
                      extra=object@extra,
                      objfun=obfunct,
                      object=object,
                      Coefs=Coefobject,
                      misc.list = object@misc,
                      everything = FALSE,
                      mu.function = slot(object@family, "inverse"))
            } else 
                stop("not implemented in S-PLUS yet")

            if(optim.control$trace) {
                if(ans$convergence == 0)
                    cat("Successful convergence\n") else 
                    cat("Unsuccessful convergence\n")
                if(exists("flush.console"))
                    flush.console()
            }
            if(ans$convergence == 0) {
                OFvalues = c(OFvalues, ans$value)
                OFpar = rbind(OFpar, ans$par)
            }
        }
        if(length(OFpar)) {
            index = if(minimize.obfunct)
                    (1:nrow(OFpar))[OFvalues==min(OFvalues)] else
                    (1:nrow(OFpar))[OFvalues==max(OFvalues)]
            if(length(index) > 1) {
                warning(paste("multiple solutions found for observation ", i1,
                              ". Choosing one randomly.", sep=""))
                index = sample(index, size=1)
            } else if(length(index) == 0)
                stop("length(index) is zero")
            BestOFpar = rbind(BestOFpar, OFpar[index,])
            BestOFvalues = c(BestOFvalues, OFvalues[index])
        } else {
            BestOFpar = rbind(BestOFpar, rep(as.numeric(NA), len=Rank))
            BestOFvalues = c(BestOFvalues, NA)
        }
    }

    pretty = function(BestOFpar, newdata, Rank) {
        if(Rank==1) {
            BestOFpar = c(BestOFpar) 
            names(BestOFpar) = dimnames(newdata)[[1]]
        } else
            dimnames(BestOFpar) = list(dimnames(newdata)[[1]],
                if(Rank==1) "lv" else paste("lv", 1:Rank, sep=""))
        BestOFpar
    }

    if(type=="lv") {
        BestOFpar = pretty(BestOFpar, newdata, Rank)
        attr(BestOFpar,"objectiveFunction")=pretty(BestOFvalues,newdata,Rank=1)
        BestOFpar
    } else {
        etaValues = muValues = NULL   #
        if(Quadratic)
            vcValues = array(0, c(Rank,Rank,nn))
        for(i1 in 1:nn) {
            ans = if(Quadratic) .my.calib.objfunction.qrrvglm(BestOFpar[i1, ],
                          y=newdata[i1,],
                          extra=object@extra,
                          objfun=obfunct,
                          Coefs=Coefobject,
                          misc.list = object@misc,
                          everything = TRUE,
                          mu.function = slot(object@family, "inverse")) else
                  .my.calib.objfunction.cao(BestOFpar[i1, ],
                          y=newdata[i1,],
                          extra=object@extra,
                          objfun=obfunct,
                          object=object,
                          Coefs=Coefobject,
                          misc.list = object@misc,
                          everything = TRUE,
                          mu.function = slot(object@family, "inverse"))
            muValues = rbind(muValues, matrix(ans$mu, nrow=1))
            etaValues = rbind(etaValues, matrix(ans$eta, nrow=1))
            if(Quadratic)
                vcValues[,,i1] = ans$vcmat  # Can be NULL for "cao" objects
        }
        if(type=="response") {
             dimnames(muValues) = dimnames(newdata)
             muValues
        } else if(type=="predictors") {
             dimnames(etaValues) = list(dimnames(newdata)[[1]],
                                        dimnames(object@predictors)[[2]])
             etaValues
        } else if(type=="vcov") {
             if(Quadratic)
             dimnames(vcValues) = list(as.character(1:Rank), 
                                       as.character(1:Rank),
                                       dimnames(newdata)[[1]])
             vcValues
        } else if(type=="all3or4") {
             if(Quadratic)
             dimnames(vcValues) = list(as.character(1:Rank), 
                                       as.character(1:Rank),
                                       dimnames(newdata)[[1]])
             dimnames(muValues) = dimnames(newdata)
             dimnames(etaValues) = list(dimnames(newdata)[[1]],
                                        dimnames(object@predictors)[[2]])
             BestOFpar = pretty(BestOFpar, newdata, Rank)
             attr(BestOFpar,"objectiveFunction") =
                  pretty(BestOFvalues,newdata,Rank=1)
             list(lv=BestOFpar,
                  predictors=etaValues,
                  response=muValues,
                  vcov=if(Quadratic) vcValues else NULL)
        } else stop("type not matched")
    }
}
       
.my.calib.objfunction.qrrvglm = function(bnu, y, extra=NULL,
                        objfun, Coefs,
                        misc.list,
                        everything=TRUE,
                        mu.function) {

    bnumat = cbind(bnu)
    Rank = length(bnu)
    eta = cbind(c(Coefs@B1)) + Coefs@A %*% bnumat  # bix1 = intercept only
    M = misc.list$M
    for(s in 1:M) {
        temp = Coefs@D[,,s,drop=FALSE]
        dim(temp) = dim(temp)[1:2]  # c(Rank, Rank)
        eta[s,1] = eta[s,1] + t(bnumat) %*% temp %*% bnumat
    }
    eta = matrix(eta, 1, M, byrow=TRUE)
    mu = rbind(mu.function(eta, extra))  # Make sure it has one row 
    value = objfun(mu=mu, y=y,
                   w=1,  # ignore prior.weights on the object
                   residuals=FALSE, eta=eta, extra=extra)
    if(everything) {
        vcmat = matrix(0, Rank, Rank)
        for(s in 1:M) {
            vec1 = cbind(Coefs@A[s,]) + 2 *
                   matrix(Coefs@D[,,s], Rank, Rank) %*% bnumat
            vcmat = vcmat + mu[1,s] * vec1 %*% t(vec1)
        }
        vcmat = solve(vcmat)
    } else vcmat = NULL
    if(everything) list(eta=eta, mu=mu, value=value, vcmat=vcmat) else value
}

.my.calib.objfunction.cao = function(bnu, y, extra=NULL,
                        objfun, object, Coefs,
                        misc.list,
                        everything=TRUE,
                        mu.function) {
    Rank = length(bnu)
    NOS = Coefs@NOS 
    eta = matrix(as.numeric(NA), 1, NOS)
    for(j in 1:NOS) {
        eta[1,j] = predictcao(object, grid=bnu, sppno=j, 
                              Rank=Rank, deriv=0)$yvals
    }
    mu = rbind(mu.function(eta, extra))  # Make sure it has one row 
    value = objfun(mu=mu, y=y,
                   w=1,  # ignore prior.weights on the object
                   residuals=FALSE, eta=eta, extra=extra)
    vcmat = NULL  # No theory as of yet to compute the vcmat
    if(everything) list(eta=eta, mu=mu, value=value, vcmat=vcmat) else value
}


setMethod("calibrate", "qrrvglm", function(object, ...)
          calibrate.qrrvglm(object, ...))


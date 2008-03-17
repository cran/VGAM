# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.





vlabel <- function(xn, ncolBlist, M, separator=":") {

    if(length(xn) != length(ncolBlist))
        stop("length of first two arguments not equal")

    n1 <- rep(xn, ncolBlist)
    if(M==1)
        return(n1)
    n2 <- as.list(ncolBlist)
    n2 <- lapply(n2, seq)
    n2 <- unlist(n2)
    n2 <- as.character(n2)
    n2 <- paste(separator, n2, sep="")
    n3 <- rep(ncolBlist, ncolBlist)
    n2[n3==1] <- ""
    n1n2 <- paste(n1, n2, sep="")
    n1n2
}


lm2vlm.model.matrix <- function(x, Blist=NULL, assign.attributes=TRUE,
                                M=NULL, xij=NULL, Aarray=NULL, Aindex=NULL)
{
    


    if(length(Blist) != ncol(x))
        stop("length(Blist) != ncol(x)")

    if(length(xij)) {
        if(inherits(xij, "formula"))
            xij = list(xij)
        if(!is.list(xij))
            stop("xij is not a list of formulae")
    }

    if(!is.numeric(M))
        M <- nrow(Blist[[1]])

    if(length(xij)) {
        Blist.NAed = Blist 
        atx = attr(x, "assign")
        for(i in 1:length(xij)) {
            form = xij[[i]]
            if(length(form) != 3) 
                stop(paste("xij[[", i, "]] is not a formula with a response"))
            tform = terms(form)
            atform = attr(tform, "term.labels")   # doesn't include response 
            if(length(atform) != M) {
                stop(paste("xij[[", i, "]] does not contain", M, " terms"))
            }
            for(k in 1:length(atform)) {
                for(s in atx[[(atform[[k]])]]) {
                    if(length(Blist[[s]])) {
                        Blist[[s]] = ei(k, M)   # Easy for later
                        Blist.NAed[[s]] = Blist[[s]] * NA    # NA'ed 
                    }
                }
            }
        }
    }
    n <- nrow(x)
    if(all(trivial.constraints(Blist)) && !length(Aarray)) {
        xbig <- if(M > 1) kronecker(x, diag(M)) else x
        ncolBlist <- rep(M, ncol(x))
    } else {
        allB <- matrix(unlist(Blist), nrow=M)
        ncolBlist <- unlist(lapply(Blist, ncol))
        R <- sum(ncolBlist)

        X1 <- rep(c(t(x)), rep(ncolBlist,n))
        dim(X1) <- c(R, n)
        BB <- kronecker(matrix(1,n,1), allB)
        if(length(Aarray)) { 
            tmp34 = aperm(Aarray, c(1,3,2)) # c(M,n,r)
            for(ii in 1:length(Aindex))
            BB[,Aindex[[ii]]] = c(tmp34)
        }
        xbig <- kronecker(t(X1), matrix(1,M,1)) * BB
    }

    dn <- labels(x)
    yn <- dn[[1]]
    xn <- dn[[2]]
    dimnames(xbig) <- list(vlabel(yn, rep(M, n), M), 
                           vlabel(xn, ncolBlist, M))

    if(assign.attributes) {
    
        attr(xbig, "contrasts")   <- attr(x, "contrasts")
        attr(xbig, "factors")     <- attr(x, "factors")
        attr(xbig, "formula")     <- attr(x, "formula")
        attr(xbig, "class")       <- attr(x, "class")
        attr(xbig, "order")       <- attr(x, "order")
        attr(xbig, "term.labels") <- attr(x, "term.labels")
    

        nasgn <- oasgn <- attr(x, "assign")
        low <- 0
        for(i in 1:length(oasgn)) {
            len <- length(oasgn[[i]]) * ncolBlist[oasgn[[i]][1]]
            nasgn[[i]] <- (low+1):(low+len)
            low = low + len
        }
        if(low != ncol(xbig))
            stop("something gone wrong")
        attr(xbig, "assign") <- nasgn
    

        fred <- unlist(lapply(nasgn, length)) / unlist(lapply(oasgn, length))
        vasgn <- vector("list", sum(fred))
        k <- 0
        for(i in 1:length(oasgn)) {
            temp <- matrix(nasgn[[i]], ncol=length(oasgn[[i]]))
            for(j in 1:nrow(temp)) {
                k <- k + 1
                vasgn[[k]] <- temp[j,]
            }
        }
        names(vasgn) <- vlabel(names(oasgn), fred, M)
        attr(xbig, "vassign") <- vasgn


        attr(xbig, "constraints") <- Blist
    }


    xasgn <- attr(x, "assign")

    if(length(xij)) {
        rm.col.index = NULL    # Remove these columns from xbig
        for(i in 1:length(xij)) {
            form = xij[[i]]  # deparse(form1[[3]]) 
            tform = terms(form)
            atform = attr(tform, "term.labels")   # doesn't include response 
            response.name = (dimnames(attr(tform, "factors"))[[1]])[1]
    
            ptr0 = NULL 
            for(s in 1:M)
                if(length(nasgn[[atform[s]]])) {
                    ptr0 = s 
                    break 
                }
            if(!is.numeric(ptr0)) stop("no destination column indices")
            dest.col.index = nasgn[[atform[ptr0]]]

            if(M > 1)
                for(k in ((1:M)[-ptr0])) {
                    from.col.index = nasgn[[atform[k]]]   # May be NULL 
                    if(length(from.col.index)) {
                        xbig[,dest.col.index] = xbig[,dest.col.index] +
                                                xbig[,from.col.index]

                        rm.col.index = c(rm.col.index, from.col.index)

                        vasgn[[atform[k]]] = NULL    # Delete it
                    }
                }

            d2 = dimnames(xbig)[[2]]
            d2[dest.col.index] = vlabel(response.name, 
                                 length(dest.col.index), M=M, separator="")
            dimnames(xbig) = list(dimnames(xbig)[[1]], d2) 

            ptr = (1:length(names(vasgn)))[(names(vasgn)==atform[[ptr0]])]
            names(vasgn)[ptr] = response.name

        }

        if(length(rm.col.index))
            xbig = xbig[,-rm.col.index,drop=FALSE] # Delete the columns in 1 go

        if(assign.attributes) {
            attr(xbig, "constraints") <- Blist.NAed   # Not quite right
            attr(xbig, "vassign") <- vasgn
            attr(xbig, "assign") <- nasgn
            attr(xbig, "xij") <- xij
        }

    }


    xbig
}


model.matrixvlm = function(object, type=c("vlm","lm"), ...) {

    if(mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
    type <- match.arg(type, c("vlm","lm"))[1]

    x <- slot(object, "x")
    if(!length(x)) {
        data = model.frame(object, xlev=object@xlevels, ...) 

        kill.con = if(length(object@contrasts)) object@contrasts else NULL

        x <- vmodel.matrix.default(object, data=data,
                                   contrasts.arg = kill.con)
        if(is.R()) {

if(TRUE) {
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
}
            tt = terms(object)
            attr(x, "assign") <- attrassigndefault(x, tt)
        }
    }



    if(type == "lm") {
        return(x)
    } else {
        M <- object@misc$M  
        Blist <- object@constraints # Is NULL if there were no constraints?
        lm2vlm.model.matrix(x=x, Blist=Blist, xij=object@control$xij)
    }
}




setMethod("model.matrix",  "vlm", function(object, ...)
           model.matrixvlm(object, ...))






 if(is.R()) {

model.framevlm = function(object, ...) {

    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
    if(length(nargs) || !length(object@model)) {
        fcall <- object@call
        fcall$method <- "model.frame"
        fcall[[1]] <- as.name("vlm")

        fcall$smart <- FALSE
        if(length(object@smart.prediction)) {
            setup.smart("read", smart.prediction=object@smart.prediction)
        }

        fcall[names(nargs)] <- nargs
        env <- environment(object@terms$terms) # @terms or @terms$terms ??
        if (is.null(env)) 
            env <- parent.frame()
        ans = eval(fcall, env, parent.frame())

        if(length(object@smart.prediction)) {
            wrapup.smart()
        }

        ans
    } else object@model
}


if(!isGeneric("model.frame"))
    setGeneric("model.frame", function(formula, ...)
        standardGeneric("model.frame"))

setMethod("model.frame",  "vlm", function(formula, ...)
           model.framevlm(object=formula, ...))

}




vmodel.matrix.default = function (object, data = environment(object), 
    contrasts.arg = NULL, xlev = NULL, ...) {

    t <- terms(object)
    if (is.null(attr(data, "terms"))) 
        data <- model.frame(object, data, xlev = xlev) else
    {
        reorder <- match(sapply(attr(t, "variables"), deparse, 
            width.cutoff = 500)[-1], names(data))
        if (any(is.na(reorder))) 
            stop("model frame and formula mismatch in model.matrix()")
        data <- data[, reorder, drop = FALSE]
    }
    int <- attr(t, "response")
    if (length(data)) {
        contr.funs <- as.character(getOption("contrasts"))
        isF <- sapply(data, function(x) is.factor(x) || is.logical(x))
        isF[int] <- FALSE
        isOF <- sapply(data, is.ordered)
        namD <- names(data)
        for (nn in namD[isF]) if (is.null(attr(data[[nn]], "contrasts"))) 
            contrasts(data[[nn]]) <- contr.funs[1 + isOF[nn]]
        if (!is.null(contrasts.arg) && is.list(contrasts.arg)) {
            if (is.null(namC <- names(contrasts.arg))) 
                stop("invalid contrasts argument")
            for (nn in namC) {
                if (is.na(ni <- match(nn, namD))) 
                  warning(paste("Variable", nn,
                      "absent, contrast ignored")) else {
                  ca <- contrasts.arg[[nn]]
                  if (is.matrix(ca)) 
                    contrasts(data[[ni]], ncol(ca)) <- ca else
                    contrasts(data[[ni]]) <- contrasts.arg[[nn]]
                }
            }
        }
    } else {
        isF <- FALSE
        data <- list(x = rep(0, nrow(data)))
    }
    ans <- .Internal(model.matrix(t, data))
    cons <- if (any(isF)) 
        lapply(data[isF], function(x) attr(x, "contrasts")) else NULL
    attr(ans, "contrasts") <- cons
    ans
}









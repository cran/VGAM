# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.











 attrassigndefault = function(mmat, tt) {
    if (!inherits(tt, "terms"))
        stop("need terms object")
    aa = attr(mmat, "assign")
    if (is.null(aa))
        stop("argument is not really a model matrix")
    ll = attr(tt, "term.labels")
    if (attr(tt, "intercept") > 0)
        ll = c("(Intercept)", ll)
    aaa = factor(aa, labels = ll)
    split(order(aa), aaa)
}


 attrassignlm = function(object, ...)
     attrassigndefault(model.matrix(object), object@terms)



 vlabel = function(xn, ncolBlist, M, separator=":") {

    if (length(xn) != length(ncolBlist))
        stop("length of first two arguments not equal")

    n1 = rep(xn, ncolBlist)
    if (M == 1)
        return(n1)
    n2 = as.list(ncolBlist)
    n2 = lapply(n2, seq)
    n2 = unlist(n2)
    n2 = as.character(n2)
    n2 = paste(separator, n2, sep="")
    n3 = rep(ncolBlist, ncolBlist)
    n2[n3==1] = ""
    n1n2 = paste(n1, n2, sep="")
    n1n2
}


 lm2vlm.model.matrix = function(x, Blist=NULL, assign.attributes=TRUE,
                                M=NULL, xij=NULL, Xm2=NULL) {




    if (length(Blist) != ncol(x))
        stop("length(Blist) != ncol(x)")

    if (length(xij)) {
        if (inherits(xij, "formula"))
            xij = list(xij)
        if (!is.list(xij))
            stop("'xij' is not a list of formulae")
    }

    if (!is.numeric(M))
        M = nrow(Blist[[1]])

    nrow_X_lm = nrow(x)
    if (all(trivial.constraints(Blist) == 1)) {
        X_vlm = if (M > 1) kronecker(x, diag(M)) else x
        ncolBlist = rep(M, ncol(x))
    } else {
        allB = matrix(unlist(Blist), nrow=M)
        ncolBlist = unlist(lapply(Blist, ncol))
        Rsum = sum(ncolBlist)

        X1 = rep(c(t(x)), rep(ncolBlist, nrow_X_lm))
        dim(X1) = c(Rsum, nrow_X_lm)
        X_vlm = kronecker(t(X1), matrix(1, M, 1)) *
                kronecker(matrix(1, nrow_X_lm, 1), allB)
        rm(X1)
    }

    dn = labels(x)
    yn = dn[[1]]
    xn = dn[[2]]
    dimnames(X_vlm) = list(vlabel(yn, rep(M, nrow_X_lm), M), 
                           vlabel(xn, ncolBlist, M))

    if (assign.attributes) {
        attr(X_vlm, "contrasts")   = attr(x, "contrasts")
        attr(X_vlm, "factors")     = attr(x, "factors")
        attr(X_vlm, "formula")     = attr(x, "formula")
        attr(X_vlm, "class")       = attr(x, "class")
        attr(X_vlm, "order")       = attr(x, "order")
        attr(X_vlm, "term.labels") = attr(x, "term.labels")
    
        nasgn = oasgn = attr(x, "assign")
        lowind = 0
        for(ii in 1:length(oasgn)) {
            mylen = length(oasgn[[ii]]) * ncolBlist[oasgn[[ii]][1]]
            nasgn[[ii]] = (lowind+1):(lowind+mylen)
            lowind = lowind + mylen
        } # End of ii
        if (lowind != ncol(X_vlm))
            stop("something gone wrong")
        attr(X_vlm, "assign") = nasgn
    

        fred = unlist(lapply(nasgn, length)) / unlist(lapply(oasgn, length))
        vasgn = vector("list", sum(fred))
        kk = 0
        for(ii in 1:length(oasgn)) {
            temp = matrix(nasgn[[ii]], ncol=length(oasgn[[ii]]))
            for(jloc in 1:nrow(temp)) {
                kk = kk + 1
                vasgn[[kk]] = temp[jloc,]
            }
        }
        names(vasgn) = vlabel(names(oasgn), fred, M)
        attr(X_vlm, "vassign") = vasgn

        attr(X_vlm, "constraints") = Blist
    } # End of if (assign.attributes)




    if (!length(xij)) return(X_vlm)







    at.x = attr(x, "assign")
    at.vlmx = attr(X_vlm, "assign")
    at.Xm2 = attr(Xm2, "assign")

    for(ii in 1:length(xij)) {
        form.xij = xij[[ii]]
        if (length(form.xij) != 3) 
            stop("xij[[", ii, "]] is not a formula with a response")
        tform.xij = terms(form.xij)
        aterm.form = attr(tform.xij, "term.labels") # Does not include response
        if (length(aterm.form) != M)
            stop("xij[[", ii, "]] does not contain ", M, " terms")

        name.term.y = as.character(form.xij)[2]
        cols.X_vlm = at.vlmx[[name.term.y]]  # May be > 1 in length.

        x.name.term.2 = aterm.form[1]   # Choose the first one
        One.such.term = at.Xm2[[x.name.term.2]]
        for(bbb in 1:length(One.such.term)) {
            use.cols.Xm2 = NULL
            for(sss in 1:M) {
                x.name.term.2 = aterm.form[sss]
                one.such.term = at.Xm2[[x.name.term.2]]
                use.cols.Xm2 = c(use.cols.Xm2, one.such.term[bbb])
            } # End of sss

            allXk = Xm2[,use.cols.Xm2,drop=FALSE]
            cmat.no = (at.x[[name.term.y]])[1] # First one will do (all the same).
            cmat = Blist[[cmat.no]]
            Rsum.k = ncol(cmat)
            tmp44 = kronecker(matrix(1, nrow_X_lm, 1), t(cmat)) *
                    kronecker(allXk, matrix(1,ncol(cmat), 1)) # n*Rsum.k x M

            tmp44 = array(t(tmp44), c(M, Rsum.k, nrow_X_lm))
            tmp44 = aperm(tmp44, c(1,3,2)) # c(M, n, Rsum.k)
            rep.index = cols.X_vlm[((bbb-1)*Rsum.k+1):(bbb*Rsum.k)]
            X_vlm[,rep.index] = c(tmp44) 
        } # End of bbb
    } # End of for(ii in 1:length(xij))

    if (assign.attributes) {
        attr(X_vlm, "vassign") = vasgn
        attr(X_vlm, "assign") = nasgn
        attr(X_vlm, "xij") = xij
    }
    X_vlm
}






 model.matrixvlm = function(object, type=c("vlm","lm","lm2","bothlmlm2"),
                            ...) {



    if (mode(type) != "character" && mode(type) != "name")
    type = as.character(substitute(type))
    type = match.arg(type, c("vlm","lm","lm2","bothlmlm2"))[1]



    x = slot(object, "x")
    Xm2 = slot(object, "Xm2")

    if (!length(x)) {
        data = model.frame(object, xlev=object@xlevels, ...) 

        kill.con = if (length(object@contrasts)) object@contrasts else NULL

        x = vmodel.matrix.default(object, data=data,
                                  contrasts.arg = kill.con)
        tt = terms(object)
        attr(x, "assign") = attrassigndefault(x, tt)
    }

    if ((type == "lm2" || type == "bothlmlm2") && !length(Xm2)) {
        object.copy2 = object
        data = model.frame(object.copy2, xlev=object.copy2@xlevels, ...) 

        kill.con = if (length(object.copy2@contrasts))
                   object.copy2@contrasts else NULL

        Xm2 = vmodel.matrix.default(object.copy2, data=data,
                                    contrasts.arg = kill.con)
        ttXm2 = terms(object.copy2@misc$form2)
        attr(Xm2, "assign") = attrassigndefault(Xm2, ttXm2)
    }



    if (type == "lm") {
        return(x)
    } else if (type == "lm2") {
        return(Xm2)
    } else if (type == "bothlmlm2") {
        return(list(X=x, Xm2=Xm2))
    } else {
        M = object@misc$M  
        Blist = object@constraints # Is NULL if there were no constraints?
        lm2vlm.model.matrix(x=x, Blist=Blist, xij=object@control$xij, Xm2=Xm2)
    }
}




setMethod("model.matrix",  "vlm", function(object, ...)
           model.matrixvlm(object, ...))







 model.framevlm = function(object, 
                           setupsmart=TRUE, wrapupsmart=TRUE, ...) {

    dots = list(...)
    nargs = dots[match(c("data", "na.action", "subset"), names(dots), 0)]
    if (length(nargs) || !length(object@model)) {
        fcall = object@call
        fcall$method = "model.frame"
        fcall[[1]] = as.name("vlm")

        fcall$smart = FALSE
        if (setupsmart && length(object@smart.prediction)) {
            setup.smart("read", smart.prediction=object@smart.prediction)
        }

        fcall[names(nargs)] = nargs
        env = environment(object@terms$terms) # @terms or @terms$terms ??
        if (is.null(env)) 
            env = parent.frame()
        ans = eval(fcall, env, parent.frame())

        if (wrapupsmart && length(object@smart.prediction)) {
            wrapup.smart()
        }
        ans
    } else object@model
}


if (!isGeneric("model.frame"))
    setGeneric("model.frame", function(formula, ...)
        standardGeneric("model.frame"))

setMethod("model.frame",  "vlm", function(formula, ...)
           model.framevlm(object=formula, ...))





 vmodel.matrix.default = function(object, data = environment(object),
                                  contrasts.arg = NULL, xlev = NULL, ...) {
 print("20120221; in vmodel.matrix.default")
    t <- if (missing(data)) terms(object) else terms(object, data = data)
    if (is.null(attr(data, "terms")))
        data <- model.frame(object, data, xlev = xlev) else {
        reorder <- match(sapply(attr(t, "variables"), deparse,
            width.cutoff = 500)[-1], names(data))
        if (any(is.na(reorder)))
            stop("model frame and formula mismatch in model.matrix()")
        if (!identical(reorder, seq_len(ncol(data))))
            data <- data[, reorder, drop = FALSE]
    }
    int <- attr(t, "response")
    if (length(data)) {
        contr.funs <- as.character(getOption("contrasts"))
        namD <- names(data)
        for (i in namD) if (is.character(data[[i]])) {
            data[[i]] <- factor(data[[i]])
            warning(gettextf("variable '%s' converted to a factor",
                i), domain = NA)
        }
        isF <- sapply(data, function(x) is.factor(x) || is.logical(x))
        isF[int] <- FALSE
        isOF <- sapply(data, is.ordered)
        for (nn in namD[isF]) if (is.null(attr(data[[nn]], "contrasts")))
            contrasts(data[[nn]]) <- contr.funs[1 + isOF[nn]]
        if (!is.null(contrasts.arg) && is.list(contrasts.arg)) {
            if (is.null(namC <- names(contrasts.arg)))
                stop("invalid 'contrasts.arg' argument")
            for (nn in namC) {
                if (is.na(ni <- match(nn, namD)))
                  warning(gettextf(
                    "variable '%s' is absent, its contrast will be ignored",
                    nn), domain = NA) else {
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


    ans  <-          (model.matrix(t, data))




    cons <- if (any(isF))
        lapply(data[isF], function(x) attr(x, "contrasts")) else NULL
    attr(ans, "contrasts") <- cons
    ans
}




depvar.vlm <- function(object, ...) {
  object@y
}



if (!isGeneric("depvar"))
    setGeneric("depvar", function(object, ...) standardGeneric("depvar"),
               package = "VGAM")


setMethod("depvar",  "vlm", function(object, ...)
           depvar.vlm(object, ...))
setMethod("depvar",  "rrvglm", function(object, ...)
           depvar.vlm(object, ...))
setMethod("depvar",  "qrrvglm", function(object, ...)
           depvar.vlm(object, ...))
setMethod("depvar",  "cao", function(object, ...)
           depvar.vlm(object, ...))
setMethod("depvar",  "rcam", function(object, ...)
           depvar.vlm(object, ...))







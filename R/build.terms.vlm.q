# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.



    if (!isGeneric("terms"))
        setGeneric("terms", function(x, ...) standardGeneric("terms"))





terms.vlm = function(x, ...) {
    v = x@terms
    if (!length(v))
        stop("terms slot is empty")
    v = v$terms
    if (!length(v))
        stop("no terms component")
    v
}


setMethod("terms", "vlm", function(x, ...) terms.vlm(x, ...))





Build.terms.vlm = function(x, coefs, cov = NULL, assign, collapse = TRUE, M,
                           dimname=NULL, coefmat = NULL) {


    cov.true = !is.null(cov)
    if (collapse) {
        fit = matrix(x %*% coefs, ncol=M, byrow=TRUE)
        dimnames(fit) = dimname
        if (M==1)
            fit = c(fit)
        if (cov.true) {
            var = ((x %*% cov) * x) %*% rep(1, length(coefs))
            list(fitted.values = fit, se.fit = if (M==1) c(sqrt(var)) else 
                 matrix(sqrt(var), ncol=M, byrow=TRUE, dimnames=dimname))
        } else {
            fit
        }
    } else {

        constant = attr(x, "constant")
        if (!is.null(constant)) {
            constant = as.vector( t(coefmat) %*% constant )
        }
    
        if (missing(assign))
            assign = attr(x, "assign")
        if (is.null(assign))
            stop("Need an 'assign' list")
        fit = array(0, c(nrow(x), length(assign)),
                    list(dimnames(x)[[1]], names(assign)))
        if (cov.true)
            se = fit
        TL = sapply(assign, length)
        simple = TL == 1
        complex = TL > 1
        if (any(simple)) {
            asss = unlist(assign[simple])
            ones = rep(1, nrow(x))
            fit[, simple] = x[, asss] * outer(ones, coefs[asss])
            if (cov.true)
                se[,simple] = abs(x[,asss]) * outer(ones, sqrt(diag(cov))[asss])
        }
        if (any(complex)) {
            assign = assign[complex]
            for(term in names(assign)) {
                TT = assign[[term]]
                xt = x[, TT]
                fit[, term] = xt %*% coefs[TT]
                if (cov.true)
                  se[, term] = sqrt(drop(((xt %*% cov[TT, TT]) * xt) %*%
                               rep(1, length(TT))))
            }
        }
        attr(fit, "constant") = constant
    
        if (cov.true) list(fitted.values = fit, se.fit = se) else fit
    }
}




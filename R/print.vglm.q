# These functions are
# Copyright (C) 1998-2009 T.W. Yee, University of Auckland. All rights reserved.


print.vglm <- function(x, ...) {
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
    print.default(coef, ...)

    rank <- x@rank
    if(!length(rank))
        rank <- sum(!nas)
    nobs <- if(length(x@df.total)) x@df.total else length(x@residuals)
    rdf <- x@df.residual
    if(!length(rdf))
        rdf <- nobs - rank
    cat("\nDegrees of Freedom:", nobs, "Total;", rdf, "Residual\n")

    if(length(deviance(x)))
        cat("Residual Deviance:", format(deviance(x)), "\n")
    llx = logLik.vlm(object = x)

    if(length(llx))
        cat("Log-likelihood:", format(llx), "\n")

    if(length(x@criterion)) {
        ncrit <- names(x@criterion)
        for(i in ncrit)
            if(i!="loglikelihood" && i!="deviance")
                cat(paste(i, ":", sep=""), format(x@criterion[[i]]), "\n")
    }

    invisible(x)
}


print.vgam <- function(x, digits=2, ...) {

    if(!is.null(cl <- x@call)) {
        cat("Call:\n")
        dput(cl)
    }

    coef <- x@coefficients
    nas <- is.na(coef)

    rank <- x@rank
    if(is.null(rank))
        rank <- sum(!nas)
    nobs <- if(length(x@df.total)) x@df.total else length(x@residuals)
    rdf <- x@df.residual
    if(is.null(rdf))
        rdf <- nobs - rank
    cat("\nDegrees of Freedom:", nobs, "Total;",
        format(round(rdf, dig=digits)), "Residual\n")

    if(length(deviance(x)))
        cat("Residual Deviance:", format(deviance(x)), "\n")

    llx = logLik.vlm(object = x)

    if(length(llx))
        cat("Log-likelihood:", format(llx), "\n")

    criterion <- attr(terms(x), "criterion")  # 11/8/03; x@terms$terms,
    if(!is.null(criterion) && criterion!="coefficients")
        cat(paste(criterion, ":", sep=""), format(x[[criterion]]), "\n")

    invisible(x)
}




setMethod("print",  "vlm", function(x, ...)  print.vlm(x, ...))
setMethod("print", "vglm", function(x, ...) print.vglm(x, ...))
setMethod("print", "vgam", function(x, ...) print.vgam(x, ...))


    setMethod("show",  "vlm", function(object)  print.vlm(object))
    setMethod("show", "vglm", function(object) print.vglm(object))
    setMethod("show", "vgam", function(object) print.vgam(object))



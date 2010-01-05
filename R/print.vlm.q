# These functions are
# Copyright (C) 1998-2010 T.W. Yee, University of Auckland. All rights reserved.



if (!is.R()) {
setMethod("show", "vlm",
    function(object)
    print.vlm(object))
}

setMethod("print", "vlm",
    function(x, ...)
    print.vlm(x, ...))

print.vlm <- function(x, ...) {
    if (!is.null(cl <- x@call)) {
        cat("Call:\n")
        dput(cl)
    }

    coef <- x@coefficients
    cat("\nCoefficients:\n")
    print(coef, ...)

    rank <- x@rank
    if (is.null(rank))
        rank <- sum(!is.na(coef))
    n <- x@misc$n 
    M <- x@misc$M 
    nobs <- if (length(x@df.total)) x@df.total else n*M
    rdf <- x@df.residual
    if (is.null(rdf))
        rdf <- (n - rank) * M
    cat("\nDegrees of Freedom:", nobs, "Total;", rdf, "Residual\n")

    if (length(deviance(x)) && is.finite(deviance(x)))
        cat("Deviance:", format(deviance(x)), "\n")
    if (length(x@rss) && is.finite(x@rss))
        cat("Residual Sum of Squares:", format(x@rss), "\n")

    invisible(x)
}




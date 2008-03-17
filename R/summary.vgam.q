# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.





summaryvgam <- function(object, dispersion=NULL, digits=options()$digits-2)
{

    if(length(dispersion) && dispersion == 0 &&
       length(object@family@summary.dispersion) &&
       !object@family@summary.dispersion) {
        stop(paste("Can't use the general VGLM formula (based on a residual",
                   "sum of squares) for computing the dispersion parameter."))
    }

    newobject <- object 
    class(newobject) <- "vglm"
    stuff <- summaryvglm(newobject, dispersion=dispersion)
    rdf <- stuff@df[2] <- object@df.residual  # NA 

    M <- object@misc$M
    n.big <- object@misc$n.big
    rank <- if(is.null(object@qr$rank)) length(object@coefficients) else
            object@qr$rank








    useF <- object@misc$useF
    if(is.null(useF))
        useF <- FALSE

    df <- unlist(lapply(object@misc$new.assign, length))
    nldf <- object@nl.df

    if(length(df)) {
        aod <- as.matrix(round(df, 1))
        dimnames(aod) <- list(names(df), "Df")
        if(!is.null(object@nl.chisq)) {
            aod <- cbind(aod, NA, NA, NA)
            nl.chisq <- object@nl.chisq / object@dispersion

            special = abs(nldf) < 0.1  # This was the quick fix in s.vam()  
            nldf[special] = 1          # Give it a plausible value for pchisq & pf

            snames <- names(nldf)
            aod[snames, 2] <- round(nldf, 1)
            aod[snames, 3] <- if(useF) nl.chisq/nldf  else nl.chisq
            aod[snames, 4] <- if(useF) 1-pf(nl.chisq/nldf, nldf, rdf) else 
                1-pchisq(nl.chisq, nldf)

            if(any(special)) {
                aod[snames[special], 2:4] = NA 
            }

            rnames <- c("Df", "Npar Df", "Npar Chisq", "P(Chi)")
            if(useF)
                rnames[3:4] <- c("Npar F", "Pr(F)")
            dimnames(aod) <- list(names(df), rnames)
        heading <- if(useF)
        "\nDF for Terms and Approximate F-values for Nonparametric Effects\n"
        else
      "\nDF for Terms and Approximate Chi-squares for Nonparametric Effects\n"
        } else heading <- "DF for Terms\n\n"
        aod <- as.vanova(data.frame(aod, check.names=FALSE), heading)

        if(is.R()) class(aod) = "data.frame"
    }
    else aod <- if(is.R()) data.frame() else NULL

    answer <-
    new("summary.vgam",
        object,
        call=stuff@call,
        cov.unscaled=stuff@cov.unscaled,
        correlation=stuff@correlation,
        df=stuff@df,
        sigma=stuff@sigma)

    slot(answer, "coefficients") = stuff@coefficients  # Replace
    if(is.numeric(stuff@dispersion))
        slot(answer, "dispersion") = stuff@dispersion

    presid = residuals(object, type="pearson")
    if(length(presid))
        answer@pearson.resid= as.matrix(presid)

        slot(answer, "anova") = aod 

    answer
}




printsummary.vgam <- function(x, quote=TRUE, prefix="", digits=options()$digits-2)
{

    M <- x@misc$M


    cat("\nCall:\n")
    dput(x@call)

    presid <- x@pearson.resid
    rdf <- x@df[2]
    if(FALSE && !is.null(presid) && all(!is.na(presid))) {
        cat("\nPearson Residuals:\n")
        if(rdf/M > 5) {
            rq <-  apply(as.matrix(presid), 2, quantile) # 5 x M
            dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", "Max"),
                                 x@misc$predictors.names)
            print(t(rq), digits=digits)
        } else
        if(rdf > 0) {
            print(presid, digits=digits)
        }
    }

    cat("\nNumber of linear predictors:   ", M, "\n")

    if(!is.null(x@misc$predictors.names))
    if(M==1) 
        cat("\nName of linear predictor:",
            paste(x@misc$predictors.names, collapse=", "), "\n") else if(M<=5)
        cat("\nNames of linear predictors:",
            paste(x@misc$predictors.names, collapse=", "), "\n")

    prose <- ""
    if(length(x@dispersion)) {
        if(is.logical(x@misc$estimated.dispersion) &&
           x@misc$estimated.dispersion)
            prose <- "(Estimated) " else {

            if(is.numeric(x@misc$default.dispersion) &&
               x@dispersion==x@misc$default.dispersion)
                prose <- "(Default) "

            if(is.numeric(x@misc$default.dispersion) &&
               x@dispersion!=x@misc$default.dispersion)
                prose <- "(Pre-specified) "
        }
        cat(paste("\n", prose, "Dispersion Parameter for ",
            x@family@vfamily[1],
            " family:   ", format(round(x@dispersion, digits)), "\n", sep=""))
    }

    if(length(deviance(x)))
        cat("\nResidual Deviance: ", format(round(deviance(x), digits)),
            "on", format(round(rdf, 3)), "degrees of freedom\n")
    if(length(logLik(x)))
        cat("\nLog-likelihood:", format(round(logLik(x), digits)),
            "on", format(round(rdf, 3)), "degrees of freedom\n")

    if(length(x@criterion)) {
        ncrit <- names(x@criterion)
        for(i in ncrit)
            if(i!="loglikelihood" && i!="deviance")
                cat(paste(i, ":", sep=""), format(x@criterion[[i]]), "\n")
    }


    cat("\nNumber of Iterations: ", x@iter, "\n")

    if(length(x@anova)) {
        printvanova(x@anova, dig=digits)   # ".vanova" for Splus6
    }

    invisible(NULL)
}




    setMethod("summary", "vgam",
             function(object, ...)
             summaryvgam(object, ...))

    setMethod("print", "summary.vgam",
             function(x, ...)
             printsummary.vgam(x, ...))


    setMethod("show", "summary.vgam",
             function(object)
             printsummary.vgam(object))



 
 
printvanova <- function(x, digits=.Options$digits, ...)
{
    rrr <- row.names(x) 
    heading <- attr(x, "heading")
    if(!is.null(heading))
        cat(heading, sep="\n")
    attr(x, "heading") <- NULL
    for(i in 1:length(x)) {
        xx <- x[[i]]
        xna <- is.na(xx)
        xx <- format(zapsmall(xx, digits))
        xx[xna] <- ""
        x[[i]] <- xx
    }
    if(is.R()) {
        print.data.frame(as.data.frame(x, row.names=rrr))
        invisible(x)
    } else {
        print.data.frame(as.data.frame(x, row.names=rrr))
        invisible(x)
    }
}

as.vanova <- function(x, heading)
{
    if(!is.data.frame(x))
        stop("x must be a data frame")
    rrr <- row.names(x) 
    attr(x, "heading") <- heading
    if(is.R()) { 
        x <- as.data.frame(x, row.names=rrr)
    } else {
        x <- as.data.frame(x, row.names=rrr)
    }
    x
}


if(!is.R()) {

setMethod("print", "vanova",
          function(x, ...)
          printvanova(x, ...))

}




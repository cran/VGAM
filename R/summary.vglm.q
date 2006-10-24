# These functions are
# Copyright (C) 1998-2006 T.W. Yee, University of Auckland. All rights reserved.







yformat = function(x, digits=options()$digits) {
    format(ifelse(abs(x)<0.001, signif(x, digits), round(x, digits)))
}




summaryvglm <- function(object, correlation=FALSE, dispersion=NULL, digits=NULL)
{



    if(length(dispersion) && dispersion == 0 && 
       length(object@family@summary.dispersion) && 
       !object@family@summary.dispersion) {
        stop(paste("Can't use the general VGLM formula (based on a residual",
                   "sum of squares) for computing the dispersion parameter.")) 
    }

    stuff <- summaryvlm(as(object, "vlm"),
                         correlation=correlation,
                         dispersion=dispersion)



    answer <-
    new("summary.vglm",
        object,
        coef3=stuff@coef3,
        cov.unscaled=stuff@cov.unscaled,
        correlation=stuff@correlation,
        df=stuff@df,
        sigma=stuff@sigma)

    presid = resid(object, type="pearson")
    if(length(presid))
        answer@pearson.resid = as.matrix(presid)

    slot(answer, "misc") = stuff@misc  # Replace

    if(is.numeric(stuff@dispersion))
        slot(answer, "dispersion") = stuff@dispersion

    answer
}






setMethod("logLik",  "summary.vglm", function(object, ...)
    logLik.vlm(object, ...))


printsummary.vglm <- function(x, digits = NULL, quote = TRUE, prefix = "")
{

    M <- x@misc$M
    coef <- x@coef3   # icients
    correl <- x@correlation

    digits <- if(is.null(digits)) options()$digits - 2 else digits

    cat("\nCall:\n")
    dput(x@call)

    presid <- x@pearson.resid
    rdf <- x@df[2]
    if(length(presid) && all(!is.na(presid)) && is.finite(rdf))
    {
        cat("\nPearson Residuals:\n")
        if(rdf/M > 5) 
        {
            rq <-  apply(as.matrix(presid), 2, quantile) # 5 x M
            dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", "Max"),
                                 x@misc$predictors.names)
            print(t(rq), digits = digits)
        } else
        if(rdf > 0) {
            print(presid, digits = digits)
        }
    }

    cat("\nCoefficients:\n")
    print.default(coef, digits = digits)

    cat("\nNumber of linear predictors: ", M, "\n")

    if(!is.null(x@misc$predictors.names))
    if(M==1) 
        cat("\nName of linear predictor:",
            paste(x@misc$predictors.names, collapse=", "), "\n") else if(M<=5)
        cat("\nNames of linear predictors:",
            paste(x@misc$predictors.names, collapse=", "), fill=TRUE)

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
            " family:   ", yformat(x@dispersion, digits), "\n", sep=""))
    }

    if(length(deviance(x))) {
        cat("\nResidual Deviance:", yformat(deviance(x), digits))
        if(is.finite(rdf))
            cat(" on", round(rdf, digits), "degrees of freedom\n") else
            cat("\n")
    }
    if(length(logLik(x))) {
        cat("\nLog-likelihood:", yformat(logLik(x), digits))
        if(is.finite(rdf))
            cat(" on", round(rdf, digits), "degrees of freedom\n") else
            cat("\n")
    }

    if(length(x@criterion)) {
        ncrit <- names(x@criterion)
        for(i in ncrit)
            if(i!="loglikelihood" && i!="deviance")
                cat(paste(i, ":", sep=""), yformat(x@criterion[[i]], digits),
                    "\n")
    }


    cat("\nNumber of Iterations:", format(trunc(x@iter)), "\n")

    if(!is.null(correl)) 
    {
        p.big <- dim(correl)[2]
        if(p.big > 1) 
        {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits))
            correl[!ll] <- ""
            print(correl[-1,  -p.big, drop = FALSE], quote = FALSE, digits = 
                digits)
        }
    }
    invisible(NULL)
}



    setMethod("summary", "vglm",
             function(object, ...)
             summaryvglm(object, ...))

    setMethod("print", "summary.vglm",
             function(x, ...)
             invisible(printsummary.vglm(x, ...)))

    setMethod("show", "summary.vglm",
             function(object)
             invisible(printsummary.vglm(object)))







vcovdefault <- function(object, ...) {
    if(is.null(object@vcov))
        stop("no default")
    object@vcov
}

vcovvlm <- function(object, dispersion=NULL, untransform=FALSE) {
    so <- summaryvlm(object, corr=FALSE, dispersion=dispersion)
    d = if(any(slotNames(so) == "dispersion") && 
           is.Numeric(so@dispersion)) so@dispersion else 1
    answer = d * so@cov.unscaled

    if(!untransform) return(answer)

    if(!is.logical(object@misc$intercept.only))
       stop(paste("cannot determine whether the object is",
                  "an intercept-only fit, i.e., y ~ 1 is the response"))
    if(!object@misc$intercept.only)
       stop("object must be an intercept-only fit, i.e., y ~ 1 is the response")

    M = object@misc$M
    Links = object@misc$link
    if(length(Links) != M && length(Links) != 1)
       stop("cannot obtain the link functions to untransform the object")


    tvector = numeric(M)
    etavector = predict(object)[1,]   # Contains transformed parameters
    earg = object@misc$earg  # This could be a NULL
    if(!is.null(earg) && M > 1 && (!is.list(earg) || length(earg) != M))
        stop(paste("the earg component of object@misc should be of length ",
                   M, sep=""))
    for(ii in 1:M) {
        TTheta = etavector[ii]  # Transformed theta
        use.earg = if(M == 1 || is.null(earg)) earg else earg[[ii]]
        if(is.list(use.earg) && !length(use.earg))
            use.earg = NULL
        newcall = paste(Links[ii],
                        "(theta=TTheta, earg=use.earg, inverse=TRUE)", sep="")
        newcall = parse(text=newcall)[[1]]
        Theta = eval(newcall) # Theta, the untransformed parameter
        newcall = paste(Links[ii],
                        "(theta=Theta, earg=use.earg, deriv=1)", sep="")
        newcall = parse(text=newcall)[[1]]
        tvector[ii] = eval(newcall)
    }
    tvector = abs(tvector)
    answer = (cbind(tvector) %*% rbind(tvector)) * answer
    if(length(dmn2 <- names(object@misc$link)) == M)
        dimnames(answer) = list(dmn2, dmn2)
    answer
}

setMethod("vcov", "vlm",
         function(object, ...)
         vcovvlm(object, ...))

setMethod("vcov", "vglm",
         function(object, ...)
         vcovvlm(object, ...))






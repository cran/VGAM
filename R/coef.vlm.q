# These functions are
# Copyright (C) 1998-2011 T.W. Yee, University of Auckland.
# All rights reserved.






coefvlm <- function(object, matrix.out=FALSE, label=TRUE) {

    ans <- object@coefficients
    if (!label)
        names(ans) <- NULL
    if (!matrix.out)
        return(ans)

 
    ncolx <- object@misc$p   # = length(object@constraints)
    M <- object@misc$M

    Blist <- object@constraints
    if (all(trivial.constraints(Blist) == 1)) {
        Bmat <- matrix(ans, nrow=ncolx, ncol=M, byrow=TRUE)
    } else {
        Bmat <- matrix(as.numeric(NA), nrow=ncolx, ncol=M)

        if (!matrix.out)
            return(ans) 

        ncolBlist <- unlist(lapply(Blist, ncol)) 
        nasgn <- names(Blist)
        temp <- c(0, cumsum(ncolBlist))
        for(ii in 1:length(nasgn)) {
            index <- (temp[ii]+1):temp[ii+1]
            cmat <- Blist[[nasgn[ii]]]
            Bmat[ii,] <- cmat %*% ans[index]
        }
    }

    if (label) {
        d1 <- object@misc$colnames.x
        d2 = object@misc$predictors.names # Could be NULL
        dimnames(Bmat) <- list(d1, d2)
    }

    Bmat
} # end of coefvlm




setMethod("coefficients", "vlm", function(object, ...)
           coefvlm(object, ...))
setMethod("coef", "vlm", function(object, ...)
           coefvlm(object, ...))
setMethod("coefficients", "vglm", function(object, ...)
           coefvlm(object, ...))
setMethod("coef", "vglm", function(object, ...)
           coefvlm(object, ...))



  
setMethod("coefficients", "summary.vglm", function(object, ...) object@coef3)
setMethod("coef",         "summary.vglm", function(object, ...) object@coef3)









Coef.vlm <- function(object, ...) {
    LL <- length(object@family@vfamily)
    funname = paste("Coef.", object@family@vfamily[LL], sep="")
    if (exists(funname)) {
        newcall = paste("Coef.", object@family@vfamily[LL],
                        "(object, ...)", sep="")
        newcall = parse(text=newcall)[[1]]
        eval(newcall)
    } else
    if (length(tmp2 <- object@misc$link) &&
       object@misc$intercept.only &&
       trivial.constraints(object@constraints)) {

        answer = eta2theta(rbind(coefvlm(object)),
                           link = object@misc$link,
                           earg = object@misc$earg)
        answer = c(answer)
        if (length(ntmp2 <- names(tmp2)) == object@misc$M)
            names(answer) = ntmp2
        answer
    } else {
        coefvlm(object, ... )
    }
}

setMethod("Coefficients", "vlm", function(object, ...)
               Coef.vlm(object, ...))
setMethod("Coef", "vlm", function(object, ...)
               Coef.vlm(object, ...))







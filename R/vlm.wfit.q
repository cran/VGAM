# These functions are
# Copyright (C) 1998-2009 T.W. Yee, University of Auckland. All rights reserved.








vlm.wfit <- function(x, z, Blist, wz=NULL, U=NULL, 
        matrix.out=FALSE, XBIG=FALSE, rss=TRUE, qr=FALSE, x.ret=FALSE,
        offset=NULL,
        omit.these=NULL, only.rss=FALSE,
        ncolx=if(matrix.out && XBIG) stop("need this argument") else ncol(x),
        xij=NULL, Aarray=NULL, Aindex=NULL, lp.names=NULL, 
        Eta.range=NULL, ...)
{

    missing.Blist <- missing(Blist)
    z = as.matrix(z)
    n <- nrow(z)
    M <- ncol(z)
    if(!only.rss) {
        contrast.save <- attr(x, "contrasts")
        znames <- dimnames(z)[[2]]
    }

    if(length(offset))
        z <- z - offset
    if(missing(U) || !length(U)) {
        U <- vchol(wz, M=M, n=n, silent=FALSE)
    }
    dU <- dim(U)
    if(dU[2] != n)
        stop("input unconformable")

    xbig.save <- if(XBIG) {
            x 
        } else {
            if(missing.Blist || !length(Blist))
                Blist = replace.constraints(vector("list", ncol(x)), 
                                            diag(M), 1:ncol(x)) # NULL
            lm2vlm.model.matrix(x=x, Blist=Blist, M=M, assign.attributes=FALSE,
                                xij = xij, Aarray=Aarray, Aindex=Aindex)
        }
    xbig <- mux111(U, xbig.save, M=M)
    z.big <- mux22(U, z, M=M, upper=TRUE, as.mat=FALSE)


    if(length(omit.these)) {
        xbig = xbig[!omit.these,,drop=FALSE] 
        z.big = z.big[!omit.these]
    }

    ans <- if(!is.R()) lm.fit.qr(x=xbig, y=z.big, qr=qr, ...) else  
               lm.fit(xbig, z.big, ...)

    if(rss) {
        ans$rss <- sum(ans$resid^2)
        if(only.rss) return(list(rss=ans$rss))
    }


    if(length(omit.these) && any(omit.these))
        stop("code beyond here can't handle omitted observations")


    fv <- ans$fitted.values
    dim(fv) <- c(M, n)
    fv <- vbacksub(U, fv, M=M, n=n) # Have to premultiply fv by U


    if(length(Eta.range)) {
        if(length(Eta.range) != 2)
            stop("length(Eta.range) must equal 2")
        fv = ifelse(fv < Eta.range[1], Eta.range[1], fv)
        fv = ifelse(fv > Eta.range[2], Eta.range[2], fv)
    }

    ans$fitted.values <- if(M==1) c(fv) else fv
    if(M > 1)
        dimnames(ans$fitted.values) <- list(dimnames(z)[[1]], znames)
    ans$residuals <- if(M==1) c(z-fv) else z-fv
    if(M > 1)
        dimnames(ans$residuals) <- list(dimnames(ans$residuals)[[1]], znames)
    ans$misc <- list(M=M, n=n)
    ans$call <- match.call()

    ans$constraints <- Blist
    ans$contrasts <- contrast.save
    if(x.ret) 
        ans$xbig <- xbig.save

    if(!is.null(offset))
        ans$fitted.values <- ans$fitted.values + offset




    if(!matrix.out)
        return(ans)


    dx2 = if(XBIG) NULL else dimnames(x)[[2]]
    B <- matrix(as.numeric(NA), nrow=M, ncol=ncolx, dimnames=list(lp.names, dx2))
    if(is.null(Blist)) {
        Blist = replace.constraints(vector("list", ncolx), diag(M), 1:ncolx)
    }
    ncolBlist <- unlist(lapply(Blist, ncol)) 
    temp <- c(0, cumsum(ncolBlist))
    for(i in 1:ncolx) {
        index <- (temp[i]+1):temp[i+1]
        cm <- Blist[[i]]
        B[,i] <- cm %*% ans$coef[index]
    }
    ans$mat.coefficients <- t(B)
    ans
}


print.vlm.wfit <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }

    coef <- x$coefficients
    cat("\nCoefficients:\n")
    print(coef, ...)

    rank <- x$rank
    if(is.null(rank))
        rank <- sum(!is.na(coef))
    n <- x$misc$n 
    M <- x$misc$M 
    rdf <- x$df.resid
    if(is.null(rdf))
        rdf <- (n - rank) * M
    cat("\nDegrees of Freedom:", n*M, "Total;", rdf, "Residual\n")

    if(!is.null(x$rss))
        cat("Residual Sum of Squares:", format(x$rss), "\n")

    invisible(x)
}




# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.









vlm.wfit <- function(xmat, zmat, Blist, wz = NULL, U = NULL, 
    matrix.out = FALSE, is.vlmX = FALSE, rss = TRUE, qr = FALSE,
    x.ret = FALSE,
    offset = NULL,
    omit.these = NULL, only.rss = FALSE,
    ncolx = if (matrix.out && is.vlmX) {
        stop("need argument 'ncolx'") 
    } else {
            ncol(xmat)
    },
    xij = NULL,
    lp.names = NULL, Eta.range = NULL, Xm2 = NULL, ...) {
    missing.Blist <- missing(Blist)
    zmat = as.matrix(zmat)
    n <- nrow(zmat)
    M <- ncol(zmat)
    if (!only.rss) {
        contrast.save <- attr(xmat, "contrasts")
        znames <- dimnames(zmat)[[2]]
    }

    if (length(offset)) {
        zmat <- zmat - offset
    }
    if (missing(U) || !length(U)) {
        U <- vchol(wz, M = M, n = n, silent = FALSE)
    }
    dU <- dim(U)
    if (dU[2] != n) {
        stop("input unconformable")
    }

    X_vlm_save <- if (is.vlmX) {
            xmat 
        } else {
            if (missing.Blist || !length(Blist)) {
                Blist = replace.constraints(vector("list", ncol(xmat)),
                                            diag(M), 1:ncol(xmat)) # NULL
            }
            lm2vlm.model.matrix(x=xmat, Blist=Blist, M = M,
                                assign.attributes = FALSE,
                                xij = xij,
                                Xm2=Xm2)
        }
    X_vlm <- mux111(U, X_vlm_save, M = M)
    z_vlm <- mux22(U, zmat, M = M, upper = TRUE, as.matrix = FALSE)


    if (length(omit.these)) {
        X_vlm = X_vlm[!omit.these,,drop = FALSE] 
        z_vlm = z_vlm[!omit.these]
    }

    ans <- lm.fit(X_vlm, z_vlm, ...)

    if (rss) {
        ans$rss <- sum(ans$resid^2)
        if (only.rss) return(list(rss = ans$rss))
    }

    if (length(omit.these) && any(omit.these)) {
        stop("code beyond here cannot handle omitted observations")
    }


    fv <- ans$fitted.values
    dim(fv) <- c(M, n)
    fv <- vbacksub(U, fv, M = M, n = n) # Have to premultiply fv by U


    if (length(Eta.range)) {
        if (length(Eta.range) != 2) {
            stop("length(Eta.range) must equal 2")
        }
        fv = ifelse(fv < Eta.range[1], Eta.range[1], fv)
        fv = ifelse(fv > Eta.range[2], Eta.range[2], fv)
    }

    ans$fitted.values <- if (M == 1) c(fv) else fv
    if (M > 1) {
        dimnames(ans$fitted.values) <- list(dimnames(zmat)[[1]], znames)
    }
    ans$residuals <- if (M == 1) c(zmat-fv) else zmat-fv
    if (M > 1) {
        dimnames(ans$residuals) <- list(dimnames(ans$residuals)[[1]], znames)
    }
    ans$misc <- list(M = M, n = n)
    ans$call <- match.call()

    ans$constraints <- Blist
    ans$contrasts <- contrast.save
    if (x.ret) {
        ans$X_vlm <- X_vlm_save
    }

    if (!is.null(offset)) {
        ans$fitted.values <- ans$fitted.values + offset
    }




    if (!matrix.out) {
        return(ans)
    }


    dx2 = if (is.vlmX) NULL else dimnames(xmat)[[2]]
    B = matrix(as.numeric(NA),
               nrow = M, ncol = ncolx, dimnames = list(lp.names, dx2))
    if (is.null(Blist)) {
        Blist = replace.constraints(vector("list", ncolx), diag(M), 1:ncolx)
    }
    ncolBlist <- unlist(lapply(Blist, ncol)) 
    temp <- c(0, cumsum(ncolBlist))
    for(ii in 1:ncolx) {
        index <- (temp[ii]+1):temp[ii+1]
        cm <- Blist[[ii]]
        B[,ii] <- cm %*% ans$coef[index]
    }
    ans$mat.coefficients <- t(B)
    ans
}




if (FALSE)
print.vlm.wfit <- function(x, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }

  coef <- x$coefficients
  cat("\nCoefficients:\n")
  print(coef, ...)

  rank <- x$rank
  if (is.null(rank)) {
    rank <- sum(!is.na(coef))
  }
  n <- x$misc$n 
  M <- x$misc$M 
  rdf <- x$df.resid
  if (is.null(rdf)) {
    rdf <- (n - rank) * M
  }
  cat("\nDegrees of Freedom:", n*M, "Total;", rdf, "Residual\n")

  if (!is.null(x$rss)) {
    cat("Residual Sum of Squares:", format(x$rss), "\n")
  }

  invisible(x)
}




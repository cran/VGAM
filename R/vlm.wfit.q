# These functions are
# Copyright (C) 1998-2014 T.W. Yee, University of Auckland.
# All rights reserved.









vlm.wfit <-
  function(xmat, zmat, Hlist, wz = NULL, U = NULL, 
           matrix.out = FALSE, is.vlmX = FALSE, res.ss = TRUE, qr = FALSE,
           x.ret = FALSE,
           offset = NULL,
           omit.these = NULL, only.res.ss = FALSE,
           ncolx = if (matrix.out && is.vlmX) {
                     stop("need argument 'ncolx'") 
                   } else {
                     ncol(xmat)
                   },
           xij = NULL,
           lp.names = NULL, Eta.range = NULL, Xm2 = NULL, ...) {


  missing.Hlist <- missing(Hlist)
  zmat <- as.matrix(zmat)
  n <- nrow(zmat)
  M <- ncol(zmat)
  if (!only.res.ss) {
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

  X.vlm.save <- if (is.vlmX) {
        xmat 
      } else {
        if (missing.Hlist || !length(Hlist)) {
          Hlist <- replace.constraints(vector("list", ncol(xmat)),
                                       diag(M), 1:ncol(xmat))  # NULL
        }
        lm2vlm.model.matrix(x = xmat, Hlist = Hlist, M = M,
                            assign.attributes = FALSE,
                            xij = xij,
                            Xm2 = Xm2)
      }
  X.vlm <- mux111(U, X.vlm.save, M = M)
  z.vlm <- mux22(U, zmat, M = M, upper = TRUE, as.matrix = FALSE)


  if (length(omit.these)) {
    X.vlm <- X.vlm[!omit.these, , drop = FALSE] 
    z.vlm <- z.vlm[!omit.these]
  }






  ans <- lm.fit(X.vlm, y = z.vlm, ...)

  if (res.ss) {
    ans$res.ss <- sum(ans$resid^2)
    if (only.res.ss)
      return(list(res.ss = ans$res.ss))
  }

  if (length(omit.these) && any(omit.these)) {
    stop("code beyond here cannot handle omitted observations")
  }


  fv <- ans$fitted.values
  dim(fv) <- c(M, n)
  fv <- vbacksub(U, fv, M = M, n = n)  # Have to premultiply fv by U


  if (length(Eta.range)) {
    if (length(Eta.range) != 2) {
      stop("length(Eta.range) must equal 2")
    }
    fv <- ifelse(fv < Eta.range[1], Eta.range[1], fv)
    fv <- ifelse(fv > Eta.range[2], Eta.range[2], fv)
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

  ans$constraints <- Hlist
  ans$contrasts <- contrast.save
  if (x.ret) {
    ans$X.vlm <- X.vlm.save
  }

  if (!is.null(offset)) {
    ans$fitted.values <- ans$fitted.values + offset
  }




  if (!matrix.out) {
    return(ans)
  }


  dx2 <- if (is.vlmX) NULL else dimnames(xmat)[[2]]
  B <- matrix(as.numeric(NA),
             nrow = M, ncol = ncolx, dimnames = list(lp.names, dx2))
  if (is.null(Hlist)) {
    Hlist <- replace.constraints(vector("list", ncolx), diag(M), 1:ncolx)
  }
  ncolHlist <- unlist(lapply(Hlist, ncol)) 
  temp <- c(0, cumsum(ncolHlist))
  for (ii in 1:ncolx) {
    index <- (temp[ii]+1):temp[ii+1]
    cm <- Hlist[[ii]]
    B[, ii] <- cm %*% ans$coef[index]
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

  if (!is.null(x$res.ss)) {
    cat("Residual Sum of Squares:", format(x$res.ss), "\n")
  }

  invisible(x)
}




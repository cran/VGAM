# These functions are
# Copyright (C) 1998-2016 T.W. Yee, University of Auckland.
# All rights reserved.











vlm.wfit <-
  function(xmat, zmat, Hlist, wz = NULL, U = NULL, 
           matrix.out = FALSE, is.vlmX = FALSE, ResSS = TRUE, qr = FALSE,
           x.ret = FALSE,
           offset = NULL,
           omit.these = NULL, only.ResSS = FALSE,
           ncolx = if (matrix.out && is.vlmX) {
                     stop("need argument 'ncolx'") 
                   } else {
                     ncol(xmat)
                   },
           xij = NULL,
           lp.names = NULL, Eta.range = NULL, Xm2 = NULL,

           Xvlm.aug = NULL,
           ps.list = NULL,
           constraints = NULL, first.ps = FALSE,
           trace = FALSE,

           ...) {
  mgcvvgam <- length(ps.list)

  missing.Hlist <- missing(Hlist)
  zmat <- as.matrix(zmat)
  n <- nrow(zmat)
  M <- ncol(zmat)
  if (!only.ResSS) {
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










  if (mgcvvgam) {
    m.objects <- psv2magic(x.VLM = X.vlm,
                           constraints = constraints,
                           lambda.vlm = attr(Xvlm.aug, "lambda.vlm"),
                           ps.list = ps.list)
    if (FALSE && trace) {
      cat("m.objects$sp \n")
      print( m.objects$sp )
      cat("m.objects$off \n")
      print( m.objects$off )
    }

    if (first.ps) {
      m.objects$sp <- rep_len(-1, length(m.objects$sp))
    }



    magicfit <- mgcv::magic(y   = z.vlm,
                            X   = m.objects$x.VLM.new,
                            sp  = m.objects$sp,
                            S   = m.objects$S.arg,
                            off = m.objects$off,
                            gcv = FALSE)
    SP <- magicfit$sp
    if (FALSE && trace) {
      cat("SP \n")
      print( SP )
    }

    length.lambda.vlm <- sapply(attr(Xvlm.aug, "lambda.vlm"), length)  # lambda.new
    sp.opt <- vector("list", length(length.lambda.vlm))  # list()
    iioffset <- 0
    for (ii in seq_along(length.lambda.vlm)) {
      sp.opt[[ii]] <- SP[iioffset + 1:length.lambda.vlm[ii]]
      iioffset <- iioffset + length.lambda.vlm[ii]
    }
    names(sp.opt) <- names(ps.list$which.X.ps)
    if (FALSE && trace) {
      cat("sp.opt \n")
      print( sp.opt )
    }

    ps.list$lambdalist <- sp.opt
    Xvlm.aug <- Pen.psv(constraints = constraints, ps.list = ps.list)


    first.ps <- FALSE  # May have been TRUE on entry but is FALSE on exit


    X.vlm <- rbind(X.vlm, Xvlm.aug)
    z.vlm <- c(z.vlm, rep(0, nrow(Xvlm.aug)))
  }




  ans <- lm.fit(X.vlm, y = z.vlm, ...)



  if (mgcvvgam) {
    ans$residuals     <- head(ans$residuals,     n*M)
    ans$effects       <- head(ans$effects,       n*M)
    ans$fitted.values <- head(ans$fitted.values, n*M)
    ans$qr$qr         <- head(ans$qr$qr,         n*M)
  }





  if (ResSS) {
    ans$ResSS <- sum(ans$resid^2)
    if (only.ResSS)
      return(list(ResSS = ans$ResSS))
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








  if (mgcvvgam) {
    ans$first.ps <- first.ps  # Updated.
    ans$ps.list  <- ps.list   # Updated wrt "lambdalist" component.
    ans$Xvlm.aug <- Xvlm.aug  # Updated matrix.
    ans$magicfit <- magicfit  # Updated.
  }



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
  B <- matrix(NA_real_,
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
}  # vlm.wfit





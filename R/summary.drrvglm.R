# These functions are
# Copyright (C) 1998-2024 T.W. Yee, University of Auckland.
# All rights reserved.










 summary.drrvglm <-
  function(object, correlation = FALSE,
           dispersion = NULL, digits = NULL,
           numerical = TRUE,
           h.step = 0.005,
           omit123 = FALSE,
           omit13 = FALSE,   # TRUE
           fixA = FALSE,
           presid = FALSE,  # TRUE
 signif.stars = getOption("show.signif.stars"),
           nopredictors = FALSE,
           eval0 = TRUE, ...) {








  object@control$trace <- FALSE  # Suppress


  if (!is.Numeric(h.step, length.arg = 1) ||
      abs(h.step) > 1 || h.step == 0)
    stop("bad input for 'h.step'")


  if (is.null(dispersion))
    dispersion <- object@misc$dispersion


  newobject <- as(object, "vglm")
  stuff <- summaryvglm(newobject,
                       correlation = correlation,
                       dispersion = dispersion,
                       presid = presid)


  if ((length(dispersion) &&
       any(dispersion != 1)) ||
      (length(stuff@dispersion) &&
       any(stuff@dispersion != 1)))
    stop("VGAM no longer supports dispersion ",
         "parameters")
  dispersion <- 1  # For beyond here ,,,,,,,,


  answer <-
    new(Class = "summary.drrvglm",
        object,
        call = stuff@call,
        coef3 = stuff@coef3,
        cov.unscaled = stuff@cov.unscaled,
        correlation = stuff@correlation,
        df = stuff@df,
        sigma = stuff@sigma)



  if (presid && length(stuff@pearson.resid))
    slot(answer, "pearson.resid") <-
    stuff@pearson.resid


  tmp5 <-
    get.drrvglm.se1(object, omit13 = omit13,
      numerical = numerical, h.step = h.step,
      omit123 = omit123, fixA = fixA, ...)
  if (eval0 && (
      any(diag(tmp5$cov.unscaled) <= 0) ||
      any(eigen(tmp5$cov.unscaled,
          symmetric = TRUE)$value <= 0))) {
    warning("cov.unscaled is not pos-definite")
  }  # eval0
  answer@cov.unscaled <- tmp5$cov.unscaled



  answer@df[1] <- answer@df[1] + tmp5$n.elts.tA
  answer@df[2] <- answer@df[2] - tmp5$n.elts.tA

  answer@coef3 <-
    get.drrvglm.se2(answer@cov.unscaled,
              dispersion = dispersion,
              coefficients = tmp5$allcoefs)

  answer@dispersion <- dispersion
  answer@sigma <- sqrt(dispersion)


answer@misc$signif.stars <- signif.stars  # 201606
answer@misc$nopredictors <- nopredictors  # 201509

  answer
}  # summary.drrvglm







setClass("summary.drrvglm",
         representation("summary.rrvglm",
                        H.A = "list",
                        H.C = "list"))



setMethod("summary", "drrvglm",
         function(object, ...)
         summary.drrvglm(object, ...))



show.summary.drrvglm <-
  function(x, digits = NULL, quote = TRUE,
           prefix = "",
           signif.stars = NULL) {
  show(as(x, "summary.rrvglm"))

  invisible(x)
  NULL
}

 setMethod("show", "summary.drrvglm",
           function(object)
             show.summary.drrvglm(x = object))






 get.drrvglm.se1 <-
  function(object,
           omit13 = FALSE, omit123 = FALSE,
           numerical = TRUE,
           fixA = FALSE, h.step = 0.0001,
           trace.arg = FALSE,
           check.2 = FALSE,  # TRUE,
           ...) {



  F <- FALSE
  T <- TRUE
  rrcontrol <- object@control
  covun.RAvcov <- chol2inv(object@misc$RAvcov)
  covun.RCvcov <- chol2inv(object@misc$RCvcov)
  pell2.AB1 <- solve(covun.RAvcov)
  pell2.B1C <- solve(covun.RCvcov)
  if (check.2) {
    RA.ei <- eigen(covun.RAvcov, symmetric = TRUE)
    RC.ei <- eigen(covun.RCvcov, symmetric = TRUE)
    cat("covun.RAvcov is ",
        ifelse(all(RA.ei$val > 0), "", "NOT "),
        "positive-definite.\n", sep = "")
    cat("covun.RCvcov is ",
        ifelse(all(RC.ei$val > 0), "", "NOT "),
        "positive-definite.\n", sep = "")
  }  # check.2

  colx1.index <- rrcontrol$colx1.index
  colx2.index <- rrcontrol$colx2.index
  p1 <- length(colx1.index)  # May be 0
  p2 <- length(colx2.index)
  Rank <- object@control$Rank
  H.C <- if (length(object@misc$H.C))
    object@misc$H.C else object@H.C
  H.A.alt <- if (length(object@misc$H.A.alt))
    object@misc$H.A.alt else object@H.A.alt
  ncol.H.C <- sapply(H.C, ncol)
  ncol.H.A.alt <- sapply(H.A.alt, ncol)
  Hlist <- constraints(object)  # type = "term"?
  ncolHlist <- unlist(lapply(Hlist, ncol))
  ncolH.C <- unlist(lapply(H.C, ncol))

  ind3 <- seq(sum(ncol.H.A.alt))  # 1st subset: A
  pell.11   <- pell2.AB1[ ind3,  ind3, drop = F]
  B1a.pell2 <- pell2.AB1[-ind3, -ind3, drop = F]
  pell.12   <- pell2.AB1[ ind3, -ind3, drop = F]




  cs.ncolHlist <- cumsum(c(1, ncolHlist))
  ind.colx1 <- ind.colx2 <- NULL
  for (i in colx1.index)
    ind.colx1 <- c(ind.colx1,
    (cs.ncolHlist[i]):(cs.ncolHlist[i + 1] - 1))
  for (i in colx2.index)
    ind.colx2 <- c(ind.colx2,
    (cs.ncolHlist[i]):(cs.ncolHlist[i + 1] - 1))
  B1b.pell2 <- pell2.B1C[ind.colx1, ind.colx1,
                         drop = FALSE]
  C.pell2 <- pell2.B1C[ind.colx2, ind.colx2,
                       drop = FALSE]


  B1.check <- abs(max(B1a.pell2 - B1b.pell2))
  if (check.2)
    cat("B1.check (0?) is", B1.check, "\n")
  if (B1.check > 1e-5)
    warning("estimate of B1 differs substantia",
            "lly between two overlapping fits")


  pell.23 <- pell2.B1C[ind.colx1, ind.colx2,
                       drop = FALSE]
  pell.33 <- C.pell2


  X.lm <- if (length(object@x)) object@x else
    model.matrix(object, type = "lm")
  x1mat <- if (p1)
    X.lm[, colx1.index, drop = FALSE] else NULL
  x2mat <- X.lm[, colx2.index, drop = FALSE]
  Amat <- object@A.est
  Cmat <- object@C.est
  if (!length(M <- object@misc$M))
    M <- npred(object)
  str0 <- rrcontrol$str0
  Index.corner <- rrcontrol$Index.corner
  B1mat <- if (p1)
    coefvlm(object, matrix.out = TRUE)[
            colx1.index, , drop = FALSE] else NULL


  if (!omit13) {
  delct.da <- if (numerical) {
    num.deriv.drrr(object, M = M, Rank = Rank,
          x1mat = x1mat, x2mat = x2mat,
          p2 = p2, h.step = h.step,
          Index.corner = Index.corner,
          Aimat = Amat,
          B1mat = B1mat, Cimat = Cmat,
          H.A.alt = H.A.alt, H.C = H.C,
          ncol.H.A.alt = ncol.H.A.alt,
          ncol.H.C = ncol.H.C,
          xij = rrcontrol$xij,
          str0 = str0)
  } else {
    stop("dctda.fast.only() currently ",
         "unavailable for 'drrvglm' objects")
    warning("20240103; this call to ",
      "dctda.fast.only() needs work for drrvglm")
    thetA <- c(Amat[-c(Index.corner, str0), ])
    wz <- U <- zmat <- NULL  # Avoid a warning
    dctda.fast.only(theta = thetA, wz = wz,
                    U = U, zmat, M = M, 
                    r = Rank, x1mat = x1mat,
                    x2mat = x2mat, p2 = p2,
                    Index.corner, Aimat = Amat,
                    B1mat = B1mat, Cimat = Cmat,
                    xij = object@control$xij,
                    str0 = str0)
  }
  }  # !omit13


  pell.13 <- if (omit13)
    matrix(0, sum(ncol.H.A.alt),
              sum(ncolH.C)) else
    delct.da %*% (-pell.33)  # Need 2 mux by -1


  if (omit123) {
    pell.13 <- pell.13 * 0   # zero it
    if (fixA) {
      pell.12 <- pell.12 * 0   # zero it
    } else {
      pell.23 <- pell.23 * 0   # zero it
    }
  }  # omit123




  NEell2.partials <- if (fixA) {
    rbind(cbind(pell.11, pell.12, pell.13),
          cbind(rbind(t(pell.12), t(pell.13)),
                pell2.B1C))  # Huge blk @ bot RHS
  } else {  # fixC == T effectively
    rbind(cbind(pell2.AB1,  # Huge blk @ top LHS
                rbind(pell.13, pell.23)),
          cbind(t(pell.13), t(pell.23), pell.33))
  }


  cov.unscaled <- solve(NEell2.partials)


  prefx <- param.names("I(latvar.mat)", Rank, T)
  suffx <- c(sapply(ncol.H.A.alt, seq))
  Aelts.names <- if (Rank == 1) {
    paste(prefx, suffx, sep = ":")
  } else {
    iptr <- 1
    tmp5 <- rep(" ", length(suffx))
    cs.ncol.H.A.alt <- cumsum(ncol.H.A.alt)
    for (i in seq(length(suffx))) {
      tmp5[i] <- paste(prefx[iptr], suffx[i],
                       sep = ":")
      if (i >= cs.ncol.H.A.alt[iptr])
        iptr <- iptr + 1
    }
    tmp5
  }
  cnames <- c(Aelts.names,
              object@misc$colnames.X.vlm)
  dimnames(cov.unscaled) <- list(cnames, cnames)

  n.elts.tildeA <- if (is(object, "drrvglm"))
     sum(ncol.H.A.alt) else
     (M - Rank - length(str0)) * Rank
  allcoefs <- c(object@misc$Avec,
                object@coefficients)
  list(cov.unscaled  = cov.unscaled,
       allcoefs      = allcoefs,
       n.elts.tA     = n.elts.tildeA,
       ResSS         = object@ResSS)
}  # get.drrvglm.se1






get.drrvglm.se2 <-
    function(cov.unscaled, dispersion = 1,
             coefficients) {

  dn8 <-  dimnames(cov.unscaled)[[1]]
  ans <- matrix(coefficients,
                length(coefficients), 4)
  ans[, 2] <- sqrt(dispersion) *
              sqrt(diag(cov.unscaled))
  ans[, 3] <- ans[, 1] / ans[, 2]
  ans[, 4] <- pnorm(-abs(ans[, 3]))
  dimnames(ans) <-
    list(dn8, c("Estimate", "Std. Error",
                "z value", "Pr(>|z|)"))
  ans
}  # get.drrvglm.se2







 num.deriv.drrr <-
  function(object, M, Rank = 1, x1mat, x2mat,
    p2, Index.corner, Aimat, B1mat, Cimat,
    h.step = 0.0001,  # colx2.index,
    H.A.alt = list(), H.C = list(),  # new
    ncol.H.A.alt = rep(M-Rank, Rank),  # "rrvglm"
    ncol.H.C = rep(Rank, p2),  # "rrvglm"
    xij = NULL, str0 = NULL) {




  nn <- nrow(x2mat)
  if (nrow(Cimat) != p2 || ncol(Cimat) != Rank)
    stop("'Cimat' wrong shape")

  dct.da <- matrix(NA_real_,
              sum(ncol.H.A.alt), sum(ncol.H.C))

  cptr <- 1
  if (!length(B1Cvec <- object@misc$B1Cvec))
    stop("could not retrieve B1Cvec")
  for (vvv in 1:Rank) {
    for (ttt in 1:ncol.H.A.alt[vvv]) {
      small.Hlist <- vector("list", p2)
      pAmat <- Aimat
      pAmat[, vvv] <- pAmat[, vvv] + h.step *
        (H.A.alt[[vvv]])[, ttt]  # One coln
      for (ii in 1:p2)  # Only for x2mat
        small.Hlist[[ii]] <- pAmat

      offset <- if (length(object@offset))
                  object@offset else 0
      if (all(offset == 0)) offset <- 0
      neweta <- x2mat %*% Cimat %*% t(pAmat)
      if (is.numeric(x1mat))
        neweta <- neweta + x1mat %*% B1mat
      object@predictors <- neweta


      newmu <- object@family@linkinv(neweta,
                                     object@extra)
      object@fitted.values <- as.matrix(newmu)

      fred <- weights(object, type = "work",
                deriv = TRUE, ignore.slot = TRUE)
      if (!length(fred))
        stop("cannot get object@weights & @deriv")
      wz <- fred$weights
      deriv.mu <- fred$deriv

      U <- vchol(wz, M = M, n = nn, silent = TRUE)
      tvfor <- vforsub(U, as.matrix(deriv.mu),
                       M = M, n = nn)
      newzmat <- neweta - offset +
                 vbacksub(U, tvfor, M = M, n = nn)
      if (is.numeric(x1mat))
        newzmat <- newzmat - x1mat %*% B1mat
      newfit <- vlm.wfit(xmat = x2mat,
           zmat = newzmat, qr = FALSE,
           Hlist = small.Hlist, U = U, 
           matrix.out = FALSE, is.vlmX = FALSE,
           ResSS = TRUE, x.ret = FALSE,
           offset = NULL, xij = xij)
      dct.da[cptr, ] <-  # 1 elt at a time. 
        (tail(newfit$coef, sum(ncol.H.C)) -
         tail(B1Cvec, sum(ncol.H.C))) / h.step
      cptr <- cptr + 1
    }  # tt
  }  # ss

  dct.da
}  # num.deriv.drrr










setMethod("coefficients", "summary.drrvglm",
          function(object, ...)
          object@coef3)
setMethod("coef",         "summary.drrvglm",
          function(object, ...)
          object@coef3)































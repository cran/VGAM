# These functions are
# Copyright (C) 1998-2018 T.W. Yee, University of Auckland.
# All rights reserved.














hdeff.vglm <-
  function(object,
           derivative = NULL,
           se.arg = FALSE,
           subset = NULL,  # Useful for Cox model as a poissonff().
           hstep = 0.005,  # Formerly 'Delta', recycled to length 2
           fd.only = FALSE,
           ...) {


  if (is.Numeric(hstep, positive = TRUE)) {
      if (length(hstep) > 2)
        warning("length(hstep) too large; recycling to 2 values")
    hstep <- rep(hstep, length = 2)
    if (any(hstep > 0.1))
      warning("probably some values of 'hstep' are too large")
  } else {
    stop("bad input for argument 'hstep'")
  }



  type <- if (length(derivative)) {
    if (is.Numeric(derivative, length.arg = 1, positive = TRUE,
                   integer.valued = TRUE) &&
        derivative %in% 1:2)
        "derivatives" else stop("bad input for argument 'derivative'")
    } else {
      "logical"
    }



  Fam <- if (inherits(object, "vlm")) {
    object@family
  } else {
   stop("cannot get at the 'family' slot")
  }
  Fam.infos <- Fam@infos()

  if (is.logical(Fam.infos$hadof) && !Fam.infos$hadof)
    return(NULL)

  dfun <-
    if (is.logical(Fam.infos$hadof) && Fam.infos$hadof) {
      Fam@hadof
    } else {
      NULL  # This means its not implemented yet
    }



  link1parameter <- Fam.infos$link1parameter
  if (is.null(link1parameter))
    link1parameter <- TRUE  # The default really




  M <- npred(object)  # Some constraints span across responses
  ind5 <- iam(NA, NA, both = TRUE, M = M)
  MM12 <- M * (M + 1) / 2
  all.Hk <- constraints(object, matrix = TRUE)
  X.vlm <- model.matrix(object, type = "vlm")
  eta.mat <- predict(object)
  n.LM <- NROW(eta.mat)
  pwts <- weights(object, type = "prior")

  mylinks <- linkfun(object)  # Of length 1 for GLMs, char only
  wwt.0 <- weights(object, type = "working", ignore.slot = TRUE)
  if (ncol(wwt.0) < MM12)
    wwt.0 <- cbind(wwt.0, matrix(0, n.LM, MM12 - ncol(wwt.0)))
  dim.wz <- dim(wwt.0)  # Inefficient

  p.VLM <- ncol(all.Hk)
  M1 <- npred(object, type = "one.response")

  vc2 <- vcov(object)
  SE2 <- diag.ixwx <- diag(vc2)
  SE1 <- sqrt(SE2)
  cobj <- coef(object)
  SE2.deriv1 <- vec.Wald.deriv1 <- rep_len(NA_real_, p.VLM)
  names(vec.Wald.deriv1) <- names(cobj)
  if (type == "derivatives" && derivative == 2) {
    SE2.deriv2 <- vec.Wald.deriv2 <- vec.Wald.deriv1
  }


  fd.use <- is.null(dfun) || fd.only
  if ((blot.out <- !fd.use && any(colSums(all.Hk != 0) > 1)) &&
      (type == "derivatives" && derivative == 2))
    warning("2nd derivs available only when M1==1 and with ",
            "trivial constraints; ",
            "try setting 'fd.only = TRUE'; returning NAs")







  D3thetas.Detas3 <-  # May not be needed
  D2thetas.Detas2 <-
  D1thetas.Detas1 <-
  Param.mat       <- matrix(NA_real_, n.LM, M)
  if (link1parameter) {
    if (!fd.use) {
    for (jay in 1:M) {      # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      Param.mat[, jay] <-
      Param.vec <- eta2theta(eta.mat[, jay], mylinks[jay],
                             earg = object@misc$earg[[jay]])

      D1thetas.Detas1[, jay] <-
        dtheta.deta(Param.vec,
                    link = mylinks[jay],
                    earg = object@misc$earg[[jay]])
      D2thetas.Detas2[, jay] <-
        d2theta.deta2(Param.vec,
                      link = mylinks[jay],
                      earg = object@misc$earg[[jay]])


    if (type == "derivatives" && derivative == 2) {
      D3thetas.Detas3[, jay] <-
        d3theta.deta3(Param.vec,
                      link = mylinks[jay],
                      earg = object@misc$earg[[jay]])
    }


    }  # for (jay)
    wz.tet <-  D1thetas.Detas1[, ind5$row] *
               D1thetas.Detas1[, ind5$col]  # n x MM12
      Der1 <- D1thetas.Detas1
      Der2 <- D2thetas.Detas2
    }  # !fd.use
  } else {
    Param.mat <- eta2theta(eta.mat, mylinks,
                           earg = object@misc$earg)



    myearg <- object@misc$earg[[1]]
    build.list <- list(theta = Param.mat, inverse = TRUE, deriv = 1)
    build.list <- c(build.list, myearg)  # Hopefully no dups arg names
    build.list$all.derivs <- TRUE  # For multinomial, etc.
    Der1 <- do.call(what = mylinks, args = build.list)
    if (type == "derivatives" && derivative == 2) {
      build.list$deriv <- 2
      Der2 <- do.call(what = mylinks, args = build.list)
    }
  }  # if (link1parameter) and (!link1parameter)




  kvec.use <- 1:p.VLM
  if (length(subset))
    kvec.use <- kvec.use[subset]  # & !is.na(subset)




  for (kay in kvec.use) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  dwz.dbetas   <-
  d2wz.dbetas2 <- 0  # Good for the first instance of use.


  wetas.kay <- which.etas(object, kay = kay)
  for (jay in wetas.kay) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    vecTF.jay <- as.logical(eijfun(jay, M))
    bix.jk <- X.vlm[vecTF.jay, kay]  # An n-vector; allows for xij.



    if (fd.use) {
      for (dirr in 1:2) {
        temp1 <- object
        temp1@predictors[, jay] <- temp1@predictors[, jay] +
                                   (c(1, -1)[dirr]) * hstep[1]
        temp1@fitted.values <- cbind(
          temp1@family@linkinv(eta = temp1@predictors,
                               extra = temp1@extra))  # Make sure a matrix
        if (dirr == 1)
          wwt.f <- weights(temp1, type = "working", ignore.slot = TRUE)
        if (dirr == 2)
          wwt.b <- weights(temp1, type = "working", ignore.slot = TRUE)
      }  # dirr
      if (ncol(wwt.f) < MM12)
        wwt.f <- cbind(wwt.f, matrix(0, n.LM, MM12 - ncol(wwt.f)))
      if (ncol(wwt.b) < MM12)
        wwt.b <- cbind(wwt.b, matrix(0, n.LM, MM12 - ncol(wwt.b)))


      cdiff <- (wwt.f - wwt.b) / (2 * hstep[1])
      cdiff2 <- (wwt.f - 2 * wwt.0 + wwt.b) / (hstep[1]^2)

      dwz.dbetas   <- dwz.dbetas   + as.matrix(cdiff)  * bix.jk
      d2wz.dbetas2 <- d2wz.dbetas2 + as.matrix(cdiff2) * bix.jk^2



    } else {  # !fd.use
    if (link1parameter) {
      dfun1 <- dfun(eta.mat, extra = object@extra,
                    linpred.index = jay,
                    w = pwts, dim.wz = dim.wz,
                    deriv = 1)
      dfun0 <- dfun(eta.mat, extra = object@extra,
                    linpred.index = jay,
                    w = pwts, dim.wz = dim.wz,
                    deriv = 0)
      use.ncol <- NCOL(dfun0)  # Reference value really

        if (NCOL(dfun1 ) < use.ncol)
          dfun1  <- cbind(dfun1, matrix(0, n.LM, use.ncol - NCOL(dfun1)))
        if (use.ncol < NCOL(wz.tet))
          wz.tet <- wz.tet[, 1:use.ncol, drop = FALSE]

      write.into.wz <- function(jay, nxM) {
        M <- NCOL(nxM)
        wz <- matrix(0, NROW(nxM), M*(M+1)/2)
        for (uuu in 1:M)
          wz[, iam(jay, uuu, M = M)] <- (1 + (jay == uuu)) * nxM[, uuu]
        wz
      }

      wz.lhs <- write.into.wz(jay, D1thetas.Detas1)
      if (use.ncol < NCOL(wz.lhs))
        wz.lhs <- wz.lhs[, 1:use.ncol, drop = FALSE]
      dwz.dtheta.Der1 <- dfun1 * wz.tet * Der1[, jay] +
        Der2[, jay] * dfun0 * wz.lhs

      if (!is.matrix(dwz.dtheta.Der1))
        dwz.dtheta.Der1 <- as.matrix(dwz.dtheta.Der1)
    }  # else {






    if (link1parameter) {
      dwz.dbetakk <- dwz.dtheta.Der1 * bix.jk  # * Der1[, jay]
    } else {
      dwz.dbetakk <- 0
      for (uuu in 1:M) {
        dfun1 <- dfun(eta.mat, extra = object@extra,
                      linpred.index = uuu,
                      w = pwts, dim.wz = dim.wz,
                      deriv = 1)
        dwz.dbetakk <- dwz.dbetakk +
          dfun1 * Der1[, iam(uuu, jay, M = M)]
      }  # for uuu
      dwz.dbetakk <- dwz.dbetakk * bix.jk
    }

    if (!is.matrix(dwz.dbetakk))
      dwz.dbetakk <- as.matrix(dwz.dbetakk)
    dwz.dbetas <- dwz.dbetas + dwz.dbetakk  # Summed over 1:M

    }  # !fd.use



    if (type == "derivatives" && derivative == 2) {


    if (fd.use) {
  all.bix.jk.mat <- matrix(X.vlm[, kay], n.LM, M, byrow = TRUE)
  crossprod.bix.jk.mat <- all.bix.jk.mat[, ind5$row.index] *
                          all.bix.jk.mat[, ind5$col.index]

  if (length(wetas.kay) > 1) {
  for (sss in wetas.kay) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  for (ttt in wetas.kay) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
 if (sss < ttt) {
    for (dirr in 1:4) {
      temp1 <- object
      temp1@predictors[, sss] <- temp1@predictors[, sss] +
          (c(1, 1, -1, -1)[dirr]) * hstep[1]
      temp1@predictors[, ttt] <- temp1@predictors[, ttt] +
          (c(1, -1, 1, -1)[dirr]) * hstep[2]
      temp1@fitted.values <- cbind(
      temp1@family@linkinv(eta = temp1@predictors,
                           extra = temp1@extra))  # Make sure a matrix
      if (dirr == 1)
        wwt.1 <- weights(temp1, type = "working", ignore.slot = TRUE)
      if (dirr == 2)
        wwt.2 <- weights(temp1, type = "working", ignore.slot = TRUE)
      if (dirr == 3)
        wwt.3 <- weights(temp1, type = "working", ignore.slot = TRUE)
      if (dirr == 4)
        wwt.4 <- weights(temp1, type = "working", ignore.slot = TRUE)
    }  # dirr
    if (ncol(wwt.1) < MM12)
      wwt.1 <- cbind(wwt.1, matrix(0, n.LM, MM12 - ncol(wwt.1)))
    if (ncol(wwt.2) < MM12)
      wwt.2 <- cbind(wwt.2, matrix(0, n.LM, MM12 - ncol(wwt.2)))
    if (ncol(wwt.3) < MM12)
      wwt.3 <- cbind(wwt.3, matrix(0, n.LM, MM12 - ncol(wwt.3)))
    if (ncol(wwt.4) < MM12)
      wwt.4 <- cbind(wwt.4, matrix(0, n.LM, MM12 - ncol(wwt.4)))
    cdiff2 <- ((wwt.1 - wwt.2 - wwt.3 + wwt.4) / (4 *
              hstep[1]  * hstep[2]))

    d2wz.dbetas2 <- d2wz.dbetas2 + 2 *  # Twice
      as.matrix(cdiff2) *
      crossprod.bix.jk.mat[, iam(sss, ttt, M = M)]
  }  # if (sss < ttt)
  }  # ttt in wetas.kay
  }  # sss in wetas.kay
  }  # length(wetas.kay) > 1

    } else {  # !fd.use

      if (link1parameter) {
        Der3 <- D3thetas.Detas3
      }
      dfun2 <- dfun(eta.mat, extra = object@extra,
                    linpred.index = jay,
                    w = pwts, dim.wz = dim.wz,
                    deriv = 2)
      use.ncol <- if (link1parameter) NCOL(dwz.dtheta.Der1) else MM12
      if (NCOL(dfun2) < use.ncol)
        dfun2 <- cbind(dfun2, matrix(0, n.LM, use.ncol - NCOL(dfun2)))

      d2wz.dtheta2 <- if (link1parameter &&
                          M1 == 1 &&
                          length(wetas.kay) == 1) {
        dfun2 * (Der1[, jay])^4 +
        (dwz.dtheta.Der1 / Der1[, jay]) * Der2[, jay] +
        4 * dfun1 * Der2[, jay] * (Der1[, jay])^2 +
        2 * dfun0 * Der3[, jay] * Der1[, jay]
      } else {
        NA * dfun2
      }
      d2wz.dbetakk2 <- d2wz.dtheta2 * bix.jk^2
      if (!is.matrix(d2wz.dbetakk2))
        d2wz.dbetakk2 <- as.matrix(d2wz.dbetakk2)

    d2wz.dbetas2 <- d2wz.dbetas2 + d2wz.dbetakk2  # Summed over 1:M
    } # !fd.use
    }  # (type == "derivatives" && derivative == 2)
  }  # for (jay in wetas.kay) # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,








    tmp.mat <- mux111(t(dwz.dbetas), X.vlm, M = M, upper = FALSE)


    dA.dbeta <- crossprod(X.vlm, tmp.mat)  # p.VLM x p.VLM




    SE2.kay <- SE2[kay]





    small.temp1 <- dA.dbeta %*% vc2[, kay, drop = FALSE]
    small.d1ixwx.dbeta1 <- -(vc2[kay, , drop = FALSE] %*% small.temp1)
    SE2.deriv1[kay] <- small.d1ixwx.dbeta1




    vec.Wald.deriv1[kay] <- (1 - 0.5 * cobj[kay] *
                            SE2.deriv1[kay] / SE2.kay) / SE1[kay]



    if (type == "derivatives" && derivative == 2) {
      tmp.mat <- mux111(t(d2wz.dbetas2), X.vlm, M = M, upper = FALSE)
      d2A.dbeta2 <- crossprod(X.vlm, tmp.mat)  # p.VLM x p.VLM


    temp1 <- dA.dbeta %*% vc2



    small.d2ixwx.dbeta2 <-
        vc2[kay, , drop = FALSE] %*%
        (2 * temp1 %*% dA.dbeta - d2A.dbeta2) %*%
        vc2[, kay, drop = FALSE]
      SE2.deriv2[kay] <- if (blot.out) NA else small.d2ixwx.dbeta2



      vec.Wald.deriv2[kay] <- if (blot.out) NA else
        (-SE2.deriv1[kay] +
         0.5 * cobj[kay] *
         (1.5 * ((SE2.deriv1[kay])^2) / SE2.kay -
                  SE2.deriv2[kay])) / (SE2.kay^1.5)
    }  # derivative == 2
  }  # for (kay in kvec.use)  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,



  SE.deriv1 <- if (se.arg) 0.5 * SE2.deriv1 / SE1 else NULL
  ans <-
  switch(type,
         logical = vec.Wald.deriv1 < 0,  # yettodo: 2nd-deriv test later
         derivatives =
           if (derivative == 1) {
             if (se.arg)
               cbind(deriv1    = vec.Wald.deriv1,
                     SE.deriv1 = SE.deriv1) else
               vec.Wald.deriv1
           } else {
             cbind(deriv1 = vec.Wald.deriv1,
                   deriv2 = vec.Wald.deriv2,
                   SE.deriv1 = if (se.arg) SE.deriv1 else NULL,
                   SE.deriv2 = if (se.arg) 0.5 * (SE2.deriv2 / SE1
                               - 0.5 * SE2.deriv1^2  / SE2^1.5) else
                               NULL)
           })
  if (length(subset))
    ans <- if (is.matrix(ans))
      ans[kvec.use, , drop = FALSE] else
      ans[kvec.use]
  ans
}  # hdeff.vglm





 if (!isGeneric("hdeff"))
    setGeneric("hdeff", function(object, ...) standardGeneric("hdeff"))


setMethod("hdeff", "vglm", function(object, ...)
          hdeff.vglm(object, ...))







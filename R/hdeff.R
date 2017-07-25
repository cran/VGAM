# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.














hdeff.vglm <-
  function(object,
           derivative = NULL,
           se.arg = FALSE,
           ...) {


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
  dfun <-
    if (is.logical(Fam.infos$hadof) && Fam.infos$hadof) {
      Fam@hadof
    } else {
      return(NULL)  # This means its not implemented yet
      stop(gettextf(
      "hdeff() not implemented for family '%s' yet", Fam), domain = NA)
    }


  link1parameter <- Fam.infos$link1parameter
  if (is.null(link1parameter))
    link1parameter <- TRUE  # The default really




  M <- npred(object)  # Some constraints span across responses
  all.Hk <- constraints(object, matrix = TRUE)
  X.vlm <- model.matrix(object, type = "vlm")
  eta.mat <- predict(object)
  nnn <- NROW(eta.mat)
  pwts <- weights(object, type = "prior")

  mylinks <- linkfun(object)  # Of length 1 for GLMs, char only
  wwts <- weights(object, type = "working")
  dim.wz <- dim(wwts)  # Inefficient

  p.VLM <- ncol(all.Hk)
  M1 <- npred(object, type = "one.response")

  vc2 <- vcov(object)
  SE2 <- diag.ixwx <- diag(vc2)
  SE1 <- sqrt(SE2)
  cobj <- coef(object)
  se2.deriv1 <- vec.deriv1 <- rep_len(NA_real_, p.VLM)
  names(vec.deriv1) <- names(cobj)
  if (type == "derivatives" && derivative == 2) {
    se2.deriv2 <- vec.deriv2 <- vec.deriv1
  }



  D3thetas.Detas3 <-  # May not be needed
  D2thetas.Detas2 <-
  D1thetas.Detas1 <-
  Param.mat       <- matrix(NA_real_, nnn, M)
  if (link1parameter) {
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
    ind5 <- iam(NA, NA, both = TRUE, M = M)
    wz.tet <-  D1thetas.Detas1[, ind5$row] *
               D1thetas.Detas1[, ind5$col]  # n x MM12
      Der1 <- D1thetas.Detas1
      Der2 <- D2thetas.Detas2
  } else {
    MM12 <- M * (M + 1) / 2
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



  for (kay in 1:p.VLM) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  dwz.dbetas   <-
  d2wz.dbetas2 <- 0  # Good for the first instance of use.


  for (jay in 1:M) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  if (all.Hk[jay, kay] != 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    vecTF.jay <- as.logical(eijfun(jay, M))
    bix.jk <- X.vlm[vecTF.jay, kay]  # An n-vector





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
          dfun1  <- cbind(dfun1, matrix(0, nnn, use.ncol - NCOL(dfun1)))
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



    if (type == "derivatives" && derivative == 2) {
      if (link1parameter) {
        Der3 <- D3thetas.Detas3
      }
      dfun2 <- dfun(eta.mat, extra = object@extra,
                    linpred.index = jay,
                    w = pwts, dim.wz = dim.wz,
                    deriv = 2)
      use.ncol <- if (link1parameter) NCOL(dwz.dtheta.Der1) else MM12
      if (NCOL(dfun2) < use.ncol)
        dfun2 <- cbind(dfun2, matrix(0, nnn, use.ncol - NCOL(dfun2)))

      d2wz.dtheta2 <- if (link1parameter &&
                          M1 == 1) {
        dfun2 * (Der1[, jay])^4 + (dwz.dtheta.Der1 / Der1[, jay]) * Der2[, jay] +
        4 * dfun1 * Der2[, jay] * (Der1[, jay])^2 +
        2 * dfun0 * Der3[, jay] * Der1[, jay]
      } else {
        NA * dfun2
      }
      d2wz.dbetakk2 <- d2wz.dtheta2 * bix.jk^2
      if (!is.matrix(d2wz.dbetakk2))
        d2wz.dbetakk2 <- as.matrix(d2wz.dbetakk2)

    d2wz.dbetas2 <- d2wz.dbetas2 + d2wz.dbetakk2  # Summed over 1:M
    }  # (type == "derivatives" && derivative == 2)
  }  # if (all.Hk[jay, kay] != 0)
  }  # for (jay in 1:M) # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,








    tmp.mat <- mux111(t(dwz.dbetas), X.vlm, M = M, upper = FALSE)


    dA.dbeta <- crossprod(X.vlm, tmp.mat)  # p.VLM x p.VLM

    temp1 <- dA.dbeta %*% vc2
    d1ixwx.dbeta1 <- -(vc2 %*% temp1)
    SE2.kay <- SE2[kay]
    se2.deriv1[kay] <- diag(d1ixwx.dbeta1)[kay]
    vec.deriv1[kay] <- (1 - 0.5 * cobj[kay] *
                       se2.deriv1[kay] / SE2.kay) / SE1[kay]



    if (type == "derivatives" && derivative == 2) {
      tmp.mat <- mux111(t(d2wz.dbetas2), X.vlm, M = M, upper = FALSE)
      d2A.dbeta2 <- crossprod(X.vlm, tmp.mat)  # p.VLM x p.VLM
      d2ixwx.dbeta2 <-
        vc2 %*% (2 * temp1 %*% dA.dbeta - d2A.dbeta2) %*% vc2
      se2.deriv2[kay] <- diag(d2ixwx.dbeta2)[kay]
      vec.deriv2[kay] <-
        (-se2.deriv1[kay] +
         0.5 * cobj[kay] *
         (1.5 * ((se2.deriv1[kay])^2) / SE2.kay -
                  se2.deriv2[kay])) / (SE2.kay^1.5)
    }  # derivative == 2
  }  # for (kay in 1:p.VLM)  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,


  se.deriv1 <- if (se.arg) 0.5 * se2.deriv1 / SE1 else NULL
  switch(type,
         logical = vec.deriv1 < 0,  # yettodo: 2nd-deriv test later
         derivatives =
           if (derivative == 1) {
             if (se.arg)
               cbind(deriv1    = vec.deriv1,
                     se.deriv1 = se.deriv1) else
               vec.deriv1
           } else {
             cbind(deriv1 = vec.deriv1,
                   deriv2 = vec.deriv2,
                   se.deriv1 = if (se.arg) se.deriv1 else NULL,
                   se.deriv2 = if (se.arg) 0.5 * (se2.deriv2 / SE1
                               - 0.5 * se2.deriv1^2  / SE2^1.5) else
                               NULL)
           })
}  # hdeff.vglm





 if (!isGeneric("hdeff"))
    setGeneric("hdeff", function(object, ...) standardGeneric("hdeff"))


setMethod("hdeff", "vglm", function(object, ...)
          hdeff.vglm(object, ...))







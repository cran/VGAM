# These functions are
# Copyright (C) 1998-2021 T.W. Yee, University of Auckland.
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
    link1parameter <- TRUE  # The default, for ordinary 1-par links





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
    mixture.links <- is.logical(Fam.infos$mixture.links) &&
                     Fam.infos$mixture.links
    Param.mat <- eta2theta(eta.mat, mylinks,
                           earg = object@misc$earg,
                           delete.coln = !mixture.links)








    if (mixture.links) {



if (FALSE) {  ###################################################
      myrle <- rle(mylinks)
      if (length(myrle$value) != 2)
        stop("can only handle two types of links in two chunks")
      M.1 <- myrle$length[1]
      myearg1 <- object@misc$earg[[1]]
      build.list1 <- list(theta = Param.mat[, 1:(1 + M.1)],
                          inverse = TRUE, deriv = 1)
      build.list1 <- c(build.list1, myearg1)  # No dups arg names..
      build.list1$all.derivs <- TRUE  # For "multilogitlink".
      Der11 <- do.call(mylinks[1], build.list1)
      M.2 <- 1  # Corresponding to, e.g., loglink("lambda")
      lastone <- length(object@misc$earg)
      tmp5 <- ncol(Param.mat)
      myearg2 <- object@misc$earg[[lastone]]
      myearg2$theta <- Param.mat[, (2 + M.1):tmp5]
      myearg2$inverse <- TRUE
      myearg2$deriv <- 1
      Der12 <- do.call(mylinks[lastone], myearg2)
      Der1 <- wz.merge(Der11, Der12, M.1, M.2)  # Combine them
}  # if (FALSE)  ###################################################




  llink <- length(mylinks)  # length(link)
  vecTF <- mylinks == "multilogitlink"
  Ans <- NULL  # Growing matrix data structure
  M.1 <- 0
  offset.Param.mat <- 0  # Coz rowSums(Param.mat)==1 for multilogitlink
  iii <- 1
  while (iii <= llink) {
    first.index <- last.index <- iii  # Ordinary case
    special.case <- vecTF[iii]  # && sum(vecTF) < length(vecTF)

    if (special.case) {
      next.i <- iii+1
      while (next.i <= llink) {
        if (vecTF[next.i]) {
          last.index <- next.i
          next.i <- next.i + 1
        } else {
          break
        }
      }  # while
    }  # special.case

    iii <- iii + last.index - first.index + 1  # For next time

    myearg2 <- object@misc$earg[[first.index]]  # Only one will do
    if (special.case) {
      build.list2 <-  # rowSums of Param.mat subset are all == 1:
        list(theta = Param.mat[, first.index:(last.index+1)],
             inverse = TRUE, deriv = 1)
      offset.Param.mat <- offset.Param.mat + 1
      myearg2 <- c(build.list2, myearg2)  # No dups arg names..
      myearg2$all.derivs <- TRUE  # For "multilogitlink".
      myearg2$M <- last.index - first.index  + 1
    }  # special.case

    use.earg <- myearg2
    use.earg[["inverse"]] <- TRUE  # New
    if (!special.case)
      use.earg[["theta"]] <-
        Param.mat[, offset.Param.mat + (first.index:last.index)]
    use.earg$deriv <- 1
    use.function.name <- mylinks[first.index]  # e.g., "multilogitlink"

    Ans2 <- do.call(use.function.name, use.earg)
    delete.coln <- FALSE
    if (special.case && delete.coln)
      Ans2 <- Ans2[, -use.earg$refLevel]

    M.2 <- last.index - first.index + 1
    Ans <- if (length(Ans)) wz.merge(Ans, Ans2, M.1, M.2) else Ans2
    M.1 <- M.1 + M.2
  }  # while (iii <= llink)







    } else {
 # Handle multinomial, etc.
      myearg <- object@misc$earg[[1]]  # Only ONE anyway
      build.list <- list(theta = Param.mat,
                    inverse = TRUE, deriv = 1)  # This line is important
      build.list <- c(build.list, myearg)  # Hopefully no dups arg names
      build.list$all.derivs <- TRUE  # For "multilogitlink".
      Der1 <- do.call(mylinks, build.list)  # n x MM12 for multinomial
    }  # mixture.links & !mixture.links

 
    if (type == "derivatives" && derivative == 2) {
    if (mixture.links) {
      build.list1$deriv <- 2
      Der21 <- do.call(mylinks[1], build.list1)
      myearg2$deriv <- 2
      Der22 <- do.call(mylinks[lastone], myearg2)
      Der2 <- wz.merge(Der21, Der22, M.1, M.2)  # Combine them
    } else {
      build.list$deriv <- 2
      Der2 <- do.call(mylinks, build.list)  # n x M for multinomial
     } # mixture.links & !mixture.links
   }  # derivative == 2
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
      }  # write.into.wz

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
    }  # link1parameter and !link1parameter

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








hdeff.matrix <-
  function(object,
           ...) {
  if (!is.matrix(object) || nrow(object) != 2 || ncol(object) != 2) 
    stop("argument 'object' is not a 2 x 2 matrix")
  if (any(c(object) <= 0))
    stop("some cells are not positive valued")
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  if (!all(is.wholenumber(c(object))))
    stop("some cells are not integer-valued")

  N0 <- sum(object[1, ])
  N1 <- sum(object[2, ])
  p0 <- object[1, 2] / N0
  p1 <- object[2, 2] / N1
  oratio <- object[1, 1] * object[2, 2] / (object[1, 2] * object[2, 1])
  beta2 <- log(oratio)
  lhs <- 1 + N1 * p1 * (1 - p1) / (N0 * p0 * (1 - p0))
  rhs <- beta2 * (p1 - 0.5)
  1 + N1 * p1 * (1 - p1) / (N0 * p0 * (1 - p0)) < beta2 * (p1 - 0.5)
}  # hdeff.matrix




hdeff.numeric <-
  function(object,
           byrow = FALSE,
           ...) {
  if (length(c(object)) > 4)
    stop("length of argument 'object' greater than 4")
  if (length(c(object)) < 4)
    warning("length of argument 'object' less than 4, so recycling")
  hdeff(matrix(c(object), 2, 2, byrow = byrow))  # Recycles if needed
}  # hdeff.numeric







 if (!isGeneric("hdeff"))
    setGeneric("hdeff", function(object, ...) standardGeneric("hdeff"))


setMethod("hdeff", "vglm", function(object, ...)
          hdeff.vglm(object, ...))


setMethod("hdeff", "matrix", function(object, ...)
          hdeff.matrix(object, ...))

setMethod("hdeff", "numeric", function(object, ...)
          hdeff.numeric(object, ...))









hdeffsev <-
  function(x, y,
           dy, ddy,  # 1st and 2nd derivs           
           allofit = FALSE,
           tol0 = 0.1,
           severity.table = c("None", "Faint", "Weak",
                              "Moderate", "Strong", "Extreme",
                              "Undetermined")) { 







  severity <- rep_len(severity.table[7], length(x))  # Initialize
  names(severity) <- names(x)
  zeta <- x + y * dy
  dzeta.dx <- 1 + dy^2 + y * ddy
  ind.none     <- dy > 0 &
                  ifelse(0 <= x, ddy > 0, ddy < 0) &  # 20181105
                  dzeta.dx > 0
  severity[ind.none] <- severity.table[1]
  severity[abs(x) < tol0] <- severity.table[1]  # Additional condition
  ind.faint    <- dy > 0 &
                  ifelse(0 <= x, ddy <= 0, ddy >= 0) &
                  dzeta.dx > 0
  severity[ind.faint] <- severity.table[2]
  ind.weak     <- dy > 0 &
                  ifelse(0 <= x, ddy < 0, ddy > 0) &
                  dzeta.dx < 0
  severity[ind.weak] <- severity.table[3]
  ind.moderate <- dy <= 0 &  # Note <= rather than <
                  ifelse(0 <= x, ddy < 0, ddy > 0) & 
                  dzeta.dx < 0
  severity[ind.moderate] <- severity.table[4]
  ind.strong   <- dy < 0 &
                  ifelse(0 <= x, ddy < 0, ddy > 0) & 
                  dzeta.dx > 0
  severity[ind.strong] <- severity.table[5]
  ind.extreme  <- dy < 0 &
                  ifelse(0 <= x, ddy >= 0, ddy < 0) & 
                  dzeta.dx > 0
  severity[ind.extreme] <- severity.table[6]
  if (allofit) list(severity  = severity,
                    zeta      = zeta,
                    dzeta.dx  = dzeta.dx) else
               severity
}



seglines <-
  function(x, y,
           dy, ddy,  # 1st and 2nd derivs           
           lwd = 2,
           cex = 2,
           plot.it = TRUE,
           add.legend = TRUE,
           position.legend = "topleft",
           lty.table = c("solid", "dashed",
                         "solid", "dashed",
                         "solid", "dashed",
                         "solid"), 
           col.table = rainbow.sky[-5],  # Omit "yellow"
           pch.table = 7:1,
           severity.table = c("None", "Faint", "Weak",
                              "Moderate", "Strong", "Extreme",
                              "Undetermined"),
           tol0 = 0.1,
           FYI = FALSE,
           ...) {


  answer <- hdeffsev(x, y, dy, ddy,
                     severity.table = severity.table,
                     tol0 = tol0, allofit = FYI)
  severity <- if (FYI) answer$severity else answer

  if (plot.it) {
    myrle <- rle(severity)
    myrle$cslength <- cumsum(myrle$length)
    mycol <- col.table[match(severity, severity.table)]
    mylty <- lty.table[match(severity, severity.table)]
    mypch <- pch.table[match(severity, severity.table)]

    single.points <- FALSE  # Assumes all lines()
    pointsvec <- NULL  # History of points used
    for (iloc in seq(length(myrle$values))) {
      end.val <- myrle$cslength[iloc]
      start.val <- end.val + 1 - myrle$length[iloc]
      if (start.val < end.val) {
        lines(x[start.val:end.val],
              y[start.val:end.val],
              lwd = lwd, col = mycol[start.val:end.val],
              lty = mylty[start.val:end.val])
      } else {
        single.points <- TRUE
        pointsvec <- c(pointsvec, mypch[start.val])
        points(x[start.val:end.val],
               y[start.val:end.val],
               col = mycol[start.val:end.val],
               pch = mypch[start.val:end.val],
               cex = cex)
      }

      if (FYI) {
        some.val <- sample(start.val:end.val,
                           min(2, end.val-start.val+1))  
        segments(x[some.val],
                 y[some.val],
                 ifelse(x[some.val] > 0,
                         answer$zeta[some.val],
                        -answer$zeta[some.val]),
                 0,
                 col = "purple", lty = "dashed")
      }  # FYI
    }  # for

    pointsvec <- unique(pointsvec)
    use.pch.table <- pch.table
    for (ii in 1:7)
      if (single.points && !any(use.pch.table[ii] == pointsvec))
        use.pch.table[ii] <- NA


    if (add.legend) {
      if (FALSE && !any(is.element(severity,  severity.table[7]))) {
        use.pch.table <- use.pch.table[-7]
        col.table <- col.table[-7]
        severity.table <- severity.table[-7]
      }
      ind3 <- match(severity.table, severity)  # length of 1st argument
      keep.ind <- !is.na(ind3)
      use.pch.table <- use.pch.table[keep.ind]
      col.table <- col.table[keep.ind]
      severity.table <- severity.table[keep.ind]
      
      legend(position.legend, lwd = lwd, lty = lty.table,
             pch = use.pch.table, col = col.table,
             legend = severity.table)
      invisible(severity)
    }
  } else {
    return(severity)
  }
}  # seglines


















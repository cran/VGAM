# These functions are
# Copyright (C) 1998-2025 T.W. Yee, University of Auckland.
# All rights reserved.













hdeff.vglm <-
  function(object,
       derivative = NULL,
       se.arg = FALSE,
       subset = NULL,  # Useful for Cox model as a poissonff().
       theta0 = 0,  # Recycled to the necessary length 20210406
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
        "derivatives" else
        stop("bad input for argument 'derivative'")
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
    link1parameter <- TRUE  # (default) for ordinary 1-par links





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
  if (length(theta0) > p.VLM)
    warning("length of argument 'theta0' is loo long. ",
            "Truncating it.")
  theta0 <- rep_len(theta0, p.VLM)
  M1 <- npred(object, type = "one.response")

  vc2 <- vcov(object)
  SE2 <- diag.ixwx <- diag(vc2)
  SE1 <- sqrt(SE2)
  cobj <- coef(object) - theta0
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
    for (jay in 1:M) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
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
      build.list1 <- c(build.list1,
                       myearg1)  # No dups arg names..
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
}  # if (FALSE)  ##############################################




  llink <- length(mylinks)  # length(link)
  vecTF <- mylinks == "multilogitlink"
  Ans <- NULL  # Growing matrix data structure
  M.1 <- 0
  offset.Param.mat <- 0  # Coz
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
    use.function.name <-
      mylinks[first.index]  # e.g., "multilogitlink"

    Ans2 <- do.call(use.function.name, use.earg)
    delete.coln <- FALSE
    if (special.case && delete.coln)
      Ans2 <- Ans2[, -use.earg$refLevel]

    M.2 <- last.index - first.index + 1
    Ans <- if (length(Ans))
             wz.merge(Ans, Ans2, M.1, M.2) else Ans2
    M.1 <- M.1 + M.2
  }  # while (iii <= llink)







    } else {
 # Handle multinomial, etc.
      myearg <- object@misc$earg[[1]]  # Only ONE anyway
      build.list <- list(theta = Param.mat,
                         inverse = TRUE,
                         deriv = 1)  # This line is important
      build.list <-
        c(build.list, myearg)  # Hopefully no dups arg names
      build.list$all.derivs <- TRUE  # For "multilogitlink".
      Der1 <- do.call(mylinks,
                      build.list)  # n x MM12 4 multinomial
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
      Der2 <- do.call(mylinks,
                      build.list)  # n x M for multinomial
     } # mixture.links & !mixture.links
   }  # derivative == 2
  }  # if (link1parameter) and (!link1parameter)




  kvec.use <- 1:p.VLM
  if (length(subset))
    kvec.use <- kvec.use[subset]  #&!is.na(subset)




  for (kay in kvec.use) {  # ,,,,,,,,,,,,,,,,,,,,,
  dwz.dbetas   <-
  d2wz.dbetas2 <- 0  # Good 4the 1st instance of use.


  wetas.kay <- which.etas(object, kay = kay)
  for (jay in wetas.kay) {  # ,,,,,,,,,,,,,,,,,,,
    vecTF.jay <- as.logical(eijfun(jay, M))
    bix.jk <- X.vlm[vecTF.jay, kay]  # An n-vector for xij.



    if (fd.use) {
      for (dirr in 1:2) {
        temp1 <- object
        temp1@predictors[, jay] <- temp1@predictors[, jay] +
                                   (c(1, -1)[dirr]) * hstep[1]
        temp1@fitted.values <- cbind(
          temp1@family@linkinv(eta = temp1@predictors,
                  extra = temp1@extra))  # Make sure a matrix
        if (dirr == 1)
          wwt.f1 <- weights(temp1, type = "working",
                            ignore.slot = TRUE)
        if (dirr == 2)
          wwt.b1 <- weights(temp1, type = "working",
                            ignore.slot = TRUE)
      }  # dirr
      if (ncol(wwt.f1) < MM12)
        wwt.f1 <- cbind(wwt.f1,
                        matrix(0, n.LM, MM12 - ncol(wwt.f1)))
      if (ncol(wwt.b1) < MM12)
        wwt.b1 <- cbind(wwt.b1,
                        matrix(0, n.LM, MM12 - ncol(wwt.b1)))


      cdiff1 <- (wwt.f1 - wwt.b1) / (2 * hstep[1])
      cdiff2 <- (wwt.f1 - 2 * wwt.0 + wwt.b1) / (hstep[1])^2

      dwz.dbetas   <- dwz.dbetas   + as.matrix(cdiff1) * bix.jk
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
          dfun1  <- cbind(dfun1,
                          matrix(0, n.LM, use.ncol - NCOL(dfun1)))
        if (use.ncol < NCOL(wz.tet))
          wz.tet <- wz.tet[, 1:use.ncol, drop = FALSE]

      write.into.wz <- function(jay, nxM) {
        M <- NCOL(nxM)
        wz <- matrix(0, NROW(nxM), M*(M+1)/2)
        for (uuu in 1:M)
          wz[, iam(jay, uuu, M = M)] <-
            (1 + (jay == uuu)) * nxM[, uuu]
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
  for (ttt in wetas.kay) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
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
        wwt.1 <- weights(temp1, type = "working",
                         ignore.slot = TRUE)
      if (dirr == 2)
        wwt.2 <- weights(temp1, type = "working",
                         ignore.slot = TRUE)
      if (dirr == 3)
        wwt.3 <- weights(temp1, type = "working",
                         ignore.slot = TRUE)
      if (dirr == 4)
        wwt.4 <- weights(temp1, type = "working",
                         ignore.slot = TRUE)
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
      use.ncol <- if (link1parameter)
                    NCOL(dwz.dtheta.Der1) else MM12
      if (NCOL(dfun2) < use.ncol)
        dfun2 <- cbind(dfun2,
                       matrix(0, n.LM, use.ncol - NCOL(dfun2)))

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
  }  # for (jay in wetas.kay) # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,








    tmp.mat <- mux111(t(dwz.dbetas), X.vlm, M = M, upper = FALSE)


    dA.dbeta <- crossprod(X.vlm, tmp.mat)  # p.VLM x p.VLM




    SE2.kay <- SE2[kay]





    small.temp1 <- dA.dbeta %*% vc2[, kay, drop = FALSE]
    small.d1ixwx.dbeta1 <- -(vc2[kay, , drop = FALSE] %*%
                             small.temp1)
    SE2.deriv1[kay] <- small.d1ixwx.dbeta1




    vec.Wald.deriv1[kay] <- (1 - 0.5 * cobj[kay] *
                  SE2.deriv1[kay] / SE2.kay) / SE1[kay]



    if (type == "derivatives" && derivative == 2) {
      tmp.mat <- mux111(t(d2wz.dbetas2), X.vlm,
                        M = M, upper = FALSE)
      d2A.dbeta2 <- crossprod(X.vlm, tmp.mat)  # p.VLM x p.VLM


    temp1 <- dA.dbeta %*% vc2



    small.d2ixwx.dbeta2 <-
        vc2[kay, , drop = FALSE] %*%
        (2 * temp1 %*% dA.dbeta - d2A.dbeta2) %*%
        vc2[, kay, drop = FALSE]
        SE2.deriv2[kay] <- if (blot.out) NA else
                           small.d2ixwx.dbeta2



      vec.Wald.deriv2[kay] <- if (blot.out) NA else
        (-SE2.deriv1[kay] +
         0.5 * cobj[kay] *
         (1.5 * ((SE2.deriv1[kay])^2) / SE2.kay -
                  SE2.deriv2[kay])) / (SE2.kay^1.5)
    }  # derivative == 2
  }  # for (kay in kvec.use)  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,



  SE.deriv1 <- if (se.arg) 0.5 * SE2.deriv1 / SE1 else NULL
  ans <-
  switch(type,
         logical = vec.Wald.deriv1 < 0,  # yettodo
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


  if (any(is.na(ans)) && !fd.only) {
    warning("NAs detected. Setting 'fd.only = TRUE' and ",
            "making a full recursive call")
    ans <- hdeff.vglm(object, derivative = derivative,
                      se.arg = se.arg, subset = subset,
                      theta0 = theta0, hstep = hstep,
                      fd.only = TRUE, ...)
  }



  ans
}  # hdeff.vglm








hdeff.matrix <-
  function(object,
           ...) {
  if (!is.matrix(object) || nrow(object) != 2 ||
      ncol(object) != 2) 
    stop("argument 'object' is not a 2 x 2 matrix")
  if (any(c(object) <= 0))
    stop("some cells are not positive valued")
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)
      abs(x - round(x)) < tol
  if (!all(is.wholenumber(c(object))))
    stop("some cells are not integer-valued")

  N0 <- sum(object[1, ])
  N1 <- sum(object[2, ])
  p0 <- object[1, 2] / N0
  p1 <- object[2, 2] / N1
  oratio <- object[1, 1] *
            object[2, 2] / (object[1, 2] * object[2, 1])
  beta2 <- log(oratio)
  lhs <- 1 + N1 * p1 * (1 - p1) / (N0 * p0 * (1 - p0))
  rhs <- beta2 * (p1 - 0.5)
  1 + N1 * p1 * (1 - p1) / (N0 * p0 * (1 - p0)) <
  beta2 * (p1 - 0.5)
}  # hdeff.matrix




hdeff.numeric <-
  function(object,
           byrow = FALSE,
           ...) {
  if (length(c(object)) > 4)
    stop("length of argument 'object' greater than 4")
  if (length(c(object)) < 4)
    warning("length of argument 'object' < 4, so recycling")
  hdeff(matrix(c(object), 2, 2,
               byrow = byrow))  # Recycles if needed
}  # hdeff.numeric







if (!isGeneric("hdeff"))
  setGeneric("hdeff", function(object, ...)
    standardGeneric("hdeff"))


setMethod("hdeff", "vglm", function(object, ...)
          hdeff.vglm(object, ...))


setMethod("hdeff", "matrix", function(object, ...)
          hdeff.matrix(object, ...))

setMethod("hdeff", "numeric", function(object, ...)
          hdeff.numeric(object, ...))









hdeffsev0 <-
  function(x, y,
           dy, ddy,  # 1st and 2nd derivs
           allofit = FALSE,
           eta0 = 0,  # NA or NULL means dont know, may be Inf
           COPS0 = eta0,  # Assumption. May be Inf.
           severity.table =    # if (ndepends)
                            c("None",
                              "Faint",  # Retaining this
                              "Weak",
                              "Moderate",
                              "Strong",
                              "Extreme",  # ==Extreme1--
                              "Undetermined")   #  else
          ) { 



  if (is.unsorted(x))
    stop("argument 'x' must be sorted")

  if (is.na(eta0) || is.null(eta0)) {
    splinethrudata1 <- function(X)
      spline(x, y, xout = X)$y
    eta0 <- uniroot(splinethrudata1,
                    interval = range(x))$root
  }

  LLL <- length(severity.table)
  severity <- rep_len(severity.table[LLL], length(x))
  names(severity) <- names(x)

  zeta <- x + y * dy
  dzeta.dx <- 1 + dy^2 + y * ddy
  tanzeta <- x - y / dy
  dtanzeta.dx <- y * ddy / dy^2







  ind.none <- dy >= 0   # &  # Not SSD: 20220829
  severity[ind.none] <- severity.table[1]



w10.n <- w10.p <- NULL
ind.w10.n <- min(which(x < COPS0 & dy >= 0))
w10.n <- x[ind.w10.n]
ind.w10.p <- max(which(x > COPS0 & dy >= 0))
w10.p <- x[ind.w10.p]

w20.n <- w20.p <- NULL
ind.w20.n <- min(which(x < COPS0 & ddy >= 0))
w20.n <- x[ind.w20.n]
ind.w20.p <- min(which(x > COPS0 & ddy >= 0))
w20.p <- x[ind.w20.p]


if (is.infinite(COPS0)) {
  w10.p <- (-w10.n)  # Used for tanline
  w20.p <- (-w20.n)  # Used for tanline
}
if (is.infinite(COPS0) && COPS0 < 0)
  stop("cannot handle -Inf just now")








  ind.faint    <- dy >= 0 &
      ifelse(COPS0 <= x,
             w20.n <= tanzeta & tanzeta <= w10.n,
             w10.p <= tanzeta & tanzeta <= w20.p)
  severity[ind.faint] <- severity.table[2]





  ind.weak     <- dy >= 0 &
                  ifelse(COPS0 <= x,
                         tanzeta <= w20.n,
                         tanzeta >= w20.p)
  severity[ind.weak] <- severity.table[3]







  ind.moderate <- dy <= 0 &  # Note <= rather than <
                  ifelse(COPS0 <= x,
                         w10.p <= x & x <= w20.p,
                         w20.n <= x & x <= w10.n)
  severity[ind.moderate] <- severity.table[4]






  ind.strong   <- dy <= 0 &  # Note <= rather than <
                  ifelse(COPS0 <= x,
                         w20.p <= x,
                         x <= w20.n)
  severity[ind.strong] <- severity.table[5]








  if (any(ind.strong, na.rm = TRUE)) {




    w20.n.next <- tanzeta[ind.w20.n]
    w20.p.next <- tanzeta[ind.w20.p]
    ind.extreme  <- dy <= 0 &  # Note <= rather than <
                    ifelse(COPS0 <= x,
                           w20.p.next <= x,
                           x <= w20.n.next)
    severity[ind.extreme] <- severity.table[6]
  }  # Extreme done here.



  if (FALSE && !is.na(w20.n)) {
    w20.n.next <- 1  # zz
    w20.n.next <- 1  # zz

  }

  if (allofit)
    list(severity     = severity,
         zeta         = zeta,
         dzeta.dx     = dzeta.dx,
         x            = x,
         y            = y,
         tanzeta      = tanzeta,
         dtanzeta.dx  = dtanzeta.dx) else
    severity
}  # hdeffsev0








hdeffsev2 <-
  function(x, y,
           dy, ddy,  # 1st and 2nd derivs           
           allofit = FALSE,
           ndepends = FALSE,  # 20240703
           eta0 = 0,  # NA or NULL means dont know, may be Inf
           severity.table = c("None", "Faint", "Weak",
             "Moderate", "Strong",
             "Extreme",
             "Undetermined")[if (ndepends) TRUE else
             c(1, 4, 6, 7)],
           tol0 = 0.1) { 

  if ((Lx <- length(x)) != length(y) ||
       Lx != length(dy) ||  Lx != length(ddy))
    stop("Args 'x', 'y', 'dy' and 'ddy' ",
         "must have equal lengths")

  LSE <- length(severity.table)
  severity <- rep(severity.table[LSE], Lx)
  names(severity) <- names(x)
  zeta <- x + y * dy  # Normal line
  dzeta.dx <- 1 + dy^2 + y * ddy
  tanzeta <- x - y / dy
  dtanzeta.dx <- y * ddy / dy^2

  if (is.na(eta0) || is.null(eta0)) {
    splinethrudata1 <- function(X)
      spline(x, y, xout = X)$y
    eta0 <- uniroot(splinethrudata1,
                    interval = range(x))$root
  }

  if (ndepends) {
  warning("20240703; answer is sample size dependent!")
  fullans <- hdeffsev2(x, y, dy, ddy,
           allofit = FALSE,  # TRUE,
           ndepends = FALSE,  # Sparse answer
           eta0 = eta0,
           tol0 = tol0)
  return(fullans)
  ind.none <-
      dy > 0 &
      ifelse(eta0 <= x, ddy > 0, ddy < 0) &  # 20181105
      dzeta.dx > 0
  severity[ind.none] <- severity.table[1]
  severity[abs(x) < tol0] <- severity.table[1]  # Additional cond.
  ind.faint <-
      dy > 0 &
      ifelse(eta0 <= x, ddy <= 0, ddy >= 0) &
      dzeta.dx > 0
  severity[ind.faint] <- severity.table[2]
  ind.weak <-
      dy > 0 &
      ifelse(eta0 <= x, ddy < 0, ddy > 0) &
      dzeta.dx < 0
  severity[ind.weak] <- severity.table[3]
  ind.moderate <-
      dy <= 0 &  # Note <= rather than <
      ifelse(0 <= x, ddy < 0, ddy > 0) & 
      dzeta.dx < 0
  severity[ind.moderate] <- severity.table[4]
  ind.strong <-
      dy < 0 &
      ifelse(0 <= x, ddy < 0, ddy > 0) & 
      dzeta.dx > 0
  severity[ind.strong] <- severity.table[5]
  ind.extreme <-
      dy < 0 &
      ifelse(0 <= x, ddy >= 0, ddy < 0) & 
      dzeta.dx > 0
  severity[ind.extreme] <- severity.table[6]
  } else {  # ------------------------------
  vecTF.xy <- is.finite(x) & is.finite(dy) &
              is.finite(y) & is.finite(ddy)
  ind.none <- vecTF.xy & dy > 0 
  severity[ind.none] <- severity.table[1]
  ind.moderate <-
    vecTF.xy & dy <= 0 &  # <=, not <
    ifelse(eta0 <= x, ddy <= 0, ddy >= 0)
  severity[ind.moderate] <- severity.table[2]
  ind.extreme  <-
    vecTF.xy & dy <= 0 &
    ifelse(eta0 <= x, ddy >= 0, ddy <= 0)
  severity[ind.extreme] <- severity.table[3]
  }

  if (allofit)
    list(severity     = severity,
         zeta         = zeta,
         dzeta.dx     = dzeta.dx,
         x            = x,
         y            = y,
         tanzeta      = tanzeta,
         dtanzeta.dx  = dtanzeta.dx) else
    severity
}  # hdeffsev2 (was the original)







seglines <-
  function(x, y,
           dy, ddy,  # 1st and 2nd derivs
           lwd = 2,
           cex = 2,
           plot.it = TRUE,
           add.legend = TRUE,
           cex.legend = 1,   # par()$cex,
           position.legend = "topleft",
           eta0 = NA,
           COPS0 = NA,  # Using eta0 instead
           lty.table = c("solid", "dashed",
                         "solid", "dashed",
                         "solid", "dashed",
                         "solid"),
           col.table = rainbow.sky[-5],  # -yellow
           pch.table = 7:1,   # 7:1,
           severity.table = c("None",
                              "Faint",
                              "Weak",
                              "Moderate",
                              "Strong",
                              "Extreme",
                              "Undetermined"),
           FYI = FALSE,
           ...) {




  Six <- 7


  answer <- hdeffsev0(x, y, dy, ddy,
                      severity.table = severity.table,
                      eta0 = eta0,
                      COPS0 = COPS0,  # Using eta0 instead
                      allofit = FYI)
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
    for (ii in 1:Six)
      if (single.points && !any(use.pch.table[ii] == pointsvec))
        use.pch.table[ii] <- NA


    if (add.legend) {
      if (TRUE &&
          !any(is.element(severity,
                          severity.table[Six]))) {
        use.pch.table <- use.pch.table[-Six]
        col.table <- col.table[-Six]
        severity.table <- severity.table[-Six]
      }
      ind3 <- match(severity.table, severity)  # len. of 1st arg
      keep.ind <- !is.na(ind3)
      use.pch.table <- use.pch.table[keep.ind]
      col.table <- col.table[keep.ind]
      severity.table <- severity.table[keep.ind]

      legend(position.legend, lwd = lwd, lty = lty.table,
             pch = use.pch.table,
             col = col.table,
             cex = cex.legend,  # Overall size
             legend = severity.table)
      invisible(severity)
    }  # add.legend
  } else {
    return(severity)
  }
}  # seglines






copsvglm <-
  function(object,
           level = 0.999,   #1 - 1e-6,
           beta.range = confint(object, level = level),
           muxrange = 4,  # 2, or 4, or 7, or...
           iter.max = 8,  # Gridsearch  # Not 20,
           tol = 1e-5,  # 4 optim()
           subset = NULL,  # 20250219
           do1 = TRUE,  # F 4 propodds()
           slowtrain = FALSE,  # TRUE,
           ...) {
  T <- TRUE; F <- FALSE
  if (is.null(subset)) {
    subset <- TRUE
  }

  if (length(muxrange) != 1 || muxrange < 1.1)
    stop("'muxrange' must be > 1.1")
  if (!is.Numeric(iter.max, positive = TRUE,
       length.arg = 1, integer.valued = TRUE) ||
      iter.max <= 2)
    stop("'iter.max' must be an integer > 2")
  if (!isFALSE(do1) && !isTRUE(do1))
    stop("'do1' must be TRUE or FALSE")


  M <- npred(object)
  cobj <- coef(object)  # Original coeffs
  pwts <- c(weights(object, type = "prior"))
  nmcobj <- names(cobj)
  p.VLM <- length(cobj)
  ind.ok <- if (is.character(subset))
    match(subset, nmcobj) else  # Unsorted
    (1:p.VLM)[subset]  # Sorted
  Hlist <- constraints(object)
  has.intercept <- names(Hlist[1]) == "(Intercept)"
  if (has.intercept) {
    if (do1) {
      ind.ok <- unique(c(1:ncol(Hlist[[1]]), ind.ok))
    } else {
      ind.ok <- setdiff(ind.ok, 1:ncol(Hlist[[1]]))
      if (!length(ind.ok))
        stop("No coefficients! Change 'subset' or set ",
             "'do1 = TRUE'")
    }
  }  # has.intercept



  if (!is(object, "vglm")) stop("not a 'vglm' object")
  nn <- nobs(object)
  Xvlm  <- model.matrix(object, type = "vlm")
  etamat <- predict(object)  # Contains offsets


  newinfo <-  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    function(beta.try1, jkay)  {
    copy.object <- object 
    copy.object@coefficients <- cobj
    copy.object@coefficients[jkay] <- beta.try1
    newetamat <- matrix(Xvlm %*%
                        copy.object@coefficients,
                        nn, M, byrow = TRUE)
    copy.object@predictors <- newetamat
    newmu <- object@family@linkinv(newetamat,
                             extra = object@extra)
    copy.object@fitted.values <- newmu
    newwz <- weights(copy.object, type = "working")
    UU <- vchol(newwz, M = M, n = nn)  # silent = T
    UtXvlm <- mux111(cc = UU, xmat = Xvlm, M = M,
                     slowtrain = slowtrain,
                     whichj = jkay)
    total.info <- sum((UtXvlm[, jkay])^2)
    total.info
  }  # newinfo


  gridpts <- 1 + 2 * iter.max
  gmat <- matrix(NA_real_, p.VLM, gridpts)
  tmp <- 1:iter.max
  nmgmat <- as.character(sort(c(-tmp, 0, tmp)))
  dimnames(gmat) <- list(nmcobj, nmgmat)
  jmat <- gmat  # To store total info
  midpts <- gmat[, '0'] <- rowMeans(beta.range)
  gmat[,  '1'] <- beta.range[, 2]
  gmat[, '-1'] <- beta.range[, 1]
  if (iter.max > 1)  # Redundant
    for (iter in 2:iter.max) {
      muxr <- muxrange * iter  # Linear spread
      acit <- as.character( iter)
      Acit <- as.character(-iter)
      gmat[, acit] <- midpts + muxr *
                     (beta.range[, 2] - midpts)
      gmat[, Acit] <- midpts - muxr *
                     (midpts - beta.range[, 1])
    }  # for iter
    for (rr in ind.ok)   #1:nrow(gmat)
      for (cc in 1:ncol(gmat))
        jmat[rr, cc] <- newinfo(gmat[rr, cc], rr)

  initLHS <- initRHS <- numeric(p.VLM)
  for (rr in ind.ok) {  # 1:p.VLM
    wmax <- which.max(jmat[rr, ])
    if (wmax %in% c(1, ncol(jmat)))
      stop("try increasing 'iter.max' and/or 'muxrange'")
    initLHS[rr] <- gmat[rr, wmax - 1]  # Bracket the
    initRHS[rr] <- gmat[rr, wmax + 1]  # maximum
  }  # rr

  copsvec <- cobj  # Named storage
  for (rr in ind.ok) {  # 1:p.VLM
    ofit <- optimize(newinfo, maximum = TRUE,
                     lower = initLHS[rr],
                     upper = initRHS[rr],
                     tol = tol, jkay = rr)
    copsvec[rr] <- ofit$maximum
  }  # rr

  copsvec[ind.ok]
}  # copsvglm








copsvglmb <-
  function(object,
           level = 0.999,   #1 - 1e-6,
           beta.range = confint(object, level = level),
           muxrange = 4,  # 2, or 4, or 7, or...
           iter.max = 20,
           tol = 1e-5,  # 4 bisection.basic()
           nmax = NULL,  # 4 bisection.basic()
           subset = NULL,  # 20250219
           ...) {
  T <- TRUE; F <- FALSE
  if (is.null(subset)) subset <- TRUE

  M <- npred(object)
  cobj <- coef(object)  # Original coeffs
  pwts <- c(weights(object, type = "prior"))

  vfamily <- object@family@vfamily
  Fam <- if (inherits(object, "vlm")) object@family else
    stop("cannot get at the 'family' slot")
  Fam.infos <- Fam@infos()
  if (!is(object, "vglm")) stop("not a 'vglm' object")
  if (!any(vfamily == "binomialff") || M > 1 ||
      !any(linkfun(object) == "logitlink"))
    stop("not a (single) logistic regression")

  nn <- nobs(object)
  Xvlm  <- model.matrix(object, type = "vlm")
  etamat <- predict(object)  # Contains offsets
  if (F) {
  Etamat <- matrix(Xvlm %*% cobj, nn, M, byrow = T)
  offset <- if (length(object@offset))
              object@offset else matrix(0, 1, 1)
  if (any(offset != 0))
    Etamat <- Etamat + offset
 }



  dinfo.dbetak <-  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    function(beta2)  {  #, jkay, jay=1, M=1
    jay2 <- numeric(length(beta2))
    for (uu in 1:length(beta2)) {
      eta <- etamat - cobj[uu] * Xvlm[, uu]
      eta <- eta + beta2[uu] * Xvlm[, uu]
      prob <- logitlink(eta, inverse = TRUE)
      wz <- pwts * prob * (1 - prob) * (1 - 2 * prob)
      jay2[uu] <- sum((Xvlm[, uu])^3 * wz)
    }
    jay2
  }  # dinfo.dbetak 







    try.int <- beta.range  # Fixed seems bad 
    iter <- 0; muxr <- muxrange
    midpts <- rowMeans(beta.range)
    Ping <- bisection.basic(dinfo.dbetak, Ping = T,
                            a = try.int[, 1],
                            b = try.int[, 2])
    while(any(!Ping)) {  # Need to bracket the root
      iter <- iter + 1
      if (iter > iter.max + 1)
        stop("could not bracket the root within ",
             iter,  " iterations")
      muxr <- muxrange * iter  # Linear spread
      try.int[!Ping, 1] <- midpts[!Ping] -
        (midpts[!Ping] - beta.range[!Ping, 1]) * muxr
      try.int[!Ping, 2] <- midpts[!Ping] +
        (beta.range[!Ping, 2] - midpts[!Ping]) * muxr
      Ping <- bisection.basic(dinfo.dbetak, Ping = T,
                              a = try.int[, 1],
                              b = try.int[, 2])
    }  # while


        ofit <- bisection.basic(dinfo.dbetak,
                         a = try.int[, 1],
                         b = try.int[, 2],
                         tol = tol,  # jay = jay,
                         nmax = nmax)  # M = M
        copsvec <- ofit

    copsvec[subset]
}  # copsvglmb







if (!isGeneric("cops"))
  setGeneric("cops",
             function(object, ...)
                 standardGeneric("cops"),
             package = "VGAM")

setMethod("cops", "vglm",
          function(object, ...) {
  copsvglm(object, ...)
})




DDfun <- function(expr, name, order = 0) {
  if (order < 0) stop("'order' must be >= 0")
  if (order == 0) return(expr)
  if (order == 1) D(expr, name) else
  DDfun(D(expr, name), name, order - 1)
}





hdeffsev <-
  function(object,
           hdiff = 0.005,
           eta0 = 0,
           subset = NULL,
           maxderiv = 6,
           severity.table = c("None", "Faint",
             "Weak", "Moderate", "Strong",
             "ExtremeI",
             "ExtremeII",
             "ExtremeIII",
             "ExtremeIV+",
             "Undetermined"),  # Needed, last
           lookup = c(0, 0.5, 0.7, 1, 1.3, 2:5),
           tx.some = TRUE,  # log, cauchit
           wsdmvec = NULL,  # Can input here
           ...) {
  lkv <- lookup
  stab <- severity.table
  if (length(lkv) + 1 != length(stab) ||
      max(lkv) > maxderiv - 1)  # Bookkeep checks
    stop("'severity.table', 'maxderiv' and ",
         "'lookup' do not match")

  vecTF <- FALSE
  wvec <- if (length(wsdmvec)) wsdmvec else
    wsdm(object, eta0 = eta0,
         subset = subset, hdiff = hdiff,
         maxderiv = maxderiv, ...)
  ans1 <- character(length(wvec))
  names(ans1) <- names(wvec)
  ans1[is.na(wvec)] <- stab[length(stab)]

  if (tx.some && !length(wsdmvec)) {  # object inputted
    links.coef <- linkfun(object, by.var = TRUE)
    vecTF <-
      (wvec < 1 & links.coef == "loglink") |
      (wvec > 1 & links.coef == "cauchitlink")
    wvec[vecTF] <- sqrt(wvec[vecTF])
  }
  if (any(vecTF))
    warning("'loglink' and/or 'cauchitlink' d",
            "etected; setting 'tx.some = TRUE'",
            " may be more accurate")


  for (jay in seq(lkv)) {
    ans1[wvec >= lkv[jay]] <- stab[jay]
  }

  ans1
}  # hdeffsev







cops3 <- function(object, ...) UseMethod("cops3")


cops3.default <-
  function(object, ...)
    stop("there is no default method!!")


cops3.glm <-
  function(object,
           level = 0.999,   #1 - 1e-6,
           beta.range = confint.default(object,
                                        level = level),
           muxrange = 4,  # 2, or 4, or 7, or...
           iter.max = 8,  # Gridsearch  # Not 20,
           tol = 1e-5,  # 4 optim()
           subset = NULL,  # 20250219
           do1 = TRUE,  # F 4 propodds()
           ...) {
  T <- TRUE; F <- FALSE
  if (is.null(subset)) {
    subset <- TRUE
  }

  if (length(muxrange) != 1 || muxrange < 1.1)
    stop("'muxrange' must be > 1.1")
  if (!is.Numeric(iter.max, positive = TRUE,
       length.arg = 1, integer.valued = TRUE) ||
      iter.max <= 2)
    stop("'iter.max' must be an integer > 2")
  if (!isFALSE(do1) && !isTRUE(do1))
    stop("'do1' must be TRUE or FALSE")


  M <- 1   #npred(object)
  cobj <- coef(object)  # Original coeffs
  nmcobj <- names(cobj)
  p.VLM <- length(cobj)
  ind.ok <- if (is.character(subset))
    match(subset, nmcobj) else  # Unsorted
    (1:p.VLM)[subset]  # Sorted
  if (!inherits(object, "glm")) stop("not a 'glm' object")
  nn <- nobs(object)
  Xvlm  <- model.matrix(object)   #, type = "vlm"
  etamat <- predict(object)  # Contains offsets



  has.intercept <- all(Xvlm[, 1] == 1)
  tltmp <- attr(terms(object), "term.labels")
  Hlist <- vector("list", has.intercept + length(tltmp))
  for (i in 1:length(Hlist))
    Hlist[[i]] <- diag(1)
  names(Hlist) <-
    c(if (has.intercept) "(Intercept)" else NULL, tltmp)
  if (has.intercept) {
    if (do1) {
      ind.ok <- unique(c(1:ncol(Hlist[[1]]), ind.ok))
    } else {
      ind.ok <- setdiff(ind.ok, 1:ncol(Hlist[[1]]))
      if (!length(ind.ok))
        stop("No coefficients! Change 'subset' or set ",
             "'do1 = TRUE'")
    }
  }  # has.intercept



  newinfo <-  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    function(beta.try1, jkay)  {
    copy.object <- object 
    copy.object$qr <- NULL
    copy.object$coefficients <- cobj
    copy.object$coefficients[jkay] <- beta.try1
    neweta <- c(Xvlm %*% copy.object$coefficients)
    copy.object$linear.predictors <- neweta
    newmu <- object$family$linkinv(neweta)
    copy.object$fitted.values <- newmu

    variance <- object$family$variance
    mu.eta <- object$family$mu.eta
    pwt <- weights(object, type = "prior")
    newwz <- pwt * (mu.eta(neweta)^2) / variance(newmu)
    UU <- vchol(newwz, M = M, n = nn)  # silent = T
    UtXvlm <- mux111(cc = UU, xmat = Xvlm, M = M,
                     whichj = jkay)
    total.info <- sum((UtXvlm[, jkay])^2)
    total.info
  }  # newinfo


  gridpts <- 1 + 2 * iter.max
  gmat <- matrix(NA_real_, p.VLM, gridpts)
  tmp <- 1:iter.max
  nmgmat <- as.character(sort(c(-tmp, 0, tmp)))
  dimnames(gmat) <- list(nmcobj, nmgmat)
  jmat <- gmat  # To store total info
  beta.range <- cbind(beta.range)  # Wald, silent, 2-cols
  midpts <- gmat[, '0'] <- rowMeans(beta.range)
  gmat[,  '1'] <- beta.range[, 2]
  gmat[, '-1'] <- beta.range[, 1]
  if (iter.max > 1)  # Redundant
    for (iter in 2:iter.max) {
      muxr <- muxrange * iter  # Linear spread
      acit <- as.character( iter)
      Acit <- as.character(-iter)
      gmat[, acit] <- midpts + muxr *
                     (beta.range[, 2] - midpts)
      gmat[, Acit] <- midpts - muxr *
                     (midpts - beta.range[, 1])
    }  # for iter
    for (rr in ind.ok)   #1:nrow(gmat)
      for (cc in 1:ncol(gmat))
        jmat[rr, cc] <- newinfo(gmat[rr, cc], rr)

  initLHS <- initRHS <- numeric(p.VLM)
  for (rr in ind.ok) {  # 1:p.VLM
    wmax <- which.max(jmat[rr, ])
    if (wmax %in% c(1, ncol(jmat)))
      stop("try increasing 'iter.max' and/or 'muxrange'")
    initLHS[rr] <- gmat[rr, wmax - 1]  # Bracket the
    initRHS[rr] <- gmat[rr, wmax + 1]  # maximum
  }  # rr

  copsvec <- cobj  # Named storage
  for (rr in ind.ok) {  # 1:p.VLM
    ofit <- optimize(newinfo, maximum = TRUE,
                     lower = initLHS[rr],
                     upper = initRHS[rr],
                     tol = tol, jkay = rr)
    copsvec[rr] <- ofit$maximum
  }  # rr

  copsvec[ind.ok]
}  # cops3.glm





copsdvglm <-
    function(object, doffset = 0.1,
             ...) {
  if (!is.Numeric(doffset, 1L, positive = TRUE))
    stop("bad input for 'doffset'")
  cc.tmp <- abs(cops(object, ...))
  cf.tmp <- abs(coef(object))
    cc.tmp / (doffset + cf.tmp)
}


if (!isGeneric("copsd"))
  setGeneric("copsd",
             function(object, ...)
                 standardGeneric("copsd"),
             package = "VGAM")


setMethod("copsd", "vglm",
          function(object, ...) {
  copsdvglm(object, ...)
})




copsd3 <- function(object, ...) UseMethod("copsd3")
copsd3.default <-
  function(object, ...)
    stop("there is no default method!!")
copsd3.glm <-
    function(object, doffset = 0.1,
             ...) {
  cc.tmp <- abs(cops3(object, ...))
  cf.tmp <- abs(coef(object))
    cc.tmp / (doffset + cf.tmp)
}



























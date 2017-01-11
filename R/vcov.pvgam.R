# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.














vcov.pvgam <- function(object, ...) {
  vcovpvgam(object, ...)
}



 vcovpvgam <-
 function(object,
          special = FALSE,
          frequentist = FALSE, dispersion = NULL, unconditional = FALSE,
          ...) {





 if (!special) {
   return(vcovvlm(object, ...))
  }



 warning("vcovpvgam() is only 50% finished")


 print("in vcovpvgam; hi 2a")
 print("class(object)")
 print( class(object) )




  M <- npred(object)
  n <- nobs(object, type = "lm")
  wz <- weights(object, type = "working")
  X.vlm.save <- model.matrix(object, type = "vlm")
  U <- vchol(wz, M = M, n = n)
  X.vlm <- mux111(U, X.vlm.save, M = M)
  X.vlm.aug <- rbind(X.vlm,
                     model.matrix(object, type = "penalty"))


  qr1 <- qr(X.vlm.aug)
  qr2 <- qr(X.vlm)
   poststuff <-
     mgcv::magic.post.proc(X.vlm.aug,
                           object = object@ospsslot$magicfit, w = NULL)
  magicfit <- object@ospsslot$magicfit
  rV <- magicfit$rV
  Vb   <- poststuff$Vb
  Ve   <- poststuff$Ve
  hhat <- poststuff$hat
  eedf <- poststuff$edf


  scale.param <- 1  # Assumed


  vc <- if (frequentist) {

    mat1 <- solve(crossprod(qr.R(qr1)))
    scale.param *
        (mat1 %*% crossprod(qr.R(qr2)) %*% mat1)
  } else {
    Vc <- NULL  # Corrected ML or REML is not available.

    Vp  <- scale.param * tcrossprod(solve(qr.R(qr1)))

    Vp2 <- rV %*% t(rV)  # * sig2  # For checking

 print("max(abs(Vp - Vp2)); should be 0")
 print( max(abs(Vp - Vp2)) )


    if (FALSE) {
    He <- SemiParFit$fit$hessian
    He.eig <- eigen(He, symmetric=TRUE)
    Vb <- He.eig$vectors %*%
          tcrossprod(diag(1/He.eig$values),
                     He.eig$vectors)  # this could be taken from magic as well
    Vb <- (Vb + t(Vb) ) / 2


    HeSh <- He - SemiParFit$fit$S.h
    F <- Vb%*%HeSh  # diag(SemiParFit$magpp$edf)


    HeSh <- He
    Ve <- Vb
    F <- F1 <- diag(rep(1,dim(Vb)[1]))
    R <- SemiParFit$bs.mgfit$R
    }



    if (unconditional && !is.null(Vc))
      Vc else Vp
  }  # Bayesian
  if (is.null(dispersion)) {
    sig2 <- 1  # zz
    vc <- summary(object)@dispersion * vc / sig2
  } else {
    sig2 <- summary(object)@dispersion  # 1  # zz
    vc <- dispersion * vc / sig2
  }



 print("head(sort(diag(vc)))")
 print( head(sort(diag(vc))) )
 print("head(sort(diag(Ve)))")
 print( head(sort(diag(Ve))) )


 print("tail(sort(diag(vc)))")
 print( tail(sort(diag(vc))) )
 print("tail(sort(diag(Ve)))")
 print( tail(sort(diag(Ve))) )



 print("head(sort(diag(vc))) / head(sort(diag(Ve)))")
 print( head(sort(diag(vc))) / head(sort(diag(Ve))) )




 print("max(abs(sort(diag(vc)) - sort(diag(Ve))))")
 print( max(abs(sort(diag(vc)) - sort(diag(Ve)))) )





  vc
}






setMethod("vcov", "pvgam",
         function(object, ...)
         vcovpvgam(object, ...))








startstoppvgam <-
  function(object, ...) {

  which.X.sm.osps <- object@ospsslot$sm.osps.list$which.X.sm.osps
  if (!length(which.X.sm.osps))
    stop("no 'sm.os()' or 'sm.ps()' term in 'object'")
  all.ncol.Hk <- unlist(lapply(constraints(object, type = "term"), ncol))
  names.which.X.sm.osps <- names(which.X.sm.osps)
  endf <- rep_len(NA_real_, sum(all.ncol.Hk[names.which.X.sm.osps]))
  names(endf) <- vlabel(names.which.X.sm.osps,
                        all.ncol.Hk[names.which.X.sm.osps],
                        M = npred(object))
  stopstart <- NULL


  iptr <- 1
  iterm <- 1
  for (ii in names(all.ncol.Hk)) {
    if (length(which.X.sm.osps[[ii]])) {
      temp3 <- -1 + iptr + all.ncol.Hk[ii] * length(which.X.sm.osps[[ii]])
      new.index <- iptr:temp3  # Includes all component functions wrt xk
      iptr <- iptr + length(new.index)  # temp3
      mat.index <- matrix(new.index, ncol = all.ncol.Hk[ii], byrow = TRUE)
      for (jay in 1:all.ncol.Hk[ii]) {
        cf.index <- mat.index[, jay]


        stopstart <- c(stopstart, list(cf.index))


        iterm <- iterm + 1
      }  # for
    } else {
      iptr <- iptr + all.ncol.Hk[ii]
    }
  }  # ii
  names(stopstart) <- names(endf)
  stopstart
}











summarypvgam <-
  function(object, dispersion = NULL,
           digits = options()$digits-2,
           presid = TRUE) {

  stuff <- summaryvglm(object, dispersion = dispersion,
           digits = digits,
           presid = presid)

  answer <-
  new("summary.pvgam",
      object,
      call = stuff@call,
      cov.unscaled = stuff@cov.unscaled,
      correlation = stuff@correlation,
      df = stuff@df,
      sigma = stuff@sigma)

  answer@misc$nopredictors <- stuff@misc$nopredictors
  answer@ospsslot <- object@ospsslot


  slot(answer, "coefficients") <- stuff@coefficients  # Replace

  coef3 <- stuff@coef3
  aassign <- attr(model.matrix(object, type = "vlm"),  "assign")
  myterms <- names(object@ospsslot$sm.osps.list$which.X.sm.osps)
  index.exclude <- NULL
  for (ii in myterms) {
    index.exclude <- c(index.exclude, unlist(aassign[[ii]]))
  }
  slot(answer, "coef3") <- coef3[-index.exclude, , drop = FALSE]



  if (is.numeric(stuff@dispersion))
    slot(answer, "dispersion") <- stuff@dispersion

  if (presid) {
    Presid <- residuals(object, type = "pearson")
    if (length(Presid))
      answer@pearson.resid <- as.matrix(Presid)
  }







  pinv <- function(V, M, rank.tol = 1e-6) {
    D <- eigen(V, symmetric = TRUE)
    M1 <- length(D$values[D$values > rank.tol * D$values[1]])
    if (M > M1)
      M <- M1  # avoid problems with zero eigen-values

    if (M+1 <= length(D$values))
      D$values[(M+1):length(D$values)] <- 1
    D$values <- 1 / D$values
    if (M+1 <= length(D$values))
      D$values[(M+1):length(D$values)] <- 0
    res <- D$vectors %*% (D$values * t(D$vectors))  ##D$u%*%diag(D$d)%*%D$v
    attr(res, "rank") <- M
    res
  }  ## end of pinv



  startstop <- startstoppvgam(object)
  m <- length(startstop)

  df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, m)
  names(chi.sq) <- names(startstop)
  p.type <- 5  # Frequentist
  est.disp <- if (is.logical(object@misc$estimated.dispersion))
    object@misc$estimated.dispersion else FALSE

  pvgam.residual.df <- df.residual_pvgam(object)



  for (i in 1:m) {

    p <- coef(as(object, "pvgam"))[(startstop[[i]])]  # params for smooth

    endf <- endfpvgam(object, diag.all = TRUE)  # This is ENDF+1 actually


    edf1[i] <- edf[i] <- sum(endf[(startstop[[i]])])
    if (FALSE && !is.null(object$edf1))
      edf1[i] <- sum(object$edf1[(startstop[[i]])])

    V <- if (p.type == 5) {
      Ve <- vcov(object, special = FALSE)
      Ve[(startstop[[i]]), (startstop[[i]]), drop = FALSE]
    } else {
      Vp <- vcov(object, special = TRUE, frequentist = FALSE)
      Vp[(startstop[[i]]), (startstop[[i]]), drop = FALSE]
    }

    if (p.type == 5) {
      M1 <- length(startstop[[i]])  # zz
      M <- min(M1,
               ceiling(2*sum(endf[(startstop[[i]])]))
               )
      V <- pinv(V, M)
      chi.sq[i] <- t(p) %*% V %*% p
      df[i] <- attr(V, "rank")
    }



    if (p.type == 5) {
      s.pv[i] <- if (est.disp) {
        pf(chi.sq[i] / df[i], df1 = df[i], df2 = pvgam.residual.df,
           lower.tail = FALSE)
      } else {
        pchisq(chi.sq[i], df = df[i], lower.tail = FALSE)
      }
      if (df[i] < 0.1)
        s.pv[i] <- NA
    }


    if (est.disp) {
      if (p.type == 5) {
        s.table <- cbind(edf, df, chi.sq / df, s.pv)
        dimnames(s.table) <- list(names(chi.sq),
                                  c("edf", "Est.rank", "F", "p-value"))
      } else {
        s.table <- cbind(edf, df, chi.sq/df, s.pv)
        dimnames(s.table) <- list(names(chi.sq),
                                  c("edf", "Ref.df", "F", "p-value"))
      }
    } else {
      if (p.type == 5) {
 # This case is commonly executed
        s.table <- cbind(edf, df, chi.sq, s.pv)
        dimnames(s.table) <- list(names(chi.sq),
                                  c("edf", "Est.rank", "Chi.sq", "p-value"))
      } else {
        s.table <- cbind(edf, df, chi.sq, s.pv)
        dimnames(s.table) <- list(names(chi.sq),
                                  c("edf", "Ref.df", "Chi.sq", "p-value"))
      }
    }  # else
  }  # for (i)
  answer@post$s.table <- s.table





  aod <- data.frame(message = 'this does not work yet')
  slot(answer, "anova") <- aod

  answer
}  # summarypvgam()







show.summary.pvgam <-
  function(x, quote = TRUE, prefix = "",
           digits = options()$digits-2,
           signif.stars = getOption("show.signif.stars")) {




  show.summary.vglm(x, quote = quote, prefix = prefix,
                    digits = digits, top.half.only = TRUE)



  startstop <- startstoppvgam(x)
  m <- length(startstop)
  s.table <- x@post$s.table
  if (0 < m && length(s.table)) {
    cat("\nApproximate significance of smooth terms:\n")
    printCoefmat(s.table, digits = digits,
                 signif.stars = signif.stars, has.Pvalue = TRUE,
                 na.print = "NA", cs.ind = 1)
  }






  M <- x@misc$M


  Presid <- x@pearson.resid
  rdf <- x@df[2]


  cat("\nNumber of linear/additive predictors:   ", M, "\n")

  if (!is.null(x@misc$predictors.names))
  if (M == 1)
    cat("\nName of linear/additive predictor:",
        paste(x@misc$predictors.names, collapse = ", "), "\n") else
  if (M <= 5)
    cat("\nNames of linear/additive predictors:",
        paste(x@misc$predictors.names, collapse = ", "), "\n")

  prose <- ""
  if (length(x@dispersion)) {
    if (is.logical(x@misc$estimated.dispersion) &&
        x@misc$estimated.dispersion) {
      prose <- "(Estimated) "
    } else {
      if (is.numeric(x@misc$default.dispersion) &&
          x@dispersion == x@misc$default.dispersion)
        prose <- "(Default) "

      if (is.numeric(x@misc$default.dispersion) &&
          x@dispersion != x@misc$default.dispersion)
        prose <- "(Pre-specified) "
    }
    cat(paste("\n", prose, "Dispersion Parameter for ",
        x@family@vfamily[1],
        " family:   ",
        format(round(x@dispersion, digits)), "\n", sep = ""))
  }

  if (length(deviance(x)))
    cat("\nResidual deviance: ", format(round(deviance(x), digits)),
        "on", format(round(rdf, 3)), "degrees of freedom\n")

  if (length(logLik.vlm(x)))
    cat("\nLog-likelihood:", format(round(logLik.vlm(x), digits)),
        "on", format(round(rdf, 3)), "degrees of freedom\n")

  if (length(x@criterion)) {
    ncrit <- names(x@criterion)
    for (ii in ncrit)
      if (ii != "loglikelihood" && ii != "deviance")
        cat(paste(ii, ":", sep = ""), format(x@criterion[[ii]]), "\n")
  }




  if (is.Numeric(x@ospsslot$iter.outer)) {
    cat("\nNumber of outer iterations: ", x@ospsslot$iter.outer, "\n")
    cat("\nNumber of IRLS iterations at final outer iteration: ", x@iter,
        "\n")
  } else {
    cat("\nNumber of IRLS iterations: ", x@iter, "\n")
  }


  if (FALSE && length(x@anova)) {
    show.vanova(x@anova, digits = digits)   # ".vanova" for Splus6
  }

  invisible(NULL)
}  # show.summary.pvgam()








setMethod("summary", "pvgam",
          function(object, ...)
          summarypvgam(object, ...))



setMethod("show", "summary.pvgam",
          function(object)
          show.summary.pvgam(object))






psintpvgam <- function(object, ...) {
  object@ospsslot$sm.osps.list$ps.int
}


if (!isGeneric("psint"))
  setGeneric("psint", function(object, ...)
              standardGeneric("psint"),
             package = "VGAM")


setMethod("psint", "pvgam",
         function(object, ...)
         psintpvgam(object, ...))







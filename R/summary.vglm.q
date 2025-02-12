# These functions are
# Copyright (C) 1998-2025 T.W. Yee, University of Auckland.
# All rights reserved.













summaryvglm <-
  function(object, correlation = FALSE,
           dispersion = NULL, digits = NULL,
           presid = FALSE,  # TRUE,
           HDEtest = TRUE,  # Added 20180203
           hde.NA = TRUE,
           threshold.hde = 0.001,
           signif.stars = getOption("show.signif.stars"),
           nopredictors = FALSE,
           lrt0.arg = FALSE,
           score0.arg = FALSE,
           wald0.arg = FALSE,
           values0 = 0,
           subset = NULL,
           omit1s = TRUE,
           wsdm.arg = FALSE,  # 20241116
           hdiff = 0.005,  # 20250110
           retry = TRUE,   # FALSE,
           mux.hdiff = 1,
           eps.wsdm = 0.15,
           Mux.div = 3,
           doffset.wsdm = NULL,  # 20241120
           ...  # Added 20151211
          ) {


  missing.HDEtest <- missing(HDEtest)









  if (length(dispersion) &&
      dispersion == 0 &&
      length(object@family@summary.dispersion) &&
      !object@family@summary.dispersion) {
    stop("cannot use the general VGLM formula (based on a residual ",
         "sum of squares) for computing the dispersion parameter")
  }



  stuff <- summaryvlm(
                      object,

                      presid = FALSE,

                      correlation = correlation,
                      dispersion = dispersion,

                      lrt0.arg   = lrt0.arg,
                      score0.arg = score0.arg,
                      wald0.arg  = wald0.arg,
                      values0 = values0,
                      subset = subset,
                      omit1s = omit1s)





  infos.fun <- object@family@infos
  infos.list <- infos.fun()
  summary.pvalues <- if (is.logical(infos.list$summary.pvalues))
    infos.list$summary.pvalues else TRUE
  if (!summary.pvalues && ncol(stuff@coef3) == 4)
    stuff@coef3 <- stuff@coef3[, -4]  # Del pvalues coln
  if (length(idoff <- infos.list$doffset) &&
      length(doffset.wsdm) == 0 &&
      length(dlink <- object@misc$link) == 1 &&
      any(colnames(infos.list$doffset) ==
      object@misc$link)) {
    doffset.wsdm <- idoff[, dlink]
  }
  if (!length(doffset.wsdm))
    doffset.wsdm <- 1  # Final alternative





  answer <-
  new("summary.vglm",
      object,
      coef3 = stuff@coef3,
      coef4lrt0 = stuff@coef4lrt0,  # Might be an empty "matrix"
      coef4score0 = stuff@coef4score0,  # Might be an empty "matrix"
      coef4wald0 = stuff@coef4wald0,  # Might be an empty "matrix"
      cov.unscaled = stuff@cov.unscaled,
      correlation = stuff@correlation,
      df = stuff@df,
      sigma = stuff@sigma)


  if (presid) {
    Presid <- resid(object, type = "pearson")
    if (length(Presid))
      answer@pearson.resid <- as.matrix(Presid)
  }

  slot(answer, "misc") <- stuff@misc  # Replace


  answer@misc$signif.stars <- signif.stars  # 20140728
  answer@misc$nopredictors <- nopredictors  # 20150831


  if (is.numeric(stuff@dispersion))
    slot(answer, "dispersion") <- stuff@dispersion






  try.this <- findFirstMethod("summaryvglmS4VGAM",
                              object@family@vfamily)
  if (length(try.this)) {
    new.postslot <-
    summaryvglmS4VGAM(object = object,
                      VGAMff = new(try.this),
                      ...)
    answer@post <- new.postslot
  } else {
  }



  control <- object@control
  if (missing.HDEtest &&
      length(temp <- object@control$summary.HDEtest)) {
    HDEtest <- temp
  }

  if (HDEtest) {
    answer@post$hdeff <- hdeff(object, derivative = 1, se.arg = TRUE)
    answer@post$hde.NA <- hde.NA
    answer@post$threshold.hde <- threshold.hde
  }



  if (!isFALSE(wsdm.arg) && !isTRUE(wsdm.arg))
    stop("bad input for argument 'wsdm.arg'")
  answer@post$wsdm.arg <- wsdm.arg  # Save it
  if (!wsdm.arg) return(answer)  # Exit now.
  dmax.WSDM <- 4  # max.deriv; zzzz
  upsvec <-
    wsdm(object, maxderiv = dmax.WSDM,
         doffset = doffset.wsdm,
         hdiff = hdiff,  # Do via attr:
         retry = retry,
         mux.hdiff = mux.hdiff,
         subset = subset,
         eps.wsdm = eps.wsdm,
         Mux.div = Mux.div,
         warn.retry = FALSE)  # No warning here
  answer@coef3 <-  # Insert WSDM as 4th coln
    cbind(answer@coef3,
          "WSDM" = round(upsvec, digits = 2))
    answer@post$WSDM <- upsvec  # No round()
    answer@post$max.deriv.WSDM <- dmax.WSDM
    answer@post$wsdm.arg <- wsdm.arg


  answer
}  # summaryvglm





setMethod("summaryvglmS4VGAM",
          signature(VGAMff = "cumulative"),
  function(object, VGAMff, ...) {
   object@post <-
     callNextMethod(VGAMff = VGAMff,
                    object = object, ...)
  object@post$reverse <- object@misc$reverse
  cfit <- coef(object, matrix = TRUE)
  M <- ncol(cfit)
  if (rownames(cfit)[1] ==  "(Intercept)" &&
  identical(constraints(object)[[1]], diag(M)))
object@post$expcoeffs <- exp(coef(object)[-(1:M)])
  object@post
})



setMethod("showsummaryvglmS4VGAM",
          signature(VGAMff = "cumulative"),
  function(object, VGAMff, ...) {

  if (length(object@post$expcoeffs)) {
    cat("\nExponentiated coefficients:\n")
    print(object@post$expcoeffs)
  }
  if (FALSE) {
    if (object@post$reverse)
    cat("Reversed\n\n") else
    cat("Not reversed\n\n")
  }
})





setMethod("summaryvglmS4VGAM",
          signature(VGAMff = "cratio"),
  function(object, VGAMff, ...) {
   object@post <-
     callNextMethod(VGAMff = VGAMff,
                    object = object, ...)
  object@post$reverse <- object@misc$reverse
  cfit <- coef(object, matrix = TRUE)
  M <- ncol(cfit)
  if (rownames(cfit)[1] ==  "(Intercept)" &&
  identical(constraints(object)[[1]], diag(M)))
object@post$expcoeffs <- exp(coef(object)[-(1:M)])
  object@post
})

setMethod("showsummaryvglmS4VGAM",
          signature(VGAMff = "cratio"),
  function(object, VGAMff, ...) {
  if (length(object@post$expcoeffs)) {
    cat("\nExponentiated coefficients:\n")
    print(object@post$expcoeffs)
  }
})



setMethod("summaryvglmS4VGAM",
          signature(VGAMff = "sratio"),
  function(object, VGAMff, ...) {
   object@post <-
     callNextMethod(VGAMff = VGAMff,
                    object = object, ...)
  object@post$reverse <- object@misc$reverse
  cfit <- coef(object, matrix = TRUE)
  M <- ncol(cfit)
  if (rownames(cfit)[1] ==  "(Intercept)" &&
  identical(constraints(object)[[1]], diag(M)))
object@post$expcoeffs <- exp(coef(object)[-(1:M)])
  object@post
})

setMethod("showsummaryvglmS4VGAM",
          signature(VGAMff = "sratio"),
  function(object, VGAMff, ...) {
  if (length(object@post$expcoeffs)) {
    cat("\nExponentiated coefficients:\n")
    print(object@post$expcoeffs)
  }
})





setMethod("summaryvglmS4VGAM",
          signature(VGAMff = "acat"),
  function(object, VGAMff, ...) {
   object@post <-
     callNextMethod(VGAMff = VGAMff,
                    object = object, ...)
  object@post$reverse <- object@misc$reverse
  cfit <- coef(object, matrix = TRUE)
  M <- ncol(cfit)
  if (rownames(cfit)[1] ==  "(Intercept)" &&
  identical(constraints(object)[[1]], diag(M)))
object@post$expcoeffs <- exp(coef(object)[-(1:M)])
  object@post
})

setMethod("showsummaryvglmS4VGAM",
          signature(VGAMff = "acat"),
  function(object, VGAMff, ...) {
  if (length(object@post$expcoeffs)) {
    cat("\nExponentiated coefficients:\n")
    print(object@post$expcoeffs)
  }
})







setMethod("summaryvglmS4VGAM",
          signature(VGAMff = "multinomial"),
  function(object,
           VGAMff,
           ...) {
   object@post <-
     callNextMethod(VGAMff = VGAMff,
                    object = object,
                    ...)
  object@post$refLevel <- object@misc$refLevel
  object@post
})



setMethod("showsummaryvglmS4VGAM",
          signature(VGAMff = "multinomial"),
  function(object,
           VGAMff,
           ...) {
  if (length(object@post$refLevel))
    cat("\nReference group is level ",
        object@post$refLevel, " of the response\n")
  callNextMethod(VGAMff = VGAMff,
                 object = object,
                 ...)
})



setMethod("summaryvglmS4VGAM",
          signature(VGAMff = "VGAMcategorical"),
  function(object,
           VGAMff,
           ...) {
  object@post
})


setMethod("showsummaryvglmS4VGAM",
          signature(VGAMff = "VGAMcategorical"),
  function(object,
           VGAMff,
           ...) {
})









setMethod("logLik",  "summary.vglm",
          function(object, ...)
  logLik.vlm(object, ...))




show.summary.vglm <-
  function(x,
       digits = max(3L, getOption("digits") - 3L),  # Same as glm()
       quote = TRUE,
       prefix = "",
       presid = length(x@pearson.resid) > 0,  # FALSE,  # TRUE,
       HDEtest = TRUE,
       hde.NA = TRUE,
       threshold.hde = 0.001,
       signif.stars = NULL,   # Use this if logical; 20140728
       nopredictors = NULL,   # Use this if logical; 20150831
       top.half.only = FALSE,  # Added 20160803
       ...  # Added 20151214
       ) {



  M <- x@misc$M
  coef3 <- x@coef3  # icients
  correl <- x@correlation

  cn3 <- colnames(coef3)
  if (cn3[length(cn3)] == "WSDM")
    coef3 <- coef3[, -length(cn3), drop = FALSE]

  digits <- if (is.null(digits))
    options()$digits - 2 else digits

  cat("\nCall:\n",
      paste(deparse(x@call),
                    sep = "\n", collapse = "\n"),
      "\n", sep = "")




  if (HDEtest)
  if (is.logical(x@post$hde.NA) && x@post$hde.NA) {
    if (length(hado <- x@post$hdeff)) {
      HDE <- is.Numeric(hado[, "deriv1"]) &  # Could be all NAs
             hado[, "deriv1"] < 0
      if (any(HDE) && ncol(coef3) == 4) {
        HDE <- HDE & (x@post$threshold.hde < coef3[, 4])
        coef3[HDE, 3:4] <- NA  # 3:4 means WaldStat and p-value
      }
    }
  }  # is.logical(x@post$hde.NA) && x@post$hde.NA




  Presid <- x@pearson.resid
  rdf <- x@df[2]
  pearres.out <- FALSE
  if (presid &&
      length(Presid) &&
      all(!is.na(Presid)) &&
      is.finite(rdf)) {

    if (rdf/M > 5) {
      rq <- apply(as.matrix(Presid), 2, quantile)  # 5 x M
      dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", "Max"),
                           x@misc$predictors.names)
      pearres.out <- TRUE
      cat("\nPearson residuals:\n")
      print(t(rq), digits = digits)
    } else
    if (rdf > 0) {
      pearres.out <- TRUE
      cat("\nPearson residuals:\n")
      print(Presid, digits = digits)
    }
  }


  use.signif.stars <- if (is.logical(signif.stars))
    signif.stars else x@misc$signif.stars  # 20140728
  if (!is.logical(use.signif.stars))
    use.signif.stars <- getOption("show.signif.stars")


  use.nopredictors <- if (is.logical(nopredictors))
    nopredictors else x@misc$nopredictors  # 20140728
  if (!is.logical(use.nopredictors)) {
    warning("cannot determine 'nopredictors'; choosing FALSE")
    use.nopredictors <- FALSE
  }


  Situation <- -1
  how.many <- c(length(x@coef4lrt0),
                length(x@coef4score0),
                length(x@coef4wald0))

  if (length(x@coef4lrt0)) {  # && wald0.arg
    cat(if (top.half.only)
    "\nParametric coefficients:" else
    "\nLikelihood ratio test coefficients:", "\n")
    printCoefmat(x@coef4lrt0, digits = digits,
      signif.stars = use.signif.stars,
      na.print = "NA",
      P.values = TRUE, has.Pvalue = TRUE,
 signif.legend = sum(how.many[-1]) == 0)  # Last 1
    Situation <- 3
  }



  if (length(x@coef4score0)) {  # && wald0.arg
    cat(if (top.half.only)
        "\nParametric coefficients:" else
        "\nRao score test coefficients:", "\n")
    printCoefmat(x@coef4score0, digits = digits,
      signif.stars = use.signif.stars,
      na.print = "NA",
 signif.legend = sum(how.many[3]) == 0)  # Last 1
    Situation <- 4
  }



  if (length(x@coef4wald0)) {  # && wald0.arg
    cat(if (top.half.only)
      "\nParametric coefficients:" else
      "\nWald (modified by IRLS iterations) coefficients:", "\n")
    printCoefmat(x@coef4wald0, digits = digits,
                 signif.stars = use.signif.stars,
                 na.print = "NA")
    Situation <- 1
  } else
  if (length(coef3) && Situation < 0) {
    cat(if (top.half.only)
        "\nParametric coefficients:" else
        "\nCoefficients:", "\n")
    printCoefmatt(coef3, digits = digits,
                  signif.stars = use.signif.stars,
                  wsdmvec = x@post$WSDM,
                  dig.wsdm = 2,
                  na.print = "NA")
    Situation <- 2
  }  # length(coef3)







  if (top.half.only)
    return(invisible(NULL))



  if (M >= 5 && !pearres.out)
  cat("\nNumber of linear predictors: ", M, "\n")

  if (!is.null(x@misc$predictors.names) &&
      !use.nopredictors &&
      !pearres.out) {
    if (M == 1) {
      cat("\nName of linear predictor:",
          paste(x@misc$predictors.names, collapse = ", "), "\n")
    } else
    if (M <= 12) {
      LLL <- length(x@misc$predictors.names)
      cat("\nNames of linear predictors:",
          if (LLL == 1)
            x@misc$predictors.names else
        c(paste0(x@misc$predictors.names[-LLL], sep = ","),
          x@misc$predictors.names[LLL]),  fill = TRUE)
    }
  }

  prose <- ""
  if (length(x@dispersion)) {
    if (is.logical(x@misc$estimated.dispersion) &&
       x@misc$estimated.dispersion) {
      prose <- "(Estimated) "
    }  else {
      if (is.numeric(x@misc$default.dispersion) &&
          x@dispersion == x@misc$default.dispersion)
        prose <- "(Default) "

      if (is.numeric(x@misc$default.dispersion) &&
          x@dispersion != x@misc$default.dispersion)
        prose <- "(Pre-specified) "
    }
    if (any(x@dispersion != 1))
      cat(paste0("\n", prose,
          "Dispersion Parameter for ",
         x@family@vfamily[1], " family:   ",
         yformat(x@dispersion, digits), "\n"))
  }


  if (length(deviance(x))) {
    cat("\nResidual deviance:",
        yformat(deviance(x), digits))
    if (is.finite(rdf))
      cat(" on", round(rdf, digits),
          "degrees of freedom\n") else
      cat("\n")
  }


  if (length(vll <- logLik.vlm(x))) {
    cat("\nLog-likelihood:", yformat(vll, digits))
    if (is.finite(rdf))
      cat(" on", round(rdf, digits),
          "degrees of freedom\n") else
      cat("\n")
  }


  if (length(x@criterion)) {
    ncrit <- names(x@criterion)
    for (ii in ncrit)
      if (ii != "loglikelihood" && ii != "deviance")
        cat(paste0(ii, ":"),
            yformat(x@criterion[[ii]], digits), "\n")
  }


  cat("\nNumber of Fisher scoring iterations:",
      format(trunc(x@iter)), "\n\n")


  if (!is.null(correl)) {
    ncol.X.vlm <- dim(correl)[2]
    if (ncol.X.vlm > 1) {
      cat("Correlation of Coefficients:\n\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1,  -ncol.X.vlm, drop = FALSE], quote = FALSE,
            digits = digits)
      cat("\n")
    }
  }





    if (HDEtest)
    if (Situation == 2 &&
        length(hado <- x@post$hdeff)) {
    if (is.Numeric(hado[, "deriv1"]) &  # Could be all NAs
        all(hado[, "deriv1"] > 0))
      cat("No Hauck-Donner effect found in any of the estimates\n\n")
    if (is.Numeric(hado[, "deriv1"]) &  # Could be all NAs
        any(hado[, "deriv1"] < 0)) {
      cat("Warning: Hauck-Donner effect detected",
          "in the following estimate(s):\n")
      cat(paste0("'",
            rownames(hado)[hado[, "deriv1"] < 0],
                "'", collapse = ", "))
      cat("\n\n")
    }


    if (any(colnames(x@coef3) == "WSDM") &&
        isTRUE(x@post$wsdm.arg) &&
        isFALSE(x@misc$intercept.only)) {
      p.1 <- ncol(x@constraints[["(Intercept)"]])
      maxwsdm <- max(x@post$WSDM[-seq(p.1)])

      seems.okay <- attr(x@post$WSDM, "seems.okay")
      if (is.null(seems.okay))  # In case
        seems.okay <- NA
      proviso <- ", but with uncertain accuracy"
      if (isTRUE(seems.okay)) proviso <- ""
      if (isFALSE(seems.okay))
        proviso <- ", but inaccurate"
      cat("Max-WSDM (excluding intercepts",
          proviso, "): ",
          sprintf("%5.3f", maxwsdm),
      if (abs(maxwsdm + 1 - x@post$max.deriv.WSDM)
          < 0.05) "+" else "",  # Right-censored?
          "\n\n", sep = "")
    }
  }  # Situation == 2 && length(hado)





  try.this <-
      findFirstMethod("showsummaryvglmS4VGAM",
                      x@family@vfamily)
  if (length(try.this)) {
    showsummaryvglmS4VGAM(object = x,
            VGAMff = new(try.this),
            ...)
  } else {
  }



  invisible(NULL)
}  # show.summary.vglm




setMethod("summary", "vglm",
          function(object, ...)
          summaryvglm(object, ...))






setMethod("show", "summary.vglm",
          function(object)
          show.summary.vglm(object))










if (FALSE)
show.summary.binom2.or <-
  function(x,
           digits = max(3L, getOption("digits") - 3L)  # Same as glm()
          ) {

  if (length(x@post$oratio) == 1 &&
      is.numeric(x@post$oratio)) {
    cat("\nOdds ratio: ", round(x@post$oratio, digits), "\n")
  }
}




if (FALSE)
setMethod("show", "summary.binom2.or",
          function(object)
          show.summary.vglm(object))








vcovdefault <- function(object, ...) {
  if (is.null(object@vcov))
    stop("no default")
  object@vcov
}






vcov.vlm <- function(object, ...) {

  vcovvlm(object, ...)
}  # vcov.vlm







 vcovvlm <-
   function(object,
            dispersion = NULL,
            untransform = FALSE,
            complete = TRUE,
            ...    # This line added 20230309
           ) {







  so <- summaryvlm(object, correlation = FALSE,
                   presid = FALSE,
                   dispersion = dispersion)
  d <- if (any(slotNames(so) == "dispersion") &&
           is.Numeric(so@dispersion))
       so@dispersion else 1
  answer <- d * so@cov.unscaled

  if (is.logical(OKRC <- object@misc$RegCondOK) && !OKRC)
    warning("MLE regularity conditions  violated ",
    "at the final iteration of the fitted object")

  if (!untransform)
    return(answer)





  new.way <- TRUE



  if (!is.logical(object@misc$intercept.only))
    stop("cannot determine whether the object is",
         "an intercept-only fit, i.e., 'y ~ 1' is the response")
  if (!object@misc$intercept.only)
    stop("object must be an intercept-only fit,",
         " i.e., y ~ 1 is the response")

  if (!all(trivial.constraints(constraints(object)) == 1))
    stop("object must have trivial constraints")

  M <- object@misc$M




  tvector <- numeric(M)
  etavector <- predict(object)[1, ]  # Contains
  LINK <- object@misc$link
  EARG <- object@misc$earg  # This could be a NULL
  if (is.null(EARG))
    EARG <- list(theta = NULL)
  if (!is.list(EARG))
    stop("'object@misc$earg' must be a list")




  if (length(LINK) != M &&
      length(LINK) != 1)
    stop("cannot obtain the link functions ",
         "to untransform 'object'")



  if (!is.character(LINK))
    stop("the 'link' component of 'object@misc'",
         " should be a character vector")

  learg <- length(EARG)
  llink <- length(LINK)
  if (llink != learg)
    stop("the 'earg' component of 'object@misc'",
         " should be a list of length ", learg)


  level1 <- length(EARG) > 3 &&
            length(intersect(names(EARG),
                c("theta", "inverse",
                  "deriv", "short", "tag"))) > 3
  if (level1)
    EARG <- list(oneOnly = EARG)



  learg <- length(EARG)
  for (ii in 1:M) {
    TTheta <- etavector[ii]  # Transformed theta

    use.earg      <-
      if (llink == 1) EARG[[1]] else EARG[[ii]]
    function.name <-
      if (llink == 1) LINK else LINK[ii]


    if (new.way) {
      use.earg[["inverse"]] <- TRUE  # New
      use.earg[["theta"]] <- TTheta  # New
      Theta <- do.call(function.name, use.earg)


      use.earg[["inverse"]] <- TRUE  # Reset this


      use.earg[["deriv"]] <- 1  # New
      use.earg[["theta"]] <- Theta  # Renew this
      tvector[ii] <- do.call(function.name, use.earg)
    } else {
      stop("link funs handled in the new way now")

    }
  }  # of for (ii in 1:M)

  tvector <- abs(tvector)
  answer <- (cbind(tvector) %*% rbind(tvector)) * answer
  if (length(dmn2 <- names(object@misc$link)) == M)
    dimnames(answer) <- list(dmn2, dmn2)
  answer
}  # vcovvlm







setMethod("vcov", "vlm",
          function(object, ...) vcovvlm(object, ...))


setMethod("vcov", "vglm",
         function(object, ...)
         vcovvlm(object, ...))








yformat <- function(x, digits = options()$digits) {
    format(ifelse(abs(x) < 0.001,
                  signif(x, digits),
                  round(x, digits)))
}











printCoefmatt <- function(x,
  digits = max(3L, getOption("digits") - 2L),
  signif.stars = getOption("show.signif.stars"), 
  signif.legend = signif.stars,
  dig.tst = max(1L, min(5L, digits - 1L)),
  cs.ind = 1:k,
  tst.ind = k + 1, zap.ind = integer(), 
  P.values = NULL,
  has.Pvalue = nc >= 4L &&
      length(cn <- colnames(x)) && 
      substr(cn[nc], 1L, 3L) %in% c("Pr(", "p-v"),
  eps.Pvalue = .Machine$double.eps, 
  na.print = "NA", quote = FALSE, right = TRUE,
  wsdmvec = NULL, dig.wsdm = 2,
  ...) {
  if (is.null(d <- dim(x)) || length(d) != 2L) 
      stop("'x' must be coefficient matrix/data frame")
  nc <- d[2L]
  if (is.null(P.values)) {
    scp <- getOption("show.coef.Pvalues")
    if (!is.logical(scp) || is.na(scp)) {
      warning("option \"show.coef.Pvalues\" is invalid:",
              " assuming TRUE")
      scp <- TRUE
    }
    P.values <- has.Pvalue && scp
  } else if (P.values && !has.Pvalue) 
      stop("'P.values' is TRUE but 'has.Pvalue' is not")
  if (has.Pvalue && !P.values) {
    d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
    nc <- nc - 1
    has.Pvalue <- FALSE
  } else xm <- data.matrix(x)
  k <- nc - has.Pvalue - (if (missing(tst.ind)) 1 else
    length(tst.ind))
  if (!missing(cs.ind) && length(cs.ind) > k) 
      stop("wrong k / cs.ind")
  Cf <- array("", dim = d, dimnames = dimnames(xm))
  ok <- !(ina <- is.na(xm))
  for (i in zap.ind)
    xm[, i] <- zapsmall(xm[, i], digits)
  if (length(cs.ind)) {
    acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
    if (any(ia <- is.finite(acs))) {
      digmin <- 1 + if (length(acs <-
                    acs[ia & acs != 0]))
        floor(log10(range(acs[acs != 0],
              finite = TRUE))) else 0
      Cf[, cs.ind] <- format(round(coef.se,
                             max(1L, digits - 
        digmin)), digits = digits)
    }
  }
  if (length(tst.ind)) 
    Cf[, tst.ind] <- format(round(xm[, tst.ind],
                                  digits = dig.tst), 
        digits = digits)
  if (any(r.ind <- !((1L:nc) %in%
      c(cs.ind, tst.ind, if (has.Pvalue) nc))))
    for (i in which(r.ind))
      Cf[, i] <- format(xm[, i], digits = digits)
  ok[, tst.ind] <- FALSE
  okP <- if (has.Pvalue) ok[, -nc] else ok
  x1 <- Cf[okP]
  dec <- getOption("OutDec")
  if (dec != ".") 
      x1 <- chartr(dec, ".", x1)
  x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
  if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
    Cf[okP][not.both.0] <- format(xm[okP][not.both.0],
                                  digits = max(1L, 
        digits - 1L))
  }
  if (any(ina)) 
    Cf[ina] <- na.print
  if (any(inan <- is.nan(xm))) 
    Cf[inan] <- "NaN"
  if (P.values) {
    if (!is.logical(signif.stars) ||
         is.na(signif.stars)) {
      warning("option \"show.signif.stars\" is ",
              "invalid: assuming TRUE")
      signif.stars <- TRUE
    }
    if (any(okP <- ok[, nc])) {
      pv <- as.vector(xm[, nc])
      Cf[okP, nc] <- format.pval(pv[okP],
        digits = dig.tst, eps = eps.Pvalue)
signif.stars <- signif.stars && any(pv[okP] < 0.1)
      if (signif.stars) {
    Signif <- symnum(pv, corr = FALSE, na = FALSE, 
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
          symbols = c("***", "**", "*", ".", " "))
        Cf <- cbind(Cf, format(Signif))
      }
    } else signif.stars <- FALSE
  } else signif.stars <- FALSE


  do.wsdm <- length(wsdmvec) && is.Numeric(wsdmvec)
  if (do.wsdm) {
   if (signif.stars &&  # Last one is ""
       colnames(Cf)[length(colnames(Cf))] != "")
     stop("confused about colnames(Cf)")
   use.wsdm <- wsdmvec



   mydigits <- 2   # pmin(2, digits - 2)
   use.wsdm[use.wsdm < 0.5 / 10^mydigits] <- 0

   Cf <- if (signif.stars) {  # Retain stars @ RHS
     ncc <- NCOL(Cf)   # ncol(Cf)
     cbind(Cf[, seq(ncc - 1), drop = FALSE],
           "WSDM" = yformat(use.wsdm, dig.wsdm),
           Cf[, ncc, drop = FALSE])
   } else {  # Last coln is WSDM
     cbind(Cf, "WSDM" = yformat(use.wsdm, dig.wsdm))
   }
  }  # do.wsdm


  print.default(Cf, quote = quote, right = right,
                na.print = na.print,  ...)
  if (signif.stars && signif.legend) {
    if ((w <- getOption("width")) <
         nchar(sleg <- attr(Signif, 
      "legend"))) 
      sleg <- strwrap(sleg, width = w - 2,
                      prefix = "  ")
    cat("---\nSignif. codes:  ", sleg, sep = "",
        fill = w + 4 + max(nchar(sleg, "bytes") -
               nchar(sleg)))
  }
  invisible(x)
}  # printCoefmatt













wsdm <-
  function(object,
       hdiff = 0.005,  # Recycled to length >= 2
       retry = TRUE,   # FALSE,
       mux.hdiff = 1,
       maxderiv = 5,   # 0:maxderiv
       theta0 = 0,  # Recycled 20210406
       use.hdeff = FALSE,
       doffset = NULL,  # denominator offset
       subset = NULL,  # numeric, logical, chars
       derivs.out = FALSE,
       fixed.hdiff = TRUE,
       eps.wsdm = 0.15,
       Mux.div = 3,
       warn.retry = TRUE,
       with1 = TRUE, ...) {
  seems.okay <- NA  # Unsure
  hdiff <- hdiff * mux.hdiff  # Adjustment
  if (!is.Numeric(eps.wsdm, positive = TRUE,
                  length.arg = 1))
    stop("bad 'eps.wsdm'")
  if (!is.Numeric(Mux.div - 1, positive = TRUE,
                  length.arg = 1))
    stop("bad 'Mux.div'")

  T <- TRUE; F <- FALSE
  infos.fun <- object@family@infos
  infos.list <- infos.fun()
  lambdas <- if (length(doffset)) doffset else {
    if (length(idoff <- infos.list$doffset) &&
        length(dlink <- object@misc$link) == 1 &&
        any(colnames(infos.list$doffset) ==
        object@misc$link))
      idoff[, dlink] else 1
  }





  fdderiv <-
    function(which.d = 1,  # which.deriv
             mat.Wald.deriv,
             mat.Wald.tmp,
             hdeff.output = NULL,
             use.hdiff = 0.005) {
  if (which.d == 0) return(mat.Wald.deriv)  # Done
  if (which.d >= 9) stop("excessive which.d")
  if (which.d == 1)
    mat.Wald.deriv[, "1"] <- (
        mat.Wald.tmp[, "1"]
      - mat.Wald.tmp[, "0"]) / use.hdiff
  if (which.d == 2)
    mat.Wald.deriv[, "2"] <- (
      - mat.Wald.tmp[, "0"] * 2
      + mat.Wald.tmp[, "1"]
      + mat.Wald.tmp[, "2"]) / use.hdiff^2
  if (which.d == 3)
      mat.Wald.deriv[, "3"] <- (
        mat.Wald.tmp[, "0"] * 3
      - mat.Wald.tmp[, "1"] * 3
      - mat.Wald.tmp[, "2"]
      + mat.Wald.tmp[, "3"]) / use.hdiff^3
  if (which.d == 4)
      mat.Wald.deriv[, "4"] <- (
        mat.Wald.tmp[, "0"] * 6
      - mat.Wald.tmp[, "1"] * 4
      - mat.Wald.tmp[, "2"] * 4
      + mat.Wald.tmp[, "3"]
      + mat.Wald.tmp[, "4"]) / use.hdiff^4
  if (which.d == 5)  # 20241120
      mat.Wald.deriv[, "5"] <- (
      - mat.Wald.tmp[, "0"] * 10
      + mat.Wald.tmp[, "1"] * 10
      + mat.Wald.tmp[, "2"] * 5
      - mat.Wald.tmp[, "3"] * 5
      - mat.Wald.tmp[, "4"]
      + mat.Wald.tmp[, "5"]) / use.hdiff^5
  if (which.d == 6)  # 20241120
      mat.Wald.deriv[, "6"] <- (
      - mat.Wald.tmp[, "0"] * 20
      + mat.Wald.tmp[, "1"] * 15
      + mat.Wald.tmp[, "2"] * 15
      - mat.Wald.tmp[, "3"] * 6
      - mat.Wald.tmp[, "4"] * 6
      + mat.Wald.tmp[, "5"]
      + mat.Wald.tmp[, "6"]) / use.hdiff^6
  if (which.d == 7)  # 20241120
      mat.Wald.deriv[, "7"] <- (
        mat.Wald.tmp[, "0"] * 35
      - mat.Wald.tmp[, "1"] * 35
      - mat.Wald.tmp[, "2"] * 21
      + mat.Wald.tmp[, "3"] * 21
      + mat.Wald.tmp[, "4"] * 7
      - mat.Wald.tmp[, "5"] * 7
      - mat.Wald.tmp[, "6"]
      + mat.Wald.tmp[, "7"]) / use.hdiff^7
  if (which.d == 8)  # 20241120
      mat.Wald.deriv[, "8"] <- (
        mat.Wald.tmp[, "0"] * 70
      - mat.Wald.tmp[, "1"] * 56
      - mat.Wald.tmp[, "2"] * 56
      + mat.Wald.tmp[, "3"] * 28
      + mat.Wald.tmp[, "4"] * 28
      - mat.Wald.tmp[, "5"] * 8
      - mat.Wald.tmp[, "6"] * 8
      + mat.Wald.tmp[, "7"]
      + mat.Wald.tmp[, "8"]) / use.hdiff^8
    return(mat.Wald.deriv)
  }  # fdderiv
  doone <-  # One value of dirr := \pm1, \pm2, ...
    function(dirrval,  # e.g., 2 for -1; direction
             mat.coef.tmp,
             mat.Stdr.tmp) {


  for (kay in kvec.use) {  # ,,,,,,,,,,,,,,,,,,,,

    bix.jk <- X.vlm[, kay]  # n-vector for xij

    object@predictors <- etamat0 + matrix(
       byrow = TRUE, nrow = n.LM, ncol = M,
 bix.jk * vec.step[dirrval] * hdiff.use[kay])  # @@
  
    object@fitted.values <- as.matrix(
    object@family@linkinv(object@predictors,
        extra = object@extra))

    wz.fb <- weights(object, type = "work",
                     ignore.slot = T)
    dimnames(wz.fb) <- NULL  # Faster!!
    U.fb <- if (M == 1) { # 1 x (n.VLM * M)
      wz.fb <- sqrt(wz.fb)
      dim(wz.fb) <- c(1, n.LM * M)
      wz.fb
    } else vchol(wz.fb, M, n.LM)
    cp.X.vlm <- if (M == 1) {  # dim(X.vlm)
      c(U.fb) * X.vlm
    } else mux111(U.fb, X.vlm, M)

    qrR <- qr(cp.X.vlm)  # Faster; not lm.fit()
    R <- qrR$qr[1:p.VLM, 1:p.VLM, drop = F]
    R[lower.tri(R)] <- 0
    if (p.VLM < max(dim(R)))
      stop("'R' is rank deficient")
    attributes(R) <-  # Might be unnecessary
      list(dim = c(p.VLM, p.VLM),
           dimnames = list(nmcobj, nmcobj),
           rank = p.VLM)
    covun <- chol2inv(R)
    SE1.fb <- sqrt(covun[kay, kay])
    ch.d <- as.character(dirrval)  # Safer
    mat.Stdr.tmp[kay, ch.d] <- SE1.fb
    mat.coef.tmp[kay, ch.d] <- cobj[kay] +
      vec.step[dirrval] * hdiff.use[kay]  # @@
  }  # for (kay in kvec.use)  # ,,,,,,,,,,,,,,,,

  list(mat.coef.tmp = mat.coef.tmp,
       mat.Stdr.tmp = mat.Stdr.tmp)
}  # doone


  lambdas <- rep_len(lambdas, 1 + (maxderiv - 1))
  names(lambdas) <- as.character(seq(lambdas) - 1)

  if (!is.Numeric(hdiff, positive = TRUE,
                  length.arg = 1))
    stop("bad input for argument 'hdiff'")
  if (hdiff > 0.5) warning("'hdiff' too large?")

  if (!isFALSE(with1) && !isTRUE(with1))
    stop("'with1' must be TRUE or FALSE")
  if (!isFALSE(use.hdeff) && !isTRUE(use.hdeff))
    stop("'use.hdeff' must be TRUE or FALSE")
  if (use.hdeff)
    warning("'use.hdeff' is not yet implemented")

  if (is.null(subset)) subset <- TRUE

  if (use.hdeff && isFALSE(infos.list$hadof))
    stop("Cannot apply the HDE")

  M <- npred(object)  # Constraints span across ys
  X.vlm <- if (M == 1 && length(object@x))
    object@x else  # May have dimnames
    model.matrixvlm(object,
          type = if (M == 1) "lm" else "vlm",
          label.it = FALSE)  # zz doesnt work?
  dimnames(X.vlm) <- NULL  # Faster!!
  etamat0 <- object@predictors  # yettodo: offsets
  n.LM <- NROW(etamat0)
  cobj <- object@coefficients  # coef(object)
  nmcobj <- names(cobj)
  p.VLM <- length(cobj)
  hdiff.use <- if (fixed.hdiff)
    rep(hdiff, length = p.VLM) else {
    abs(cobj) * hdiff
  }
  hdiff.use[abs(hdiff.use) < 1e-10] <- hdiff
  p.VLM <- length(cobj)
  p.VLM <- length(cobj)
  if (length(theta0) > p.VLM)
    warning("Truncating theta0")
  theta0 <- rep_len(theta0, p.VLM)
  SE1 <- sqrt(diag(chol2inv(object@R)))


  mat.coef.deriv <-
  mat.Stdr.deriv <-
  mat.Wald.deriv <-
    matrix(NA, p.VLM, 1 + maxderiv, dimnames =
           list(nmcobj,
                as.character(0:maxderiv)))
  mat.coef.deriv[, "0"] <- cobj
  mat.Stdr.deriv[, "0"] <- SE1
  mat.Wald.deriv[, "0"] <- (cobj - theta0) / SE1
  mat.coef.tmp <- mat.coef.deriv  # Temporary
  mat.Stdr.tmp <- mat.Stdr.deriv
  vec.step <- head(rep(1:9, each = 2) *
                   c(1, -1), maxderiv)
  names(vec.step) <- as.character(seq(vec.step))
  kvec.use <- 1:p.VLM
  names(kvec.use) <- nmcobj  # For char subset
  if (length(subset))
    kvec.use <- kvec.use[subset]



  upsvec2 <- rep(NA_real_, length(cobj))
  names(upsvec2) <- nmcobj

  TFmat <- NULL
  for (ddd in 0:(maxderiv - 1)) {  # ++++++++++++

  doone.ans <-
    doone(dirrval = ddd + 1,  # \in 1:maxderiv
          mat.coef.tmp = mat.coef.tmp,
          mat.Stdr.tmp = mat.Stdr.tmp)
    mat.coef.tmp <- doone.ans$mat.coef.tmp
    mat.Stdr.tmp <- doone.ans$mat.Stdr.tmp

    mat.Wald.tmp <- (mat.coef.tmp - theta0) / (
                     mat.Stdr.tmp)
    mat.Wald.deriv <-
        fdderiv(which.d = ddd + 1,  # Crucial
                mat.Wald.deriv,
                mat.Wald.tmp,
                use.hdiff = hdiff.use)
    dddp0 <- as.character(ddd)
    dddp1 <- as.character(ddd + 1)
    TFmat <- cbind(TFmat,
                (((-1)^((ddd + 1) * (
                mat.Wald.deriv[, "0"] > 0))) *
                mat.Wald.deriv[, dddp0]) < 0)
    new.indd <- apply(TFmat, 1, all) & 
                (((-1)^(ddd * (
                mat.Wald.deriv[, "0"] > 0))) *
                mat.Wald.deriv[, dddp1]) > 0
    new.indd.use <- new.indd[subset]


   tmp2 <- abs(mat.Wald.deriv[new.indd.use, dddp0])
   tmp3 <- abs(mat.Wald.deriv[new.indd.use, dddp1])
    upsvec2[new.indd.use] <- ddd + tmp2 / (
      tmp2 + lambdas[dddp0] * tmp3)

    if (all(!apply(TFmat, 1, all))) break;
    if (ddd == maxderiv - 1)
      warning("Right-censored WSDM returned. ",
              "Increase 'maxderiv'?")
  }  # for ddd  # ++++++++++++++++



  if (retry) {
    seems.okay <- TRUE
    alist <- vector("list", 3)
    alist[[1]] <- upsvec2[kvec.use]  # kay == 1
    for (kay in 2:3) {  # Compute WSDM thrice
      hdiff.o <- if (kay == 2) hdiff * Mux.div else
                 hdiff / Mux.div
      alist[[kay]] <- ans5 <-
        Recall(object, hdiff = hdiff.o,
               mux.hdiff = 1,  # Adjusted already
               maxderiv = maxderiv,
               theta0 = theta0,
               use.hdeff = use.hdeff,
               doffset = doffset,
               subset = subset,
               derivs.out = FALSE,   # derivs.out
               fixed.hdiff = fixed.hdiff,
               retry = FALSE,  # Nonrecursive!!
               eps.wsdm = eps.wsdm,
               Mux.div = Mux.div,
               warn.retry = FALSE, ...)
      if (any(abs(ans5 - upsvec2[kvec.use]) > eps.wsdm,
              na.rm = TRUE)) {  # Discordant
        if (warn.retry)
        warning("another solution quite diffe",
        "rent... best to try another 'hdiff' ",
        "value; returning the original solution")
        seems.okay <- FALSE
        break
      }
    }  # kay
  }  # retry
      

  if (any(is.na(upsvec2[kvec.use])))
    warning("Some NAs are returned")
  ans8 <- upsvec2[kvec.use]  # No lost attributes
  attr(ans8, "seems.okay") <- seems.okay


  if (!with1) {
    H1mat <- constraints(object)[["(Intercept)"]]
    if (!length(H1mat) || !is.matrix(H1mat))
      stop("'object' has no intercepts!")
    ans8 <- ans8[-seq(ncol(H1mat))]
    kvec.use <- setdiff(kvec.use, 1:ncol(H1mat))
  }


  if (derivs.out)
    list(WSDM   = ans8,
         derivs = mat.Wald.deriv[kvec.use, ]) else
    ans8
}  # wsdm




































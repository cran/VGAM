# These functions are
# Copyright (C) 1998-2023 T.W. Yee, University of Auckland.
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
    stuff@coef3 <- stuff@coef3[, -4]  # Delete the pvalues column









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



  answer
}  # summary.vglm








setMethod("summaryvglmS4VGAM",  signature(VGAMff = "cumulative"),
  function(object,
           VGAMff,
           ...) {
   object@post <-
     callNextMethod(VGAMff = VGAMff,
                    object = object,
                    ...)
  object@post$reverse <- object@misc$reverse


  cfit <- coef(object, matrix = TRUE)
  M <- ncol(cfit)
  if (rownames(cfit)[1] ==  "(Intercept)")
    object@post$expcoeffs <- exp(coef(object)[-(1:M)])


  object@post
})



setMethod("showsummaryvglmS4VGAM",  signature(VGAMff = "cumulative"),
  function(object,
           VGAMff,
           ...) {

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






setMethod("summaryvglmS4VGAM",  signature(VGAMff = "multinomial"),
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



setMethod("showsummaryvglmS4VGAM",  signature(VGAMff = "multinomial"),
  function(object,
           VGAMff,
           ...) {
  cat("\nReference group is level ", object@post$refLevel,
      " of the response\n")
  callNextMethod(VGAMff = VGAMff,
                 object = object,
                 ...)
})



setMethod("summaryvglmS4VGAM",  signature(VGAMff = "VGAMcategorical"),
  function(object,
           VGAMff,
           ...) {
  object@post
})


setMethod("showsummaryvglmS4VGAM", signature(VGAMff = "VGAMcategorical"),
  function(object,
           VGAMff,
           ...) {
})









setMethod("logLik",  "summary.vglm", function(object, ...)
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

  digits <- if (is.null(digits)) options()$digits - 2 else digits

  cat("\nCall:\n", paste(deparse(x@call), sep = "\n", collapse = "\n"),
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
  if (presid &&
      length(Presid) &&
      all(!is.na(Presid)) &&
      is.finite(rdf)) {

    if (rdf/M > 5) {
      rq <- apply(as.matrix(Presid), 2, quantile)  # 5 x M
      dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", "Max"),
                           x@misc$predictors.names)
      cat("\nPearson residuals:\n")
      print(t(rq), digits = digits)
    } else
    if (rdf > 0) {
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
    cat(if (top.half.only) "\nParametric coefficients:" else
        "\nLikelihood ratio test coefficients:", "\n")
    printCoefmat(x@coef4lrt0, digits = digits,
                 signif.stars = use.signif.stars,
                 na.print = "NA",
                 P.values = TRUE, has.Pvalue = TRUE,
                 signif.legend = sum(how.many[-1]) == 0)  # Last one
    Situation <- 3
  }



  if (length(x@coef4score0)) {  # && wald0.arg
    cat(if (top.half.only) "\nParametric coefficients:" else
        "\nRao score test coefficients:", "\n")
    printCoefmat(x@coef4score0, digits = digits,
                 signif.stars = use.signif.stars,
                 na.print = "NA",
                 signif.legend = sum(how.many[3]) == 0)  # Last one
    Situation <- 4
  }



  if (length(x@coef4wald0)) {  # && wald0.arg
    cat(if (top.half.only) "\nParametric coefficients:" else
        "\nWald (modified by IRLS iterations) coefficients:", "\n")
    printCoefmat(x@coef4wald0, digits = digits,
                 signif.stars = use.signif.stars,
                 na.print = "NA")
    Situation <- 1
  } else
  if (length(coef3) && Situation < 0) {
    cat(if (top.half.only) "\nParametric coefficients:" else
        "\nCoefficients:", "\n")
    printCoefmat(coef3, digits = digits,
                 signif.stars = use.signif.stars,
                 na.print = "NA")
    Situation <- 2
  }  # length(coef3)







  if (top.half.only)
    return(invisible(NULL))



  if (M >= 5)
  cat("\nNumber of linear predictors: ", M, "\n")

  if (!is.null(x@misc$predictors.names) && !use.nopredictors) {
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
    cat(paste("\n", prose, "Dispersion Parameter for ",
              x@family@vfamily[1],
              " family:   ", yformat(x@dispersion, digits), "\n",
              sep = ""))
  }


  if (length(deviance(x))) {
    cat("\nResidual deviance:", yformat(deviance(x), digits))
    if (is.finite(rdf))
      cat(" on", round(rdf, digits), "degrees of freedom\n") else
      cat("\n")
  }


  if (length(vll <- logLik.vlm(x))) {
    cat("\nLog-likelihood:", yformat(vll, digits))
    if (is.finite(rdf))
      cat(" on", round(rdf, digits), "degrees of freedom\n") else
      cat("\n")
  }


  if (length(x@criterion)) {
    ncrit <- names(x@criterion)
    for (ii in ncrit)
      if (ii != "loglikelihood" && ii != "deviance")
        cat(paste(ii, ":", sep = ""),
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
      cat("Warning: Hauck-Donner effect detected in the",
            "following estimate(s):\n")
      cat(paste("'", rownames(hado)[hado[, "deriv1"] < 0],
                "'", collapse = ", ", sep = ""))
      cat("\n\n")
    }
  }  # Situation == 2 && length(hado)




  try.this <- findFirstMethod("showsummaryvglmS4VGAM", x@family@vfamily)
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
            dispersion = NULL, untransform = FALSE,
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
    warning("MLE regularity conditions were violated ",
            "at the final iteration of the fitted object")

  if (!untransform)
    return(answer)





  new.way <- TRUE



  if (!is.logical(object@misc$intercept.only))
    stop("cannot determine whether the object is",
         "an intercept-only fit, i.e., 'y ~ 1' is the response")
  if (!object@misc$intercept.only)
    stop("object must be an intercept-only fit, i.e., ",
         "y ~ 1 is the response")

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
    stop("the 'earg' component of 'object@misc' must be a list")




  if (length(LINK) != M &&
      length(LINK) != 1)
    stop("cannot obtain the link functions ",
         "to untransform 'object'")



  if (!is.character(LINK))
    stop("the 'link' component of 'object@misc' should ",
         "be a character vector")

  learg <- length(EARG)
  llink <- length(LINK)
  if (llink != learg)
    stop("the 'earg' component of 'object@misc' should ",
         "be a list of length ", learg)


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
      stop("link functions handled in the new way now")

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








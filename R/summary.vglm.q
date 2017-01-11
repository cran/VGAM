# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.













yformat <- function(x, digits = options()$digits) {
  format(ifelse(abs(x) < 0.001, signif(x, digits), round(x, digits)))
}






summaryvglm <-
  function(object, correlation = FALSE,
           dispersion = NULL, digits = NULL,
           presid = TRUE,
           signif.stars = getOption("show.signif.stars"),
           nopredictors = FALSE,
           ...  # Added 20151211
          ) {







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
                      dispersion = dispersion)





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






  try.this <- findFirstMethod("summaryvglmS4VGAM", object@family@vfamily)
  if (length(try.this)) {
    new.postslot <-
    summaryvglmS4VGAM(object = object,
                      VGAMff = new(try.this),
                      ...)
    answer@post <- new.postslot
  } else {
  }



  answer
}








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


setMethod("showsummaryvglmS4VGAM",  signature(VGAMff = "VGAMcategorical"),
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
           presid = TRUE,
           signif.stars = NULL,   # Use this if logical; 20140728
           nopredictors = NULL,   # Use this if logical; 20150831
           top.half.only = FALSE,  # Added 20160803
           ...  # Added 20151214
           ) {



  M <- x@misc$M
  coef <- x@coef3  # icients
  correl <- x@correlation

  digits <- if (is.null(digits)) options()$digits - 2 else digits

  cat("\nCall:\n", paste(deparse(x@call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")



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



  if (length(coef)) {
    cat(if (top.half.only) "\nParametric coefficients:" else
        "\nCoefficients:", "\n")
    printCoefmat(coef, digits = digits,
                 signif.stars = use.signif.stars,
                 na.print = "NA")
  }

  if (top.half.only)
    return(invisible(NULL))



  cat("\nNumber of linear predictors: ", M, "\n")

  if (!is.null(x@misc$predictors.names) && !use.nopredictors) {
    if (M == 1) {
      cat("\nName of linear predictor:",
          paste(x@misc$predictors.names, collapse = ", "), "\n")
    } else
    if (M <= 5) {
      cat("\nNames of linear predictors:",
        paste(x@misc$predictors.names, collapse = ", "), fill = TRUE)
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
        cat(paste(ii, ":", sep = ""), yformat(x@criterion[[ii]], digits),
            "\n")
  }


  cat("\nNumber of iterations:", format(trunc(x@iter)), "\n")


  if (!is.null(correl)) {
    ncol.X.vlm <- dim(correl)[2]
    if (ncol.X.vlm > 1) {
      cat("\nCorrelation of Coefficients:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1,  -ncol.X.vlm, drop = FALSE], quote = FALSE,
            digits = digits)
    }
  }





  try.this <- findFirstMethod("showsummaryvglmS4VGAM", x@family@vfamily)
  if (length(try.this)) {
    showsummaryvglmS4VGAM(object = x,
            VGAMff = new(try.this),
            ...)
  } else {
  }



  invisible(NULL)
}




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
}


 vcovvlm <-
function(object, dispersion = NULL, untransform = FALSE) {



  so <- summaryvlm(object, correlation = FALSE,
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
  etavector <- predict(object)[1, ]  # Contains transformed parameters
  LINK <- object@misc$link
  EARG <- object@misc$earg  # This could be a NULL
  if (is.null(EARG))
    EARG <- list(theta = NULL)
  if (!is.list(EARG))
    stop("the 'earg' component of 'object@misc' must be a list")




  if (length(LINK) != M &&
      length(LINK) != 1)
    stop("cannot obtain the link functions to untransform 'object'")



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
              c("theta", "inverse", "deriv", "short", "tag"))) > 3
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
}







setMethod("vcov", "vlm",
         function(object, ...)
         vcovvlm(object, ...))


setMethod("vcov", "vglm",
         function(object, ...)
         vcovvlm(object, ...))






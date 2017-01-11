# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.










endfpvgam <- function(object,
                      nonlinear.edf = TRUE,
                      diag.all = FALSE,
                      return.endf = TRUE, ...) {





  M <- npred(object)
  n <- nobs(object, type = "lm")
  wz <- weights(object, type = "working")
  X.vlm.save <- model.matrix(object, type = "vlm")
  U <- vchol(wz, M = M, n = n)
  X.vlm <- mux111(U, X.vlm.save, M = M)
  X.vlm.aug <- rbind(X.vlm,
                 model.matrix(object, type = "penalty"))


  poststuff <-
    mgcv::magic.post.proc(X.vlm.aug,
                          object = object@ospsslot$magicfit, w = NULL)


  if (!return.endf)
    return(poststuff)


  which.X.sm.osps <- object@ospsslot$sm.osps.list$which.X.sm.osps
  all.ncol.Hk <- unlist(lapply(constraints(object, type = "term"), ncol))
  names.which.X.sm.osps <- names(which.X.sm.osps)
  endf <- rep_len(NA_real_, sum(all.ncol.Hk[names.which.X.sm.osps]))
  names(endf) <- vlabel(names.which.X.sm.osps,
                        all.ncol.Hk[names.which.X.sm.osps],
                        M = npred(object))
  use.index <- NULL




  qr1 <- qr(X.vlm.aug)
  qr2 <- qr(X.vlm)
  endf.all <-  diag(solve(crossprod(qr.R(qr1)), crossprod(qr.R(qr2))))
  if (diag.all)
    return(endf.all)




  startstop <- startstoppvgam(object)
  for (iterm in 1:length(startstop)) {
    endf[iterm] <- sum(endf.all[(startstop[[iterm]])])
  }
  endf[endf < 1] <- 1  # Cannot be smoother than linear


  if (nonlinear.edf) endf - 1 else endf
}  # endfpvgam()





show.pvgam <- function(object) {

  digits <- 3


  if (!is.null(cl <- object@call)) {
    cat("\nCall:\n", paste(deparse(cl), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  }



  magicfit <- object@ospsslot$magicfit




  if (FALSE) {
  XX <- model.matrix(object, type = "vlm")
  poststuff <-
    mgcv::magic.post.proc(XX,
                          object = object@ospsslot$magicfit, w = NULL)
  }



  if (FALSE) {
    edf <- rep_len(NA_real_, n.smooth)
    cat("\nEstimated degrees of freedom:\n")
    for (i in 1:n.smooth)
      edf[i] <- sum(x$edf[x$smooth[[i]]$first.para:x$smooth[[i]]$last.para])
    edf.str <- format(round(edf, digits = 4), digits = 3, scientific = FALSE)
    for (i in 1:n.smooth) {
      cat(edf.str[i], " ", sep = "")
      if (i%%7 == 0)
        cat("\n")
    }
    cat(" total =", round(sum(poststuff$edf), digits = 2), "\n")
  }


  endf <- endfpvgam(object)
  cat("\nEstimated nonlinear degrees of freedom:\n")  # based on endfpvgam()
  print(round(endf, digits = digits + 2), digits = digits, scientific = FALSE)

  if (length(endf) > 1)
  cat("Total:",
      format(sum(endf), digits = digits), "\n")

  object@post$endf <- endf  # Good to save this on the object


  if (FALSE)
  cat("\nEstimated degrees of freedom based on poststuff:",
      format(poststuff$edf, digits = digits),
      "\nTotal:",
      format(round(sum(poststuff$edf), digits = digits)), "\n")


  cat("\nUBRE score:", format(magicfit$score, digits = digits + 1), "\n\n")


  if (length(deviance(object)))
    cat("Residual deviance:", format(deviance(object)), "\n")


  llx <- logLik.vlm(object = object)
  if (length(llx))
    cat("Log-likelihood:", format(llx), "\n")




  invisible(object)
}



setMethod("show", "pvgam", function(object) show.pvgam(object))







if (!isGeneric("endf"))
    setGeneric("endf", function(object, ...)
    standardGeneric("endf"))


setMethod("endf", "pvgam", function(object, ...)
          endfpvgam(object, ...))

setMethod("endf", "summary.pvgam", function(object, ...)
          endfpvgam(object, ...))








show.vglm <- function(object) {

  if (!is.null(cl <- object@call)) {
    cat("\nCall:\n", paste(deparse(cl), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  }

  coef <- object@coefficients
  if (any(nas <- is.na(coef))) {
    if (is.null(names(coef)))
      names(coef) <- paste("b", seq_along(coef), sep = "")
    cat("\nCoefficients: (", sum(nas),
        " not defined because of singularities)\n", sep = "")
  } else {
    cat("\nCoefficients:\n")
  }
  print(coef)

  rank <- object@rank
  if (!length(rank))
    rank <- sum(!nas)
  nobs <- if (length(object@df.total)) object@df.total else
          length(object@residuals)
  rdf <- object@df.residual
  if (!length(rdf))
    rdf <- nobs - rank
  cat("\nDegrees of Freedom:", nobs, "Total;", rdf, "Residual\n")

  if (length(deviance(object)))
    cat("Residual deviance:", format(deviance(object)), "\n")
  llx <- logLik.vlm(object = object)

  if (length(llx))
    cat("Log-likelihood:", format(llx), "\n")

  if (length(object@criterion)) {
    ncrit <- names(object@criterion)
    for (ii in ncrit)
      if (ii != "loglikelihood" &&
          ii != "deviance")
          cat(paste(ii, ":", sep = ""),
              format(object@criterion[[ii]]), "\n")
  }





  try.this <- findFirstMethod("showvglmS4VGAM", object@family@vfamily)
  if (length(try.this)) {
    showvglmS4VGAM(object = object,
                   VGAMff = new(try.this))
  } else {
  }




  invisible(object)
}






show.vgam <- function(object) {

  digits <- 2


  if (!is.null(cl <- object@call)) {
    cat("\nCall:\n", paste(deparse(cl), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  }

  coef <- object@coefficients
  nas <- is.na(coef)

  rank <- object@rank
  if (is.null(rank))
      rank <- sum(!nas)
  nobs <- if (length(object@df.total)) object@df.total else
          length(object@residuals)
  rdf <- object@df.residual
  if (is.null(rdf))
    rdf <- nobs - rank
  cat("\nDegrees of Freedom:", nobs, "Total;",
      format(round(rdf, digits = digits)), "Residual\n")

  if (length(deviance(object)))
    cat("Residual deviance:", format(deviance(object)), "\n")

  llx <- logLik.vlm(object = object)

  if (length(llx))
    cat("Log-likelihood:", format(llx), "\n")

  criterion <- attr(terms(object), "criterion")
  if (!is.null(criterion) &&
      criterion != "coefficients")
    cat(paste(criterion, ":", sep = ""),
        format(object[[criterion]]), "\n")




  try.this <- findFirstMethod("showvgamS4VGAM", object@family@vfamily)
  if (length(try.this)) {
    showvgamS4VGAM(object = object,
                   VGAMff = new(try.this))
  } else {
  }



  invisible(object)
}




setMethod("show",  "vlm", function(object) show.vlm (object))
setMethod("show", "vglm", function(object) show.vglm(object))
setMethod("show", "vgam", function(object) show.vgam(object))










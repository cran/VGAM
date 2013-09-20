# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.







show.vglm <- function(object) {
  if (!is.null(cl <- object@call)) {
    cat("Call:\n")
    dput(cl)
  }

  coef <- object@coefficients
  if (any(nas <- is.na(coef))) {
    if (is.null(names(coef)))
      names(coef) <- paste("b", 1:length(coef), sep = "")  
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

  invisible(object)
}










show.vgam <- function(object) {

  digits <- 2


  if (!is.null(cl <- object@call)) {
    cat("Call:\n")
    dput(cl)
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

  invisible(object)
}




setMethod("show",  "vlm", function(object) show.vlm (object))
setMethod("show", "vglm", function(object) show.vglm(object))
setMethod("show", "vgam", function(object) show.vgam(object))








 if (FALSE)
print.vglm <- function(x, ...) {
  if (!is.null(cl <- x@call)) {
    cat("Call:\n")
    dput(cl)
  }

  coef <- x@coefficients
  if (any(nas <- is.na(coef))) {
    if (is.null(names(coef)))
      names(coef) <- paste("b", 1:length(coef), sep = "")  
    cat("\nCoefficients: (", sum(nas), 
        " not defined because of singularities)\n", sep = "")
  } else {
    cat("\nCoefficients:\n")
  }
  print.default(coef, ...)

  rank <- x@rank
  if (!length(rank))
    rank <- sum(!nas)
  nobs <- if (length(x@df.total)) x@df.total else
          length(x@residuals)
  rdf <- x@df.residual
  if (!length(rdf))
    rdf <- nobs - rank
  cat("\nDegrees of Freedom:", nobs, "Total;", rdf, "Residual\n")

  if (length(deviance(x)))
    cat("Residual deviance:", format(deviance(x)), "\n")
  llx <- logLik.vlm(object = x)

  if (length(llx))
    cat("Log-likelihood:", format(llx), "\n")

  if (length(x@criterion)) {
    ncrit <- names(x@criterion)
    for (ii in ncrit)
      if (ii != "loglikelihood" && ii != "deviance")
          cat(paste(ii, ":", sep = ""),
              format(x@criterion[[ii]]), "\n")
  }

  invisible(x)
}



 if (FALSE)
print.vgam <- function(x, digits = 2, ...) {

  if (!is.null(cl <- x@call)) {
    cat("Call:\n")
    dput(cl)
  }

  coef <- x@coefficients
  nas <- is.na(coef)

  rank <- x@rank
  if (is.null(rank))
      rank <- sum(!nas)
  nobs <- if (length(x@df.total)) x@df.total else
          length(x@residuals)
  rdf <- x@df.residual
  if (is.null(rdf))
    rdf <- nobs - rank
  cat("\nDegrees of Freedom:", nobs, "Total;",
      format(round(rdf, dig = digits)), "Residual\n")

  if (length(deviance(x)))
    cat("Residual deviance:", format(deviance(x)), "\n")

  llx <- logLik.vlm(object = x)

  if (length(llx))
    cat("Log-likelihood:", format(llx), "\n")

  criterion <- attr(terms(x), "criterion")  # 11/8/03; x@terms$terms,
  if (!is.null(criterion) &&
      criterion != "coefficients")
    cat(paste(criterion, ":", sep = ""), format(x[[criterion]]), "\n")

  invisible(x)
}



 if (FALSE) {

setMethod("print",  "vlm", function(x, ...)  print.vlm(x, ...))
setMethod("print", "vglm", function(x, ...) print.vglm(x, ...))
setMethod("print", "vgam", function(x, ...) print.vgam(x, ...))



    setMethod("show",  "vlm", function(object)  print.vlm(object))
    setMethod("show", "vglm", function(object) print.vglm(object))
    setMethod("show", "vgam", function(object) print.vgam(object))

}




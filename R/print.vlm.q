# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.








show.vlm <- function(object) {
  if (!is.null(cl <- object@call)) {
    cat("Call:\n")
    dput(cl)
  }

  coef <- object@coefficients
  cat("\nCoefficients:\n")
  print(coef)

  rank <- object@rank
  if (is.null(rank))
    rank <- sum(!is.na(coef))
  n <- object@misc$n
  M <- object@misc$M
  nobs <- if (length(object@df.total)) object@df.total else n * M
  rdf <- object@df.residual
  if (is.null(rdf))
    rdf <- (n - rank) * M
  cat("\nDegrees of Freedom:", nobs, "Total;",
      rdf, "Residual\n")

  if (length(deviance(object)) &&
      is.finite(deviance(object)))
    cat("Deviance:", format(deviance(object)), "\n")
  if (length(object@ResSS) &&
      is.finite(object@ResSS))
    cat("Residual Sum of Squares:", format(object@ResSS), "\n")

  invisible(object)
}



setMethod("show", "vlm",
    function(object)
    show.vlm(object))







if (FALSE)
print.vlm <- function(x, ...) {
  if (!is.null(cl <- x@call)) {
    cat("Call:\n")
    dput(cl)
  }

  coef <- x@coefficients
  cat("\nCoefficients:\n")
  print(coef, ...)

  rank <- x@rank
  if (is.null(rank))
    rank <- sum(!is.na(coef))
  n <- x@misc$n
  M <- x@misc$M
  nobs <- if (length(x@df.total)) x@df.total else n * M
  rdf <- x@df.residual
  if (is.null(rdf))
    rdf <- (n - rank) * M
  cat("\nDegrees of Freedom:", nobs, "Total;",
      rdf, "Residual\n")

  if (length(deviance(x)) &&
      is.finite(deviance(x)))
    cat("Deviance:", format(deviance(x)), "\n")
  if (length(x@ResSS) &&
      is.finite(x@ResSS))
    cat("Residual Sum of Squares:", format(x@ResSS), "\n")

  invisible(x)
}




if (!is.R()) {
setMethod("show", "vlm",
    function(object)
    print.vlm(object))
}



if (FALSE)
setMethod("print", "vlm",
    function(x, ...)
    print.vlm(x, ...))




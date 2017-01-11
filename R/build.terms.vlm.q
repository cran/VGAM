# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.





if (!isGeneric("terms"))
  setGeneric("terms", function(x, ...) standardGeneric("terms"))



terms.vlm <- function(x, ...) {
  termsvlm(x, ...)
}



termsvlm <- function(x, form.number = 1, ...) {
  if (!is.Numeric(form.number, integer.valued = TRUE,
                  length.arg = 1, positive = TRUE) ||
      form.number > 2)
    stop("argument 'form.number' must be 1 or 2")


  v <- if (form.number == 1) {
    v <- x@terms
    if (!length(v))
      stop("terms slot is empty")
    v$terms
  } else if (form.number == 2) {
    x@misc$Terms2
  }
  if (length(v)) {
    v
  } else {
    warning("no terms component; returning a NULL")
    NULL
  }
}





setMethod("terms", "vlm", function(x, ...) terms.vlm(x, ...))





Build.terms.vlm <-
  function(x, coefs, cov = NULL, assign, collapse = TRUE, M,
           dimname = NULL, coefmat = NULL) {


  cov.true <- !is.null(cov)
  if (collapse) {
    fit <- matrix(x %*% coefs, ncol = M, byrow = TRUE)
    dimnames(fit) <- dimname
    if (M == 1)
      fit <- c(fit)
    if (cov.true) {
      var <- rowSums((x %*% cov) * x)
      list(fitted.values = fit,
           se.fit = if (M == 1) c(sqrt(var)) else
                    matrix(sqrt(var), ncol = M,
                           byrow = TRUE, dimnames = dimname))
    } else {
      fit
    }
  } else {

    constant <- attr(x, "constant")
    if (!is.null(constant)) {
      constant <- as.vector(t(coefmat) %*% constant)
    }

    if (missing(assign))
      assign <- attr(x, "assign")
    if (is.null(assign))
      stop("Need an 'assign' list")
    fit <- array(0, c(nrow(x), length(assign)),
                 list(dimnames(x)[[1]], names(assign)))
    if (cov.true)
      se <- fit
    TL <- sapply(assign, length)
    simple <- (TL == 1)
    complex <- (TL > 1)
    if (any(simple)) {
      asss <- unlist(assign[simple])
      ones <- rep_len(1, nrow(x))
      fit[, simple] <- x[, asss] * outer(ones, coefs[asss])
      if (cov.true)
        se[, simple] <- abs(x[, asss]) * outer(ones, sqrt(diag(cov))[asss])
    }
    if (any(complex)) {
      assign <- assign[complex]
      for (term in names(assign)) {
        TT <- assign[[term]]
        xt <- x[, TT]
        fit[, term] <- xt %*% coefs[TT]
        if (cov.true) {
          se[, term] <- sqrt(rowSums((xt %*% cov[TT, TT]) * xt))
        }
      }
    }
    attr(fit, "constant") <- constant

    if (cov.true)
      list(fitted.values = fit,
           se.fit        = se) else
      fit
  }
}  # Build.terms.vlm()




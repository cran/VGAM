# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.









formula.vlm <- function(x, ...)
  formulavlm(x, ...)








formulavlm <- function(x, form.number = 1, ...) {
  if (!is.Numeric(form.number, integer.valued = TRUE,
                  length.arg = 1, positive = TRUE) ||
      form.number > 2)
    stop("argument 'form.number' must be 1 or 2")

  if (!any(slotNames(x) == "misc"))
    stop("cannot find slot 'misc'")

  if (form.number == 1) x@misc$formula else x@misc$form2
}



formulaNA.VGAM <- function(x, ...) {
  stop("a formula does not make sense for object 'x'")
}





setMethod("formula", "vlm",
          function(x, ...)
              formulavlm(x = x, ...))

setMethod("formula", "vglm",
          function(x, ...)
              formulavlm(x = x, ...))



setMethod("formula", "vgam",
          function(x, ...)
              formulavlm(x = x, ...))

setMethod("formula", "rrvglm",
          function(x, ...)
              formulavlm(x = x, ...))

setMethod("formula", "qrrvglm",
          function(x, ...)
              formulavlm(x = x, ...))

setMethod("formula", "grc",
          function(x, ...)
              formulavlm(x = x, ...))










variable.namesvlm <- function(object, full = FALSE, ...) {
  qrslot <- object@qr
  if (!length(qrslot$qr)) {
    use.this <- object@x
    if (!length(use.this))
      stop("argument 'object' has empty 'qr' and 'x' slots.")
  } else {
    use.this <- qrslot$qr
  }
  if (full) dimnames(use.this)[[2]] else
  if (object@rank) dimnames(use.this)[[2]][seq_len(object@rank)] else
  character(0)
}




variable.namesrrvglm <- function(object, ...) {

  qrslot <- object@qr
  if (!length(qrslot$qr)) {
    use.this <- object@x
    if (!length(use.this))
      stop("argument 'object' has empty 'qr' and 'x' slots.")
  } else {
    use.this <- qrslot$qr
  }
  dimnames(use.this)[[2]]
}







case.namesvlm <- function(object, full = FALSE, ...) {
  w <- weights(object, type="prior")
  use.this <- residuals(object, type = "working")
  if (!length(use.this))
    use.this <- object@x
  if (!length(use.this))
    use.this <- object@y
  if (!length(use.this))
    stop("argument 'object' has empty 'x' and 'y' slots.")
  dn <- dimnames(use.this)[[1]]
  if (full || is.null(w) || ncol(cbind(w)) != 1)
    dn else dn[w != 0]
}


setMethod("variable.names", "vlm",
          function(object, ...)
              variable.namesvlm(object = object, ...))

setMethod("variable.names", "vglm",
          function(object, ...)
              variable.namesvlm(object = object, ...))

setMethod("variable.names", "vgam",
          function(object, ...)
              variable.namesvlm(object = object, ...))

setMethod("variable.names", "rrvglm",
          function(object, ...)
              variable.namesrrvglm(object = object, ...))

setMethod("variable.names", "qrrvglm",
          function(object, ...)
              variable.namesvlm(object = object, ...))

setMethod("variable.names", "grc",
          function(object, ...)
              variable.namesvlm(object = object, ...))






setMethod("case.names", "vlm",
          function(object, ...)
              case.namesvlm(object = object, ...))

setMethod("case.names", "vglm",
          function(object, ...)
              case.namesvlm(object = object, ...))

setMethod("case.names", "vgam",
          function(object, ...)
              case.namesvlm(object = object, ...))

setMethod("case.names", "rrvglm",
          function(object, ...)
              case.namesvlm(object = object, ...))

setMethod("case.names", "qrrvglm",
          function(object, ...)
              case.namesvlm(object = object, ...))

setMethod("case.names", "grc",
          function(object, ...)
              case.namesvlm(object = object, ...))








has.interceptvlm <- function(object, form.number = 1, ...) {
  if (!is.Numeric(form.number, integer.valued = TRUE,
                  length.arg = 1, positive = TRUE) ||
      form.number > 2)
    stop("argument 'form.number' must be 1 or 2")


  if (form.number == 1) {
    if (is.numeric(aa <- attr(terms(object), "intercept")))
      as.logical(aa) else
      FALSE
  } else if (form.number == 2) {
    if (is.numeric(aa <- attr(terms(object, form.number = 2), "intercept")))
      as.logical(aa) else
      FALSE
  }
}


if (!isGeneric("has.intercept"))
    setGeneric("has.intercept", function(object, ...)
               standardGeneric("has.intercept"),
               package = "VGAM")


setMethod("has.intercept",  "vlm", function(object, ...)
           has.interceptvlm(object, ...))







term.namesvlm <- function(model, form.number = 1, ...) {
  if (!is.Numeric(form.number, integer.valued = TRUE,
                  length.arg = 1, positive = TRUE) ||
      form.number > 2)
    stop("argument 'form.number' must be 1 or 2")

    aa <- if (has.intercept(model, form.number = form.number))
          "(Intercept)" else NULL
    bb <- attr(terms(model, form.number = form.number), "term.labels")
    c(aa, bb)
}


if (!isGeneric("term.names"))
    setGeneric("term.names", function(model, ...)
               standardGeneric("term.names"),
               package = "VGAM")


setMethod("term.names",  "vlm", function(model, ...)
           term.namesvlm(model, ...))







responseNamevlm <- function(model, form.number = 1, ...) {
  TERMS.MODEL <-terms(model, form.number = form.number)
  if (length(aa <- attr(TERMS.MODEL, "dataClasses")) &&
      length(bb <- attr(TERMS.MODEL, "response"   )) &&
      bb == 1) {
    names(aa)[1]
  } else {
    NULL
  }
}


if (!isGeneric("responseName"))
  setGeneric("responseName", function(model, ...)
             standardGeneric("responseName"),
             package = "VGAM")


setMethod("responseName",  "vlm", function(model, ...)
           responseNamevlm(model, ...))






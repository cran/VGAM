# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.







formulavlm <- function(x, fnumber = 1, ...) {
  if (!is.Numeric(fnumber, integer.valued = TRUE,
                  length.arg = 1, positive = TRUE) ||
      fnumber > 2)
    stop("argument 'fnumber' must be 1 or 2")

  if (!any(slotNames(x) == "misc"))
    stop("cannot find slot 'misc'")

  if (fnumber == 1) x@misc$formula else x@misc$form2
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



















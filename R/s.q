# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.




s <- function(x, df = 4, spar = 0, ...) {

  xs <- substitute(x)
  ans <- as.character(xs)
  if (length(ans) > 1)
    stop("argument 'x' must be of length one")

  call <- deparse(sys.call())

  if (ncol(as.matrix(x)) > 1)
    stop("argument 'x' must be a vector")
  if (!is.null(levels(x))) {
    x <- if (is.ordered(x)) {
      as.vector(x)
    } else
      stop("unordered factors cannot be used as smoothing variables")
  }
  attr(x, "spar") <- spar
  attr(x, "df") <- df
  attr(x, "call") <- call
  attr(x, "class") <- "smooth"
  attr(x, "s.xargument") <- ans   # Needed for prediction and constraints


  a <- is.na(x)
  if (any(a))
    attr(x, "NAs") <- seq(along = x)[a]

  x
}




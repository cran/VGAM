# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.





effects.vlm <- function(object, ...) {
  cat("Sorry, this function has not been written yet. Returning a NULL.\n")
  invisible(NULL)
}


if (!isGeneric("effects"))
  setGeneric("effects", function(object, ...)
             standardGeneric("effects"))


if (is.R()) {
  setMethod("effects",  "vlm", function(object, ...)
            effects.vlm(object, ...))
} else {
  setMethod("effects",  "vlm", function(object, ...)
            effects.vlm(object, ...))
}



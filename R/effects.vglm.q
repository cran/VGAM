# These functions are
# Copyright (C) 1998-2024 T.W. Yee, University of Auckland.
# All rights reserved.








effects.vlm <- function(object, ...) {
  warning("Sorry, this function has not been written yet. ",
          "Returning a NULL.")
  invisible(NULL)
}


if (!isGeneric("effects"))
  setGeneric("effects", function(object, ...)
             standardGeneric("effects"))


  setMethod("effects",  "vlm", function(object, ...)
            effects.vlm(object, ...))






Influence.vglm <-
  function(object, weighted = TRUE, ...) {

  dl.deta <- weights(object, deriv = TRUE, type = "working")$deriv
  if (!is.matrix(dl.deta))
    dl.deta <- cbind(dl.deta)

      
  if (!weighted)
    stop("currently only the weighted version is returned")
    
  X.vlm <- model.matrix(object, type = "vlm")
  nn <- nobs(object)
  p.vlm <- ncol(X.vlm)
  M <- npred(object)
  dl.dbeta.vlm <- matrix(0, nn, p.vlm)
  for (jay in 1:M) {
    vecTF <- rep(FALSE, M)
    vecTF[jay] <- TRUE  # Recycling
    dl.dbeta.vlm <- dl.dbeta.vlm +
        X.vlm[vecTF, , drop = FALSE] * dl.deta[, jay]
  }

  inv.info <- vcov(object)
  inffuns <- dl.dbeta.vlm %*% inv.info
  if (M > 1) {
    rns <- unlist(strsplit(rownames(inffuns), ":"))
    rownames(inffuns) <- rns[c(TRUE, FALSE)]
  }
  inffuns
}  # Influence.vglm





if (!isGeneric("Influence"))
  setGeneric("Influence", function(object, ...)
             standardGeneric("Influence"))


setMethod("Influence",  "vglm", function(object, ...)
          Influence.vglm(object, ...))


setMethod("Influence",  "vgam", function(object, ...)
          stop("This methods function has not been written"))



setMethod("Influence",  "rrvglm", function(object, ...)
          stop("This methods function has not been written"))





























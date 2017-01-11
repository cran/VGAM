# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.







plotvlm <- function(object, residuals = NULL, rugplot= FALSE, ...) {
  stop("sorry, this function hasn't been written yet")
}






plotvglm <-
  function(x,
           which = "(All)",



           ...) {


  show <- rep(FALSE, 10000)
  if (is.character(which) && which == "(All)") {
    show[TRUE] <- TRUE
  } else {
    show[which] <- TRUE
  }



  presid <- resid(x, type = "pearson")
  if (!is.matrix(presid) == 1)
    presid <- as.matrix(presid)
  lapred <- predict(x)
  M <- ncol(lapred)
  for (jay in 1:M) {
    if (show[jay]) {
      use.x <- lapred[, jay]
      if (one.x <- diff(range(use.x)) < 1e-10)
        use.x[TRUE] <- jitter(mean(use.x))
      plot(use.x, presid[, jay],
           ylab = "Pearson residuals",
           xlab = paste(if (one.x) "Jittered l" else
                        "L", "inear predictor ", jay, sep = ""),
           ...)
    }
  }



  hvmat <- hatvalues(x)
  for (jay in 1:M) {
    if (show[M + jay]) {
      use.x <- hvmat[, jay]
      if (one.x <- diff(range(use.x)) < 1e-10)
        use.x[TRUE] <- jitter(mean(use.x))
      plot(use.x, presid[, jay],
           ylab = "Pearson residuals",
           xlab = paste(if (one.x) "Jittered h" else
                        "H", "at values for linear predictor ", jay,
                        sep = ""),
           ...)
    }
  }






  invisible(x)
}












setMethod("plot", "vlm",
           function(x, y, ...) {
           if (!missing(y))
             stop("cannot process the 'y' argument")
           invisible(plotvlm(x, y, ...))})
setMethod("plot", "vglm",
           function(x, y, ...) {
           if (!missing(y))
             stop("cannot process the 'y' argument")
           invisible(plotvglm(x = x, ...))})








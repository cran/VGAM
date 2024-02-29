# These functions are
# Copyright (C) 1998-2024 T.W. Yee, University of Auckland.
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










 spikeplot <-
  function(x,
           freq = FALSE,  # Default are proportions
           as.table = FALSE,
           col = par("col"),  #"pink2",  # Parent distribution
           lty = par("lty"),  #
           lwd = par("lwd"),  #
           lend = par("lend"),  # "butt", "round", etc.
           type = "h",  #
           xlab = deparse1(substitute(x)),  # NULL,
           ylab = NULL,
           capped = FALSE, cex =  sqrt(lwd) / 2, pch = 19,
           pcol = col,
           scol = NULL, slty = NULL, slwd = NULL,
           new.plot = TRUE, offset.x = 0,  # 20211123
           ymux = 1,  # 20211129
           ...) {  # ... allows many graphical params, e.g., xlim

      deparse1 <- function(expr, collapse = " ",
                           width.cutoff = 500L, ...)
    paste(deparse(expr, width.cutoff, ...), collapse = collapse)
  xlabel <- xlab
  ylabel <- if (length(ylab)) ylab else
             ifelse(freq, "Frequency", "Proportion")
  if (!is.numeric(x))
    stop("argument 'x' is not numeric")
  tx <- table(x)  # exclude, useNA
  ntx <- names(tx)
  x.use <- x.use2 <- as.numeric(ntx)
  if (as.table) {
    y.use <- y.use2 <- ymux * (if (freq) tx else tx / sum(tx))
    if (new.plot)
      plot(y.use, col = col, xlab = xlabel, ylab = ylabel,
           type = type, lwd = lwd, lty = lty, lend = lend, ...) else
      points(y.use, col = col,  #  xlab = xlabel, ylab = ylabel,
             type = type, lwd = lwd, lty = lty, lend = lend, ...)
  } else {
      y.use <- ymux * (if (freq) as.vector(tx) else
                       as.vector(tx / sum(tx)))



    specialvals <- NULL  # None 
    if ((length.sargs <- length(scol) + length(slty) + length(slwd))) {
      combo.args <- c(scol, slty, slwd)
      specialvals <- unique(sort(unlist(combo.args)))

 
      ooo <- match(x.use, specialvals)   # combo.vals
      x.use2 <- x.use[is.na(ooo)]
      y.use2 <- y.use[is.na(ooo)]
    } else {
      x.use2 <- x.use
      y.use2 <- y.use
    }

    if (new.plot)
      plot(x.use2 + offset.x, y.use2,
           type = type, xlab = xlabel, ylab = ylabel,
           col = col, lty = lty, lwd = lwd, lend = lend, ...) else
    points(x.use2 + offset.x, y.use2,
           type = type,  # xlab = xlabel, ylab = ylabel,
           col = col, lty = lty, lwd = lwd, lend = lend, ...)

      
    if (length.sargs) {
      if (length(scol)) {
        vec_scol <- unlist(scol, use.names = FALSE)
        rep_scol <- rep(names(scol), times = sapply(scol, length))
        names(vec_scol) <- rep_scol
      }  # length(scol)
      if (length(slty)) {
        vec_slty <- unlist(slty, use.names = FALSE)
        rep_slty <- rep(names(slty), times = sapply(slty, length))
        names(vec_slty) <- rep_slty
      }  # length(slty)
      if (length(slwd)) {
        vec_slwd <- unlist(slwd, use.names = FALSE)
        rep_slwd <- rep(names(slwd), times = sapply(slwd, length))
        names(vec_slwd) <- rep_slwd
      }  # length(slwd)


      for (xx in specialvals) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        use_scol <- col[1]
        use_slty <- lty[1]
        use_slwd <- lwd[1]
        if (length(scol) && any(vec_scol == xx)) {
          use_scol <- names(vec_scol[vec_scol == xx])
        }
        if (length(slty) && any(vec_slty == xx)) {
          use_slty <- names(vec_slty[vec_slty == xx])
          bits <- substring(use_slty, 1:nchar(use_slty),
                            1:nchar(use_slty))
          if (all(bits %in% as.character(0:9)))
           use_slty <- as.numeric(use_slty)
        }
        if (length(slwd) && any(vec_slwd == xx)) {
          use_slwd <- as.numeric(names(vec_slwd[vec_slwd == xx]))
        }

        points(xx + offset.x, y.use[x.use == xx], type = type,
               col = use_scol, lty = use_slty,
               lwd = use_slwd, lend = lend)

        if (capped)
          points(xx + offset.x, y.use[x.use = xx],
                 cex = cex, pch = pch, col = use_scol)
      }  # for ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    }  # length.sargs
  }  # !as.table

  if (capped)
    points(x.use2 + offset.x, y.use2, cex = cex, pch = pch, col = pcol)

  invisible(tx)
}  # spikeplot




















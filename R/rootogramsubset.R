 # This is a subset of rootogram.R from \pkg{countreg},
 # for use by \pkg{VGAM}.
 # The functions have been renamed to rootogram0.











rootogram0 <- function(object, ...) {
  UseMethod("rootogram0")
}








rootogram0.default <-
  function(object, fitted, breaks = NULL,
           style = c("hanging", "standing", "suspended"),
           scale = c("sqrt", "raw"), plot = TRUE,
           width = NULL, xlab = NULL, ylab = NULL,
           main = NULL, lowsup = 0L, ...) {


  ## rectangle style
  scale <- match.arg(scale, c("sqrt", "raw"))[1]
  style <- match.arg(style, c("hanging", "standing", "suspended"))[1]

  ## default annotation
  if (is.null(xlab)) {
    xlab <- if (is.null(names(dimnames(object)))) {
      deparse(substitute(object))
    } else {
      names(dimnames(object))[1L]
    }
  }
  if(is.null(ylab)) {
    ylab <- if(scale == "raw") "Frequency" else "sqrt(Frequency)" 
  }
  if(is.null(main)) main <- deparse(substitute(fitted))
  
  ## breaks, midpoints, widths
  if (is.null(breaks)) {
    x <- as.numeric(names(object)) + lowsup
    if(length(x) < 1L) x <- lowsup:(length(object) - 1L)
    breaks <- (head(x, -1L) + tail(x, -1L))/2
    breaks <- c(2 * head(x, 1L) - head(breaks, 1L), breaks,
      2 * tail(x, 1L) - tail(breaks, 1L))
    if(is.null(width)) width <- 0.9
  } else {
    x <- (head(breaks, -1L) + tail(breaks, -1L)) / 2
    if(is.null(width)) width <- 1
  }

  ## raw vs. sqrt scale
  if (scale == "sqrt") {
    obsrvd <- sqrt(as.vector(object))
    expctd <- sqrt(as.vector(fitted))
  } else {
    obsrvd <- as.vector(object)
    expctd <- as.vector(fitted)
  }

  ## height/position of rectangles
  y <- if(style == "hanging") expctd - obsrvd else 0
  height <- if(style == "suspended") expctd - obsrvd else obsrvd

  ## collect everything as data.frame
    rval <- data.frame(observed = as.vector(object),
                       expected = as.vector(fitted),
    x = x, y = y, width = diff(breaks) * width, height = height,
    line = expctd)
  attr(rval, "style") <- style
  attr(rval, "scale") <- scale
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  class(rval) <- c("rootogram0", "data.frame")
  
  ## also plot by default
  if (plot) plot(rval, ...)
  
  ## return invisibly
  invisible(rval)
}  # rootogram0.default






plot.rootogram0 <- function(x,
  xlim = NULL, ylim = NULL, xlab = NULL,
  ylab = NULL, main = NULL,
  border = "black", fill = "lightgray",
  col = "#B61A51",
  lwd = 2, pch = 19, lty = 1, max = NULL,
  type = NULL, axes = TRUE, ...) {




 

  if(is.null(x$group)) x$group <- 1L
  n <- max(x$group)
  if(is.null(type))
    type <- ifelse(any(table(x$group) > 20L), "l", "b")

  ## annotation
  if(is.null(xlab)) xlab <- TRUE
  if(is.null(ylab)) ylab <- TRUE
  if(is.null(main)) main <- TRUE
  xlab <- rep(xlab, length.out = n)
  ylab <- rep(ylab, length.out = n)
  main <- rep(main, length.out = n)
  if(is.logical(xlab)) xlab <- ifelse(xlab, attr(x, "xlab"), "")
  if(is.logical(ylab)) ylab <- ifelse(ylab, attr(x, "ylab"), "")
  if(is.logical(main)) main <- ifelse(main, attr(x, "main"), "")

  ## plotting function
  rootogram1 <- function(d, ...) {
    ## rect elements
    xleft <- d$x - d$width/2
    xright <- d$x + d$width/2
    ybottom <- d$y
    ytop <- d$y + d$height
    j <- unique(d$group)
    
    ## defaults
    if(is.null(xlim)) xlim <- range(c(xleft, xright))
    if(is.null(ylim)) ylim <- range(c(ybottom, ytop, d$line))

    ## draw rootogram
    plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
         xlab = xlab[j], ylab = ylab[j], main = main[j],
         axes = FALSE, ...)
    if (axes) {
      axis(1, ...)
      axis(2, ...)
    }
    rect(xleft, ybottom, xright, ytop, border = border, col = fill)
    abline(h = 0, col = border)
    lines(d$x, d$line,
      col = col, pch = pch, type = type, lty = lty, lwd = lwd)
   }
   
   ## draw plots
   if(n > 1L) par(mfrow = n2mfrow(n))
   for(i in 1L:n)
     rootogram1(x[x$group == i, ], ...)
}  # plot.rootogram0 













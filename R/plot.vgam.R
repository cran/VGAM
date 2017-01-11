# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.













plotvgam <-
plot.vgam <-
  function(x, newdata = NULL, y = NULL, residuals = NULL, rugplot = TRUE,
           se = FALSE, scale = 0,
           raw = TRUE, offset.arg = 0, deriv.arg = 0, overlay = FALSE,
           type.residuals = c("deviance", "working", "pearson", "response"),
           plot.arg = TRUE, which.term = NULL, which.cf = NULL,
           control = plotvgam.control(...),
           varxij = 1, ...) {

  missing.control <- missing(control)

  na.act <- x@na.action
  x@na.action <- list()

  if (!is.Numeric(varxij, integer.valued = TRUE,
                  length.arg = 1, positive = TRUE))
      stop("bad input for the 'varxij' argument")
  if (any(slotNames(x) == "control")) {
    x@control$varxij <- varxij
  }


  missing.type.residuals <- missing(type.residuals)
  if (mode(type.residuals) != "character" &&
      mode(type.residuals) != "name")
    type.residuals <- as.character(substitute(type.residuals))
  if (!missing.type.residuals)
    type.residuals <- match.arg(type.residuals,
        c("deviance", "working", "pearson", "response"))[1]


  if (!is.Numeric(deriv.arg, integer.valued = TRUE,
                  length.arg = 1) ||
      deriv.arg < 0)
    stop("bad input for the 'deriv' argument")

  if (se && deriv.arg > 0) {
    warning("standard errors not available with derivatives. ",
            "Setting 'se = FALSE'")
    se <- FALSE
  }

  preplot.object <- x@preplot
  if (!length(preplot.object)) {
    preplot.object <- preplotvgam(x, newdata = newdata,
                                  raw = raw,
                                  deriv.arg = deriv.arg, se = se,
                                  varxij = varxij)
  }

  x@preplot <- preplot.object


  if (!is.null(residuals) && length(residuals) == 1) {
    if (residuals) {
      if (missing.type.residuals) {
        for (rtype in type.residuals)
          if (!is.null(residuals <- resid(x, type = rtype))) break
      } else {
       residuals = resid(x, type = type.residuals)
        if (!length(residuals))
          warning("residuals are NULL. Ignoring 'residuals = TRUE'")
      }
    } else {
      residuals <- NULL
    }
  }

  if (!missing.control) {
    control <- c(plotvgam.control( .include.dots = FALSE, ...),
                control, plotvgam.control(...))
  }

  x@post$plotvgam.control <- control  # Add it to the object

  if (plot.arg)
    plotpreplotvgam(preplot.object, residuals = residuals,
                    rugplot = rugplot, scale = scale, se = se,
                    offset.arg = offset.arg,
                    deriv.arg = deriv.arg,
                    overlay = overlay,
                    which.term = which.term, which.cf = which.cf,
                    control = control)

  x@na.action <- na.act  # Restore its original value
  invisible(x)
}




ylim.scale <- function(ylim, scale = 0) {
  if (length(ylim) != 2 ||
      ylim[2] < ylim[1])
    stop("error in 'ylim'")
  try <- ylim[2] - ylim[1]
  if (try > scale) ylim else
    c(ylim[1] + ylim[2] - scale,
      ylim[1] + ylim[2] + scale) / 2
}





getallresponses <- function(xij) {
  if (!is.list(xij))
    return("")

  allterms <- lapply(xij, terms)
  allres <- NULL
  for (ii in seq_along(xij))
    allres <- c(allres,
                as.character(attr(allterms[[ii]], "variables"))[2])
  allres
}



 headpreplotvgam <-
  function(object, newdata = NULL,
           terms = attr((object@terms)$terms, "term.labels"),
           raw = TRUE, deriv.arg = deriv.arg, se = FALSE,
           varxij = 1) {
  Terms <- terms(object)  # 20030811; object@terms$terms
  aa <- attributes(Terms)
  all.terms <- labels(Terms)
  xvars <- parse(text = all.terms)




  names(xvars) <- all.terms
  terms <- sapply(terms, match.arg, all.terms)

  Interactions <- aa$order > 1
  if (any(Interactions)) {
    stop("cannot handle interactions")
  }

  xvars <- xvars[terms]
  xnames <- as.list(terms)
  names(xnames) <- terms
  modes <- sapply(xvars, mode)
  for (term in terms[modes != "name"]) {
    evars <- all.names(xvars[term], functions = FALSE, unique = TRUE)
    if (!length(evars))
      next
    xnames[[term]] <- evars
    evars <- parse(text=evars)
    if (length(evars) == 1) {
      evars <- evars[[1]]
    } else
    if (length(evars) > 1 &&
        length(intersect(getallresponses(object@control$xij), names(xnames)))
       ) {




      evars <- evars[[varxij]]
    } else {
      evars <- c(as.name("list"), evars)
      mode(evars) <- "call"
    }
    xvars[[term]] <- evars
  }


  xvars <- c(as.name("list"), xvars)
  mode(xvars) <- "call"
  if (length(newdata)) {
    xvars <- eval(xvars, newdata)
  } else {
    Call <- object@call
    if (!is.null(Call$subset) | !is.null(Call$na.action) |
        !is.null(options("na.action")[[1]])) {
      Rownames <- names(fitted(object))
      if (!(Rl <- length(Rownames)))
        Rownames <- dimnames(fitted(object))[[1]]

      if (length(object@x) && !(Rl <- length(Rownames)))
        Rownames <- (dimnames(object@x))[[1]]
      if (length(object@y) && !(Rl <- length(Rownames)))
        Rownames <- (dimnames(object@y))[[1]]

      if (!(Rl <- length(Rownames)))
        stop("need to have names for fitted.values ",
             "when call has a 'subset' or 'na.action' argument")

      form <- paste("~", unlist(xnames), collapse = "+")
      Mcall <- c(as.name("model.frame"), list(formula =
                 terms(as.formula(form)),
                 subset = Rownames, na.action = function(x) x))
      mode(Mcall) <- "call"
      Mcall$data <- Call$data
      xvars <- eval(xvars, eval(Mcall))
    } else {
      ecall <- substitute(eval(expression(xvars)))
      ecall$local <- Call$data
      xvars <- eval(ecall)
    }
  }
  list(xnames = xnames, xvars = xvars)
}








preplotvgam <-
  function(object, newdata = NULL,
           terms = attr((object@terms)$terms, "term.labels"),
           raw = TRUE, deriv.arg = deriv.arg, se = FALSE,
           varxij = 1) {

  result1 <- headpreplotvgam(object, newdata = newdata, terms = terms,
                             raw = raw, deriv.arg = deriv.arg, se = se,
                             varxij = varxij)

  xvars  <- result1$xvars
  xnames <- result1$xnames



  if (FALSE && !is.null(object@control$jix)) {




    myxij <- object@control$xij
    if (length(myxij)) {
    }
  }

  pred <- if (length(newdata)) {
    predict(object, newdata, type = "terms",
            raw = raw, se.fit = se, deriv.arg = deriv.arg)
  } else {
    predict(object, type = "terms",
            raw = raw, se.fit = se, deriv.arg = deriv.arg)
  }

  fits <- if (is.atomic(pred)) NULL else pred$fit
  se.fit <- if (is.atomic(pred)) NULL else pred$se.fit

  if (is.null(fits))
    fits <- pred
  fred <- attr(fits, "vterm.assign")   # NULL for M==1
  Constant <- attr(fits, "constant")  # NULL if se = TRUE

  gamplot <- xnames

  loop.var <- names(fred)
  for (term in loop.var) {
    .VGAM.x <- xvars[[term]]

    myylab <- if (all(substring(term, 1:nchar(term),
                                      1:nchar(term)) != "("))
              paste("partial for", term) else term

    TT <- list(x = .VGAM.x,
               y = fits[, (if (is.null(fred)) term else fred[[term]])],
               se.y = if (is.null(se.fit)) NULL else
                     se.fit[, (if (is.null(fred)) term else fred[[term]])],
               xlab = xnames[[term]],
               ylab = myylab)
    class(TT) <- "preplotvgam"
    gamplot[[term]] <- TT
  }
  attr(gamplot, "Constant") <- Constant
  invisible(gamplot)
}






 plotpreplotvgam <-
  function(x, y = NULL, residuals = NULL,
           rugplot = TRUE, se = FALSE, scale = 0,
           offset.arg = 0, deriv.arg = 0, overlay = FALSE,
           which.term = NULL, which.cf = NULL,
           control = NULL) {
  listof <- inherits(x[[1]], "preplotvgam")
  if (listof) {
    TT <- names(x)
    if (is.null(which.term))
      which.term <- TT  # Plot them all
    plot.no <- 0
    for (ii in TT) {
      plot.no <- plot.no + 1
      if ((is.character(which.term) && any(which.term == ii)) ||
          (is.numeric(which.term) && any(which.term == plot.no)))
        plotpreplotvgam(x[[ii]], y = NULL,
                        residuals, rugplot = rugplot, se = se,
                        scale = scale,
                        offset.arg = offset.arg,
                        deriv.arg = deriv.arg, overlay = overlay,
                        which.cf = which.cf,
                        control = control)
    }
  } else {
    dummy <- function(residuals = NULL, rugplot = TRUE,
                      se = FALSE, scale = 0,
                      offset.arg = 0, deriv.arg = 0, overlay = FALSE,
                      which.cf = NULL, control = plotvgam.control())
     c(list(residuals = residuals, rugplot = rugplot,
            se = se, scale = scale,
            offset.arg = offset.arg, deriv.arg = deriv.arg,
            overlay = overlay, which.cf = which.cf), control)

    dd <- dummy(residuals = residuals, rugplot = rugplot,
                se = se, scale = scale,
                offset.arg = offset.arg, deriv.arg = deriv.arg,
                overlay = overlay,
                which.cf = which.cf,
                control = control)

    uniq.comps <- unique(c(names(x), names(dd)))
    Call <- c(as.name("vplot"), c(dd, x)[uniq.comps])
    mode(Call) <- "call"
    invisible(eval(Call))
  }
}


vplot.default <- function(x, y, se.y = NULL, xlab = "", ylab = "",
                          residuals = NULL, rugplot = FALSE,
                          scale = 0, se = FALSE,
                          offset.arg = 0, deriv.arg = 0, overlay = FALSE,
                          which.cf = NULL, ...) {
  switch(data.class(x)[1],
         logical = vplot.factor(factor(x), y, se.y, xlab, ylab, residuals,
                                rugplot, scale, se,
                                offset.arg = offset.arg,
                                overlay = overlay, ...),
         if (is.numeric(x)) {
           vplot.numeric(as.vector(x), y, se.y, xlab, ylab,
                         residuals, rugplot, scale, se,
                         offset.arg = offset.arg, overlay = overlay, ...)
         } else {
           warning("The 'x' component of '", ylab, "' has class '",
                   class(x), "'; no vplot() methods available")
         }
        )  # End of switch
}



vplot.list <-
  function(x, y, se.y = NULL, xlab, ylab,
           residuals = NULL, rugplot = FALSE, scale = 0, se = FALSE,
           offset.arg = 0, deriv.arg = 0, overlay = FALSE,
           which.cf = NULL, ...) {

  if (is.numeric(x[[1]])) {
    vplot.numeric(x[[1]], y, se.y, xlab, ylab,
                  residuals, rugplot, scale, se,
                  offset.arg = offset.arg, deriv.arg = deriv.arg,
                  overlay = overlay, ...)
  } else {
    stop("this function has not been written yet")
  }
}




 plotvgam.control <-
  function(which.cf = NULL,
           xlim = NULL, ylim = NULL,
           llty = par()$lty,
           slty = "dashed",
           pcex = par()$cex,
           pch = par()$pch,
           pcol = par()$col,
           lcol = par()$col,
           rcol = par()$col,
           scol = par()$col,
           llwd = par()$lwd,
           slwd = par()$lwd,
           add.arg = FALSE,
           one.at.a.time = FALSE,
           .include.dots = TRUE,
           noxmean = FALSE,
           shade = FALSE, shcol = "gray80",
           ...) {


  ans <-
  list(which.cf = which.cf,
       xlim = xlim, ylim = ylim,
       llty = llty, slty = slty,
       pcex = pcex, pch = pch,
       pcol = pcol, lcol = lcol, rcol = rcol, scol = scol,
       llwd = llwd, slwd = slwd,
       add.arg = add.arg,
       noxmean = noxmean,
       one.at.a.time = one.at.a.time,
       shade = shade, shcol = shcol)

  if (.include.dots) {
    c(list(...), ans)
  } else {
    default.vals <- plotvgam.control()
    return.list <- list()
    for (ii in names(default.vals)) {
      replace.val <- !((length(ans[[ii]]) == length(default.vals[[ii]])) &&
            (length(default.vals[[ii]]) > 0) &&
            identical(ans[[ii]], default.vals[[ii]]))

      if (replace.val)
        return.list[[ii]] <- ans[[ii]]
    }
    if (length(return.list)) {
      names(return.list) <- names(return.list)
      return.list
    } else {
      NULL
    }
  }
}






vplot.numeric <-
  function(x, y, se.y = NULL, xlab, ylab,
           residuals = NULL, rugplot = FALSE,
           se = FALSE, scale = 0,
           offset.arg = 0, deriv.arg = 0, overlay = FALSE,
           which.cf = NULL,
           xlim = NULL, ylim = NULL,
           llty = par()$lty,
           slty = "dashed",
           pcex = par()$cex,
           pch = par()$pch,
           pcol = par()$col,
           lcol = par()$col,
           rcol = par()$col,
           scol = par()$col,
           llwd = par()$lwd,
           slwd = par()$lwd,
           add.arg = FALSE,
           one.at.a.time = FALSE,
           noxmean = FALSE,
           separator = ":",
           shade = FALSE, shcol = "gray80",
           ...) {




    ylim0 <- ylim

    if (length(y)/length(x)  != round(length(y)/length(x)))
      stop("length of 'x' and 'y' do not seem to match")
    y <- as.matrix(y)
    if (!length(which.cf))
      which.cf <- 1:ncol(y)  # Added 20040807

    if (!is.null(se.y))
      se.y <- as.matrix(se.y)
    if (!is.null(se.y) && anyNA(se.y))
      se.y <- NULL

    if (!is.null(residuals))  {
      residuals <- as.matrix(residuals)
      if (ncol(residuals) != ncol(y)) {
        warning("ncol(residuals) != ncol(y) so residuals are not plotted")
        residuals <- NULL
      }
    }

    offset.arg <- matrix(offset.arg, nrow(y), ncol(y), byrow = TRUE)
    y <- y + offset.arg

    ylab <- add.hookey(ylab, deriv.arg)


    if (xmeanAdded <-
       (se && !is.null(se.y) && !noxmean &&
        all(substring(ylab, 1:nchar(ylab), 1:nchar(ylab)) != "("))) {
      x <- c(x, mean(x))
      y <- rbind(y, 0 * y[1, ])
      se.y <- rbind(se.y, 0 * se.y[1, ])
      if (!is.null(residuals))
        residuals <- rbind(residuals, NA*residuals[1, ])  # NAs not plotted
    }

    ux <- unique(sort(x))
    ooo <- match(ux, x)
    uy <- y[ooo, , drop = FALSE]


    xlim.orig <- xlim
    ylim.orig <- ylim
    xlim <- range(if (length(xlim)) NULL else ux, xlim, na.rm = TRUE)
    ylim <- range(if (length(ylim)) NULL else uy[, which.cf],
                  ylim, na.rm = TRUE)


    if (rugplot) {
      usex <- if (xmeanAdded) x[-length(x)] else x
      jx <- jitter(usex[!is.na(usex)])
      xlim <- range(if (length(xlim.orig)) NULL else jx,
                    xlim.orig, na.rm = TRUE)
    }

    if (se && !is.null(se.y)) {
      se.upper <- uy + 2 * se.y[ooo, , drop = FALSE]
      se.lower <- uy - 2 * se.y[ooo, , drop = FALSE]

      ylim <- if (length(ylim.orig)) range(ylim.orig) else
              range(c(ylim, se.upper[, which.cf], se.lower[, which.cf]))
    }

    if (!is.null(residuals)) {
      if (length(residuals) == length(y)) {
        residuals <- as.matrix(y + residuals)
        ylim <- if (length(ylim.orig)) range(ylim.orig) else
                range(c(ylim, residuals[, which.cf]), na.rm = TRUE)
      } else {
        residuals <- NULL
        warning("Residuals do not match 'x' in \"", ylab,
                "\" preplot object")
      }
    }


  all.missingy <- all(is.na(y))

  if (all.missingy)
    return()

  if (!length(ylim.orig))
    ylim <- ylim.scale(ylim, scale)

  if (overlay) {
    if (!length(which.cf))
      which.cf <- 1:ncol(uy)  # Added 20040807
    if (!add.arg) {
      matplot(ux, uy[, which.cf], type = "n",
              xlim = xlim, ylim = ylim,
              xlab = xlab, ylab = ylab, ...)
    }
    matlines(ux, uy[, which.cf],
             lwd = llwd, col = lcol, lty = llty)
    if (!is.null(residuals)) {
      if (ncol(y) == 1) {
        points(x, residuals, pch = pch, col = pcol, cex = pcex)
      } else {
        matpoints(x, residuals[, which.cf],
                  pch = pch, col = pcol, cex = pcex)  # add.arg = TRUE,
      }
    }
    if (rugplot)
      rug(jx, col = rcol)
    if (se && !is.null(se.y)) {
     matlines(ux, se.upper[, which.cf], lty =  slty, lwd = slwd, col = scol)
     matlines(ux, se.lower[, which.cf], lty =  slty, lwd = slwd, col = scol)
    }
  } else {
    YLAB <- ylab

    pcex <- rep_len(pcex,  ncol(uy))
    pch  <- rep_len(pch ,  ncol(uy))
    pcol <- rep_len(pcol,  ncol(uy))
    lcol <- rep_len(lcol,  ncol(uy))
    llty <- rep_len(llty,  ncol(uy))
    llwd <- rep_len(llwd,  ncol(uy))
    slty <- rep_len(slty,  ncol(uy))
    rcol <- rep_len(rcol,  ncol(uy))
    scol <- rep_len(scol,  ncol(uy))
    slwd <- rep_len(slwd,  ncol(uy))

    for (ii in 1:ncol(uy)) {
      if (!length(which.cf) ||
         ( length(which.cf) && any(which.cf == ii))) {

        if (is.Numeric(ylim0, length.arg = 2)) {
          ylim <- ylim0
        } else {
          ylim <- range(ylim0, uy[, ii], na.rm = TRUE)
          if (se && !is.null(se.y))
            ylim <- range(ylim0, se.lower[, ii], se.upper[, ii],
                          na.rm = TRUE)
          if (!is.null(residuals))
            ylim <- range(c(ylim, residuals[, ii]), na.rm = TRUE)
          ylim <- ylim.scale(ylim, scale)
        }
        if (ncol(uy) > 1 && length(separator))
          YLAB <- paste(ylab, separator, ii, sep = "")

        if (!add.arg) {
          if (one.at.a.time) {
            readline("Hit return for the next plot ")
          }
           plot(ux, uy[, ii], type = "n",
                xlim = xlim, ylim = ylim,
                xlab = xlab, ylab = YLAB, ...)
        }

        lines(ux, uy[, ii], lwd = llwd[ii], col = lcol[ii], lty = llty[ii])
        if (!is.null(residuals))
          points(x, residuals[, ii], pch = pch[ii],
                 col = pcol[ii], cex = pcex[ii])
        if (rugplot)
          rug(jx, col = rcol[ii])

        if (se && !is.null(se.y)) {
          if (shade) {
            polygon(c(ux, rev(ux), ux[1]),
                    c(se.upper[, ii], rev(se.lower[, ii]), se.upper[1, ii]),
                    col = shcol, border = NA)
            lines(ux, uy[, ii], lwd = llwd[ii], col = lcol[ii], lty = llty[ii])
          } else {
            lines(ux, se.upper[, ii], lty = slty[ii], lwd = slwd[ii],
                  col = scol[ii])
            lines(ux, se.lower[, ii], lty = slty[ii], lwd = slwd[ii],
                  col = scol[ii])
          }  # !shade
        }  # se && !is.null(se.y))
      }
    }  # for()
  }  # overlay
}  # vplot.numeric()



vplot.matrix <-
  function(x, y, se.y = NULL, xlab, ylab,
           residuals = NULL, rugplot = FALSE, scale = 0, se = FALSE,
           offset.arg = 0, deriv.arg = 0, overlay = FALSE,
           which.cf = NULL, ...) {
  stop("You shouldn't ever call this function!")
}


add.hookey <- function(ch, deriv.arg = 0) {

  if (!is.Numeric(deriv.arg, integer.valued = TRUE,
                  length.arg = 1) ||
      deriv.arg < 0)
      stop("bad input for the 'deriv' argument")

  if (deriv.arg == 0)
    return(ch)

  hookey <- switch(deriv.arg, "'", "''", "'''", "''''",
                              "'''''", stop("too high a derivative"))
  nc <- nchar(ch)
  sub <- substring(ch, 1:nc, 1:nc)
  if (nc >= 2 && sub[1] == "s" && sub[2] == "(") {
    paste("s", hookey, substring(ch, 2, nc), sep = "", coll = "")
  } else {
    paste(ch, hookey, sep = "", collapse = "")
  }
}



vplot.factor <-
  function(x, y, se.y = NULL, xlab, ylab,
           residuals = NULL, rugplot = FALSE, scale = 0,
           se = FALSE, xlim = NULL, ylim = NULL,
           offset.arg = 0, deriv.arg = 0, overlay = FALSE,
           which.cf = NULL, ...) {
  if (deriv.arg > 0)
    return(NULL)

  if (length(y)/length(x)  != round(length(y)/length(x)))
    stop("length of 'x' and 'y' do not seem to match")
  y <- as.matrix(y)

  if (!is.null(se.y))
    se.y <- as.matrix(se.y)
  if (!is.null(se.y) && anyNA(se.y))
    se.y <- NULL

  if (!is.null(residuals))  {
    residuals <- as.matrix(residuals)
    if (ncol(residuals) != ncol(y)) {
      warning("ncol(residuals) != ncol(y) so residuals are not plotted")
      residuals <- NULL
    }
  }

    if (overlay) {
      vvplot.factor(x, y,
                    se.y = if (is.null(se.y)) NULL else se.y,
                    xlab = xlab, ylab = ylab,
                    residuals = residuals,
                    rugplot = rugplot, scale = scale,
                    se = se, xlim = xlim, ylim = ylim, ...)
  } else {
    for (ii in 1:ncol(y)) {
      ylab <- rep_len(ylab, ncol(y))
      if (ncol(y) > 1)
        ylab <- dimnames(y)[[2]]
      vvplot.factor(x, y[, ii,drop = FALSE],
                    se.y = if (is.null(se.y)) NULL else
                           se.y[, ii,drop = FALSE],
                    xlab = xlab, ylab = ylab[ii],
                    residuals = if (is.null(residuals))
                        NULL else residuals[, ii,drop = FALSE],
                    rugplot = rugplot, scale = scale,
                    se = se, xlim = xlim, ylim = ylim, ...)

    }
  }
  invisible(NULL)
}





vvplot.factor <-
  function(x, y, se.y = NULL, xlab, ylab,
           residuals = NULL, rugplot = FALSE, scale = 0,
           se = FALSE, xlim = NULL, ylim = NULL,
           ...) {

  M <- ncol(y)
  nn <- as.numeric(table(x))
  codex <- as.numeric(x)
  ucodex <- seq(nn)[nn > 0]
  ooo <- match(ucodex, codex, 0)

  uy <- y[ooo, , drop = FALSE]
  ylim <- range(ylim, uy)
  xlim <- range(c(0, sum(nn), xlim))
  rightx <- cumsum(nn)
  leftx <- c(0, rightx[ -length(nn)])
  ux <- (leftx + rightx)/2
  delta <- (rightx - leftx)/8

  jx <- runif(length(codex), (ux - delta)[codex], (ux + delta)[codex])
  nnajx <- jx[!is.na(jx)]

  if (rugplot)
    xlim <- range(c(xlim, nnajx))
  if (se && !is.null(se.y)) {
    se.upper <- uy + 2 * se.y[ooo, , drop = FALSE]
    se.lower <- uy - 2 * se.y[ooo, , drop = FALSE]
    ylim <- range(c(ylim, se.upper, se.lower))
  }
  if (!is.null(residuals)) {
    if (length(residuals) == length(y)) {
      residuals <- y + residuals
      ylim <- range(c(ylim, residuals))
    } else {
      residuals <- NULL
      warning("Residuals do not match 'x' in \"", ylab,
              "\" preplot object")
    }
  }
  ylim <- ylim.scale(ylim, scale)
  Levels <- levels(x)
  if (!all(nn)) {
    keep <- nn > 0
    nn <- nn[keep]
    ux <- ux[keep]
    delta <- delta[keep]
    leftx <- leftx[keep]
    rightx <- rightx[keep]
    Levels <- Levels[keep]
  }


  about <- function(ux, M, Delta = 1 / M) {
    if (M == 1) return(cbind(ux))
    ans <- matrix(NA_real_, length(ux), M)
    grid <- seq(-Delta, Delta, len = M)
    for (ii in 1:M) {
      ans[, ii] <- ux + grid[ii]
    }
    ans
  }

  uxx <- about(ux, M, Delta = min(delta))
  xlim <- range(c(xlim, uxx))

  matplot(ux, uy, ylim = ylim, xlim = xlim, xlab = "", type = "n",
          ylab = ylab, axes = FALSE, frame.plot = TRUE, ...)
  mtext(xlab, 1, 2, adj = 0.5)
  axis(side = 2)
  lpos <- par("mar")[3]
  mtext(Levels, side = 3, line = lpos/2, at = ux, adj = 0.5, srt = 45)

  for (ii in 1:M)
    segments(uxx[, ii] - 1.0 * delta, uy[, ii],
             uxx[, ii] + 1.0 * delta, uy[, ii])
  if (!is.null(residuals)) {
    for (ii in 1:M) {
      jux <- uxx[, ii]
      jux <- jux[codex]
      jux <- jux + runif(length(jux), -0.7*min(delta), 0.7*min(delta))
      if (M == 1) points(jux, residuals[, ii]) else
                  points(jux, residuals[, ii], pch = as.character(ii))
    }
  }
  if (rugplot)
    rug(nnajx)
  if (se) {
    for (ii in 1:M) {
      segments(uxx[, ii] + 0.5*delta, se.upper[, ii],
               uxx[, ii] - 0.5*delta, se.upper[, ii])
      segments(uxx[, ii] + 0.5*delta, se.lower[, ii],
               uxx[, ii] - 0.5*delta, se.lower[, ii])
      segments(uxx[, ii], se.lower[, ii],
               uxx[, ii], se.upper[, ii], lty = 2)
    }
  }
  invisible(diff(ylim))
}


if (!isGeneric("vplot"))
  setGeneric("vplot", function(x, ...) standardGeneric("vplot"))


setMethod("vplot", "factor", function(x, ...)
           vplot.factor(x, ...))
setMethod("vplot", "list", function(x, ...)
           vplot.list(x, ...))
setMethod("vplot", "matrix", function(x, ...)
           vplot.matrix(x, ...))
setMethod("vplot", "numeric", function(x, ...)
           vplot.numeric(x, ...))




setMethod("plot", "vgam",
           function(x, y, ...) {
           if (!missing(y))
             stop("cannot process the 'y' argument")
           invisible(plot.vgam(x = x, y = y, ...))})






plotqrrvglm <- function(object,
               rtype = c("response", "pearson", "deviance", "working"),
               ask = FALSE,
               main = paste(Rtype, "residuals vs latent variable(s)"),
               xlab = "Latent Variable",
               I.tolerances = object@control$eq.tolerances,
               ...) {
  M <- object@misc$M
  n <- object@misc$n
  Rank <- object@control$Rank
  Coef.object <- Coef(object, I.tolerances = I.tolerances)
  rtype <- match.arg(rtype,
                     c("response", "pearson", "deviance", "working"))[1]
  res <- resid(object, type = rtype)

  my.ylab <- if (length(object@misc$ynames)) object@misc$ynames else
             rep_len(" ", M)
  Rtype <- switch(rtype, pearson = "Pearson", response = "Response",
                  deviance = "Deviance", working = "Working")

  done <- 0
  for (rr in 1:Rank)
    for (ii in 1:M) {
      plot(Coef.object@latvar[, rr],
           res[, ii],
           xlab = paste(xlab, if (Rank == 1) "" else rr, sep = ""),
           ylab = my.ylab[ii],
           main = main, ...)
      done <- done + 1
      if (done >= prod(par()$mfrow) && ask && done != Rank*M) {
          done <- 0
          readline("Hit return for the next plot: ")
      }
    }
  object
}


setMethod("plot", "qrrvglm", function(x, y, ...)
         invisible(plotqrrvglm(object = x, ...)))







put.caption <- function(text.arg = "(a)",
                        w.x = c(0.50, 0.50),
                        w.y = c(0.07, 0.93), ...) {
  text(text.arg,
       x = weighted.mean(par()$usr[1:2], w = w.x),
       y = weighted.mean(par()$usr[3:4], w = w.y), ...)
}










setMethod("plot", "pvgam",
           function(x, y, ...) {
           if (!missing(y))
             stop("cannot process the 'y' argument")
           invisible(plot.vgam(x = x, y = y, ...))})








# These functions are
# Copyright (C) 1998-2011 T.W. Yee, University of Auckland.
# All rights reserved.















 rcam <- function(y,
         family = poissonff,
         Rank = 0,
         Musual = NULL,
         Index.corner = if (!Rank) NULL else 1 + Musual * (1:Rank),
         rprefix = "Row.",
         cprefix = "Col.",
         szero = if (!Rank) NULL else
                           { if (Musual == 1) 1 else
                                  setdiff(1:(Musual*ncol(y)),
                                          c( # 1:Musual,
                                            1 + (1:ncol(y)) * Musual,
                                            Index.corner))},
                  summary.arg = FALSE, h.step = 0.0001,
                  rbaseline = 1, cbaseline = 1, ...) {
                           






  if (!is.character(rprefix))
    stop("argument 'rprefix' must be character")
  if (!is.character(cprefix))
    stop("argument 'cprefix' must be character")

  if (is.character(family))
      family <- get(family)
  if (is.function(family))
      family <- ((family)())
  if (!inherits(family, "vglmff")) {
      stop("'family = ", family, "' is not a VGAM family function")
  }
  efamily = family


  if (!is.Numeric(Musual)) {
    iefamily <- efamily@infos
    if (is.function(iefamily))
      Musual <- (iefamily())$Musual
  }
  if (!is.Numeric(Musual)) {
    warning("cannot determine the value of 'Musual'.",
            "Assuming the value one.")
    Musual <- 1
  }


  object.save <- y
  y <- if (is(y, "rrvglm")) {
    object.save@y
  } else {
    as(as.matrix(y), "matrix")
  }
  if (length(dim(y)) != 2 || nrow(y) < 3 || ncol(y) < 3)
    stop("argument 'y' must be a matrix with >= 3 rows & columns, or ",
         "a rrvglm() object")



  eifun <- function(i, n) diag(n)[, i, drop = FALSE]

  .rcam.df <- data.frame("Row.2" = eifun(2, nrow(y)))
  colnames( .rcam.df )<- paste(rprefix, "2", sep = "")  # Overwrite "Row.2"



  yn1 <- if (length(dimnames(y)[[1]])) dimnames(y)[[1]] else
            paste("X2.", 1:nrow(y), sep = "")
  warn.save = options()$warn
  options(warn = -3)    # Suppress the warnings (hopefully, temporarily)
  if (any(!is.na(as.numeric(substring(yn1, 1, 1)))))
      yn1 <- paste("X2.", 1:nrow(y), sep = "")
  options(warn = warn.save)


  nrprefix <- as.name(rprefix)
  ncprefix <- as.name(cprefix)


  assign(rprefix, factor(1:nrow(y)))
  modmat.row <- substitute(
           model.matrix( ~ .rprefix ), list( .rprefix = nrprefix ))
  assign(cprefix, factor(1:ncol(y)))
  modmat.col <- substitute(
           model.matrix( ~ .cprefix ), list( .cprefix = ncprefix ))
  modmat.row <- eval( modmat.row )
  modmat.col <- eval( modmat.col )








  Hlist <- list("(Intercept)" = matrix(1, ncol(y), 1))

  for(ii in 2:nrow(y)) {
    Hlist[[   paste(rprefix, ii, sep = "")]] <- matrix(1, ncol(y), 1)


    .rcam.df[[paste(rprefix, ii, sep = "")]] <- modmat.row[, ii]
  }


  for(ii in 2:ncol(y)) {


    Hlist[[   paste(cprefix, ii, sep = "")]] <- modmat.col[, ii, drop = FALSE]
    .rcam.df[[paste(cprefix, ii, sep = "")]] <- rep(1, nrow(y))
  }

  if (Rank > 0) {
    for(ii in 2:nrow(y)) {
      Hlist[[yn1[ii]]] <- diag(ncol(y))
      .rcam.df[[yn1[ii]]] <- eifun(ii, nrow(y))
    }
  }


  dimnames(.rcam.df) <- list(if (length(dimnames(y)[[1]]))
                           dimnames(y)[[1]] else
                           as.character(1:nrow(y)),
                           dimnames(.rcam.df)[[2]])

  str1 <- paste("~ ", rprefix, "2", sep = "")

  if (nrow(y) > 2) 
    for(ii in 3:nrow(y)) {
      str1 <- paste(str1, paste(rprefix, ii, sep = ""), sep = " + ")
    }



  for(ii in 2:ncol(y)) {
    str1 <- paste(str1, paste(cprefix, ii, sep = ""), sep = " + ")
  }


  str2 <- paste("y ", str1)
  if (Rank > 0) {
    for(ii in 2:nrow(y))
      str2 <- paste(str2, yn1[ii], sep = " + ")
  }


  controlfun <- if (Rank == 0) rrvglm.control else rrvglm.control
  controlfun <- if (Rank == 0)   vglm.control else rrvglm.control  # orig.


  mycontrol <- controlfun(Rank = Rank,
                          Index.corner = Index.corner,
                          szero = szero, ...)

  if (mycontrol$trace) {
  }



  if ((mindim <- min(nrow(y), ncol(y))) <= Rank) {
    stop("argument 'Rank' is too high. Must be a value from 0 ",
         "to ", mindim - 1, " inclusive")
  }



  if (Rank > 0)
    mycontrol$Norrr <- as.formula(str1)  # Overwrite this

  assign(".rcam.df", .rcam.df, envir = VGAM:::VGAMenv)

  warn.save <- options()$warn
  options(warn = -3)  # Suppress the warnings (hopefully, temporarily)

  if (mycontrol$trace) {
  }


  if (Musual > 1) {
    orig.Hlist <- Hlist
    for (ii in 1:length(Hlist))
      Hlist[[ii]] <- kronecker(Hlist[[ii]], rbind(1, 0))
    Hlist[["(Intercept)"]] <-
      cbind(Hlist[["(Intercept)"]],
            kronecker(matrix(1, nrow(orig.Hlist[[1]]), 1), rbind(0, 1)))



    if (mycontrol$trace) {
    }

  }





  answer <- if (Rank > 0) {
    if (is(object.save, "rrvglm")) object.save else 
      rrvglm(as.formula(str2),
             family = family,
             constraints = Hlist,
             control = mycontrol, data = .rcam.df)
  } else {
    if (is(object.save, "vglm")) object.save else 
        vglm(as.formula(str2),
             family = family,
             constraints = Hlist,
             control = mycontrol, data = .rcam.df)
  }

  options(warn = warn.save)


  answer <- if (summary.arg) {
    if (Rank > 0) {
      summary.rrvglm(as(answer, "rrvglm"), h.step = h.step)
    } else { 
      summary(answer)
    }
  } else { 
    as(answer, ifelse(Rank > 0, "rcam", "rcam0"))
  }


  answer@misc$rbaseline <- rbaseline
  answer@misc$cbaseline <- cbaseline

  answer
}








summaryrcam = function(object, ...) {
    rcam(object, summary.arg = TRUE, ...)
}









 setClass("rcam0", representation(not.needed = "numeric"),
          contains = "vglm")  # Added 20110506

 setClass("rcam", representation(not.needed = "numeric"),
          contains = "rrvglm")


setMethod("summary", "rcam0",
          function(object, ...)
          summaryrcam(object, ...))


setMethod("summary", "rcam",
          function(object, ...)
          summaryrcam(object, ...))










 Rcam <- function (mat, rbaseline = 1, cbaseline = 1) {

  mat <- as.matrix(mat)
  RRR <- dim(mat)[1]
  CCC <- dim(mat)[2]
    
  if (is.null(rownames(mat))) 
    rnames <- paste("X", 1:RRR, sep = "") else  
                 rnames  <- rownames(mat)

  if (is.null(colnames(mat))) 
    cnames <- paste("Y", 1:CCC, sep = "") else  
                 cnames  <- colnames(mat)

  r.index <- if (is.character(rbaseline))  
               which(rownames(mat) == rbaseline) else
                     if (is.numeric(rbaseline)) rbaseline else
                         stop("argement 'rbaseline' must be numeric", 
                               "or character of the level of row")
 
  c.index <- if (is.character(cbaseline))  
               which(colnames(mat) == cbaseline) else
                     if (is.numeric(cbaseline)) cbaseline else
                         stop("argement 'cbaseline' must be numeric",
                               "or character of the level of row")

  if (length(r.index) != 1)
    stop("Could not match with argument 'rbaseline'")

  if (length(c.index) != 1)
    stop("Could not match with argument 'cbaseline'")


  yswap <- rbind(mat[r.index:RRR, ],
                 if (r.index > 1) mat[1:(r.index - 1),] else NULL)
  yswap <- cbind(yswap[, c.index:CCC],
                 if (c.index > 1) yswap[, 1:(c.index - 1)] else  NULL)

  new.rnames <- rnames[c(r.index:RRR,
                         if (r.index > 1) 1:(r.index - 1) else NULL)]
  new.cnames <- cnames[c(c.index:CCC, 
                         if (c.index > 1) 1:(c.index - 1) else NULL)]
  colnames(yswap) <- new.cnames
  rownames(yswap) <- new.rnames
  
  yswap
}













 plotrcam0  <- function (object,
     centered = TRUE, whichplots = c(1, 2),
     hline0 = TRUE, hlty = "dashed", hcol = par()$col, hlwd = par()$lwd,
                         rfirst = 1, cfirst = 1,
                         rtype = "h", ctype = "h",
                         rcex.lab = 1, rcex.axis = 1, # rlabels = FALSE,
                         rtick = FALSE,
                         ccex.lab = 1, ccex.axis = 1, # clabels = FALSE,
                         ctick = FALSE,
                         rmain = "Row effects", rsub = "",
                         rxlab = "", rylab = "Row effects",
                         cmain = "Column effects", csub = "",
                         cxlab = "", cylab = "Column effects",
                         rcol = par()$col, ccol = par()$col,
                         ...) {

 
  nparff <- if (is.numeric(object@family@infos()$Musual)) {
    object@family@infos()$Musual
  } else {
    1
  }



  if (is.numeric(object@control$Rank) && object@control$Rank != 0)
    warning("argument 'object' is not Rank-0")


  n_lm  = nrow(object@y)

  cobj <- coefficients(object)

  upperbound = if (!is.numeric(object@control$Rank) ||
                   object@control$Rank == 0) length(cobj) else
               length(object@control$colx1.index)

  orig.roweff <- c("Row.1" = 0, cobj[(nparff + 1) : (nparff + n_lm - 1)])
  orig.coleff <- c("Col.1" = 0, cobj[(nparff + n_lm) : upperbound])
  last.r <- length(orig.roweff)
  last.c <- length(orig.coleff)


  orig.raxisl  <- rownames(object@y)
  orig.caxisl  <- colnames(object@y) 
    
  roweff.orig <- 
  roweff <- orig.roweff[c(rfirst:last.r,
                          if (rfirst > 1) 1:(rfirst-1) else NULL)]
  coleff.orig <- 
  coleff <- orig.coleff[c(cfirst:last.c,
                          if (cfirst > 1) 1:(cfirst-1) else NULL)]

  if (centered) {
    roweff = scale(roweff, scale = FALSE)  # Center it only
    coleff = scale(coleff, scale = FALSE)  # Center it only
  }

  raxisl <- orig.raxisl[c(rfirst:last.r,
                          if (rfirst > 1) 1:(rfirst-1) else NULL)]

  caxisl <- orig.caxisl[c(cfirst:last.c, 
                          if (cfirst > 1) 1:(cfirst-1) else NULL)]


  if (any(whichplots == 1, na.rm = TRUE)) {
    plot(roweff, type = rtype, 
         axes = FALSE, col = rcol, main = rmain,
         sub  = rsub, xlab = rxlab, ylab = rylab, ...)

    axis(1, at = 1:length(raxisl),
         cex.lab = rcex.lab,  
         cex.axis = rcex.axis,
         label = raxisl)
    axis(2, cex.lab = rcex.lab, ...)  # las = rlas)

    if (hline0)
      abline(h = 0, lty = hlty, col = hcol, lwd = hlwd)
  }


  if (any(whichplots == 2, na.rm = TRUE)) {
    plot(coleff, type = ctype, 
         axes = FALSE, col = ccol, main = cmain, # lwd = 2, xpd = FALSE,
         sub  = csub, xlab = cxlab, ylab = cylab, ...)

    axis(1, at = 1:length(caxisl),
         cex.lab = ccex.lab,
         cex.axis = ccex.axis,
         label = caxisl)
    axis(2, cex.lab = ccex.lab, ...)  # las = clas)
    
    if (hline0)
      abline(h = 0, lty = hlty, col = hcol, lwd = hlwd)
  }




  object@post$row.effects = roweff
  object@post$col.effects = coleff
  object@post$raw.row.effects = roweff.orig
  object@post$raw.col.effects = coleff.orig

  invisible(object)
}





setMethod("plot", "rcam0",
          function(x, y, ...)
          plotrcam0(object = x, ...))


setMethod("plot", "rcam",
          function(x, y, ...)
          plotrcam0(object = x, ...))















moffset <- function (mat, roffset = 0, coffset = 0, postfix = "") {





  if ((is.numeric(roffset) && (roffset == 0)) &&
      (is.numeric(coffset) && (coffset == 0)))
    return(mat)


  vecmat = c(unlist(mat))
  ind1 <- if (is.character(roffset))
                 which(rownames(mat) == roffset) else
                       if (is.numeric(roffset)) roffset + 1 else
                           stop("argument 'roffset' not matched (character). ",
                                 "It must be numeric, ",
                                 "else character and match the ",
                                 "row names of the response")
  ind2 <- if (is.character(coffset))
                 which(colnames(mat) == coffset) else
                       if (is.numeric(coffset)) coffset + 1 else
                           stop("argument 'coffset' not matched (character). ",
                                 "It must be numeric, ",
                                 "else character and match the ",
                                 "column names of the response")

  if (!is.Numeric(ind1, positive = TRUE, integ = TRUE, allow = 1) ||
      !is.Numeric(ind2, positive = TRUE, integ = TRUE, allow = 1))
    stop("bad input for arguments 'roffset' and/or 'coffset'")
  if (ind1 > nrow(mat))
    stop("too large a value for argument 'roffset'")
  if (ind2 > ncol(mat))
    stop("too large a value for argument 'coffset'")


  start.ind = (ind2 - 1)* nrow(mat) + ind1


  svecmat = vecmat[c(start.ind:(nrow(mat) * ncol(mat)),
                     0:(start.ind - 1))]

  rownames.mat = rownames(mat)
  if (length(rownames.mat) != nrow(mat))
    rownames.mat = paste("Row.", 1:nrow(mat), sep = "")

  colnames.mat = colnames(mat)
  if (length(colnames.mat) != ncol(mat))
    colnames.mat = paste("Col.", 1:ncol(mat), sep = "")


  newrn = if (roffset > 0)
            c(rownames.mat[c(ind1:nrow(mat))],
              paste(rownames.mat[0:(ind1-1)], postfix, sep = "")) else
           rownames.mat

  newcn = c(colnames.mat[c(ind2:ncol(mat), 0:(ind2 - 1))])
  if (roffset > 0)
    newcn = paste(newcn, postfix, sep = "")

  newmat = matrix(svecmat, nrow(mat), ncol(mat),
                  dimnames = list(newrn, newcn))
  newmat
}



















confint_rrnb <- function(rrnb2) {

  if (class(rrnb2) != "rrvglm")
    stop("argument 'rrnb2' does not appear to be a rrvglm() object")

  if (!any(rrnb2@family@vfamily == "negbinomial"))
    stop("argument 'rrnb2' does not appear to be a negbinomial() fit")

  if (rrnb2@control$Rank != 1)
    stop("argument 'rrnb2' is not Rank-1")

  if (rrnb2@misc$M != 2)
    stop("argument 'rrnb2' does not have M = 2")

  if (!all(rrnb2@misc$link == "loge"))
    stop("argument 'rrnb2' does not have log links for both parameters")

  a21.hat <- (Coef(rrnb2)@A)["log(size)", 1]
  beta11.hat <- Coef(rrnb2)@B1["(Intercept)", "log(mu)"]
  beta21.hat <- Coef(rrnb2)@B1["(Intercept)", "log(size)"]
  delta1.hat <- exp(a21.hat * beta11.hat - beta21.hat)
  delta2.hat <- 2 - a21.hat

  se.a21.hat <- sqrt(vcovrrvglm(rrnb2)["I(lv.mat)", "I(lv.mat)"])


  ci.a21 <- a21.hat +  c(-1, 1) * 1.96 * se.a21.hat
  (ci.delta2 <- 2 - rev(ci.a21))  # The 95 percent confidence interval

  list(a21.hat    = a21.hat,
       beta11.hat = beta11.hat,
       beta21.hat = beta21.hat,
       ci.delta2  = ci.delta2,
       delta1     = delta1.hat,
       delta2     = delta2.hat,
       se.a21.hat = se.a21.hat)
}



confint_nb1 <- function(nb1) {

  if (class(nb1) != "vglm")
    stop("argument 'nb1' does not appear to be a vglm() object")

  if (!any(nb1@family@vfamily == "negbinomial"))
    stop("argument 'nb1' does not appear to be a negbinomial() fit")

  if (!all(unlist(constraints(nb1)[-1]) == 1))
    stop("argument 'nb1' does not appear to have parallel = TRUE")

  if (!all(unlist(constraints(nb1)[1]) == c(diag(nb1@misc$M))))
    stop("argument 'nb1' does not have parallel = FALSE for the intercept")

  if (nb1@misc$M != 2)
    stop("argument 'nb1' does not have M = 2")

  if (!all(nb1@misc$link == "loge"))
    stop("argument 'nb1' does not have log links for both parameters")

  cnb1 <- coefficients(as(nb1, "vglm"), matrix = TRUE)
  mydiff <- (cnb1["(Intercept)", "log(size)"] - cnb1["(Intercept)", "log(mu)"])
  delta0.hat <- exp(mydiff)
  (phi0.hat <- 1 + 1 / delta0.hat)  # MLE of phi0

  myvcov <- vcovvlm(as(nb1, "vglm"))  # Not great; improve this!
  myvec <- cbind(c(-1, 1, rep(0, len = nrow(myvcov) - 2)))
  (se.mydiff <- sqrt(t(myvec) %*%  myvcov %*%  myvec))
  ci.mydiff <- mydiff + c(-1.96, 1.96) * se.mydiff
  ci.delta0 <- ci.exp.mydiff <- exp(ci.mydiff)
  (ci.phi0 <- 1 + 1 / rev(ci.delta0)) # The 95 percent conf int. for phi0

  list(ci.phi0    = ci.phi0,
       delta0     = delta0.hat,
       phi0       = phi0.hat)
}












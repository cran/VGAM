# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.















 rcam <- function(y,
         family = poissonff,
         Rank = 0,
         Musual = NULL,
         weights = NULL,
         which.lp = 1,
         Index.corner = if (!Rank) NULL else 1 + Musual * (1:Rank),
         rprefix = "Row.",
         cprefix = "Col.",
         offset = 0,
         szero = if (!Rank) NULL else
                           { if (Musual == 1) 1 else
                                  setdiff(1:(Musual*ncol(y)),
                                          c( # 1:Musual,
                                            1 + (1:ncol(y)) * Musual,
                                            Index.corner))},
         summary.arg = FALSE, h.step = 0.0001,
         rbaseline = 1, cbaseline = 1,
         ...) {
                           




  noroweffects = FALSE
  nocoleffects = FALSE

  if (!is.Numeric(which.lp, allowable.length = 1, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'which.lp'")

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

  .rcam.df <-
    if (!noroweffects) data.frame("Row.2" = eifun(2, nrow(y))) else
    if (!nocoleffects) data.frame("Col.2" = eifun(2, nrow(y))) else
    stop("at least one of 'noroweffects' and 'nocoleffects' must be FALSE")
  colnames( .rcam.df ) <- paste(rprefix, "2", sep = "") # Overwrite "Row.2"



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

  if (!noroweffects)
  for(ii in 2:nrow(y)) {
    Hlist[[   paste(rprefix, ii, sep = "")]] <- matrix(1, ncol(y), 1)


    .rcam.df[[paste(rprefix, ii, sep = "")]] <- modmat.row[, ii]
  }


  if (!nocoleffects)
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
  controlfun <- if (Rank == 0)   vglm.control else rrvglm.control # orig.


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

  assign(".rcam.df", .rcam.df, envir = VGAM::VGAMenv)

  warn.save <- options()$warn
  options(warn = -3)  # Suppress the warnings (hopefully, temporarily)

  if (mycontrol$trace) {
  }


  if (Musual > 1) {
    orig.Hlist <- Hlist
    kmat1 = rbind(1, 0)
    kmat0 = rbind(0, 1)

    kmat1 = matrix(0, nrow = Musual, ncol = 1)
    kmat1[which.lp, 1] = 1
    kmat0 = matrix(1, nrow = Musual, ncol = 1)
    kmat0[which.lp, 1] = 0

    for (ii in 1:length(Hlist)) {
      Hlist[[ii]] <- kronecker(Hlist[[ii]],
                               kmat1)
    }
    Hlist[["(Intercept)"]] <-
      cbind(Hlist[["(Intercept)"]],
            kronecker(matrix(1, nrow(orig.Hlist[[1]]), 1),
                      kmat0))



    if (mycontrol$trace) {
    }

  }



  offset.matrix = matrix(offset, nrow = nrow(y),
                                 ncol = ncol(y) * Musual) # byrow = TRUE

  answer <- if (Rank > 0) {
    if (is(object.save, "rrvglm")) object.save else
      rrvglm(as.formula(str2),
             family = family,
             constraints = Hlist,
             offset = offset.matrix,
             weights = if (length(weights)) weights else rep(1, length = nrow(y)),
             ...,
             control = mycontrol, data = .rcam.df)
  } else {
    if (is(object.save, "vglm")) object.save else
        vglm(as.formula(str2),
             family = family,
             constraints = Hlist,
             offset = offset.matrix,
             weights = if (length(weights)) weights else rep(1, length = nrow(y)),
             ...,
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
  answer@misc$which.lp  <- which.lp
  answer@misc$offset    <- offset.matrix

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
  if (is.null(orig.raxisl))
    orig.raxisl = as.character(1:nrow(object@y))
  if (is.null(orig.caxisl))
    orig.caxisl = as.character(1:ncol(object@y))
    
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
         labels = raxisl)
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
         labels = caxisl)
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

  if (!is.Numeric(ind1, positive = TRUE, integer.valued = TRUE, allowable.length = 1) ||
      !is.Numeric(ind2, positive = TRUE, integer.valued = TRUE, allowable.length = 1))
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



















confint_rrnb <- function(rrnb2, level = 0.95) {

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


  ci.a21 <- a21.hat +  c(-1, 1) * qnorm(1 - (1 - level)/2) * se.a21.hat
  (ci.delta2 <- 2 - rev(ci.a21)) # e.g., the 95 percent CI

  list(a21.hat    = a21.hat,
       beta11.hat = beta11.hat,
       beta21.hat = beta21.hat,
       CI.a21     = ci.a21,
       ci.delta2  = ci.delta2,
       delta1     = delta1.hat,
       delta2     = delta2.hat,
       se.a21.hat = se.a21.hat)
}





confint_nb1 <- function(nb1, level = 0.95) {




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

  ci.mydiff <- mydiff + c(-1, 1) * qnorm(1 - (1 - level)/2) * se.mydiff

  ci.delta0 <- ci.exp.mydiff <- exp(ci.mydiff)
  (ci.phi0 <- 1 + 1 / rev(ci.delta0)) # e.g., the 95 percent CI for phi0

  list(ci.phi0    = ci.phi0,
       CI.delta0  = ci.delta0,
       delta0     = delta0.hat,
       phi0       = phi0.hat)
}




plota21 <- function(rrvglm2, plot.it = TRUE, nseq.a21 = 31,
                    se.eachway = c(5, 5), # == c(LHS, RHS),
                    trace.arg = TRUE, ...) {





  if (class(rrvglm2) != "rrvglm")
    stop("argument 'rrvglm2' does not appear to be a rrvglm() object")

  if (rrvglm2@control$Rank != 1)
    stop("argument 'rrvglm2' is not Rank-1")

  if (rrvglm2@misc$M != 2)
    stop("argument 'rrvglm2' does not have M = 2")


  loglik.orig <- logLik(rrvglm2)
  temp1 <- confint_rrnb(rrvglm2) # zz

  a21.hat <- (Coef(rrvglm2)@A)[2, 1]
  se.a21.hat <- temp1$se.a21.hat


  se.a21.hat <- sqrt(vcov(rrvglm2)["I(lv.mat)", "I(lv.mat)"])


  big.ci.a21 <- a21.hat +  c(-1, 1) * se.eachway * se.a21.hat
  seq.a21 <- seq(big.ci.a21[1], big.ci.a21[2], length = nseq.a21)
  Hlist.orig <- constraints.vlm(rrvglm2, type = "lm")


  alreadyComputed <- !is.null(rrvglm2@post$a21.matrix)


  a21.matrix <- if (alreadyComputed) rrvglm2@post$a21.matrix else
                cbind(a21 = seq.a21, loglikelihood = 0)
  prev.etastart <- predict(rrvglm2) # Halves the computing time
  funname <- "vglm"
  listcall <- as.list(rrvglm2@call)


  if (!alreadyComputed)
  for (ii in 1:nseq.a21) {
    if (trace.arg)
      print(ii)

       argslist <- vector("list", length(listcall) - 1)
       for (kay in 2:(length(listcall)))
         argslist[[kay - 1]] <- listcall[[kay]]

       names(argslist) <- c(names(listcall)[-1])

       argslist$trace       = trace.arg
       argslist$etastart    = prev.etastart
       argslist$constraints = Hlist.orig


       for (kay in 2:length(argslist[["constraints"]])) {
         argslist[["constraints"]][[kay]] <- rbind(1, a21.matrix[ii, 1])
       }


       fitnew = do.call(what = funname, args = argslist)

       a21.matrix[ii, 2] <- logLik(fitnew)

       prev.etastart <- predict(fitnew)
  }



  if (plot.it) {
    plot(a21.matrix[ ,1], a21.matrix[ ,2], type = "l",
            col = "blue",
            xlab = expression(a[21]), ylab = "Log-likelihood") # ...

    abline(v = (Hlist.orig[[length(Hlist.orig)]])[2, 1],
           col = "darkorange", lty = "dashed")

    abline(h = loglik.orig,
           col = "darkorange", lty = "dashed")

    abline(h = loglik.orig - qchisq(0.95, df = 1),
           col = "darkorange", lty = "dashed")

    abline(v = a21.hat +  c(-1, 1) * 1.96 * se.a21.hat,
           col = "gray50", lty = "dashed", lwd = 2.0)

  } # End of (plot.it)

  rrvglm2@post <- list(a21.matrix = a21.matrix)
  invisible(rrvglm2)
}











if (FALSE)
Qvar <- function(object, factor.name = NULL,
                 level1.name = "level1",
                 ...) {




  object.xlevels = if (is(object, "vglm")) {
    object@xlevels
  } else {
    factor(rownames(object))
  }

  myvcov = if (is(object, "vglm")) {
    if (length(object.xlevels) == 0)
      stop("no factors amongst the model.matrix of 'object'")

    if (is.null(factor.name)) {
      if (length(object.xlevels) > 1)
        stop("more than one factor in the model.matrix of 'object'")

      factor.name = names(object.xlevels)
      object.xlevels = object@xlevels[[1]]
    } else {
      object.xlevels = object.xlevels[[factor.name]]
    }





    colptr = attr(model.matrix(object), "vassign")
    colptr = colptr[[factor.name]]
    vcov(object)[colptr, colptr, drop = FALSE]
  } else if (is.matrix(object)) {
    object
  } else {
    stop("argument 'object' is not a vglm() object or a matrix")
  }




  myvcov = rbind(0, cbind(0, myvcov))

  LL = nrow(myvcov)
  if (LL <= 3)
    stop("the factor must have at least three levels")


  vcov0 = myvcov
  for (ilocal in 1:LL)
    for (jlocal in ilocal:LL)
      myvcov[ilocal, jlocal] =
      myvcov[jlocal, ilocal] = vcov0[ilocal, ilocal] + vcov0[jlocal, jlocal] -
                               2 * vcov0[ilocal, jlocal]

  allvcov = myvcov
  rownames(allvcov) =
    c(paste(if (is.matrix(object)) level1.name else factor.name,
            if (is.matrix(object)) NULL else object.xlevels[1], sep = ""),
            rownames(vcov0)[-1])
  colnames(allvcov) = rownames(allvcov)









  diag(allvcov) = rep(1,             len = LL) # Any positive value should do





  Allvcov = allvcov




  wmat   = matrix(1, LL, LL)
  diag(wmat) = sqrt( .Machine$double.eps )


  logAllvcov = log(Allvcov)
  attr(logAllvcov, "Prior.Weights") = wmat
  logAllvcov
}












Qvar <- function(object, factorname = NULL, coef.indices = NULL,
                 labels = NULL, dispersion = NULL,
                 reference.name = "(reference)",
                 estimates = NULL
                ) {




  coef.indices.saved <- coef.indices
  if (!is.matrix(object)) {
      model <- object
      if (is.null(factorname) && is.null(coef.indices)) {
        stop("arguments \"factorname\" and \"coef.indices\" are both NULL")
      }

      if (is.null(coef.indices)) {
          tmodel <- terms(model)
          modelmat <- if (is.matrix(model@x)) model@x else
                         model.matrix(tmodel,
                                      data = model@model)




          colptr = attr(model.matrix(object), "vassign")
          colptr = colptr[[factorname]]
          coef.indices <- colptr


          contmat <- if (length(model@xlevels[[factorname]]) ==
                         length(coef.indices)) {
              diag(length(coef.indices))
          } else {
              eval(call(model@contrasts[[factorname]],
                        model@xlevels[[factorname]]))
          }
          rownames(contmat) <- model@xlevels[[factorname]]

          if (is.null(estimates))
            estimates <- contmat %*% coefvlm(model)[coef.indices]

          covmat <- vcovvlm(model, dispersion = dispersion)
          covmat <- covmat[coef.indices, coef.indices, drop = FALSE]
          covmat <- contmat %*% covmat %*% t(contmat)
      } else {
          k <- length(coef.indices)
          refPos <- numeric(0)
          if (0 %in% coef.indices) {
              refPos <- which(coef.indices == 0)
              coef.indices <- coef.indices[-refPos]
          }
          covmat <- vcovvlm(model, dispersion = dispersion)
          covmat <- covmat[coef.indices, coef.indices, drop = FALSE]

          if (is.null(estimates))
            estimates <- coefvlm(model)[coef.indices]

          if (length(refPos) == 1) {
              if (length(estimates) != k)
                estimates <- c(0, estimates)
              covmat <- rbind(0, cbind(0, covmat))
              names(estimates)[1] <-
              rownames(covmat)[1] <-
              colnames(covmat)[1] <- reference.name
              if (refPos != 1) {
                perm <- if (refPos == k) c(2:k, 1) else
                        c(2:refPos, 1, (refPos + 1):k)
                  estimates <- estimates[perm]
                  covmat <- covmat[perm, perm, drop = FALSE]
              }
          }
      }


      return(Recall(covmat,
                    factorname = factorname,
                    coef.indices = coef.indices.saved,
                    labels = labels,
                    dispersion = dispersion,
                    estimates = estimates
                    )
             )
  } else {

      covmat <- object
      if (length(labels))
        rownames(covmat) <- colnames(covmat) <- labels
      if ((LL <- dim(covmat)[1]) <= 2)
        stop("This function works only for factors with 3 or more levels")
  }




  allvcov = covmat
  for (ilocal in 1:LL)
    for (jlocal in ilocal:LL)
      allvcov[ilocal, jlocal] =
      allvcov[jlocal, ilocal] = covmat[ilocal, ilocal] +
                                covmat[jlocal, jlocal] -
                                covmat[ilocal, jlocal] * 2

  diag(allvcov) = rep(1.0, len = LL) # Any positive value should do


  wmat   = matrix(1.0, LL, LL)
  diag(wmat) = sqrt( .Machine$double.eps )

  logAllvcov = log(allvcov)
  attr(logAllvcov, "Prior.Weights") = wmat
  attr(logAllvcov, "estimates") = estimates
  attr(logAllvcov, "coef.indices") = coef.indices
  attr(logAllvcov, "factorname") = factorname
  attr(logAllvcov, "regularVar") = diag(covmat)

  logAllvcov
}






WorstErrors <- function(qv.object) {
  stop("20110729; does not work")

    reducedForm <- function(covmat, qvmat){
        nlevels <- dim(covmat)[1]
     firstRow <- covmat[1, ]
     ones <- rep(1, nlevels)
     J <- outer(ones, ones)
     notzero <- 2:nlevels
     r.covmat <- covmat + (firstRow[1]*J) -
                      outer(firstRow, ones) -
                     outer(ones, firstRow)
     r.covmat <- r.covmat[notzero, notzero]
     qv1 <- qvmat[1, 1]
     r.qvmat <- (qvmat + qv1*J)[notzero, notzero]
     list(r.covmat, r.qvmat)}
    covmat <- qv.object$covmat
    qvmat <- diag(qv.object$qvframe$quasiVar)
    r.form <- reducedForm(covmat, qvmat)
    r.covmat <- r.form[[1]]
    r.qvmat <- r.form[[2]]
    inverse.sqrt <- solve(chol(r.covmat))
    evalues <- eigen(t(inverse.sqrt) %*% r.qvmat %*% inverse.sqrt,
                    symmetric = TRUE)$values
    sqrt(c(min(evalues), max(evalues))) - 1
}




IndentPrint <- function(object, indent = 4, ...){
  stop("20110729; does not work")

    zz <- ""
    tc <- textConnection("zz", "w", local = TRUE)
    sink(tc)
    try(print(object, ...))
    sink()
    close(tc)
    indent <- paste(rep(" ", indent), sep = "", collapse = "")
    cat(paste(indent, zz, sep = ""), sep = "\n")}



Print.qv <- function(x, ...){
  stop("20110729; does not work")

}



summary.qvar <- function(object, ...) {


  relerrs = 1 - sqrt(exp(residuals(object, type = "response")))
  diag(relerrs) = NA

    minErrSimple <- round(100 * min(relerrs, na.rm = TRUE), 1)
    maxErrSimple <- round(100 * max(relerrs, na.rm = TRUE), 1)



  estimates = c(object@extra$attributes.y$estimates)
  if (!length(names(estimates)) &&
      is.matrix(object@extra$attributes.y$estimates))
    names( estimates) = rownames(object@extra$attributes.y$estimates)
  if (!length(names(estimates)))
    names( estimates) = paste("Level", 1:length(estimates), sep = "")


  regularVar = c(object@extra$attributes.y$regularVar)
  QuasiVar <- exp(diag(fitted(object))) / 2
  QuasiSE  <- sqrt(QuasiVar)


  structure(list(estimate = estimates,
                 SE            = sqrt(regularVar),  # zz dispersion parameter??
                 minErrSimple  = minErrSimple,
                 maxErrSimple  = maxErrSimple,
                 quasiSE  = QuasiSE,
                 object   = object,
                 quasiVar = QuasiVar),
            class = "summary.qvar")
}



print.summary.qvar <- function(x, ...) {

  object = x$object
  minErrSimple  = x$minErrSimple
  maxErrSimple  = x$maxErrSimple

  x$minErrSimple = NULL
  x$maxErrSimple = NULL
  x$object = NULL


    if (length(cl <- object@call)) {
        cat("Call:\n")
        dput(cl)
    }


    facname = c(object@extra$attributes.y$factorname)
    if (length(facname))
      cat("Factor name: ", facname, "\n")


    if (length(object@dispersion))
        cat("\nDispersion: ", object@dispersion, "\n\n")

  x = as.data.frame(c(x))
  print.data.frame(x)


        cat("\nWorst relative errors in SEs of simple contrasts (%): ",
            minErrSimple, ", ", maxErrSimple, "\n")

  invisible(x)
}



plotqvar <- function(object,
                     intervalWidth = 2,
                     ylab = "Estimate",
                     xlab = NULL, # x$factorname,
                     ylim = NULL,
                     main = "",
                     levelNames = NULL,
                     conf.level = 0.95,
                     warn.ratio = 10,
                     border = "transparent",  # None
                     points.arg = TRUE,
                     length.arrows = 0.25, angle = 30,
                     scol = par()$col,
                     slwd = par()$lwd,
                     slty = par()$lty,
                     ...) {




    if (!is.numeric(intervalWidth) &&
        !is.numeric(conf.level))
      stop("at least one of arguments 'intervalWidth' and 'conf.level' ",
            "should be numeric")





  if (!any("normal1" %in% object@family@vfamily))
    stop("argument 'object' dos not appear to be a ",
         "rcam(, normal1) object")

  estimates = c(object@extra$attributes.y$estimates)
  if (!length(names(estimates)) &&
      is.matrix(object@extra$attributes.y$estimates))
    names( estimates) = rownames(object@extra$attributes.y$estimates)
  if (!length(names(estimates)))
    names( estimates) = paste("Level", 1:length(estimates), sep = "")




  QuasiVar <- exp(diag(fitted(object))) / 2
  QuasiSE  <- sqrt(QuasiVar)

    if (!is.numeric(estimates))
      stop("Cannot plot, because there are no 'proper' parameter estimates")
    if (!is.numeric(QuasiSE))
      stop("Cannot plot, because there are no quasi standard errors")



    faclevels <- factor(names(estimates), levels = names(estimates))


    xvalues <- seq(along = faclevels)
    tops  <- estimates + intervalWidth * QuasiSE
    tails <- estimates - intervalWidth * QuasiSE




    if (is.numeric(conf.level)) {
      zedd = abs(qnorm((1 - conf.level) / 2))
      lsd.tops  <- estimates + zedd * QuasiSE / sqrt(2)
      lsd.tails <- estimates - zedd * QuasiSE / sqrt(2)
      if (max(QuasiSE) / min(QuasiSE) > warn.ratio)
        warning("Quasi SEs appear to be quite different... the ",
                "LSD intervals may not be very accurate")
    } else {
      lsd.tops  <- NULL
      lsd.tails <- NULL
    }




    if (is.null(ylim))
      ylim <- range(c(tails, tops, lsd.tails, lsd.tops), na.rm = TRUE)

    if (is.null(xlab))
      xlab <- "Factor level"

    plot(faclevels, estimates,
         border = border,
         ylim = ylim, xlab = xlab, ylab = ylab,
         main = main, ...)


    if (points.arg)
      points(estimates, ...)


    if (is.numeric(intervalWidth)) {
      segments(xvalues, tails, xvalues, tops,
               col = scol, lty = slty, lwd = slwd)
    }


    if (is.numeric(conf.level)) {
      arrows(xvalues, lsd.tails, xvalues, lsd.tops,
             col = scol, lty = slty, lwd = slwd, code = 3,
             length = length.arrows, angle = angle)

    }




  if (any(slotNames(object) == "post")) {
    object@post$estimates  = estimates 
    object@post$xvalues    = xvalues  
    if (is.numeric(intervalWidth)) {
      object@post$tails = tails
      object@post$tops  = tops
    }
    if (is.numeric(conf.level)) {
      object@post$lsd.tails = lsd.tails
      object@post$lsd.tops  = lsd.tops
    }
  }


    invisible(object)
}











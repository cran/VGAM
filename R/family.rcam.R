# These functions are
# Copyright (C) 1998-2011 T.W. Yee, University of Auckland.
# All rights reserved.

















 rcam <- function(y, Rank = 0,
         family = poissonff,
         Musual = NULL,
         Index.corner = if (!Rank) NULL else 1 + Musual * (1:Rank),
         rprefix = "Row.",
         cprefix = "Col.",
         szero = if (!Rank) NULL else
                           {if (Musual == 1) 1 else
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
  controlfun <- if (Rank == 0) vglm.control else rrvglm.control  # orig.


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
    as(answer, ifelse(Rank > 0, "rrvglm", "vglm"))
  }


  answer@misc$rbaseline <- rbaseline
  answer@misc$cbaseline <- cbaseline

  answer
}








summaryrcam = function(object, ...) {
    rcam(object, summary.arg = TRUE, ...)
}









 setClass("rcam", representation(not.needed = "numeric"),
          contains = "rrvglm")


setMethod("summary", "rcam",
          function(object, ...)
          summaryrcam(object, ...))














  Rcam <- function (mat, rbaseline = 1, cbaseline = 1) {


  mat <- as.matrix(mat)
  RRR <- dim(mat)[1]
  CCC <- dim(mat)[2]
    

  if (is.null(rownames(mat))) 
    rnames <- paste("X", 1:RRR, sep="") else  
                 rnames  <- rownames(mat)

  if (is.null(colnames(mat))) 
    cnames <- paste("Y", 1:CCC, sep="") else  
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
    

  yswap <- rbind(mat[r.index:RRR, ], 
                  if (r.index > 1) mat[1:(r.index - 1),] else NULL)


  if (length(r.index) != 1) 
    stop("Could not match with argument 'rbaseline'")

  if (length(c.index) != 1) 
    stop("Could not match with argument 'cbaseline'")

  yswap <- cbind(yswap[, c.index:CCC], 
               if (c.index > 1) yswap[, 1:(c.index - 1)] else  NULL)

  new.rnames <- rnames[ c(r.index:RRR, 
                if (r.index > 1) 1:(r.index - 1) else NULL)]
  
  new.cnames <- cnames[ c(c.index:CCC, 
                if (c.index > 1) 1:(c.index - 1) else NULL)]

  colnames(yswap) <- new.cnames
  rownames(yswap) <- new.rnames  
  
  yswap
}





 plotrcam0  <- function (object, rfirst = 1, cfirst = 1,
                         rtype = "h", ctype = "h",
                         rlas = 1, rcex.lab = 1, 
                         rcex.axis = 1, rlabels = FALSE,
                         rtick = FALSE, clas = 1, ccex.lab = 1,
                         ccex.axis = 1, clabels = FALSE, ctick = FALSE,
                         rmain = "Row effects", rsub = "",
                         rxlabel = "", rylabel = "Row effects",
                         cmain = "Column effects", csub = "", cxlabel= "",
                         cylabel = "Column effects",
                         rcol = par()$col, ccol = par()$col,
                         ...) {

 
  if (object@family@infos()$Musual == 1) nparff <- 1 else nparff <- 2


  orig.roweff <- c(0, coefficients(object)[(nparff+1): (nparff+nrow(object@y)-1)])
 
  orig.coleff <- c(0, coefficients(object)[(nparff+nrow(object@y)):
                                    (length(coefficients(object)))])
  rlast <- length(orig.roweff)
  clast <- length(orig.coleff)

  orig.raxisl  <- rownames(object@y)
  orig.caxisl  <- colnames(object@y) 
    
  roweff <- orig.roweff[c(rfirst:rlast, 
             if (rfirst > 1) 1:(rfirst-1) else NULL)]

  coleff <- orig.coleff[c(cfirst:clast, 
             if (cfirst > 1) 1:(cfirst-1) else NULL)]

  raxisl <- orig.raxisl[ c(rfirst:rlast, 
             if (rfirst > 1) 1:(rfirst-1) else NULL)]

  caxisl <- orig.caxisl[ c(cfirst:clast, 
             if (cfirst > 1) 1:(cfirst-1) else NULL)]


   plot(roweff, type = rtype, 
        axes = FALSE, col = rcol,
        main = rmain,
        sub  = rsub,
        xlab = rxlabel, ylab = rylabel, ...)

    
    axis(1, at = 1:length(raxisl), cex.lab = rcex.lab,  
            cex.axis = rcex.axis, label = raxisl)
    axis(2, cex.lab = rcex.lab, las = rlas)
    axis(3:4, labels = rlabels, tick = rtick)



    plot(coleff, type = ctype, col = ccol, # lwd=2, xpd=F,
         axes = FALSE, main = cmain,
         sub  = csub, xlab = cxlabel, ylab = cylabel, ...)

    axis(1, at = 1:length(caxisl), cex.lab = ccex.lab,
            cex.axis = ccex.axis, label = caxisl)
    axis(2, cex.lab= ccex.lab, las = clas)
    
    
  invisible(object)
}











moffset <- function (mat, roffset=1, coffset=1){

   y <- mat

   
    rowoffset <- function(y, roffset=1) {
      y <- as.matrix(y)  

      roffset <- if (is.character(roffset))  
                 which(rownames(y) == roffset) else
                       if (is.numeric(roffset)) roffset else
                           stop("argument rstart/cstart must be numeric ",
                                 "or character of the level of row/column")

      if (roffset == 1) ye <- y else {   
          ye <- y[1:roffset-1,,drop = FALSE]
          ye <- rbind(y[(roffset):nrow(y),,drop = FALSE],
                   cbind(ye[,2:ncol(y),drop = FALSE], ye[,1,drop = FALSE]))
          ye
      }
    }        
    
    if (((coffset >= 1) && (coffset <= ncol(y))) ||
        ((roffset >=1) && (roffset <= nrow(y)))) {
         y <- rowoffset(y, roffset)
         y <- t(rowoffset(t(y), coffset))
         y
         } else
           stop ("Error argument in 'rstart' or 'cstart'.",
                 "It must be numeric or chacarter argument of row or column of",
                 "'mat' matrix input")
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

  a21.hat <- (Coef(rrnb2)@A)["log(k)", 1]
  beta11.hat <- Coef(rrnb2)@B1["(Intercept)", "log(mu)"]
  beta21.hat <- Coef(rrnb2)@B1["(Intercept)", "log(k)"]
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
  mydiff <- (cnb1["(Intercept)", "log(k)"] - cnb1["(Intercept)", "log(mu)"])
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











  # ref: Kus, section 4.1, pg 4500
  # updated on 22/15/2010
dexppois <- function(x, lambda, betave = 1, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  N <- max(length(x), length(lambda), length(betave))
  x <- rep(x, len = N); lambda = rep(lambda, len = N);
  betave <- rep(betave, len = N)

  logdensity <- rep(log(0), len = N)
  xok <- (0 < x)
 # logdensity[xok] <- log(lambda[xok]) + log(betave[xok]) - lambda[xok] -
  #                   betave[xok] * x[xok] + lambda[xok] * exp(-betave[xok] *
   #                  x[xok]) - log(expm1(lambda[xok])) + 1
  
 logdensity[xok] <- log(lambda[xok]) + log(betave[xok]) -
                    log1p(-exp(-lambda[xok])) - lambda[xok] - betave[xok] *
                    x[xok] + lambda[xok] * exp(-betave[xok] * x[xok])
   
  logdensity[lambda <= 0] <- NaN
  logdensity[betave <= 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}

  # ref: calculated from F(x) from Kus, pg 4499
  # Not working 13/12/10
  # updated and working on 22/15/2010
rexppois <- function(n, lambda, betave = 1) {
  # ans <- log(lambda/(log((exp(lambda) + 1) * runif(n))))/betave
  ans <- -log(log(runif(n) * -(expm1(lambda)) +
         exp(lambda)) / lambda) / betave
  ans[(lambda <= 0) | (betave <= 0)] <- NaN
  ans
}


  # ref: calculated from F(x) from Kus, pg 4499
  # Not working 13/12/10
  # updated and working on 22/15/2010
qexppois<- function(p, lambda, betave = 1) {
  # ans <- log(lambda/(log((exp(lambda) + 1) * p)))/betave
  ans <- -log(log(p * -(expm1(lambda)) +
         exp(lambda)) / lambda) / betave
  ans[(lambda <= 0) | (betave <= 0)] = NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans
}



  # ref: Kus, eqn 2, pg 4499
  # Updated on 22/12/2010
pexppois<- function(q, lambda, betave = 1) {
  #ans <- -(exp(lambda * exp(-betave * q)) - exp(lambda))/expm1(lambda)
  ans <-(exp(lambda * exp(-betave * q)) - exp(lambda)) / -expm1(lambda)  
  ans[(lambda <= 0) | (betave <= 0)] <- NaN
  ans
}




 exppoisson = function (llambda = "loge", lbetave = "loge",
                        elambda = list(), ebetave = list(),
                        ilambda = 1.1,    ibetave = 1.5, 
                        zero = NULL) {

  if (mode(llambda) != "character" && mode(llambda) != "name") 
    llambda = as.character(substitute(llambda))
  if (mode(lbetave) != "character" && mode(lbetave) != "name") 
    lbetave = as.character(substitute(lbetave))

  if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE)) 
    stop("bad input for argument 'zero'")

  if (!is.Numeric(ilambda, posit = TRUE)) 
    stop("bad input for argument 'ilambda'")
  if (length(ibetave) && !is.Numeric(ibetave, posit = TRUE)) 
    stop("bad input for argument 'ibetave'")

  ilambda[ilambda == 1] = 1.1

  if (!is.list(ebetave)) 
    ebetave = list()
  if (!is.list(lambda)) 
    elambda = list()

  new("vglmff",
  blurb = c("Exponential Poisson Distribution \n \n", 
            "Links:    ",
             namesof("lambda", llambda, earg = elambda), ", ",
             namesof("betave", lbetave, earg = ebetave), "\n",
            "Mean:     ",
            "(lambda/(expm1(lambda) * betave)) *",
            "genhypergeo(c(1,1),c(2,2),lambda)"),

  # genhypergeo() from package: hypergeo
  # ref = mean from Kus pg 4499

  constraints = eval(substitute(expression({
    constraints = cm.zero.vgam(constraints, x, .zero , M)
    }), list( .zero = zero))),

  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")

    predictors.names = c(
      namesof("lambda", .llambda, earg = .elambda, short = TRUE),
      namesof("betave", .lbetave, earg = .ebetave, short = TRUE))
    if (!length(etastart)) {

      lambda.init = if (!is.Numeric( .ilambda , posit = TRUE))
                      stop("argument 'ilambda' must be positive") else
                     rep( .ilambda , len = n)
      betave.init = if (length( .ibetave )) 
                      rep( .ibetave , len = n) else
                      stop("zz need to fix this code")
                   ## (lambda.init/(expm1(lambda.init) * (y + 1/8))) *
                   ##  genhypergeo(c(1,1),c(2,2),lambda.init)


      betave.init = rep(weighted.mean(betave.init, w = w), len = n)
      etastart = cbind(theta2eta(lambda.init, .llambda , earg = .elambda ), 
                       theta2eta(betave.init, .lbetave , earg = .ebetave ))
    }
   }), list( .llambda = llambda, .lbetave = lbetave, 
             .ilambda = ilambda, .ibetave = ibetave, 
             .elambda = elambda, .ebetave = ebetave))), 

  inverse = eval(substitute(function(eta, extra = NULL) {
    lambda = eta2theta(eta[, 1], .llambda , earg = .elambda )
    betave = eta2theta(eta[, 2], .lbetave , earg = .ebetave )
  warning("returning dud means")
    runif(nrow(eta))
  }, list( .llambda = llambda, .lbetave = lbetave, 
           .elambda = elambda, .ebetave = ebetave))), 

  last = eval(substitute(expression({
    misc$link =    c(lambda = .llambda , betave = .lbetave )
    misc$earg = list(lambda = .elambda , betave = .ebetave )
    misc$expected = TRUE
    
  }), list( .llambda = llambda, .lbetave = lbetave,
            .elambda = elambda, .ebetave = ebetave))), 

  loglikelihood = eval(substitute(function(mu, y, w, 
                  residuals = FALSE, eta, extra = NULL) {
    lambda = eta2theta(eta[, 1], .llambda , earg = .elambda )
    betave = eta2theta(eta[, 2], .lbetave , earg = .ebetave )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
    sum(w * dexppois(x = y, lambda = lambda, betave = betave, log = TRUE))
    }
  }, list( .lbetave = lbetave , .llambda = llambda , 
           .elambda = elambda , .ebetave = ebetave ))), 

  vfamily = c("exppoisson"),

  deriv = eval(substitute(expression({
    lambda = eta2theta(eta[, 1], .llambda , earg = .elambda )
    betave = eta2theta(eta[, 2], .lbetave , earg = .ebetave )

    dl.dbetave = 1/betave - y - y * lambda * exp(-betave * y)
    dl.dlambda = 1/lambda - 1/expm1(lambda) - 1 + exp(-betave * y)

    dbetave.deta = dtheta.deta(betave, .lbetave , earg = .ebetave )
    dlambda.deta = dtheta.deta(lambda, .llambda , earg = .elambda )

    w * cbind(dl.dlambda * dlambda.deta,
              dl.dbetave * dbetave.deta)
  }), list( .llambda = llambda , .lbetave = lbetave,
            .elambda = elambda, .ebetave = ebetave ))), 

  weight = eval(substitute(expression({
    
    temp1 = -expm1(-lambda)
    
    ed2l.dlambda2 = (1 + exp(2 * lambda) - lambda^2 * exp(lambda) - 2 * 
                    exp(lambda)) / (lambda * temp1)^2


    ed2l.dbetave2 = 1 / betave^2 - (lambda^2 * exp(-lambda) / (4 * betave^2 *
                    temp1)) * genhypergeo(c(2,2,2),c(3,3,3),lambda) 

    ed2l.dbetavelambda = (lambda * exp(-lambda) / (4 * betave * temp1)) *
                         genhypergeo(c(2,2),c(3,3),lambda)   

    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M)] = ed2l.dlambda2 * dlambda.deta^2
    wz[, iam(2, 2, M)] = ed2l.dbetave2 * dbetave.deta^2
    wz[, iam(1, 2, M)] = dbetave.deta * dlambda.deta * ed2l.dbetavelambda
    w * wz
  }), list( .zero = zero ))))
}





dgenray <- function(x, shape, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  N <- max(length(x), length(shape), length(scale))
  x <- rep(x, len = N)
  shape <- rep(shape, len = N)
  scale <- rep(scale, len = N)

  logdensity <- rep(log(0), len = N)
  if (any(xok <- (x > 0))) {
    temp1 <- x[xok] / scale[xok]
    logdensity[xok] <- log(2) + log(shape[xok]) + log(x[xok]) -
                       2 * log(scale[xok]) - temp1^2  +
                       (shape[xok] - 1) * log1p(-exp(-temp1^2))
  }
  logdensity[(shape <= 0) | (scale <= 0)] <- NaN
  if (log.arg) {
    logdensity
  } else {
     exp(logdensity)
  }
}


pgenray <- function(q, shape, scale = 1) {
  ans <- (-expm1(-(q/scale)^2))^shape
  ans[q <= 0] <- 0
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}



rgenray <- function(n, shape, scale = 1) {
  ans <- qgenray(runif(n), shape = shape, scale = scale)
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}


qgenray <- function(p, shape, scale = 1) {
  ans <- scale * sqrt(-log1p(-(p^(1/shape))))
  ans[(shape <= 0) | (scale <= 0)] = NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans
}



      




genrayleigh.control <- function(save.weight = TRUE, ...)
{
    list(save.weight = save.weight)
}





 genrayleigh = function (lshape = "loge", lscale = "loge",
                         eshape = list(), escale = list(),
                         ishape = NULL,   iscale = NULL,
                         tol12 = 1.0e-05, 
                         nsimEIM = 300, zero = 1) {

  if (mode(lshape) != "character" && mode(lshape) != "name")
    lshape = as.character(substitute(lshape))
  if (mode(lscale) != "character" && mode(lscale) != "name")
    lscale = as.character(substitute(lscale))

  if (length(ishape) && !is.Numeric(ishape, posit = TRUE))
    stop("bad input for argument 'ishape'")
  if (length(iscale) && !is.Numeric(iscale, posit = TRUE)) 
    stop("bad input for argument 'iscale'")

  if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
    stop("bad input for argument 'zero'")
    if (!is.Numeric(nsimEIM, allow = 1, integ = TRUE) || nsimEIM <= 50)
        stop("'nsimEIM' should be an integer greater than 50")

  if (!is.list(escale))
    escale = list()
  if (!is.list(eshape))
    eshape = list()


  new("vglmff",
  blurb = c("Generalized Rayleigh distribution\n",
            "Links:    ",
            namesof("shape", lshape, earg = eshape), ", ",
            namesof("scale", lscale, earg = escale), "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1) 
      stop("response must be a vector or a one-column matrix")
       
    predictors.names = c(
      namesof("shape", .lshape , earg = .eshape , short = TRUE),
      namesof("scale", .lscale , earg = .escale , short = TRUE))

    if (!length(etastart)) {
      genrayleigh.Loglikfun = function(scale, y, x, w, extraargs) {
        temp1 <- y / scale
        shape = -1 / weighted.mean(log1p(-exp(-temp1^2)), w = w)

        ans <-
        sum(w * (log(2) + log(shape) + log(y) - 2 * log(scale) -
                 temp1^2  + (shape - 1) * log1p(-exp(-temp1^2))))
        ans
      }
      scale.grid = seq(0.2 * sd(y), 5 * sd(y), len = 29)
      scale.init = if (length( .iscale )) .iscale else
                   getMaxMin(scale.grid, objfun = genrayleigh.Loglikfun,
                             y = y, x = x, w = w)
      scale.init = rep(scale.init, length = length(y))
 
      shape.init = if (length( .ishape )) .ishape else
                   -1 / weighted.mean(log1p(-exp(-(y/scale.init)^2)), w = w)
      shape.init = rep(shape.init, length = length(y))

      etastart = cbind(theta2eta(shape.init, .lshape, earg = .eshape),
                       theta2eta(scale.init, .lscale, earg = .escale))
        }
    }), list( .lscale = lscale, .lshape = lshape,
              .iscale = iscale, .ishape = ishape,
              .escale = escale, .eshape = eshape))), 

  inverse = eval(substitute(function(eta, extra = NULL) {
    shape = eta2theta(eta[, 1], .lshape , earg = .eshape )
    Scale = eta2theta(eta[, 2], .lscale , earg = .escale )
    qgenray(p = 0.5, shape = shape, scale = Scale)
  }, list( .lshape = lshape, .lscale = lscale, 
           .eshape = eshape, .escale = escale ))),

  last = eval(substitute(expression({
    misc$link =    c(shape = .lshape , scale = .lscale )
    misc$earg = list(shape = .eshape , scale = .escale )
    misc$expected = TRUE
    misc$nsimEIM = .nsimEIM
  }), list( .lshape = lshape, .lscale = lscale,
            .eshape = eshape, .escale = escale,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(function(mu, y, w, 
                  residuals = FALSE, eta, extra = NULL) {

    shape = eta2theta(eta[, 1], .lshape , earg = .eshape )
    Scale = eta2theta(eta[, 2], .lscale , earg = .escale )

    if (residuals) stop("loglikelihood residuals",
                        "not implemented yet") else {
      sum(w * dgenray(x = y, shape = shape, scale = Scale, log = TRUE))
    }
  }, list( .lshape = lshape , .lscale = lscale , 
           .eshape = eshape , .escale = escale ))), 
      
  vfamily = c("genrayleigh"), 

  deriv = eval(substitute(expression({
    shape = eta2theta(eta[, 1], .lshape , earg = .eshape )
    Scale = eta2theta(eta[, 2], .lscale , earg = .escale )
    dshape.deta = dtheta.deta(shape, .lshape , earg = .eshape )
    dscale.deta = dtheta.deta(Scale, .lscale , earg = .escale )
    dthetas.detas = cbind(dshape.deta, dscale.deta)

    temp1 <- y / Scale
    temp2 <- exp(-temp1^2)
    temp3 <- temp1^2 / Scale
    AAA   <- 2 * temp1^2 / Scale  # 2 * y^2 / Scale^3
    BBB   <- -expm1(-temp1^2)     # denominator
    dl.dshape = 1/shape + log1p(-temp2)
    dl.dscale = -2 / Scale + AAA * (1 - (shape - 1) * temp2 / BBB)

    dl.dshape[!is.finite(dl.dshape)] = max(dl.dshape[is.finite(dl.dshape)])

    answer <- w * cbind(dl.dshape, dl.dscale) * dthetas.detas
    answer
  }), list( .lshape = lshape , .lscale = lscale,
            .eshape = eshape,  .escale = escale ))),

  weight = eval(substitute(expression({


    run.varcov = 0
    ind1 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    for(ii in 1:( .nsimEIM )) {
        ysim = rgenray(n = n, shape = shape, scale = Scale)

        temp1 <- ysim / Scale
        temp2 <- exp(-temp1^2)  # May be 1 if ysim is very close to 0.
        temp3 <- temp1^2 / Scale
        AAA   <- 2 * temp1^2 / Scale  # 2 * y^2 / Scale^3
        BBB   <- -expm1(-temp1^2)     # denominator
        dl.dshape = 1/shape + log1p(-temp2)
        dl.dscale = -2 / Scale + AAA * (1 - (shape - 1) * temp2 / BBB)

        dl.dshape[!is.finite(dl.dshape)] = max(dl.dshape[is.finite(dl.dshape)])

        temp3 = cbind(dl.dshape, dl.dscale)
        run.varcov = run.varcov + temp3[, ind1$row.index] *
                                  temp3[, ind1$col.index]
    }
    run.varcov = run.varcov / .nsimEIM

    wz = if (intercept.only)
        matrix(colMeans(run.varcov, na.rm = FALSE),
               n, ncol(run.varcov), byrow = TRUE) else run.varcov
    wz = wz * dthetas.detas[, ind1$row] * dthetas.detas[, ind1$col]
    w * wz
  }), list( .lshape = lshape , .lscale = lscale,
            .eshape = eshape,  .escale = escale,
            .tol12 = tol12, .nsimEIM = nsimEIM ))))
}




      




  # Ref: Mudholker pg 293
  # Updated and working: 06/01/11
desnorm <- function(x, location = 0, Scale = 1, epsilon = 0, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)
  
  N <- max(length(x), length(location), length(Scale), length(epsilon))
  x <- rep(x, len = N)
  location <- rep(location, len = N)
  Scale <- rep(Scale, len = N)
  epsilon <- rep(epsilon, len = N)
  zedd <- (x - location)/Scale
  
  logdensity <- rep(log(0), len = N)
  xneg <- (zedd < 0)
  xpos <- (zedd >= 0)
  logdensity[xneg] <- 1/2 * log(2 * pi) - zedd[xneg]^2/(2 * (1 + epsilon[xneg])^2)
  logdensity[xpos] <- 1/2 * log(2 * pi) - zedd[xpos]^2/(2 * (1 - epsilon[xpos])^2)
  logdensity[(epsilon < -1) | (epsilon > 1)] <- NaN  
  
  if (log.arg)
    logdensity
  else exp(logdensity)
  
}
      

  # Ref: Mudholker pg 293
  # Updated and working: 06/01/11
pesnorm <- function(q, location = 0, Scale = 1, epsilon = 0) {

  N <- max(length(q), length(location), length(Scale), length(epsilon))
  q <- rep(q, len = N)
  location <- rep(location, len = N)
  Scale <- rep(Scale, len = N)
  epsilon <- rep(epsilon, len = N)
  zedd <- (q - location)/Scale
  
  qneg <- (zedd < 0)
  qpos <- (zedd >= 0)
  ans <- rep(0, len = length(q))
  
  ans[qneg] <- (1 + epsilon[qneg]) * pnorm(q = zedd[qneg]/(1 + epsilon[qneg]),
                mean = 0, sd = 1)
  ans[qpos] <- epsilon[qpos] + (1 - epsilon[qpos]) * pnorm(q = zedd[qpos]/(1 - 
               epsilon[qpos]),mean = 0, sd = 1)
  
  ans
}











dexpgeom <- function(x, scale = 1, shape, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  N <- max(length(x), length(scale), length(shape))
  x <- rep(x, len = N)
  scale <- rep(scale, len = N)
  shape <- rep(shape, len = N)

  logdensity <- rep(log(0), len = N)
  if (any(xok <- (x > 0))) {
    temp1 <- (-x[xok]) * scale[xok]
    logdensity[xok] <- log(scale[xok]) + log1p(-shape[xok]) + 
                       temp1 - 2 * log1p(-shape[xok] * exp(temp1))
  }

  logdensity[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  if (log.arg) {
    logdensity
  } else {
     exp(logdensity)
  }
}



pexpgeom <- function(q, scale = 1, shape) {
  temp1 <- (-q) * scale
  ans <- -expm1(temp1) / (1 - shape * exp(temp1))
  ans[q <= 0] <- 0
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}


rexpgeom <- function(n, scale = 1, shape) {
  ans <- qexpgeom(runif(n), shape = shape, scale = scale)
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}


 
qexpgeom <- function(p, scale = 1, shape) {
  ans <- (-1/scale) * log((p - 1) / (p * shape - 1))
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans
}





expgeometric.control <- function(save.weight = TRUE, ...)
{
    list(save.weight = save.weight)
}



 expgeometric = function (lscale = "loge", lshape = "logit",
                          escale = list(), eshape = list(),
                          iscale = NULL,   ishape = NULL, 
                          zero = 1, nsimEIM = 400) {


  if (mode(lshape) != "character" && mode(lshape) != "name")
    lshape = as.character(substitute(lshape))
  if (mode(lscale) != "character" && mode(lscale) != "name")
    lscale = as.character(substitute(lscale))

  if (length(ishape))
    if (!is.Numeric(ishape, posit = TRUE) || any(ishape >= 1))
      stop("bad input for argument 'ishape'")

  if (length(iscale))
    if (!is.Numeric(iscale, posit = TRUE))
    stop("bad input for argument 'iscale'")

  if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
    stop("bad input for argument 'zero'")

  if (!is.list(escale))
    escale = list()
  if (!is.list(eshape))
    eshape = list()

  if (!is.Numeric(nsimEIM, allow = 1, integ = TRUE))
      stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
      stop("'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Exponential geometric distribution\n\n",
            "Links:    ",
            namesof("Scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     ", "(shape - 1) * log(1 - ",
            "shape) / (Scale * shape)"), 
                           
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
 


  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")

    predictors.names = c(
      namesof("Scale", .lscale , earg = .escale , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {

      scale.init = if (is.Numeric( .iscale , posit = TRUE)) {
                     rep( .iscale , len = n)
                   } else {
                      1 / sd(y)  # The papers scale parameter beta
                   }

      shape.init = if (is.Numeric( .ishape , posit = TRUE)) {
                     rep( .ishape , len = n)
                   } else {
                      rep(2 - exp(scale.init * median(y)), len = n)
                   }
      shape.init[shape.init >= 0.95] = 0.95
      shape.init[shape.init <= 0.05] = 0.05

      etastart = cbind(theta2eta(scale.init, .lscale , earg = .escale ),
                       theta2eta(shape.init, .lshape , earg = .eshape ))

    }
   }), list( .lscale = lscale, .lshape = lshape, 
             .iscale = iscale, .ishape = ishape, 
             .escale = escale, .eshape = eshape))), 

  inverse = eval(substitute(function(eta, extra = NULL) {
    Scale = eta2theta(eta[, 1], .lscale , earg = .escale )
    shape = eta2theta(eta[, 2], .lshape , earg = .eshape )
    
    (shape - 1) * log1p(-shape) / (Scale * shape)

  }, list( .lscale = lscale, .lshape = lshape, 
           .escale = escale, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link =    c(Scale = .lscale , shape = .lshape )
    misc$earg = list(Scale = .escale , shape = .eshape )
    misc$expected = TRUE
    misc$nsimEIM = .nsimEIM
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(function(mu, y, w, 
                  residuals = FALSE, eta, extra = NULL) {

    Scale = eta2theta(eta[, 1], .lscale , earg = .escale )
    shape = eta2theta(eta[, 2], .lshape , earg = .eshape )

    if (residuals) stop("loglikelihood residuals",
                        "not implemented yet") else {
      sum(w * dexpgeom(x = y, shape = shape, scale = Scale, log = TRUE))
    }
  }, list( .lscale = lscale , .lshape = lshape , 
           .escale = escale , .eshape = eshape ))), 
      
  vfamily = c("expgeometric"), 

  deriv = eval(substitute(expression({
    Scale = eta2theta(eta[, 1], .lscale , earg = .escale )
    shape = eta2theta(eta[, 2], .lshape , earg = .eshape )

     temp2 <- exp(-Scale * y)
     temp3 <- shape * temp2
     dl.dscale =  1 / Scale - y   - 2 * y * temp3 / (1 - temp3)
     dl.dshape = -1 / (1 - shape) + 2 *     temp2 / (1 - temp3)

    dscale.deta = dtheta.deta(Scale, .lscale , earg = .escale )            
    dshape.deta = dtheta.deta(shape, .lshape , earg = .eshape )
    dthetas.detas = cbind(dscale.deta, dshape.deta)

    answer <- w * cbind(dl.dscale, dl.dshape) * dthetas.detas
    answer
  }), list( .lscale = lscale , .lshape = lshape,
            .escale = escale,  .eshape = eshape ))),

  weight = eval(substitute(expression({
  








        run.varcov = 0
        ind1 = iam(NA, NA, M=M, both = TRUE, diag = TRUE)

        if (length( .nsimEIM )) {
            for(ii in 1:( .nsimEIM )) {
                ysim = rexpgeom(n, scale=Scale, shape=shape)

                temp2 <- exp(-Scale * ysim)
                temp3 <- shape * temp2
                dl.dscale = 1 / Scale - ysim - 2 * ysim *
                            temp3 / (1 - temp3)
                dl.dshape = -1 / (1 - shape) + 2 * temp2 / (1 - temp3)

                temp6 = cbind(dl.dscale, dl.dshape)
                run.varcov = run.varcov +
                           temp6[,ind1$row.index] * temp6[,ind1$col.index]
            }

            run.varcov = run.varcov / .nsimEIM

            wz = if (intercept.only)
                matrix(colMeans(run.varcov),
                       n, ncol(run.varcov), byrow = TRUE) else run.varcov

            wz = wz * dthetas.detas[, ind1$row] *
                      dthetas.detas[, ind1$col]
        }

    w * wz
  }), list( .nsimEIM = nsimEIM ))))
}










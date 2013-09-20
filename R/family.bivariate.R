# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.












dbiclaytoncop <- function(x1, x2, alpha = 0, log = FALSE){
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  A <- x1^(-alpha) + x2^(-alpha) - 1
  logdensity <- log1p(alpha) -
                (1 + alpha) * (log(x1) + log(x2)) - 
                (2 + 1 / alpha) * log(abs(A))  # Avoid warning

  out.square <- (x1 < 0) | (x1 > 1) | (x2 < 0) | (x2 > 1)
  logdensity[out.square] <- log(0.0)


  index0 <- (rep(alpha, length = length(A)) < sqrt(.Machine$double.eps))
  if (any(index0))
    logdensity[index0] <- log(1.0)


  index1 <- (rep(alpha, length = length(A)) < 0.0) | (A < 0.0)
  if (any(index1))
    logdensity[index1] <- NaN








  if (log.arg) logdensity else exp(logdensity)
}



rbiclaytoncop <- function(n, alpha = 0) {
  if (any(alpha < 0))
    stop("argument 'alpha' must be greater or equal to 0")

  u1 <- runif(n = n)
  v2 <- runif(n = n)

  u2 <- (u1^(-alpha) *
        (v2^(-alpha / (1 + alpha)) - 1) + 1)^(-1 / alpha)


  index0 <- (rep(alpha, length = length(u1)) < sqrt(.Machine$double.eps))
  if (any(index0))
    u2[index0] <- runif(sum(index0))

  cbind(u1, u2)
}



 biclaytoncop <- function(lalpha    = "loge",
                          ialpha    = NULL,
                          imethod   = 1,
                          parallel  = FALSE,
                          zero = NULL) {
  
  apply.parint <- TRUE


  lalpha <- as.list(substitute(lalpha))
  ealpha <- link2list(lalpha)
  lalpha <- attr(ealpha, "function.name")


  if (length(ialpha) && any(ialpha <= 0))
    stop("argument 'ialpha' must have values in (0, Inf)")



  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) || imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")

  new("vglmff",
  blurb = c(" bivariate clayton copula distribution)\n","Links:    ",
                namesof("alpha", lalpha, earg = ealpha)),

  constraints = eval(substitute(expression({
        constraints <- cm.vgam(matrix(1, M, 1), x = x,
                               bool = .parallel , 
                               constraints = constraints,
                               apply.int = .apply.parint )

        dotzero <- .zero
        Musual <- 1
        Yusual <- 2
        eval(negzero.expression)
  }), list( .zero = zero,
            .apply.parint = apply.parint,
            .parallel = parallel ))),

  infos = eval(substitute(function(...) {
        list(Musual = 1,
             Yusual = 2,
             apply.parint = .apply.parint,
             parallel = .parallel,
             zero = .zero )
      }, list( .zero = zero,
               .apply.parint = apply.parint, 
               .parallel = parallel ))),

  initialize = eval(substitute(expression({
    Musual <- 1
    Yusual <- 2

    temp5 <-
      w.y.check(w = w, y = y,
                Is.positive.y = TRUE,
                ncol.w.max = Inf,
                ncol.y.max = Inf,
                ncol.y.min = Yusual,
                out.wy = TRUE,
                colsyperw = Yusual,
                maximize = TRUE)

    w <- temp5$w
    y <- temp5$y


    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    extra$Yusual <- Yusual
    M <- Musual * (ncoly / Yusual)
    mynames1 <- paste("alpha", if (M / Musual > 1) 1:(M / Musual) else "",
                      sep = "")
    predictors.names <- c(
      namesof(mynames1, .lalpha , earg = .ealpha , short = TRUE))


    extra$dimnamesy1 <- dimnames(y)[[1]]
    if (length(dimnames(y)))
      extra$dimnamesy2 <- dimnames(y)[[2]]
    
    if (!length(etastart)) {
      
      alpha.init <- matrix(if (length( .ialpha )) .ialpha else 0 + NA,
                           n, M / Musual, byrow = TRUE)

      if (!length( .ialpha ))
        for (spp. in 1:(M / Musual)) {
          ymatj <- y[, (Yusual * spp. - 1):(Yusual * spp.)]

              
              
          alpha.init0 <- if ( .imethod == 1) {
            k.tau <- kendall.tau(ymatj[, 1], ymatj[, 2], exact = FALSE,
                                 max.n = 500)

            max(0.1, 2 * k.tau / (1 - k.tau))  # Must be positive
          } else if ( .imethod == 2) {
            spearman.rho <-  max(0.05, cor(ymatj[, 1],
                                           ymatj[, 2], meth = "spearman"))
            rhobit(spearman.rho)
          } else {
            pearson.rho <- max(0.05, cor(ymatj[, 1], ymatj[, 2]))
            rhobit(pearson.rho)
          }




          if (any(is.na(alpha.init[, spp.])))
            alpha.init[, spp.] <- alpha.init0
        }
          
      etastart <- theta2eta(alpha.init, .lalpha , earg = .ealpha )
    }
  }), list( .imethod = imethod,
                .lalpha = lalpha,
                .ealpha = ealpha,
                .ialpha = ialpha ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {

        eta <- as.matrix(eta)
        fv.matrix <- matrix(0.5, nrow(eta), extra$ncoly)
        
        
        if (length(extra$dimnamesy2))
          dimnames(fv.matrix) <- list(extra$dimnamesy1,
                                      extra$dimnamesy2)
        fv.matrix
      }  , list( .lalpha = lalpha,
                 .ealpha = ealpha ))),

  last = eval(substitute(expression({
        
    Musual <- extra$Musual
    Yusual <- extra$Yusual
    misc$link <- rep( .lalpha , length = M)
    temp.names <- mynames1
    names(misc$link) <- temp.names
    
    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:M) {
      misc$earg[[ii]] <- .ealpha
    }

    misc$Musual <- Musual
    misc$Yusual <- Yusual
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$parallel  <- .parallel
    misc$apply.parint <- .apply.parint
    misc$multipleResponses <- TRUE

  }) , list( .imethod = imethod,
             .parallel = parallel, .apply.parint = apply.parint,
             .lalpha = lalpha,
             .ealpha = ealpha ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      Alpha <- eta2theta(eta, .lalpha , earg = .ealpha )

      if (residuals) stop("loglikelihood residuals not ", 
                          "implemented yet") else {

      sum(c(w) * dbiclaytoncop(x1  = c(y[, c(TRUE, FALSE)]),
                               x2  = c(y[, c(FALSE, TRUE)]),
                               alpha = c(Alpha), log = TRUE))
      }
    } , list( .lalpha = lalpha,
              .ealpha = ealpha,
              .imethod = imethod ))),
  vfamily = c("biclaytoncop"),
  deriv = eval(substitute(expression({
    Alpha <- eta2theta(eta, .lalpha , earg = .ealpha )
    Yindex1 <- extra$Yusual * (1:(extra$ncoly/extra$Yusual)) - 1
    Yindex2 <- extra$Yusual * (1:(extra$ncoly/extra$Yusual))




    
    AA <- y[, Yindex1]^(-Alpha) + y[, Yindex2]^(-Alpha) - 1
    dAA.dalpha <- -y[, Yindex1]^(-Alpha) * log(y[, Yindex1]) -
                   y[, Yindex2]^(-Alpha) * log(y[, Yindex2])
    dl.dalpha <- 1 / (1 + Alpha) - log(y[, Yindex1] * y[, Yindex2]) -
                 dAA.dalpha / AA * (2 + 1 / Alpha ) + log(AA) / Alpha^2
   


    dalpha.deta <- dtheta.deta(Alpha, .lalpha , earg = .ealpha )

    dl.deta <- c(w) * cbind(dl.dalpha) * dalpha.deta
    dl.deta
  }), list( .lalpha = lalpha,
            .ealpha = ealpha,
            .imethod = imethod ))),

  weight = eval(substitute(expression({


    par <- Alpha +1 #20130808
    denom1 <- (3 * par -2) * (2 * par - 1)
    denom2 <- 2 * (par - 1)
    v1 <- trigamma(1 / (denom2))
    v2 <- trigamma(par / (denom2))
    v3 <- trigamma((2 * par - 1) / (denom2))
    Rho. <- 1 / denom1 * (1 + par / (denom2) * (v1 - v2) +
            1 / (denom2) * (v2 - v3))
    
    out <- 1 / par^2 + 2 / (par * (par - 1) * (2 * par - 1)) +
           4 * par / (3 * par - 2) - 2 * (2 * par - 1) * Rho. / (par - 1)
    ned2l.dalpha  <- out


    wz <- ned2l.dalpha * dalpha.deta^2
 c(w) * wz
  }), list( .lalpha = lalpha,
            .ealpha = ealpha,
            .imethod = imethod ))))
}

















dbistudentt <- function(x1, x2, df, rho = 0, log = FALSE) {




  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  logdensity <-
    -(df/2 + 1) * log1p(
    (x1^2 + x2^2 - 2 * rho * x1 * x2) / (df * (1 - rho^2))) -
    log(2 * pi) - 0.5 * log1p(-rho^2)  # -

  logdensity[df <= 0] <- NaN  # Not picked up by dt().

  if (log.arg) logdensity else exp(logdensity)
}




if (FALSE)
bistudent.deriv.dof <-  function(u, v, nu, rho) {

  
  t1 <- qt(u, nu, 1, 0)
  t2 <- qt(v, nu, 1, 0)
  t3 <- -(nu + 2.0) / 2.0
  t10 <- nu * (1.0 - rho * rho)
  t4 <- -2.0 * t1 * t2 / t10
  t11 <- (t1 * t1 + t2 * t2 - 2.0 * rho * t1 * t2)
  t5 <- 2.0 * t11 * rho / t10 / (1.0 - rho * rho)
  t6 <- 1.0 + (t11 / t10)
  t7 <- rho / (1.0 - rho * rho)
  out <- (t3 * (t4 + t5) / t6  +  t7)
}







 bistudentt <-
   function(ldf     = "loglog",
            lrho    = "rhobit",
            idf     = NULL,
            irho    = NULL,
            imethod = 1,
            parallel = FALSE,
            zero = -1) {




  apply.parint <- TRUE

  ldof <- as.list(substitute(ldf))
  edof <- link2list(ldof)
  ldof <- attr(edof, "function.name")

  lrho <- as.list(substitute(lrho))
  erho <- link2list(lrho)
  lrho <- attr(erho, "function.name")


  idof <- idf
  if (length(idof) &&
      any(idof <= 1))
    stop("argument 'idf' must have values in (1,Inf)")


  if (length(irho) &&
      any(abs(irho) >= 1))
    stop("argument 'irho' must have values in (-1,1)")



  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  new("vglmff",
  blurb = c("Bivariate student-t distribution\n",
            "Links:    ",
            namesof("df",  ldof, earg = edof), ", ",
            namesof("rho", lrho, earg = erho)),

  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(1, M, 1), x = x,
                           bool = .parallel , 
                           constraints = constraints,
                           apply.int = .apply.parint )

    dotzero <- .zero
    Musual <- 2
    Yusual <- 2
    eval(negzero.expression)
  }), list( .zero = zero,
            .apply.parint = apply.parint,
            .parallel = parallel ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
         Yusual = 2,
         apply.parint = .apply.parint ,
         parallel = .parallel ,
         zero = .zero )
  }, list( .zero = zero,
           .apply.parint = apply.parint, 
           .parallel = parallel ))),

  initialize = eval(substitute(expression({
    Musual <- 2
    Yusual <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              ncol.y.min = Yusual,
              out.wy = TRUE,
              colsyperw = Yusual,
              maximize = TRUE)

    w <- temp5$w
    y <- temp5$y


    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    extra$Yusual <- Yusual
    M <- Musual * (ncoly / Yusual)
    mynames1 <- paste("df",  if (M / Musual > 1) 1:(M / Musual) else "",
                      sep = "")
    mynames2 <- paste("rho", if (M / Musual > 1) 1:(M / Musual) else "",
                      sep = "")
    predictors.names <- c(
      namesof(mynames1, .ldof , earg = .edof , short = TRUE),
      namesof(mynames2, .lrho , earg = .erho , short = TRUE))[
              interleave.VGAM(M, M = Musual)]


    extra$dimnamesy1 <- dimnames(y)[[1]]
    if (length(dimnames(y)))
      extra$dimnamesy2 <- dimnames(y)[[2]]

    if (!length(etastart)) {

      dof.init <- matrix(if (length( .idof )) .idof else 0 + NA,
                         n, M / Musual, byrow = TRUE)
      rho.init <- matrix(if (length( .irho )) .irho else 0 + NA,
                         n, M / Musual, byrow = TRUE)

      if (!length( .idof ) || !length( .irho ))
      for (spp. in 1:(M / Musual)) {
        ymatj <- y[, (Musual * spp. - 1):(Musual * spp.)]


        dof.init0 <- if ( .imethod == 1) {


          2 + rexp(n = 1, rate = 0.1)
        } else {
          10
        }

        if (any(is.na(dof.init[, spp.])))
          dof.init[, spp.] <- dof.init0


        rho.init0 <- if ( .imethod == 2) {
          runif(n, min = -1 + 0.1, max = 1 - 0.1)
        } else {
          cor(ymatj[, 1], ymatj[, 2])
        }

        if (any(is.na(rho.init[, spp.])))
          rho.init[, spp.] <- rho.init0

      }

      etastart <-
        cbind(theta2eta(dof.init, .ldof , earg = .edof ),
              theta2eta(rho.init, .lrho , earg = .erho ))

      etastart <- etastart[, interleave.VGAM(M, M = Musual)]

    }
  }), list( .imethod = imethod,
            .lrho = lrho, .ldof = ldof,
            .erho = erho, .edof = edof,
            .idof = idof, .irho = irho ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {

    eta <- as.matrix(eta)
    fv.matrix <- matrix(0.0, nrow(eta), extra$ncoly)


    if (length(extra$dimnamesy2))
      dimnames(fv.matrix) <- list(extra$dimnamesy1,
                                  extra$dimnamesy2)
    fv.matrix
  }  , list( .lrho = lrho, .ldof = ldof,
             .erho = erho, .edof = edof ))),

  last = eval(substitute(expression({

    Musual <- extra$Musual
    Yusual <- extra$Yusual
    misc$link <-
      c(rep( .ldof , length = M / Musual),
        rep( .lrho , length = M / Musual))[
                       interleave.VGAM(M, M = Musual)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = Musual)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:(M / Musual)) {
      misc$earg[[Musual*ii-1]] <- .edof
      misc$earg[[Musual*ii  ]] <- .erho
    }

    misc$Musual <- Musual
    misc$Yusual <- Yusual
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$parallel  <- .parallel
    misc$apply.parint <- .apply.parint
    misc$multipleResponses <- TRUE

  }) , list( .imethod = imethod,
             .parallel = parallel,
             .apply.parint = apply.parint,
             .lrho = lrho, .ldof = ldof,
             .erho = erho, .edof = edof ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Dof <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                     .ldof , earg = .edof )
    Rho <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                     .lrho , earg = .erho )

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      Yindex1 <- extra$Yusual * (1:(extra$ncoly/extra$Yusual)) - 1
      Yindex2 <- extra$Yusual * (1:(extra$ncoly/extra$Yusual))
      sum(c(w) * dbistudentt(x1  = y[, Yindex1, drop = FALSE],
                             x2  = y[, Yindex2, drop = FALSE],
                             df  = Dof,
                             rho = Rho, log = TRUE))
    }
  } , list( .lrho = lrho, .ldof = ldof,
            .erho = erho, .edof = edof,
            .imethod = imethod ))),
  vfamily = c("bistudentt"),
  deriv = eval(substitute(expression({
    Musual <- Yusual <- 2
    Dof <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                     .ldof , earg = .edof )
    Rho <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                     .lrho , earg = .erho )
    Yindex1 <- extra$Yusual * (1:(extra$ncoly/extra$Yusual)) - 1
    Yindex2 <- extra$Yusual * (1:(extra$ncoly/extra$Yusual))


    x1 <- c(y[, Yindex1])  # Convert into a vector
    x2 <- c(y[, Yindex2])

    dee3 <- deriv3( ~
        -(Dof/2 + 1) * log(1 +
        (x1^2 + x2^2 - 2 * Rho * x1 * x2) / (Dof * (1 - Rho^2))) -
        log(2 * pi) - 0.5 * log(1 - Rho^2),
        namevec = c("Dof", "Rho"), hessian = FALSE)
    eval.d3 <- eval(dee3)

    dl.dthetas <-  attr(eval.d3, "gradient")
   
    dl.ddof <- matrix(dl.dthetas[, "Dof"], n, length(Yindex1))
    dl.drho <- matrix(dl.dthetas[, "Rho"], n, length(Yindex2))

  
  if (FALSE) {
    dd <- cbind(y, Rho, Dof)
    pp <- apply(dd, 1, function(x)
                BiCopPDF(x[1], x[2], family = 2, x[3], x[4]))
    alt.dl.ddof <- apply(dd, 1, function(x)
                     BiCopDeriv(x[1], x[2], family = 2,
                                x[3], x[4], "par2")) / pp
    alt.dl.drho <- apply(dd, 1, function(x)
                     BiCopDeriv(x[1], x[2], family = 2,
                                x[3], x[4], "par")) / pp

 print("head(dl.ddof)")
 print( head(dl.ddof) )
 print("head(alt.dl.ddof)")
 print( head(alt.dl.ddof) )

 print("max(abs(alt.dl.drho - dl.drho))")
 print( max(abs(alt.dl.drho - dl.drho)) )
 print("max(abs(alt.dl.ddof - dl.ddof))")
 print( max(abs(alt.dl.ddof - dl.ddof)) )
    
  }





    ddof.deta <- dtheta.deta(Dof, .ldof , earg = .edof )
    drho.deta <- dtheta.deta(Rho, .lrho , earg = .erho )

    ans <- c(w) * cbind(dl.ddof * ddof.deta,
                        dl.drho * drho.deta)
    ans <- ans[, interleave.VGAM(M, M = Musual)]
    ans
  }), list( .lrho = lrho, .ldof = ldof,
            .erho = erho, .edof = edof,
            .imethod = imethod ))),

  weight = eval(substitute(expression({
    wz11 <- beta(2, Dof / 2) / Dof -
            beta(3, Dof / 2) * (Dof + 2) / (4 * Dof)
    wz12 <- -Rho / (2 * (1 - Rho^2)) * (beta(2, Dof / 2) -
            beta(3, Dof / 2) * (Dof + 2) / 2)
    wz22 <- (1 + Rho^2) / (1 - Rho^2)^2 +
            (Dof^2 + 2 * Dof) * Rho^2 *
             beta(3, Dof / 2) / (4 * (1 - Rho^2)^2)
    wz22 <- wz22 + (Dof^2 + 2 * Dof) * (2 - 3 * Rho^2 + Rho^6) *   
            beta(3, Dof / 2) / (16 * (1 - Rho^2)^4)
    wz22 <- wz22 + (Dof^2 + 2 * Dof) * (1 + Rho^2) *    # Replace - by + ???
            beta(2, Dof / 2) / (4 * (1 - Rho^2)^2) # denom == 4 or 2 ???
    ned2l.ddof2   <- wz11
    ned2l.ddofrho <- wz12
    ned2l.drho2   <- wz22

    wz <- array(c(c(w) * ned2l.ddof2 * ddof.deta^2,
                  c(w) * ned2l.drho2 * drho.deta^2,
                  c(w) * ned2l.ddofrho * ddof.deta * drho.deta),
                dim = c(n, M / Musual, 3))
    wz <- arwz2wz(wz, M = M, Musual = Musual)
    wz
  }), list( .lrho = lrho, .ldof = ldof,
            .erho = erho, .edof = edof,
            .imethod = imethod ))))
}




  


dbinormcop <- function(x1, x2, rho = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  x1 <- qnorm(x1)
  x2 <- qnorm(x2)

  logdensity <- (2 * rho * x1 * x2 -
                 rho^2 * (x1^2 + x2^2)) / (2 * (1 - rho^2)) -
                0.5 * log1p(-rho^2)

  if (log.arg) logdensity else exp(logdensity)
}




pbinormcop <- function(q1, q2, rho = 0) {

  if (!is.Numeric(q1, positive = TRUE) ||
      any(q1 >= 1))
    stop("bad input for argument 'q1'")
  if (!is.Numeric(q2, positive = TRUE) ||
      any(q2 >= 1))
    stop("bad input for argument 'q2'")
  if (!is.Numeric(rho) ||
      any(abs(rho) >= 1))
    stop("bad input for argument 'rho'")

  pnorm2(x1 = qnorm(q1),
         x2 = qnorm(q2),
         cov12 = rho)
}


rbinormcop <- function(n, rho = 0) {

  ymat <- rbinorm(n = n, cov12 = rho)
  cbind(y1 = pnorm(ymat[, 1]),
        y2 = pnorm(ymat[, 2]))
}





 binormalcop <- function(lrho    = "rhobit",
                         irho    = NULL,
                         imethod = 1,
                         parallel = FALSE,
                         zero = NULL) {



  apply.parint <- TRUE


  lrho <- as.list(substitute(lrho))
  erho <- link2list(lrho)
  lrho <- attr(erho, "function.name")


  if (length(irho) &&
      any(abs(irho) >= 1))
    stop("argument 'irho' must have values in (-1,1)")



  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")

  new("vglmff",
  blurb = c("Gaussian copula (based on the bivariate normal distribution)\n",
            "Links:    ",
            namesof("rho", lrho, earg = erho)),

  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(1, M, 1), x = x,
                           bool = .parallel , 
                           constraints = constraints,
                           apply.int = .apply.parint )

    dotzero <- .zero
    Musual <- 1
    Yusual <- 2
    eval(negzero.expression)
  }), list( .zero = zero,
            .apply.parint = apply.parint,
            .parallel = parallel ))),

  infos = eval(substitute(function(...) {
    list(Musual = 1,
         Yusual = 2,
         apply.parint = .apply.parint ,
         parallel = .parallel ,
         zero = .zero )
  }, list( .zero = zero,
           .apply.parint = apply.parint, 
           .parallel = parallel ))),

  initialize = eval(substitute(expression({
    Musual <- 1
    Yusual <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              ncol.y.min = Yusual,
              out.wy = TRUE,
              colsyperw = Yusual,
              maximize = TRUE)

    w <- temp5$w
    y <- temp5$y


    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    extra$Yusual <- Yusual
    M <- Musual * (ncoly / Yusual)
    mynames1 <- paste("rho", if (M / Musual > 1) 1:(M / Musual) else "",
                      sep = "")
    predictors.names <- c(
      namesof(mynames1, .lrho , earg = .erho , short = TRUE))


    extra$dimnamesy1 <- dimnames(y)[[1]]
    if (length(dimnames(y)))
      extra$dimnamesy2 <- dimnames(y)[[2]]

    if (!length(etastart)) {

      rho.init <- matrix(if (length( .irho )) .irho else 0 + NA,
                         n, M / Musual, byrow = TRUE)

      if (!length( .irho ))
      for (spp. in 1:(M / Musual)) {
        ymatj <- y[, (Yusual * spp. - 1):(Yusual * spp.)]


        rho.init0 <- if ( .imethod == 1) {
          sin(kendall.tau(ymatj[, 1], ymatj[, 2],
                          exact = FALSE,
                          max.n = 200) * pi / 2)
        } else if ( .imethod == 2) {
          sin(cor(ymatj[, 1], ymatj[, 2],
                  method = "spearman") * pi / 6) * 2
        } else {
          cor(ymatj[, 1], ymatj[, 2])
        }





        if (any(is.na(rho.init[, spp.])))
          rho.init[, spp.] <- rho.init0
      }

      etastart <- theta2eta(rho.init, .lrho , earg = .erho )
    }
  }), list( .imethod = imethod,
            .lrho = lrho,
            .erho = erho,
            .irho = irho ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {

    eta <- as.matrix(eta)
    fv.matrix <- matrix(0.5, nrow(eta), extra$ncoly)


    if (length(extra$dimnamesy2))
      dimnames(fv.matrix) <- list(extra$dimnamesy1,
                                  extra$dimnamesy2)
    fv.matrix
  }  , list( .lrho = lrho,
             .erho = erho ))),

  last = eval(substitute(expression({

    Musual <- extra$Musual
    Yusual <- extra$Yusual
    misc$link <- rep( .lrho , length = M)
    temp.names <- mynames1
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:M) {
      misc$earg[[ii]] <- .erho
    }

    misc$Musual <- Musual
    misc$Yusual <- Yusual
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$parallel  <- .parallel
    misc$apply.parint <- .apply.parint
    misc$multipleResponses <- TRUE

  }) , list( .imethod = imethod,
             .parallel = parallel,
             .apply.parint = apply.parint,
             .lrho = lrho,
             .erho = erho ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Rho <- eta2theta(eta, .lrho , earg = .erho )

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      Yindex1 <- extra$Yusual * (1:(extra$ncoly/extra$Yusual)) - 1
      Yindex2 <- extra$Yusual * (1:(extra$ncoly/extra$Yusual))
      sum(c(w) * dbinormcop(x1  = y[, Yindex1, drop = FALSE],
                            x2  = y[, Yindex2, drop = FALSE],
                            rho = Rho, log = TRUE))
    }
  } , list( .lrho = lrho,
            .erho = erho,
            .imethod = imethod ))),
  vfamily = c("binormalcop"),
  deriv = eval(substitute(expression({
    Rho <- eta2theta(eta, .lrho , earg = .erho )
    Yindex1 <- extra$Yusual * (1:(extra$ncoly/extra$Yusual)) - 1
    Yindex2 <- extra$Yusual * (1:(extra$ncoly/extra$Yusual))

    temp7 <- 1 - Rho^2
    q.y <- qnorm(y)

    dl.drho <- ((1 + Rho^2) * q.y[, Yindex1] * q.y[, Yindex2] -
                Rho * (q.y[, Yindex1]^2 + q.y[, Yindex2]^2)) / temp7^2 +
                Rho / temp7

    drho.deta <- dtheta.deta(Rho, .lrho , earg = .erho )

    c(w) * cbind(dl.drho) * drho.deta
  }), list( .lrho = lrho,
            .erho = erho,
            .imethod = imethod ))),

  weight = eval(substitute(expression({
    ned2l.drho  <- (1 + Rho^2) / temp7^2
    wz <- ned2l.drho * drho.deta^2
    c(w) * wz
  }), list( .lrho = lrho,
            .erho = erho,
            .imethod = imethod ))))
}






bilogistic4.control <- function(save.weight = TRUE, ...) {
  list(save.weight = save.weight)
}


 bilogistic4 <- function(llocation = "identity",
                         lscale = "loge",
                         iloc1 = NULL, iscale1 = NULL,
                         iloc2 = NULL, iscale2 = NULL,
                         imethod = 1, zero = NULL) {

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")




  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2) stop("argument 'imethod' must be 1 or 2")

  new("vglmff",
  blurb = c("Bivariate logistic distribution\n\n",
            "Link:    ",
            namesof("location1", llocat, elocat), ", ",
            namesof("scale1",    lscale, escale), ", ",
            namesof("location2", llocat, elocat), ", ",
            namesof("scale2",    lscale, escale),
            "\n", "\n",
            "Means:     location1, location2"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      c(namesof("location1", .llocat, .elocat , tag = FALSE),
        namesof("scale1",    .lscale, .escale , tag = FALSE),
        namesof("location2", .llocat, .elocat , tag = FALSE),
        namesof("scale2",    .lscale, .escale , tag = FALSE))

    if (!length(etastart)) {
      if ( .imethod == 1) {
        locat.init1 <- y[, 1]
        scale.init1 <- sqrt(3) * sd(y[, 1]) / pi
        locat.init2 <- y[, 2]
        scale.init2 <- sqrt(3) * sd(y[, 2]) / pi
      } else {
        locat.init1 <- median(rep(y[, 1], w))
        locat.init2 <- median(rep(y[, 2], w))
        const4 <- sqrt(3) / (sum(w) * pi)
        scale.init1 <- const4 * sum(c(w) *(y[, 1] - locat.init1)^2)
        scale.init2 <- const4 * sum(c(w) *(y[, 2] - locat.init2)^2)
      }
      loc1.init <- if (length( .iloc1 ))
                   rep( .iloc1, length.out = n) else
                   rep(locat.init1, length.out = n)
      loc2.init <- if (length( .iloc2 ))
                   rep( .iloc2, length.out = n) else
                   rep(locat.init2, length.out = n)
      scale1.init <- if (length( .iscale1 ))
                     rep( .iscale1, length.out = n) else
                     rep(1, length.out = n)
      scale2.init <- if (length( .iscale2 ))
                     rep( .iscale2, length.out = n) else
                     rep(1, length.out = n)

      if ( .llocat == "loge")
        locat.init1 <- abs(locat.init1) + 0.001
      if ( .llocat == "loge")
        locat.init2 <- abs(locat.init2) + 0.001

      etastart <-
        cbind(theta2eta(locat.init1, .llocat , .elocat ),
              theta2eta(scale1.init, .lscale , .escale ),
              theta2eta(locat.init2, .llocat , .elocat ),
              theta2eta(scale2.init, .lscale , .escale ))
    }
  }), list(.imethod = imethod,
           .iloc1 = iloc1, .iloc2 = iloc2,
           .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale,
           .iscale1 = iscale1, .iscale2 = iscale2))),
  linkinv = function(eta, extra = NULL) {
    cbind(eta[, 1], eta[, 2])
  },
  last = eval(substitute(expression({
    misc$link <-    c(location1 = .llocat, scale1 = .lscale,
                      location2 = .llocat, scale2 = .lscale)

    misc$earg <- list(location1 = .elocat, scale1 = .escale,
                      location2 = .elocat, scale2 = .escale)

    misc$expected <- FALSE
    misc$BFGS <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    locat1 <- eta2theta(eta[, 1], .llocat , .elocat )
    Scale1 <- eta2theta(eta[, 2], .lscale , .escale )
    locat2 <- eta2theta(eta[, 3], .llocat , .elocat )
    Scale2 <- eta2theta(eta[, 4], .lscale , .escale )

    zedd1 <- (y[, 1]-locat1) / Scale1
    zedd2 <- (y[, 2]-locat2) / Scale2

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
    sum(c(w) * (-zedd1 - zedd2 -
                3 * log1p(exp(-zedd1) + exp(-zedd2)) -
                log(Scale1) - log(Scale2)))
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale ))),
  vfamily = c("bilogistic4"),
  deriv = eval(substitute(expression({
    locat1 <- eta2theta(eta[, 1], .llocat , .elocat )
    Scale1 <- eta2theta(eta[, 2], .lscale , .escale )
    locat2 <- eta2theta(eta[, 3], .llocat , .elocat )
    Scale2 <- eta2theta(eta[, 4], .lscale , .escale )

    zedd1 <- (y[, 1]-locat1) / Scale1
    zedd2 <- (y[, 2]-locat2) / Scale2
    ezedd1 <- exp(-zedd1)
    ezedd2 <- exp(-zedd2)
    denom <- 1 + ezedd1 + ezedd2

    dl.dlocat1 <- (1 - 3 * ezedd1 / denom) / Scale1
    dl.dlocat2 <- (1 - 3 * ezedd2 / denom) / Scale2
    dl.dscale1 <- (zedd1 - 1 - 3 * ezedd1 * zedd1 / denom) / Scale1
    dl.dscale2 <- (zedd2 - 1 - 3 * ezedd2 * zedd2 / denom) / Scale2

    dlocat1.deta <- dtheta.deta(locat1, .llocat , .elocat )
    dlocat2.deta <- dtheta.deta(locat2, .llocat , .elocat )
    dscale1.deta <- dtheta.deta(Scale1, .lscale , .escale )
    dscale2.deta <- dtheta.deta(Scale2, .lscale , .escale )

    if (iter == 1) {
        etanew <- eta
    } else {
        derivold <- derivnew
        etaold <- etanew
        etanew <- eta
    }
    derivnew <- c(w) * cbind(dl.dlocat1 * dlocat1.deta,
                             dl.dscale1 * dscale1.deta,
                             dl.dlocat2 * dlocat2.deta,
                             dl.dscale2 * dscale2.deta)
    derivnew
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale ))),
  weight = eval(substitute(expression({
    if (iter == 1) {
      wznew <- cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
    } else {
      wzold <- wznew
      wznew <- qnupdate(w = w, wzold=wzold, dderiv=(derivold - derivnew),
                       deta=etanew-etaold, M = M,
                       trace=trace)  # weights incorporated in args
    }
    wznew
  }), list( .lscale = lscale,
            .escale = escale,
            .llocat = llocat))))
}






dbilogis4 <- function(x1, x2, loc1 = 0, scale1 = 1,
                      loc2 = 0, scale2 = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)




    L <- max(length(x1), length(x2),
             length(loc1), length(loc2),
             length(scale1), length(scale2))
    if (length(x1    ) != L) x1     <- rep(x1,     length.out = L)
    if (length(x2    ) != L) x2     <- rep(x2,     length.out = L)
    if (length(loc1  ) != L) loc1   <- rep(loc1,   length.out = L)
    if (length(loc2  ) != L) loc2   <- rep(loc2,   length.out = L)
    if (length(scale1) != L) scale1 <- rep(scale1, length.out = L)
    if (length(scale2) != L) scale2 <- rep(scale2, length.out = L)
    zedd1 <- (-(x1 - loc1) / scale1)
    zedd2 <- (-(x2 - loc2) / scale2)
    logdensity <- log(2) + log(zedd1) + log(zedd2) - log(scale1) - 
                  log(scale1) - 3 * log1p(exp(zedd1) + exp(zedd2))
    if (log.arg) logdensity else exp(logdensity)
}



pbilogis4 <-
  function(q1, q2, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1) {

  ans <- 1 / (1 + exp(-(q1-loc1)/scale1) + exp(-(q2-loc2)/scale2))
  ans[scale1 <= 0] <- NA
  ans[scale2 <= 0] <- NA
  ans
}



rbilogis4 <- function(n, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1) {
  y1 <- rlogis(n = n, location = loc1, scale = scale1)
  ezedd1 <- exp(-(y1-loc1)/scale1)
  y2 <- loc2 - scale2 * log(1/sqrt(runif(n) / (1 + ezedd1)^2) - 1 - ezedd1)
  ans <- cbind(y1, y2)
  ans[scale2 <= 0, ] <- NA
  ans
}



 freund61 <- function(la  = "loge",
                      lap = "loge",
                      lb  = "loge",
                      lbp = "loge",
                      ia = NULL, iap = NULL, ib = NULL, ibp = NULL,
                      independent = FALSE,
                      zero = NULL) {
  la <- as.list(substitute(la))
  ea <- link2list(la)
  la <- attr(ea, "function.name")

  lap <- as.list(substitute(lap))
  eap <- link2list(lap)
  lap <- attr(eap, "function.name")

  lb <- as.list(substitute(lb))
  eb <- link2list(lb)
  lb <- attr(eb, "function.name")


  lbp <- as.list(substitute(lbp))
  ebp <- link2list(lbp)
  lbp <- attr(ebp, "function.name")



  new("vglmff",
  blurb = c("Freund (1961) bivariate exponential distribution\n",
            "Links:    ",
            namesof("a",  la,  earg = ea ), ", ",
            namesof("ap", lap, earg = eap), ", ",
            namesof("b",  lb,  earg = eb ), ", ",
            namesof("bp", lbp, earg = ebp)),
  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(c(1, 1,0,0, 0,0, 1, 1), M, 2), x = x,
                           bool = .independent ,
                           constraints = constraints,
                           apply.int = TRUE)
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list(.independent = independent, .zero = zero))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      c(namesof("a",  .la,  earg = .ea , short = TRUE), 
        namesof("ap", .lap, earg = .eap, short = TRUE), 
        namesof("b",  .lb,  earg = .eb , short = TRUE), 
        namesof("bp", .lbp, earg = .ebp, short = TRUE))
    extra$y1.lt.y2 = y[, 1] < y[, 2]

    if (!(arr <- sum(extra$y1.lt.y2)) || arr == n)
        stop("identifiability problem: either all y1<y2 or y2<y1")

    if (!length(etastart)) {
      sumx  <- sum(y[ extra$y1.lt.y2, 1]);
      sumxp <- sum(y[!extra$y1.lt.y2, 1])
      sumy  <- sum(y[ extra$y1.lt.y2, 2]);
      sumyp <- sum(y[!extra$y1.lt.y2, 2])

      if (FALSE) { # Noise:
        arr <- min(arr + n/10, n*0.95)
        sumx <- sumx * 1.1; sumxp <- sumxp * 1.2;
        sumy <- sumy * 1.2; sumyp <- sumyp * 1.3;
      }
      ainit  <- if (length(.ia))  rep(.ia, length.out = n) else
         arr / (sumx + sumyp)
      apinit <- if (length(.iap)) rep(.iap,length.out = n) else
         (n-arr)/(sumxp-sumyp)
      binit  <- if (length(.ib))  rep(.ib, length.out = n) else
         (n-arr)/(sumx +sumyp)
      bpinit <- if (length(.ib))  rep(.ibp,length.out = n) else
         arr / (sumy - sumx)

      etastart <-
        cbind(theta2eta(rep(ainit,  length.out = n), .la,  earg = .ea  ),
              theta2eta(rep(apinit, length.out = n), .lap, earg = .eap ),
              theta2eta(rep(binit,  length.out = n), .lb,  earg = .eb  ),
              theta2eta(rep(bpinit, length.out = n), .lbp, earg = .ebp ))
    }
  }), list( .la = la, .lap = lap, .lb = lb, .lbp = lbp,
            .ea = ea, .eap = eap, .eb = eb, .ebp = ebp,
            .ia = ia, .iap = iap, .ib = ib, .ibp = ibp))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha  <- eta2theta(eta[, 1], .la,  earg = .ea  )
    alphap <- eta2theta(eta[, 2], .lap, earg = .eap )
    beta   <- eta2theta(eta[, 3], .lb,  earg = .eb  )
    betap  <- eta2theta(eta[, 4], .lbp, earg = .ebp )
    cbind((alphap + beta) / (alphap * (alpha + beta)),
          (alpha + betap) / (betap * (alpha + beta)))
  }, list( .la = la, .lap = lap, .lb = lb, .lbp = lbp,
           .ea = ea, .eap = eap, .eb = eb, .ebp = ebp ))),
  last = eval(substitute(expression({
    misc$link <-    c("a" = .la, "ap" = .lap, "b" = .lb, "bp" = .lbp)

    misc$earg <- list("a" = .ea, "ap" = .eap, "b" = .eb, "bp" = .ebp)

    misc$multipleResponses <- FALSE
  }), list( .la = la, .lap = lap, .lb = lb, .lbp = lbp,
            .ea = ea, .eap = eap, .eb = eb, .ebp = ebp ))),
  loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    alpha  <- eta2theta(eta[, 1], .la,  earg = .ea  )
    alphap <- eta2theta(eta[, 2], .lap, earg = .eap )
    beta   <- eta2theta(eta[, 3], .lb,  earg = .eb  )
    betap  <- eta2theta(eta[, 4], .lbp, earg = .ebp )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      tmp88 <- extra$y1.lt.y2
      ell1 <- log(alpha[tmp88]) + log(betap[tmp88]) -
             betap[tmp88] * y[tmp88, 2] -
             (alpha+beta-betap)[tmp88] * y[tmp88, 1]
      ell2 <- log(beta[!tmp88]) + log(alphap[!tmp88]) -
             alphap[!tmp88] * y[!tmp88, 1] -
             (alpha+beta-alphap)[!tmp88] * y[!tmp88, 2]
    sum(w[tmp88] * ell1) + sum(w[!tmp88] * ell2) }
  }, list( .la = la, .lap = lap, .lb = lb, .lbp = lbp,
           .ea = ea, .eap = eap, .eb = eb, .ebp = ebp ))),
  vfamily = c("freund61"),
  deriv = eval(substitute(expression({
    tmp88  <- extra$y1.lt.y2
    alpha  <- eta2theta(eta[, 1], .la,  earg = .ea  )
    alphap <- eta2theta(eta[, 2], .lap, earg = .eap )
    beta   <- eta2theta(eta[, 3], .lb,  earg = .eb  )
    betap  <- eta2theta(eta[, 4], .lbp, earg = .ebp )

    dalpha.deta  <- dtheta.deta(alpha,  .la,  earg = .ea  )
    dalphap.deta <- dtheta.deta(alphap, .lap, earg = .eap )
    dbeta.deta   <- dtheta.deta(beta,   .lb,  earg = .eb  )
    dbetap.deta  <- dtheta.deta(betap,  .lbp, earg = .ebp )

    d1 <- 1/alpha - y[, 1]
    d1[!tmp88] <- -y[!tmp88, 2]
    d2 <- 0 * alphap
    d2[!tmp88] <- 1/alphap[!tmp88] - y[!tmp88, 1] + y[!tmp88, 2]
    d3 <- -y[, 1]
    d3[!tmp88] <- 1/beta[!tmp88] - y[!tmp88, 2]
    d4 <- 1/betap - y[, 2] + y[, 1]
    d4[!tmp88] <- 0

    c(w) * cbind(d1 * dalpha.deta,
                 d2 * dalphap.deta,
                 d3 * dbeta.deta,
                 d4 * dbetap.deta)
  }), list( .la = la, .lap = lap, .lb = lb, .lbp = lbp,
            .ea = ea, .eap = eap, .eb = eb, .ebp = ebp ))),
  weight = eval(substitute(expression({
    py1.lt.y2 <- alpha / (alpha+beta)
    d11 <- py1.lt.y2 / alpha^2
    d22 <- (1-py1.lt.y2) / alphap^2
    d33 <- (1-py1.lt.y2) / beta^2
    d44 <- py1.lt.y2 / betap^2

    wz <- matrix(0, n, M) # diagonal
    wz[, iam(1, 1, M)] <- dalpha.deta^2  * d11
    wz[, iam(2, 2, M)] <- dalphap.deta^2 * d22
    wz[, iam(3, 3, M)] <- dbeta.deta^2   * d33
    wz[, iam(4, 4, M)] <- dbetap.deta^2  * d44

    c(w) * wz
  }), list( .la = la, .lap = lap, .lb = lb, .lbp = lbp,
            .ea = ea, .eap = eap, .eb = eb, .ebp = ebp ))))
}








 bigamma.mckay <- function(lscale = "loge",
                           lshape1 = "loge",
                           lshape2 = "loge",
                           iscale = NULL,
                           ishape1 = NULL,
                           ishape2 = NULL,
                           imethod = 1,
                           zero = 1) {
  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")

  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")


  if (!is.null(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
      stop("argument 'iscale' must be positive or NULL")
  if (!is.null(ishape1))
    if (!is.Numeric(ishape1, positive = TRUE))
      stop("argument 'ishape1' must be positive or NULL")
  if (!is.null(ishape2))
    if (!is.Numeric(ishape2, positive = TRUE))
      stop("argument 'ishape2' must be positive or NULL")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2.5)
    stop("argument 'imethod' must be 1 or 2")



  new("vglmff",
  blurb = c("Bivariate gamma: McKay's distribution\n",
            "Links:    ",
            namesof("scale",  lscale), ", ",
            namesof("shape1", lshape1), ", ",
            namesof("shape2", lshape2)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    if (any(y[, 1] >= y[, 2]))
      stop("the second column minus the first column must be a vector ",
           "of positive values")


    predictors.names <-
      c(namesof("scale",  .lscale,  .escale,  short = TRUE), 
        namesof("shape1", .lshape1, .eshape1, short = TRUE), 
        namesof("shape2", .lshape2, .eshape2, short = TRUE))

    if (!length(etastart)) {
      momentsY <- if ( .imethod == 1) {
        cbind(median(y[, 1]),  # This may not be monotonic
              median(y[, 2])) + 0.01
      } else {
        cbind(weighted.mean(y[, 1], w),
              weighted.mean(y[, 2], w))
      }

      mcg2.loglik <- function(thetaval, y, x, w, extraargs) {
        ainit <- a <- thetaval
        momentsY <- extraargs$momentsY
          p <- (1/a) * abs(momentsY[1]) + 0.01
          q <- (1/a) * abs(momentsY[2] - momentsY[1]) + 0.01
            sum(c(w) * (-(p+q)*log(a) - lgamma(p) - lgamma(q) +
                 (p - 1)*log(y[, 1]) +
                 (q - 1)*log(y[, 2]-y[, 1]) - y[, 2] / a ))
      }

      a.grid <- if (length( .iscale )) c( .iscale ) else
         c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
      extraargs <- list(momentsY = momentsY)
      ainit <- getMaxMin(a.grid, objfun = mcg2.loglik,
                         y = y, x = x, w = w, maximize = TRUE,
                         extraargs = extraargs)
      ainit <- rep(if(is.Numeric( .iscale )) .iscale else ainit,
                   length.out = n)
      pinit <- (1/ainit) * abs(momentsY[1]) + 0.01
      qinit <- (1/ainit) * abs(momentsY[2] - momentsY[1]) + 0.01

      pinit <- rep(if(is.Numeric( .ishape1 )) .ishape1 else pinit,
                   length.out = n)
      qinit <- rep(if(is.Numeric( .ishape2 )) .ishape2 else qinit,
                   length.out = n)

      etastart <-
        cbind(theta2eta(ainit, .lscale),
              theta2eta(pinit, .lshape1),
              theta2eta(qinit, .lshape2))
    }
  }), list( .lscale = lscale, .lshape1 = lshape1, .lshape2 = lshape2,
            .escale = escale, .eshape1 = eshape1, .eshape2 = eshape2,
            .iscale = iscale, .ishape1 = ishape1, .ishape2 = ishape2,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    a <- eta2theta(eta[, 1], .lscale  ,  .escale )
    p <- eta2theta(eta[, 2], .lshape1 , .eshape1 )
    q <- eta2theta(eta[, 3], .lshape2 , .eshape2 )
    cbind("y1" = p*a,
          "y2" = (p+q)*a)
  }, list( .lscale = lscale, .lshape1 = lshape1, .lshape2 = lshape2,
           .escale = escale, .eshape1 = eshape1, .eshape2 = eshape2 ))),
  last = eval(substitute(expression({
    misc$link <-    c("scale"  = .lscale ,
                      "shape1" = .lshape1 ,
                      "shape2" = .lshape2 )

    misc$earg <- list("scale"  = .escale ,
                      "shape1" = .eshape1 ,
                      "shape2" = .eshape2 )

    misc$ishape1 <- .ishape1
    misc$ishape2 <- .ishape2
    misc$iscale <- .iscale
    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lscale = lscale, .lshape1 = lshape1, .lshape2 = lshape2,
            .escale = escale, .eshape1 = eshape1, .eshape2 = eshape2,
            .iscale = iscale, .ishape1 = ishape1, .ishape2 = ishape2,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    a <- eta2theta(eta[, 1], .lscale  ,  .escale )
    p <- eta2theta(eta[, 2], .lshape1 , .eshape1 )
    q <- eta2theta(eta[, 3], .lshape2 , .eshape2 )

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
    sum(c(w) * (-(p+q)*log(a) - lgamma(p) - lgamma(q) +
               (p - 1)*log(y[, 1]) + (q - 1)*log(y[, 2]-y[, 1]) -
               y[, 2] / a))
  }, list( .lscale = lscale, .lshape1 = lshape1, .lshape2 = lshape2,
           .escale = escale, .eshape1 = eshape1, .eshape2 = eshape2 ))),
  vfamily = c("bigamma.mckay"),
  deriv = eval(substitute(expression({
    aparam <- eta2theta(eta[, 1], .lscale  ,  .escale )
    shape1 <- eta2theta(eta[, 2], .lshape1 , .eshape1 )
    shape2 <- eta2theta(eta[, 3], .lshape2 , .eshape2 )

    dl.da <- (-(shape1+shape2) + y[, 2] / aparam) / aparam
    dl.dshape1 <- -log(aparam) - digamma(shape1) + log(y[, 1])
    dl.dshape2 <- -log(aparam) - digamma(shape2) + log(y[, 2]-y[, 1])

    c(w) * cbind(dl.da      * dtheta.deta(aparam, .lscale),
                 dl.dshape1 * dtheta.deta(shape1, .lshape1),
                 dl.dshape2 * dtheta.deta(shape2, .lshape2))
  }), list( .lscale = lscale, .lshape1 = lshape1, .lshape2 = lshape2,
            .escale = escale, .eshape1 = eshape1, .eshape2 = eshape2 ))),
  weight = eval(substitute(expression({
    d11 <- (shape1+shape2) / aparam^2
    d22 <- trigamma(shape1)
    d33 <- trigamma(shape2)
    d12 <- 1 / aparam
    d13 <- 1 / aparam
    d23 <- 0

    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M)] <- dtheta.deta(aparam, .lscale)^2 * d11
    wz[, iam(2, 2, M)] <- dtheta.deta(shape1, .lshape1)^2 * d22
    wz[, iam(3, 3, M)] <- dtheta.deta(shape2, .lshape2)^2 * d33
    wz[, iam(1, 2, M)] <- dtheta.deta(aparam, .lscale) *
                          dtheta.deta(shape1, .lshape1) * d12
    wz[, iam(1, 3, M)] <- dtheta.deta(aparam, .lscale) *
                          dtheta.deta(shape2, .lshape2) * d13
    wz[, iam(2, 3, M)] <- dtheta.deta(shape1, .lshape1) *
                          dtheta.deta(shape2, .lshape2) * d23

    c(w) * wz
  }), list( .lscale = lscale, .lshape1 = lshape1,
                              .lshape2 = lshape2 ))))
}











rbifrankcop <- function(n, alpha) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n
  if (!is.Numeric(alpha, positive = TRUE))
    stop("bad input for argument 'alpha'")
  alpha <- rep(alpha, length.out = use.n)
  U <- runif(use.n)
  V <- runif(use.n)

  T <- alpha^U + (alpha - alpha^U) * V
  X <- U
  index <- abs(alpha - 1) < .Machine$double.eps
  Y <- U
  if (any(!index))
    Y[!index] <- logb(T[!index] / (T[!index] +
                      (1 - alpha[!index]) * V[!index]),
                      base = alpha[!index])
  ans <- matrix(c(X, Y), nrow = use.n, ncol = 2)
  if (any(index)) {
    ans[index, 1] <- runif(sum(index))  # Uniform density for alpha == 1
    ans[index, 2] <- runif(sum(index))
  }
  ans
}


pbifrankcop <- function(q1, q2, alpha) {
  if (!is.Numeric(q1))                     stop("bad input for 'q1'")
  if (!is.Numeric(q2))                     stop("bad input for 'q2'")
  if (!is.Numeric(alpha, positive = TRUE)) stop("bad input for 'alpha'")

  L <- max(length(q1), length(q2), length(alpha))
  if (length(alpha) != L) alpha <- rep(alpha, length.out = L)
  if (length(q1   ) != L) q1    <- rep(q1,    length.out = L)
  if (length(q2   ) != L) q2    <- rep(q2,    length.out = L)

  x <- q1; y <- q2
  index <- (x >= 1 & y <  1) | (y >= 1 & x <  1) |
           (x <= 0 | y <= 0) | (x >= 1 & y >= 1) |
           (abs(alpha - 1) < .Machine$double.eps)
  ans <- as.numeric(index)
  if (any(!index))
  ans[!index] <- logb(1 + ((alpha[!index])^(x[!index]) - 1)*
                 ((alpha[!index])^(y[!index]) - 1)/(alpha[!index] - 1), 
                 base = alpha[!index])
  ind2 <- (abs(alpha - 1) < .Machine$double.eps)
  ans[ind2] <- x[ind2] * y[ind2]
  ans[x >= 1 & y <  1] <- y[x >= 1 & y < 1]  # P(Y2 < q2) = q2
  ans[y >= 1 & x <  1] <- x[y >= 1 & x < 1]  # P(Y1 < q1) = q1
  ans[x <= 0 | y <= 0] <- 0
  ans[x >= 1 & y >= 1] <- 1
  ans
}


dbifrankcop <- function(x1, x2, alpha, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if (!is.Numeric(x1))                     stop("bad input for 'x1'")
  if (!is.Numeric(x2))                     stop("bad input for 'x2'")
  if (!is.Numeric(alpha, positive = TRUE)) stop("bad input for 'alpha'")

  L <- max(length(x1), length(x2), length(alpha))
  if (length(alpha) != L) alpha <- rep(alpha, length.out = L)
  if (length(x1   ) != L) x1    <- rep(x1,    length.out = L)
  if (length(x2   ) != L) x2    <- rep(x2,    length.out = L)

  if (log.arg) {
    denom <- alpha-1 + (alpha^x1  - 1) * (alpha^x2  - 1)
    denom <- abs(denom)
    log((alpha - 1) * log(alpha)) + (x1+x2)*log(alpha) - 2 * log(denom)
  } else {
    temp <- (alpha - 1) + (alpha^x1 - 1) * (alpha^x2 - 1)
    index <- (abs(alpha - 1) < .Machine$double.eps)
    ans <- x1
    if (any(!index))
      ans[!index] <- (alpha[!index] - 1) * log(alpha[!index]) *
                     (alpha[!index])^(x1[!index] +
                                      x2[!index]) / (temp[!index])^2
    ans[x1 <= 0 | x2 <= 0 | x1 >= 1 | x2 >= 1] <- 0
    ans[index] <- 1
    ans
  }
}




bifrankcop.control <- function(save.weight = TRUE, ...) {
  list(save.weight = save.weight)
}





 bifrankcop <- function(lapar = "loge", iapar = 2, nsimEIM = 250) {

  lapar <- as.list(substitute(lapar))
  eapar <- link2list(lapar)
  lapar <- attr(eapar, "function.name")


  if (!is.Numeric(iapar, positive = TRUE))
    stop("argument 'iapar' must be positive")


  if (length(nsimEIM) &&
     (!is.Numeric(nsimEIM, allowable.length = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50))
    stop("argument 'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Frank's bivariate copula\n",
            "Links:    ",
            namesof("apar", lapar, earg = eapar )),
  initialize = eval(substitute(expression({

    if (any(y <= 0) || any(y >= 1))
      stop("the response must have values between 0 and 1") 

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      c(namesof("apar", .lapar , earg = .eapar, short = TRUE))

    if (length(dimnames(y)))
      extra$dimnamesy2 <- dimnames(y)[[2]]

    if (!length(etastart)) {
      apar.init <- rep(.iapar, length.out = n)
      etastart <- cbind(theta2eta(apar.init, .lapar , earg = .eapar ))
    }
  }), list( .lapar = lapar, .eapar = eapar, .iapar = iapar))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    apar <- eta2theta(eta, .lapar , earg = .eapar )
    fv.matrix <- matrix(0.5, length(apar), 2)
    if (length(extra$dimnamesy2))
      dimnames(fv.matrix) <- list(names(eta), extra$dimnamesy2)
    fv.matrix
  }, list( .lapar = lapar, .eapar = eapar ))),
  last = eval(substitute(expression({
    misc$link <-    c("apar" = .lapar )

    misc$earg <- list("apar" = .eapar )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$pooled.weight <- pooled.weight
    misc$multipleResponses <- FALSE
  }), list( .lapar = lapar, .eapar = eapar, .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    apar <- eta2theta(eta, .lapar , earg = .eapar )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
        sum(c(w) * dbifrankcop(x1 = y[, 1], x2 = y[, 2],
                               alpha = apar, log = TRUE))
    }
  }, list( .lapar = lapar, .eapar = eapar ))),
  vfamily = c("bifrankcop"),
  deriv = eval(substitute(expression({
    apar <- eta2theta(eta, .lapar , earg = .eapar )
    dapar.deta <- dtheta.deta(apar, .lapar , earg = .eapar )

    de3 <- deriv3(~ (log((apar - 1) * log(apar)) + (y1+y2)*log(apar) -
                      2 * log(apar-1 + (apar^y1  - 1) * (apar^y2  - 1))),
                    name = "apar", hessian = TRUE)

    denom <- apar-1 + (apar^y[, 1]  - 1) * (apar^y[, 2]  - 1)
    tmp700 <- 2*apar^(y[, 1]+y[, 2]) - apar^y[, 1] - apar^y[, 2]
    numerator <- 1 + y[, 1] * apar^(y[, 1] - 1) * (apar^y[, 2]  - 1) + 
                     y[, 2] * apar^(y[, 2] - 1) * (apar^y[, 1]  - 1)
    Dl.dapar <- 1/(apar - 1) + 1/(apar*log(apar)) +
                (y[, 1]+y[, 2])/apar - 2 * numerator / denom
    c(w) * Dl.dapar * dapar.deta
  }), list( .lapar = lapar,
            .eapar = eapar, .nsimEIM = nsimEIM ))),
  weight = eval(substitute(expression({
  if ( is.Numeric( .nsimEIM)) {

    pooled.weight <- FALSE  # For @last


    run.mean <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- rbifrankcop(n, alpha = apar)
        y1 <- ysim[, 1]; y2 <- ysim[, 2];
        eval.de3 <- eval(de3)
        d2l.dthetas2 <-  attr(eval.de3, "hessian")
        rm(ysim)
        temp3 <- -d2l.dthetas2[, 1, 1]   # M = 1
        run.mean <- ((ii - 1) * run.mean + temp3) / ii
    }
    wz <- if (intercept.only)
        matrix(mean(run.mean), n, dimm(M)) else run.mean

    wz <- wz * dapar.deta^2
    c(w) * wz
  } else {
      nump <- apar^(y[, 1]+y[, 2]-2) * (2 * y[, 1] * y[, 2] +
                    y[, 1]*(y[, 1] - 1) + y[, 2]*(y[, 2] - 1)) - 
                    y[, 1]*(y[, 1] - 1) * apar^(y[, 1]-2) - 
                    y[, 2]*(y[, 2] - 1) * apar^(y[, 2]-2)
      D2l.dapar2 <- 1/(apar - 1)^2 + (1+log(apar))/(apar*log(apar))^2 +
                    (y[, 1]+y[, 2])/apar^2 + 2 *
                    (nump / denom - (numerator/denom)^2)
      d2apar.deta2 <- d2theta.deta2(apar, .lapar , earg = .eapar )
      wz <- c(w) * (dapar.deta^2 * D2l.dapar2 - Dl.dapar * d2apar.deta2)
      if (TRUE && intercept.only) {
        wz <- cbind(wz)
        sumw <- sum(w)
        for (iii in 1:ncol(wz))
          wz[,iii] <- sum(wz[, iii]) / sumw
        pooled.weight <- TRUE
        wz <- c(w) * wz   # Put back the weights
      } else {
        pooled.weight <- FALSE
      }
    wz
  }
  }), list( .lapar = lapar,
            .eapar = eapar, .nsimEIM = nsimEIM ))))
}





 gammahyp <- function(ltheta = "loge", itheta = NULL, expected = FALSE) {

  ltheta <- as.list(substitute(ltheta))
  etheta <- link2list(ltheta)
  ltheta <- attr(etheta, "function.name")

  if (!is.logical(expected) || length(expected) != 1)
      stop("argument 'expected' must be a single logical")


  new("vglmff",
  blurb = c("Gamma hyperbola bivariate distribution\n",
            "Links:    ",
            namesof("theta", ltheta, etheta)),
  initialize = eval(substitute(expression({
    if (any(y[, 1] <= 0) || any(y[, 2] <= 1))
      stop("the response has values that are out of range") 

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      c(namesof("theta", .ltheta , .etheta , short = TRUE))

    if (!length(etastart)) {
      theta.init <- if (length( .itheta)) {
        rep( .itheta , length.out = n) 
      } else {
        1 / (y[, 2] - 1 + 0.01)
      }
      etastart <-
        cbind(theta2eta(theta.init, .ltheta , .etheta ))
    }
  }), list( .ltheta = ltheta, .etheta = etheta, .itheta = itheta))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    theta <- eta2theta(eta, .ltheta , .etheta )
    cbind(theta*exp(theta), 1+1/theta)
  }, list( .ltheta = ltheta, .etheta = etheta ))),
  last = eval(substitute(expression({
    misc$link <-    c("theta" = .ltheta )

    misc$earg <- list("theta" = .etheta )

    misc$expected <- .expected 
    misc$multipleResponses <- FALSE
  }), list( .ltheta = ltheta,
            .etheta = etheta, .expected = expected ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    theta <- eta2theta(eta, .ltheta , .etheta )
    if (residuals) stop("loglikelihood residuals not ",
                      "implemented yet") else {
      sum(c(w) * (-exp(-theta) * y[, 1] / theta - theta * y[, 2]))
    }
  }, list( .ltheta = ltheta, .etheta = etheta ))),
  vfamily = c("gammahyp"),
  deriv = eval(substitute(expression({
    theta <- eta2theta(eta, .ltheta , .etheta )
    Dl.dtheta <- exp(-theta) * y[, 1] * (1+theta) / theta^2 - y[, 2]
    DTHETA.deta <- dtheta.deta(theta, .ltheta , .etheta )
    c(w) * Dl.dtheta * DTHETA.deta
  }), list( .ltheta = ltheta, .etheta = etheta ))),
  weight = eval(substitute(expression({
    temp300 <- 2 + theta * (2 + theta)
    if ( .expected ) {
      D2l.dtheta2 <- temp300 / theta^2
      wz <- c(w) * DTHETA.deta^2 * D2l.dtheta2
    } else {
      D2l.dtheta2 <- temp300 * y[, 1] * exp(-theta) / theta^3
      D2theta.deta2 <- d2theta.deta2(theta, .ltheta )
      wz <- c(w) * (DTHETA.deta^2 * D2l.dtheta2 -
                    Dl.dtheta * D2theta.deta2)
    }
    wz
  }), list( .ltheta = ltheta,
            .etheta = etheta, .expected = expected ))))
}



 morgenstern <- function(lapar = "rhobit",
                         iapar = NULL, tola0 = 0.01,
                         imethod = 1) {
  lapar <- as.list(substitute(lapar))
  earg  <- link2list(lapar)
  lapar <- attr(earg, "function.name")

  if (length(iapar) &&
     (!is.Numeric(iapar, allowable.length = 1) ||
      abs(iapar) >= 1))
    stop("argument 'iapar' must be a single number between -1 and 1")

  if (!is.Numeric(tola0, allowable.length = 1, positive = TRUE))
      stop("argument 'tola0' must be a single positive number")

  if (length(iapar) && abs(iapar) <= tola0)
      stop("argument 'iapar' must not be between -tola0 and tola0")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2.5)
      stop("argument 'imethod' must be 1 or 2")


  new("vglmff",
  blurb = c("Morgenstern's bivariate exponential distribution\n",
            "Links:    ",
            namesof("apar", lapar, earg = earg )),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <-
      c(namesof("apar", .lapar , earg = .earg , short = TRUE))

    if (length(dimnames(y)))
      extra$dimnamesy2 = dimnames(y)[[2]]

    if (!length(etastart)) {
      ainit  <- if (length(.iapar))  rep( .iapar , length.out = n) else {
        mean1 <- if ( .imethod == 1) median(y[, 1]) else mean(y[, 1])
        mean2 <- if ( .imethod == 1) median(y[, 2]) else mean(y[, 2])
        Finit <- 0.01 + mean(y[, 1] <= mean1 & y[, 2] <= mean2)
            ((Finit+expm1(-mean1)+exp(-mean2)) / exp(-mean1-mean2) - 1) / (
             expm1(-mean1) * expm1(-mean2))
          }
        etastart <-
          theta2eta(rep(ainit, length.out = n), .lapar , earg = .earg )
      }
  }), list( .iapar = iapar, .lapar = lapar, .earg = earg,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha <- eta2theta(eta, .lapar , earg = .earg )
    fv.matrix <- matrix(1, length(alpha), 2)
    if (length(extra$dimnamesy2))
        dimnames(fv.matrix) = list(names(eta), extra$dimnamesy2)
    fv.matrix
  }, list( .lapar = lapar, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c("apar" = .lapar )

    misc$earg <- list("apar" = .earg  )

    misc$expected <- FALSE
    misc$pooled.weight <- pooled.weight
    misc$multipleResponses <- FALSE
  }), list( .lapar = lapar, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      alpha  <- eta2theta(eta, .lapar , earg = .earg )
      alpha[abs(alpha) < .tola0 ] <- .tola0
      if (residuals) stop("loglikelihood residuals not ",
                          "implemented yet") else {
      denom <- (1 + alpha - 2*alpha*(exp(-y[, 1]) + exp(-y[, 2])) +
               4*alpha*exp(-y[, 1] - y[, 2]))
      sum(c(w) * (-y[, 1] - y[, 2] + log(denom)))
    }
  }, list( .lapar = lapar, .earg = earg, .tola0=tola0 ))),
  vfamily = c("morgenstern"),
  deriv = eval(substitute(expression({
    alpha  <- eta2theta(eta, .lapar , earg = .earg )
    alpha[abs(alpha) < .tola0 ] <- .tola0
    numerator <- 1 - 2*(exp(-y[, 1]) + exp(-y[, 2])) +
                 4*exp(-y[, 1] - y[, 2])
    denom <- (1 + alpha - 2*alpha*(exp(-y[, 1]) + exp(-y[, 2])) +
             4 *alpha*exp(-y[, 1] - y[, 2]))
    dl.dalpha <- numerator / denom

    dalpha.deta <- dtheta.deta(alpha,  .lapar , earg = .earg )

    c(w) * cbind(dl.dalpha * dalpha.deta)
  }), list( .lapar = lapar, .earg = earg, .tola0=tola0 ))),
  weight = eval(substitute(expression({
    d2l.dalpha2 <- dl.dalpha^2
    d2alpha.deta2 <- d2theta.deta2(alpha,  .lapar , earg = .earg )
    wz <- c(w) * (dalpha.deta^2 * d2l.dalpha2 - d2alpha.deta2 * dl.dalpha)
    if (TRUE &&
       intercept.only) {
        wz <- cbind(wz)
      sumw <- sum(w)
      for (iii in 1:ncol(wz))
        wz[,iii] <- sum(wz[, iii]) / sumw
      pooled.weight <- TRUE
      wz <- c(w) * wz   # Put back the weights
    } else {
      pooled.weight <- FALSE
    }
    wz
  }), list( .lapar = lapar, .earg = earg ))))
}




rfgm <- function(n, alpha) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (!is.Numeric(alpha))
    stop("bad input for argument 'alpha'")
  if (any(abs(alpha) > 1))
    stop("argument 'alpha' has values out of range")

  y1 <- V1 <- runif(use.n)
  V2 <- runif(use.n)
  temp <- 2*y1 - 1
  A <- alpha * temp - 1
  B <- sqrt(1 - 2 * alpha * temp + (alpha*temp)^2 + 4 * alpha * V2 * temp)
  y2 <- 2 * V2 / (B - A)
  matrix(c(y1, y2), nrow = use.n, ncol = 2)
}



dfgm <- function(x1, x2, alpha, log = FALSE) {
  if (!is.logical(log.arg <- log) ||
      length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(alpha))
    stop("bad input for 'alpha'")
  if (any(abs(alpha) > 1))
    stop("'alpha' values out of range")
  if ( !is.logical( log.arg ) ||
       length( log.arg ) != 1 )
    stop("bad input for argument 'log'")

  L <- max(length(x1), length(x2), length(alpha))
  if (length(x1)    != L)  x1    <- rep(x1,    length.out = L)
  if (length(x2)    != L)  x2    <- rep(x2,    length.out = L)
  if (length(alpha) != L)  alpha <- rep(alpha, length.out = L)
  ans <- 0 * x1
  xnok <- (x1 <= 0) | (x1 >= 1) | (x2 <= 0) | (x2 >= 1)
  if ( log.arg ) {
    ans[!xnok] <- log1p(alpha[!xnok] * (1-2*x1[!xnok]) * (1-2*x2[!xnok]))
    ans[xnok] <- log(0)
  } else {
    ans[!xnok] <-   1 + alpha[!xnok] * (1-2*x1[!xnok]) * (1-2*x2[!xnok])
    ans[xnok] <- 0
    if (any(ans < 0))
      stop("negative values in the density (alpha out of range)")
  }
  ans
}


pfgm <- function(q1, q2, alpha) {
  if (!is.Numeric(q1))     stop("bad input for 'q1'")
  if (!is.Numeric(q2))     stop("bad input for 'q2'")
  if (!is.Numeric(alpha))  stop("bad input for 'alpha'")
  if (any(abs(alpha) > 1)) stop("'alpha' values out of range")

  L <- max(length(q1), length(q2), length(alpha))
  if (length(q1)    != L)     q1 <- rep(q1,    length.out = L)
  if (length(q2)    != L)     q2 <- rep(q2,    length.out = L)
  if (length(alpha) != L)  alpha <- rep(alpha, length.out = L)

  x <- q1
  y <- q2
  index <- (x >= 1 & y <  1) |
           (y >= 1 & x <  1) |
           (x <= 0 | y <= 0) |
           (x >= 1 & y >= 1)
  ans <- as.numeric(index)
  if (any(!index)) {
    ans[!index] <-    q1[!index] *   q2[!index] * (1 + alpha[!index] *
                   (1-q1[!index])*(1-q2[!index]))
  }
  ans[x >= 1 & y<1] <- y[x >= 1 & y<1]  # P(Y2 < q2) = q2
  ans[y >= 1 & x<1] <- x[y >= 1 & x<1]  # P(Y1 < q1) = q1
  ans[x <= 0 | y <= 0] <- 0
  ans[x >= 1 & y >= 1] <- 1
  ans
}






 fgm <- function(lapar = "rhobit", iapar = NULL,
                 imethod = 1) {

  lapar <- as.list(substitute(lapar))
  earg  <- link2list(lapar)
  lapar <- attr(earg, "function.name")


  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3.5)
    stop("argument 'imethod' must be 1 or 2 or 3")

  if (length(iapar) &&
     (abs(iapar) >= 1))
    stop("'iapar' should be less than 1 in absolute value")


  new("vglmff",
  blurb = c("Farlie-Gumbel-Morgenstern distribution\n",
            "Links:    ",
            namesof("apar", lapar, earg = earg )),
  initialize = eval(substitute(expression({
    if (any(y < 0) || any(y > 1))
      stop("the response must have values in the unit square")

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      namesof("apar", .lapar , earg = .earg , short = TRUE)

    if (length(dimnames(y)))
      extra$dimnamesy2 <- dimnames(y)[[2]]

    if (!length(etastart)) {
      ainit  <- if (length( .iapar ))  .iapar else {


      if ( .imethod == 1) {
        3 * cor(y[, 1], y[, 2], method = "spearman")
      } else if ( .imethod == 2) {
        9 * kendall.tau(y[, 1], y[, 2]) / 2
      } else {
        mean1 <- if ( .imethod == 1) weighted.mean(y[, 1], w) else
                 median(y[, 1])
        mean2 <- if ( .imethod == 1) weighted.mean(y[, 2], w) else
                 median(y[, 2])
        Finit <- weighted.mean(y[, 1] <= mean1 & y[, 2] <= mean2, w)
        (Finit / (mean1 * mean2) - 1) / ((1 - mean1) * (1 - mean2))
      }
    }

    ainit <- min(0.95, max(ainit, -0.95))

    etastart <-
      theta2eta(rep(ainit, length.out = n), .lapar , earg = .earg )
    }
  }), list( .iapar = iapar, .lapar = lapar, .earg = earg,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha <- eta2theta(eta, .lapar , earg = .earg )
    fv.matrix <- matrix(0.5, length(alpha), 2)
    if (length(extra$dimnamesy2))
      dimnames(fv.matrix) <- list(names(eta), extra$dimnamesy2)
    fv.matrix
  }, list( .lapar = lapar, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c("apar" = .lapar )

    misc$earg <- list("apar" = .earg  )

    misc$expected <- FALSE
    misc$multipleResponses <- FALSE
  }), list( .lapar = lapar, .earg = earg))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    alpha <- eta2theta(eta, .lapar , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dfgm(x1 = y[, 1],
                      x2 = y[, 2], alpha = alpha, log = TRUE))
    }
  }, list( .lapar = lapar, .earg = earg ))),
  vfamily = c("fgm"),
  deriv = eval(substitute(expression({
    alpha  <- eta2theta(eta, .lapar , earg = .earg )

    dalpha.deta <- dtheta.deta(alpha, .lapar , earg = .earg )

    numerator <- (1 - 2 * y[, 1])  * (1 - 2 * y[, 2])
    denom <- 1 + alpha * numerator

    mytolerance <- .Machine$double.eps
    bad <- (denom <= mytolerance)   # Range violation
    if (any(bad)) {
      cat("There are some range violations in @deriv\n")
      flush.console()
      denom[bad] <- 2 * mytolerance
    }
    dl.dalpha <- numerator / denom
    c(w) * cbind(dl.dalpha * dalpha.deta)
  }), list( .lapar = lapar, .earg = earg))),

  weight = eval(substitute(expression({
  wz <- lerch(alpha^2, 2, 1.5) / 4  # Checked and correct
  wz <- wz * dalpha.deta^2
    c(w) * wz
  }), list( .lapar = lapar, .earg = earg))))
}




 bigumbelI <- function(lapar = "identity", iapar = NULL, imethod = 1) {

  lapar <- as.list(substitute(lapar))
  earg  <- link2list(lapar)
  lapar <- attr(earg, "function.name")


  if (length(iapar) &&
      !is.Numeric(iapar, allowable.length = 1))
    stop("'iapar' must be a single number")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2.5)
    stop("argument 'imethod' must be 1 or 2")


  new("vglmff",
  blurb = c("Gumbel's Type I bivariate distribution\n",
            "Links:    ",
            namesof("apar", lapar, earg = earg )),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <-
      c(namesof("apar", .lapar , earg = .earg , short = TRUE))

    if (!length(etastart)) {
      ainit  <- if (length( .iapar ))  rep( .iapar, length.out = n) else {
        mean1 <- if ( .imethod == 1) median(y[, 1]) else mean(y[, 1])
        mean2 <- if ( .imethod == 1) median(y[, 2]) else mean(y[, 2])
        Finit <- 0.01 + mean(y[, 1] <= mean1 & y[, 2] <= mean2)
        (log(Finit+expm1(-mean1)+exp(-mean2))+mean1+mean2)/(mean1*mean2)
      }
      etastart <-
        theta2eta(rep(ainit,  length.out = n), .lapar , earg = .earg )
      }
  }), list( .iapar = iapar, .lapar = lapar, .earg = earg,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha <- eta2theta(eta, .lapar , earg = .earg )
    cbind(rep(1, len = length(alpha)),
          rep(1, len = length(alpha)))
  }, list( .lapar = lapar, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c("apar" = .lapar )

    misc$earg <- list("apar" = .earg  )

    misc$expected <- FALSE
    misc$pooled.weight <- pooled.weight
    misc$multipleResponses <- FALSE
  }), list( .lapar = lapar, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    alpha  <- eta2theta(eta, .lapar , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      denom <- (alpha*y[, 1] - 1) * (alpha*y[, 2] - 1) + alpha
      mytolerance <- .Machine$double.xmin
      bad <- (denom <= mytolerance)  # Range violation
      if (any(bad)) {
        cat("There are some range violations in @deriv\n")
        flush.console()
      }
      sum(bad) * (-1.0e10) + 
      sum(w[!bad] * (-y[!bad, 1] - y[!bad, 2] +
          alpha[!bad] * y[!bad, 1] * y[!bad, 2] + log(denom[!bad])))
    }
  }, list( .lapar = lapar, .earg = earg ))),
  vfamily = c("bigumbelI"),
  deriv = eval(substitute(expression({
    alpha  <- eta2theta(eta, .lapar , earg = .earg )
    numerator <- (alpha * y[, 1] - 1) * y[, 2] +
                 (alpha * y[, 2] - 1) * y[, 1] + 1
    denom <- (alpha * y[, 1] - 1) * (alpha * y[, 2] - 1) + alpha
    denom <- abs(denom)

    dl.dalpha <- numerator / denom + y[, 1] * y[, 2]

    dalpha.deta <- dtheta.deta(alpha,  .lapar , earg = .earg )

    c(w) * cbind(dl.dalpha * dalpha.deta)
  }), list( .lapar = lapar, .earg = earg ))),
  weight = eval(substitute(expression({
    d2l.dalpha2 <- (numerator/denom)^2 - 2*y[, 1]*y[, 2] / denom
    d2alpha.deta2 <- d2theta.deta2(alpha, .lapar , earg = .earg )
    wz <- c(w) * (dalpha.deta^2 * d2l.dalpha2 - d2alpha.deta2 * dl.dalpha)
    if (TRUE &&
           intercept.only) {
            wz <- cbind(wz)
      sumw <- sum(w)
      for (iii in 1:ncol(wz))
        wz[, iii] <- sum(wz[, iii]) / sumw
      pooled.weight <- TRUE
      wz <- c(w) * wz   # Put back the weights
    } else {
      pooled.weight <- FALSE
    }
    wz
  }), list( .lapar = lapar, .earg = earg ))))
}







pplack <- function(q1, q2, oratio) {
  if (!is.Numeric(q1)) stop("bad input for 'q1'")
  if (!is.Numeric(q2)) stop("bad input for 'q2'")
  if (!is.Numeric(oratio, positive = TRUE)) stop("bad input for 'oratio'")

  L <- max(length(q1), length(q2), length(oratio))
  if (length(q1) != L)  q1 <- rep(q1, length.out = L)
  if (length(q2) != L)  q2 <- rep(q2, length.out = L)
  if (length(oratio) != L)  oratio <- rep(oratio, length.out = L)

  x <- q1; y <- q2
  index <- (x >= 1 & y <  1) | (y >= 1 & x <  1) |
           (x <= 0 | y <= 0) | (x >= 1 & y >= 1) |
           (abs(oratio - 1) < 1.0e-6)  #  .Machine$double.eps
  ans <- as.numeric(index)
  if (any(!index)) {
    temp1 <- 1 + (oratio[!index]  - 1) * (q1[!index] + q2[!index])
    temp2 <- temp1 - sqrt(temp1^2 - 4 * oratio[!index] *
             (oratio[!index] - 1) * q1[!index] * q2[!index])
    ans[!index] <- 0.5 * temp2 / (oratio[!index] - 1)
  }

  ind2 <- (abs(oratio - 1) < 1.0e-6) # .Machine$double.eps
  ans[ind2] <- x[ind2] * y[ind2]
  ans[x >= 1 & y<1] <- y[x >= 1 & y<1] # P(Y2 < q2) = q2
  ans[y >= 1 & x<1] <- x[y >= 1 & x<1] # P(Y1 < q1) = q1
  ans[x <= 0 | y <= 0] <- 0
  ans[x >= 1 & y >= 1] <- 1
  ans
}



rplack <- function(n, oratio) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n


  y1 <- U <- runif(use.n)
  V <- runif(use.n)
  Z <- V * (1-V)
  y2 <- (2*Z*(y1*oratio^2 + 1 - y1) + oratio * (1 - 2 * Z) -
        (1 - 2 * V) *
        sqrt(oratio * (oratio + 4*Z*y1*(1-y1)*(1-oratio)^2))) / (oratio +
        Z*(1-oratio)^2)
  matrix(c(y1, 0.5 * y2), nrow = use.n, ncol = 2)
}



dplack <- function(x1, x2, oratio, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  ans <- log(oratio) + log1p((oratio - 1) *
         (x1+x2 - 2*x1*x2)) - 1.5 *
         log((1 + (x1+x2)*(oratio - 1))^2 -
             4 * oratio * (oratio - 1)*x1*x2)
  ans[ # !is.na(x1) & !is.na(x2) & !is.na(oratio) &
     ((x1 < 0) | (x1 > 1) | (x2 < 0) | (x2 > 1))] <- log(0)


  if (log.arg) ans else exp(ans)
}



plackett.control <- function(save.weight = TRUE, ...) {
  list(save.weight = save.weight)
}



 plackett <- function(link = "loge", ioratio = NULL,
                      imethod = 1, nsimEIM = 200) {

  link <- as.list(substitute(link))
  earg  <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(ioratio) && (!is.Numeric(ioratio, positive = TRUE)))
    stop("'ioratio' must be positive")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  new("vglmff",
  blurb = c("Plackett distribution\n",
            "Links:    ",
            namesof("oratio", link, earg = earg )),
  initialize = eval(substitute(expression({
    if (any(y < 0) || any(y > 1))
      stop("the response must have values in the unit square")

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      namesof("oratio", .link , earg = .earg, short = TRUE)

    if (length(dimnames(y)))
      extra$dimnamesy2 <- dimnames(y)[[2]]

    if (!length(etastart)) {
      orinit <- if (length( .ioratio ))  .ioratio else {
          if ( .imethod == 2) {
            scorp <- cor(y)[1, 2]
            if (abs(scorp) <= 0.1) 1 else
            if (abs(scorp) <= 0.3) 3^sign(scorp) else
            if (abs(scorp) <= 0.6) 5^sign(scorp) else
            if (abs(scorp) <= 0.8) 20^sign(scorp) else 40^sign(scorp)
          } else {
            y10 <- weighted.mean(y[, 1], w)
            y20 <- weighted.mean(y[, 2], w)
            (0.5 + sum(w[(y[, 1] <  y10) & (y[, 2] <  y20)])) *
            (0.5 + sum(w[(y[, 1] >= y10) & (y[, 2] >= y20)])) / (
            ((0.5 + sum(w[(y[, 1] <  y10) & (y[, 2] >= y20)])) *
             (0.5 + sum(w[(y[, 1] >= y10) & (y[, 2] <  y20)]))))
          }
        }
        etastart <-
          theta2eta(rep(orinit, length.out = n),
                    .link , earg = .earg )
    }
  }), list( .ioratio = ioratio, .link = link, .earg = earg,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    oratio <- eta2theta(eta, .link , earg = .earg )
    fv.matrix <- matrix(0.5, length(oratio), 2)
    if (length(extra$dimnamesy2))
        dimnames(fv.matrix) <- list(dimnames(eta)[[1]], extra$dimnamesy2)
    fv.matrix
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c(oratio = .link)

    misc$earg <- list(oratio = .earg)

    misc$expected <- FALSE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .link = link, .earg = earg,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    oratio <- eta2theta(eta, .link , earg = .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
        sum(c(w) * dplack(x1 = y[, 1], x2 = y[, 2],
                          oratio = oratio, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("plackett"),
  deriv = eval(substitute(expression({
    oratio  <- eta2theta(eta, .link , earg = .earg )
    doratio.deta <- dtheta.deta(oratio, .link , earg = .earg )
    y1 <- y[, 1]
    y2 <- y[, 2]
    de3 <- deriv3(~ (log(oratio) + log(1+(oratio - 1) *
                 (y1+y2-2*y1*y2)) - 1.5 *
                 log((1 + (y1+y2)*(oratio - 1))^2 -
                 4 * oratio * (oratio - 1)*y1*y2)),
                 name = "oratio", hessian = FALSE)
    eval.de3 <- eval(de3)

    dl.doratio <-  attr(eval.de3, "gradient")

    c(w) * dl.doratio * doratio.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    sd3 <- deriv3(~ (log(oratio) + log(1+(oratio - 1) *
          (y1sim+y2sim-2*y1sim*y2sim)) - 1.5 *
          log((1 + (y1sim+y2sim)*(oratio - 1))^2 -
          4 * oratio * (oratio - 1)*y1sim*y2sim)),
                    name = "oratio", hessian = FALSE)
    run.var <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- rplack(n, oratio=oratio)
      y1sim <- ysim[, 1]
      y2sim <- ysim[, 1]
        eval.sd3 <- eval(sd3)
        dl.doratio <-  attr(eval.sd3, "gradient")
        rm(ysim, y1sim, y2sim)
        temp3 <- dl.doratio
        run.var <- ((ii - 1) * run.var + temp3^2) / ii
    }
    wz <- if (intercept.only)
        matrix(colMeans(cbind(run.var)),
               n, dimm(M), byrow = TRUE) else cbind(run.var)

    wz <- wz * doratio.deta^2
    c(w) * wz
  }), list( .link = link, .earg = earg, .nsimEIM = nsimEIM ))))
}




damh <- function(x1, x2, alpha, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  L <- max(length(x1), length(x2), length(alpha))
  alpha <- rep(alpha,  length.out = L)
  x1    <- rep(x1,     length.out = L)
  x2    <- rep(x2,     length.out = L)
  temp <- 1 - alpha*(1-x1)*(1-x2)

  if (log.arg) {
    ans <- log1p(-alpha+2*alpha*x1*x2/temp) - 2*log(temp)
    ans[(x1 <= 0) | (x1 >= 1) | (x2 <= 0) | (x2 >= 1)] <- log(0)
  } else {
    ans <- (1-alpha+2*alpha*x1*x2/temp) / (temp^2)
    ans[(x1 <= 0) | (x1 >= 1) | (x2 <= 0) | (x2 >= 1)] <- 0
  }
  ans[abs(alpha) > 1] <- NA
  ans
}


pamh <- function(q1, q2, alpha) {
  if (!is.Numeric(q1)) stop("bad input for 'q1'")
  if (!is.Numeric(q2)) stop("bad input for 'q2'")
  if (!is.Numeric(alpha)) stop("bad input for 'alpha'")

  L <- max(length(q1), length(q2), length(alpha))
  if (length(q1) != L)  q1 <- rep(q1, length.out = L)
  if (length(q2) != L)  q2 <- rep(q2, length.out = L)
  if (length(alpha) != L)  alpha <- rep(alpha, length.out = L)

  x <- q1; y <- q2
  index <- (x >= 1 & y < 1) | (y >= 1 & x <  1) |
           (x <= 0 | y<= 0) | (x >= 1 & y >= 1)
  ans <- as.numeric(index)
  if (any(!index)) {
      ans[!index] <- (q1[!index]*q2[!index]) / (1 -
                     alpha[!index]*(1-q1[!index])*(1-q2[!index]))
  }
  ans[x >= 1 & y <  1] <- y[x >= 1 & y < 1] # P(Y2 < q2) = q2
  ans[y >= 1 & x <  1] <- x[y >= 1 & x < 1] # P(Y1 < q1) = q1
  ans[x <= 0 | y <= 0] <- 0
  ans[x >= 1 & y >= 1] <- 1
  ans[abs(alpha) > 1] <- NA
  ans
}


ramh <- function(n, alpha) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (any(abs(alpha) > 1))
    stop("'alpha' values out of range")

  U1 <- V1 <- runif(use.n)
  V2 <- runif(use.n)
  b <- 1-V1
  A <- -alpha*(2*b*V2+1)+2*alpha^2*b^2*V2+1
  B <- alpha^2*(4*b^2*V2-4*b*V2+1)+alpha*(4*V2-4*b*V2-2)+1
  U2 <- (2*V2*(alpha*b - 1)^2)/(A+sqrt(B))
  matrix(c(U1, U2), nrow = use.n, ncol = 2)
}


amh.control <- function(save.weight = TRUE, ...) {
  list(save.weight = save.weight)
}


 amh <- function(lalpha = "rhobit", ialpha = NULL,
                 imethod = 1, nsimEIM = 250) {
  lalpha <- as.list(substitute(lalpha))
  ealpha <- link2list(lalpha)
  lalpha <- attr(ealpha, "function.name")



  if (length(ialpha) && (abs(ialpha) > 1))
    stop("'ialpha' should be less than or equal to 1 in absolute value")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
    imethod > 2)
    stop("imethod must be 1 or 2")

  if (length(nsimEIM) &&
    (!is.Numeric(nsimEIM, allowable.length = 1,
                  integer.valued = TRUE) ||
     nsimEIM <= 50))
  stop("'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Ali-Mikhail-Haq distribution\n",
            "Links:    ",
            namesof("alpha", lalpha, earg = ealpha )),
  initialize = eval(substitute(expression({
    if (any(y < 0) || any(y > 1))
        stop("the response must have values in the unit square")

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      c(namesof("alpha", .lalpha, earg = .ealpha, short = TRUE))

    if (length(dimnames(y)))
      extra$dimnamesy2 <- dimnames(y)[[2]]

    if (!length(etastart)) {
      ainit  <- if (length( .ialpha ))  .ialpha else {
          mean1 <- if ( .imethod == 1) weighted.mean(y[, 1], w) else
                   median(y[, 1])
          mean2 <- if ( .imethod == 1) weighted.mean(y[, 2], w) else
                   median(y[, 2])
          Finit <- weighted.mean(y[, 1] <= mean1 & y[, 2] <= mean2, w)
          (1 - (mean1 * mean2 / Finit)) / ((1-mean1) * (1-mean2))
      }
      ainit <- min(0.95, max(ainit, -0.95))
      etastart <-
        theta2eta(rep(ainit, length.out = n), .lalpha, earg = .ealpha )
    }
  }), list( .lalpha = lalpha, .ealpha = ealpha, .ialpha = ialpha,
            .imethod = imethod))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha <- eta2theta(eta, .lalpha, earg = .ealpha )
    fv.matrix <- matrix(0.5, length(alpha), 2)
    if (length(extra$dimnamesy2))
        dimnames(fv.matrix) <- list(names(eta), extra$dimnamesy2)
    fv.matrix
  }, list( .lalpha = lalpha, .ealpha = ealpha ))),
  last = eval(substitute(expression({
    misc$link <-    c("alpha" = .lalpha )

    misc$earg <- list("alpha" = .ealpha )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lalpha = lalpha,
            .ealpha = ealpha, .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    alpha <- eta2theta(eta, .lalpha, earg = .ealpha )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * damh(x1 = y[, 1], x2 = y[, 2],
                      alpha = alpha, log = TRUE))
    }
  }, list( .lalpha = lalpha, .ealpha = ealpha ))),
  vfamily = c("amh"),
  deriv = eval(substitute(expression({
    alpha <- eta2theta(eta, .lalpha, earg = .ealpha )

    dalpha.deta <- dtheta.deta(alpha, .lalpha, earg = .ealpha )

    y1 <- y[, 1]
    y2 <- y[, 2]
    de3 <- deriv3(~ (log(1 - alpha+
                        (2 * alpha*y1*y2/(1-alpha*(1-y1)*(1-y2)))) -
                    2 * log(1 - alpha*(1-y1)*(1-y2))) ,
                    name = "alpha", hessian = FALSE)
    eval.de3 <- eval(de3)

    dl.dalpha <-  attr(eval.de3, "gradient")

    c(w) * dl.dalpha * dalpha.deta
  }), list( .lalpha = lalpha, .ealpha = ealpha ))),
  weight = eval(substitute(expression({
    sd3 <- deriv3(~ (log(1 - alpha +
                        (2 * alpha * y1sim * y2sim / (1 - alpha *
                         (1 - y1sim) * (1-y2sim)))) -
                     2 * log(1-alpha*(1-y1sim)*(1-y2sim))),
                     name = "alpha", hessian = FALSE)
    run.var <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- ramh(n, alpha = alpha)
      y1sim <- ysim[, 1]
      y2sim <- ysim[, 1]
      eval.sd3 <- eval(sd3)
      dl.alpha <-  attr(eval.sd3, "gradient")
      rm(ysim, y1sim, y2sim)
      temp3 <- dl.dalpha
      run.var <- ((ii - 1) * run.var + temp3^2) / ii
    }

    wz <- if (intercept.only)
        matrix(colMeans(cbind(run.var)),
               n, dimm(M), byrow = TRUE) else cbind(run.var)

    wz <- wz * dalpha.deta^2

    c(w) * wz
  }), list( .lalpha = lalpha,
            .ealpha = ealpha, .nsimEIM = nsimEIM ))))
}















dbinorm <- function(x1, x2, mean1 = 0, mean2 = 0,
                    var1 = 1, var2 = 1, cov12 = 0,
                    log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  sd1 <- sqrt(var1)
  sd2 <- sqrt(var2)
  rho <- cov12 / (sd1 * sd2)


  temp5 <- 1 - rho^2
  zedd1 <- (x1 - mean1) / sd1
  zedd2 <- (x2 - mean2) / sd2
  logpdf <- -log(2 * pi) - log(sd1) - log(sd2) -
              0.5 * log1p(-rho^2) +
            -(0.5 / temp5)  * (zedd1^2 + (-2 * rho * zedd1 + zedd2) * zedd2)
  if (log.arg) logpdf else exp(logpdf)
}



rbinorm <- function(n, mean1 = 0, mean2 = 0,
                    var1 = 1, var2 = 1, cov12 = 0) {

  Y1 <- rnorm(n)
  Y2 <- rnorm(n)
  X1 <- sqrt(var1) * Y1 + mean1
  delta <- sqrt(var2 - (cov12^2) / var1)
  X2 <- cov12 * Y1 / sqrt(var1) + delta * Y2 + mean2

  ans <- cbind(X1, X2)
  ans[is.na(delta), ] <- NA

  ans
}




 binormal <- function(lmean1 = "identity",
                      lmean2 = "identity",
                      lsd1   = "loge",
                      lsd2   = "loge",
                      lrho   = "rhobit",
                      imean1 = NULL,       imean2 = NULL,
                      isd1   = NULL,       isd2   = NULL,
                      irho   = NULL,       imethod = 1,
                      eq.mean = FALSE,     eq.sd = FALSE,
                      zero = 3:5) {

  lmean1 <- as.list(substitute(lmean1))
  emean1 <- link2list(lmean1)
  lmean1 <- attr(emean1, "function.name")

  lmean2 <- as.list(substitute(lmean2))
  emean2 <- link2list(lmean2)
  lmean2 <- attr(emean2, "function.name")

  lsd1 <- as.list(substitute(lsd1))
  esd1 <- link2list(lsd1)
  lsd1 <- attr(esd1, "function.name")

  lsd2 <- as.list(substitute(lsd2))
  esd2 <- link2list(lsd2)
  lsd2 <- attr(esd2, "function.name")

  lrho <- as.list(substitute(lrho))
  erho <- link2list(lrho)
  lrho <- attr(erho, "function.name")




  trivial1 <- is.logical(eq.mean) && length(eq.mean) == 1 && !eq.mean
  trivial2 <- is.logical(eq.sd  ) && length(eq.sd  ) == 1 && !eq.sd
  if (!trivial1 && !trivial2)
    stop("only one of 'eq.mean' and 'eq.sd' can be assigned a value")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  new("vglmff",
  blurb = c("Bivariate normal distribution\n",
            "Links:    ",
            namesof("mean1", lmean1, earg = emean1 ), ", ",
            namesof("mean2", lmean2, earg = emean2 ), ", ",
            namesof("sd1",   lsd1,   earg = esd1   ), ", ",
            namesof("sd2",   lsd2,   earg = esd2   ), ", ",
            namesof("rho",   lrho,   earg = erho   )),
  constraints = eval(substitute(expression({
    temp8.m <- diag(5)[, -2]
    temp8.m[2, 1] <- 1
    temp8.s <- diag(5)[, -4]
    temp8.s[4, 3] <- 1
    constraints <- cm.vgam(temp8.m, x = x,
                           bool = .eq.mean ,
                           constraints = constraints, apply.int = TRUE)
    constraints <- cm.vgam(temp8.s, x = x,
                           bool = .eq.sd ,
                           constraints = constraints, apply.int = TRUE)
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero,
            .eq.sd   = eq.sd,
            .eq.mean = eq.mean ))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <- c(
      namesof("mean1", .lmean1 , earg = .emean1 , short = TRUE),
      namesof("mean2", .lmean2 , earg = .emean2 , short = TRUE),
      namesof("sd1",   .lsd1 ,   earg = .esd1 ,   short = TRUE),
      namesof("sd2",   .lsd2 ,   earg = .esd2 ,   short = TRUE),
      namesof("rho",   .lrho ,   earg = .erho ,   short = TRUE))

    if (length(dimnames(y)))
      extra$dimnamesy2 <- dimnames(y)[[2]]

    if (!length(etastart)) {
      imean1 <- rep(if (length( .imean1 )) .imean1 else
                   weighted.mean(y[, 1], w = w), length.out = n)
      imean2 <- rep(if (length( .imean2 )) .imean2 else
                   weighted.mean(y[, 2], w = w), length.out = n)
      isd1   <- rep(if (length( .isd1 )) .isd1 else  sd(y[, 1]),
                   length.out = n)
      isd2   <- rep(if (length( .isd2 )) .isd2 else  sd(y[, 2]),
                   length.out = n)
      irho   <- rep(if (length( .irho )) .irho else cor(y[, 1], y[, 2]),
                   length.out = n)

      if ( .imethod == 2) {
        imean1 <- abs(imean1) + 0.01
        imean2 <- abs(imean2) + 0.01
      }
      etastart <-
        cbind(theta2eta(imean1, .lmean1 , earg = .emean1 ),
              theta2eta(imean2, .lmean2 , earg = .emean2 ),
              theta2eta(isd1,   .lsd1 ,   earg = .esd1 ),
              theta2eta(isd2,   .lsd2 ,   earg = .esd2 ),
              theta2eta(irho,   .lrho ,   earg = .erho ))
    }
  }), list( .lmean1 = lmean1, .lmean2 = lmean2,
            .emean1 = emean1, .emean2 = emean2,
            .lsd1   = lsd1  , .lsd2   = lsd2  , .lrho = lrho,
            .esd1   = esd1  , .esd2   = esd2  , .erho = erho,
            .imethod = imethod,
            .imean1 = imean1, .imean2 = imean2,
            .isd1   = isd1,   .isd2   = isd2,
            .irho   = irho ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mean1 <- eta2theta(eta[, 1], .lmean1, earg = .emean1)
    mean2 <- eta2theta(eta[, 2], .lmean2, earg = .emean2)
    fv.matrix <- cbind(mean1, mean2)
    if (length(extra$dimnamesy2))
      dimnames(fv.matrix) <- list(names(eta), extra$dimnamesy2)
    fv.matrix
  }  , list( .lmean1 = lmean1, .lmean2 = lmean2,
             .emean1 = emean1, .emean2 = emean2,
             .lsd1   = lsd1  , .lsd2   = lsd2  , .lrho = lrho,
             .esd1   = esd1  , .esd2   = esd2  , .erho = erho ))),

  last = eval(substitute(expression({
    misc$link <-    c("mean1" = .lmean1,
                      "mean2" = .lmean2,
                      "sd1"   = .lsd1,
                      "sd2"   = .lsd2,
                      "rho"   = .lrho )

    misc$earg <- list("mean1" = .emean1,
                      "mean2" = .emean2, 
                      "sd1"   = .esd1,
                      "sd2"   = .esd2,
                      "rho"   = .erho )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }) , list( .lmean1 = lmean1, .lmean2 = lmean2,
             .emean1 = emean1, .emean2 = emean2,
             .lsd1   = lsd1  , .lsd2   = lsd2  , .lrho = lrho,
             .esd1   = esd1  , .esd2   = esd2  , .erho = erho ))),
  loglikelihood = eval(substitute(
          function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    mean1 <- eta2theta(eta[, 1], .lmean1, earg = .emean1)
    mean2 <- eta2theta(eta[, 2], .lmean2, earg = .emean2)
    sd1   <- eta2theta(eta[, 3], .lsd1  , earg = .esd1  )
    sd2   <- eta2theta(eta[, 4], .lsd2  , earg = .esd2  )
    Rho   <- eta2theta(eta[, 5], .lrho  , earg = .erho  )

    if (residuals) stop("loglikelihood residuals not ",
                          "implemented yet") else {
      sum(c(w) * dbinorm(x1 = y[, 1], x2 = y[, 2],
                         mean1 = mean1, mean2 = mean2,
                         var1 = sd1^2, var2 = sd2^2, cov12 = Rho *sd1*sd2,
                         log = TRUE))
    }
  } , list( .lmean1 = lmean1, .lmean2 = lmean2,
            .emean1 = emean1, .emean2 = emean2,
            .lsd1   = lsd1  , .lsd2   = lsd2  , .lrho = lrho,
            .esd1   = esd1  , .esd2   = esd2  , .erho = erho,
            .imethod = imethod ))),
  vfamily = c("binormal"),
  deriv = eval(substitute(expression({
    mean1 <- eta2theta(eta[, 1], .lmean1, earg = .emean1)
    mean2 <- eta2theta(eta[, 2], .lmean2, earg = .emean2)
    sd1   <- eta2theta(eta[, 3], .lsd1  , earg = .esd1  )
    sd2   <- eta2theta(eta[, 4], .lsd2  , earg = .esd2  )
    Rho   <- eta2theta(eta[, 5], .lrho  , earg = .erho  )

    zedd1 <- (y[, 1] - mean1) / sd1
    zedd2 <- (y[, 2] - mean2) / sd2
    temp5 <- 1 - Rho^2

    SigmaInv <- matrix(0, n, dimm(2))
    SigmaInv[, iam(1, 1, M = 2)] <- 1 / ((sd1^2) * temp5)
    SigmaInv[, iam(2, 2, M = 2)] <- 1 / ((sd2^2) * temp5)
    SigmaInv[, iam(1, 2, M = 2)] <- -Rho / (sd1 * sd2 * temp5)
    dl.dmeans <- mux22(t(SigmaInv), y - cbind(mean1, mean2), M = 2,
                       as.matrix = TRUE)
    dl.dsd1   <- -1 / sd1 + zedd1 * (zedd1 - Rho * zedd2) / (sd1 * temp5)
    dl.dsd2   <- -1 / sd2 + zedd2 * (zedd2 - Rho * zedd1) / (sd2 * temp5)
    dl.drho   <- -Rho * (zedd1^2 - 2 * Rho * zedd1 * zedd2 +
                        zedd2^2) / temp5^2 +
                zedd1 * zedd2 / temp5 +
                Rho / temp5

    dmean1.deta <- dtheta.deta(mean1, .lmean1) 
    dmean2.deta <- dtheta.deta(mean2, .lmean2) 
    dsd1.deta   <- dtheta.deta(sd1  , .lsd1  ) 
    dsd2.deta   <- dtheta.deta(sd2  , .lsd2  ) 
    drho.deta   <- dtheta.deta(Rho  , .lrho  ) 
    dthetas.detas  <- cbind(dmean1.deta,
                           dmean2.deta,
                           dsd1.deta,
                           dsd2.deta,
                           drho.deta)

    c(w) * cbind(dl.dmeans[, 1],
                 dl.dmeans[, 2],
                 dl.dsd1,
                 dl.dsd2,
                 dl.drho) * dthetas.detas
  }), list( .lmean1 = lmean1, .lmean2 = lmean2,
            .emean1 = emean1, .emean2 = emean2,
            .lsd1   = lsd1  , .lsd2   = lsd2  , .lrho = lrho,
            .esd1   = esd1  , .esd2   = esd2  , .erho = erho,
            .imethod = imethod ))),

  weight = eval(substitute(expression({
    wz <- matrix(0.0, n, dimm(M))
    wz[, iam(1, 1, M)] <- SigmaInv[, iam(1, 1, M = 2)]
    wz[, iam(2, 2, M)] <- SigmaInv[, iam(2, 2, M = 2)]
    wz[, iam(1, 2, M)] <- SigmaInv[, iam(1, 2, M = 2)]
    wz[, iam(3, 3, M)] <- (1 + 1 / temp5) / sd1^2
    wz[, iam(4, 4, M)] <- (1 + 1 / temp5) / sd2^2
    wz[, iam(3, 4, M)] <- -(Rho^2) / (temp5 * sd1 * sd2)
    wz[, iam(5, 5, M)] <- (1 + Rho^2) / temp5^2
    wz[, iam(3, 5, M)] <- -Rho / (sd1 * temp5)
    wz[, iam(4, 5, M)] <- -Rho / (sd2 * temp5)
    for (ilocal in 1:M)
      for (jlocal in ilocal:M)
        wz[, iam(ilocal, jlocal, M)] <- wz[, iam(ilocal, jlocal, M)] *
                                        dthetas.detas[, ilocal] *
                                        dthetas.detas[, jlocal]
      c(w) * wz
  }), list( .lmean1 = lmean1, .lmean2 = lmean2,
            .emean1 = emean1, .emean2 = emean2,
            .lsd1   = lsd1  , .lsd2   = lsd2  , .lrho = lrho,
            .esd1   = esd1  , .esd2   = esd2  , .erho = erho,
            .imethod = imethod ))))
}









gumbelI <-
  function(la = "identity", earg = list(), ia = NULL, imethod = 1) {

  la <- as.list(substitute(la))
  earg  <- link2list(la)
  la <- attr(earg, "function.name")



  if (length(ia) && !is.Numeric(ia, allowable.length = 1))
      stop("'ia' must be a single number")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2.5)
      stop("argument 'imethod' must be 1 or 2")


  new("vglmff",
  blurb = c("Gumbel's Type I Bivariate Distribution\n",
          "Links:    ",
          namesof("a", la, earg =  earg )),
  initialize = eval(substitute(expression({
    if (!is.matrix(y) || ncol(y) != 2)
        stop("the response must be a 2 column matrix") 

    if (any(y < 0))
        stop("the response must have non-negative values only")

    predictors.names <-
      c(namesof("a", .la, earg =  .earg , short = TRUE))
    if (!length(etastart)) {
        ainit  <- if (length( .ia ))  rep( .ia, len = n) else {
            mean1 <- if ( .imethod == 1) median(y[,1]) else mean(y[,1])
            mean2 <- if ( .imethod == 1) median(y[,2]) else mean(y[,2])
                Finit <- 0.01 + mean(y[,1] <= mean1 & y[,2] <= mean2)
      (log(Finit+expm1(-mean1)+exp(-mean2))+mean1+mean2)/(mean1*mean2)
            }
            etastart <-
               theta2eta(rep(ainit,  len = n), .la, earg =  .earg )
      }
  }), list( .ia=ia, .la = la, .earg = earg, .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha <- eta2theta(eta, .la, earg =  .earg )
    cbind(rep(1, len = length(alpha)),
          rep(1, len = length(alpha)))
  }, list( .la = la ))),
  last = eval(substitute(expression({
    misc$link <-    c("a" = .la )
    misc$earg <- list("a" = .earg )

    misc$expected <- FALSE
    misc$pooled.weight <- pooled.weight
  }), list( .la = la, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    alpha  <- eta2theta(eta, .la, earg =  .earg )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
    denom <- (alpha*y[,1] - 1) * (alpha*y[,2] - 1) + alpha
    mytolerance <- .Machine$double.xmin
    bad <- (denom <= mytolerance)   # Range violation
    if (any(bad)) {
      cat("There are some range violations in @deriv\n")
          flush.console()
          denom[bad] <- 2 * mytolerance
      }
      sum(w * (-y[,1] - y[,2] + alpha*y[,1]*y[,2] + log(denom)))
      }
  }, list( .la = la, .earg = earg ))),
  vfamily = c("gumbelI"),
  deriv = eval(substitute(expression({
      alpha  <- eta2theta(eta, .la, earg =  .earg )
      numerator <- (alpha*y[,1] - 1)*y[,2] + (alpha*y[,2] - 1)*y[,1] + 1
      denom <- (alpha*y[,1] - 1) * (alpha*y[,2] - 1) + alpha
      denom <- abs(denom)
      dl.dalpha <- numerator / denom + y[,1]*y[,2]
      dalpha.deta <- dtheta.deta(alpha,  .la, earg =  .earg )
      c(w) * cbind(dl.dalpha * dalpha.deta)
  }), list( .la = la, .earg = earg ))),
  weight = eval(substitute(expression({
    d2l.dalpha2 <- (numerator/denom)^2 - 2*y[,1]*y[,2] / denom
    d2alpha.deta2 <- d2theta.deta2(alpha, .la, earg =  .earg )
    wz <- w * (dalpha.deta^2 * d2l.dalpha2 - d2alpha.deta2 * dl.dalpha)
    if (TRUE &&
        intercept.only) {
        wz <- cbind(wz)
        sumw <- sum(w)
        for (iii in 1:ncol(wz))
            wz[,iii] <- sum(wz[,iii]) / sumw
        pooled.weight <- TRUE
        wz <- c(w) * wz   # Put back the weights
    } else
        pooled.weight <- FALSE
    wz
  }), list( .la = la, .earg = earg ))))
}




kendall.tau <- function(x, y, exact = FALSE, max.n = 3000) {

  if ((N <- length(x)) != length(y))
    stop("arguments 'x' and 'y' do not have equal lengths")

  NN <- if (!exact && N > max.n) {
    cindex <- sample.int(n = N, size = max.n, replace = FALSE)
    x <- x[cindex] 
    y <- y[cindex] 
    max.n
  } else {
    N
  }


  ans3 <-
    c(.C("VGAM_C_kend_tau", as.double(x), as.double(y),
         as.integer(NN), ans = double(3),
         NAOK = TRUE, DUP = TRUE, PACKAGE = "VGAM")$ans)

  con <- ans3[1] + ans3[2] / 2  # Ties put half and half
  dis <- ans3[3] + ans3[2] / 2
  (con - dis) / (con + dis)
}




if (FALSE)
kendall.tau <- function(x, y, exact = TRUE, max.n = 1000) {

  if ((N <- length(x)) != length(y))
    stop("arguments 'x' and 'y' do not have equal lengths")
  index <- iam(NA, NA, M = N, both = TRUE)

  index$row.index <- index$row.index[-(1:N)] 
  index$col.index <- index$col.index[-(1:N)] 

  NN <- if (!exact && N > max.n) {
    cindex <- sample.int(n = N, size = max.n, replace = FALSE)
    index$row.index <- index$row.index[cindex] 
    index$col.index <- index$col.index[cindex] 
    max.n
  } else{
    choose(N, 2)
  }

  con <- sum((x[index$row.index] - x[index$col.index]) *
             (y[index$row.index] - y[index$col.index]) > 0)
  dis <- NN - con
  (con - dis) / (con + dis)
}




dbistudenttcop <- function(x1, x2, df, rho = 0, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  u1 <- qt(x1, df = df)
  u2 <- qt(x2, df = df)

  logdensity <-
    -(df/2 + 1) * log1p(
    (u1^2 + u2^2 - 2 * rho * u1 * u2) / (df * (1 - rho^2))) -
    log(2*pi) - 0.5 * log1p(-rho^2) -
  dt(u1, df = df, log = TRUE) -
  dt(u2, df = df, log = TRUE)

  if (log.arg) logdensity else exp(logdensity)
}





